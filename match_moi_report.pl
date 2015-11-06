#!/usr/bin/perl
# Parse a VCF file and generate a table of MOIs and aMOIs
#
# 2/12/2014 - D Sims
###################################################################################################
use warnings;
use strict;
use autodie;
use feature "switch";
# Need to turn of smartmatch warning, but only if the version of perl is newer
no if $] > 5.018, 'warnings', 'experimental::smartmatch'; 

use constant DEBUG => 0;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Sort::Versions;
use Term::ANSIColor;
use Data::Dump;

# Remove when in prod.
print "\n";
print colored("*" x 50, 'bold yellow on_black'), "\n";
print colored("      DEVELOPMENT VERSION OF MATCH_MOI_REPORT\n", 'bold yellow on_black');
print colored("*" x 50, 'bold yellow on_black');
print "\n\n";

my $scriptname = basename($0);
my $version = "v2.9.9_110615-dev";
my $description = <<"EOT";
Program to parse an IR VCF file to generate a list of NCI-MATCH MOIs and aMOIs.  This program requires 
the use of `convert_vcf.py` from ThermoFisher to run as it does the bulk of the file parsing.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <VCF>
    -f, --freq      Don't report SNVs / Indels below this allele frequency in decimal (0.1=10%). DEFAULT: 0.05 
    -c, --cn        Don't report CNVs below this copy number threshold.  DEFAULT: 7
    -o, --output    Send output to custom file.  Default is STDOUT.
    -v, --version   Version information
    -h, --help      Print this help information
EOT

my $help;
my $ver_info;
my $outfile;
my $freq_cutoff = 0.05;
my $cn_cutoff = 7;

GetOptions( "freq|f=f"      => \$freq_cutoff,
            "cn|c=i"        => \$cn_cutoff,
            "output|o=s"    => \$outfile,
            "version|v"     => \$ver_info,
            "help|h"        => \$help )
        or die $usage;

sub help {
	printf "%s - %s\n\n%s\n\n%s\n", $scriptname, $version, $description, $usage;
	exit;
}

sub version {
	printf "%s - %s\n", $scriptname, $version;
	exit;
}

help if $help;
version if $ver_info;

# Make sure enough args passed to script
if ( scalar( @ARGV ) < 1 ) {
    print "ERROR: No VCF file passed to script\n\n!";
    print "$usage\n";
    exit 1;
}

# Write output to either indicated file or STDOUT
my $out_fh;
if ( $outfile ) {
    #print "Writing results to $outfile...\n";
    open( $out_fh, ">", $outfile ) || die "Can't open the output file '$outfile' for writing: $!";
}

if (DEBUG) {
    print "======================================  DEBUG  ======================================\n";
    print "Params as passed into script:\n";
    print "\tCNV Threshold  => $cn_cutoff\n";
    print "\tVAF Threshold  => $freq_cutoff\n";
    print "\tOutput File    => ";
    ($outfile) ? print " $outfile\n" : print "\n";
    print "=====================================================================================\n\n";
}

########------------------------------ END ARG Parsing ---------------------------------#########
my $vcf_file = shift;
my ($gender, $cellularity, $mapd, $tot_rna_reads);
my $variant_data = read_vcf(\$vcf_file );

# Parse TSV file based on variant type and load up a hash of results.
my (%snv_indel_data, %fusion_data, %cnv_data, %ipc_data);
for my $variant (@$variant_data) {
    my $type = $$variant{'rowtype'};
    given ($type) {
        when (/([ms]np|ins|del|complex)/) { proc_snv_indel( $variant ) }
        when (/Fusion/)                   { proc_fusion( $variant ) }
        when (/CNV/)                      { proc_cnv( $variant ) }
        when (/ExprControl/)              { proc_ipc( $variant) }
        default                           {next}
    }
}

my $expr_control_sum = 0;
$expr_control_sum += $ipc_data{$_} for keys %ipc_data;
$fusion_data{'EXPR_CTRL'} = $expr_control_sum;

# Print out the combined report
gen_report( \%snv_indel_data, \%fusion_data, \%cnv_data );

sub read_vcf {
    # Read in VCF file, proc with 'convert_vcf.py' and load a data struct.
    my $input_file = shift;
    (my $converted_vcf = $$input_file) =~ s/vcf$/tsv/;
    my @variant_results;

    # Print out a nice title for the report based on the DNA and RNA sample name to make it nicer
    my ($dna_name, $rna_name) = $$input_file =~ /^(?:.*\/)?(.*?)_v\d+_(.*?)_RNA_v\d+\.vcf/;
    print_msg(sprintf("%s\n",'-'x150), 'bold ansi15');
    print_msg("NCI-MATCH MOI Report for ",'bold ansi15');
    (! $dna_name || ! $rna_name) ? print_msg("$$input_file\n",'bold ansi15') : print_msg("$dna_name DNA / $rna_name RNA\n",'bold ansi15'); 
    print_msg(sprintf("%s\n", '-'x150), 'bold ansi15');

    # Want to have MAPD, Gender, and Cellularity in the output.  So, going to have to read the VCF twice it looks like
    open( my $vcf_fh, "<", $$input_file );
    while (<$vcf_fh>) {
        if ( /^##/ ) {
            if ( /sampleGender=(\w+)/ ) {
                $gender = $1; 
                next;
            } 
            elsif ( /CellularityAsAFractionBetween0-1=(.*)$/ ) {
                $cellularity = sprintf( "%d%%", 100 * $1); 
                next;
            }
            elsif ( /mapd=(\d\.\d+)/ ) {
                $mapd = $1; 
                next;
            }
            elsif ( /TotalMappedFusionPanelReads=(\d+)/ ) {
                $tot_rna_reads = $1;
                next;
            }
        }
    }
    close $vcf_fh;

    #print "Coverting VCF file into TSV file for parsing...\n";
    qx( convert_vcf.py --force -i $$input_file -o $converted_vcf );

    open( my $fh, "<", $converted_vcf );
    chomp(my $header = <$fh>);
    # for some reason v1.3.1 of vcf_coverter is outputting quotes; remove them!
    my @header_elems = map { s/"//g; $_ } split(/\t/, $header);
    while (<$fh>) {
        $_ =~ s/"//g; 
        my %raw_data;
        chomp( my @elems = split(/\t/) );
        @raw_data{@header_elems} = @elems;
        my %filtered = filter_raw_data( \%raw_data );
        push(@variant_results, \%filtered);
    }
    close $fh;

    unlink $converted_vcf;
    return \@variant_results;
}

sub filter_raw_data {
    # Pass in the raw data from the 'read_vcf()' sub and cull out what we don't need.
    my $raw_data = shift;
    my %filtered_data;

    # Use hash instead of array of keys in order to prevent undef key warning
    my @wanted_vcf_elems = qw( rowtype call ALT CHROM FILTER ID INFO...FR INFO...OALT INFO...OID INFO...OMAPALT 
                               INFO...OPOS INFO...OREF INFO...READ_COUNT INFO.1.FRO INFO.1.NUMTILES INFO.A.AF 
                               INFO.A.FAO POS REF FORMAT.1.CN INFO...CI INFO.1.RO INFO.A.AO FUNC1.coding FUNC1.exon 
                               FUNC1.gene FUNC1.normalizedAlt FUNC1.normalizedRef FUNC1.normalizedPos 
                               FUNC1.oncomineGeneClass FUNC1.oncomineVariantClass FUNC1.protein FUNC1.transcript 
                               FUNC1.location FUNC1.function INFO.1.ANNOTATION );
    my %wanted_keys = map { $_ => '' } @wanted_vcf_elems;
    
    @filtered_data{keys %wanted_keys} = @$raw_data{keys %wanted_keys};
    return %filtered_data;
}

sub proc_snv_indel {
    # Follow rules to generate a list of MOIs and aMOIs
    my $variant_info = shift;
    my @oncomine_vc = qw( Deleterious Hotspot );
    my @undesired_locations = qw( unknown intergenic intronic utr_5 utr_3 splicesite_5 splicesite_3 upstream
                                  downstream intronic_nc ncRNA nonCoding);

    # Blacklisted variants based on Oncomine's list of recurrent SNPs and high error rate mutations.
    my @blacklisted_variants = qw(
        chr2:16082320:C:CCG
        chr2:209108317:C:T
        chr4:55593464:A:C
        chr4:55964925:G:A
        chr4:106196819:G:T
        chr5:112170746:AT:A
        chr5:112175651:A:G
        chr5:112175951:G:GA
        chr5:149449827:C:T
        chr7:116339642:G:T
        chr7:116340262:A:G
        chr7:116411990:C:T
        chr9:98209628:T:TG
        chr9:139391975:GC:G
        chr9:139399132:C:T
        chr10:89685288:T:TA
        chr10:123247514:C:CT
        chr10:123247514:CT:C
        chr11:320606:G:T
        chr11:108123551:C:T
        chr11:108138003:T:C
        chr11:108175462:G:A
        chr11:108175463:A:T
        chr13:28623587:C:T
        chr13:32972626:A:T
        chr16:68855966:G:A
        chr17:7578212:GA:G
        chr17:7579471:G:GC
        chr17:7579472:G:C
        chr17:29508455:TTA:T
        chr17:29553538:G:A
        chr19:1223125:C:G
        chr20:36030940:G:C
    );

    # Map these variables to make typing easier and the code cleaner downstream
    my $info_af  = $$variant_info{'INFO.A.AF'};
    my $ro       = $$variant_info{'INFO.1.RO'};
    my $ao       = $$variant_info{'INFO.A.AO'};
    my $function = $$variant_info{'FUNC1.function'};
    my $location = $$variant_info{'FUNC1.location'};
    my $exon     = $$variant_info{'FUNC1.exon'};
    my $gene     = $$variant_info{'FUNC1.gene'};
    my $ocp_vc   = $$variant_info{'FUNC1.oncomineVariantClass'};

    return if ( grep { $location eq $_ } @undesired_locations );
    return if $function eq 'synonymous';

    my $id = join( ':', $$variant_info{'CHROM'}, $$variant_info{'INFO...OPOS'}, $$variant_info{'INFO...OREF'}, $$variant_info{'INFO...OALT'} );
    
    # Added to prevent missing long indel assembler calls.
    my $vaf;
    if ( $info_af !~ /\d+/ ) { 
        $vaf = $ao / ($ro + $ao);
    } else {
        $vaf = $info_af;
    }

    # Substitute for the computed one to be used for the rest of the script
    $$variant_info{'VAF'} = $vaf;

    # Remove Blacklisted SNPs
    return if grep { $id eq $_ } @blacklisted_variants;

    # Get some debugging messages if there's an issue
    local $SIG{__WARN__} = sub {
        print "**** Error in parsing line ****\n\n";
        my $message = shift;
        print $message;
        print "variant id  => $id\n";
        print "vaf         => $vaf\n";
        dd $variant_info;
        print "*******************************\n";
        exit;
    };

    # Need to define the hash entry or we'll get an issue later on.
    $ocp_vc //= '---';

    if ( $vaf >= $freq_cutoff ) {
        # Anything that's a hotspot
        if ( $$variant_info{'INFO...OID'} ne '.' ) {
            # bin NOCALLs for now
            return if ( $$variant_info{'call'} eq 'NOCALL' );
            gen_var_entry( $variant_info, \$id, 'Hotspot Variant' );
        }
        # De Novo TSG frameshift calls
        elsif ( grep {$ocp_vc eq $_} @oncomine_vc ) {
            gen_var_entry( $variant_info, \$id, 'Deleterious in TSG' );
        }
        # EGFR nonframeshiftDeletion and nonframeshiftInsertion in Exon 19, 20 rule for Arms A & C
        elsif ( $gene eq 'EGFR' ) { 
            if ( $exon == 19 && $function eq 'nonframeshiftDeletion' ) {
                gen_var_entry( $variant_info, \$id, 'nonframeshiftDeletion in Exon 20' );
            }
            elsif ($exon == 20 && $function eq 'nonframeshiftInsertion') {
                gen_var_entry( $variant_info, \$id, 'nonframeshiftInsertion in Exon 19' );
            }
        }
        # ERBB2 nonframeshiftInsertion in Exon20 rule for Arm B
        elsif ( $gene eq 'ERBB2' && $exon eq '20' && $function eq 'nonframeshiftInsertion' ) {
            gen_var_entry( $variant_info, \$id, 'nonframeshiftInsertion in Exon 20' );
        }

        # KIT Exon 9 / 11 nonframeshiftInsertion and nonframeshiftDeletion rule for Arm V
        elsif ( $gene eq 'KIT' && (grep $exon == $_, (9,11)) && $function =~ /nonframeshift.*/ ) {
            gen_var_entry( $variant_info, \$id, 'nonframeshiftIndel in Exon 9 or 11' );
        }
    }
    return;
}

sub gen_var_entry {
    my $data = shift;
    my $id = shift;
    my $rule = shift;

    my $coord = "$$data{'CHROM'}:$$data{'INFO...OPOS'}";
    my $ref = $$data{'INFO...OREF'};
    my $alt = $$data{'INFO...OALT'};
    my $filter = $$data{'FILTER'};
    (my $fr = $$data{'INFO...FR'}) =~ s/^\.,//;
    my $vaf = sprintf( "%0.2f", ($$data{'VAF'} * 100) );

    my ($rcov, $acov, $tcov);
    if ( ! $$data{'INFO.A.FAO'} || ! $$data{'INFO.1.FRO'} ) {
        $rcov = $$data{'INFO.1.RO'};
        $acov = $$data{'INFO.A.AO'}; 
    } else {
        $rcov = $$data{'INFO.1.FRO'};
        $acov = $$data{'INFO.A.FAO'}; 
    }
    $tcov = $rcov + $acov;
    my $varid    = $$data{'ID'};
    my $gene     = $$data{'FUNC1.gene'};
    my $function = $$data{'FUNC1.function'};
    my $location = $$data{'FUNC1.location'};
    my $exon     = $$data{'FUNC1.exon'};
    my $hgvs     = $$data{'FUNC1.coding'};
    my $tscript  = $$data{'FUNC1.transcript'};
    my $protein  = $$data{'FUNC1.protein'};

    my ($om_gc, $om_vc);
    ($$data{'FUNC1.oncomineGeneClass'}) ? ($om_gc = $$data{'FUNC1.oncomineGeneClass'}) : ($om_gc = '---');
    ($$data{'FUNC1.oncomineVariantClass'}) ? ($om_vc = $$data{'FUNC1.oncomineVariantClass'}) : ($om_vc = '---');
    $snv_indel_data{$$id} = [$coord, $ref, $alt, $vaf, $tcov, $rcov, $acov, $varid, $gene, $tscript, $hgvs, $protein, $om_gc, $om_vc, $rule];
    return;
}

sub proc_fusion {
    # Generate Hash of fusion data for output later.
    my $variant_info = shift;
    my @drivers = qw( ABL1 AKT3 ALK AXL BRAF CDK4 EGFR ERBB2 ERG ETV1 ETV4 ETV5 FGFR1 FGFR2 FGFR3 NTRK1 NTRK3 PDGFRA PPARG
                      RAF1 RET ROS);

    if ( $$variant_info{'call'} eq 'POS' ) { 
        my ($name, $elem) = $$variant_info{'ID'} =~ /(.*)_([12])$/;
        my ($pair, $junct, $id) = split( /\./, $name );
        $id //= '-';

        #print "$$variant_info{'ID'}  => \n";
        #print "\tname:  $name\n";
        #print "\tpair:  $pair\n";
        #print "\tjunct: $junct\n";
        #print "\tID:    $$variant_info{'INFO.1.ANNOTATION'}\n";
        #print "\toID:   $id\n";
        #print '-'x50;
        #print "\n";

        # Get rid of Fusions that are below our thresholds
        if ( $id eq 'DelPositive' ) {
            return if $$variant_info{'INFO...READ_COUNT'} < 1000;
        } else {
            return if $$variant_info{'INFO...READ_COUNT'} < 25;
        }

        my $fid = join( '|', $pair, $junct, $id );
        
        if ( $pair eq 'MET-MET' || $pair eq 'EGFR-EGFR' ) {
            $fusion_data{$fid}->{'DRIVER'} = $fusion_data{$fid}->{'PARTNER'} = $$variant_info{'FUNC1.gene'};
        }
        elsif ( grep{ $_ eq $$variant_info{'FUNC1.gene'} } @drivers ) {
            $fusion_data{$fid}->{'DRIVER'} = $$variant_info{'FUNC1.gene'};
        } 
        else {
            $fusion_data{$fid}->{'PARTNER'} = $$variant_info{'FUNC1.gene'};
        }

        $fusion_data{$fid}->{'COUNT'} = $$variant_info{'INFO...READ_COUNT'};
    }
    return;
}

sub proc_ipc {
    # Get the RNA panel expression control data for output.
    my $variant_info = shift;
    $ipc_data{$$variant_info{'FUNC1.gene'}} = $$variant_info{'INFO...READ_COUNT'};
    return;
}

sub proc_cnv {
    my $variant_info = shift;

    if ( $$variant_info{'FORMAT.1.CN'} >= $cn_cutoff ) { 
        my ($ci_05, $ci_95) = $$variant_info{'INFO...CI'} =~ /0\.05:(.*?),0\.95:(.*)/;
        $cnv_data{$$variant_info{'FUNC1.gene'}} = [$$variant_info{'CHROM'}, $$variant_info{'INFO.1.NUMTILES'}, $ci_05, $$variant_info{'FORMAT.1.CN'}, $ci_95];
    }
    return;
}

sub field_width {
    # Get the longest field width for formatting later.
    my $data_ref = shift;
    my $type = shift;

    my $ref_width = 0;
    my $var_width = 0;
    my $hgvs_width = 0;
    my $filter_width= 0;
    my $fusion_id_width = 0;

    if ( $type eq 'snv' ) {
        for my $variant ( keys %$data_ref ) {
            my $ref_len = length( $$data_ref{$variant}[1] );
            my $alt_len = length( $$data_ref{$variant}[2] );
            my $filter_len = length( $$data_ref{$variant}[4] );
            my $hgvs_len = length( $$data_ref{$variant}[10] );
            $ref_width = $ref_len if ( $ref_len > $ref_width );
            $var_width = $alt_len if ( $alt_len > $var_width );
            $filter_width = $filter_len if ( $filter_len > $filter_width );
            $hgvs_width = $hgvs_len if ( $hgvs_len > $hgvs_width );
        }

        ( $filter_width > 13 ) ? ($filter_width += 4) : ($filter_width = 17);
        return ( $ref_width + 4, $var_width + 4, $filter_width, $hgvs_width + 4);
    } 
    elsif ( $type eq 'fusion' ) {
        for my $variant ( keys %$data_ref ) {
            my @id_elems = split( /\|/, $variant );
            my $length = length( "$id_elems[0].$id_elems[1]" );
            $fusion_id_width = $length if $length > $fusion_id_width;
        }
        return ( $fusion_id_width + 4 );
    }
    return;
}

sub print_msg {
    # Kludgy way to simulate tee, but with and without colored output (ANSI escape codes in .txt files not helpful!
    my ($msg, $format) = @_;

    # To outfile
    print {$out_fh} $msg if $outfile;

    # To STDOUT
    ($format) ? print colored($msg, $format) : print $msg;
    return;
}

sub gen_report {
    # Print out the final MOI Report
    my ($snv_indels, $fusion_data, $cnv_data) = @_;
    my ($w1, $w2, $w3, $w4);


    #########################
    ## SNV / Indel Output  ##
    #########################
    print_msg("::: MATCH Reportable SNVs and Indels (VAF >= $freq_cutoff) :::\n",'ansi3');
    ($w1, $w2, $w3, $w4) = field_width( $snv_indels, 'snv' );
    my @snv_indel_header = qw( Chrom:Pos Ref Alt VAF TotCov RefCov AltCov VARID Gene Transcript HGVS Protein oncomineGeneClass 
                               oncomineVariantClass Functional_Rule );
    my $snv_indel_format = "%-17s %-${w1}s %-${w2}s %-8s %-7s %-7s %-7s %-14s %-10s %-16s %-${w4}s %-16s %-21s %-22s %-21s\n";

    print_msg(sprintf($snv_indel_format, @snv_indel_header));
    if ( %$snv_indels ) {
        for my $variant ( sort{ versioncmp( $a, $b ) } keys %$snv_indels ) {
            print_msg(sprintf($snv_indel_format, @{$$snv_indels{$variant}}));
        }
    } else {
        print_msg(">>>>  No Reportable SNVs or Indels Found in Sample  <<<<\n", "red on_black");
    }
    print_msg("\n");


    ########################
    ##  CNV Result Ouput  ##
    ########################
    my @formatted_mapd;
    ($mapd >= 0.9) ? 
        (@formatted_mapd = ("**$mapd**", 'bold red on_black')) : 
        (@formatted_mapd = ($mapd,'ansi3'));

    print_msg("::: MATCH Reportable CNVs (Gender: $gender, Cellularity: $cellularity, MAPD: ", 'ansi3');
    print_msg(@formatted_mapd);
    print_msg( ", CN >= $cn_cutoff) :::\n", "ansi3");

    my $cnv_format = "%-9s %-10s %-6s %-10.3f %-10.1f %-10.3f\n";
    my @cnv_header = qw( Chr Gene Tiles CI_05 CN CI_95 );
    print_msg(sprintf("%-9s %-10s %-6s %-10s %-10s %-10s\n", @cnv_header));
    if ( %$cnv_data ) {
        for my $cnv ( sort{ versioncmp( $cnv_data{$a}->[0], $cnv_data{$b}->[0] ) } keys %$cnv_data ) {
            print_msg(sprintf($cnv_format, $cnv_data{$cnv}->[0], $cnv, @{$$cnv_data{$cnv}}[1..4]));
        }
    } else {
        print_msg(">>>>  No Reportable CNVs Found in Sample  <<<<\n", "red on_black");
    }
    print_msg("\n");

    #########################
    ##   Fusions Output    ##
    #########################
    my @read_count;
    ($tot_rna_reads < 100000) ? 
        (@read_count = ("**$tot_rna_reads**", 'bold red on_black')) : 
        (@read_count = ($tot_rna_reads,'ansi3')); 

    my $ipc_reads = $$fusion_data{'EXPR_CTRL'};
    delete $$fusion_data{'EXPR_CTRL'};

    my @ipc_output;
    ($ipc_reads < 20000) ? 
        (@ipc_output = ("**$ipc_reads**", 'bold red on_black')) : 
        (@ipc_output = ($ipc_reads, 'ansi3'));

    print_msg("::: MATCH Reportable Fusions (Total Reads: ",'ansi3');
    print_msg(@read_count);
    print_msg(', Sum Expression Control Reads: ','ansi3');
    print_msg(@ipc_output);
    print_msg( ") :::\n",'ansi3');

    ($w1) = field_width( $fusion_data, 'fusion' );
    my $fusion_format = "%-${w1}s %-12s %-12s %-15s %-15s\n";
    my @fusion_header = qw( Fusion ID Read_Count Driver_Gene Partner_Gene );
    print_msg(sprintf($fusion_format, @fusion_header));
    if ( %$fusion_data ) {
        for ( sort { versioncmp( $a, $b ) } keys %$fusion_data ) {
            my ($fusion, $junct, $id) = split( /\|/ );
            print_msg(sprintf($fusion_format, "$fusion.$junct", $id, $fusion_data{$_}->{'COUNT'}, $fusion_data{$_}->{'DRIVER'}, $fusion_data{$_}->{'PARTNER'}));
        }
    } else {
        print_msg(">>>>  No Reportable Fusions found in Sample  <<<<\n", "red on_black");
    }
    return;
}
