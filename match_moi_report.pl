#!/usr/bin/perl
# Parse a VCF file and generate a table of MOIs and aMOIs
# TODO:
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

my $scriptname = basename($0);
my $version = "v2.5.0_070615";
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
	open( $out_fh, ">", $outfile ) || die "Can't open the output file '$outfile' for writing: $!";
} else {
	$out_fh = \*STDOUT;
}

if (DEBUG) {
    print "======================================  DEBUG  ======================================\n";
    print "Params as passed into script:\n";
    print "\tCNV Threshold  => $cn_cutoff\n";
    print "\tVAF Threshold  => $freq_cutoff\n";
    print "\tOutput File    => $out_fh\n";
    print "=====================================================================================\n\n";
}

########------------------------------ END ARG Parsing ---------------------------------#########
my $vcf_file = shift;
my ($gender, $cellularity, $mapd);
my $variant_data = read_vcf(\$vcf_file );

# Parse TSV file based on variant type and load up a hash of results.
my (%snv_indel_data, %fusion_data, %cnv_data);
for my $variant (@$variant_data) {
    my $type = $$variant{'rowtype'};
    given ($type) {
        when (/([ms]np|ins|del|complex)/) { proc_snv_indel( $variant ) }
        when (/Fusion/)                   { proc_fusion( $variant ) }
        when (/CNV/)                      { proc_cnv( $variant ) }
        default                           {next}
    }
}

# Print out the combined report
gen_report( \%snv_indel_data, \%fusion_data, \%cnv_data );

sub read_vcf {
    # Read in VCF file, proc with 'convert_vcf.py' and load a data struct.
    my $input_file = shift;
    (my $converted_vcf = $$input_file) =~ s/vcf$/tsv/;
    my @variant_results;

    # Want to have MAPD, Gender, and Cellularity in the output.  So, going to have to 
    # read the VCF twice it looks like
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
        }
    }
    close $vcf_fh;

    print "Coverting VCF file into TSV file for parsing...\n";
    qx( convert_vcf.py --force -i $$input_file -o $converted_vcf );
    print "Done!\n";

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
    # Pass in the raw data from the 'read_vcf()' funct and cull out what we don't need.
    my $raw_data = shift;
    my %filtered_data;

    # Use hash instead of array of keys in order to prevent undef key warning
    my %wanted_keys = ( 
    'rowtype'  => '',
    'call'    => '',
    'ALT'   => '',
    'CHROM'   => '',
    'FILTER'   => '',
    'ID' => '',
    'INFO...FR'   => '',
    'INFO...OALT'   => '',
    'INFO...OID'   => '',
    'INFO...OMAPALT'   => '',
    'INFO...OPOS'   => '',
    'INFO...OREF'   => '',
    'INFO...READ_COUNT'   => '',
    'INFO.1.FRO'   => '',
    'INFO.1.NUMTILES' => '',
    'INFO.A.AF'   => '',
    'INFO.A.FAO'   => '',
    'POS'   => '',
    'REF'   => '',
    'FORMAT.1.CN'   => '',
    'INFO...CI'   => '',
    'INFO.1.RO'   => '',
    'INFO.A.AO'   => '',
    'FUNC1.coding'   => '',
    'FUNC1.exon'   => '',
    'FUNC1.gene'   => '',
    'FUNC1.normalizedAlt'   => '',
    'FUNC1.normalizedPos'   => '',
    'FUNC1.normalizedRef'   => '',
    'FUNC1.oncomineGeneClass'   => '',
    'FUNC1.oncomineVariantClass'  => '',
    'FUNC1.protein'   => '',
    'FUNC1.transcript'   => '',
    'FUNC1.location'   => '',
    'FUNC1.function'  => '',
    );
    
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

    return if ( grep { $$variant_info{'FUNC1.location'} eq $_ } @undesired_locations );
    return if $$variant_info{'FUNC1.function'} eq 'synonymous';

    my $id = join( ':', $$variant_info{'CHROM'}, $$variant_info{'INFO...OPOS'}, $$variant_info{'INFO...OREF'}, $$variant_info{'INFO...OALT'} );
    
    # Remove Blacklisted SNPs
    return if grep { $id eq $_ } @blacklisted_variants;

    # Added to prevent missing long indel assembler calls.
    my $vaf;
    if ( $$variant_info{'INFO.A.AF'} eq '.'  ) {
        $vaf = $$variant_info{'INFO.A.AO'} / ($$variant_info{'INFO.1.RO'} + $$variant_info{'INFO.A.AO'});
    } else {
        $vaf = $$variant_info{'INFO.A.AF'};
    }
    
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
    $$variant_info{'FUNC1.oncomineVariantClass'} //= '---';

    if ( $vaf >= $freq_cutoff ) {
        # Anything that's a hotspot
        if ( $$variant_info{'INFO...OID'} ne '.' ) {
            # bin NOCALLs for now
            return if ( $$variant_info{'call'} eq 'NOCALL' );
            gen_var_entry( $variant_info, \$id );
        }
        # De Novo TSG frameshift calls
        elsif ( grep { $$variant_info{'FUNC1.oncomineVariantClass'} eq $_ } @oncomine_vc ) {
            gen_var_entry( $variant_info, \$id );
        }
        # EGFR nonframeshiftDeletion in Exon 19 rule
        # TODO: Get some test cases to try this.  
        elsif ( $$variant_info{'FUNC1.gene'} eq 'EGFR' 
            && $$variant_info{'FUNC1.exon'} eq '19' 
            # Jason saying that now we need Deletion, Insertion, and BlockSub here.
            && $$variant_info{'FUNC1.function'} =~ /nonframeshift.*/ )
        {
                gen_var_entry( $variant_info, \$id );
        }

        # ERBB2 nonframeshiftInsertion in Exon20 rule
        elsif ( $$variant_info{'FUNC1.gene'} eq 'ERBB2' 
            && $$variant_info{'FUNC1.exon'} eq '20' 
            && $$variant_info{'FUNC1.function'} eq 'nonframeshiftInsertion' ) 
            {
                gen_var_entry( $variant_info, \$id );
            }

        # TODO add in KIT Exon 9 Indel rule....check with jason for specifics
        # TODO: Get some test cases to try this.  
        elsif ( $$variant_info{'FUNC1.gene'} eq 'KIT' 
            && grep { $$variant_info{'FUNC1.exon'} == $_ } [9,11]
            && $$variant_info{'FUNC1.function'} =~ /nonframeshift(Deletion|Insertion)/ ) 
            {
                gen_var_entry( $variant_info, \$id );
            }
    }
    return;
}

sub gen_var_entry {
    my $data = shift;
    my $id = shift;
    
    my $coord = "$$data{'CHROM'}:$$data{'INFO...OPOS'}";
    my $ref = $$data{'INFO...OREF'};
    my $alt = $$data{'INFO...OALT'};
    my $filter = $$data{'FILTER'};
    (my $fr = $$data{'INFO...FR'}) =~ s/^\.,//;
    my $vaf = sprintf( "%0.2f", ($$data{'INFO.A.AF'} * 100) );
    my ($rcov, $acov, $tcov);
    if ( $$data{'INFO.A.FAO'} eq '.' ) {
        $rcov = $$data{'INFO.1.RO'};
        $acov = $$data{'INFO.A.AO'}; 
    } else {
        $rcov = $$data{'INFO.1.FRO'};
        $acov = $$data{'INFO.A.FAO'}; 
    }
    $tcov = $rcov + $acov;
    my $varid = $$data{'ID'};
    my $gene = $$data{'FUNC1.gene'};
    my ($om_gc, $om_vc);
    ($$data{'FUNC1.oncomineGeneClass'}) ? ($om_gc = $$data{'FUNC1.oncomineGeneClass'}) : ($om_gc = '---');
    ($$data{'FUNC1.oncomineVariantClass'}) ? ($om_vc = $$data{'FUNC1.oncomineVariantClass'}) : ($om_vc = '---');

    # Remove duplicate variant entries that span VCF blocks
    #if ($varid ne '.' && exists $snv_indel_data{$$id}) {

        #delete $snv_indel_data{$$id};
    #} 

    $snv_indel_data{$$id} = [$coord, $ref, $alt, $filter, $fr, $vaf, $tcov, $rcov, $acov, $varid, $gene, $om_gc, $om_vc];

    return;
}

sub proc_fusion {
    # Generate Hash of fusion data for output later.
    my $variant_info = shift;
    if ( $$variant_info{'call'} eq 'POS' ) { 
        my ($name, $elem) = $$variant_info{'ID'} =~ /(.*?)_([12])/;
        my ($pair, $junct, $id) = split( /\./, $name );
        $id = 'NOVEL' unless $id;

        # Get rid of Fusions that are below our thresholds
        if ( $id eq 'DelPositive' ) {
            return if $$variant_info{'INFO...READ_COUNT'} < 1000;
        } else {
            return if $$variant_info{'INFO...READ_COUNT'} < 25;
        }

        my $fid = join( '|', $pair, $junct, $id );
        
        $fusion_data{$fid}->{'COUNT'} = $$variant_info{'INFO...READ_COUNT'};
        if ( $elem == 2 ) {
            $fusion_data{$fid}->{'DRIVER'} = $$variant_info{'FUNC1.gene'};
        } else {
            $fusion_data{$fid}->{'PARTNER'} = $$variant_info{'FUNC1.gene'};
        }
    }
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
    my $filter_width= 0;
    my $fusion_id_width = 0;

    if ( $type eq 'snv' ) {
        for my $variant ( keys %$data_ref ) {
            my $ref_len = length( $$data_ref{$variant}[1] );
            my $alt_len = length( $$data_ref{$variant}[2] );
            my $filter_len = length( $$data_ref{$variant}[4] );
            $ref_width = $ref_len if ( $ref_len > $ref_width );
            $var_width = $alt_len if ( $alt_len > $var_width );
            $filter_width = $filter_len if ( $filter_len > $filter_width );
        }

        ( $filter_width > 13 ) ? ($filter_width += 4) : ($filter_width = 17);
        return ( $ref_width + 4, $var_width + 4, $filter_width);
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

sub gen_report {
    # Print out the final MOI Report
    my ($snv_indels, $fusion_data, $cnv_data) = @_;

    my ($w1, $w2, $w3);

    select $out_fh;

    print colored("::: MATCH Reportable SNVs and Indels (VAF >= $freq_cutoff) :::\n", "green on_black");
    ($w1, $w2, $w3) = field_width( $snv_indels, 'snv' );
    my $snv_indel_format = "%-17s %-${w1}s %-${w2}s %-10s %-${w3}s %-8s %-10s %-10s %-10s %-14s %-10s %-21s %s\n";
    my @snv_indel_header = qw( Chrom:Pos Ref Alt Filter Filter_Reason VAF TotCov RefCov AltCov VARID Gene oncomineGeneClass oncomineVariantClass );
    printf $snv_indel_format, @snv_indel_header;
    if ( %$snv_indels ) {
        for my $variant ( sort{ versioncmp( $a, $b ) } keys %$snv_indels ) {
            printf $snv_indel_format, @{$$snv_indels{$variant}};
        }
    } else {
        print colored(">>>>  No Reportable SNVs or Indels Found in Sample  <<<<\n", "red on_black");
    }
    print "\n";

    print colored("::: MATCH Reportable Fusions :::\n", "green on_black");
    ($w1) = field_width( $fusion_data, 'fusion' );
    my $fusion_format = "%-${w1}s %-12s %-12s %-15s %-15s\n";
    my @fusion_header = qw( Fusion ID Read_Count Driver_Gene Partner_Gene );
    printf $fusion_format, @fusion_header;
    if ( %$fusion_data ) {
        for ( sort { versioncmp( $a, $b ) } keys %$fusion_data ) {
            my ($fusion, $junct, $id) = split( /\|/ );
            printf $fusion_format, "$fusion.$junct", $id, $fusion_data{$_}->{'COUNT'}, $fusion_data{$_}->{'DRIVER'}, $fusion_data{$_}->{'PARTNER'};
        }
    } else {
        print colored(">>>>  No Reportable Fusions found in Sample  <<<<\n", "red on_black");
    }
    print "\n";

    print colored("::: MATCH Reportable CNVs (Gender: $gender, Cellularity: $cellularity, MAPD: $mapd, CN >= $cn_cutoff) :::\n", "green on_black");
    my $cnv_format = "%-9s %-10s %-6s %-10.3f %-10.1f %-10.3f\n";
    my @cnv_header = qw( Chr Gene Tiles CI_05 CN CI_95 );
    printf "%-9s %-10s %-6s %-10s %-10s %-10s\n", @cnv_header;
    if ( %$cnv_data ) {
        for my $cnv ( sort{ versioncmp( $cnv_data{$a}->[0], $cnv_data{$b}->[0] ) } keys %$cnv_data ) {
            printf $cnv_format, $cnv_data{$cnv}->[0], $cnv, @{$$cnv_data{$cnv}}[1..4];
        }
    } else {
        print colored(">>>>  No Reportable CNVs Found in Sample  <<<<\n", "red on_black");
    }
    return;
}
