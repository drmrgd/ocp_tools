#!/usr/bin/perl
# Parse a VCF file and generate a table of MOIs and aMOIs
#
# 2/12/2014 - D Sims
###################################################################################################
use warnings;
use strict;
use autodie;

use constant DEBUG => 0;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Sort::Versions;
use Term::ANSIColor;
use Data::Dump;

# Remove when in prod.
#print "\n";
#print colored("*" x 50, 'bold yellow on_black'), "\n";
#print colored("      DEVELOPMENT VERSION OF MATCH_MOI_REPORT\n", 'bold yellow on_black');
#print colored("      WARN:  \n\tCNV Threshold set to 5% CI >= 4;\n\tMETe14 threshold set to >1000\n", 'bold yellow on_black');
#print colored("*" x 50, 'bold yellow on_black');
#print "\n\n";

my $scriptname = basename($0);
my $version = "v4.0.0_031016";
my $description = <<"EOT";
Program to parse an IR VCF file to generate a list of NCI-MATCH MOIs and aMOIs.  This program requires 
the NCI-MATCH CNV Report, Fusion Report, IPC Report, and vcfExtractor scripts to be in your path prior to running.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <VCF>
    -f, --freq   INT   Don't report SNVs / Indels below this allele frequency INT (DEFAULT: 5%)
    -c, --cn     INT   Don't report CNVs below this copy number threshold.  DEFAULT: 5% CI >= 4
    -o, --output STR   Send output to custom file.  Default is STDOUT.
    -r, --raw          Output raw data rather than pretty printed report that can be parsed with other tools
    -v, --version      Version information
    -h, --help         Print this help information
EOT

my $help;
my $ver_info;
my $outfile;
my $freq_cutoff = 5;
my $cn_cutoff = 7;
my $raw_output;
my $nofilter; # Retain filtered out calls.

GetOptions( "freq|f=f"      => \$freq_cutoff,
            "cn|c=i"        => \$cn_cutoff,
            "output|o=s"    => \$outfile,
            "raw|r"         => \$raw_output,
            "n|nofilter"    => \$nofilter,
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
    print "ERROR: No VCF file passed to script!\n\n";
    print "$usage\n";
    exit 1;
}

# Write output to either indicated file or STDOUT
my $out_fh;
if ( $outfile ) {
    print "Writing results to $outfile...\n";
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

# Check ENV for required programs
my @required_programs = qw( vcfExtractor.pl ocp_cnv_report.pl ocp_control_summary.pl ocp_fusion_report.pl );
for my $prog (@required_programs) {
    die "ERROR: '$prog' is required, but not found in your path!\n" unless `which $prog`;
}

########------------------------------ END ARG Parsing ---------------------------------#########
my $vcf_file = shift;

my $snv_indel_data          = proc_snv_indel(\$vcf_file);
my $cnv_data                = proc_cnv(\$vcf_file);
my $fusion_data             = proc_fusion(\$vcf_file);
$$fusion_data{'EXPR_CTRL'}  = proc_ipc(\$vcf_file);

# Print out the combined report
($raw_output) ? raw_output( $snv_indel_data, $fusion_data, $cnv_data ) : gen_report( $vcf_file, $snv_indel_data, $fusion_data, $cnv_data );

sub proc_snv_indel {
    # use new VCF extractor to handle SNV and Indel calling
    my $vcf = shift;
    my %results;

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
        chr9:139391438:TG:T
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
        chr17:7578373:T:C
        chr17:7579471:G:GC
        chr17:7579472:G:C
        chr17:29508455:TTA:T
        chr17:29553538:G:A
        chr19:1223125:C:G
        chr20:36030940:G:C
    );
    my @oncomine_vc = qw( Deleterious Hotspot );

    open(my $vcf_data, "-|", "vcfExtractor.pl -Nna $$vcf") or die "ERROR: can't parse VCF";
    while (<$vcf_data>) {
        next unless /^chr/;
        my @fields = split;
        my $id = join(':', @fields[0..2]);

        # Output blacklisted variants but with a flag if we want to.
        # XXXX
        #next if grep {$id eq $_} @blacklisted_variants;
        if ($nofilter) {
            if (grep {$id eq $_} @blacklisted_variants) {
                $id .= "_BL";
            }
        } else {
            next if grep {$id eq $_} @blacklisted_variants;
        }

        # Map these variables to make typing easier and the code cleaner downstream
        my $vaf        = $fields[3];
        my $gene       = $fields[8];
        my $ocp_vc     = $fields[15];
        my $hotspot_id = $fields[7];
        my $function   = $fields[13];
        (my $exon = $fields[12]) =~ s/Exon//;
        next if $exon eq 'intronic';
        
        if ( $vaf >= $freq_cutoff ) {
            # Anything that's a hotspot
            if ( $hotspot_id ne '.' ) {
                $results{$id} = gen_var_entry(\@fields, 'Hotspot Variant');
            }
            # De Novo TSG frameshift calls
            elsif ( grep {$ocp_vc eq $_} @oncomine_vc ) {
                $results{$id} = gen_var_entry(\@fields, 'Deleterious in TSG');
            }
            # EGFR nonframeshiftDeletion and nonframeshiftInsertion in Exon 19, 20 rule for Arms A & C
            elsif ( $gene eq 'EGFR' ) { 
                if ( $exon == 19 && $function eq 'nonframeshiftDeletion' ) {
                    $results{$id} = gen_var_entry(\@fields, 'nonframeshiftDeletion in Exon 20');
                }
                elsif ($exon == 20 && $function eq 'nonframeshiftInsertion') {
                    $results{$id} = gen_var_entry(\@fields, 'nonframeshiftInsertion in Exon 19');
                }
            }
            # ERBB2 nonframeshiftInsertion in Exon20 rule for Arm B
            elsif ( $gene eq 'ERBB2' && $exon eq '20' && $function eq 'nonframeshiftInsertion' ) {
                $results{$id} = gen_var_entry(\@fields, 'nonframeshiftInsertion in Exon 20');
            }
            # KIT Exon 9 / 11 nonframeshiftInsertion and nonframeshiftDeletion rule for Arm V
            elsif ( $gene eq 'KIT' && (grep $exon == $_, (9,11)) && $function =~ /nonframeshift.*/ ) {
                $results{$id} = gen_var_entry(\@fields, 'nonframeshiftIndel in Exon 9 or 11 of KIT');
            }
        }
    }
    #dd \%results;
    #exit;
    return \%results;
}
sub gen_var_entry {
    my ($data, $rule) = @_;
    return [@$data[0..11,13,14,15],$rule];
}

sub proc_fusion {
    # Generate Hash of fusion data using ocp_fusion_report.pl
    my $vcf_file = shift;
    my %results;

    open(my $vcf_data, '-|', "ocp_fusion_report.pl $$vcf_file") or die "ERROR: Can't parse VCF file for fusions!";
    while (<$vcf_data>) {
        # Skip header and blank lines
        next if $. < 4 || $_ =~ /^\s*$/;
        my @fields = split;

        # Get rid of Fusions that are below our thresholds
        if ( $fields[1] eq 'DelPositive' || $fields[0] eq 'MET-MET.M13M15' ) {
            next if $fields[2] < 1000;
        } else {
            next if $fields[2] < 25;
        }

        #my $fid = join( '|', $pair, $junct, $id );
        $fields[0] =~ s/\./|/;
        $results{"$fields[0]|$fields[1]"} = {
            'DRIVER'   => $fields[3],
            'PARTNER'  => $fields[4],
            'COUNT'    => $fields[2]
        };
    }
    # Get the total mapped RNA reads
    open(my $fh, "<", $$vcf_file);
    my ($mapped_reads) = map { /^##TotalMappedFusionPanelReads=(\d+)/ }<$fh>;
    close $fh;
    $results{'MAPPED_RNA'} = $mapped_reads;

    return \%results;
}

sub proc_ipc {
    # Get the RNA panel expression control sum 
    my $vcf_file = shift;

    open( my $vcf_data, '-|', "ocp_control_summary.pl -f $$vcf_file" ) || die "ERROR: Can't parse the IPC data";
    my ($expr_sum) = map {/\s+(\d+)\s*$/} <$vcf_data>; # trailing whitespace in ocp_control_summary output.
    $expr_sum //= 0;
    return $expr_sum;
}

sub proc_cnv {
    my $vcf_file = shift;
    my %results;
    my ($gender, $cellularity, $mapd);

    open(my $vcf_data, '-|', "ocp_cnv_report.pl $$vcf_file");
    while (<$vcf_data>) {
        if (/^:::/) {
            ($gender, $cellularity, $mapd) = $_ =~ /.*?Gender: (\w+), Cellularity: (.*?), MAPD: (\d\.\d+)\) :::$/;
        }
        elsif (/^chr/) {
            my @fields = split;
            my $ci_05 = $fields[8];
            my $cn = $fields[10];
            # XXX: Set CNV cutoff here with either 5% CI val and threshold or CN val and threshold
            #next unless $ci_05 >= $cn_cutoff;
            next unless $cn >= $cn_cutoff;
            $results{$fields[1]} = [@fields[0,5,8,10,9]];
        }
    }
    $results{'META'} = [$gender, $cellularity, $mapd];
    return \%results;
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
    print {$out_fh} $msg if $outfile;
    ($format) ? print colored($msg, $format) : print $msg;
    return;
}

sub raw_output {
    # Generate a raw data dump so that we can import this data easily into another tool for further parsing.
    my ($snv_indels, $fusion_data, $cnv_data) = @_;
    my $mapd = $$cnv_data{'META'}[2];

    for my $var (sort{ versioncmp( $a, $b ) } keys %$snv_indels) {
        print join(',', 'SNV', @{$$snv_indels{$var}}), "\n";
    }

    for my $var (sort{ versioncmp($$cnv_data{$a}->[0], $$cnv_data{$b}->[0])} keys %$cnv_data) {
        next if $var eq 'META';
        print join(',', 'CNV', $var, @{$$cnv_data{$var}}, $mapd), "\n";
    }

    for my $var ( sort { versioncmp( $a, $b ) } keys %$fusion_data ) {
        next if $var eq 'EXPR_CTRL' || $var eq 'MAPPED_RNA';
        my ($fusion, $junct, $id) = split( /\|/, $var );
        print join(',', 'Fusion', "$fusion.$junct", $id, $$fusion_data{$var}->{'COUNT'}, 
            $$fusion_data{$var}->{'DRIVER'}, $$fusion_data{$var}->{'PARTNER'}), "\n";
    }
    return;
}

sub gen_report {
    # Print out the final MOI Report
    my ($vcf_filename, $snv_indels, $fusion_data, $cnv_data) = @_;
    my ($w1, $w2, $w3, $w4);

    #########################
    ##    Report Header    ##
    #########################
    my ($dna_name, $rna_name) = $vcf_filename =~ /^(?:.*\/)?(.*?)_v\d+_(.*?)_RNA_v\d+\.vcf/;
    print_msg(sprintf("%s\n",'-'x150), 'bold ansi15');
    print_msg('NCI-MATCH MOI Report for ', 'bold ansi15');
    (! $dna_name || ! $rna_name) ? print_msg("$vcf_filename\n", 'bold ansi15') : print_msg("$dna_name DNA / $rna_name RNA\n", 'bold ansi15');
    print_msg(sprintf("%s\n",'-'x150), 'bold ansi15');

    #########################
    ## SNV / Indel Output  ##
    #########################
    print_msg("::: MATCH Reportable SNVs and Indels (VAF >= $freq_cutoff) :::\n",'ansi3');
    ($w1, $w2, $w3, $w4) = field_width( $snv_indels, 'snv' );
    my @snv_indel_header = qw( Chrom:Pos Ref Alt VAF TotCov RefCov AltCov VARID Gene Transcript HGVS Protein Function oncomineGeneClass 
                               oncomineVariantClass Functional_Rule );
    my $snv_indel_format = "%-17s %-${w1}s %-${w2}s %-8s %-7s %-7s %-7s %-14s %-10s %-16s %-${w4}s %-16s %-23s %-21s %-22s %-21s\n";

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
    my ($gender, $cellularity, $mapd) = @{$$cnv_data{'META'}};
    delete $$cnv_data{'META'};
    
    ($mapd >= 0.9 || $mapd == 0) ? 
        (@formatted_mapd = ("**$mapd**", 'bold red on_black')) : 
        (@formatted_mapd = ($mapd,'ansi3'));

    print_msg("::: MATCH Reportable CNVs (Gender: $gender, Cellularity: $cellularity, MAPD: ", 'ansi3');
    print_msg(@formatted_mapd);
    print_msg( ", 5% CI>= $cn_cutoff) :::\n", "ansi3");

    my $cnv_format = "%-9s %-10s %-6s %-10.3f %-10.1f %-10.3f\n";
    my @cnv_header = qw( Chr Gene Tiles CI_05 CN CI_95 );
    print_msg(sprintf("%-9s %-10s %-6s %-10s %-10s %-10s\n", @cnv_header));
    if ( %$cnv_data ) {
        for my $cnv ( sort{ versioncmp( $$cnv_data{$a}->[0], $$cnv_data{$b}->[0] ) } keys %$cnv_data ) {
            print_msg(sprintf($cnv_format, $$cnv_data{$cnv}->[0], $cnv, @{$$cnv_data{$cnv}}[1..4]));
        }
    } else {
        print_msg(">>>>  No Reportable CNVs Found in Sample  <<<<\n", "red on_black");
    }
    print_msg("\n");

    #########################
    ##   Fusions Output    ##
    #########################
    my $ipc_reads = $$fusion_data{'EXPR_CTRL'};  
    delete $$fusion_data{'EXPR_CTRL'};
    my $tot_rna_reads = $$fusion_data{'MAPPED_RNA'};
    delete $$fusion_data{'MAPPED_RNA'};

    my @read_count;
    ($tot_rna_reads < 100000) ? 
        (@read_count = ("**$tot_rna_reads**", 'bold red on_black')) : 
        (@read_count = ($tot_rna_reads,'ansi3')); 

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
            print_msg(sprintf($fusion_format, "$fusion.$junct", $id, $$fusion_data{$_}->{'COUNT'}, $$fusion_data{$_}->{'DRIVER'}, $$fusion_data{$_}->{'PARTNER'}));
        }
    } else {
        print_msg(">>>>  No Reportable Fusions found in Sample  <<<<\n", "red on_black");
    }
    return;
}
