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

my $scriptname = basename($0);
my $version = "v5.13.031618";

# Remove when in prod.
#print "\n";
#print colored("*" x 75, 'bold yellow on_black'), "\n";
#print colored("      DEVELOPMENT VERSION OF MATCH_MOI_REPORT (ver: $version)\n", 
#    'bold yellow on_black');
#print colored("*" x 75, 'bold yellow on_black');
#print "\n\n";

my $help;
my $ver_info;
my $outfile;
my $freq_cutoff = 5;
my $cn_cutoff; # if set to 4, will use 5% CI.  Else will use CN as the threshold.  No need to specify.
my $cn_upper_cutoff = 4; # Configure to capture upper and lower bound CNs in an attempt to get both amps and dels
my $cn_lower_cutoff = 1; # Configure to capture upper and lower bound CNs in an attempt to get both amps and dels
my $read_count = 100;
my $raw_output;
my $nocall;
my $blood;
my $ped_match;

my $description = <<"EOT";
Program to parse an IR VCF file to generate a list of NCI-MATCH MOIs and aMOIs.  
This program requires the NCI-MATCH CNV Report, Fusion Report, IPC Report, and 
vcfExtractor scripts to be in your path prior to running.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <VCF>
    -f, --freq   INT   Don't report SNVs / Indels below this allele frequency 
                       (DEFAULT: $freq_cutoff%)
    --cu         INT   Set upper bound for amplifications (DEFAULT:
                       5% CI >= $cn_upper_cutoff) 
    --cl         INT   Set lower bound for deletions (DEFAULT:
                       95% CI <= $cn_lower_cutoff) 
    -c, --cn     INT   Don't report CNVs below this copy number threshold.  
                       DEFAULT: off. 
    -r, --reads  INT   Don't report Fusions below this read count. 
                       DEFAULT: $read_count reads.
    -o, --output STR   Send output to custom file.  Default is STDOUT.
    -b, --blood        Report is from a Pediatric MATCH Blood specimen and no 
                       fusion panel run.
    -p, --ped_match    Data is from the Pediatric MATCH study. Different 
                       non-hotspot rules are run, and so the reporting is a 
                       little different.
    -R, --Raw          Output raw data rather than pretty printed report that 
                       can be parsed with other tools
    -n, --nocall       Do not report NOCALL variants in Fusion and CNV space. 
                       Due to noise NOCALL is always on in SNV / Indel space.
    -v, --version      Version information
    -h, --help         Print this help information
EOT

GetOptions( "freq|f=f"      => \$freq_cutoff,
            "cn|c=i"        => \$cn_cutoff,
            "cu=i"          => \$cn_upper_cutoff,
            "cl=i"          => \$cn_lower_cutoff,
            "output|o=s"    => \$outfile,
            "Raw|R"         => \$raw_output,
            "blood|b"       => \$blood,
            "ped_match|p"   => \$ped_match,
            "reads|r=i"     => \$read_count,
            "nocall|n"      => \$nocall,
            "version|v"     => \$ver_info,
            "help|h"        => \$help )
        or die $usage;

sub help {
	printf "%s - %s\n\n%s\n\n%s\n", $scriptname, $version, $description, $usage;
	exit;
}

sub print_version {
	printf "%s - %s\n", $scriptname, $version;
	exit;
}

help if $help;
print_version if $ver_info;

# Need to disable the CI cutoffs if we want to use a raw CN cutoff like in 
# MATCH prod.
($cn_upper_cutoff, $cn_lower_cutoff) = '' if $cn_cutoff;

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
    open( $out_fh, ">", $outfile ) 
        || die "Can't open the output file '$outfile' for writing: $!";
} else {
    $out_fh = \*STDOUT;
}

# Now that we have more than one study with this assay, need to configure some 
# specific params.
my $study;
($ped_match or $blood) ? ($study = 'pediatric') : ($study = 'adult');

if (DEBUG) {
    print "============================  DEBUG  ============================\n";
    print "Params as passed into script:\n";
    print "\tCNV Threshold    => $cn_cutoff\n";
    print "\tVAF Threshold    => $freq_cutoff\n";
    print "\tFusion Threshold => $read_count\n";
    print "\tOutput File      => ";
    print "\tStudy            => $study\n";
    print "\tPediatric Blood  => $blood\n";
    ($outfile) ? print " $outfile\n" : print "\n";
    print "=================================================================\n\n";
}

# Check ENV for required programs
my @required_programs = qw( vcfExtractor.pl ocp_cnv_report.pl 
    ocp_control_summary.pl ocp_fusion_report.pl match_rna_qc.pl );
for my $prog (@required_programs) {
    die "ERROR: '$prog' is required, but not found in your path!\n" unless qx(which $prog);
}

########--------------------- END ARG Parsing ------------------------#########
my $vcf_file = shift;
die "ERROR: '$vcf_file' does not exist or is not a valid VCF file!\n" unless -e $vcf_file;

my $current_version = version->parse('2.3');
my $assay_version = version->parse( vcf_version_check(\$vcf_file) );
print "[INFO]: OVAT version: $assay_version\n" if DEBUG;

my $snv_indel_data          = proc_snv_indel(\$vcf_file);
my $cnv_data                = proc_cnv(\$vcf_file);
my $fusion_data             = proc_fusion(\$vcf_file) unless $blood;  

unless ($blood) {
    if ($assay_version >= $current_version) {
        my $rna_control_data = rna_qc(\$vcf_file);
        @$fusion_data{qw(P1_SUM P2_SUM)} = @$rna_control_data{qw(pool1_total pool2_total)};
    } else {
        $$fusion_data{'EXPR_CTRL'}  = proc_ipc(\$vcf_file, \$assay_version);
    }
}

# Print out the combined report
($raw_output) 
    ? raw_output($snv_indel_data,$fusion_data,$cnv_data) 
    : gen_report($vcf_file,$snv_indel_data,$fusion_data,$cnv_data);

sub vcf_version_check {
    # Check the VCF version, as well as, determining if we loaded a DNA only 
    # file and asked for both DNA and RNA data.
    my $vcf = shift;
    open(my $vcf_fh, "<", $$vcf);
    my @header = grep{ /^#/ } <$vcf_fh>;
    die "ERROR: The input file '$$vcf' does not appear to be a valid VCF file!\n" unless @header;
    die "ERROR: You have tried to load a VCF file witout fusion data and without ",
        "selecting the DNA only option!\n" unless grep {/Fusion/} @header or $blood;
    my ($ovat_ver) = map { /^##OncomineVariantAnnotationToolVersion=(\d+\.\d+)\.\d+/ } @header;
    return $ovat_ver;
}

sub proc_snv_indel {
    # use new VCF extractor to handle SNV and Indel calling
    my $vcf = shift;
    my %results;

    # Blacklisted variants based on Oncomine's list of recurrent SNPs and high 
    # error rate mutations.
    my @blacklisted_variants = qw(
        chr2:16082320:C:CCG
        chr2:209108317:C:T
        chr3:10183734:C:T
        chr4:1806131:T:C
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
        chr7:116411923:C:T
        chr9:98209628:T:TG
        chr9:139391437:TG:T
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

    open(my $vcf_data, "-|", "vcfExtractor.pl -Nna $$vcf") 
        or die "ERROR: can't parse VCF";

    # Check to see if we're using the dev version VCF extractor and issue 
    # warning.
    my $first_line = <$vcf_data>;
    if ($first_line =~ /[\*]+/) {
        warn "Using dev version of vcfExtractor!!\n";
        print $first_line;
        while ($. < 4) {
            my $line = <$vcf_data>;
            print $line;
        }
    }
    while (<$vcf_data>) {
        next unless /^chr/;
        my @fields = split;
        my $id = join(':', @fields[0..2]);
        next if grep {$id eq $_} @blacklisted_variants;

        # Map these variables to make typing easier and the code cleaner 
        # downstream
        my $vaf        = $fields[3];
        my $gene       = $fields[8];
        my $ocp_vc     = $fields[14];
        my $hotspot_id = $fields[7];
        my $aa_change  = $fields[11];
        my $function;
        ($fields[13] eq '---') 
            ? ($function = $fields[12] and $fields[13] = $fields[12]) 
            : ($function = $fields[13]);

        # Skip anything that does not map to an exon for now..might want to get 
        # utr vars later, though
        (my $exon = $fields[12]) =~ s/Exon// if $fields[12] =~ /^Exon/;
        $exon //= '-';
        
        #next unless $gene eq 'ERBB2';
        if ( $vaf >= $freq_cutoff ) {
            # Anything that's a hotspot
            if ( $hotspot_id ne '.' or $ocp_vc eq 'Hotspot' ) {
                $results{$id} = gen_var_entry(\@fields, 'Hotspot Variant');
            }
            # De Novo TSG frameshift calls
            elsif ( $ocp_vc eq 'Deleterious' ) {
                $results{$id} = gen_var_entry(\@fields, 'Deleterious in TSG');
            }
            # EGFR nonframeshiftDeletion and nonframeshiftInsertion in Exon 19, 
            # 20 rule for Arms A & C
            elsif ( $gene eq 'EGFR' ) { 
                next if $study eq 'pediatric';
                if ( $exon eq '19' && $function eq 'nonframeshiftDeletion' ) {
                    $results{$id} = gen_var_entry(\@fields, 
                        'EGFR in-frame deletion in Exon 19');
                }
                elsif ($exon eq '20' && $function eq 'nonframeshiftInsertion') {
                    $results{$id} = gen_var_entry(\@fields, 
                        'EGFR in-frame insertion in Exon 20');
                }
            }
            # ERBB2 nonframeshiftInsertion in Exon20 rule for Arm B
            elsif ( $gene eq 'ERBB2' 
                    && $exon eq '20' 
                    && $function eq 'nonframeshiftInsertion' 
                ) {
                next if $study eq 'pediatric';
                $results{$id} = gen_var_entry(\@fields, 
                    'ERBB2 in-frame insertion in Exon 20');
            }
            # KIT Exon 9 / 11 nonframeshiftInsertion and nonframeshiftDeletion 
            # rule for Arm V
            elsif ( $gene eq 'KIT' 
                    && (grep $exon eq $_, ('9','11','13','14')) 
                    && ($function =~ /nonframeshift.*/ || $function =~ /missense/) 
                ) {
                next if $study eq 'pediatric';
                $results{$id} = gen_var_entry(\@fields, 
                    'KIT in-frame indel in Exons 9, 11, 13, or 14');
            }
        }
    }
    return \%results;
}

sub gen_var_entry {
    my ($data, $rule) = @_;
    return [@$data[0..11,13,14],$rule];
}

sub proc_fusion {
    # Generate Hash of fusion data using ocp_fusion_report.pl
    my $vcf_file = shift;
    my %results;

    my $cmd;
    ($nocall) 
        ? ($cmd = qq(ocp_fusion_report.pl -Nn $$vcf_file)) 
        : ($cmd = qq(ocp_fusion_report.pl -n $$vcf_file));
    open(my $vcf_data, '-|', $cmd) 
        or die "ERROR: Can't parse VCF file for fusions!";
    while (<$vcf_data>) {
        # Skip header and blank lines
        next if $. < 4 || $_ =~ /^\s*$/;
        my @fields = split;

        # Unifying the fusion threshold for both inter and intra-genic 
        # fusions now.
        next if $fields[2] < $read_count;

        $fields[0] =~ s/\./|/;
        $results{"$fields[0]|$fields[1]"} = {
            'DRIVER'   => $fields[3],
            'PARTNER'  => $fields[4],
            'COUNT'    => $fields[2]
        };
    }
    # Get the total mapped RNA reads
    open(my $fh, "<", $$vcf_file);
    my ($mapped_reads) = map { /^##TotalMappedFusionPanelReads=(\d+)/ } <$fh>;
    close $fh;
    $results{'MAPPED_RNA'} = $mapped_reads;
    return \%results;
}

sub rna_qc {
    my $vcf = shift;
    my %results;
    open(my $rna_qc_data, "-|", "match_rna_qc.pl $$vcf");
    chomp(my $header = <$rna_qc_data>);
    chomp(my $data   = <$rna_qc_data>);
    @results{split(/,/,$header)} = split(/,/,$data);
    return \%results;
}

sub proc_ipc {
    # Get the RNA panel expression control sum 
    my ($vcf_file,$assay_version) = @_;
    my %panel = ( '2.0' => '1', '2.2' => '2', '2.3' => '3' );
    my $cmd = "ocp_control_summary.pl -m $panel{$$assay_version} $$vcf_file";
    open( my $vcf_data, '-|', $cmd ) 
        or die "ERROR: Can't parse the expression control data";
    # Trailing whitespace in ocp_control_summary output.
    my ($expr_sum) = map {/\s+(\d+)\s*$/} <$vcf_data>; 
    $expr_sum //= 0;
    return $expr_sum;
}

sub proc_cnv {
    my $vcf_file = shift;
    my %results;
    my ($gender, $cellularity, $mapd);

    my $cmd; 
    ($nocall) 
        ? ($cmd = qq(ocp_cnv_report.pl -N $$vcf_file)) 
        : ($cmd = qq(ocp_cnv_report.pl $$vcf_file));

    open(my $vcf_data, '-|', $cmd);
    while (<$vcf_data>) {
        if (/^:::/) {
            ($gender, $cellularity, $mapd) = $_ =~ /.*?Gender: (\w+), Cellularity: (.*?), MAPD: (\d\.\d+)\) :::$/;
        }
        elsif (/^chr/) {
            my @fields = split;
            my $ci_05 = $fields[8];
            my $ci_95 = $fields[9];
            my $cn = $fields[10];
            if ($cn_upper_cutoff && $cn_lower_cutoff) {
                if ($ci_05 >= $cn_upper_cutoff || $ci_95 <= $cn_lower_cutoff) {
                    $results{$fields[1]} = [@fields[0,5,8,10,9]];
                }
            } else {
                if ($cn_cutoff == 4) {
                    next unless $ci_05 >= $cn_cutoff;
                } else {
                    next unless $cn >= $cn_cutoff;
                }
                $results{$fields[1]} = [@fields[0,5,8,10,9]];
            }
        }
    }
    $results{'META'} = [$gender, $cellularity, $mapd];
    return \%results;
}

sub field_width {
    # Get the longest field width for formatting later.
    my ($data_ref, $type) = @_;

    my $ref_width = 0;
    my $alt_width = 0;
    my $cds_width = 0;
    my $aa_width  = 0;
    my $func_width = 0;

    my $fusion_id_width = 0;

    if ( $type eq 'snv' ) {
        while (my ($variant, $data) = each %$data_ref ) {
            my $ref_len = length( $data->[1] ) + 2;
            my $alt_len = length( $data->[2] ) + 2;
            my $cds_len = length( $data->[10] ) + 2;
            my $aa_len  = length( $data->[11] ) + 2;
            my $func_len = length( $data->[12] ) + 2;

            $ref_width = $ref_len if ( $ref_len > $ref_width );
            $alt_width = $alt_len if ( $alt_len > $alt_width );
            $cds_width = $cds_len if ( $cds_len > $cds_width );
            $aa_width = $aa_len if ( $aa_len > $aa_width );
            $func_width = $func_len if ( $func_len > $func_width );
        }
        return ( $ref_width, $alt_width, $cds_width, $aa_width, $func_width );
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
    # Kludgy way to simulate tee, but with and without colored output (ANSI 
    # escape codes in .txt files not helpful!
    my ($msg, $format) = @_;
    print {$out_fh} $msg if $outfile;
    ($format) ? print colored($msg, $format) : print $msg;
    return;
}

sub raw_output {
    # Generate a raw data dump so that we can import this data easily into 
    # another tool for further parsing.
    my ($snv_indels, $fusion_data, $cnv_data) = @_;
    my $mapd = $$cnv_data{'META'}[2];
    select $out_fh;

    for my $var (sort{ versioncmp( $a, $b ) } keys %$snv_indels) {
        print join(',', 'SNV', @{$$snv_indels{$var}}), "\n";
    }

    for my $var (sort{ versioncmp($$cnv_data{$a}->[0], $$cnv_data{$b}->[0])} keys %$cnv_data) {
        next if $var eq 'META';
        print join(',', 'CNV', $var, @{$$cnv_data{$var}}, $mapd), "\n";
    }

    # Do not output fusion results if blood since no fusion panel run.
    return if $blood;

    my @skipped_keys = qw(EXPR_CTRL MAPPED_RNA P1_SUM P2_SUM);
    for my $var ( sort { versioncmp( $a, $b ) } keys %$fusion_data ) {
        next if grep {$var eq $_} @skipped_keys;
        my ($fusion, $junct, $id) = split( /\|/, $var );
        print join(',', 'Fusion', "$fusion.$junct", $id, $$fusion_data{$var}->{'COUNT'}, 
            $$fusion_data{$var}->{'DRIVER'}, $$fusion_data{$var}->{'PARTNER'}), "\n";
    }
    return;
}

sub gen_report {
    # Print out the final MOI Report
    my ($vcf_filename, $snv_indels, $fusion_data, $cnv_data) = @_;
    my $format_string;

    #########################
    ##    Report Header    ##
    #########################
    my ($dna_name, $rna_name) = $vcf_filename =~ /^(?:.*\/)?(.*?)_v\d+_(.*?)_RNA_v\d+\.vcf/;
    my $study_title;
    ($study eq 'pediatric') 
        ? ($study_title = 'Pediatric NCI-MATCH') 
        : ($study_title = 'Adult NCI-MATCH');

    print_msg(sprintf("%s\n",'-'x150), 'bold ansi15');
    print_msg("$study_title MOI Report for ", 'bold ansi15');

    (! $dna_name || ! $rna_name) 
        ? print_msg("$vcf_filename\n", 'bold ansi15') 
        : print_msg("$dna_name DNA / $rna_name RNA\n", 'bold ansi15');

    print_msg(sprintf("%s\n",'-'x150), 'bold ansi15');

    #########################
    ## SNV / Indel Output  ##
    #########################
    print_msg("::: MATCH Reportable SNVs and Indels (VAF >= $freq_cutoff) :::\n",
        'ansi3');

    my ($refwidth, $altwidth, $cdswidth, $aawidth, 
        $funcwidth) = field_width( $snv_indels, 'snv' ) if $snv_indels;

    # Need some defaults if no data.
    $refwidth  //= 5; 
    $altwidth  //= 5;
    $cdswidth  //= 5; 
    $aawidth   //= 9; 
    $funcwidth //= 10;

    my %formatter = (
        'Chrom:Pos'            => "%-16s",
        'Ref'                  => "%-${refwidth}s",
        'Alt'                  => "%-${altwidth}s",
        'VAF'                  => '%-7s',
        'TotCov'               => '%-7s',
        'RefCov'               => '%-7s',
        'AltCov'               => '%-7s',
        'VARID'                => '%-12s',
        'Gene'                 => '%-8s',
        'Transcript'           => '%-16s',
        'CDS'                  => "%-${cdswidth}s",
        'Protein'              => "%-${aawidth}s",
        'Function'             => "%-${funcwidth}s",
        'oncomineVariantClass' => '%-21s',
        'FunctionalRule'       => '%s',
    );
    #dd \%formatter;
    #exit;

    my @snv_indel_header = qw( Chrom:Pos Ref Alt VAF TotCov RefCov AltCov VARID
        Gene Transcript CDS Protein Function oncomineVariantClass 
        FunctionalRule );

    my $snv_indel_format = join(' ', @formatter{@snv_indel_header}) . "\n";

    print_msg(sprintf($snv_indel_format, @snv_indel_header));
    if ( %$snv_indels ) {
        for my $variant ( sort{ versioncmp( $a, $b ) } keys %$snv_indels ) {
            print_msg(sprintf($snv_indel_format, @{$$snv_indels{$variant}}));
        }
    } else {
        print_msg(">>>>  No Reportable SNVs or Indels Found in Sample  <<<<\n", 
            "red on_black");
    }
    print_msg("\n");

    ########################
    ##  CNV Result Ouput  ##
    ########################
    my @formatted_mapd;
    my ($gender, $cellularity, $mapd) = @{$$cnv_data{'META'}};
    delete $$cnv_data{'META'};
    
    print_msg("::: MATCH Reportable CNVs (Gender: $gender, Cellularity: $cellularity, MAPD: ", 'ansi3');
    $format_string = format_string($mapd, '>', 0.5);
    print_msg(@$format_string);

    if ($cn_upper_cutoff) {
        print_msg( ", 5% CI >= $cn_upper_cutoff, 95% CI <= $cn_lower_cutoff) :::\n", "ansi3");
    } else {
        my $cnv_param_string; 
        ($cn_cutoff == 4) 
            ? ($cnv_param_string = ", 5% CI >=" ) 
            : ($cnv_param_string = ", CN >=");
        print_msg( "$cnv_param_string $cn_cutoff) :::\n", "ansi3");
    }

    my @cnv_header = qw( Chr Gene Tiles CI_05 CN CI_95 );
    print_msg(sprintf("%-9s %-10s %-6s %-10s %-10s %-10s\n", @cnv_header));

    if ( %$cnv_data ) {
        for my $cnv ( sort{ versioncmp( $$cnv_data{$a}->[0], $$cnv_data{$b}->[0] ) } keys %$cnv_data ) {
            print_msg(sprintf('%-9s %-10s %-6s %-10.2f ', 
                    $$cnv_data{$cnv}->[0], $cnv, @{$$cnv_data{$cnv}}[1,2]));
            my @formatted_copy_number = sprintf('%-10.2f ', $$cnv_data{$cnv}[3]);

            ($$cnv_data{$cnv}[3] < 1) 
                ? push(@formatted_copy_number,'bold red on_black') 
                : push(@formatted_copy_number,'bold green on_black');
            print_msg(@formatted_copy_number);
            print_msg(sprintf("%-10.2f\n", $$cnv_data{$cnv}[4]));
        }
    } else {
        print_msg(">>>>  No Reportable CNVs Found in Sample  <<<<\n",
            "red on_black");
    }
    print_msg("\n");

    # Do not output fusion data results if we have run on blood since no 
    # fusion panel run.
    return if $blood;

    #########################
    ##   Fusions Output    ##
    #########################
    my $tot_rna_reads = $$fusion_data{'MAPPED_RNA'};
    delete $$fusion_data{'MAPPED_RNA'};

    my $ipc_reads;
    if ($$fusion_data{'EXPR_CTRL'}) {
        $ipc_reads = $$fusion_data{'EXPR_CTRL'} and delete $$fusion_data{'EXPR_CTRL'};
    }
    
    my $pool1_sum;
    if ($$fusion_data{'P1_SUM'}) {
        $pool1_sum = int($$fusion_data{'P1_SUM'}) and delete $$fusion_data{'P1_SUM'};
    }

    my $pool2_sum;
    if ($$fusion_data{'P2_SUM'}) {
        $pool2_sum = int($$fusion_data{'P2_SUM'}) and delete $$fusion_data{'P2_SUM'};
    }

    print_msg("::: MATCH Reportable Fusions (Total Mapped Reads: ",'ansi3');
    # No good version information here.  Use the presence / absence of pool 
    # specific information to deduce what version and therefore which threshold 
    # to use.
    my $rna_reads_threshold;
    # if we have a pool1_sum value, then it's the v3 assay and we need a higher 
    # threshold.  If not, use the v2 data.
    ($pool1_sum) 
        ? ($rna_reads_threshold = 500000) 
        : ($rna_reads_threshold = 100000);
    $format_string = format_string($tot_rna_reads, '<', $rna_reads_threshold);
    print_msg(@$format_string);

    if ($ipc_reads) {
        print_msg('; Expression Control Sum: ','ansi3');
        $format_string = format_string($ipc_reads, '<', 20000);
        print_msg(@$format_string);

    } else {
        print_msg('; Pool1 Expression Reads: ','ansi3');
        $format_string = format_string($pool1_sum, '<', 100000);
        print_msg(@$format_string);

        print_msg('; Pool2 Expression Reads: ','ansi3');
        $format_string = format_string($pool2_sum, '<', 100000);
        print_msg(@$format_string);
    }

    $read_count = commify($read_count);
    print_msg( "; Threshold: $read_count) :::\n",'ansi3');

    my ($w1) = field_width( $fusion_data, 'fusion' );
    my $fusion_format = "%-${w1}s %-12s %-12s %-15s %-15s\n";
    my @fusion_header = qw( Fusion ID Read_Count Driver_Gene Partner_Gene );
    print_msg(sprintf($fusion_format, @fusion_header));
    if ( %$fusion_data ) {
        for ( sort { versioncmp( $a, $b ) } keys %$fusion_data ) {
            my ($fusion, $junct, $id) = split( /\|/ );
            print_msg(sprintf($fusion_format, "$fusion.$junct", $id, 
                    $$fusion_data{$_}->{'COUNT'}, $$fusion_data{$_}->{'DRIVER'}, 
                    $$fusion_data{$_}->{'PARTNER'}));
        }
    } else {
        print_msg(">>>>  No Reportable Fusions found in Sample  <<<<\n", "red on_black");
    }
    return;
}

sub format_string {
    my ($string, $cmp, $threshold) = @_;
    my $commified_string = commify($string);
    (eval "$string $cmp $threshold")
        ? return ["***$commified_string***", 'bold red on_black']
        : return [$commified_string,'ansi3'];
}

sub commify {
    my $number = shift;
    my ($integer, $decimal) = split(/\./, $number);
    my @groups = unpack '(A3)*', reverse $integer;
    my $commified_int = join(',', map {scalar reverse $_} reverse @groups);
    ($decimal) ? return "$commified_int.$decimal" : return $commified_int;
}


sub __exit__ {
    my ($line, $msg) = @_;
    print "\n\n";
    print colored("Got exit message at line $line with message: $msg", 
        'bold white on_green');
    print "\n";
    exit;
}
