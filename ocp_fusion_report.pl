#!/usr/bin/env perl
# Generate a fusion report from VCF files derived from and IR analysis.  Can 
# output either annotated fusions only or both annotated and novel fusions. Can
# also choose to output ref calls in addtion to variant calls to make a more 
# complete report.
#
# 6/9/2014 - D Sims
#######################################################################################################
use warnings;
use strict;
use autodie;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Parallel::ForkManager;
use Data::Dump;
use Sort::Versions;

my $scriptname = basename($0);
my $version = "v4.1.053119-dev";
my $description = <<"EOT";
Print out a summary table of fusions detected by the OCP Fusion Workflow VCF 
files. Can choose to output anything seen, or just limit to annotated fusions.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <vcf_file(s)>
    -R, --Ref       Include reference variants too (DEFAULT: OFF).
    -n, --novel     Include 'Non-Targeted' fusions in the output (DEFAULT: OFF).
    -t, --threshold Only report fusions above this threshold (DEFAULT: 25).
    -g, --gene      Only output data for a specific driver gene or genes separated by a comma. 
    -f, --forks     Number of processes to run at once (DEFAULT: 24).
    -N, --NOCALL    Don't report NOCALL or FAIL Fusions.
    -r, --raw       Raw output rather that pretty printed file.
    -o, --output    Write output to file.
    -v, --version   Display version information.
    -h, --help      Display this help text.
EOT

my $help;
my $ver_info;
my $novel_filter;
my $outfile;
my $ref_calls;
my $gene;
my $threshold = 25;
my $raw_output;
my $nocall;
my $num_procs = 24;

GetOptions( 
    "Ref|R"         => \$ref_calls,
    "novel|n"       => \$novel_filter,
    "threshold|t=s" => \$threshold,
    "gene|g=s"      => \$gene,
    "output|o=s"    => \$outfile,
    "raw|r"         => \$raw_output,
    "NOCALL|N"      => \$nocall,
    "forks|f=i"     => \$num_procs,
    "help|h"        => \$help,
    "version|v"     => \$ver_info,
);

sub help { 
    printf "%s - %s\n%s\n\n%s\n", $scriptname, $version, $description, $usage;
    exit;
}

sub version_info {
    printf "%s - %s\n", $scriptname, $version;
    exit;
}

help if $help;
version_info if $ver_info;

# Check we have some files to process
if ( @ARGV < 1 ) {
    print "ERROR: You must load at least one VCF file\n";
    print $usage;
    exit 1;
}

# Set up output
my $out_fh;
if ( $outfile ) {
    open( $out_fh, ">", $outfile ) 
        || die "Can't open the output file '$outfile' for writing: $!";
} else {
    $out_fh = \*STDOUT;
}

my $novel;
($novel_filter) ? ($novel = 1) : ($novel = 0);

my @files = @ARGV;
my @genes_list = map{uc($_)} split(/,/, $gene) if $gene;

#######======================  END ARG Parsing  ======================#######
my %results;
my $fwidth=0;

my $pm = Parallel::ForkManager->new($num_procs);
$pm->run_on_finish(
    sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, 
            $data_structure_reference) = @_;
        my $vcf  = $data_structure_reference->{input};
        my $name = $data_structure_reference->{id};
        $name //= basename($vcf);
        $results{$$name} = $data_structure_reference->{result};
    }
);

for my $input_file ( @files ) {
    $pm->start and next;
    my ($return_data, $sample_id) = proc_vcf(\$input_file);
    $pm->finish(0,
        {
            result => $return_data,
            input  => $input_file,
            id     => $sample_id,
        }
    );
}
$pm->wait_all_children;

# Generate and print out the final results table(s)
select $out_fh;
for my $sample ( sort keys %results ) {
    unless ($raw_output) {
        print "::: ";
        print join(',', @genes_list) if @genes_list;
        print " Fusions in $sample :::\n\n";
    }
    
    if ( $results{$sample} ) {
        my $fusion_format = "%-${fwidth}s %-12s %-12s %-15s %-15s\n";
        my @fusion_header = qw (Fusion ID Read_Count Driver_Gene Partner_Gene);
        printf $fusion_format, @fusion_header unless $raw_output;
        for my $entry (sort { versioncmp($a,$b) } keys %{$results{$sample}}) {
            if ($gene) {
                next unless grep{ $results{$sample}->{$entry}->{'DRIVER'} eq $_ } @genes_list;
            }
            my ($fusion, $junct, $id ) = split(/\|/, $entry);
            next if $results{$sample}->{$entry}->{'COUNT'} < $threshold && ! $ref_calls;
            print_data(\$sample, "$fusion.$junct", 
                \$id, $results{$sample}->{$entry}, \$fusion_format);
        }
    } else {
        print "\t\t\t<<< No Fusions Detected >>>\n\n" unless $raw_output;
    }
}

sub proc_vcf {
    my $vcf = shift;
    # Version 1,2, and 3 drivers. Not all exist in the current version, but 
    # keep all for backward compatibility.
    my @drivers = qw(ABL1 AKT2 AKT3 ALK AR AXL BRAF BRCA1 BRCA2 CDKN2A EGFR 
        ERBB2 ERBB4 ERG ESR1 ETV1 ETV1a ETV1b ETV4 ETV4a ETV5 ETV5a ETV5d FGFR1
        FGFR2 FGFR3 FGR FLT3 JAK2 KRAS MDM4 MET MYB MYBL1 NF1 NOTCH1 NOTCH4 NRG1
        NTRK1 NTRK2 NTRK3 NUTM1 PDGFRA PDGFRB PIK3CA PPARG PRKACA PRKACB PTEN 
        RAD51B RAF1 RB1 RELA RET ROS1 RSPO2 RSPO3 TERT
    );

    my %results;

    (my $sample_name = $$vcf) =~ s/(:?_Fusion_filtered)?\.vcf$//i;
    $sample_name =~ s/_RNA//;

    open( my $in_fh, "<", $$vcf );
    while (<$in_fh>) {
        next if /^#/;
        my @data = split;

        # Get rid of FAIL and NOCALL calls to be more compatible with MATCHBox 
        # output.
        next if $nocall && ($data[6] eq 'FAIL' || $data[6] eq 'NOCALL');

        if ( $data[7] =~ /SVTYPE=Fusion/ ) {
            my ( $name, $elem ) = $data[2] =~ /(.*?)_([12])$/;
            my ($count) = map { /READ_COUNT=(\d+)/ } @data;
            my ($pair, $junct, $id) = split(/\./, $name);
            $id //= '-';
            if ($id eq 'Non-Targeted' || $id eq 'Novel') {
                next unless $novel;
            }
            my ($gene1, $gene2) = split(/-/, $pair);

            # Filter out ref calls if we don't want to view them
            if ( $count == 0 ) { next unless ( $ref_calls ) }
            my $fid = join('|', $pair, $junct, $id);

            if ( $pair eq 'MET-MET' || $pair eq 'EGFR-EGFR' ) {
                $results{$fid}->{'DRIVER'} = $results{$fid}->{'PARTNER'} = $gene1;
            }
            elsif (grep {$_ eq $gene1} @drivers) {
                $results{$fid}->{'DRIVER'} = $gene1;
                $results{$fid}->{'PARTNER'} = $gene2;
            }
            elsif (grep {$_ eq $gene2} @drivers) {
                $results{$fid}->{'DRIVER'} = $gene2;
                $results{$fid}->{'PARTNER'} = $gene1;
            }
            else {
                $results{$fid}->{'DRIVER'} = 'UNKNOWN';
                $results{$fid}->{'PARTNER'} = "$gene1,$gene2";
            }
            $results{$fid}->{'COUNT'} = $count;

            # Get some field width formatting.
            $fwidth = (length($name)+4) if ( length($name) > $fwidth );
        }
    }
    return \%results, \$sample_name;
}

sub print_data {
    my ($sample_name,$fusion_name,$id,$data,$format) = @_;

    if ($raw_output) {
        print join(',', $$sample_name, $fusion_name, $$id, $$data{'COUNT'}, 
            $$data{'DRIVER'}), "\n";
    } else {
        printf $$format, $fusion_name, $$id, $$data{'COUNT'}, $$data{'DRIVER'}, 
        $$data{'PARTNER'};
    }
    return;
}
