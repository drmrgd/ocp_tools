#!/usr/bin/env perl
# Program to read OCAv3 VCF and perform a control summary analysis.  Requires 
# having an OCAv3 fusion panel JSON derived from the panel BED file (generally 
# using rna_panel_2_json.py), which is not freely distributed at this time.
#
# 2/27/2017 - D Sims
################################################################################
use warnings;
use strict;
use autodie;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Data::Dump;
use JSON;
use List::Util qw(sum);

my $scriptname = basename($0);
my $version = "v2.0.110718";
my $description = <<"EOT";
Read in one or more OCAv3 VCF files and output gene expression and expression 
control data.  Can either output a summary table of reads binned by type, or a 
raw CSV of all expression data.  Requires a panel JSON file not part of this 
package.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <VCF_file(s)>
    -r, --raw       Out raw gene expression and expression control data rather 
                    than a summary table.
    -j, --json      Use a custom panel JSON rather than the default OCAv3 file.  
                    Experimental at this time!
    -a, --all       Output the pool reads for all reads in a pool, not just the 
                    expression control reads.
    -o, --output    Send output to custom file.  Default is STDOUT.
    -v, --version   Version information
    -h, --help      Print this help information
EOT

my $help;
my $ver_info;
my $outfile;
my $raw_output;
my $custom_json;
my $all_pool_reads;

GetOptions( 
    "json|j=s"      => \$custom_json, 
    "raw|r"         => \$raw_output,
    "output|o=s"    => \$outfile,
    "version|v"     => \$ver_info,
    "all|a"         => \$all_pool_reads,
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
    print "ERROR: Not enough arguments passed to script!\n\n";
    print "$usage\n";
    exit 1;
}
my @vcfs = @ARGV;

# Test for and load fusion panel json
my $panel_json;
($custom_json) 
    ? ($panel_json = $custom_json) 
    : ($panel_json = dirname($0) . '/fusion_panel.json');

die "ERROR: Can not find the JSON file '$panel_json', that describes the panel",
    " contents. Must be present in script source dir!\n" unless -f $panel_json;

my $json_text = do {
    open(my $json_fh, "<", $panel_json);
    local $/;
    <$json_fh>;
};

my $json = JSON->new();
my $panel_data = $json->decode($json_text);

# Write output to either indicated file or STDOUT
my $out_fh;
if ( $outfile ) {
    print "Writing output to '$outfile'.\n";
	open( $out_fh, ">", $outfile ) 
       || die "Can't open the output file '$outfile' for writing: $!";
} else {
	$out_fh = \*STDOUT;
}
#########---------------------- END ARG Parsing -----------------------#########
my (%summary_data, %raw_data);
for my $vcf (@vcfs) {
    my ($control_data, $mapped_reads, $sample_name, $ec1, $ec2, $ge1, $ge2, 
        $fus1, $fus2) = read_vcf($vcf, $panel_data);

    $summary_data{$sample_name} = [$mapped_reads, $ec1, $ec2, $ge1, $ge2, $fus1,
        $fus2];
    $summary_data{$sample_name} = {
        'mapped_reads'       => $mapped_reads,
        'pool1_ec_reads'     => $ec1,
        'pool2_ec_reads'     => $ec2,
        'pool1_ge_reads'     => $ge1,
        'pool2_ge_reads'     => $ge2,
        'pool1_fusion_reads' => $fus1,
        'pool2_fusion_reads' => $fus2
    };
    $raw_data{$sample_name} = $control_data;
}

# Print out either the raw data (i.e. data by genes) or the summary data;
select $out_fh;
if ($raw_output) {
    print join(',', qw(sample_name type gene pool counts)), "\n";
    for my $sample (sort keys %raw_data) {
        for my $type (keys %{$raw_data{$sample}}) {
            for my $gene_info (@{$raw_data{$sample}->{$type}}) {
                my @data = split(/:/,$gene_info);
                print join(',', $sample,$type,@data), "\n";
            }
        }
    }
} else {
    my @fields = qw(mapped_reads pool1_ec_reads pool2_ec_reads pool1_ge_reads 
        pool2_ge_reads);
    print join(',', 'sample_name', @fields, 'pool1_total', 'pool2_total'), "\n";

    for my $sample (sort keys %summary_data) {
        my @pool1_vals = qw(pool1_ec_reads pool1_ge_reads);
        my @pool2_vals = qw(pool2_ec_reads pool2_ge_reads);

        if ($all_pool_reads) {
            push(@pool1_vals, 'pool1_fusion_reads');
            push(@pool2_vals, 'pool2_fusion_reads');
        }

        my $pool1_total = sum(@{$summary_data{$sample}}{@pool1_vals});
        my $pool2_total = sum(@{$summary_data{$sample}}{@pool2_vals});

        print join(',', $sample, @{$summary_data{$sample}}{@fields}, $pool1_total, 
            $pool2_total), "\n";
    }
}

sub read_vcf {
    my ($vcf, $panel) = @_;
    my ($mapped_reads, $sample_name, %control_data, %summary_data);

    open(my $vcf_fh, "<", $vcf);
    while (<$vcf_fh>) {
        $mapped_reads = $1 if (/^##TotalMappedFusionPanelReads=(\d+)/);
        if (/^#CHROM/) {
            my @header = split;
            $sample_name = $header[-1] and next;
        } else {
            next unless /SVTYPE=(ExprControl|Fusion|GeneExpression)/;
        }
        my ($gene, $counts, $type, $pool) = parse_line(trim_line($_), $panel);
        $pool =~ s/1,2/1-2/; # Get comma out of field for easier parsing.
        push(@{$control_data{$type}}, "$gene:$pool:$counts");
    }
    my ($ec1, $ec2) = get_sum(\%control_data, 'ExprControl');
    my ($ge1, $ge2) = get_sum(\%control_data, 'GeneExpression');
    my ($fus1, $fus2) = get_sum(\%control_data, 'Fusion');

    # Since we have two lines for each fusion, need to divide by 2. Probably 
    # need a better hash structure to store these data!
    $fus1 = $fus1/2;
    $fus2 = $fus2/2;

    return \%control_data, $mapped_reads, $sample_name, $ec1, $ec2, $ge1, $ge2, 
        $fus1, $fus2;
}

sub get_sum {
    my ($data, $type) = @_;
    my ($sum1, $sum2);
    
    for (@{$$data{$type}}) {
        my ($gene, $pool, $count) = split(/:/);
        if ($pool eq 'pool1') {
            $sum1 += $count;
        }
        elsif ($pool eq 'pool2') {
            $sum2 += $count;
        }
        elsif ($pool eq 'pool1-2') {
            $sum1 += ($count/2);
            $sum2 += ($count/2);
        }
    }
    return $sum1, $sum2;
}

sub trim_line {
    my $line = shift;
    my @elems = split(/\t/, $line);
    (my $assay_id = $elems[2]) =~ s/_[12]$//;
    return "ID=$assay_id;$elems[7]";
}

sub parse_line {
    my ($line, $panel) = @_;
    my ($gene, $assay_type, $counts) = $line =~ /^ID=(.*?);SVTYPE=(.*?);.*?READ_COUNT=(\d+).*/;
    my $pool = $$panel{$assay_type}->{$gene};
    return $gene, $counts, $assay_type, $pool;
}
