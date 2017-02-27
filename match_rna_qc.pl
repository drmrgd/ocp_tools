#!/usr/bin/perl
use warnings;
use strict;
use autodie;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Data::Dump;
use JSON;

my $scriptname = basename($0);
my $version = "v0.9.0_022717";
my $description = <<"EOT";
Read in an OCAv3 VCF file and output RNA panel control data.  
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <input_file>
    -r, --raw       Send out raw gene data rather than full summary dataset.
    -o, --output    Send output to custom file.  Default is STDOUT.
    -v, --version   Version information
    -h, --help      Print this help information
EOT

my $help;
my $ver_info;
my $outfile;
my $raw_output;

GetOptions( "raw|r"         => \$raw_output,
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
    print "ERROR: Not enough arguments passed to script!\n\n";
    print "$usage\n";
    exit 1;
}
my @vcfs = @ARGV;

#Test for and load fusion panel json
my $panel_json = dirname($0) . '/fusion_panel.json';
die "ERROR: Can not find the JSON file 'fusion_panel.json' that describes the panel contents.  Must be present in script source dir!\n" unless -f $panel_json;
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
	open( $out_fh, ">", $outfile ) || die "Can't open the output file '$outfile' for writing: $!";
} else {
	$out_fh = \*STDOUT;
}
#########------------------------------ END ARG Parsing ---------------------------------#########
my (%summary_data, %raw_data);
for my $vcf (@vcfs) {
    my ($control_data,$mapped_reads,$sample_name,$ec1,$ec2,$ge1,$ge2) = read_vcf($vcf, $panel_data);
    $summary_data{$sample_name} = [$mapped_reads, $ec1, $ec2, $ge1, $ge2];
    $raw_data{$sample_name} = $control_data;
}

# Print out either the raw data (i.e. data by genes) or the summary data;
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
    print join(',', qw(sample_name mapped_reads pool1_ec_reads pool2_ec_reads pool1_ge_reads pool2_ge_reads pool1_total pool2_total)), "\n";
    for my $sample (sort keys %summary_data) {
        my $total1 = $summary_data{$sample}->[1] + $summary_data{$sample}->[3];
        my $total2 = $summary_data{$sample}->[2] + $summary_data{$sample}->[4];
        print join(',',$sample,@{$summary_data{$sample}}, $total1, $total2), "\n";
    }
}

sub read_vcf {
    my ($vcf, $panel) = @_;
    my ($mapped_reads,$sample_name,%control_data,%summary_data);

    open(my $vcf_fh, "<", $vcf);
    while (<$vcf_fh>) {
        $mapped_reads = $1 if (/^##TotalMappedFusionPanelReads=(\d+)/);
        if (/^#CHROM/) {
            my @header = split;
            $sample_name = $header[-1] and next;
        } else {
            next unless /SVTYPE=.*?Expr.*/;
        }
        my ($gene,$counts,$type,$pool) = parse_line(trim_line($_), $panel);
        push(@{$control_data{$type}}, "$gene:$pool:$counts");
        }
    my ($ec1, $ec2) = get_sum(\%control_data, 'ExprControl');
    my ($ge1, $ge2) = get_sum(\%control_data, 'GeneExpression');
    return \%control_data,$mapped_reads,$sample_name,$ec1,$ec2,$ge1,$ge2;
}

sub get_sum {
    my ($data,$type) = @_;
    my ($sum1,$sum2);
    
    for (@{$$data{$type}}) {
        my ($gene,$pool,$count) = split(/:/,$_);
        if ($pool eq 'pool1') {
            $sum1 += $count;
        }
        elsif ($pool eq 'pool2') {
            $sum2 += $count;
        }
        elsif ($pool eq 'pool1,2') {
            $sum1 += ($count/2);
            $sum2 += ($count/2);
        }
    }
    return $sum1, $sum2;
}

sub trim_line {
    my $line = shift;
    my @elems = split(/\t/, $line);
    return "ID=$elems[2];$elems[7]";
}

sub parse_line {
    my ($line, $panel) = @_;
    my ($gene,$assay_type,$counts) = $line =~ /^ID=(.*?);SVTYPE=(.*?);.*?READ_COUNT=(\d+).*/;
    my $pool = $$panel{$assay_type}->{$gene};
    return $gene,$counts,$assay_type,$pool;
}
