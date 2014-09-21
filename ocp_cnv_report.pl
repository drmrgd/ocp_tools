#!/usr/bin/perl
# Read in VCF file from IR and output called CNVs
# 9/20/2014 - D Sims
#################################################################################################
use warnings;
use strict;
use autodie;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Sort::Versions;
use Data::Dump;

my $scriptname = basename($0);
my $version = "v0.2.0_092014";
my $description = <<"EOT";
Input one more more VCF files from IR output and generate a report of called CNVs. Can print anything
called a CNV, or filter based on gene name, copy number, number of tiles, or hotspot calls.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <VCF_file(s)>
    Filter Options
    -n, --novel       Print non-HS CNVs (Default = OFF)
    -c, --copies      Only print CNVs with at least this copy number 
    -g, --gene        Print out results for this gene only
    -t, --tiles       Only print out results for CNVs with at least this many tiles.

    -o, --output      Send output to custom file.  Default is STDOUT.
    -v, --version     Version information
    -h, --help        Print this help information
EOT

my $help;
my $ver_info;
my $outfile;
my $novel;
my $threshold;
my $geneid;
my $tiles;

GetOptions( "novel|n"       => \$novel,
            "copies|c=i" => \$threshold,
            "gene|g=s"      => \$geneid,
            "tiles|t=i"     => \$tiles,
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
    print "ERROR: No VCF files passed to script!\n\n"; 
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

#########------------------------------ END ARG Parsing ---------------------------------#########
my %cnv_data;
my @vcfs = @ARGV;
for my $input_file (@vcfs) {
    my $sample_name;
    open( my $vcf_fh, "<", $input_file );
    while (<$vcf_fh>) {
        next if /^##/;
        my @data = split;
        if ( $data[0] =~ /^#/ ) {
            $sample_name = $data[-1];
            next;
        }
        next unless $data[4] eq '<CNV>';
        if ( $data[4] eq '<CNV>' ) {
            my $varid = join( ':', @data[0..3] );
            
            # Kludgy, but need to deal with HS field; not like others!
            $data[7] =~ s/HS/HS=Yes/;
            my @format = split( /;/, $data[7] );
            my ($cn) = $data[9] =~ /:([^:]+)$/;
            push( @format, "CN=$cn" );

            %{$cnv_data{$sample_name}->{$varid}} = map { split /=/ } @format;
        }
    }
}

#dd \%cnv_data;
#exit;

# Set up for printing output
my @outfields = qw( END LEN NUMTILES RAW_CN REF_CN CN HS );
my @header = qw( Chr Gene Start End Length Tiles Raw_CN Ref_CN CN );
my $format = "%-8s %-8s %-11s %-11s %-11s %-8s %-8s %-8s %-8s\n";

select $out_fh;

# Print out each sample's data
for my $sample ( keys %cnv_data ) {
    my $count;
    print "::: CNVs Found in $sample :::\n";
    printf $format, @header;

    for my $cnv ( sort { versioncmp ( $a, $b ) } keys %{$cnv_data{$sample}} ) {
        my ($ci_5, $ci_95) = $cnv_data{$sample}->{$cnv}->{'CI'} =~ /0\.05:(\d\.\d+),0\.95:(\d\.\d+)/;
        my ($chr, $start, $gene, undef) = split( /:/, $cnv );
        my ($end, $length, $numtiles, $raw_cn, $ref_cn, $cn, $hs) = map { $cnv_data{$sample}->{$cnv}->{$_} } @outfields;

        # Filter data
        next if ( ! $novel && ( $gene eq '.' || ! defined $hs ) );
        if ( $threshold ) {
            next unless ( $cn >= $threshold );
        }
        if ( $tiles ) {
            next unless ( $numtiles >= $tiles );
        }
        if ( $geneid ) {
            next unless ( $gene eq $geneid );
        }
        $count++;

        printf $format, $chr, $gene, $start, $end, $length, $numtiles, $raw_cn, $ref_cn, $cn;
    }
    unless ( defined $count ) {
        print ">>> No CNVs found with the applied filters! <<<\n";
    }
    print "\n";
}