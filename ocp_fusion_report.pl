#!/usr/bin/perl
# Generate a fusion report from VCF files derived from and IR analysis.  Can output either annotated
# fusions only or both annotated and novel fusions.  Can also choose to output ref calls in addtion
# to variant calls to make a more complete report.
#
# 6/9/2014 - D Sims
#######################################################################################################
use warnings;
use strict;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Data::Dump;

my $scriptname = basename($0);
my $version = "v0.8.1_010815";
my $description = <<"EOT";
Print out a summary table of fusions detected by the OCP Fusion Workflow VCF files. Can choose to output
anything seen, or just limit to annotated fusions.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <vcf_file(s)>
    -r, --ref       Include reference variants too (DEFAULT: off)
    -n, --novel     Include 'Novel' fusions in the output (DEFAULT: on)
    -o, --output    Write output to file <default =  STDOUT>
    -v, --version   Display version information
    -h, --help      Display this help text
EOT

my $help;
my $ver_info;
my $novel=1;
my $outfile;
my $ref_calls;

GetOptions( 
    "ref|r"        => \$ref_calls,
    "novel|n"      => \$novel,
    "output|o=s"   => \$outfile,
    "help|h"       => \$help,
    "version|v"    => \$ver_info,
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
    open( $out_fh, ">", $outfile ) || die "Can't open the output file '$outfile' for writing: $!";
} else {
    $out_fh = \*STDOUT;
}

my @files = @ARGV;

#######===========================  END ARG Parsing  #######=========================== 

my %results;
my $fwidth=0;
for my $input_file ( @files ) {
    (my $name = $input_file) =~ s/(:?_Fusion_filtered)?\.vcf$//i;
    $name =~ s/_RNA//;

    open( my $in_fh, "<", $input_file ) || die "Can't open the file '$input_file' for reading: $!";
    while (<$in_fh>) {
        next if /^#/;
        my @data = split;
        if ( grep { /Fusion/ } @data ) {
            my ( $fname, $felem ) = $data[2] =~ /(.*?)_([12])$/;
            my ($count, $gene) = map { /READ_COUNT=(\d+);GENE_NAME=(.*?);/ } @data;

            # Filter out ref calls if we don't want to view them
            if ( $count == 0 ) { next unless ( $ref_calls ) }

            my ($annot) = map { /ANNOTATION=(.*);FUNC/ } @data;
            $annot //= 'NOVEL';
            my $fid = join( '|', $fname, $annot );

            if ( $felem == 2 ) {
                $results{$name}->{$fid}->{'DRIVER'} = [$gene, $count];
            } else {
                $results{$name}->{$fid}->{'PARTNER'} = [$gene, $count];
            }

            # Get dynamic field width for fusion ID for the final output table
            $fwidth = (length($fname)+4) if ( length($fname) > $fwidth );
        }
    }
    # At least generate a hash entry for a sample with no fusions for the final output table below
    if ( ! $results{$name} ) {
        $results{$name} = undef;
    }
}

#dd \%results;
#exit;

# Generate and print out the final results table(s)
select $out_fh;
for my $sample ( sort keys %results ) {
    print "::: Fusions in $sample :::\n\n";
    if ( $results{$sample} ) {
        printf "%-${fwidth}s %12s %24s\n", ' ', 'Driver', 'Partner'; 
        printf "%-${fwidth}s%-15s %-8s %-15s %-8s %-8s\n", qw( Fusion Gene Count Gene Count Annot );
        for my $fusion ( sort keys %{$results{$sample}} ) {
            my ($fus_id, $annot) = split( /\|/, $fusion );
            next if ( $annot eq 'NOVEL' && ! $novel );
            printf "%-${fwidth}s", $fus_id;
            printf "%-15s %-8s %-15s %-8s %-8s\n", @{$results{$sample}->{$fusion}->{'DRIVER'}}, @{$results{$sample}->{$fusion}->{'PARTNER'}}, $annot;
        }
        print "\n";
    } else {
        print "\t\t\t<<< No Fusions Detected >>>\n\n";
    }
}
