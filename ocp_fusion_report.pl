#!/usr/bin/perl
# Print out Novel fusions as found in the OCP Fusion Pipeline.  Will expand later, but for now just 
# print out the novel ones.  Can add some filters and stuff later if wanted.
# 6/9/2014 - D Sims
#######################################################################################################
use warnings;
use strict;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Data::Dump;

my $scriptname = basename($0);
my $version = "v0.6.0_082614";
my $description = <<"EOT";
Print out a summary table of fusions detected by the OCP Fusion Workflow VCF files. Can choose to output
anything seen, or just limit to annotated fusions.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <vcf_file(s)>
    -r, --ref       Include reference variants too (DEFAULT = off)
    -n, --novel     Include 'Novel' fusions in the output (DEFAULT = off)
    -o, --output    Write output to file <default =  STDOUT>
    -v, --version   Display version information
    -h, --help      Display this help text
EOT

my $help;
my $ver_info;
my $novel;
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
        my ( $count, $gene, $annot );
        if ( grep { /Fusion/ } @data ) {
            my $fusion = $data[2];
            ($count, $gene) = map { /READ_COUNT=(\d+);GENE_NAME=(.*?);/ } @data;
            # TODO: Not sure this is accurate.  how are novel fusions reported.  Are those without IDs novel?
            #( grep { /NOVEL/ } @data ) ? ($annot = "NOVEL") : ( ($annot) = map { /ANNOTATION=(.*?);FUNC/ } @data );
            ( ! grep { /ANNOTATION/ } @data ) ? ($annot = "NOVEL") : ( ($annot) = map { /ANNOTATION=(.*?);FUNC/ } @data );
            $fwidth = length($fusion) if ( length($fusion) > $fwidth );
            $results{$name}->{$fusion} = [$gene,$count,$annot];
        }
    }
}
#dd \%results;
#exit;

$fwidth += 4; # Give a little extra padding

select $out_fh;
for my $sample ( sort keys %results ) {
    my $count;
    printf "%-20s\n", $sample;
    printf "\t%-${fwidth}s%-15s %-8s %-8s\n", "Fusion", "Gene", "Count", "Annot";
    for my $fusion ( sort keys %{$results{$sample}} ) {
        (my $fname = $fusion) =~ s/(.*?\.(:?\S+))\..*_\d+$/$1/;
        if ( ! $novel ) {
            next if ( ${$results{$sample}->{$fusion}}[2] eq 'NOVEL' );
        }
        next if ( ${$results{$sample}->{$fusion}}[1] == 0 && ! $ref_calls );
        printf "\t%-${fwidth}s", $fname;
        printf "%-15s %-8s %-8s\n", @{$results{$sample}->{$fusion}};
        $count++;
    }
    print "\t\t\t<<< No Fusions Detected >>>\n" unless $count;
    print "\n";
}
