#!/usr/bin/perl
# Print out Novel fusions as found in the OCP Fusion Pipeline.  Will expand later, but for now just 
# print out the novel ones.  Can add some filters and stuff later if wanted.
# 6/9/2014 - D Sims
#######################################################################################################
use warnings;
use strict;

use Getopt::Long;
use File::Basename;
use Data::Dump;

my $scriptname = basename($0);
my $version = "v0.4.1_061614";
my $description = <<"EOT";
Print out a summary table of fusions detected by the OCP Fusion Workflow VCF files. Can choose to output
anything seen, or just limit to annotated fusions.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <vcf_file(s)>
    -r, --ref       Include reference variants too (DEFAULT = off)
    -n, --novel     Include 'Novel' fusions in the output
    -o, --output    Write output to file <default is STDOUT>
    -v, --version   Display version information
    -h, --help      Display this help text
EOT

my $help;
my $ver_info;
my $novel;
my $outfile;
my $ref_calls;

GetOptions( 
    "ref|r"      => \$ref_calls,
    "novel"      => \$novel,
    "output=s"   => \$outfile,
    "help"       => \$help,
    "version"    => \$ver_info,
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
my @fusions;
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
            $results{$name}->{$fusion} = [$gene,$count,$annot];
        }
    }
}
#dd \%results;
#exit;

select $out_fh;
for my $sample ( sort keys %results ) {
    printf "%-20s\n", $sample;
    printf "\t%-34s %-8s %-8s %-8s\n", "Fusion", "Gene", "Count", "Annot";
    for my $fusion ( sort keys %{$results{$sample}} ) {
        (my $fname = $fusion) =~ s/(.*?\.(:?\S+))\..*_\d+$/$1/;
        if ( ! $novel ) {
            next if ( grep { /NOVEL/ } @{$results{$sample}->{$fusion}} );
        }
        printf "\t%-35s", $fname;
        printf "%-8s %-8s %-8s\n", @{$results{$sample}->{$fusion}};
    }
    print "\n\n";
}

sub width_format {
    my $data = shift;
    my $width = 0;

    for my $key ( %$data ) {
        for my $fusion ( %{$$data{$key}} ) {
            $width = length($fusion) if ( length($fusion) > $width );
        }
    }
    return ( $width + 4 );
}
