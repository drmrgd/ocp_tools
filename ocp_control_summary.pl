#!/usr/bin/perl
# Get expression control data from IR Fusion output
use warnings;
use strict;

use Getopt::Long;
use File::Basename;
use Data::Dump;

my $scriptname = basename($0);
my $version = "v0.3.0_061014";
my $description = <<"EOT";
Program to pull out control data from VCF files generated from the OCP fusion pipeline.  For now, I'm 
running this on files that have already been filtered for reference calls.  This can be added if need 
be.  Also the sample naming is based on one file name type only, and will probably not work great for 
others.  
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <vcf_file(s)>
    -FP, --FP       Output the 5P3P assay data with the HKG control data
    -o, --output    Write output to file <default is STDOUT>
    -v, --version   Display version information
    -h, --help      Display this help text
EOT

my $help;
my $ver_info;
my $five_to_three;
my $outfile;

GetOptions( 
    "FP"         => \$five_to_three,
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

my @files = @ARGV;
my %results;
#my @controls = qw( LMNA TBP MYC HMBS ITGB7 NTRK1_5p3p ALK_5p3p ROS1_5p3p RET_5p3p );
my @controls = qw( LMNA TBP MYC HMBS ITGB7 NTRK1_5p3p ALK_5p3p ROS1_5p3p RET_5p3p );

if ( ! $five_to_three ) {
    splice( @controls, 5 );
}

#dd \@controls;
#exit;

my $out_fh;
if ( $outfile ) {
    open( $out_fh, ">", $outfile ) || die "Can't open  the file 'outfile' for writing: $!";
} else {
    $out_fh = \*STDOUT;
}

####### ------------------------------------- END ARG PARSING ---------------------------- #######

for my $input_file ( @files ) {
    # TODO Figure out a good naming scheme.  Might have to standardize file naming initially for now
    #      do this in a few steps to trim out everything we know we don't want for sure.
    (my $name = $input_file) =~ s/(:?_Fusion_filtered)?\.vcf$//i;
    $name =~ s/_RNA//;
    open( my $in_fh, "<", $input_file ) || die "Can't open the file '$input_file' for reading: $!";
    my %parsed_data;
    while (<$in_fh>) {
        next if ( /^#/ || /READ_COUNT=0/ );
        if ( /ExprControl/ ) {
            my @data = split;
            my ($gene, $count) = map { /GENE_NAME=(.*?);READ_COUNT=(\d+);.*/ } @data;
            $parsed_data{$gene} = $count;
        }
        # Add in the 5P3P Data
        if ( /5p3pAssays/ ) {
            my @data = split;
            my ($gene,$count, $ratio) = map { /GENE_NAME=(.*?);READ_COUNT=(\d+,\d+);5P_3P_ASSAYS=(.*?);/ } @data;
            $parsed_data{"${gene}_5p3p"} = $count;
        }
    }

    # Summarize and tie together results for output.
    for my $control ( @controls ) {
        ( exists $parsed_data{$control} ) ? ( $results{$name}->{$control} = $parsed_data{$control} ) : ($results{$name}->{$control} = "---");
    }
}
#dd \%results;
#exit;

# Create header
select $out_fh;
my $padding = "%-12s " x scalar(@controls);
printf "%-19s $padding\n", "Samples", @controls;

for my $sample ( sort keys %results ) {
    printf "%-20s", $sample;
    for my $control ( @controls ) {
        printf "%-13s", $results{$sample}->{$control};
    }
    print "\n";
}
