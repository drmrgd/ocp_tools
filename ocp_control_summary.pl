#!/usr/bin/perl
# Get expression control data from IR Fusion output
use warnings;
use strict;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use List::Util qw(max);
use Data::Dump;

my $scriptname = basename($0);
my $version = "v1.1.0_080515";
my $description = <<"EOT";
Program to pull out control data from VCF files generated from the OCP fusion pipeline on IR. Can
report both the internal expression control data and the 5'3'Assay data.  
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <vcf_file(s)>
    -F, --FP        Output the 5'3' Assay data.
    -f, --focus     Data is from the Oncomine Focus Panel (LMNA replaced with LRP1)
    -o, --output    Write output to file <default is STDOUT>
    -v, --version   Display version information
    -h, --help      Display this help text
EOT

my $help;
my $ver_info;
my $five_to_three;
my $outfile;
my $focus_panel;

GetOptions( 
    "FP|F"         => \$five_to_three,
    "focus|f"      => \$focus_panel,
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

my @files = @ARGV;
my %results;
my @expr_controls;

# Need to specify which depending on the panel run.
($focus_panel) ? 
    (@expr_controls = qw( LRP1 TBP MYC HMBS ITGB7 )) : 
    (@expr_controls = qw( LMNA TBP MYC HMBS ITGB7 ));

my @fp_controls = qw( NTRK1_5p3p ALK_5p3p ROS1_5p3p RET_5p3p );

my $out_fh;
if ( $outfile ) {
    open( $out_fh, ">", $outfile ) || die "Can't open  the file 'outfile' for writing: $!";
} else {
    $out_fh = \*STDOUT;
}

####### ------------------------------------- END ARG PARSING ---------------------------- #######

my @controls;
for my $input_file ( @files ) {
    (my $name = basename($input_file)) =~ s/\.vcf//;
    open( my $in_fh, "<", $input_file ) || die "Can't open the file '$input_file' for reading: $!";
    my %parsed_data;
    while (<$in_fh>) {
        next if /^#/;
        if ( /ExprControl/ ) {
            my @data = split;
            my ($gene, $count) = map { /GENE_NAME=(.*?);READ_COUNT=(\d+);.*/ } @data;
            $parsed_data{expr}->{$gene} = $count;
        }
        # Add in the 5P3P Data
        if ( /5p3pAssays/ ) {
            my @data = split;
            my ($gene, $count, $ratio) = map { /GENE_NAME=(.*?);READ_COUNT=(\d+,\d+);5P_3P_ASSAYS=(.*?);/ } @data;
            $parsed_data{fptp}->{"${gene}_5p3p"} = [$count,sprintf("%.4g", $ratio)];
        }
    }

    # If no data collected, this might be a DNA only sample and not run through fusion pipeline
    if ( ! %parsed_data ) {
        print "ERROR: No control data found in the VCF file '$input_file'. Is this data from an RNA sample run through the fusion pipeline?\n";
        exit 1;
    }

    # Convert zeros in output to '---' to make the table a little cleaner
    for my $control ( @expr_controls ) {
        ($parsed_data{expr}->{$control} == 0) ? 
        ($results{$name}->{expr}{$control} = ' ---') : 
        ($results{$name}->{expr}{$control} = $parsed_data{expr}->{$control} );
    }
    
    for my $control ( @fp_controls ) {
        ($parsed_data{fptp}->{$control}[0] eq '0,0') ? 
        ($results{$name}->{fptp}{$control} = [' ---', ' ---']) : 
        ($results{$name}->{fptp}{$control} = $parsed_data{fptp}->{$control} );
    }
}

# Get the longest sample name width
my ($width) = max( map { length($_)+4 } keys %results );
my $top_pad = ($width+54);

# Create header
select $out_fh;
my $epad = "%-10s" x scalar(@expr_controls);
my ($fpad, $sub_header);
if ($five_to_three) {
    $fpad = "%-25s " x scalar(@fp_controls);
    $sub_header= (sprintf("%-12s  %-12s", 'Counts', 'Imb')) x 4;
}

printf "%${top_pad}s $fpad\n", '', @fp_controls if $five_to_three;
printf "%-${width}s $epad", 'Samples', @expr_controls;
($five_to_three) ? print "$sub_header\n" : print "\n";

# Print out all control data;
for my $sample ( sort keys %results ) {
    printf "%-${width}s", $sample;
    for my $control ( @expr_controls ) {
        printf "%-10s", $results{$sample}->{expr}{$control};
    }
    if ( $five_to_three ) {
        for my $control ( @fp_controls ) {
            my ($count, $ratio) = @{$results{$sample}->{fptp}{$control}};
            $count =~ s/,/ : /;
            printf "%-13s %-12s", $count, $ratio;
        }
    }
    print "\n";
}
