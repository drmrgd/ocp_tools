#!/usr/bin/perl
# Get expression control data from IR Fusion output
use warnings;
use strict;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use List::Util qw(max sum);
use Data::Dump;

my $scriptname = basename($0);
my $version = "v2.1.0_090116";
my $description = <<"EOT";
Program to pull out control data from VCF files generated from the OCP fusion pipeline on IR. Can
report both the internal expression control data and the 5'3'Assay data.  
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <vcf_file(s)>
    -F, --FP                  Output the 5'3' Assay data.
    -m, --match_version  INT  Version of the NCI-MATCH panel to define which controls used.
    -o, --output              Write output to file <default is STDOUT>
    -v, --version             Display version information
    -h, --help                Display this help text
EOT

my $help;
my $ver_info;
my $five_to_three;
my $outfile;
my $match_version;

GetOptions( 
    "FP|F"               => \$five_to_three,
    "match_version|m=i"  => \$match_version,
    "output|o=s"         => \$outfile,
    "help|h"             => \$help,
    "version|v"          => \$ver_info,
) or die $usage;

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

# Make sure enough args passed to script
if (@ARGV < 1) {
    print "ERROR: No VCF files passed to script!\n\n"; 
    print "$usage\n";
    exit 1;
}
my @files = @ARGV;

my $out_fh;
if ( $outfile ) {
    open( $out_fh, ">", $outfile ) || die "Can't open  the file 'outfile' for writing: $!";
} else {
    $out_fh = \*STDOUT;
}

####### ------------------------------------- END ARG PARSING ---------------------------- #######
my (@expr_controls, @fp_controls, @used_controls);
my %results;
for my $input_file ( @files ) {
    # Need to specify which depending on the panel run.
    @expr_controls = select_expr_controls(\$input_file, $match_version);
    for my $control (@expr_controls) {
        push(@used_controls, $control) unless grep {$_ eq $control} @used_controls;
    }
    @fp_controls = qw( NTRK1_5p3p ALK_5p3p ROS1_5p3p RET_5p3p );
    (my $name = basename($input_file)) =~ s/\.vcf//;
    open( my $in_fh, "<", $input_file ) || die "Can't open the file '$input_file' for reading: $!";
    my %parsed_data;
    my $sum = 0;
    my $sample;
    while (<$in_fh>) {
        if (/^#CHROM/) {
            my @elems = split(/\s+/, $_);
            $sample = $elems[-1];
            $sample //= $name;
            $results{$name}->{sample} = $sample;
        }

        next if /^#/;
        if ( /SVTYPE=ExprControl/ ) {
            my @data = split;
            my ($gene, $count) = map { /GENE_NAME=(.*?);READ_COUNT=(\d+);.*/ } @data;
            $parsed_data{expr}->{$gene} = $count;
            $sum += $count;
        }
        
        # Add in the expression control sum data
        $parsed_data{expr}->{'Total'} = $sum;

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
        ($parsed_data{expr}->{$control} == 0)
        ?  ($results{$name}->{expr}{$control} = ' ---')
        :  ($results{$name}->{expr}{$control} = $parsed_data{expr}->{$control} );
    }
    
    for my $control ( @fp_controls ) {
        ($parsed_data{fptp}->{$control}[0] eq '0,0')
        ?  ($results{$name}->{fptp}{$control} = [' ---', ' ---'])
        :  ($results{$name}->{fptp}{$control} = $parsed_data{fptp}->{$control} );
    }
}

# Get the longest sample name width
#my ($width) = max( map { length($_)+4 } keys %results );
my $width = 0;
for my $file (keys %results) {
    my $len = length($results{$file}->{sample});
    $width = $len if $len > $width;
}
$width += 4;
my $top_pad = ($width+62);

# Create header
select $out_fh;
my $epad = "%-10s" x scalar(@used_controls);
my ($fpad, $sub_header);
if ($five_to_three) {
    $fpad = "%-25s " x scalar(@fp_controls);
    $sub_header= (sprintf("%-12s  %-12s", 'Counts', 'Imb')) x 4;
}

printf "%${top_pad}s $fpad\n", '', @fp_controls if $five_to_three;
printf "%-${width}s$epad", 'Samples', sort @used_controls;
($five_to_three) ? print "$sub_header\n" : print "\n";

# Print out all control data;
for my $sample ( sort keys %results ) {
    #printf "%-${width}s", $sample;
    printf "%-${width}s", $results{$sample}->{sample};
    for my $control (sort @used_controls) {
        $results{$sample}->{expr}{$control} //= 'N/A';
        printf "%-10s", $results{$sample}->{expr}{$control}
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

sub select_expr_controls {
    my ($vcf,$version) = @_;
    my %expr_controls = (
        1  => [qw( LMNA TBP MYC HMBS ITGB7 Total )],
        2  => [qw( LRP1 TBP MYC HMBS ITGB7 Total )],
    );
    return @{$expr_controls{$version}} if $version;

    open(my $fh, "<", $$vcf);
    my ($ovat_version) = map { /^##OncomineVariantAnnotationToolVersion=(\d\.\d)\.\d$/ } <$fh>;
    close $fh;
    ($ovat_version eq '2.0') ? return @{$expr_controls{1}} : return @{$expr_controls{2}};
}
