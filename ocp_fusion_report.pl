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
use Sort::Versions;

my $scriptname = basename($0);
my $version = "v1.4.0_031416";
my $description = <<"EOT";
Print out a summary table of fusions detected by the OCP Fusion Workflow VCF files. Can choose to output
anything seen, or just limit to annotated fusions.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <vcf_file(s)>
    -r, --ref       Include reference variants too (DEFAULT: off)
    -n, --novel     Include 'Novel' fusions in the output (DEFAULT: on)
    -g, --gene      Only output data for a specific driver gene.
    -o, --output    Write output to file <default =  STDOUT>
    -v, --version   Display version information
    -h, --help      Display this help text
EOT

my $help;
my $ver_info;
my $novel=1;
my $outfile;
my $ref_calls;
my $gene;

GetOptions( 
    "ref|r"        => \$ref_calls,
    "novel|n"      => \$novel,
    "gene|g=s"     => \$gene,
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
my @drivers = qw( ABL1 AKT3 ALK AXL BRAF CDK4 EGFR ERBB2 ERG ETV1 ETV1a ETV1b ETV4 ETV4a ETV5 ETV5a ETV5b ETV5d 
                  FGFR1 FGFR2 FGFR3 MET NTRK1 NTRK2 NTRK3 PDGFRA PPARG RAF1 RET ROS1);

for my $input_file ( @files ) {
    (my $sample_name = $input_file) =~ s/(:?_Fusion_filtered)?\.vcf$//i;
    $sample_name =~ s/_RNA//;

    open( my $in_fh, "<", $input_file ) || die "Can't open the file '$input_file' for reading: $!";
    while (<$in_fh>) {
        next if /^#/;
        my @data = split;
        if ( grep { /Fusion/ } @data ) {
            my ( $name, $elem ) = $data[2] =~ /(.*?)_([12])$/;
            #my ($count, $gene) = map { /READ_COUNT=(\d+);GENE_NAME=(.*?);/ } @data;
            my ($count) = map { /READ_COUNT=(\d+)/ } @data;
            my ($pair, $junct, $id) = split(/\./, $name);
            $id //= '-';

            my ($gene1, $gene2) = split(/-/, $pair);

            # Filter out ref calls if we don't want to view them
            if ( $count == 0 ) { next unless ( $ref_calls ) }
            my $fid = join('|', $pair, $junct, $id);

            #print "$data[2]  => \n";
            #print "\tname:  $name\n";
            #print "\tgene:  $gene\n";
            #print "\tgene1: $gene1\n";
            #print "\tgene2: $gene2\n";
            #print "\tpair:  $pair\n";
            #print "\tjunct: $junct\n";
            #print "\tID:   $id\n";
            #print '-'x50;
            #print "\n";

            if ( $pair eq 'MET-MET' || $pair eq 'EGFR-EGFR' ) {
                $results{$sample_name}->{$fid}->{'DRIVER'} = $results{$sample_name}->{$fid}->{'PARTNER'} = $gene1;
            }
            elsif (grep {$_ eq $gene1} @drivers) {
                $results{$sample_name}->{$fid}->{'DRIVER'} = $gene1;
                $results{$sample_name}->{$fid}->{'PARTNER'} = $gene2;
            }
            elsif (grep {$_ eq $gene2} @drivers) {
                $results{$sample_name}->{$fid}->{'DRIVER'} = $gene2;
                $results{$sample_name}->{$fid}->{'PARTNER'} = $gene1;
            }
            else {
                $results{$sample_name}->{$fid}->{'DRIVER'} = 'UNKNOWN';
                $results{$sample_name}->{$fid}->{'PARTNER'} = "$gene1,$gene2";
            }
            $results{$sample_name}->{$fid}->{'COUNT'} = $count;

            # Get some field width formatting.
            $fwidth = (length($name)+4) if ( length($name) > $fwidth );
        }
    }

    # At least generate a hash entry for a sample with no fusions for the final output table below
    if ( ! $results{$sample_name} ) {
        $results{$sample_name} = undef;
    }
}

#dd \%results;
#exit;

# Generate and print out the final results table(s)
select $out_fh;
for my $sample ( sort keys %results ) {
    print "::: ";
    print uc($gene) if $gene;
    print " Fusions in $sample :::\n\n";

    if ( $results{$sample} ) {
        my $fusion_format = "%-${fwidth}s %-12s %-12s %-15s %-15s\n";
        my @fusion_header = qw (Fusion ID Read_Count Driver_Gene Partner_Gene);
        printf $fusion_format, @fusion_header;
        for (sort { versioncmp($a,$b) } keys %{$results{$sample}} ) {
            if ($gene) {
                next unless $results{$sample}->{$_}->{'DRIVER'} eq uc($gene);
            }
            my ($fusion, $junct, $id ) = split(/\|/);
            printf $fusion_format, "$fusion.$junct", $id, $results{$sample}->{$_}->{'COUNT'}, $results{$sample}->{$_}->{'DRIVER'}, $results{$sample}->{$_}->{'PARTNER'};
        }
        print "\n";
    } else {
        print "\t\t\t<<< No Fusions Detected >>>\n\n";
    }
}
