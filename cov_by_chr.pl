#!/usr/bin/perl
# Get coverage metrics by chromosome for multiple BAM files using samtools idxstats

use warnings;
use strict;
use autodie;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Sort::Versions;
use IPC::Cmd 'can_run';
use List::Util 'sum';
use Data::Dump;

use constant DEBUG => 0;

my $scriptname = basename($0);
my $version = "v0.7.0_082214";
my $description = <<"EOT";
Using samtools idxstats, print out a table comparing read depth per chromosome for a set of BAM files.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <bam_file(s)>
    -c, --csv       Generate a CSV output instead of plain table
    -o, --output    Send output to custom file.  Default is STDOUT.
    -v, --version   Version information
    -h, --help      Print this help information
EOT

my $help;
my $ver_info;
my $outfile;
my $csv_out;

GetOptions( "csv|c"         => \$csv_out,
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
    print "ERROR: Need to input at least 1 BAM file to process\n"; 
    print $usage;
    exit 1;
}

if ( ! can_run('samtools') ) {
    print "ERROR: It does not look like `samtools` is installed on this system.  You must install samtools to use this script.\n";
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
my @bam_files = @ARGV;
verify_bams( \@bam_files );

my %data;
my @samples;

for my $bam_file ( @bam_files ) {
    # split it all out...useful at some point later?
    my ($sample, $type, $site, $barcode) = $bam_file =~ /(.*?)[-_]([DR]NA)_(\w+)_(IonXpress_\d+)/;
    $site = uc($site);

    if ( DEBUG ) {
        print "############  DEBUG  ############\n";
        print "\tsample   => $sample\n";
        print "\ttype     => $type\n";
        print "\tsite     => $site\n";
        print "\tbarcode  => $barcode\n";
        print "#################################\n";
    }

    my $sample_string = "${sample}_${type}";
    push( @samples, $sample_string );
    get_stats( \$bam_file, \$sample_string, \$site, \%data );
}


my $col_format= "  %-13s" x scalar(@samples);

for my $site ( keys %data ) {
    my %sample_totals;
    print  {$out_fh} "::: Reads by Chromosome for $site :::\n";
    if ( $csv_out ) {
        print {$out_fh} join( ',', @samples ), "\n";
    } else {
        printf {$out_fh} "      $col_format\n", @samples;
    }

    for my $chr ( sort { versioncmp( $a, $b ) } keys %{$data{$site}} ) { 
        my @chr_counts;;
        #print {$out_fh} "$chr\t";
        ($csv_out) ? print {$out_fh} "$chr," : print {$out_fh} "$chr\t";
        for my $elem ( @{$data{$site}->{$chr}} ) {
            for my $sample ( keys %$elem ) {
                push( @chr_counts, $$elem{$sample} );
                push( @{$sample_totals{$sample}}, $$elem{$sample} );
            }
        }
        if ( $csv_out ) {
            print {$out_fh} join( ',', @chr_counts ), "\n";
        } else {
            printf {$out_fh} "$col_format\n", @chr_counts; 
        }
    }
    #print {$out_fh} "Total:    ";
    ($csv_out) ? print {$out_fh} "Total:" : print {$out_fh} "Total:    ";
    for my $sample (@samples) {
        my $total = sum(@{$sample_totals{$sample}});
        if ( $csv_out ) {
            print {$out_fh} ",$total";
        } else {
            printf {$out_fh} "%-13s  ", $total;
        }
    }
    print "\n";
}
        

sub get_stats {
    my ( $file, $sample, $site, $output_data ) = @_;
    my %parsed_data;
    my $total;

    print "proc $$sample...\n\n" if DEBUG;

    open( my $idxstats, "-|", "samtools idxstats $$file" ) || die "ERROR: Can't open the stream: $!";

    while (<$idxstats>) {
        next unless ( /chr[\d+YX]/ ); 
        my ($chr, $length, $counts, undef) = split;
        my $holder = { $$sample => $counts };
        push( @{$$output_data{$$site}->{$chr}}, { $$sample => $counts } );
        $total += $counts;
    }
}

sub verify_bams {
    my $bfiles = shift;
    my $counter = 0;

    for my $file ( @$bfiles ) {
        if ( ! -e $file ) {
            print "ERROR: File '$file' does not exist. Skipping...\n";
            next;
        }
        elsif ( ! -e "$file.bai" ) {
            print "ERROR: There is no bam index for '$file'.  You must run 'samtools index' on '$file' first.  Skipping...\n";
            next;
        }
        else {
            $counter++;
        }
    }
    if ( ! $counter ) {
        print "ERROR: There were no valid BAM files to be processed\n";
        exit 1;
    }
}
