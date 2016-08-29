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
use JSON -support_by_pp;
use Parallel::ForkManager;
use Data::Dump;

use constant DEBUG => 0;

my $scriptname = basename($0);
my $version = "v2.2.0_082916";
my $description = <<"EOT";
Input one more more VCF files from IR output and generate a report of called CNVs. Can print anything
called a CNV, or filter based on gene name, copy number, number of tiles, or hotspot calls.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <VCF_file(s)>
    Filter Options
    -n, --novel       Print non-HS CNVs (Default = OFF)
    -c, --copies      Only print CNVs with at least this copy number 
    -g, --gene        Print out results for this gene only. Can also input a list of comma separated gene names to search 
    -t, --tiles       Only print out results for CNVs with at least this many tiles.
    -a, --annot       Only print CNVs with Oncomine Annotations.
    -N, --NOCALL      Do not output NOCALL results (Default: OFF)

    Output Options
    -o, --output      Send output to custom file.  Default is STDOUT.
    -r, --raw         Output raw CSV data, rather than pretty printed report.
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
my $annot;
my $nocall;
my $raw_output;

GetOptions( "novel|n"       => \$novel,
            "copies|c=i"    => \$threshold,
            "gene|g=s"      => \$geneid,
            "tiles|t=i"     => \$tiles,
            "annot|a"       => \$annot,
            "output|o=s"    => \$outfile,
            "NOCALL|N"      => \$nocall,
            "version|v"     => \$ver_info,
            "raw|r"         => \$raw_output,
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
my $cnv_data;
my @vcfs = @ARGV;

my $pm = new Parallel::ForkManager(48);
$pm->run_on_finish(
    sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
        my $vcf = $data_structure_reference->{input};
        my $name = $data_structure_reference->{sample};
        $name //= basename($vcf);
        #$cnv_data{$name} = $data_structure_reference->{result};
        $cnv_data = $data_structure_reference->{result};
    }
);

for my $input_file (@vcfs) {
    $pm->start and next;
    my ($return_data, $sample_name, $cellularity, $gender, $mapd)  = proc_vcf(\$input_file);
    $pm->finish(0, 
        { result      => $return_data, 
          sample      => $$sample_name, 
          cellularity => $cellularity, 
          gender      => $gender, 
          mapd        => $mapd,
          input       => $input_file, 
        }
     );
}
$pm->wait_all_children;

#dd \%cnv_data;
dd $cnv_data;
exit;


# Set up for printing output
my @outfields = qw( END LEN NUMTILES RAW_CN REF_CN CN HS FUNC );
my @header = qw( Chr Gene Start End Length Tiles Raw_CN Ref_CN CI_05 CI_95 CN Annot );
my $format = "%-8s %-8s %-11s %-11s %-11s %-8s %-8s %-8s %-8s %-8s %-8s %-18s\n";

select $out_fh;

# Print out each sample's data
for my $sample ( keys %cnv_data ) {
    my ($id, $gender, $mapd, $cell) = split( /:/, $sample );
    my $count;
    print "::: CNVs Data For $id (Gender: $gender, Cellularity: $cell, MAPD: $mapd) :::\n";
    printf $format, @header;

    for my $cnv ( sort { versioncmp ( $a, $b ) } keys %{$cnv_data{$sample}} ) {
        last if $cnv eq 'NONE';
        # Seems to be a bug in the same the CI are reported for deletions.  Solution correctly reports the value
        # in the VCF, but it's not so informative.  This will give a better set of data.
        my ($ci_5, $ci_95) = $cnv_data{$sample}->{$cnv}->{'CI'} =~ /0\.05:(.*?),0\.95:(.*)$/; 
        my ($chr, $start, $gene, undef) = split( /:/, $cnv );
        my ($end, $length, $numtiles, $raw_cn, $ref_cn, $cn, $hs, $func) = map { $cnv_data{$sample}->{$cnv}->{$_} } @outfields;
        $hs //= 'No';

        # Filter data
        next if ( ! $novel && ($gene eq '.' || $hs eq 'No') );
        if ( $threshold ) {
            next unless ( $cn >= $threshold );
        }
        if ( $tiles ) {
            next unless ( $numtiles >= $tiles );
        }
        if ( $geneid ) {
            my @genelist = split(/,/, $geneid);
            # Allow for case insensitive searching...I'm too lazy for the shift key!
            next unless ( grep { $gene eq uc($_) } @genelist );
        }

        # Need to add this to fix null 5% and 95% CI values if they don't exist
        $ci_5 //= 0;
        $ci_95 //= 0;

        # Get OVAT Annot Data
        my ($gene_class, $variant_class);
        if ( $func && $func =~ /oncomine/ ) {
            my $json_annot = JSON->new->allow_singlequote->decode($func);
            my $parsed_annot = $$json_annot[0];
            $gene_class = $$parsed_annot{'oncomineGeneClass'};
            $variant_class = $$parsed_annot{'oncomineVariantClass'};
        } else {
            $gene_class = $variant_class = '---';
        }

        # Filter out non-oncomine CNVs
        if ( $annot ) {
            next unless $gene_class ne '---';
        }
        printf $format, $chr, $gene, $start, $end, $length, $numtiles, $raw_cn, $ref_cn, $ci_5, $ci_95, $cn, $gene_class;
        $count++;
    }
    unless ( defined $count ) {
        print "\t\t>>> No CNVs found with the applied filters! <<<\n";
    }
    print "\n";
}

sub proc_vcf {
    my $vcf = shift;
    my ($gender, $mapd, $cellularity, $sample_name);
    my %results;

    open( my $vcf_fh, "<", $$vcf);
    while (<$vcf_fh>) {
        if ( /^##/ ) {
            if ( $_ =~ /sampleGender=(\w+)/ ) {
                $gender = $1 and next;
            }
            # Need to add to accomodate the new CNV plugin; may not have the same field as the normal IR data.
            if ($_ =~ /AssumedGender=([mf])/) {
                ($1 eq 'm') ? ($gender='Male') : ($gender='Female');
                next;
            }
            elsif ( $_ =~ /mapd=(\d\.\d+)/ ) {
                $mapd = $1 and next;
            }
            elsif ( $_ =~ /CellularityAsAFractionBetween0-1=(.*)$/ ) {
                $cellularity = $1 and next;
            }
        } 

        my @data = split;
        if ( $data[0] =~ /^#/ ) {
            $sample_name = $data[-1] and next;
        }
        next unless $data[4] eq '<CNV>';

        my $sample_id = join( ':', $sample_name, $gender, $mapd, $cellularity );

        # Let's handle NOCALLs for MATCHBox compatibility (prefer to filter on my own though).
        if ($nocall && $data[6] eq 'NOCALL') {
            ${$cnv_data{$sample_id}->{'NONE'}} = '';
            next;
        }

        my $varid = join( ':', @data[0..3] );
        
        # Kludgy, but need to deal with hotspots (HS) field; not like others!
        $data[7] =~ s/HS/HS=Yes/;
        $data[7] =~ s/SD;/SD=NA;/; # sometimes data in this field and sometimes not.  

        my @format = split( /;/, $data[7] );
        my ($cn) = $data[9] =~ /:([^:]+)$/;
        push( @format, "CN=$cn" );

        %{$results{$sample_id}->{$varid}} = map { split /=/ } @format;
    }
    if (DEBUG) {
        print "="x40, "  DEBUG  ", "="x40, "\n";
        print "\tSample Name:  $sample_name\n";
        print "\tCellularity:  $cellularity\n";
        print "\tGender:       $gender\n";
        print "\tMAPD:         $mapd\n";
        print "="x89, "\n";
    }
    return \%results, \$sample_name, \$cellularity, \$gender, \$mapd;
}
