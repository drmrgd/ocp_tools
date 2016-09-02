#!/usr/bin/perl
# Read in VCF file from IR and output called CNVs
# 9/20/2014 - D Sims
#################################################################################################
#use warnings;
#use strict;
use autodie;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Sort::Versions;
use JSON -support_by_pp;
use Parallel::ForkManager;
use Data::Dump;

use constant DEBUG => 1;

# Remove when in prod.
print "\n";
print colored("*" x 50, 'bold yellow on_black'), "\n";
print colored("      DEVELOPMENT VERSION OF OCP_CNV_REPORT\n", 'bold yellow on_black');
print colored("*" x 50, 'bold yellow on_black');
print "\n\n";

my $scriptname = basename($0);
my $version = "v2.5.0_090116-dev";
my $description = <<"EOT";
Input one more more VCF files from IR output and generate a report of called CNVs. Can print anything
called a CNV, or filter based on gene name, copy number, number of tiles, or hotspot calls.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <VCF_file(s)>
    Filter Options
    --cn  INT      Only report amplifications above this threshold (DEFAULT: CN >= 0)
    --cu  INT      Upper bound for amplifications based on 5% CI (DEFAULT: 5% CI >= 4)
    --cl  INT      Lower bound for deletions based on 95% CI (DEFAULT: 95% CI <= 1)
    -n, --novel    Print non-HS CNVs (Default = OFF)
    -g, --gene     Print out results for this gene only. Can also input a list of comma separated gene names to search 
    -t, --tiles    Only print out results for CNVs with at least this many tiles.
    -a, --annot    Only print CNVs with Oncomine Annotations.
    -N, --NOCALL   Do not output NOCALL results (Default: OFF)

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
my $copy_number = 0;
my $cu;
my $cl;
my $geneid;
my $tiles;
my $annot;
my $nocall;
my $raw_output;

GetOptions( "novel|n"             => \$novel,
            "copy-number|cn=f"    => \$copy_number,
            "upper-copies|cu=f"   => \$cu,
            "lower-copies|cl=f"   => \$cl,
            "gene|g=s"            => \$geneid,
            "tiles|t=i"           => \$tiles,
            "annot|a"             => \$annot,
            "output|o=s"          => \$outfile,
            "NOCALL|N"            => \$nocall,
            "version|v"           => \$ver_info,
            "raw|r"               => \$raw_output,
            "help|h"              => \$help )
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

# We don't need both cn and cu or cl
if ($cl or $cu) {
    undef $copy_number;
}
my @genelist = split(/,/, $geneid) if $geneid;

my %filters = (
    'cn'    => $copy_number,
    'cu'    => $cu,
    'cl'    => $cl,
    'gene'  => [@genelist],
    'tiles' => $tiles,
    'annot' => $annot,
    'novel' => $novel,
);

if (DEBUG) {
    print '='x35, '  DEBUG  ', '='x35, "\n";
    print "Filters being employed\n";
    while (my ($keys, $values) = each %filters) {
        $values //= 'undef';
        printf "\t%-7s => %s\n",$keys,$values;
    }
    print '='x79, "\n";
}

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

my $pm = new Parallel::ForkManager(48);
$pm->run_on_finish(
    sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
        my $vcf = $data_structure_reference->{input};
        my $name = $data_structure_reference->{id};
        $name //= basename($vcf);
        $cnv_data{$$name} = $data_structure_reference->{result};
    }
);

for my $input_file (@vcfs) {
    $pm->start and next;
    my ($return_data, $sample_id)  = proc_vcf(\$input_file);

    $pm->finish(0, 
        { 
          result  => $return_data, 
          input   => $input_file, 
          id      => $sample_id,
        }
     );
}
$pm->wait_all_children;

# Set up for printing output
my @outfields = qw( END LEN NUMTILES RAW_CN REF_CN CN HS FUNC );
my @header = qw( Chr Gene Start End Length Tiles Raw_CN Ref_CN CI_05 CI_95 CN Annot );
my $format = "%-8s %-8s %-11s %-11s %-11s %-8s %-8s %-8s %-8s %-8s %-8s %-18s\n";

select $out_fh;

# Print out each sample's data
for my $sample ( keys %cnv_data ) {
    my ($id, $gender, $mapd, $cell) = split( /:/, $sample );
    my $count;
    print "::: CNV Data For $id (Gender: $gender, Cellularity: $cell, MAPD: $mapd) :::\n";
    printf $format, @header;

    for my $cnv ( sort { versioncmp ( $a, $b ) } keys %{$cnv_data{$sample}} ) {
        my %mapped_cnv_data;
        last if $cnv eq 'NONE';
        # Seems to be a bug in the same the CI are reported for deletions.  Solution correctly reports the value
        # in the VCF, but it's not so informative.  This will give a better set of data.
        my ($ci_5, $ci_95) = $cnv_data{$sample}->{$cnv}->{'CI'} =~ /0\.05:(.*?),0\.95:(.*)$/; 
        my ($chr, $start, $gene, undef) = split( /:/, $cnv );
        #my ($end, $length, $numtiles, $raw_cn, $ref_cn, $cn, $hs, $func) = map { $cnv_data{$sample}->{$cnv}->{$_} } @outfields;

        %mapped_cnv_data = map{ $_ => $cnv_data{$sample}->{$cnv}->{$_} } @outfields;
        @mapped_cnv_data{qw(ci_05 ci_95)} = ($ci_5,$ci_95);
        @mapped_cnv_data{qw(chr start gene undef)} = split(/:/, $cnv);

        #$hs //= 'No';
        #$ci_5 //= 0;
        #$ci_95 //= 0;
        $mapped_cnv_data{HS} //= 'No';
        $mapped_cnv_data{ci_05} //= 0;
        $mapped_cnv_data{ci_95} //= 0;

        # Get OVAT Annot Data
        my ($gene_class, $variant_class);
        my $func = $mapped_cnv_data{FUNC};
        #if ( $mapped_cnv_data{$func} && $mapped_cnv_data{$func} =~ /oncomine/ ) {
        if ( $func && $func =~ /oncomine/ ) {
            my $json_annot = JSON->new->allow_singlequote->decode($func);
            my $parsed_annot = $$json_annot[0];
            $gene_class = $$parsed_annot{'oncomineGeneClass'};
            $variant_class = $$parsed_annot{'oncomineVariantClass'};
        } else {
            $gene_class = $variant_class = '---';
        }
        $mapped_cnv_data{GC} = $gene_class;
        $mapped_cnv_data{VC} = $variant_class;

        my @filtered_data = filter_results(\%mapped_cnv_data, \%filters);
        if (@filtered_data) {
            dd \@filtered_data;
            print '-'x150, "\n";
        }

        #printf $format, $chr, $gene, $start, $end, $length, $numtiles, $raw_cn, $ref_cn, $ci_5, $ci_95, $cn, $gene_class;
        #$count++;
    #}
    #unless ( defined $count ) {
        #print "\t\t>>> No CNVs found with the applied filters! <<<\n";
    #}
    #print "\n";
    }
}

sub return_data {
    my $data = shift;
    my @fields = qw( chr gene start END LEN NUMTILES RAW_CN REF_CN ci_05 ci_95 CN GC );
    return @$data{@fields};
}

sub copy_number_filter {
    # if there is a 5% / 95% CI filter, use that, otherwise if there's a cn filter use that, and 
    # if there's nothing, just return it all.
    my ($data, $threshold) = @_;
    my ($cn, $cu, $cl) = @$threshold;

    if ($cu and $cl) {
        ($$data{ci_05} >= $cu || $$data{ci_95} <= $cl) ? return 1 : return 0;
    } 
    elsif ($cn) {
        ($$data{cn} >= $cn) ? return 1 : return 0;
    } else {
        return 1;
    }
}

sub filter_results {
    my ($data, $filters) = @_;
    my @cn_thresholds = $$filters{qw(cn cu cl)};

    # Filter non-hotspots and novel if we don't want them.
    #unless ($$filters{novel} && ($$data{gene} eq '.' or $$data{HS})) return;


    # Gene level filter
    if ($$filters{gene}) {
        if ( grep { $$data{gene} eq uc($_) } @{$$filters{gene}} ) {
            return return_data($data) if copy_number_filter($data, \@cn_thresholds);
        } 
    }

    # OVAT Filter
    if ($$filters{annot} and $$data{GC} ne '---') {
        return return_data($data) if copy_number_filter($data, \@cn_thresholds);
    }

    # Number of tiles filter
    if ($$filters{tiles} and $$data{NUMTILES} > $$filters{tiles}) {
        return return_data($data) if copy_number_filter($data, \@cn_thresholds);
    }

    return;
    return return_data($data);
}

sub proc_vcf {
    my $vcf = shift;
    my ($sample_id, $gender, $mapd, $cellularity, $sample_name);
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

        #my $sample_id = join( ':', $sample_name, $gender, $mapd, $cellularity );
        $sample_id = join( ':', $sample_name, $gender, $mapd, $cellularity );

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

        %{$results{$varid}} = map { split /=/ } @format;
    }
    if (DEBUG) {
        print "="x35, "  DEBUG  ", "="x35, "\n";
        print "\tSample Name:  $sample_name\n";
        print "\tCellularity:  $cellularity\n";
        print "\tGender:       $gender\n";
        print "\tMAPD:         $mapd\n";
        print "="x79, "\n";
    }
    #return \%results, \$sample_name, \$cellularity, \$gender, \$mapd;
    return \%results, \$sample_id;
}
