#!/usr/bin/env perl
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
use Term::ANSIColor;

use constant DEBUG => 0;

my $scriptname = basename($0);
my $version = "v3.4.053119-dev";

# Remove when in prod.
#print "\n";
#print colored("*" x 75, 'bold yellow on_black'), "\n";
#print colored("      DEVELOPMENT VERSION OF $scriptname\n", 'bold yellow on_black');
#print colored("*" x 75, 'bold yellow on_black');
#print "\n\n";

my $description = <<"EOT";
Input one more more VCF files from IR output and generate a report of called 
CNVs. Can print anything called a CNV, or filter based on gene name, copy 
number, number of tiles, or hotspot calls. Also can output a raw formatted table
that is easily imported into Excel and R for downstream analysis.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <VCF_file(s)>
    Filter Options
    --cn  INT        Only report amplifications above this threshold (DEFAULT: 
                     CN >= 0)
    --cu  INT        Upper bound for amplifications based on 5% CI (DEFAULT: 
                     5% CI >= 4)
    --cl  INT        Lower bound for deletions based on 95% CI (DEFAULT: 
                     95% CI <= 1)
    -n, --novel      Print non-HS CNVs (Default = OFF)
    -g, --gene       Print out results for this gene only. Can also input a list
                     of comma separated gene names to search 
    -t, --tiles      Only print out results for CNVs with at least this many 
                     tiles.
    -a, --annot      Only print CNVs with Oncomine Annotations.
    -p, --procs      Number of processor cores to use to process files in 
                     parallel. (Default: 48).
    -N, --NOCALL     Do not output NOCALL results (Default: OFF)

    Output Options
    -o, --output      Send output to custom file.  Default is STDOUT.
    -f, --format      Format to use for output.  Can choose 'csv', 'tsv', or 
                      pretty print (as 'pp').  DEFAULT: pp. 
                      This format still retains the header and other data.
    -r, --raw         Output as a raw CSV format that can be input into Excel.
                      The header output is not in columns as 
                      a part of the whole.
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
my $format = 'pp';
my $raw_output;
my $num_procs = 48;

GetOptions( "novel|n"             => \$novel,
            "copy-number|cn=f"    => \$copy_number,
            "upper-copies|cu=f"   => \$cu,
            "lower-copies|cl=f"   => \$cl,
            "gene|g=s"            => \$geneid,
            "tiles|t=i"           => \$tiles,
            "annot|a"             => \$annot,
            "output|o=s"          => \$outfile,
            "format|f=s"          => \$format,
            "raw|r"               => \$raw_output,
            "NOCALL|N"            => \$nocall,
            "procs|p=i"       => \$num_procs,
            "version|v"           => \$ver_info,
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
        if ($keys eq 'gene') {
            printf "\t%-7s => %s\n",$keys,join(',',@$values);
        } else {
            printf "\t%-7s => %s\n",$keys,$values;
        }
    }
    print '='x79, "\n";
}

my %formats = (
    'csv'   => ',',
    'tsv'   => "\t",
    'pp'    => '',
);

# Set the output format delimiter
my $delimiter;
unless (defined $formats{$format}) {
    die "ERROR: '$format' is not a valid option as a delimiter!\n";
}
($format) ? ($delimiter = $formats{$format}) : ($delimiter = $formats{pp});

# Make sure enough args passed to script
if ( scalar( @ARGV ) < 1 ) {
    print "ERROR: No VCF files passed to script!\n\n"; 
    print "$usage\n";
    exit 1;
}

# Write output to either indicated file or STDOUT
my $out_fh;
if ( $outfile ) {
	open( $out_fh, ">", $outfile );
} else {
	$out_fh = \*STDOUT;
}

#########--------------------- END ARG Parsing -----------------------#########
my %cnv_data;
my @vcfs = @ARGV;

my $pm = Parallel::ForkManager->new($num_procs);
$pm->run_on_finish(
    sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, 
            $data_structure_reference) = @_;
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

#dd \%cnv_data;
# exit;

my %results;
for my $sample ( keys %cnv_data ) {
    $results{$sample} = [];
    my @outfields = qw( END LEN NUMTILES RAW_CN REF_CN CN HS FUNC );
    for my $cnv ( sort { versioncmp ( $a, $b ) } keys %{$cnv_data{$sample}} ) {
        my %mapped_cnv_data;
        last if $cnv eq 'NONE';

        # Seems to be a bug in the same the CI are reported for deletions.
        # Solution correctly reports the value in the VCF, but it's not so 
        # informative.  This will give a better set of data.
        my ($ci_5, $ci_95) = $cnv_data{$sample}->{$cnv}->{'CI'} =~ /0\.05:(.*?),0\.95:(.*)$/; 
        my ($chr, $start, $gene, undef) = split( /:/, $cnv );
        %mapped_cnv_data = map{ $_ => $cnv_data{$sample}->{$cnv}->{$_} } @outfields;
        @mapped_cnv_data{qw(ci_05 ci_95)} = (sprintf("%.2f",$ci_5),
            sprintf("%.2f",$ci_95));
        @mapped_cnv_data{qw(chr start gene undef)} = split(/:/, $cnv);

        $mapped_cnv_data{HS}    //= 'No';
        $mapped_cnv_data{ci_05} //= 0;
        $mapped_cnv_data{ci_95} //= 0;

        # Get OVAT Annot Data
        my ($gene_class, $variant_class);
        my $func = $mapped_cnv_data{FUNC};
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
        push(@{$results{$sample}}, \@filtered_data) if @filtered_data;
    }
}
print_results(\%results, $delimiter);

sub print_results {
    my ($data, $delimiter) = @_;
    my @header = qw( Chr Gene Start End Length Tiles Raw_CN Ref_CN CI_05 CI_95
        CN Annot );
    my $pp_format = "%-8s %-8s %-11s %-11s %-11s %-8s %-8s %-8s %-8s %-8s %-8s 
        %-18s\n";

    select $out_fh;

    # Print out comma separated dataset for easy import into Excel and whatnot.
    raw_output($data) if $raw_output;

    for my $sample (keys %$data) {
        my ($id, $gender, $mapd, $cellularity) = split( /:/, $sample );
        print "::: CNV Data For $id (Gender: $gender, Cellularity: ", 
            "$cellularity, MAPD: $mapd) :::\n";
        ($delimiter) 
            ? print join($delimiter, @header),"\n" 
            : printf $pp_format, @header;
        if ( ! @{$$data{$sample}} ) {
            print ">>>>  No Reportable CNVs Found in Sample  <<<<\n"; 
        } else {
            for my $cnv (@{$$data{$sample}}) {
                ($delimiter) 
                    ? print join($delimiter, @$cnv), "\n" 
                    : printf $pp_format, @$cnv;
            }
        }
        print "\n";
    }
    return;
}

sub raw_output {
    my $data = shift;
    my @header = qw( Sample Gender Cellularity MAPD Chr Gene Start End Length 
        Tiles Raw_CN Ref_CN CI_05 CI_95 CN Annot );
    select $out_fh;
    print join(',', @header), "\n";

    for my $sample (keys %$data) {
        my @elems = split(/:/, $sample);
        for my $cnv (@{$$data{$sample}}) {
            print join(',', @elems, @$cnv), "\n";
        }
    }
    exit;
}

sub filter_results {
    # Filter out CNV data prior to printing it all out.
    my ($data, $filters) = @_;
    my @cn_thresholds = @$filters{qw(cn cu cl)};

    # Gene level filter
    return if (@{$filters{gene}}) and ! grep {$$data{gene} eq $_} @{$filters{gene}};

    # Filter non-hotspots and novel if we don't want them.
    return if ! $$filters{novel} and ($$data{HS} eq 'No' || $$data{gene} eq '.');

    # Number of tiles filter
    return if ($$filters{tiles} and $$data{NUMTILES} < $$filters{tiles});

    # OVAT Filter
    return if ($$filters{annot} and $$data{GC} eq '---');
    
    # We made it the whole way through; check for copy number thresholds
    (copy_number_filter($data, \@cn_thresholds)) 
        ? return return_data($data) 
        : return;
}

sub return_data {
    my $data = shift;
    my @fields = qw(chr gene start END LEN NUMTILES RAW_CN REF_CN ci_05 ci_95 
        CN GC);
    return @$data{@fields};
}

sub copy_number_filter {
    # if there is a 5% / 95% CI filter, use that, otherwise if there's a cn 
    # filter use that, and if there's nothing, just return it all.
    my ($data, $threshold) = @_;
    my ($cn, $cu, $cl) = @$threshold;

    if ($cn) {
        return 1 if ($$data{CN} >= $cn);
    }
    elsif ($cu) {
        return 1 if $cl and ($$data{ci_05} >= $cu || $$data{ci_95} <= $cl);
        return 1 if ($$data{ci_05} >= $cu);
    }
    elsif ($cl) {
        return 1 if $cu and ($$data{ci_05} >= $cu || $$data{ci_95} <= $cl);
        return 1 if ($$data{ci_95} <= $cl);
    } else {
        # Return everything if there are no filters.
        return 1;
    }
    return 0;
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
            # Need to add to accomodate the new CNV plugin; may not have the 
            # same field as the normal IR data.
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
        $sample_id = join( ':', $sample_name, $gender, $mapd, $cellularity );

        next unless $data[4] eq '<CNV>';
        # Let's handle NOCALLs for MATCHBox compatibility (prefer to filter on 
        # my own though).
        if ($nocall && $data[6] eq 'NOCALL') {
            ${$cnv_data{$sample_id}->{'NONE'}} = '';
            next;
        }

        my $varid = join( ':', @data[0..3] );
        
        # Will either get a HS entry in the line, or nothing at all. Kludgy, but
        # need to deal with hotspots (HS) field; not like others (i.e. not a 
        # key=val struct)!
        $data[7] =~ s/HS/HS=Yes/;

        # New v5.2 standard deviation value; sometimes data in this field and 
        # sometimes not. 
        $data[7] =~ s/SD;/SD=NA;/; 

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
    return \%results, \$sample_id;
}

sub __exit__ {
    my ($line, $msg) = @_;
    print "\n\n";
    print colored("Got exit message at line: $line with message:\n$msg", 
        'bold white on_green');
    print "\n";
    exit;
}
