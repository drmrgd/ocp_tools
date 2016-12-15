#!/usr/bin/perl
# Script to read in one or more VCF files from a MATCH Control run and output a report.  Eventually I'll link this with R so that we
# can generate a nice box and whisker plot
# D Sims - 8/15/16
####################################################################################################################################
use warnings;
use strict;
use autodie;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Data::Dump;
use Parallel::ForkManager;
use Term::ANSIColor;
use Sort::Versions;

my $scriptname = basename($0);
my $version = "v1.4.2_101216";
my $description = <<"EOT";
Generate a summary MATCH control report.  Need to input a list of VCF files and the version of the MATCH control used.  By
default we'll use version 2.  Also, can output as a pretty printed report, or can output as a CSV for importation into 
other tools like R.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <vcf_file(s)>
    -l, --lookup    Version of lookup table to use (1 or 2). DEFAULT: 2.
    -s, --site      Manually chose the MATCH site rather than deducing it from the file data.
    -f, --format    Method to format the report output (pp: pretty print, csv: CSV file).  DEFAULT: 'pp'.
    -o, --output    Send output to custom file.  Default is STDOUT.
    -v, --version   Version information
    -h, --help      Print this help information
EOT

my $help;
my $ver_info;
my $outfile;
my $lookup_table = 2;
my $format = 'pp';
my $sequencing_site;

GetOptions( "lookup|l=s"    => \$lookup_table, 
            "format|f=s"    => \$format,
            "output|o=s"    => \$outfile,
            "site|s=s"      => \$sequencing_site,
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
    print "ERROR: Not enough arguments passed to script!\n\n";
    print "$usage\n";
    exit 1;
}

die "ERROR: You must choose either 'pp' or 'csv' for report format with the '-f' option!\n" if ($format ne 'pp' && $format ne 'csv');

# Validate site input
my @valid_sites = qw(NCI MDA YSM MGH);
if ($sequencing_site) {
    $sequencing_site = uc($sequencing_site);
    die "ERROR: Invalid site '$sequencing_site'!\n" unless grep { $sequencing_site eq $_ } @valid_sites;
}

# Write output to either indicated file or STDOUT
my $out_fh;
if ( $outfile ) {
	open( $out_fh, ">", $outfile ) || die "Can't open the output file '$outfile' for writing: $!";
    print "Writing output to '$outfile'\n";
} else {
	$out_fh = \*STDOUT;
}

#########------------------------------ END ARG Parsing ---------------------------------#########
my @vcfs = @ARGV;

my %control_data;
my $pm = new Parallel::ForkManager(48);
$pm->run_on_finish(
    sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
        my $vcf = $data_structure_reference->{input};
        my ($name) = $vcf =~ /^.*?(SampleControl_.*?_\d+)(?:_v[0-9]+)?.*/;
        $name //= basename($vcf);
        $control_data{$name} = $data_structure_reference->{result};
    }
);

for my $vcf (@vcfs) {
    $pm ->start and next;
    my $return_data = proc_vcf(\$vcf);
    $pm->finish(0, { result => $return_data, input => $vcf });
}
$pm->wait_all_children;

print "INFO: Using lookup table for MATCH control v$lookup_table...\n";
check_results(\%control_data, $lookup_table);
generate_report(\%control_data, $format);

sub check_results {
    my ($data,$lookup_table) = @_;

    my %v2_lookup_table = (
        'chr3:178916946:G:C:PIK3CA'  => '',
        'chr7:140453136:A:T:BRAF'    => '',
        'chr10:89717715:T:TA:PTEN'   => '',
        'chr13:32968850:C:A:BRCA2'   => '',
        'chr13:48916815:CACTT:C:RB1' => '',
        'chr17:7574002:CG:C:TP53'    => '',
        'chr17:ERBB2'                => '',
        'chr17:RPS6KB1'              => '',
        'chr20:ZNF217'               => '',
        'ALK-PTPN3.A11P3'            => '',
        'EML4-ALK.E6aA20'            => '',
        'MET-MET.M13M15'             => '',
    );

    my %v1_lookup_table  = (
        'chr3:178916946:G:C:PIK3CA'  => '',
        'chr7:140453136:A:T:BRAF'    => '',
        'chr10:89717716:T:TA:PTEN'   => '',
        'chr13:32968850:C:A:BRCA2'   => '',
        'chr13:48916815:CACTT:C:RB1' => '',
        'chr17:7574002:CG:C:TP53'    => '',
        'chr17:ERBB2'                => '',
        'chr17:RPS6KB1'              => '',
        'chr20:ZNF217'               => '',
        'EML4-ALK.E6aA20'            => '',
        'EML4-ALK.E6bA20'            => '',
    );

    my %lookup_tables = (
        1  => \%v1_lookup_table,
        2  => \%v2_lookup_table,
    );

    unless (defined $lookup_tables{$lookup_table}) {
        print "ERROR: Lookup table '$lookup_table' is not a valid lookup table!  Valid tables are:\n";
        #print "\t$_  => $lookup_tables{$_}\n" for keys %lookup_tables;
        print "\t$_  => v$_\n" for sort keys %lookup_tables;
        exit 1;
    }

    my %results;
    my %truth_table = %{$lookup_tables{$lookup_table}};

    for my $sample (keys %$data) {
        for my $variant (keys %{$$data{$sample}}) {
            if ( exists $truth_table{$variant} ) {
                $truth_table{$variant} = 1;
                push(@{$$data{$sample}->{$variant}}, 'POS');
            } else {
                push(@{$$data{$sample}->{$variant}}, 'FP');
            }
        }

        # Look for FN calls too.
        for my $control (keys %truth_table) {
            my $num_elems = () = $control =~ /:/g;
            my @var_string;
            if ($num_elems == 4) {
                my ($pos, $ref, $alt, $gene) = $control =~ /(^chr.*):(\w+):(\w+):(.*?)$/;
                @var_string = ('SNV_Indel', $pos, $ref, $alt, 'vaf', 'totcov', 'refcov', 'altcov', 'varid', $gene);
                push(@var_string, ('xxx')x7);
            }
            elsif ($num_elems == 1) {
                my @var_elems = split(/:/, $control);
                @var_string = ('CNV', $var_elems[1], $var_elems[0], 'numtiles', '5%ci', 'cn', '95%ci', 'mapd');
            }
            else {
                (my $driver) = map { $_ =~ /(ALK|MET)/ } split(/:/,$control);
                @var_string = ('Fusion', $control, 'id', 'counts', $driver, 'partner');
            }
            push(@{$$data{$sample}->{$control}}, @var_string, 'NEG') unless $$data{$sample}->{$control};
        }
    }
    return;
}

sub generate_report {
    my ($data,$format) = @_;
    my %want_fields = (
        'SNV_Indel'    => [qw(0 9 1 2 3 4 5 17)],
        'CNV'    => [qw(0 1 2 5 8)],
        'Fusion' => [qw(0 4 1 3 6)],
    );

    my @output_strings;
    my $sample_width = get_width([keys %$data]);
    my @header = qw(Sample Site VarID Type Gene Position Ref Alt VAF_CN Cov_Reads Measurement Call);
    push(@output_strings, format_report_string(\$format, $sample_width, \@header));

    for my $sample ( sort{versioncmp($a, $b)} keys %$data ) {
        my $site;
        ($sequencing_site)
            ? ($site = $sequencing_site)
            : ($site = get_site_name($sample));

        # If no input and no deduced site, then just print a placeholder string
        warn "WARN: No sequencing site info available!\n" if $site eq '---';

        for my $variant (sort{ $$data{$sample}->{$b}[0] cmp $$data{$sample}->{$a}[0] } keys $$data{$sample}) {
            my $type = $$data{$sample}->{$variant}[0];
            my @wanted_indices = @{$want_fields{$type}};
            my @var_data = ($sample, $site, $variant, @{$$data{$sample}->{$variant}}[@wanted_indices]);

            if ($type eq 'Fusion') {
                splice(@var_data, 5, 0, '---','---');
                splice(@var_data, 8, 0, '---');
            }
            elsif ($type eq 'CNV') {
                splice(@var_data, 6, 0, '---', '---');
                splice(@var_data, 9, 0, '---');
            }
            # Handle negative results.
            map {$_ = '0'} @var_data[8,9] if $var_data[10] eq 'NEG';
            ($var_data[3] =~ /[SC]NV/) ? splice(@var_data,10,0,$var_data[8]) : splice(@var_data,10,0,$var_data[9]);
            push(@output_strings, format_report_string(\$format, $sample_width, \@var_data));
        }
    }
    print {$out_fh} join("\n", @output_strings), "\n";
}

sub get_site_name {
    my $sample_string = shift;
    my %sites = (
        'MoCha'  => 'NCI',
        'MDACC'  => 'MDA',
        'MGH',   => 'MGH',
        'Yale'   => 'YSM',
        'NA'     => '---',
    );
    (my $site) = $sample_string =~ /^SampleControl_(\w+?)_\d+/;
    $site //= 'NA';
    return $sites{$site};
}

sub proc_vcf {
    my ($vcf, $name) = @_;
    my %results;

    # Don't output these calls since we want them filtered anyway.
    my @filtered_variants = qw( chr17:7579473:G:C:TP53 );
    #push(@filtered_variants, 'EML4-ALK.E6bA20') if $lookup_table == 2;
    push(@filtered_variants, 'EML4-ALK.E6bA20') if $lookup_table =~ /[23]/; 
    my $cmd = qq(match_moi_report.pl -n -r -R1000 -c7 $$vcf) ;

    open( my $moi_report_pipe, '-|', $cmd);
    while (<$moi_report_pipe>) {
        chomp;
        #next unless /^chr/ || /^[-\w\d]+\.[^\.]+\b/;
        next if /NOTE/;
        my @data = split(/,/);
        my $varid; 
        if ($data[0] eq 'SNV') {
            $data[0] = 'SNV_Indel';
            next if $data[7] < 25; # Routbort rule
            $varid = join(':', @data[1..3,9]);
        }
        elsif ($data[0] eq 'CNV') {
            next if $data[-1] > 0.5; # Get rid of failed panel due to MAPD
            $varid = join(':', @data[2,1]);
        }
        elsif ($data[0] eq 'Fusion') {
            $varid = $data[1];
        }
        $results{$varid} = [@data] unless grep { $varid eq $_ } @filtered_variants;
    }
    return \%results;
}

sub get_width {
    my $sample_names = shift;
    my $width = 0;
    for (@$sample_names) {
        $width = length($_) if length($_) > $width;
    }
    return ($width+2);
}

sub format_report_string {
    my ($format, $width, $data) = @_;
    # Sample Site VarID Type Gene Position Ref Alt VAF_CN Cov_Reads Measurement Call
    my $pp_format = "%-${width}s %-5s %-27s %-11s %-10s %-20s %-7s %-20s %-8s %-11s %-11s %-6s";
    ($$format eq 'pp') ? return sprintf($pp_format, @$data) : return join(',', @$data);
}

sub __exit__ {
    my ($line, $msg) = @_;
    print "\n\n";
    print colored("Got exit message at line: $line with message:\n$msg", 'bold white on_green');
    print "\n";
    exit;
}
