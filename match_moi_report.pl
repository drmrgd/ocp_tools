#!/usr/bin/perl
# Parse a VCF file and generate a table of MOIs and aMOIs
# TODO:
#   - Capture NOCALLS into a hash to print out later (maybe CLI opt to output?)
#   - Set up CLI opts to set VAF, CN, and Fusion counts as thresholds (fusions use the EGFRvIII vs Others rule)
#   - Set up CLI opts to either print a fancy report, or dump to CSV with a column for variant type
#   - Add field to SNV and Indel section to id what rule caused the variant to make the report (Hotspot, EGFR Exon 19, etc)
#
# 2/12/2014 - D Sims
##################################################################################################################################################
use warnings;
use strict;
use autodie;
use feature "switch";
no if $] > 5.018, 'warnings', 'experimental::smartmatch'; # Need to turn of smartmatch warning, but only if the version of perl is newer

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Sort::Versions;
use Data::Dump;

my $scriptname = basename($0);
my $version = "v1.0.1_022015";
my $description = <<"EOT";
Program to parse an IR VCF file to generate a list of NCI-MATCH MOIs and aMOIs.  This program requires the use of `convert_vcf.py` from 
ThermoFisher to run as it does the bulk of the file parsing.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <VCF>
    -o, --output    Send output to custom file.  Default is STDOUT.
    -v, --version   Version information
    -h, --help      Print this help information
EOT

my $help;
my $ver_info;
my $outfile;

GetOptions( "output|o=s"    => \$outfile,
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
    print "ERROR: No VCF file passed to script\n\n!";
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
my $vcf_file = shift;
my ($gender, $cellularity, $mapd);
my $variant_data = read_vcf(\$vcf_file );

my (%snv_indel_data, %fusion_data, %cnv_data);
for my $variant (@$variant_data) {
    my $type = $$variant{'rowtype'};
    given ($type) {
        when (/([ms]np|ins|del|complex)/) { proc_snv_indel( $variant ) }
        when (/Fusion/)                   { proc_fusion( $variant ) }
        when (/CNV/)                      { proc_cnv( $variant ) }
        default                           {next}
    }
}

gen_report( \%snv_indel_data, \%fusion_data, \%cnv_data );


sub read_vcf {
    # Read in VCF file, proc with 'convert_vcf.py' and load a data struct.
    my $input_file = shift;
    (my $converted_vcf = $$input_file) =~ s/vcf$/tsv/;
    my @variant_results;

    # Want to have MAPD, Gender, and Cellularity in the output.  So, going to have to read the VCF twice it looks like
    open( my $vcf_fh, "<", $$input_file );
    while (<$vcf_fh>) {
        if ( /^##/ ) {
            if ( /sampleGender=(\w+)/ ) {
                $gender = $1; 
                next;
            } 
            elsif ( /CellularityAsAFractionBetween0-1=(.*)$/ ) {
                $cellularity = sprintf( "%d%%", 100 * $1); 
                next;
            }
            elsif ( /mapd=(\d\.\d+)/ ) {
                $mapd = $1; 
                next;
            }
        }
    }
    close $vcf_fh;

    qx( convert_vcf.py --force -i $$input_file -o $converted_vcf );
    open( my $fh, "<", $converted_vcf );
    chomp( my @header = split(/\t/, <$fh>) );
    while (<$fh>) {
        my %raw_data;
        chomp( my @elems = split(/\t/) );
        @raw_data{@header} = @elems;
        my %filtered = filter_raw_data( \%raw_data );
        push(@variant_results, \%filtered);
    }
    close $fh;
    unlink $converted_vcf;
    return \@variant_results;
}

sub filter_raw_data {
    # Pass in the raw data from the 'read_vcf()' funct and cull out what we don't need.
    my $raw_data = shift;
    my %filtered_data;

    # Use hash instead of array of keys in order to prevent undef key warning
    my %wanted_keys = ( 
    'rowtype'  => '',
    'call'    => '',
    'ALT'   => '',
    'CHROM'   => '',
    'FILTER'   => '',
    'ID' => '',
    'INFO...FR'   => '',
    'INFO...OALT'   => '',
    'INFO...OID'   => '',
    'INFO...OMAPALT'   => '',
    'INFO...OPOS'   => '',
    'INFO...OREF'   => '',
    'INFO...READ_COUNT'   => '',
    'INFO.1.FRO'   => '',
    'INFO.1.NUMTILES' => '',
    'INFO.A.AF'   => '',
    'INFO.A.FAO'   => '',
    'POS'   => '',
    'REF'   => '',
    'FORMAT.1.CN'   => '',
    'INFO...CI'   => '',
    'INFO.1.RO'   => '',
    'INFO.A.AO'   => '',
    'FUNC1.coding'   => '',
    'FUNC1.exon'   => '',
    'FUNC1.gene'   => '',
    'FUNC1.normalizedAlt'   => '',
    'FUNC1.normalizedPos'   => '',
    'FUNC1.normalizedRef'   => '',
    'FUNC1.oncomineGeneClass'   => '',
    'FUNC1.oncomineVariantClass'  => '',
    'FUNC1.protein'   => '',
    'FUNC1.transcript'   => '',
    'FUNC1.location'   => '',
    'FUNC1.function'  => '',
    );
    
    @filtered_data{keys %wanted_keys} = @$raw_data{keys %wanted_keys};
    return %filtered_data;
}

sub proc_snv_indel {
    # Follow rules to generate a list of MOIs and aMOIs
    my $variant_info = shift;

    return unless ( $$variant_info{'FUNC1.location'} eq 'exonic' || $$variant_info{'FUNC1.function'} eq 'synonymous' );
    my $id = join( ':', $$variant_info{'CHROM'}, $$variant_info{'INFO...OPOS'}, $$variant_info{'INFO...OREF'}, $$variant_info{'INFO...OALT'} ); 

    if ( $$variant_info{'INFO.A.AF'} > 0 ) {
        # Anything that's a hotspot
        if ( $$variant_info{'INFO...OID'} ne '.' ) {
            # bin NOCALLs for now
            return if ( $$variant_info{'call'} eq 'NOCALL' );
            gen_var_entry( $variant_info, \$id );
        }
        # De Novo TSG frameshift calls
        elsif ( $$variant_info{'FUNC1.oncomineGeneClass'} ) {
            gen_var_entry( $variant_info, \$id );
        }
        # EGFR nonframeshiftDeletion in Exon 19 rule
        elsif ( $$variant_info{'FUNC1.gene'} eq 'EGFR' && $$variant_info{'FUNC1.exon'} eq '19' && $$variant_info{'FUNC1.function'} eq 'nonframeshiftDeletion' ) {
            gen_var_entry( $variant_info, \$id );
        }
        # ERBB2 nonframeshiftInsertion in Exon20 rule
        elsif ( $$variant_info{'FUNC1.gene'} eq 'ERBB2' && $$variant_info{'FUNC1.exon'} eq '20' && $$variant_info{'FUNC1.function'} eq 'nonframeshiftInsertion' ) {
            gen_var_entry( $variant_info, \$id );
        }
    }
    return;
}

sub gen_var_entry {
    my $data = shift;
    my $id = shift;
    
    my $coord = "$$data{'CHROM'}:$$data{'INFO...OPOS'}";
    my $ref = $$data{'INFO...OREF'};
    my $alt = $$data{'INFO...OALT'};
    my $filter = $$data{'FILTER'};
    (my $fr = $$data{'INFO...FR'}) =~ s/^\.,//;
    my $vaf = sprintf( "%0.2f", ($$data{'INFO.A.AF'} * 100) );
    my ($rcov, $acov, $tcov);
    if ( $$data{'INFO.A.FAO'} eq '.' ) {
        $rcov = $$data{'INFO.1.RO'};
        $acov = $$data{'INFO.A.AO'}; 
    } else {
        $rcov = $$data{'INFO.1.FRO'};
        $acov = $$data{'INFO.A.FAO'}; 
    }
    $tcov = $rcov + $acov;
    my $varid = $$data{'ID'};
    my $gene = $$data{'FUNC1.gene'};
    my ($om_gc, $om_vc);
    ($$data{'FUNC1.oncomineGeneClass'}) ? ($om_gc = $$data{'FUNC1.oncomineGeneClass'}) : ($om_gc = '---');
    ($$data{'FUNC1.oncomineVariantClass'}) ? ($om_vc = $$data{'FUNC1.oncomineVariantClass'}) : ($om_vc = '---');

    # Remove duplicate variant entries that span VCF blocks
    if ($varid ne '.' && exists $snv_indel_data{$$id}) {
        delete $snv_indel_data{$$id};
    } 

    $snv_indel_data{$$id} = [$coord, $ref, $alt, $filter, $fr, $vaf, $tcov, $rcov, $acov, $varid, $gene, $om_gc, $om_vc];
    return;
}

sub proc_fusion {
    # Generate Hash of fusion data for output later.
    my $variant_info = shift;
    if ( $$variant_info{'call'} eq 'POS' ) { 
        my ($name, $elem) = $$variant_info{'ID'} =~ /(.*?)_([12])/;
        my ($pair, $junct, $id) = split( /\./, $name );
        $id = 'NOVEL' unless $id;
        my $fid = join( '|', $pair, $junct, $id );
        
        $fusion_data{$fid}->{'COUNT'} = $$variant_info{'INFO...READ_COUNT'};
        if ( $elem == 2 ) {
            $fusion_data{$fid}->{'DRIVER'} = $$variant_info{'FUNC1.gene'};
        } else {
            $fusion_data{$fid}->{'PARTNER'} = $$variant_info{'FUNC1.gene'};
        }
    }
    return;
}

sub proc_cnv {
    my $variant_info = shift;

    if ( $$variant_info{'FORMAT.1.CN'} >= 7 ) { 
        my ($ci_05, $ci_95) = $$variant_info{'INFO...CI'} =~ /0\.05:(.*?),0\.95:(.*)/;
        $cnv_data{$$variant_info{'FUNC1.gene'}} = [$$variant_info{'CHROM'}, $$variant_info{'INFO.1.NUMTILES'}, $ci_05, $$variant_info{'FORMAT.1.CN'}, $ci_95];
    }
    return;
}

sub field_width {
    # Get the longest field width for formatting later.
    my $data_ref = shift;
    my $type = shift;

    my $ref_width = 0;
    my $var_width = 0;
    my $filter_width= 0;
    my $fusion_id_width = 0;

    if ( $type eq 'snv' ) {
        for my $variant ( keys %$data_ref ) {
            my $ref_len = length( $$data_ref{$variant}[1] );
            my $alt_len = length( $$data_ref{$variant}[2] );
            my $filter_len = length( $$data_ref{$variant}[4] );
            $ref_width = $ref_len if ( $ref_len > $ref_width );
            $var_width = $alt_len if ( $alt_len > $var_width );
            $filter_width = $filter_len if ( $filter_len > $filter_width );
        }

        ( $filter_width > 13 ) ? ($filter_width += 4) : ($filter_width = 17);
        return ( $ref_width + 4, $var_width + 4, $filter_width);
    } 
    elsif ( $type eq 'fusion' ) {
        for my $variant ( keys %$data_ref ) {
            my @id_elems = split( /\|/, $variant );
            my $length = length( "$id_elems[0].$id_elems[1]" );
            $fusion_id_width = $length if $length > $fusion_id_width;
        }
        return ( $fusion_id_width + 4 );
    }
    return;
}

sub gen_report {
    # Print out the final MOI Report
    my ($snv_indels, $fusion_data, $cnv_data) = @_;

    my ($w1, $w2, $w3);

    select $out_fh;

    print "::: MATCH Reportable SNVs and Indels :::\n";
    ($w1, $w2, $w3) = field_width( $snv_indels, 'snv' );
    my $snv_indel_format = "%-17s %-${w1}s %-${w2}s %-10s %-${w3}s %-8s %-10s %-10s %-10s %-14s %-10s %-21s %s\n";
    my @snv_indel_header = qw( Chrom:Pos Ref Alt Filter Filter_Reason VAF TotCov RefCov AltCov VARID Gene oncomineGeneClass oncomineVariantClass );
    printf $snv_indel_format, @snv_indel_header;
    if ( %$snv_indels ) {
        for my $variant ( sort{ versioncmp( $a, $b ) } keys %$snv_indels ) {
            printf $snv_indel_format, @{$$snv_indels{$variant}};
        }
    } else {
        print ">>>>  No Reportable SNVs or Indels Found in Sample  <<<<\n";
    }
    print "\n";

    print "::: MATCH Reportable Fusions :::\n";
    ($w1) = field_width( $fusion_data, 'fusion' );
    my $fusion_format = "%-${w1}s %-12s %-12s %-15s %-15s\n";
    my @fusion_header = qw( Fusion ID Read_Count Driver_Gene Partner_Gene );
    printf $fusion_format, @fusion_header;
    if ( %$fusion_data ) {
        for ( sort { versioncmp( $a, $b ) } keys %$fusion_data ) {
            my ($fusion, $junct, $id) = split( /\|/ );
            printf $fusion_format, "$fusion.$junct", $id, $fusion_data{$_}->{'COUNT'}, $fusion_data{$_}->{'DRIVER'}, $fusion_data{$_}->{'PARTNER'};
        }
    } else {
        print ">>>>  No Reportable Fusions found in Sample  <<<<\n";
    }
    print "\n";

    print "::: MATCH Reportable CNVs (Gender: $gender, Cellularity: $cellularity, MAPD: $mapd) :::\n";
    my $cnv_format = "%-9s %-10s %-6s %-10.3f %-10.1f %-10.3f\n";
    my @cnv_header = qw( Chr Gene Tiles CI_05 CN CI_95 );
    printf "%-9s %-10s %-6s %-10s %-10s %-10s\n", @cnv_header;
    if ( %$cnv_data ) {
        for my $cnv ( sort{ versioncmp( $a, $b ) } keys %$cnv_data ) {
            printf $cnv_format, $cnv_data{$cnv}->[0], $cnv, @{$$cnv_data{$cnv}}[1..4];
        }
    } else {
        print ">>>>  No Reportable CNVs Found in Sample  <<<<\n";
    }
    print "\n";
    return;
}
