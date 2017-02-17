#!/usr/bin/perl
# Load up a JSON of genes and categories for the RNA Fusions panel (the GeneExpression assays), along with a VCF
# file, and do out parsing and outputting based on those data.
# 2/9/2017 - D Sims
###############################################################################################################
use warnings;
use strict;
use autodie;

use File::Basename;
use Data::Dump;
use JSON;

my $version = '0.1.0_020917';

sub load_json {
    my $json_file = dirname($0) . '/gene_expression_assays.json';
    #print "json: $json_file\n";
    #exit;
    die "ERROR: Can't find the 'gene_expression_assays.json' file!\n" unless -f $json_file;

    my $jtext = do {
        open( my $json_fh, "<", $json_file);
        local $\;
        <$json_fh>;
    };
    my $jobj = JSON->new;
    my $jdata = $jobj->decode($jtext);
    return $jdata;
}

# Get JSON of genes and load up for easy lookup later
#my $genes_json = shift;
my $genes = load_json();
# dd $genes;

# load and parse the VCF
my $vcf = shift;
(my $sample_name = $vcf) =~ s/\.vcf//;
die "error: you must load at least one VCF file!\n" unless $vcf;

my %results;
open(my $fh, "<", $vcf);
while (<$fh>) {
    my @data = split;
    next unless /SVTYPE=GeneExpression/;
    my ($counts) = $_ =~ /READ_COUNT=(\d+)/;
    my @fields = split(/[-.]/,$data[2]);
    push(@{$results{$fields[0]}}, {$data[2] => $counts});
}

# Dump it all out depending on gene.  will add gene expression vs WT vs splice later based on JSON input.
#my $query_gene = 'EGFR';
#my @genes = qw(ALK AR BRAF BRCA1 BRCA2 CDKN2A EGFR ERBB2 MDM4 MET NTRK1 RB1);
my @genes = qw(EGFR);

my $count = 0;
for my $query_gene (@genes) {
    for my $assay (@{$results{$query_gene}}) {
        while ( my ($k,$v) = each %$assay ) {
            #print join(": ", $k, $v ), "\n";
            print join(",", $sample_name, $k, $v ), "\n";
            $count += $v;
        }
    }
}
print "$sample_name,total,$count\n";
