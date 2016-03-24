#!/usr/bin/perl

=head1 NAME

gff2pep.pl - Convert GFF3 and nucleotide Fasta file to translated a.a. sequences

=head1 SYNOPSIS

perl gff2pep.pl -i <fasta_nucl> -g <GFF> -t <genetic_code> -o <output fasta>

perl gff2pep.pl --help

=head1 DESCRIPTION

Read features from GFF3 file, get nucleotide sequences from Fasta file and
convert to a.a. sequences. Requires BioPerl

=head1 ARGUMENTS

=over 8

=item --fasta_in|-i <file>

Nucleotide sequences in Fasta format

=item --gff|-g <file>

Feature table in GFF3 format

=item --fasta_out|-o <file>

Name for output file (Fasta)

=item --transl_table|-t <integer>

Translation table (genetic code) to use. (Default: 4)

=item --help|-h

This help message

=back

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2015, Brandon Seah

... GPL-3 ...

=cut

# convert GFF and Fasta nucleotide file to translated peptide sequences

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

use Bio::DB::Fasta;
use Bio::Seq;
use Bio::SeqIO;

my $fasta_file;
my $gff_file;
my $output_fasta;
my $transl_table=4;
my %feature_hash;

GetOptions ('fasta_in|i=s'=>\$fasta_file,
            'gff|g=s' => \$gff_file,
            'fasta_out|o=s'=>\$output_fasta,
            'transl_table|t=i'=>\$transl_table,
            'help|h' => sub { pod2usage(-exitstatus=>2, -verbose=>2); },
            'man|m' => sub { pod2usage(-exitstatus=>0, -verbose=>2); }
            ) or pod2usage(-message=>"Invalid input", -exitstatus=>2);

if (!defined $fasta_file || !defined $gff_file) {
    pod2usage (-message=>"Invalid input", -exitstatus=>2);
}

read_gff_file();
test();


sub read_gff_file {
    open(GFF_IN, "<", $gff_file) or die ("$!\n");
    while (<GFF_IN>) {
        if ($_ =~ /^\#/) { next; } # skip comment lines
        chomp;
        my @theline = split "\t", $_;
        (my $feature_id) = $theline[8] =~ m/ID=(\S+)/;  # Extract ID field from feature table notes field
        $feature_hash{$feature_id}{"contig"} = $theline[0];
        $feature_hash{$feature_id}{"type"} = $theline[2];
        $feature_hash{$feature_id}{"start"} = $theline[3];
        $feature_hash{$feature_id}{"end"} = $theline[4];
        if ($theline[6] eq "+") { $feature_hash{$feature_id}{"strand"} = 1; }
        elsif ($theline[6] eq "-") { $feature_hash{$feature_id}{"strand"} = -1; }
    }
    close (GFF_IN);
}

sub test {
    my $db = Bio::DB::Fasta->new($fasta_file);
    my $fasta_out = Bio::SeqIO->new(-file=> ">$output_fasta", -format=>'fasta');

    foreach my $thekey (keys %feature_hash) {
        if ($feature_hash{$thekey}{"type"} ne "CDS") {
            next;
        }
        my $contig = $feature_hash{$thekey}{"contig"};
        my $start = $feature_hash{$thekey}{"start"};
        my $end = $feature_hash{$thekey}{"end"};
        my $strand = $feature_hash{$thekey}{"strand"};
        my $seq = $db->seq("$contig:$start..$end/$strand");
        my $nt_out = Bio::Seq->new(-seq=>$seq,
                                   -display_id=>$thekey,
                                   -alphabet=>'dna');
        my $aa_out = $nt_out->translate(-codontable_id=>$transl_table);
        $fasta_out->write_seq($aa_out);
        #print $seq;
        #print $contig."\t".$start."\t".$end;
        #print "\n";
    }
}