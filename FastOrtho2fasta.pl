#!/usr/bin/env perl

# Parse FastOrtho results file (.end) to get list of clusters found only once per genome
# Align each cluster with Muscle

use strict;
use warnings;

use Bio::Align::Utilities qw(:all);
use Bio::LocatableSeq;
use Bio::DB::Fasta;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Getopt::Long;
use Cwd;

my $fastortho_result;
my $fasta_db;
my $concat_file;
my $outpath = getcwd();

GetOptions ('fastortho_end|e=s'=>\$fastortho_result,    # FastOrtho .end file
            'fasta|f=s'=>\$fasta_db,                    # Concatenated file of aa sequences
            'output|o=s'=>\$outpath                     # Output path (if not current folder)
            ) or die ("$!\n");

## MAIN #######################################################################
my %cluster_hash;           # Store the protein IDs that correspond to each cluster
my %clusters_with_paralogs; # Clusters that contain members from same genome (i.e. putative paralogs)
my %alignment_hash;
my %genomes_all;            # Store the names of all the genomes, as read from the FastOrtho.end file

read_FastOrtho2();

#print STDERR "Clusters with putative paralogs: ".scalar (keys %clusters_with_paralogs)."\n";
read_write_Fasta();
align_Fasta();

## SUBS #######################################################################

sub read_FastOrtho2 {
    # Read and parse FastOrtho output for ALL clusters
    open(FASTORTHO, "<", $fastortho_result) or die ("$!\n");
        while (<FASTORTHO>) {
        chomp;
        my @theline = split "\t", $_;
        if ($theline[0] =~ /^(ORTHOMCL\d+) /) {
            my $current_cluster = $1;
            #print $current_cluster."\n";
            my @theproteins = split /\s+/, $theline[1];
            foreach my $theentry (@theproteins) {
                if ($theentry =~ /^(\S+)\((\S+)\)/) {
                    # Array reference - store product ID as value of hash with genome name as key
                    push @{$cluster_hash{$current_cluster}{$2}}, $1;
                    # Record all the genomes that occur in these clusters
                    $genomes_all{$2}=1;
                }
            }
            #print "\n";
        }
    }
    close (FASTORTHO);
}

sub read_write_Fasta {
    # For each cluster, write a separate Fasta file, with the genome names as headers 
    my $db = Bio::DB::Fasta->new($fasta_db);
    my @ids = $db->get_all_primary_ids;
    #print STDERR scalar @ids."\t sequences read from Fasta file\n";    # Report number of sequences
    foreach my $thecluster (sort keys %cluster_hash) {
        my $seqio_obj = Bio::SeqIO->new(-file=>">$outpath/$thecluster.fasta",-format=>'fasta');
        for my $thegenome (sort keys %{$cluster_hash{$thecluster}}) {
            foreach my $theprotein (@{$cluster_hash{$thecluster}{$thegenome}}) {
                my $seq = $db->get_Seq_by_id($theprotein);
                my $newseq=Bio::Seq->new(-seq=>$seq->seq(), -id=>$theprotein);
                $seqio_obj->write_seq($newseq);
            }
        }
    }
}

sub align_Fasta {
    # For each cluster's Fasta file, align with Muscle
    foreach my $thecluster (sort keys %cluster_hash) {
        system ("muscle -in $outpath/$thecluster.fasta -out $outpath/$thecluster.muscle.aln 2\> /dev/null");
        my $in=Bio::AlignIO->new(-file=>"$outpath/$thecluster.muscle.aln",-format=>'fasta');
        my $thein=$in->next_aln();
        $alignment_hash{$thecluster}=$thein;
    }
}