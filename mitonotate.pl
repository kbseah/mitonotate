#!/usr/bin/env perl

=head1 NAME

mitonotate.pl - Annotation pipeline for ciliate mitochondrial genomes

=head1 SYNOPSIS

    perl mitonotate.pl -o <output dir> \
                       -p <output prefix> \
                       -f <fasta table> \
                       -d <db table> \
                       -r <flags table>

    perl mitonotate.pl --help [Full help message]

    perl mitonotate.pl --man [Manual page]

=head1 DESCRIPTION

Annotation pipeline for ciliate mitochondrial genomes. May also be useful
for other small organellar genomes. Arguments are supplied via config
files (use templates to create your own). Input is a set of Fasta files
containing genomic sequences of the organellar genomes. Mitonotate will
call third-party tools to annotate features (RNA genes, ORFs), and for
functional annotation.

If you use Mitonotate, please acknowledge the author and also cite the
dependencies (full list in documentation).

=head1 ARGUMENTS

=over 8

=item --outdir|-o <path>

Path to output folder where all results files will be written.

=item --prefix|-p <string>

File name prefix for output files (not a path).

=item --fasta|-f <file>

Genomic Fasta files for annotation. Input is a tab-delimited file, with
short alias for the genome in first column, and path to Fasta file in
second column

=item --db|-d <file>

Databases used by Mitonotate. Input is a tab-delimited file, with alias
for database in first column, and path to database in second column. Use
the template supplied with Mitonotate to create your own (aliases are
hard coded in the pipeline).

=item --run|-r <file>

List of annotation tasks to perform. Input is a tab-delimited file, with
task name in first column, and either '1' (perform task) or '0' (skip task)
in second column. Use template supplied to create your own.

=item --gencode|-g <integer>

Genetic code to use for Prodigal ORF prediction (NCBI translation table no.)
Default: 4

=item --cpus|-c <integer>

Number of CPUs for commands that are parallelized.
Default: 8

=item --help|-h

Full help message

=item --man|-m

Manual page

=back

=head1 COPYRIGHT AND LICENSE

Copyright 2016, Brandon Seah (kbseah@mpi-bremen.de)

LICENSE
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

## FRONTMATTER ################################################################

# Modules
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Path qw(make_path);
use File::Spec;
use File::Temp;
use Cwd;
use FindBin qw($Bin $Script);
use threads;
use Log::Message::Simple;
use IPC::Cmd qw(can_run);
use Bio::SeqIO;
use Bio::Tools::SeqStats;

# Redirect MSG output to STDERR
local $Log::Message::Simple::MSG_FH = \*STDERR;

my $version = "0.1";
my $author = "Brandon Seah";

# Input parameters
my $fasta_list_file;    # Table of genome shortnames and fasta files
my $db_config_file;     # Table of paths to databases
my $run_file;           # Table of flags for subroutines
my $outdir;             # Output folder
my $prefix;             # Prefix for output files
my $gencode = 4;        # Genetic code (NCBI translation table no.)
my $CPUs = 8;           # Number CPUs for subroutines where applicable

# Lists
# Dependencies to check
my @depends = qw (perl cat sed makeblastdb blastp prodigal parallel muscle nhmmer hmmscan FastOrtho tRNAscan-SE tmhmm hhblits addss.pl reformat.pl multithread.pl);
# Outputs to parse (keys to %output_files hash)
my @toparse = qw (gravy sprot mfannot pfam tmhmm);
# Flags to check in config file
my @flags_check = qw (prodigal mitossu mitolsu tRNAscan tmhmm concat_pep FastOrtho FastOrtho_cluster pfam makeblastdb reciprocal gravy sprot mfannot mitocogs hhpred read writetbl advertising);
# DBs to check in config file
my @db_check = qw (MFANNOT_BLASTDB SWISSPROT_BLASTDB SWISSPROT_FASTA HHPRED_DB PFAM_A MITOCOGS_BLASTDB MITOCOGS_SEQINFO MITOCOGS_COGINFO MTSSU_MODEL MTLSU_MODEL GCODE_FILE);
# Fields to report (keys to %{pep{...}})
my @toreport = qw (genome startpos endpos dir length tmhmm gravy pfam sprot_name sprot_eval mfannot_name mfannot_eval mitocogs_hit mitocogs_eval mitocogs_cog mitocogs_func ortho);

# Hashes for inputs and files
my %dep_paths;      # Paths to dependencies
my %fasta_list;     # List of genomic fasta files
my %db_config;      # List of databases
my %output_paths;
my %output_files;
my %flags;          # Hash of flags for which tasks to perform
my %ortholist;      # List of ortholog cluster names

# Hashes for results
my %pep;    # Results keyed by ORF name
my %clust;  # Results keyed by cluster name

# Get options or show help
GetOptions ('fasta|f=s'=>\$fasta_list_file,
            'db|d=s'=>\$db_config_file,
            'run|r=s'=>\$run_file,
            'outdir|o=s'=>\$outdir,
            'prefix|p=s'=>\$prefix,
            'gencode|g=i'=>\$gencode,
            'cpus|c=i'=>\$CPUs,
            'help|h'=> sub { pod2usage( -exitstatus => 2, -verbose => 2); },
            'man|m'=> sub { pod2usage ( -exitstatus => 0, -verbose => 2) },
            );
if (!defined $run_file || !defined $fasta_list_file || !defined $db_config_file) {
    pod2usage(-message => "Insufficient options were supplied", -existatus => 2);
}

## MAIN #######################################################################

msg ("Starting Mitonotate $version by $author",1);
msg ("Script location: $Script", 1);

# Read input and config files #############################
%fasta_list = hash_TSV ($fasta_list_file);
%db_config = hash_TSV ($db_config_file);
%flags = hash_TSV ($run_file);

# Check dependencies
check_env();

# Check config files
check_hash_keys_vs_array("Run flags",\@flags_check,\%flags); 
check_hash_keys_vs_array("DB paths",\@db_check,\%db_config);
  # Pass by reference because mixing array and hash

# Make output folder
$outdir = File::Spec->rel2abs ($outdir);
msg("Writing output to folder: $outdir",1);
my $outfix = $outdir."/".$prefix;
make_path ($outdir);

# Run structural annotations ##############################
foreach my $thekey (keys %fasta_list) {
    my $prodout = $outdir."/".$thekey;
    my $thefile = File::Spec->rel2abs($fasta_list{$thekey});
    msg ("Structural annotate genome $thekey",1);
    
    # Predict ORFs
    if ($flags{"prodigal"} == 1) {
        msg("Predict ORFs with Prodigal",1);
        do_prodigal($thefile,
                    $prodout,
                    $thekey
                    );
    } else {
        msg ("Skip predict ORFs with Prodigal", 0);
    }
    $output_files{$thekey}{"prodigal.pep"} = "$prodout.prodigal_single.pep";
    $output_files{$thekey}{"prodigal.gff"} = "$prodout.prodigal_single.gff";
    
    # Predict rRNAs with custom mito models
    if ($flags{"mitossu"} == 1) {
        msg("Predict mtSSU rRNA wtih nhmmer",1);
        do_nhmmer($thefile,
                  $db_config{"MTSSU_MODEL"},
                  "$prodout.nhmmer.mitossu",
                  $thekey,"mitossu"
                  );
    } else {
        msg ("Skip predict mtSSU rRNA with nhmmer",0);
    }
    $output_files{$thekey}{"nhmmer.mitossu"} = "$prodout.nhmmer.mitossu.tblout";
    if ($flags{"mitolsu"} == 1) {
        msg ("Predict mtLSU rRNA with nhmmer",1);
        do_nhmmer($thefile,
                  $db_config{"MTLSU_MODEL"},
                  "$prodout.nhmmer.mitolsu",
                  $thekey,
                  "mitolsu"
                  );
    } else {
        msg ("Skip predict mtLSU rRNA with nhmmer",0);
    }
    $output_files{$thekey}{"nhmmer.mitolsu"} = "$prodout.nhmmer.mitolsu.tblout";
    
    ## Predict tRNAs with tRNAscan-SE # Linear version; now parallelized below
    #if ($flags{"tRNAscan"} == 1) {
    #    msg("Predict tRNAs with tRNAscan-SE",1);
    #    do_tRNAscan($thefile,
    #                $prodout,
    #                $db_config{"GCODE_FILE"},
    #                $thekey
    #                );
    #} else {
    #    msg ("Skip predict tRNAs with tRNAscan-SE",0);
    #}
    #$output_files{$thekey}{"tRNAscan"} = "$prodout.tRNAscan.out";

    # Concatenate predicted ORFs
    if ($flags{"concat_pep"} == 1) {
        msg ("Concatenate prodigal ORFs",1);
        my $pepfile = $output_files{$thekey}{"prodigal.pep"};
        msg ("Running subcommand: cat $pepfile \>\> $outfix.concat.pep", 0);
        system ("cat $pepfile \>\> $outfix.concat.pep");
    } else {
        msg ("Skip concatenate prodigal ORFs",0);
    }
    $output_files{"concat.pep"} = "$outfix.concat.pep";
}

# Paralellize tRNAscan-SE #################################
# Adapted from http://chicken.genouest.org/perl/multi-threading-with-perl/

if ($flags{"tRNAscan"} == 1) {
    msg ("Start $CPUs threads for tRNAscan", 1);
    my @jobs_arr = (keys %fasta_list);  # List of jobs to be run
    my $total_jobs = scalar @jobs_arr;  # Count total no of jobs
    my @running;                        # List of running jobs
    my @Threads;                        # List of all threads
    my %jobs2name;                      # Hash for job names
    
    while (scalar @Threads < $total_jobs) { # While there are still jobs to run
        @running = threads->list(threads::running); # List running jobs
        if (scalar @running < $CPUs && scalar @jobs_arr > 0) {
            # If slot available and still jobs left to do
            # Get next job params
            my $thekey = pop @jobs_arr; 
            my $prodout = $outdir."/".$thekey; 
            my $thefile = File::Spec->rel2abs($fasta_list{$thekey});
            msg ("Predict tRNAs for $thekey with tRNAscan-SE",1);
            # Start thread
            my $thread = threads->new(\&do_tRNAscan,
                                      $thefile,
                                      $prodout,
                                      $db_config{"GCODE_FILE"},
                                      $thekey);
            # Hash thread to job name
            $jobs2name{$thread} = $thekey;
            $output_files{$thekey}{"tRNAscan"} = "$prodout.tRNAscan.out";
            # Add thread to list of started threads
            push @Threads, $thread;
        }
        @running = threads->list(threads::running); # Update list of running jobs
        foreach my $thr (@Threads) { # For all started threads
            if ($thr->is_joinable()) { # Check which are joinable
                $thr->join(); # If so, join it and report
                msg ("tRNAscan job for $jobs2name{$thr} complete",1);
            }
        }
        @running = threads->list(threads::running); # Update list of running jobs
        sleep(1); # Delay to prevent incessant looping
    }
    while (scalar @running > 0) { # Sweep up remaining started jobs
        foreach my $thr (@Threads) {
            if ($thr->is_joinable()) {
                $thr->join();
            }
        }
        @running = threads->list(threads::running);
        sleep(1); # Delay to prevent incessant looping
    }
} else {
    msg ("Skip predict tRNAs with tRNAscan-SE",0);
}

# Prepare Blastdb of concatenated ORFs ####################
if ($flags{"makeblastdb"} == 1) {
    msg ("Make Blastdb of concatenated ORFs", 1);
    do_makeblastdb($output_files{"concat.pep"},"prot");
} else {
    msg ("Skip make blastdb of concatenated ORFs", 0);
}
$output_files{"concat.pep.blastdb"} = $output_files{"concat.pep"};

# Run functional annotations ##############################
msg ("Run functional annotations",1);

# Predict transmembrane helices with TMHMM
if ($flags{"tmhmm"} == 1) {
    msg ("Predict TM helices with tmhmm",1);
    do_tmhmm($output_files{"concat.pep"},$outfix);
} else {
    msg ("Skip predict TM helices with tmhmm",0);
}
$output_files{"tmhmm"} = "$outfix.tmhmm.out";

# Perform reciprocal Blastp for clustering
if ($flags{"reciprocal"} == 1) {
    msg ("All-vs-all reciprocal blastp", 1);
    do_parallel_blastp($output_files{"concat.pep"},
                       $output_files{"concat.pep.blastdb"},
                       $outfix,
                       "reciprocal",
                       500);
} else {
    msg ("Skip all-vs-all reciprocal blastp",0);
}
$output_files{"concat.pep.reciprocal.blastp"} = "$outfix.reciprocal.blastp.out6";

# FastOrtho clustering of orthologs from reciprocal Blastp results
if ($flags{"FastOrtho"} == 1) {
    msg ("FastOrtho clustering from reciprocal Blastp results", 1);
    make_FastOrtho_options($output_files{"concat.pep.reciprocal.blastp"},
                           "$outfix.FastOrtho.options");
    do_FastOrtho("$outfix.FastOrtho.options");
} else {
    msg ("Skip FastOrtho clustering from reciprocal Blastp results", 0);
}
$output_files{"concat.pep.FastOrtho"} = "$outfix.FastOrtho.end";

# Create Fasta files and alignments for FastOrtho clusters
if ($flags{"FastOrtho_cluster"} == 1) {
    msg ("Align FastOrtho clusters with Muscle",1);
    run_prog ("perl",
              "$Bin/FastOrtho2fasta.pl",
              "-e",         # FastOrtho .end file
              $output_files{"concat.pep.FastOrtho"},
              "-f",         # Concatenated aa sequences
              $output_files{"concat.pep"},
              "-o $outdir"  # Output folder
              );
} else {
    msg ("Skip align FastOrtho clusters with Muscle",0);
}

# Calculate GRAVY scores for peptide sequences
if ($flags{"gravy"} == 1) {
    msg ("Calculate GRAVY scores", 1);
    do_gravy($output_files{"concat.pep"},"$outfix.gravy.out");
} else {
    msg ("Skip calculate GRAVY scores",0);
}
$output_files{"gravy"} = "$outfix.gravy.out";

# Blastp vs Swissprot mitochondrial
if ($flags{"sprot"} == 1) {
    msg ("Blastp vs Sprot database", 1);
    do_parallel_blastp($output_files{"concat.pep"},
                       $db_config{"SWISSPROT_BLASTDB"},
                       $outfix,
                       "vs_swissprot",
                       1);
} else {
    msg ("Skip blastp vs Sprot database", 0);
}
$output_files{"sprot"} = "$outfix.vs_swissprot.blastp.out6";

# Blastp vs MFannot results 
if ($flags{"mfannot"} == 1) {
    msg ("Blastp vs MFannot results", 1);
    do_parallel_blastp($output_files{"concat.pep"},
                       $db_config{"MFANNOT_BLASTDB"},
                       $outfix,
                       "vs_mfannot",
                       1);
} else {
    msg ("Skip blastp vs MFannot results", 0);
}
$output_files{"mfannot"} = "$outfix.vs_mfannot.blastp.out6";

# Blastp vs MitoCOGs DB
if ($flags{"mitocogs"} == 1) {
    msg ("Blastp vs MitoCOGs database", 1);
    do_parallel_blastp($output_files{"concat.pep"},
                       $db_config{"MITOCOGS_BLASTDB"},
                       $outfix,
                       "vs_mitocogs",
                       1);
} else {
    msg ("Skip blastp vs MitoCOGs database",0);
}
$output_files{"mitocogs"} = "$outfix.vs_mitocogs.blastp.out6";

# HMMer vs Pfam-A
if ($flags{"pfam"} == 1) {
    msg ("HMMer vs Pfam-A", 1);
    do_hmmscan($output_files{"concat.pep"},
               $db_config{"PFAM_A"},
               $outfix,
               "vs_pfam_A");
} else {
    msg ("Skip HMMer vs Pfam-A",0);
}
$output_files{"pfam"} = "$outfix.vs_pfam_A.tblout";

# Prepare alignments for HHBlits
# Get names of Ortholog clusters
open(ORTHOEND, "<", $output_files{"concat.pep.FastOrtho"})
    or error ("Cannot open FastOrtho cluster output. $!", 1);
while (<ORTHOEND>) {
    chomp;
    my @splitline = split " ";
    $ortholist{$splitline[0]} = 1;
}
close(ORTHOEND);
# Reformat Fasta alignments and perform HHblits
if ($flags{"hhpred"} == 1) {
    msg ("Reformat ortholog clusters for HHBlits",1);
    foreach my $theortho (keys %ortholist) {
        run_prog("reformat.pl",
                 "fas",
                 "a3m",
                 "$outdir/$theortho.muscle.aln",
                 "$outdir/$theortho.muscle.a3m",
                 "2\> /dev/null");
        run_prog("addss.pl",
                 "$outdir/$theortho.muscle.a3m",
                 "2\> /dev/null");
    }
    do_hhblits($outdir,$db_config{"HHPRED_DB"});
} else {
    msg ("Skip reformat ortholog clusters for HHBlits",0);
}

# Read and hash results ###################################

if ($flags{"read"} == 1) {
    msg ("Read annotation results into memory", 1);
    # Hash prodigal protein names from pep files
    foreach my $thegenome (keys %fasta_list) {
        read_prodigal_pep($output_files{$thegenome}{"prodigal.pep"},
                          $thegenome);
    }
    
    # Hash pfam tblout results
    read_pfam_tblout($output_files{"pfam"});
    
    # Hash FastOrtho cluster memmbership
    read_FastOrtho($output_files{"concat.pep.FastOrtho"});
    
    # Hash TMHMM results
    read_tmhmm($output_files{"tmhmm"});
    
    # Hash GRAVY scores and peptide lengths
    read_gravy($output_files{"gravy"});
    
    # Hash Sprot blastp hits
    read_sprot_blastp($db_config{"SWISSPROT_FASTA"},
                      $output_files{"sprot"});
    
    # Hash MFannot blastp hits
    read_mfannot_blastp($output_files{"mfannot"});
    
    # Hash MitoCOGs blastp hits
    read_mitocogs_blastp($output_files{"mitocogs"},
                         $db_config{"MITOCOGS_SEQINFO"},
                         $db_config{"MITOCOGS_COGINFO"});
    
    # Hash HHpred results for FastOrtho clusters
    foreach my $theclust (keys %ortholist) {
        # Get array ref for top hhpred hit
        $clust{$theclust}{"hhpred"} = read_hhpred_hhr("$outdir/$theclust.muscle.hhr");
    }
    
    } else {
        msg ("Skip read annotation results into memory", 0);
}

# Write report table ######################################

if ($flags {"writetbl"} == 1) {
    msg ("Write annotation to file $outfix.mitonotate.out",1 );
    
    open(FINALOUT, "> $outfix.mitonotate.out")
      or error ("Cannot open $outfix.mitonotate.out for writing $!",1);
    # Header row
    my @clust_heads = qw (hhpred_id hhpred_name hhpred_prob hhpred_eval hhpred_score);
      # Names for cluster annotation fields
    print FINALOUT "\#";
    print FINALOUT join "\t", ("id",@toreport,@clust_heads);
    print FINALOUT "\n";
    
    # Populate table
    foreach my $pepid (sort {$a cmp $b} keys %pep) {
        my @line;
        push @line, $pepid;
        # Peptide annotations
        foreach my $field (@toreport) {
            if (!defined $pep{$pepid}{$field}) {
                push @line, "NA";
            } else {
                if ($field eq "pfam") { # For fields that are arrays
                    my $fieldjoin = join ",", @{$pep{$pepid}{$field}};
                    push @line, $fieldjoin;
                } else {
                    push @line, $pep{$pepid}{$field};
                }
            }
        }
        # Cluster annotations
        if (defined $pep{$pepid}{"ortho"}) {
            my $theclust = $pep{$pepid}{"ortho"};
            push @line, @{$clust{$theclust}{"hhpred"}};
        } else {
            push @line, qw(NA NA NA NA NA); # hacky
        }
        
        print FINALOUT join "\t", @line;
        print FINALOUT "\n";
    }
    close (FINALOUT);
    } else {
        msg ("Skip write annotation to file",0);
}
# Clean up and write log file #############################

msg ("Pipeline complete, results in $outdir",1);
msg ("Write log file to $outfix.log", 1);

# Dump message stack into log file
my $msgs  = Log::Message::Simple->stack_as_string;
open(LOGOUT, ">>", "$outfix.log")
    or error ("Cannot open log file for writing. $!", 1);
print LOGOUT $msgs;
print LOGOUT "\n\n";
close(LOGOUT);

msg ("$author thanks you for using Mitonotate $version", 1);
if ($flags{"advertising"} == 1) {
    # Shameless advertising
    msg ("...",1);
    msg ("While I have your attention... ", 1);
    msg ("...", 1);
    msg ("perhaps you might be interested in another of our fine products:", 1);
    msg ("gbtools - an R package for metagenome binning visualization", 1);
}

## SUBS #######################################################################

sub hash_TSV {
    my $infile = shift @_;
    my %outhash;
    open(IN, "<", $infile)
        or error ("Cannot read input file. $!", 1);
    while (<IN>) {
        chomp;
        if (m/^#/) { # Skip comment lines
            next();
        }
        my @splitline = split /\s+/; # Split by any whitespace
        if (scalar @splitline >= 2) {
            # Hash column 2 by column 1, if table
            $outhash{$splitline[0]} = $splitline[1];
            # Ignore column 3+
        } elsif (scalar @splitline == 1) {
            # Dummy value if simple list
            $outhash{$splitline[0]} = "1";
        } else {
            next();
        }
    }
    close(IN);
    return %outhash;
}

sub check_hash_keys_vs_array {
    my ($name, $arrayinref, $hashinref) = @_;
    my @arrayin = @$arrayinref;
    my %hashin = %$hashinref;
    my %arrhash;
    foreach my $thekey (@arrayin) {
        $arrhash{$thekey} = 1;
    }
    foreach my $thekey2 (keys %hashin) {
        if (!defined $arrhash{$thekey2}) {
            error("Element $thekey2 in $name not defined",1);
        }
    }
}

sub do_prodigal {
    my ($infile, $outpath, $genome) = @_;
    run_prog("prodigal",
             "-q",
             "-g $gencode",
             "-c",
             "-m",
             "-f gff",
             "-p single",
             "-a $outpath.prodigal_single.pep",
             "-o $outpath.prodigal_single.gff",
             "-d $outpath.prodigal_single.orf",
             "-i $infile",
             "2\> /dev/null");
}

sub do_nhmmer {
    my ($infile, $hmm, $outpath, $genome, $modelname) = @_;
    run_prog("nhmmer",
             "--tblout $outpath.tblout",
             "$hmm",
             "$infile",
             "\> /dev/null");
}

sub do_tRNAscan {
    my ($infile,$outpath,$gcodefile,$genome) = @_;
    run_prog("tRNAscan-SE",
             "-g $gcodefile",
             "-O",
             "-o $outpath.tRNAscan.out",
             "$infile",
             "2\> /dev/null");
}

sub do_tmhmm {
    my ($infile, $outpath) = @_;
    run_prog ("cat",
              $infile,
              "\|",
              "sed",
              "'s/\\*//g'",
              "\|",
              "tmhmm",
              "\>",
              "$outpath.tmhmm.out");
}

sub do_gravy {
    my ($infile, $outfile) = @_;
    open(GRAVYOUT, "> $outfile")
        or error ("Cannot open file for writing Gravy results. $!", 1);
    my $seqio = Bio::SeqIO->new(-file=>$infile,
                                -format=>"fasta");
    my $seq;
    while ($seq = $seqio->next_seq()) {
        my $newseq;
        my $theseq = $seq->seq();
        if ($theseq =~ m/[\*]/) {
            # Strip STOP codon character
            $theseq =~ s/\*//; 
        }
        $seq->seq($theseq);
        # Col 1 - Sequence name
        print GRAVYOUT $seq->display_id();
        print GRAVYOUT "\t";
        # Col 2 - GRAVY hydropathicity score
        if ($theseq =~ m/X/) {
            print GRAVYOUT "NA";
        } else {
            # Skip if includes ambiguous amino acid
            print GRAVYOUT Bio::Tools::SeqStats->hydropathicity($seq);
        }
        print GRAVYOUT "\t";
        # Col 3 - Sequence length
        print GRAVYOUT $seq->length(); 
        print GRAVYOUT "\n";
    }
    close (GRAVYOUT);
}

sub do_makeblastdb {
    my ($infile, $dbtype) = @_; 
    run_prog ("makeblastdb",
              "-dbtype $dbtype",
              "-in $infile",
              "-out $infile");
}

sub do_parallel_blastp {
    my ($infile, $db, $outpath, $jobname, $num_targets) = @_;
    run_prog ("cat",
              "$infile",
              "\|",
              "parallel",
              "--no-notice --gnu -j $CPUs --recstart '\>' -N 20 --pipe",
              "blastp",
              "-query -",
              "-db $db",
              "-max_target_seqs $num_targets",
              "-outfmt 6",
              "\>",
              "$outpath.$jobname.blastp.out6");
}

sub do_hmmscan {
    my ($infile, $db, $outpath, $jobname) = @_;
    run_prog ("hmmscan",
              "--cpu $CPUs",
              "--tblout $outpath.$jobname.tblout",
              "-o $outpath.$jobname.hmmscan.out",
              "$db",
              "$infile");
}

sub make_FastOrtho_options {
    my ($blast_tblout,$optfilename) = @_;
    # Make "options" configuration file for FastOrtho to run
    my @options = ("--working_directory $outdir",
                   "--project_name $prefix.FastOrtho",
                    "--blast_file $blast_tblout",
                    "--query_start_index 6",
                    "--query_end_index 7",
                    "--subject_start_index 8",
                    "--subject_end_index 9",
                    "--alignment_length_index 3",
                    "--query_index 0",
                    "--subject_index 1",
                    "--e_value_index 10",
                    "--percent_idenity_index 2",
                    "--use_tab_split",
                    "--pv_cutoff 1e-5",
                    "--pi_cutoff 0.000000",
                    "--pmatch_cutoff 0.000000",
                    "--maximum_weight 316.000000",
                    "--mcl_path mcl",
                    "--inflation 1.5",
                    "--result_file $outfix.FastOrtho.end");
    foreach my $thekey (sort {$a cmp $b} keys %fasta_list){
        my $thefile = $output_files{$thekey}{"prodigal.pep"};
        push @options, "--single_genome_fasta $thefile";
    }
    open(OPTOUT, ">", $optfilename)
        or error ("Cannot open FastOrtho configuration file for writing $!", 1);
    foreach my $optline (@options){
        print OPTOUT $optline."\n";
    }
    close(OPTOUT);
    msg ("FastOrtho options file written to $optfilename",0);
}

sub do_FastOrtho {
    my ($optfilename) = @_;
    run_prog ("FastOrtho",
              "--option_file $optfilename",
              "2\> /dev/null");
}

sub do_hhblits {
    my ($outpath, $db) = @_;
    run_prog ("multithread.pl",
              "'$outpath/\*.a3m'",
              "'hhblits",
              "-i \$file",
              "-d $db",
              "-o \$name.hhr",
              "-n 2'",
              "-cpu $CPUs",
              "2\> /dev/null");
}

sub read_prodigal_pep {
    my ($pepfile, $genome) = @_;
    open(PEP, "<", $pepfile)
      or error ("Cannot open $pepfile for reading $!", 1);
    # Get accession numbers of all the protein seq to be annotated
    while (<PEP>) {
        chomp;
        if ($_ =~ /^>(\w+) # (\d+) # (\d+) # (\S+) #/) {
            $pep{$1}{"genome"} = $genome;
            $pep{$1}{"startpos"} = $2;
            $pep{$1}{"endpos"} = $3;
            $pep{$1}{"dir"} = $4;
        }
    }
    close(PEP);
}

sub read_pfam_tblout {
    my ($infile) = @_;
    open(PFAMRES, "<", $infile)
      or error ("Cannot open $infile for reading $!", 1);
    while (<PFAMRES>) {
        chomp;
        if ($_ =~ /^#/) { next; }   # Skip comment lines
        my @splitline = split /\s+/, $_;
        my $domain = $splitline[0];
        my $query = $splitline[2];
        my $inc = $splitline[17];
        if ($inc > 0) {
            push @{$pep{$query}{"pfam"}}, $domain;
        }
    }
    close (PFAMRES);
}

sub read_FastOrtho {
    my ($infile) = @_;
    open(FASTORTHO, "<", $infile)
      or error ("Cannot open $infile for reading: $!", 1);
    while (<FASTORTHO>) {
        chomp;
        if ($_ =~ /(ORTHOMCL\d+) \(.+\):\t(.*)/) {
            my $ortho_id = $1;
            my @members = split / /, $2;
            foreach my $themember (@members) {
                my ($theseq) = $themember =~ /(\S+)\(\S+\)/;
                $pep{$theseq}{"ortho"} = $ortho_id unless (!defined $theseq);
            }
        }
    }
    close(FASTORTHO);
}

sub read_tmhmm {
    my ($infile) = @_;
    open(TMHMM, "<", $infile)
      or error ("Cannot open $infile for reading $!", 1);
    while (<TMHMM>) {
        chomp;
        if ($_ =~ /# (\S+) Number of predicted TMHs:\s+(\d+)/) {
            #print $1."\t".$2."\n";
            $pep{$1}{"tmhmm"} = $2;
    }
}
close(TMHMM);
}

sub read_gravy {
    my ($infile) = @_;
    open(GRAVYIN, "<", $infile)
      or error ("Cannot open $infile for reading $!", 1);
    while (<GRAVYIN>) {
        chomp;
        my @splitline = split "\t";
        $pep{$splitline[0]}{"gravy"} = $splitline[1];
        $pep{$splitline[0]}{"length"} = $splitline[2];
    }
    
    close (GRAVYIN);
}
sub read_hhpred_hhr {
    my ($infile) = @_;
    my $hitnum=0;
    my %hitdex;
    my @fields = qw(ID name prob eval score);
    open(HHRIN, "<", $infile)
      or error ("Cannot open $infile for reading $!",1);
    while (<HHRIN>) {
        chomp;
        if ($_ =~ m/^No (\d+)/) {
            # New record entry
            $hitnum = $1; # Update record number
        } elsif ($_ =~ m/^>(\S+) (.+?) OX=\d+/) {
            $hitdex{$hitnum}{"ID"} = $1;
            $hitdex{$hitnum}{"name"} = $2;
        } elsif ($_ =~ m/Probab=(\S+?)\s+E-value=(\S+?)\s+Score=(\S+?)/) {
            $hitdex{$hitnum}{"prob"} = $1;
            $hitdex{$hitnum}{"eval"} = $2;
            $hitdex{$hitnum}{"score"} = $3;
        }
    }
    close(HHRIN);
    #foreach my $thekey (keys %hitdex) {
    #    print $thekey;
    #    foreach my $thefield (@fields) {
    #        print "\t";
    #        print $hitdex{$thekey}{$thefield};
    #    }
    #    print "\n";
    #}
    my @output;
    foreach my $thefield (@fields) {
        push @output, $hitdex{"1"}{$thefield};
    }
    return(\@output);
}

sub read_sprot_blastp {
    my ($dbfile, $blastout) = @_;
    my %sprot_prod;
    # Get accession numbers and correpsonding product names for Sprot database
    open (SPROTDB, "<", $dbfile) or error ("Cannot open $dbfile for reading $!", 1);
    while (<SPROTDB>) {
        chomp;
        #if ($_ =~ /^>(\w+) .*~~~(.+)/) {
        if ($_ =~ /^>(\w+) .*~~~(.*)~~~(.*)/) {
        #if ($_ =~ /^>(\w+) /) {
            #$sprot_gen_hash{$1} = $2;
            $sprot_prod{$1} = $3;
        }
    }
    close (SPROTDB);
    # Read Blastp outfmt 6 output
    open(SPBLAST, "<", $blastout) or die ("$!\n");
    while (<SPBLAST>) {
        chomp;
        my @splitline = split "\t";
        my ($thepep, $sprotname, $sproteval) = ($splitline[0],$splitline[1],$splitline[10]);
        if (!defined $pep{$thepep}{"sprot_name"}) {
            $pep{$thepep}{"sprot_name"} = $sprot_prod{$sprotname};
            $pep{$thepep}{"sprot_eval"} = $sproteval;
        } else {
            if ($pep{$thepep}{"sprot_name"} ne $sprot_prod{$sprotname}) {
                $pep{$thepep}{"sprot_name"} = "XXX";
            }
        }
    }
    close(SPBLAST);
}

sub read_mfannot_blastp {
    my ($blastout) = @_;
    open(MFBLAST, "<", $blastout)
      or error ("Cannot open file $blastout for reading $!", 1);
    # Get consensus MFannot annotation result
    while (<MFBLAST>) {
        chomp;
        my @splitline = split "\t", $_;
        my @splithit = split ";", $splitline[1];
        my @thehit = split /_/, $splithit[0];
        my $theprod;
        if ($thehit[1] eq "rearrange" || $thehit[1] eq "mito" || $thehit[1] eq "015981") { # hacky
            $theprod = $thehit[2];
        } else { $theprod=$thehit[1]; }
        my ($query, $eval) = ($splitline[0],$splitline[10]);
        if (!defined $pep{$query}{"mfannot_name"}) {
            if ($theprod =~ m/orf/) {
                # If putative hypothetical, give name 'orf'
                $pep{$query}{"mfannot_name"} = "orf";
            } else {
                # Else give putative product name
                $pep{$query}{"mfannot_name"} = $theprod;
            }
            # Record Blast hit evalue
            $pep{$query}{"mfannot_eval"} = $eval;
        } else {
            if ($pep{$query}{"mfannot_name"} ne $theprod) {
                # If conflicting product name exists, mark 'XXX'
                $pep{$query}{"mfannot_name"} = "XXX" unless $theprod =~ m/orf/;
            }
            
        }
    }
    close(MFBLAST);
}

sub read_mitocogs_blastp {
    my ($blastout, $seqinfo, $coginfo) = @_;
    my %seq2cog;
    my %cog2func;
    open (SEQINFO, "<", $seqinfo)
      or error ("Cannot open file $seqinfo for reading $!", 1);
    while (<SEQINFO>) {
        chomp;
        my @splitline = split ",";
        $seq2cog{$splitline[0]} = $splitline[3];
    }
    close (SEQINFO);
    open(COGINFO, "<", $coginfo)
      or error ("Cannot open file $coginfo for reading $!", 1);
    while (<COGINFO>) {
        chomp;
        my @splitline = split ",";
        $cog2func{$splitline[0]} = $splitline[1];
    }
    close(COGINFO);
    open (BLASTIN, "<", $blastout)
      or error ("Cannot open file $blastout for reading $!", 1);
    while (<BLASTIN>) {
        chomp;
        my @splitline = split /\s+/;
        my @splitsubj = split /\|/, $splitline[1];
        $pep{$splitline[0]}{"mitocogs_hit"} = $splitsubj[1];
        $pep{$splitline[0]}{"mitocogs_eval"} = $splitline[10];
        $pep{$splitline[0]}{"mitocogs_cog"} = $seq2cog{$splitsubj[1]};
        $pep{$splitline[0]}{"mitocogs_func"} = $cog2func{$seq2cog{$splitsubj[1]}};
          # Inception...
    }
    close (BLASTIN);
}

sub run_prog {
    # Adapted from phyloFlash code by Elmar Pruesse
    my $prog = shift @_;
    my @args = @_;
    my $cmd = join " ", ($prog, @args);
    msg("Run subcommand: $cmd", 0);
    system($cmd) == 0
        or error("Tool execution failed! Error was '$!' and return code '$?'",1);
}

sub test_output {
    open(OUT, ">", "$outdir/$prefix.log")
        or error ("Cannot open test log for writing. $!", 1);
    foreach my $thekey (sort {$a cmp $b} keys %fasta_list) {
        print OUT $thekey."\t".$fasta_list{$thekey}."\n";
    }
    close(OUT);
}

sub check_env {
    # Adapted from phyloFlash code by Elmar Pruesse
    foreach my $dep (@depends) {
        if ($dep_paths{$dep} = can_run($dep)) {
            msg("Dependency $dep found at $dep_paths{$dep}",1);
        } else {
            error("Dependency $dep not found",1);
        }
    }
}