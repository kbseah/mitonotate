#!/usr/bin/perl

=head1 NAME

synplot.pl - Make synteny plots for small genomes

=head1 SYNOPSIS

    perl synplot.pl -f <fasta1>,<fasta2>,<fasta3> \
                    -g <gff1>,<gff2>,<gff3> \
                    -o <output_prefix>

    perl synplot.pl --help

=head1 DESCRIPTION

Make synteny plots for a set of contigs or small genomes. Given a set of Fasta
files, each representing a single genome, and GFF3 feature tables with CDS
features for those Fasta files, perform Blastp between each adjacent pair of
genomes, and make synteny plots.

Requires the accompanying R script synplot.R in the same folder as the perl
script, and command Rscript in path.

=head1 ARGUMENTS

=over 8

=item --fasta|-f <file>,<file>,<file>

List of file names separated by commas; nucleotide fasta files containing
contigs to be compared, in order that they will appear on synteny plot.

=item --gff|-g <file>,<file>,<file>

List of file names separated by commas; GFF files containing predicted CDS to
be compared by Blastp against each other. Must be in same order as the list of
Fasta files given to -f parameter (above).

=item --out|-o <string>

Prefix for output file names. (Default: test)

=item --bidir

Flag: Take bidirectional reciprocal best blast hits only. (Default: Off)

=item --plotonly

Flag: Draw plot only. Requires precomputed intermediate files with same file
name prefix as supplied to -o parameter above.

=item --gencode <integer>

Genetic code for a.a. sequence translation. (Default: 4)

=item --help

This help message.

=item --cds|-c <string>

How to depict CDS direction. Allowed values: color, arrow. No quotation marks.
(Default: color)

=item --color_id <string>

Color corresponding to max value in color scale, used to show percentage ID for
two CDSs connected by a stripe in the synteny plot.

=item --color_cds_f <string>

Color for forward-directed CDSs (if "--cds color" specified), or color of arrow
(if "--cds arrow" specified).

=item --color_cds_r <string>
Color for reverse-directed CDSs (if "--cds color" specified, otherwise ignored)

=back

=head1 OUTPUT

All output files have output prefix as given to --out.

Synteny plot is in PDF format at <output_prefix>.synteny.pdf.

Intermediate files: Blastp output in tabular format (outfmt 6), and tables
containing coordinates of features and connecting polygons for plotting (suffix
.tab).

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

use warnings;
use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::DB::Fasta;
use FindBin qw($Bin $Script);
use Pod::Usage;

my $input_fasta;            # List of Fasta files containing inputs
my $input_gff;              # List of GFF files containing gene predictions
my $gencode=4;              # Genetic code
my $outfix="test";          # Output prefix
my $cdstype="color";        # How to indicate CDS direction
my $colmax="red";           # Color for maximum %ID in color scale
my $colcds1="darkblue";     # Color for CDS 
my $colcds2="darkgreen";
my @valid_cds = qw(color arrow); # List of valid CDS drawing types
my %contig_length_hash;     # Hash of hash of contig lengths
my %total_length_hash;      # Hash of total Fasta file lengths
my %x1_hash; # Don't ask
my %x0_hash;
my $plotonly;               # Flag - do not run Blastp, only call synplot.R to draw plot
my $bidir;                  # Flag - Bidrectional best hits only

if (! @ARGV) {
    pod2usage (-message=>"No input parameters supplied", -exitstatus=>2);
}

GetOptions (
    'fasta|f=s' => \$input_fasta,   # Input should be list of Fasta files, comma-separated
    'gff|g=s' => \$input_gff,       # Input should be list of GFF files, comma-separated
    'plotonly|p'=>\$plotonly,       #
    'out|o=s' =>\$outfix,           # Output prefix
    'gencode=i' =>\$gencode,      # Genetic code
    'cds|c=s' =>\$cdstype,          # How to indicate CDS direction
    'color_id=s' =>\$colmax,        # Color for maximum %id in synteny plot
    'color_cds_f=s' =>\$colcds1,
    'color_cds_r=s' =>\$colcds2,
    'bidir' => \$bidir,             # Bidirectional best hits only
    'help|h' => sub { pod2usage(-exitstatus=>2, verbose=>2); },
    'man|m' => sub { pod2usage(-exitstatus=>0, verbose=>2); },
);

## MAIN #######################################################################

my @input_fasta_list = split /,/, $input_fasta;
my @input_gff_list = split /,/, $input_gff;

# Names for output files
my $output_tab_0 = $outfix.".output_0.tab";
my $output_tab_1 = $outfix.".output_1.tab";  # Table to store summary information for plotting
my $output_tab_2 = $outfix.".output_2.tab";  # Table of feature information for plotting
my $output_tab_3 = $outfix.".output_3.tab";  # Table of polygons for plotting

if (!$plotonly) {
    # Headers for output tables
    open(OUTPUT0, ">", "$output_tab_0") or die ("$!\n");
    print OUTPUT0 join ("\t", qw(genome y label)). "\n";
    close(OUTPUT0);
    open(OUTPUT1, ">", "$output_tab_1") or die ("$!\n");
    print OUTPUT1 join("\t", qw(genome contig length start stop y))."\n";
    close(OUTPUT1);
    open(OUTPUT2, ">", "$output_tab_2") or die ("$!\n");
    print OUTPUT2 join ("\t", qw(genome gene start stop cumulstart cumulstop color))."\n";
    close(OUTPUT2);
    open(OUTPUT3, ">", "$output_tab_3") or die ("$!\n");
    print OUTPUT3 join ("\t", qw (genome1 gene1 genome2 gene2 g1x0 g1x1 g2x1 g2x0 g1x0 g1y0 g1y1 g2y1 g2y0 g1y0 pid))."\n";
    close(OUTPUT3);
    
    # Run functions
    parse_fasta_gff();
    if ($bidir) {
        run_blast_pairs_bidir();
    } else {
        run_blast_pairs();
    }
}

# Call R script to generate plot
system ("Rscript $Bin/synplot.R --args $output_tab_0 $output_tab_1 $output_tab_2 $output_tab_3 $outfix.synteny.pdf $cdstype $colmax");

## SUBROUTINES ################################################################

sub parse_fasta_gff {
    # Some ideas and code from here: https://www.biostars.org/p/46281/
    ## Record contig lengths for each Fasta file
    my %contig_zero_position;
    my %running_total;
    
    for my $i (0 .. scalar(@input_fasta_list)-1) {
        my $the_fasta = $input_fasta_list[$i];
        $the_fasta =~ /(.*)\.fasta/;
        my $label = $1;
        open(OUTPUT0, ">>", "$output_tab_0") or die ("$!\n");;
        print OUTPUT0 $the_fasta."\t".$i."\t".$label."\n";
        close(OUTPUT0);
        $running_total{$the_fasta} = 0;
        ## Predict ORFs with Prodigal if not already supplied in GFF file
        # system ("prodigal -m -c -g 4 -a $the_fasta.prodigal.pep -q -p single -f gff -o $the_fasta.prodigal.gff");  
        my $the_fasta_object = Bio::SeqIO->new(-file => $the_fasta);
        my $seq_object;
        open(OUTPUT1, ">>", "$output_tab_1") or die ("$!\n");
        while ($seq_object = $the_fasta_object->next_seq) {
            $contig_length_hash{$the_fasta}{$seq_object->display_id} = $seq_object->length;
            $total_length_hash{$the_fasta} += $seq_object->length;
            $contig_zero_position{$the_fasta}{$seq_object->display_id} = $running_total{$the_fasta};
            $running_total{$the_fasta} += $seq_object->length;      # update running total
            #Tabulate contig lengths and positions along concatenated plot of genome
            print OUTPUT1 join ("\t",
                                $the_fasta,
                                $seq_object->display_id,
                                $seq_object->length,
                                $contig_zero_position{$the_fasta}{$seq_object->display_id},
                                $running_total{$the_fasta},
                                $i
                               )."\n"; 
        }
        #print $the_fasta."\t".$total_length_hash{$the_fasta}."\n";
        close (OUTPUT1);
    ## Index Fasta files and parse corresponding GFF files
        # Load fasta sequences to memory
        my $db = Bio::DB::Fasta->new($the_fasta);
        #my %CDS;
        # Output file for translated CDS sequences
        my $outfile_pep = Bio::SeqIO->new(-format=>'fasta', -file=> ">$the_fasta.pep"); 
        # Open GFF file
        open(GFF, "<", $input_gff_list[$i]) or die ("$!\n");
        open (OUTPUT2, ">>", "$output_tab_2") or die ("$!\n");
        while (<GFF>) {
            # For each feature
            chomp;
            if (!/^\#/) {
                # If not a comment line
                # Split GFF line into fields
                my @array = split("\t",$_);
                # Split notes field into elements
                my @attrs = split(";",$array[8]);
                $attrs[0] =~ s/ID=//; 
                # Gene name parsed from ID field
                my $gene_name = $attrs[0];            
                my $start;
                my $stop;
                # What type of feature
                my $type = $array[2];
                my $current_contig = $array[0];
                # Extract gene sequence and ID from loaded Fasta file using GFF feature table
                my $gene_seq = $db->seq($array[0],
                                        $array[3],
                                        $array[4]
                                        ); 
                #print $db->seq($array[0],$array[3],$array[4])."\n";
                # Output sequence object
                my $output_gene = Bio::Seq->new(
                    -seq => $gene_seq,
                    -id => $gene_name,
                    -display_id => $gene_name,
                    -alphabet => 'dna',
                );
                if ($array[6] eq '+') {
                    $start = $array[3];
                    $stop = $array[4];
                } elsif ($array[6] eq '-') {
                    # Reverse complement if feature is '-'
                    $output_gene=$output_gene->revcom();
                    $start = $array[4];
                    $stop=$array[3];
                } 
                if ($type eq "CDS") {
                    # If CDS, write translation to file
                    # Translation table 4 (protozoan mitochondrial)
                    my $output_pep = $output_gene->translate(-codontable_id=>$gencode);    
                    $outfile_pep->write_seq($output_pep);
                }
                my $zerostart = $start + $contig_zero_position{$the_fasta}{$current_contig};
                my $zerostop = $stop + $contig_zero_position{$the_fasta}{$current_contig};
                # Define color for CDS depending on transcription direction
                my $cds_color;  
                if ($zerostart < $zerostop) {
                    $cds_color=$colcds1;
                }
                elsif ($zerostart >= $zerostop) {
                    $cds_color=$colcds2;
                }
                print OUTPUT2 join ("\t",
                                    $the_fasta,
                                    $gene_name,
                                    $start,
                                    $stop,
                                    $zerostart,
                                    $zerostop,
                                    $cds_color) . "\n";
                $x0_hash{$the_fasta}{$gene_name} = $zerostart;
                $x1_hash{$the_fasta}{$gene_name} = $zerostop;
            }
        }
        close (OUTPUT2);
        close(GFF);
    }
}

sub run_blast_pairs {
    for my $i (0 .. scalar(@input_fasta_list)-2) {
        my $j = $i + 1;
        my $blastfile1 = "$input_fasta_list[$i].pep";
        my $blastfile2 = "$input_fasta_list[$i+1].pep";
        system ("blastp -subject $blastfile1 -query $blastfile2 -evalue 1e-3 -outfmt 6 -max_target_seqs 1 -out $outfix.blastout.$i.out6");
        open(OUTPUT3, ">>", "$output_tab_3") or die ("$!\n");
        open(HITS, "<", "$outfix.blastout.$i.out6") or die ("$!\n");
        while (<HITS>) {
            # Convert Blast hit results (pairs of genes with best hits to each other) to polygons for drawing synteny diagrams
            chomp;
            my @splitline = split("\t",$_);
            my ($query,$subject,$pid) = ($splitline[0],$splitline[1],$splitline[2]);
            print OUTPUT3 $input_fasta_list[$i]."\t".$query."\t".$input_fasta_list[$i+1]."\t".$subject."\t";
            print OUTPUT3 join("\t",
                               $x0_hash{$input_fasta_list[$i+1]}{$query},
                               $x1_hash{$input_fasta_list[$i+1]}{$query},
                               $x1_hash{$input_fasta_list[$i]}{$subject},
                               $x0_hash{$input_fasta_list[$i]}{$subject},
                               $x0_hash{$input_fasta_list[$i+1]}{$query},
                               $j,
                               $j,
                               $i,
                               $i,
                               $j,
                               $pid
                               )."\n";
        }
        close(HITS);
        close(OUTPUT3);
    }
}

sub run_blast_pairs_bidir {
    # Take bidirectional best hits only
    for my $i (0 .. scalar(@input_fasta_list)-2) {
        my $j = $i + 1;
        my $blastfile1 = "$input_fasta_list[$i].pep";
        my $blastfile2 = "$input_fasta_list[$i+1].pep";
        my %hit_ji;
        system ("blastp -subject $blastfile1 -query $blastfile2 -evalue 1e-3 -outfmt 6 -max_target_seqs 1 -out $outfix.blastout.$i.out6");
        system ("blastp -subject $blastfile2 -query $blastfile1 -evalue 1e-3 -outfmt 6 -max_target_seqs 1 -out $outfix.blastout.$i.r.out6");
        open(HITJI, "<", "$outfix.blastout.$i.r.out6") or die ("$!\n");
        while (<HITJI>) {
            chomp;
            my @splitline = split "\t";
            $hit_ji{$splitline[0]} = $splitline[1];
        }
        close(HITJI);
        open(OUTPUT3, ">>", "$output_tab_3") or die ("$!\n");
        open(HITIJ, "<", "$outfix.blastout.$i.out6") or die ("$!\n");
            while (<HITIJ>) {
                chomp;
                my @splitline = split "\t";
                if (defined $hit_ji{$splitline[1]} && $hit_ji{$splitline[1]} eq $splitline[0]) {
                    # Only print hits which are bidirectional best hits
                    my ($query, $subject, $pid) = ($splitline[0], $splitline[1], $splitline[2]);
                        print OUTPUT3 $input_fasta_list[$i]."\t".$query."\t".$input_fasta_list[$i+1]."\t".$subject."\t";
                        print OUTPUT3 join("\t",
                                           $x0_hash{$input_fasta_list[$i+1]}{$query},
                                           $x1_hash{$input_fasta_list[$i+1]}{$query},
                                           $x1_hash{$input_fasta_list[$i]}{$subject},
                                           $x0_hash{$input_fasta_list[$i]}{$subject},
                                           $x0_hash{$input_fasta_list[$i+1]}{$query},
                                           $j,
                                           $j,
                                           $i,
                                           $i,
                                           $j,
                                           $pid
                                           )."\n";
                }
            }
        close(HITIJ);
        close(OUTPUT3);
    }
}