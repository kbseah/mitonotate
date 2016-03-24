#!/bin/bash

for i in *.fasta
do
prodigal -m -c -g 4 -q -p single -f gff -o $i.gff -i $i
done