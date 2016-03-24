# MITONOTATE

Annotation pipeline for ciliate mitochondrial genomes.

This software package supplements the following publication:

> Seah et al. XXX. ...

Please refer to built-in help message for more information.

```bash
 $ perl mitonotate.pl --help
```

## Description of the pipeline

### Structural annotation

* tRNA genes
  * `tRNAscan-SE` with genetic code 4
* rRNA genes
  * `nhmmer` with custom HMM models for SSU and LSU genes
* Protein-coding genes
  * `Prodigal` in single mode, genetic code 4

### Functional annotation

Protein-coding genes (a.a. sequences) are compared against the following databases:

* Swissprot mitochondrial sequences (curated set supplied with Prokka), by `blastp`
* MitoCOGs database (Kannan et al. 2014), by `blastp`
* MFannot-annotated ciliate mitochondrial genomes, by `blastp`
* Pfam-A domains, by `hmmscan`

Sequence statistics calculated:

* No. of transmembrane domains, by `tmhmm`
* Hydropathicity score ("GRAVY" score), with BioPerl
* Sequence length (a.a.)

Ortholog clustering is performed with FastOrtho:

* Reciprocal all-vs-all `blastp`
* `FastOrtho` implementation of the OrthoMCL algorithm (using `mcl`)
* Protein sequences for each ortholog cluster extracted, aligned with `muscle`
* HMM-HMM search against UniProt20 database with `hhblits`

## Output

Results are combined into a single tab-separated table, sorted by the protein sequence ID. Actual assignment of a product name has to be done manually, especially for putative distant homologs.

Output from intermediate steps are also saved, to allow for troubleshooting.

## Configuring the annotation run

`Mitonotate` requires three configuration files to tell it where to look for input, databases, and which steps in the pipeline to run.

Descriptions of the input files are given in the built-in help message.

Templates for each input file are supplied with the script. It is recommended to build your run by modifying the supplied templates, and not by writing your own configuration files from scratch.

## Software dependencies

`Mitonotate` checks your path for the required dependencies on startup. It addition, it also requires BioPerl to calculate sequence statistics. Please remember to cite the dependencies if you use them!

### Known issues

* MitoCOGs files downloaded from NCBI have DOS-style line endings. When running in Linux/Unix, the files must first be converted with `dos2unix` or a similiar utility.

## Support

We regret that no support for the use of this pipeline or the installation of dependencies can be offered for the foreseeable future. Caveat emptor!