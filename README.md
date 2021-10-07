# rmsFinder

A tool for finding Type II restriction-modification systems (RMS) in bacterial genomes.

Given a fasta file of protein sequences rmsFinder will search these proteins against the REBASE database with the aim of identifying any Type II RMS enzymes that are present. REBASE is a comprehensive and curated database of RMS that is regularly updated. It contains different categories of enzyme sequences, ranging from those with experimental support for their restriction site ('gold') to those that have only been predicted bioinformatically based on similarity to known enzymes ('putative').

REBASE citation:

Roberts, R.J., Vincze, T., Posfai, J., Macelis, D.
REBASE-a database for DNA restriction and modification: enzymes, genes and genomes.
*Nucleic Acids Research* 43: D298-D299 (2015).
doi: 10.1093/nar/gku1046

Type II RMS typically consist of a methytransferase (MTase) and a restriction enzyme (REase) with the same target specificity, usually located next to each other in the genome.

Other RMS types are more complex - they may be added one day.

## Usage

Before running rmsFinder for the first time you will need to download the protein database files from REBASE.

```
python updateDB.py
```

This should take <2 minutes and requires ~100 MB of storage in the `data` directory after it has run. REBASE is frequently updated with new proteins, so you can use this command to update the database files periodically. Bear in mind your results could change if you update the database between analyses.

Normal usage of rmsFinder - to predict cognate Type II systems in a genome - looks like

```
python rmsFinder.py --genbank {YOUR_GENBANK}.gbk --mode MT,RE --output {YOUR_OUTPUT} --db all
```

The `--db` option allows you to search against three categories of REBASE protein sequences:
* `gold` - the highest standard . See [http://rebase.neb.com/cgi-bin/rebgoldlist]
* `nonputative` - only those sequences which are not putative
* `all` - includes putative sequences predicted by REBASE based on similarity to existing sequences, but for which no experimental validation is known. In many cases if you are searching an NCBI genome, you will find a 100% match to a putative prediction for a protein because REBASE have already identified it against the existing proteins.

## Workflow

Briefly, rmsFinder takes as input a genbank file. It searches for putative Type II MTases using HMMs produced and curated by Pedro Oliveira, available at [https://github.com/oliveira-lab/RMS/tree/master/RM_HMMs]. It then searches for putative Type II REases using blastp (e<0.001) against REBASE sequences (because Type II REases diverge rapidly, producing poor MSAs that are not good for building HMMs).

rmsFinder then matches these putative MTases or REases to existing REBASE proteins. If the protein similarity is greater than some user-defined thresholds (defaults: 55% for MTases, 50% for REases), the recognition sequence (RS) of the putative protein is predicted to be the same as the matched protein. Finally, if there exists within the genome a Type II MTase and REase recognising the same RS within 4 or fewer genes of each other, rmsFinder predicts that a Type II RMS recognising that RS is present in the genome. Output files include a list of all putative MTases, REases, and RMS.  

This workflow is inspired by the methodology in Oliveira, Touchon and Rocha. *PNAS* (2016) doi: 10.1073/pnas.1603257113.
