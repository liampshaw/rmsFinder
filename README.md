# rmsFinder

A tool for finding Type II restriction-modification systems (RMS) in bacterial genomes and predicting their target sequences.

Given a fasta file of protein sequences, `rmsFinder` will search these proteins against the REBASE database to identify Type II RMS enzymes that are present, then try to predict their target sequences. This can be done at different stringencies of evidence.

Type II RMS usually consist of a methytransferase (MTase) and a restriction enzyme (REase) with the same target specificity, usually located next to each other in the genome. Homologs can have relatively low similarities (~50\%) and still recognise the same target sequence.

Other RMS types are more complex - they may be added one day.

## Workflow

The workflow is inspired by the methodology in the following paper:

P. H. Oliveira, M. Touchon, E. P. C. Rocha  
Regulation of genetic flux between bacteria by restriction–modification systems  
*PNAS* 113 (20) 5658-5663 (2016)  
doi: [10.1073/pnas.1603257113](https://doi.org/10.1073/pnas.1603257113)

Briefly, `rmsFinder` takes a file of protein sequences as input. It searches for putative Type II MTases using HMMs. One can specify either the HMMs produced and curated by Pedro Oliveira, available [here](https://github.com/oliveira-lab/RMS/tree/master/RM_HMMs) (`hmm oliveira`) or those produced by Tesson et al. for [DefenseFinder](https://github.com/mdmparis/defense-finder-models/tree/master/profiles) (`--hmm tesson`: using files for Type II [MTases](https://github.com/mdmparis/defense-finder-models/blob/master/profiles/RM_Type_II__Type_II_MTases.hmm), [REases](https://github.com/mdmparis/defense-finder-models/blob/master/profiles/RM_Type_II__Type_II_REases.hmm), and [IIG](https://raw.githubusercontent.com/mdmparis/defense-finder-models/master/profiles/RM_Type_IIG__Type_IIG.hmm); these HMMs are provided along with rmsFinder (in `data` folder) but all credit for them goes to those athors). It then searches for putative Type II REases using blastp (e<0.001) against REBASE sequences (because Type II REases diverge rapidly, producing poor MSAs that are not good for building HMMs).

`rmsFinder` then matches these putative MTases or REases to existing REBASE proteins. If the protein identity is greater than some user-defined thresholds (defaults: 55% for MTases, 50% for REases - see below for a note on thresholds), the recognition sequence (RS) of the putative protein is predicted to be the same as the matched protein. Finally, if there exists within the genome a Type II MTase and REase recognising the same RS within 4 or fewer genes of each other, `rmsFinder` predicts that a Type II RMS recognising that RS is present in the genome. Output files include a list of all putative MTases, REases, and RMS.  

**A note on thresholds.** Oliveira et al. (2016) defined their thresholds for protein-protein matches on an empirical basis. There was not a lot of data to inform an exact threshold, but I believe that their choice of 55% for MTases and 50% for REases was in terms of pairwise identity calculated from blast in a local alignment. In contrast, `rmsFinder` uses a global alignment. In my experience these values are typically similar. But be aware of this when using `rmsFinder`, particularly if using it to analyse proteins with extra domains (or spuriously long N-terminus) because this may affect the results of the global alignment. That is, two proteins might be considered as not equivalent under the global identity threshold when they might meet the Oliveira et al. (2016) criteria under a local alignment. On the other hand, using a global alignment is exact rather than heuristic and is probably better at excluding pseudogenes that are still sufficiently large to be annotated. A future update may have the option to use either global identity or local similarity (which is already calculated because blast is used to select the initial matches that are globally aligned). 

## Running rmsFinder

Normal usage looks something like

```
python rmsFinder.py genome.gbk --genbank --output results_genome
```

More detailed arguments can be specified - see below.

### Setup

First, install dependencies using conda:

```
conda create env -f environment.yml
```

Before running `rmsFinder` for the first time you will need to download the protein database files from REBASE.

```
python updateDB.py
```

This should take several minutes and requires ~100 MB of storage in the `data` directory after it has run. REBASE is frequently updated with new proteins, so you can use this command to update the database files periodically. Bear in mind results could change if you update the database between analyses!

### Usage

```
usage: rmsFinder [-h] --output OUTPUT (--genbank | --fasta | --panacotafasta) [--db DB]
                 [--mode MODE] [--dontcollapse] [--hmm HMM]
                 input [input ...]

Predict presence of Type II restriction-modification systems in a genome.

positional arguments:
  input            Input protein file (genbank or fasta).

optional arguments:
  -h, --help       show this help message and exit
  --genbank        Input file is genbank format
  --fasta          Input file is fasta format
  --panacotafasta  Input file is protein fasta output from panacota
  --db DB          Which database to use: gold, regular, all (default: gold)
  --mode MODE      Mode of running: RMS, MT, RE, MT+RE (default: RMS)
  --dontcollapse   Whether to keep all blast hits for proteins rather than just their
                   top hit (default: False)
  --hmm HMM        Which HMM to use

required arguments:
  --output OUTPUT  Output prefix
```

`input`: The input file. `rmsFinder` takes either genbank or fasta files as input, which you should flag with either `--genbank` or `--fasta` flag. You can also specify `--panacotafasta` where the protein fasta is as outputted by `PanACoTA annotate` on a nucleotide fasta.

Example genbank files can be downloaded from NCBI with ```ncbi-acc-download``` (available [here](https://github.com/kblin/ncbi-acc-download/)), for example:

```
ncbi-acc-download NZ_LR025099
```

The benefit of providing a genbank file is that it includes information on the relative positions of the proteins in the genome, including contig localisation. This allows `rmsFinder` to use a positional threshold to predict the presence of a given RMS.

Alternatively, you can provide a protein fasta. In this case, `rmsFinder` will predict RMS but will not use any positional information i.e. the presence of an MTase and REase recognising the same target sequence is sufficient. (to do: add support for gff input)

`--output`: The prefix for output files. Under normal running, three output files are created: `output_MT.csv`, `output_RE.csv` and `output_RMS.csv`.

`--db`: This allows you to search against three different categories of enzyme sequences in REBASE.  ranging from those with experimental support for their restriction site ('gold') to those that have only been predicted bioinformatically based on similarity to known enzymes ('putative').
* ```gold``` - the highest standard, with experimental support for their restriction site. See [http://rebase.neb.com/cgi-bin/rebgoldlist](here).
* ```nonputative``` - only those sequences which are not putative
* ```all``` - includes putative sequences predicted bioinformatically by REBASE based on similarity to existing sequences, but for which no experimental validation is known. In many cases if you are searching an NCBI genome using `rmsFinder`, you will find a 100% match to a putative prediction for a protein. This is because the REBASE team have already identified it against the existing proteins.  

`--mode`: The default is to search for MTases and REases and then predict RMS. You can also only search for MTases (`MT`), REases (`RE`) or both but without RMS prediction (`MT+RE`).

`--dontcollapse`: The default is to keep only the top blast hit for a protein. However, it is sometimes useful to inspect all the hits that meet the threshold for a given protein.

`--hmm`: by default uses Oliveira HMMs. However, more recent HMMs can be used with [https://github.com/mdmparis/defense-finder-models/tree/master/profiles](DefenseFinder) HMMs for Type II REases and MTases can be used with `--hmm tesson`, which may recover more systems for some clades.

### Thresholds

The usual mode for `rmsFinder` is to search for MTases and/or REases at the same time (```--mode MT,RE```) and then attempt to predict the presence of Type II RMS i.e. a MTase and REase with matching specificity, close in the genome. Results for MTases and REases will be stored even if no Type II RMS are found.

There are two sets of thresholds. First, the thresholds required to identify a protein as a putative MTase or REase. These are an e-value of 0.001 and a lateral coverage of >50%.

Second, there are the thresholds for predicting the target sequence recognised by the protein. Default thresholds for the required global similarity required between a protein sequence and its best hit for the prediction of its target specificity to be accurate are those found by Oliveira et al. (2016):

* MTase: 55\%
* REase: 50\%

N.B. These global similarities are computed with `Align.PairwiseAligner` from `SeqIO`.

The default threshold for proximity in the genome is <5 genes apart.

These thresholds should be suitable for most uses but can be altered within ```rmsFinder.py``` (to do: add the option to pass a parameter file with thresholds).  


## REBASE

`rmsFinder` is only possible because of REBASE: a comprehensive and well-maintained database of known restriction enzymes and other associated enzymes. It is regularly updated by its curators who generously make it both available to all and free to distribute. `rmsFinder` is an independent piece of software and is not affiliated with REBASE.

REBASE citation:

R. J. Roberts, T. Vincze, J. Posfai, D. Macelis  
REBASE-a database for DNA restriction and modification: enzymes, genes and genomes.  
*Nucleic Acids Research* 43: D298-D299 (2015).  
doi: [10.1093/nar/gku1046](http://doi.org/10.1093/nar/gku1046)
