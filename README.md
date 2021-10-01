# rmsFinder

Scripts for finding restriction-modification systems in bacterial genomes.

Currently in development for Type II systems. Other systems are more complex - they may be added one day.

## Usage

Before running rmsFinder for the first time you will need to download the data files from REBASE:

```
python rmsFinder.py --updatedb
```

or equivalently

```
python updateDB.py
```

This should take <2 minutes and requires ~100 MB of storage in the `data` directory after it has run. REBASE is frequently updated, so you can use this command to update the database files periodically. Bear in mind your results could change if you update the database between analyses.

Normal usage of rmsFinder - to predict cognate Type II systems in a genome - looks like

```
python rmsFinder.py --genbank {YOUR_GENBANK}.gbk --mode MT,RE --output {YOUR_OUTPUT} --db all
```

The `--db` option allows you to search against three categories of REBASE protein sequences:
* `gold` - the highest standard . See [http://rebase.neb.com/cgi-bin/rebgoldlist]
* `nonputative` - only those sequences which are not putative
* `all` - includes putative sequences predicted by REBASE based on similarity to existing sequences, but for which no experimental validation is known. In many cases if you are searching an NCBI genome, you will find a 100% match to a putative prediction for a protein because REBASE have already identified it against the existing proteins.
