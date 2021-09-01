from Bio import SeqIO
from Bio import Align
import re
import subprocess
import logging
import wget
import glob, os
import datetime


def convertRebaseToFasta(rebase_file, output_fasta):
    '''Converts downloaded Rebase files to fasta.'''
    headerFlag = True
    with open(output_fasta, 'w') as output:
        for line in open(rebase_file, 'r').readlines():
            if line.startswith('>'):
                headerFlag = False
            if headerFlag==False:
                line = line.strip('\n')
                if line.startswith('>'):
                    line = re.sub(' \(', '(', line) # to ensure unique name in fasta for proteins which have multiple subunits
                    output.write('\n%s\n' % line)
                else:
                    line = re.sub(' ', '', line)
                    line = re.sub('<>', '', line)
                    _ = output.write(line)

def extractEnzymesSpecifiedType(fasta_file, type, output_fasta):
    '''Extracts only protein sequences of a particular type from a converted
    REBASE fasta file (assumes this information is in the fasta descriptions).
    Only stores those with a known RecSeq.'''
    rebase_seqs = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    type_names = [seq for seq in rebase_seqs if type in rebase_seqs[seq].description and 'RecSeq' in rebase_seqs[seq].description]
    type_seqs = [rebase_seqs[record] for record in type_names]
    with open(output_fasta, 'w') as output:
        for record in type_seqs:
            _ = output.write('>%s\n%s\n' % (str(record.description), str(record.seq)))

def makeBlastDB(db_fasta):
    '''Makes blast database from a fasta file.'''
    makeblastdb_command = ['makeblastdb',
                        '-in', db_fasta,
                        '-dbtype', 'prot']
    makeblastdb_process = subprocess.Popen(makeblastdb_command,
                        stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    makeblastdb_out, _ = makeblastdb_process.communicate() # Read the output from stdout
    return

def downloadFromREBASE(file, output):
    '''Downloads a file from REBASE.'''
    rebase_url = 'ftp://ftp.neb.com/pub/rebase/'
    attempts = 1
    while attempts<4:
        try:
            wget.download(rebase_url+file, out=output)
            break
        except Exception as e:
            print("\n  Seemed to encounter error while downloading file: ", file)
            print("  Error code:", e)
            print("  (if the file downloaded, then you can ignore this)")
            attempts += 1
    return

def writeLookupTables(gold_fasta, nonputative_fasta, all_fasta, output):
    '''Writes lookup table.'''
    proteins_gold = SeqIO.to_dict(SeqIO.parse(gold_fasta, 'fasta')).keys()
    proteins_nonputative = SeqIO.to_dict(SeqIO.parse(nonputative_fasta, 'fasta')).keys()
    proteins_all = SeqIO.to_dict(SeqIO.parse(all_fasta, 'fasta')).keys()
    proteins_putative = [x for x in proteins_all if not x in proteins_nonputative]
    proteins_nonputative = [x for x in proteins_nonputative if x not in proteins_gold]
    proteins_lookup_dict = {x: 'gold' for x in proteins_gold}
    proteins_lookup_dict.update({x: 'nonputative' for x in proteins_nonputative})
    proteins_lookup_dict.update({x: 'putative' for x in proteins_putative})
    with open(output, 'w') as f:
        for k, v in proteins_lookup_dict.items():
            _ = f.write('%s %s\n' % (k, v))
    return

def main():
    # Logger details
    level = logging.INFO
    format = '%(message)s'
    handlers = [logging.StreamHandler()]
    logging.basicConfig(level = level, format = format, handlers = handlers)
    logging.info('Started updating local REBASE databases.')

    # Downloading
    logging.info('Downloading all REBASE protein sequences.')
    downloadFromREBASE('protein_seqs.txt',
                output='data/All.txt')

    logging.info('\nDownloading all Gold protein sequences.')
    downloadFromREBASE('protein_gold_seqs.txt',
                output='data/Gold.txt')

    # Converting
    logging.info('\n\nConverting REBASE files to fasta format...')
    if os.path.exists('data/All.txt'):
        convertRebaseToFasta('data/All.txt', 'data/All.faa')
    if os.path.exists('data/Gold.txt'):
        convertRebaseToFasta('data/Gold.txt', 'data/Gold.faa')
    logging.info('Done!')

    # Extracting Type II
    logging.info('\nExtracting Type II sequences...')
    extractEnzymesSpecifiedType('data/All.faa', ':Type II methyltransferase', 'data/Type_II_MT_nonputative.faa') # excludes putative
    extractEnzymesSpecifiedType('data/All.faa', ':Type II restriction enzyme', 'data/Type_II_RE_nonputative.faa') # excludes putative
    extractEnzymesSpecifiedType('data/All.faa', 'Type II methyltransferase', 'data/Type_II_MT_all.faa')
    extractEnzymesSpecifiedType('data/All.faa', 'Type II restriction enzyme', 'data/Type_II_RE_all.faa')
    logging.info('Done!')
    # Extracting Gold Type II
    logging.info('\nExtracting Gold Type II sequences...')
    extractEnzymesSpecifiedType('data/Gold.faa', 'Type II methyltransferase', 'data/Type_II_MT_gold.faa')
    extractEnzymesSpecifiedType('data/Gold.faa', 'Type II restriction enzyme', 'data/Type_II_RE_gold.faa')
    logging.info('Done!')

    # Making lookup tables
    writeLookupTables('data/Type_II_MT_gold.faa', 'data/Type_II_MT_nonputative.faa', 'data/Type_II_MT_all.faa', 'data/Type_II_MT_dict.txt')
    writeLookupTables('data/Type_II_RE_gold.faa', 'data/Type_II_RE_nonputative.faa', 'data/Type_II_RE_all.faa', 'data/Type_II_RE_dict.txt')

    # Make blast databases
    logging.info('\nMaking blast databases...')
    makeBlastDB('data/Type_II_MT_all.faa')
    makeBlastDB('data/Type_II_RE_all.faa')
    makeBlastDB('data/Type_II_MT_nonputative.faa')
    makeBlastDB('data/Type_II_RE_nonputative.faa')
    makeBlastDB('data/Type_II_MT_gold.faa')
    makeBlastDB('data/Type_II_RE_gold.faa')
    logging.info('Done!')

    # Removing REBASE files
    os.remove('data/All.txt')
    os.remove('data/All.faa')
    os.remove('data/Gold.txt')
    os.remove('data/Gold.faa')
    for f in glob.glob('data/*.tmp'):
        os.remove(f)

    # Record when files where downloaded
    with open('data/download.log', 'w') as f:
        now = datetime.datetime.now()
        format = "%d/%m/%Y %H:%M:%S"
        f.write('Databases last downloaded on {}\n'.format(now.strftime(format)))


if __name__ == "__main__":
    main()
