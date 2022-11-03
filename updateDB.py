from Bio import SeqIO
from Bio import Align
import re
import subprocess
import logging
import wget
import glob, os
import datetime
import rmsFinder as rf
import argparse

def get_options():
    parser = argparse.ArgumentParser(description='Download REBASE sequences and set up databases for Type II RM systems.',
                                     prog='updateDB')
    parser.add_argument('--force', help='Force overwrite.', required=False, action='store_true')
    parser.add_argument('--recompile', help='Only make the blast databases again (without redownloading files).', required=False, action='store_true')
    parser.add_argument('--keepall',
                    help='Only make the blast databases again (without redownloading files).',
                    required=False, action='store_true')

    return parser.parse_args()

def removeFile(file_str):
    '''Removes a file (after checking it exists first.)'''
    if os.path.exists(file_str):
        os.remove(file_str)
        return



def convertRebaseToFasta(rebase_file, output_fasta):
    '''Converts downloaded Rebase files to fasta.
    Args:
        rebase_file (str)
            The downloaded filename.
        output_fasta (str)
            The filename of the output fasta.
    Returns:
        None
    '''
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
    return

def extractEnzymesSpecifiedType(fasta_file, type, output_fasta):
    '''Extracts only protein sequences of a particular type from a converted
    REBASE fasta file (assumes this information is in the fasta descriptions).
    Only stores those with a known RecSeq.
    Args:
        fasta_file (str)
            Fasta filename
        type (str)
            The desired type of sequence to subset e.g. "Type II methyltransferase"
        output_fasta (str)
            The output fasta filename
    Returns:
        None
    '''
    rebase_seqs = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    type_names = [seq for seq in rebase_seqs if type in rebase_seqs[seq].description and 'RecSeq' in rebase_seqs[seq].description]
    type_seqs = [rebase_seqs[record] for record in type_names]
    with open(output_fasta, 'w') as output:
        for record in type_seqs:
            _ = output.write('>%s\n%s\n' % (str(record.description), str(record.seq)))
    return

def makeBlastDB(db_fasta, outdir=rf.get_data('db')):
    '''Makes blast database from a fasta file.
    Args:
        db_fasta (str)
            The fasta filename
    Returns:
        None
    '''
    makeblastdb_command = ['makeblastdb',
                        '-in', db_fasta,
                        '-dbtype', 'prot',
                        '-out', outdir+'/'+re.sub('.*\/', '', db_fasta)]
    makeblastdb_process = subprocess.Popen(makeblastdb_command,
                        stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    makeblastdb_out, _ = makeblastdb_process.communicate() # Read the output from stdout
    return

def downloadFromREBASE(file, output):
    '''Downloads a file from REBASE.
    Args:
        file (str)
            The filename on REBASE (full path is added)
        output (str)
            Where to save the output.
    Returns:
        None
    '''
    rebase_url = 'ftp://ftp.neb.com/pub/rebase/'
    attempts = 1
    while attempts<4:
        try:
            wget.download(rebase_url+file, out=output)
            break
        except Exception as e:
            print("\n  Warning: seemed to encounter error while downloading file: ", file)
            print("  Error code:", e)
            print("  Will try again. (If the file downloaded based on the byte count, then you can ignore this.)")
            attempts += 1
    return

def writeLookupTables(gold_fasta, nonputative_fasta, all_fasta, output):
    '''Writes lookup table for the different categories of REBASE sequence.
    Args:
        gold_fasta (str)
            Fasta file with the Gold category proteins.
        nonputative_fasta (str)
            Fasta file with the nonputative category proteins.
        all_fasta (str)
            Fasta file with all proteins.
        output (str)
            Filename to save the text output in.
    Returns:
        None
    '''
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
    args = get_options()
    if args.force==False:
        if os.path.exists(rf.get_data('download.log')):
            last_downloaded = open(rf.get_data('download.log'), 'r').readlines()[0]
            print(last_downloaded)
            print('Do you still want to continue?\nIf you want to download REBASE databases again and overwrite these older versions, use --force.')
            return
    # Logger details
    level = logging.INFO
    format = '%(message)s'
    handlers = [logging.StreamHandler()]
    logging.basicConfig(level = level, format = format, handlers = handlers)
    logging.info('Started updating local REBASE databases.')
    if os.path.exists(rf.get_data('download.log')):
        last_downloaded = open(rf.get_data('download.log'), 'r').readlines()[0]
        logging.info(last_downloaded+'. \nThese databases are being overwritten!\n\n')
    if not os.path.exists(rf.get_data('db')):
        os.mkdir(rf.get_data('db'))

    # Downloading
    if not args.recompile:
        logging.info('Downloading all REBASE protein sequences.')
        downloadFromREBASE('protein_seqs.txt',
                    output=rf.get_data('All.txt'))

        logging.info('\nDownloading all Gold protein sequences.')
        downloadFromREBASE('protein_gold_seqs.txt',
                    output=rf.get_data('Gold.txt'))

        # Converting
        logging.info('\n\nConverting REBASE files to fasta format...')
        if os.path.exists(rf.get_data('All.txt')):
            convertRebaseToFasta(rf.get_data('All.txt'), rf.get_data('All.faa'))
        if os.path.exists(rf.get_data('Gold.txt')):
            convertRebaseToFasta(rf.get_data('Gold.txt'), rf.get_data('Gold.faa'))
        logging.info('Done!')

        # Extracting Type II
        logging.info('\nExtracting Type II sequences...')
        extractEnzymesSpecifiedType(rf.get_data('All.faa'), ':Type II methyltransferase', rf.get_data('Type_II_MT_nonputative.faa')) # excludes putative
        extractEnzymesSpecifiedType(rf.get_data('All.faa'), ':Type II restriction enzyme', rf.get_data('Type_II_RE_nonputative.faa')) # excludes putative
        extractEnzymesSpecifiedType(rf.get_data('All.faa'), 'Type II methyltransferase', rf.get_data('Type_II_MT_all.faa'))
        extractEnzymesSpecifiedType(rf.get_data('All.faa'), 'Type II restriction enzyme', rf.get_data('Type_II_RE_all.faa'))
        logging.info('Done!')
        # Extracting Gold Type II
        logging.info('\nExtracting Gold Type II sequences...')
        extractEnzymesSpecifiedType(rf.get_data('Gold.faa'), 'Type II methyltransferase', rf.get_data('Type_II_MT_gold.faa'))
        extractEnzymesSpecifiedType(rf.get_data('Gold.faa'), 'Type II restriction enzyme', rf.get_data('Type_II_RE_gold.faa'))
        logging.info('Done!')

        # Making lookup tables
        writeLookupTables(rf.get_data('Type_II_MT_gold.faa'), rf.get_data('Type_II_MT_nonputative.faa'), rf.get_data('Type_II_MT_all.faa'), rf.get_data('Type_II_MT_dict.txt'))
        writeLookupTables(rf.get_data('Type_II_RE_gold.faa'), rf.get_data('Type_II_RE_nonputative.faa'), rf.get_data('Type_II_RE_all.faa'), rf.get_data('Type_II_RE_dict.txt'))

    # Make blast databases
    logging.info('\nMaking blast databases...')
    makeBlastDB(rf.get_data('Type_II_MT_all.faa'))
    makeBlastDB(rf.get_data('Type_II_RE_all.faa'))
    makeBlastDB(rf.get_data('Type_II_MT_nonputative.faa'))
    makeBlastDB(rf.get_data('Type_II_RE_nonputative.faa'))
    makeBlastDB(rf.get_data('Type_II_MT_gold.faa'))
    makeBlastDB(rf.get_data('Type_II_RE_gold.faa'))
    logging.info('Done!')

    # Removing REBASE files
    if not args.keepall:
        removeFile(rf.get_data('All.txt'))
        removeFile(rf.get_data('All.faa'))
        removeFile(rf.get_data('Gold.txt'))
        removeFile(rf.get_data('Gold.faa'))
    for f in glob.glob(rf.get_data('*.tmp')):
        os.remove(f)

    logging.info('\nDatabases all downloaded, rmsFinder is ready to run. \n\nRemember to cite REBASE.\n\nRoberts, R.J., Vincze, T., Posfai, J., Macelis, D.\nREBASE-a database for DNA restriction and modification: enzymes, genes and genomes.\nNucleic Acids Res. 43: D298-D299 (2015).\ndoi: 10.1093/nar/gku1046\nOfficial REBASE web site - http://rebase.neb.com')

    # Record when files where downloaded
    with open(rf.get_data('download.log'), 'w') as f:
        now = datetime.datetime.now()
        format = "%d/%m/%Y %H:%M:%S"
        f.write('Databases last downloaded on {}\n'.format(now.strftime(format)))


if __name__ == "__main__":
    main()
