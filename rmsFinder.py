import pandas as pd
from Bio import SeqIO
from Bio import Align
import itertools as iter
import re
import subprocess
import os
import numpy as np
import argparse
import logging

# For absolute paths
_ROOT = os.path.abspath(os.path.dirname(__file__))

def get_options():
    parser = argparse.ArgumentParser(description='Predict presence of Type II restriction-modification systems in a genome.',
                                     prog='rmsFinder')
    input_group = parser.add_mutually_exclusive_group(required=True) # mutually exclusive group
    input_group.add_argument('--genbank', help='Genbank file') # either genbank or fasta, but not both.
    input_group.add_argument('--fasta', help='Alternatively: a fasta file (protein)')
    parser.add_argument('--output', help='Output prefix', required=True)
    parser.add_argument('--mode', help='Mode', required=True)
    parser.add_argument('--collapse', help='Whether to collapse output to best hit')
    return parser.parse_args()


def get_data(path):
    '''Returns the absolute path for a data file.
    Args:
        path (str)
            The path to the data file
    Returns:
        abs_path (str)
            The absolute path
    '''
    return os.path.join(_ROOT, 'data', path)


def searchHMM(query_protein_file, hmm_file):
    '''Searches a proteome against a specified hmm file.
    Args:
        query_protein_file (str)
            Filename of the query fasta
        hmm_file (str)
            Filename of the HMM profile
    Returns:
        df (dict)
            Dict of proteins with their top hit in the HMM profile
    '''
    # Run hmmsearch with hmmer
    tmp_file = query_protein_file+'tmp.out'
    hmmsearch_command = ['hmmsearch', '--noali', '--tblout', tmp_file,
                        hmm_file, query_protein_file]
    hmmsearch_process = subprocess.Popen(hmmsearch_command,
                        stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    # Read the output from stdout
    hmmsearch_out, _ = hmmsearch_process.communicate()
    # Parse the output
    with open(tmp_file, 'r') as f:
        df = dict()
        for line in f.readlines():
            if not line.startswith('#'):
                line = line.split()
                if line[0] in df.keys():
                     # if E-value of hit is lower than existing entry, replace
                    if float(line[4]) < float(df[line[0]][3]):
                        df[line[0]] = line[1:]
                    else: # otherwise no (keep top hit)
                        pass
                else: # if not in dict, add the hit
                    df[line[0]] = line[1:]
    os.remove(tmp_file)
    return df

def subsetFasta(input_fasta, seq_names, output_fasta):
    '''Subsets a fasta file and creates a new fasta.
    Args:
        input_fasta (str)
            The original fasta file
        seq_names (list)
            The list of sequence names to be subsetted out
        output_fasta (str)
            The output file to write the subset to
    Returns:
        None
    '''
    seqs = SeqIO.to_dict(SeqIO.parse(input_fasta, 'fasta'))
    subset_seqs = [seqs[record] for record in seq_names]
    with open(output_fasta, 'w') as output_file:
        for record in subset_seqs:
            output_file.write('>%s\n%s\n' % (record.id, str(record.seq)))
    return


def blastpAgainstDB(query_fasta, db_fasta, db_built=True, evalue_threshold=0.001, format_string='qseqid sseqid pident length qlen evalue', top_hits_only=True):
    '''Blasts a fasta of query proteins against a database and returns the hits.
    Args:
        query_fasta (str)
            Filename of query fasta
        db_fasta (str)
            Filename of DB fasta
        db_built (Bool)
            Indicates whether the database is already constructed.
        evalue_threshold (float)
            Threshold to use with blastp
        format_string (str)
            How to format the blast output (outfmt 6)
    Returns:
        blast_results_top (pd dataframe)
            Dataframe of top hits for each protein in query fasta (can be multiple for each if there are multiple top hits)
    '''
    if db_built==False: # Build blastdb if have been told it doesn't exist
        makeblastdb_command = ['makeblastdb',
                            '-in', db_fasta,
                            '-dbtype', 'prot']
        makeblastdb_process = subprocess.Popen(makeblastdb_command,
                            stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        makeblastdb_out, _ = makeblastdb_process.communicate() # Read the output from stdout
    # Run blastp
    blastp_command = ['blastp',
                    '-query', query_fasta,
                    '-db', db_fasta,
                    '-outfmt', '6 '+format_string,
                    '-evalue', str(evalue_threshold)]
    blastp_process = subprocess.Popen(blastp_command,
                        stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    blastp_out, _ = blastp_process.communicate() # Read the output from stdout
    blastp_output = re.split('\n|\t',blastp_out.decode()) # Decode
    # Length of format string
    N_columns = len(format_string.split())


    if len(blastp_output)>N_columns: # should be N_column entries plus empty entry at end of output, so more implies at least one hit
        blastp_output.pop(-1) # remove empty last entry
        blast_results = pd.DataFrame(np.reshape(blastp_output, newshape=(int(np.floor(len(blastp_output)/N_columns)), N_columns))) # reformat
        blast_results.columns = format_string.split() # rename columns
        blast_results = blast_results.astype({'qseqid': 'str', 'sseqid': 'str',
                                                'pident': 'float64','length': 'int64',
                                                'qlen': 'int64', 'evalue':'str' })
                                                # Need the scientific format, so store evalue as string for safety if need later
        # change dtypes
        # max hit dict - for each protein
        if top_hits_only==True:
            max_hit_dict = {k: max(blast_results[blast_results['qseqid']==k]['pident']) for k in set(blast_results['qseqid'])}
            # Keep only top hits for each protein
            blast_results_top = pd.concat([blast_results[[blast_results['qseqid'][i]==k and blast_results['pident'][i]==max_hit_dict[k] for i in range(len(blast_results))]] for k in max_hit_dict.keys()])
            return(blast_results_top)
        else:
            return(blast_results)
    else:
        return(None)

def parseGenBank(genbank_file, genbank2fasta_output):
    '''Parses a GenBank file into a fasta of proteins with useful information.
    Args:
        genbank_file (str)
            File in genbank format. As downloaded with e.g.
                $ ncbi-acc-download NZ_LR025099
            If in gzipped format, parseGenbank will gunzip then gzip at end.
        genbank2fasta_output (str)
            Output fasta file
    Returns:
        None
    '''
    protein_dict = {}
    counter = 0
    gzipFlag = False
    if genbank_file.endswith('.gz'): # Gunzip if we need to
        gzipFlag = True
        gunzip_command = ['gunzip', genbank_file]
        gunzip_process = subprocess.Popen(gunzip_command,
                            stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        gunzip_process.wait()
        genbank_file = genbank_file[:-3]
    with open(genbank2fasta_output, 'w') as f:
        observed_proteins = []
        for record in SeqIO.parse(genbank_file, 'genbank'): # read assumes only one entry
            for feature in record.features:
                if feature.type=='CDS':
                    counter += 1
                    if 'translation' in feature.qualifiers.keys():
                        protein_id = feature.qualifiers['protein_id'][0]
                        if protein_id in observed_proteins:
                            pass
                        else:
                            f.write('>%s %s product="%s"\n%s\n' % (protein_id, counter, feature.qualifiers['product'][0], feature.qualifiers['translation'][0]))
                            observed_proteins.append(protein_id)
    if gzipFlag==True: # Gzip if we gunzipped
        gzip_command = ['gzip', genbank_file]
        gzip_process = subprocess.Popen(gzip_command,stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        gzip_process.wait()
    return


def getRS(queries, fasta_file):
    '''Returns the RecSeq for a set of queries in a (REBASE) fasta file.'''
    seqs = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    description_queries = [seqs[query].description.split('\t') for query in queries]
    rec_seqs = [[x for x in description_query if 'RecSeq' in x][0] for description_query in description_queries]
    rec_seqs = [rec_seq.split(':')[1] for rec_seq in rec_seqs]
    return(rec_seqs)

def globalSimilarity(seq_a, seq_b):
    '''Returns the optimal global alignment score for two sequences.
    Args:
        seq_a, seq_b (str)
            The sequences to align

    Returns:
        pident (float)
            The % similarity score (using the length of the shortest sequence)
    '''
    aligner = Align.PairwiseAligner()
    global_score = aligner.score(seq_a, seq_b)
    min_length = min([len(seq_a), len(seq_b)])
    global_pident = global_score/min_length * 100
    return global_pident


def collapseBestHits(best_hits):
    '''Collapses a table of best hits to just a single entry for each protein.
    Args:
        best_hits (DataFrame)
            Pandas dataframe of best hits, from blastpAgainstDB
    Returns:
        best_hits_collapse (DataFrame)
            Dataframe with one row per protein but all information stored
    '''
    best_hits_collapse = pd.concat([best_hits[best_hits['qseqid']==x].head(1) for x in set(best_hits['qseqid'])])
    # Add alternative seq_ids
    best_hits_collapse = best_hits_collapse.assign(other_REBASE_hits = [','.join(list(best_hits[best_hits['qseqid']==x]['sseqid'])[1:]) for x in set(best_hits['qseqid'])])
    best_hits_collapse = best_hits_collapse.assign(n_REBASE_hits = [len(list(best_hits[best_hits['qseqid']==x]['sseqid'])) for x in set(best_hits['qseqid'])])
    # Check if all target sequences are the same
    best_hits_collapse = best_hits_collapse.assign(identicalTarget = [len(set(best_hits[best_hits['qseqid']==x]['target']))==1 for x in best_hits_collapse['qseqid']])
    return(best_hits_collapse)

def parseCounterPreparedFasta(input_fasta):
    '''Gets counters for fasta prepared from parseGenBank.'''
    counter_dict = {}
    for line in open(input_fasta,'r').readlines():
        if line.startswith('>'):
            entries = line.split()
            protein_id = re.sub('>', '', entries[0])
            counter_dict[protein_id] = int(entries[1])
    return(counter_dict)

def parseCDSFromGenomic(input_fasta):
    '''Takes a _cds_from_genomic.fna file downloaded from NCBI and produces a dict of the proteins with their locations.
    Args:
        input_fasta (str)
            Filename of _cds_from_genomic.fna file
    Returns:
        counters (dict)
            Dictionary matching protein_id to counter number (for genomic position)
    '''
    counters = {}
    for line in open(input_fasta, 'r').readlines():
        if line.startswith('>'):
            entries = line.split()
            where_is_protein = ['protein_id' in x for x in entries]
            if any(where_is_protein)==True:
                protein_id = entries[where_is_protein.index(True)]
                protein_id =  re.sub('\]', '', re.sub('\[_=', '', re.sub('[a-z]', '', protein_id)))
                counter = int(re.sub('.*_', '', entries[0]))
                counters[protein_id] = counter
    return(counters)

def predictRMS(hits_MT, hits_RE, position_threshold=5, mt_threshold=55, re_threshold=50):
    '''Predicts RM system based on tables of hits to MTases and REases.
    Args:
        hits_MT (DataFrame)
            Hits to MTases
        hits_RE (DataFrame)
            Hits to REases
        position_threshold (int)
            Proteins need to be < this threshold to count as system present
        Optional: mt_threshold, re_threshold (int)
            Similarity thresholds required to rely on prediction of target sequence.
            Default are the values 55% (MTase) and 50% (REase) from Oliveira 2016
    Returns:
        predicted_rms (list)
            Target sequences as keys with MT and RE proteins and positions stored as values
    '''
    # Filter hits based on Oliveira thresholds
    hits_MT = hits_MT[hits_MT['similarity']>mt_threshold]
    hits_RE = hits_RE[hits_RE['similarity']>re_threshold]
    # Check for any intersection of targets
    target_overlap = set(hits_MT['target']).intersection(set(hits_RE['target']))
    if len(target_overlap) > 0:
        predicted_rms = []
        for t in target_overlap:
            MT_positions = list(hits_MT[hits_MT['target']==t]['position'])
            RE_positions = list(hits_RE[hits_RE['target']==t]['position'])
            # Want all pairwise combinations of these
            separations = [(x, y, abs(x-y)) for x in MT_positions for y in RE_positions] # List of tuples storing position of MT and RE and separation

            for s in separations:
                if s[2]<position_threshold:
                    rms_entry = [t, s[0], s[1], list((hits_MT[hits_MT['position']==s[0]]['qseqid']))[0], list((hits_RE[hits_RE['position']==s[1]]['qseqid']))[0]]
                    if rms_entry[3]==rms_entry[4]: # If the hit is for the same protein, don't include it (only want paired RE and MT)
                        pass
                    else:
                        predicted_rms.append(rms_entry)
        logger.info(predicted_rms)
        if len(predicted_rms)!=0:
            rms_results = pd.DataFrame(predicted_rms, columns=['sequence', 'pos_MT', 'pos_RE', 'prot_MT', 'prot_RE'])
            logger.info(rms_results)
            # Add similarity scores and best hit
            rms_results['sim_MT'] = rms_results.apply(lambda row : hits_MT[hits_MT['qseqid']==row['prot_MT']]['similarity'], axis=1)
            rms_results['hit_MT'] = rms_results.apply(lambda row : hits_MT[hits_MT['qseqid']==row['prot_MT']]['sseqid'], axis=1)

            rms_results['sim_RE'] = rms_results.apply(lambda row : hits_RE[hits_RE['qseqid']==row['prot_RE']]['similarity'], axis=1)
            rms_results['hit_RE'] = rms_results.apply(lambda row : hits_RE[hits_RE['qseqid']==row['prot_RE']]['sseqid'], axis=1)
            return(rms_results)
        else:
            return(None)
    else:
        return(None)

def searchMTasesTypeII(proteome_fasta, cds_from_genomic_fasta=False, evalue_threshold=0.001, coverage_threshold=0.5, collapse=True, MTase_db='protein_seqs_Type_II_MTases.faa'):
    '''Searches for Type II MTases.
    Args:
        proteome_fasta (str)
            Fasta file of proteins.
        cds_from_genomic_fasta (str)
            Fasta file of untranslated nucleotides of cds (contains counter info.)
        evalue_threshold (float)
            Threshold to keep hits at. Default 0.001 as in Oliveira 2016
        coverage_threshold (float)
            Threshold of coverage. Default: 0.5 (i.e. 50%) as in Oliveira 2016
    Returns:
        blast_hits_collapse (DataFrame)
            DataFrame of best hits, collapsed to one row per protein
    '''
    MTase_db_file = get_data('protein_seqs_Type_II_MTases.faa')

    # Using Oliveira Type II MTase HMM profiles to search
    hmm_dict_MT = searchHMM(proteome_fasta, get_data('Type_II_MTases.hmm'))
    logging.info('Found %d raw hits for MTases.' % len(hmm_dict_MT))
    #print(hmm_dict_MT)

    # Filter hits
    hits_MT_filt = {k:v for k,v in hmm_dict_MT.items() if float(v[3])<evalue_threshold}
    logging.info('Found %d filtered hits for MTases.' % len(hits_MT_filt))
    #print(hits_MT_filt)

    # Subset only the hits out from the proteome
    tmp_fasta = proteome_fasta+'tmp_MT.faa'
    subsetFasta(proteome_fasta, list(hits_MT_filt.keys()), tmp_fasta)

    # Blast these hits against all Type II MTases to find best matches
    blast_hits_MT = blastpAgainstDB(tmp_fasta, MTase_db_file, db_built=True)
    # Store the sequences for global alignment
    protein_seqs = SeqIO.to_dict(SeqIO.parse(tmp_fasta, 'fasta'))
    rebase_seqs = SeqIO.to_dict(SeqIO.parse(MTase_db_file, 'fasta'))
    # Remove tmp fasta file
    os.remove(tmp_fasta)
    logging.info('Found %d best matches for MTases.' % len(blast_hits_MT))
    #print(blast_hits_MT)

    # If no hits?
    if blast_hits_MT is None:
        return
    else:
        # Filter coverage threshold - add
        blast_hits_MT = blast_hits_MT.assign(coverage_threshold_met=list(blast_hits_MT['length'] > coverage_threshold*blast_hits_MT['qlen'])) # Condition of 50% coverage as in Oliveira 2016

        # Get the recognition sites of the best hits
        rs_MT = getRS(blast_hits_MT['sseqid'], MTase_db_file)
        blast_hits_MT = blast_hits_MT.assign(target=rs_MT)

        # Add genomic position if requested
        if cds_from_genomic_fasta==True:
            counter_dict = parseCounterPreparedFasta(proteome_fasta)
            blast_hits_MT = blast_hits_MT.assign(position=[counter_dict[x] for x in blast_hits_MT['qseqid']])

        # Add the global similarity of the best hit. Need to have the sequences available
        blast_hits_MT['similarity'] = blast_hits_MT.apply(lambda row : globalSimilarity(str(protein_seqs[row['qseqid']].seq),
                     str(rebase_seqs[row['sseqid']].seq)), axis = 1)

        # Collapse the table to best hits
        if collapse==True:
            blast_hits_collapse = collapseBestHits(blast_hits_MT)
            return(blast_hits_collapse)
        else:
            return(blast_hits_MT)

def searchREasesTypeII(proteome_fasta, cds_from_genomic_fasta=False, evalue_threshold=0.001, coverage_threshold=0.5, collapse=True, REase_db='protein_seqs_Type_II_REases.faa'):
    '''Searches a file of proteins against all known REases.
    Args:
        proteome_fasta (str)
            Fasta file with proteins
        cds_from_genomic_fasta (str)
            Fasta file of untranslated nucleotides of cds
        evalue_threshold (float)
            Threshold to filter blast hits at. Default: 0.001 as in Oliveira 2016
        coverage_threshold (float)
            Threshold of coverage. Default: 0.5 (i.e. 50%) as in Oliveira 2016
    Returns:
        blast_hits_collapse (DataFrame)
            DataFrame of best hits, one row per protein
    '''
    # Blasting for REases
    REase_db_file = get_data(REase_db)
    blast_hits_RE = blastpAgainstDB(proteome_fasta, REase_db_file, db_built=True)
    # Store the sequences for global alignment
    protein_seqs = SeqIO.to_dict(SeqIO.parse(proteome_fasta, 'fasta'))
    rebase_seqs = SeqIO.to_dict(SeqIO.parse(REase_db_file, 'fasta'))

    # Filter out hits
    blast_hits_RE = blast_hits_RE.assign(coverage_threshold_met=list(blast_hits_RE['length'] > coverage_threshold*blast_hits_RE['qlen'])) # Condition of 50% coverage as in Oliveira 2016
    blast_hits_RE_filt = blast_hits_RE[blast_hits_RE['coverage_threshold_met']==True]

    # Add genomic position, if requested
    if cds_from_genomic_fasta==True:
        counter_dict = parseCounterPreparedFasta(proteome_fasta)
        blast_hits_RE_filt = blast_hits_RE_filt.assign(position=[counter_dict[x] for x in blast_hits_RE_filt['qseqid']])

    # Add the recognition sequences
    blast_hits_RE_filt = blast_hits_RE_filt.assign(target=getRS(blast_hits_RE_filt['sseqid'], REase_db_file))

    # Add the global similarity of the best hit
    blast_hits_RE_filt['similarity'] = blast_hits_RE_filt.apply(lambda row : globalSimilarity(str(protein_seqs[row['qseqid']].seq),
                 str(rebase_seqs[row['sseqid']].seq)), axis = 1)

    # Collapse the table to best hits
    if collapse==True:
        blast_hits_collapse = collapseBestHits(blast_hits_RE_filt)
        return(blast_hits_collapse)
    else:
        return(blast_hits_RE_filt)


def main():
    args = get_options()
    output = args.output
    mode = args.mode
    collapse_hits = True
    if args.collapse=='F':
        collapse_hits = False

    # Logger
    level = logging.INFO
    format = '  %(message)s'
    handlers = [logging.StreamHandler()]
    logging.basicConfig(level = level, format = format, handlers = handlers)
    logging.info('Started running.')

    if args.genbank is not None:
        genbank_file = args.genbank
        proteome_fasta = genbank_file+'.tmp.faa'
        parseGenBank(genbank_file, proteome_fasta) # Make fasta file the way we like it
    elif args.fasta is not None:
        proteome_fasta = args.fasta

    if 'MT' in mode: # Search for MTases
        MT_hits = searchMTasesTypeII(proteome_fasta, True, collapse=collapse_hits)
        if MT_hits is not None:
            MT_hits.to_csv(output+'_MT.csv', index=False, float_format="%.3f")
        else:
            logging.info('No MTase hits.')
        logging.info('Finished searching for MTases.')
    if 'RE' in mode: # Search for REases
        RE_hits = searchREasesTypeII(proteome_fasta, True)
        RE_hits.to_csv(output+'_RE.csv', index=False, float_format="%.3f")
        logging.info('Finished searching for REases.')
    if 'RE' in mode and 'MT' in mode: # Predict R-M systems if both searched for
        rms_predictions = predictRMS(MT_hits, RE_hits)
        #print(rms_predictions)
        if rms_predictions is not None:
            logging.info('Predicted presence of %d Type II R-M systems.' % len(rms_predictions))
            rms_predictions.to_csv(output+'_RMS.csv', index=False, float_format="%.3f")
        else:
            logging.info('Predicted no Type II R-M systems.')
    os.remove(proteome_fasta) # Remove the proteome fasta we made


if __name__ == "__main__":
    main()
