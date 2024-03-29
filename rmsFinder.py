import pandas as pd
from Bio import SeqIO
from Bio import Align
from Bio.Align import substitution_matrices
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
    parser.add_argument('input', nargs='+', help='Input protein file (genbank or fasta).')
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('--output', help='Output prefix', required=True)
    input_group = parser.add_mutually_exclusive_group(required=True) # mutually exclusive group
    input_group.add_argument('--genbank', help='Input file is genbank format', action='store_true', default=False) # either genbank, fasta, or panacotafasta, but not both.
    input_group.add_argument('--fasta', help='Input file is fasta format', action='store_true', default=False)
    input_group.add_argument('--panacotafasta', help='Input file is protein fasta output from panacota', action='store_true', default=False)
    input_group.add_argument('--transcdsfasta', help='Input file is translated CDS fasta from NCBI', action='store_true', default=False)
    parser.add_argument('--db', help='Which database to use: gold, regular, all (default: gold)', required=False, default='gold')
    parser.add_argument('--mode', help='Mode of running: RMS, MT, RE, MT+RE, IIG (default: RMS)', required=False, default='RMS')
    parser.add_argument('--dontcollapse', help='Whether to keep all blast hits for proteins rather than just their top hit (default: False)', action='store_true')
    parser.add_argument('--hmm', help='Which HMM to use', required=False, default='oliveira')
    parser.add_argument('--forceprocessing', help='Forces processing of the input', required=False, action='store_true')
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

def makeTmpFile(file_path, suffix, prefix='TMP'):
    '''Makes a TMP_ file from a given file descriptor, taking path into account
    so that TMP_ file will be in same directory.'''
    if '/' in file_path:
        file_str = re.sub('.*/', '', file_path)
        preamble_str = file_path[:-len(file_str)]
        tmp_path = preamble_str+'TMP_'+file_str+'.'+suffix
    else:
        tmp_path = 'TMP_'+file_path+'.'+suffix
    return tmp_path

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
    tmp_file = makeTmpFile(query_protein_file, 'out')
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
        for record in SeqIO.parse(genbank_file, 'genbank'): # there can be multiple records
            record_name = record.name
            record_description = record.description
            record_type = 'unknown'
            if 'chromosome' in record_description:
                record_type = 'chromosome'
            elif 'plasmid' in record_description:
                record_type = 'plasmid'
            counter = 0
            for feature in record.features:
                if feature.type=='CDS':
                    counter += 1
                    if 'translation' in feature.qualifiers.keys():
                        protein_id = feature.qualifiers['protein_id'][0]
                        if protein_id in observed_proteins: # don't write if already written
                            pass
                        else:
                            f.write('>%s %s %s %s location=%s product="%s"\n%s\n' % (protein_id, record_name, counter, record_type,
                            str(feature.location), feature.qualifiers['product'][0],
                            feature.qualifiers['translation'][0]))
                            observed_proteins.append(protein_id)
    if gzipFlag==True: # Gzip if we gunzipped
        gzip_command = ['gzip', genbank_file]
        gzip_process = subprocess.Popen(gzip_command,stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        gzip_process.wait()
    return

def parsePanacota(panacota_fasta_file, panacota2fasta_output):
    '''Parses a Panacota fasta file into a fasta of proteins formatted the expected way for rmsFinder.
    Args:
        panacota_fasta_file (str)
            Panacota protein file
        panacota2fasta_output (str)
            Output fasta file
    Returns:
        None
    '''
    protein_dict = {}
    with open(panacota2fasta_output, 'w') as f:
        observed_proteins = []
        for record in SeqIO.parse(panacota_fasta_file, 'fasta'): # there can be multiple records
            protein_id = record.name
            protein_description = record.description
            record_type = 'unknown'
            if 'chromosome' in protein_description:
                record_type = 'chromosome'
            elif 'plasmid' in protein_description:
                record_type = 'plasmid'
            counter = int(re.sub('.*_', '', protein_id))
            record_name = re.sub('_.*', '', protein_id)[:-1]
            f.write('>%s %s %s %s location=%s product="%s"\n%s\n' % (protein_id, record_name, counter, record_type,'unknown', 'unknown', str(record.seq)))
            observed_proteins.append(protein_id)
    return

def parseNCBITranslatedCDSFasta(translated_cds_fasta_file, ncbi2fasta_output):
    '''Parses a NCBI translated CDS fasta file into a fasta of proteins formatted the expected way for rmsFinder.
    Args:
        translated_cds_fasta_file (str)
            NCBI protein file of translated CDS (_translated_cds.faa)
        ncbi2fasta_output (str)
            Output fasta file
    Returns:
        None
    '''
    protein_dict = {}
    with open(ncbi2fasta_output, 'w') as f:
        observed_proteins = []
        for record in SeqIO.parse(translated_cds_fasta_file, 'fasta'): # there can be multiple records
            protein_id = record.name
            protein_description = record.description
            record_type = 'unknown'
            if 'chromosome' in protein_description:
                record_type = 'chromosome'
            elif 'plasmid' in protein_description:
                record_type = 'plasmid'
            counter = int(re.sub('.*_', '', protein_id))
            record_name = re.sub('_prot.*', '', protein_id)
            f.write('>%s %s %s %s location=%s product="%s"\n%s\n' % (protein_id, record_name, counter, record_type,'unknown', 'unknown', str(record.seq)))
            observed_proteins.append(protein_id)
    return


def getRS(queries, fasta_file):
    '''Returns the RecSeq for a set of queries in a (REBASE) fasta file.'''
    seqs = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    description_queries = [seqs[query].description.split('\t') for query in queries]
    rec_seqs = [[x for x in description_query if 'RecSeq' in x][0] for description_query in description_queries]
    rec_seqs = [rec_seq.split(':')[1] for rec_seq in rec_seqs]
    return(rec_seqs)

def globalPercIdentity(seq_a, seq_b):
    '''Returns the optimal global percentage identity for two sequences.
    Args:
        seq_a, seq_b (str)
            The sequences to align

    Returns:
        global_pident (float)
            The % identity score (using the length of the shortest sequence)
    '''
    blosum62 = substitution_matrices.load('BLOSUM62')
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = blosum62
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -5
    aligner.mode = 'global'

    best_global_alignment = aligner.align(seq_a, seq_b)[0]
    matches = 0
    # Iterate over the aligned sequences and count matches
    for (a, b) in zip(best_global_alignment.aligned[0], best_global_alignment.aligned[1]):
        # In each tuple (a, b), 'a' and 'b' are (start, stop) pairs indicating the aligned regions
        start_a, stop_a = a
        start_b, stop_b = b
        matches += sum(1 for i in range(stop_a - start_a) if seq_a[start_a + i] == seq_b[start_b + i])
    min_length = min([len(seq_a), len(seq_b)]) # use shortest sequence length for pident calculation
    global_pident = matches/min_length * 100
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
    best_hits = best_hits.sort_values(['qseqid', 'hit_type'], ascending=[False, True]) # hit_type goes gold, nonputative, putative so ensures best hit is on top
    best_hits_collapse = pd.concat([best_hits[best_hits['qseqid']==x].head(1) for x in set(best_hits['qseqid'])])
    # CURRENTLY NOT DONE: Add alternative seq_ids (NOT STORED because it's redundant information)
    #best_hits_collapse = best_hits_collapse.assign(other_REBASE_hits = [','.join(list(best_hits[best_hits['qseqid']==x]['sseqid'])[1:]) for x in set(best_hits['qseqid'])])
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
            counter_dict[protein_id] = [entries[1], int(entries[2]), re.sub('location', '', str(entries[3])), entries[4]]
    return(counter_dict)

def predictRMS(hits_MT, hits_RE, with_position, position_threshold=5, mt_threshold=55, re_threshold=50):
    '''Predicts RM system based on tables of hits to MTases and REases.
    Args:
        hits_MT (DataFrame)
            Hits to MTases
        hits_RE (DataFrame)
            Hits to REases
        with_position (Bool)
            Whether to include positional information in making prediction.
        Default arguments: mt_threshold, re_threshold, position_threshold (int)
            Similarity thresholds required to rely on prediction of target sequence.
            Default are the values 55% (MTase) and 50% (REase) from Oliveira 2016. And position_threshold is <5.
    Returns:
        predicted_rms (list)
            Target sequences as keys with MT and RE proteins and positions stored as values
    '''
    # Filter hits based on Oliveira thresholds
    hits_MT = hits_MT[hits_MT['identity']>mt_threshold]
    hits_RE = hits_RE[hits_RE['identity']>re_threshold]
    # Add index
    hits_MT.index = hits_MT['qseqid']
    hits_RE.index = hits_RE['qseqid']
    # Subset to similarities for adding to dataframe


    # Check for any intersection of targets
    if with_position==True:
        # Subset to relevant columns (with genomic location - assumes all passed in correctly which is a somewhat dangerous assumption)
        hits_MT_subset = hits_MT[['qseqid', 'sseqid', 'identity', 'genomic_location']]
        hits_MT_subset.columns = ['qseqid', 'hit_MT', 'sim_MT', 'loc_MT']
        hits_RE_subset = hits_RE[['qseqid', 'sseqid', 'identity', 'genomic_location']]
        hits_RE_subset.columns = ['qseqid', 'hit_RE', 'sim_RE', 'loc_RE']

        target_overlap = set(hits_MT['target']).intersection(set(hits_RE['target']))
        if len(target_overlap) > 0:
            predicted_rms = []

            for t in target_overlap:
                MT_positions = list(hits_MT[hits_MT['target']==t]['position'])
                MT_contigs = list(hits_MT[hits_MT['target']==t]['contig'])
                RE_positions = list(hits_RE[hits_RE['target']==t]['position'])
                RE_contigs = list(hits_RE[hits_RE['target']==t]['contig'])
                MT_contig_descriptions = list(hits_MT[hits_MT['target']==t]['contig_description'])
                # Want all pairwise combinations of these
                separations = [(x, y, abs(x-y), MT_contigs[MT_positions.index(x)]==RE_contigs[RE_positions.index(y)], MT_contigs[MT_positions.index(x)], MT_contig_descriptions[MT_positions.index(x)]) for x in MT_positions for y in RE_positions] # List of tuples storing position of MT and RE and separation
                for s in separations:
                    if s[2]<position_threshold:
                        # Check if on same contig - think this can return errors
                        if s[3]==True:
                            rms_entry = [t, s[4], s[5], s[0], s[1],
                            list((hits_MT[hits_MT['position']==s[0]]['qseqid']))[0], list((hits_RE[hits_RE['position']==s[1]]['qseqid']))[0]]
                            if rms_entry[3]==rms_entry[4]: # If the hit is for the same protein, don't include it (only want paired RE and MT)
                                pass
                            else:
                                predicted_rms.append(rms_entry)
            if len(predicted_rms)!=0:
                rms_results = pd.DataFrame(predicted_rms, columns=['sequence', 'contig', 'contig_description', 'pos_MT', 'pos_RE', 'prot_MT', 'prot_RE'])
                # Add locations, similarity scores and best hit
                rms_results = rms_results.join(hits_MT_subset.set_index('qseqid'), on='prot_MT')
                rms_results = rms_results.join(hits_RE_subset.set_index('qseqid'), on='prot_RE')
                logging.info('  Predicted the following R-M systems:')
                for rms in predicted_rms:
                            logging.info(rms)
                return(rms_results)
            else:
                return(None)
        else:
            return(None)
    elif with_position==False:
        # Subset to relevant columns
        hits_MT_subset = hits_MT[['qseqid', 'sseqid', 'identity']]
        hits_MT_subset.columns = ['qseqid', 'hit_MT', 'sim_MT']
        hits_RE_subset = hits_RE[['qseqid', 'sseqid', 'identity']]
        hits_RE_subset.columns = ['qseqid', 'hit_RE', 'sim_RE']
        target_overlap = set(hits_MT['target']).intersection(set(hits_RE['target']))
        if len(target_overlap) > 0:
            predicted_rms = []
            for t in target_overlap:
                MT_matches = ','.join(hits_MT[hits_MT['target']==t]['qseqid'])
                RE_matches = ''.join(hits_RE[hits_RE['target']==t]['qseqid'])
                rms_entry = [t, 'unknown', 'unknown', 'unknown', 'unknown', MT_matches, RE_matches]
                predicted_rms.append(rms_entry)

            rms_results = pd.DataFrame(predicted_rms, columns=['sequence', 'contig', 'contig_description', 'pos_MT', 'pos_RE', 'prot_MT', 'prot_RE'])
            # Add similarity scores and best hit
            rms_results = rms_results.join(hits_MT_subset.set_index('qseqid'), on='prot_MT')
            rms_results = rms_results.join(hits_RE_subset.set_index('qseqid'), on='prot_RE')
            return(rms_results)

        else:
            return(None)


def readLookupDict(dict_file):
    '''Reads in lookup dictionary file and returns dict.'''
    # Read in MT lookup dict
    lookup_dict = {}
    for line in open(get_data(dict_file), 'r').readlines():
        line = line.strip('\n').split()
        lookup_dict[line[0]] = line[1]
    return lookup_dict

def searchMTasesTypeII(proteome_fasta, with_position=False, evalue_threshold=0.001, coverage_threshold=0.5, collapse=True, MTase_db='Type_II_MT_all.faa', MT_lookup='Type_II_MT_dict.txt', hmm=get_data('Type_II_MTases.hmm')):
    '''Searches for Type II MTases.
    Args:
        proteome_fasta (str)
            Fasta file of proteins.
        with_position (Bool)
            Whether positional info is included in the fasta or not, in headers of the form ">protein contig counter".
        evalue_threshold (float)
            Threshold to keep hits at. Default 0.001 as in Oliveira 2016
        coverage_threshold (float)
            Threshold of coverage. Default: 0.5 (i.e. 50%) as in Oliveira 2016
    Returns:
        blast_hits_collapse (DataFrame)
            DataFrame of best hits, collapsed to one row per protein
    '''
    MTase_fasta = get_data(MTase_db)
    MTase_blastdb = get_data('db/'+MTase_db)

    # Using Oliveira Type II MTase HMM profiles to search
    hmm_dict_MT = searchHMM(proteome_fasta, hmm)
    logging.info('  (hmm_raw) %d proteins matched MTases.' % len(hmm_dict_MT))

    # Filter hits
    hits_MT_filt = {k:v for k,v in hmm_dict_MT.items() if float(v[3])<evalue_threshold}
    logging.info('  (hmm_filtered) %d proteins matched MTases.' % len(hits_MT_filt))

    # Subset only the hits out from the proteome
    tmp_fasta = makeTmpFile(proteome_fasta,'_MT.faa')
    subsetFasta(proteome_fasta, list(hits_MT_filt.keys()), tmp_fasta)

    # Blast these hits against all Type II MTases to find best matches
    blast_hits_MT = blastpAgainstDB(tmp_fasta, MTase_blastdb, db_built=True)
    # Store the sequences for global alignment
    protein_seqs = SeqIO.to_dict(SeqIO.parse(tmp_fasta, 'fasta'))
    rebase_seqs = SeqIO.to_dict(SeqIO.parse(MTase_fasta, 'fasta'))
    # Remove tmp fasta file
    os.remove(tmp_fasta)


    # If no hits?
    if blast_hits_MT is None:
        # Log the blast hits
        logging.info('  (blast) 0 MTase-protein hits.')
        return
    else:
        # Log the blast hits
        logging.info('  (blast) %d MTase-protein hits.' % len(blast_hits_MT))
        # Filter coverage threshold
        blast_hits_MT = blast_hits_MT.assign(coverage_threshold_met=list(blast_hits_MT['length'] > coverage_threshold*blast_hits_MT['qlen'])) # Condition of 50% coverage as in Oliveira 2016. # TODO - decide if this additional threshold is needed
        logging.info('  (blast_cov_filtered) %d protein-REBASE hits.' % len(blast_hits_MT))
        if blast_hits_MT is None:
            return

        # Get the recognition sites of the best hits
        rs_MT = getRS(blast_hits_MT['sseqid'], MTase_fasta)
        blast_hits_MT = blast_hits_MT.assign(target=rs_MT)

        # Add genomic position if requested
        if with_position==True:
            counter_dict = parseCounterPreparedFasta(proteome_fasta)
            blast_hits_MT = blast_hits_MT.assign(contig=[counter_dict[x][0] for x in blast_hits_MT['qseqid']],
                                                position=[counter_dict[x][1] for x in blast_hits_MT['qseqid']],
                                                contig_description=[counter_dict[x][2] for x in blast_hits_MT['qseqid']],
                                                genomic_location=[re.sub('location\\=', '', counter_dict[x][3]) for x in blast_hits_MT['qseqid']])

        # Add the global similarity of the best hit. Need to have the sequences available
        blast_hits_MT['identity'] = blast_hits_MT.apply(lambda row : globalPercIdentity(str(protein_seqs[row['qseqid']].seq),
                     str(rebase_seqs[row['sseqid']].seq)), axis = 1)

        # Add the quality of the hit
        MT_lookup_dict = readLookupDict(MT_lookup)
        blast_hits_MT['hit_type'] = blast_hits_MT.apply(lambda row : MT_lookup_dict[row['sseqid']], axis=1)

        # Collapse the table to best hits
        if collapse==True:
            blast_hits_collapse = collapseBestHits(blast_hits_MT)
            logging.info('  (blast) %d putative MTases.' % len(blast_hits_collapse))
            return(blast_hits_collapse)
        else:
            return(blast_hits_MT)

def searchREasesTypeII(proteome_fasta, with_position=False, evalue_threshold=0.001, coverage_threshold=0.5, collapse=True, REase_db='protein_seqs_Type_II_REases.faa', RE_lookup='Type_II_RE_dict.txt', hmm=False, IIG=False):
    '''Searches a file of proteins against all known REases.
    Args:
        proteome_fasta (str)
            Fasta file with proteins
        with_position (Bool)
            Whether positional info is included in the fasta or not, in headers of the form ">protein contig counter".
        evalue_threshold (float)
            Threshold to filter blast hits at. Default: 0.001 as in Oliveira 2016
        coverage_threshold (float)
            Threshold of coverage. Default: 0.5 (i.e. 50%) as in Oliveira 2016
    Returns:
        blast_hits_collapse (DataFrame)
            DataFrame of best hits, one row per protein
    '''
    enzyme_name = 'REases'
    if IIG==True:
        enzyme_name = 'IIG RE/MTases'
    # Blasting for REases
    REase_fasta = get_data(REase_db)
    REase_blastdb = get_data('db/'+REase_db)

    if hmm is not False:
        hmm_dict_RE = searchHMM(proteome_fasta, hmm)
        logging.info('  (hmm_raw) %d proteins matched %s.' % (len(hmm_dict_RE), enzyme_name))

        # Filter hits
        hits_RE_filt = {k:v for k,v in hmm_dict_RE.items() if float(v[3])<evalue_threshold}
        logging.info('  (hmm_evalue_filtered) %d proteins matched %s.' % (len(hits_RE_filt), enzyme_name))
        # Subset only the hits out from the proteome
        tmp_fasta = makeTmpFile(proteome_fasta,'_RE.faa')
        subsetFasta(proteome_fasta, list(hits_RE_filt.keys()), tmp_fasta)

        # Blast these hits against all Type II REases to find best matches
        blast_hits_RE = blastpAgainstDB(tmp_fasta, REase_blastdb, db_built=True)
        # Store the sequences for global alignment
        protein_seqs = SeqIO.to_dict(SeqIO.parse(tmp_fasta, 'fasta'))
        rebase_seqs = SeqIO.to_dict(SeqIO.parse(REase_fasta, 'fasta'))
        # Remove tmp fasta file
        os.remove(tmp_fasta)

    else:
        blast_hits_RE = blastpAgainstDB(proteome_fasta, REase_blastdb, db_built=True)
        # Store the sequences for global alignment
        protein_seqs = SeqIO.to_dict(SeqIO.parse(proteome_fasta, 'fasta'))
        rebase_seqs = SeqIO.to_dict(SeqIO.parse(REase_fasta, 'fasta'))
    # Check if no hits!
    if blast_hits_RE is None:
        return(blast_hits_RE)
    # Filter out hits
    blast_hits_RE = blast_hits_RE.assign(coverage_threshold_met=list(blast_hits_RE['length'] > coverage_threshold*blast_hits_RE['qlen'])) # Condition of 50% coverage as in Oliveira 2016. # TODO - decide if this additional threshold is needed
    logging.info('  (blast_raw) %d protein-REBASE hits.' % len(blast_hits_RE))
    if blast_hits_RE is None:
        return(blast_hits_RE)
    blast_hits_RE_filt = blast_hits_RE[blast_hits_RE['coverage_threshold_met']==True]
    logging.info('  (blast_cov_filtered) %d protein-REBASE hits.' % len(blast_hits_RE_filt))

    # Add genomic position, if requested
    if with_position==True and len(blast_hits_RE_filt)>0:
        counter_dict = parseCounterPreparedFasta(proteome_fasta)
        blast_hits_RE_filt = blast_hits_RE_filt.assign(contig=[counter_dict[x][0] for x in blast_hits_RE_filt['qseqid']],
                                                position=[counter_dict[x][1] for x in blast_hits_RE_filt['qseqid']],
                                                contig_description=[counter_dict[x][2] for x in blast_hits_RE_filt['qseqid']],
                                                genomic_location=[re.sub('location\\=', '', counter_dict[x][3]) for x in blast_hits_RE_filt['qseqid']])

    # Add the recognition sequences
    blast_hits_RE = blast_hits_RE_filt.assign(target=getRS(blast_hits_RE_filt['sseqid'], REase_fasta))


    # Add the global similarity of the best hit. Need to have the sequences available
    if len(blast_hits_RE)>0:
        blast_hits_RE['identity'] = blast_hits_RE.apply(lambda row : globalPercIdentity(str(protein_seqs[row['qseqid']].seq),
                     str(rebase_seqs[row['sseqid']].seq)), axis = 1)

        # Add the quality of the hit
        RE_lookup_dict = readLookupDict(RE_lookup)
        blast_hits_RE['hit_type'] = blast_hits_RE.apply(lambda row : RE_lookup_dict[row['sseqid']], axis=1)

    # Collapse the table to best hits
    if collapse==True and len(blast_hits_RE)>0:
        blast_hits_collapse = collapseBestHits(blast_hits_RE)
        logging.info('  (blast_filtered) %d putative %s.' % (len(blast_hits_collapse), enzyme_name))

        return(blast_hits_collapse)
    else:
        return(blast_hits_RE)


def main():
    args = get_options()
    output = args.output
    mode = args.mode
    collapse_hits = not args.dontcollapse

    # Logger details
    level = logging.INFO
    format = '%(message)s'
    handlers = [logging.StreamHandler()]
    logging.basicConfig(level = level, format = format, handlers = handlers)
    logging.info('Started running rmsFinder.')

    if args.hmm=='oliveira':
        MTase_hmm = get_data('Type_II_MTases.hmm')
        logging.info('HMM: Oliveira.')
    elif args.hmm=='tesson':
        MTase_hmm = get_data('defense-finder/Type_II_MTases.hmm')
        logging.info('HMM: Tesson (defenseFinder)')
    else:
        MTase_hmm = get_data('Type_II_MTases.hmm')
        logging.info('HMM not recognised. Using HMM: Oliveira.')

    if args.db=='gold':
        logging.info('REBASE database: gold.')
        MT_db = 'Type_II_MT_gold.faa'
        RE_db = 'Type_II_RE_gold.faa'
        IIG_db = 'Type_IIG_gold.faa'
    elif args.db=='nonputative':
        logging.info('REBASE database: nonputative.')
        MT_db = 'Type_II_MT_nonputative.faa'
        RE_db = 'Type_II_RE_nonputative.faa'
        IIG_db = 'Type_IIG_nonputative.faa'
    elif args.db=='all':
        logging.info('REBASE database: all.')
        MT_db = 'Type_II_MT_all.faa'
        RE_db = 'Type_II_RE_all.faa'
        IIG_db = 'Type_IIG_all.faa'
    else:
        logging.info('ERROR: did not recognise db argument. Choose from: gold, nonputative, all')
        return

    if args.genbank==True:
        genbank_file = str(args.input[0])
        proteome_fasta = makeTmpFile(genbank_file,'faa')
        parseGenBank(genbank_file, proteome_fasta) # Make fasta file the way we like it
        include_position = True
    elif args.fasta==True:
        proteome_fasta = str(args.input[0])
        include_position = False
    elif args.panacotafasta==True:
        panacota_file = str(args.input[0])
        proteome_fasta = makeTmpFile(panacota_file,'faa')
        parsePanacota(panacota_file, proteome_fasta) # Make fasta file the way we like it
        include_position = True
    elif args.transcdsfasta==True:
        cds_fasta = str(args.input[0])
        proteome_fasta = makeTmpFile(cds_fasta,'faa')
        parseNCBITranslatedCDSFasta(cds_fasta, proteome_fasta)
        include_position = True

    # Check proteome fasta
    with open(proteome_fasta, 'r') as f:
        first_line = f.readline()
        if not first_line.startswith('>'):
            logging.info('ERROR: protein fasta file does not appear to be a fasta.')
            return
        else:
            file_characters = [x.upper() for x in list(f.readline().strip('\n'))]
            nt_characters = {'A', 'T', 'C', 'G', 'N', '-', 'X'} # Possible characters in nucleotide fasta
            sum_of_nt_characters = sum([file_characters.count(character) for character in nt_characters])
            if sum_of_nt_characters==len(file_characters):
                logging.info('ERROR:\tyou appear to have passed in a nucleotide fasta  (second line of file contains only {ATCGXN-}).\n\trmsFinder needs a protein fasta file as input.\n\tIf you want rmsFinder to process this file anyway then rerun with --forceprocessing.')
                if args.forceprocessing:
                    logging.info('WARNING: proceeding as requested with analysis anyway.')
                else:
                    return


    if 'RMS' in mode or 'MT' in mode: # Search for MTases
        logging.info('\nSearching for MTases...')
        MT_hits = searchMTasesTypeII(proteome_fasta, include_position, collapse=collapse_hits, MTase_db=MT_db, MT_lookup='Type_II_MT_dict.txt', hmm=MTase_hmm)
        if MT_hits is not None:
            MT_hits.to_csv(output+'_MT.csv', index=False, float_format="%.3f")
        else:
            pd.DataFrame(None).to_csv(output+'_MT.csv', index=False)
            logging.info('  No MTase hits.')
        logging.info('Finished searching for MTases.')

    if 'RMS' in mode or 'RE' in mode: # Search for REases
        logging.info('\nSearching for REases...')
        if args.hmm=='tesson':
            RE_hmm = get_data('defense-finder/Type_II_REases.hmm')
        else:
            RE_hmm = False
        RE_hits = searchREasesTypeII(proteome_fasta, with_position=include_position, collapse=collapse_hits, REase_db=RE_db, RE_lookup='Type_II_RE_dict.txt', hmm=RE_hmm)
        if RE_hits is not None:
            RE_hits.to_csv(output+'_RE.csv', index=False, float_format="%.3f")
        else:
            pd.DataFrame(None).to_csv(output+'_RE.csv', index=False)
            logging.info('  No REase hits.')
        logging.info('Finished searching for REases.')

    if 'RMS' in mode: # Predict RMS
        logging.info('\nPredicting RMS based on MTase and REase presence...')
        rms_predictions = None
        if (MT_hits is None):
            pd.DataFrame(None).to_csv(output+'_RMS.csv', index=False)
            logging.info('Predicted no Type II R-M systems.\n')
        elif (RE_hits is None):
            pd.DataFrame(None).to_csv(output+'_RMS.csv', index=False)
            logging.info('Predicted no Type II R-M systems.\n')
        else:
            if len(MT_hits)==0:
                pd.DataFrame(None).to_csv(output+'_RMS.csv', index=False)
                #logging.info('Predicted no Type II R-M systems.')
            elif len(RE_hits)==0:
                pd.DataFrame(None).to_csv(output+'_RMS.csv', index=False)
                #logging.info('Predicted no Type II R-M systems.')
            if len(MT_hits)>0 and len(RE_hits)>0: # check some hits exist
                rms_predictions = predictRMS(MT_hits, RE_hits, with_position=include_position)
            if rms_predictions is not None:
                logging.info('Predicted presence of %d Type II R-M systems.' % len(rms_predictions))
                if len(rms_predictions)>0:
                    rms_predictions.to_csv(output+'_RMS.csv', index=False, float_format="%.3f")
                else:
                    pd.DataFrame(None).to_csv(output+'_RMS.csv', index=False)
                    logging.info('Predicted no Type II R-M systems.\n')
            else:
                pd.DataFrame(None).to_csv(output+'_RMS.csv', index=False)
                logging.info('Predicted no Type II R-M systems.')

    if 'IIG' in mode:
        # Run in IIG Mode
        if args.hmm=='tesson':
            IIG_hmm = get_data('defense-finder/Type_IIG.hmm')
        else:
            IIG_hmm = False
        IIG_hits = searchREasesTypeII(proteome_fasta, with_position=include_position,
                                        collapse=collapse_hits, REase_db=IIG_db,
                                        RE_lookup='Type_IIG_dict.txt', hmm=IIG_hmm,
                                        IIG=True)
        if IIG_hits is not None:
            IIG_hits.to_csv(output+'_IIG.csv', index=False, float_format="%.3f")
        logging.info('Finished searching for IIG RE/MTases.')


    if os.path.exists(proteome_fasta) and args.fasta!=True:
        os.remove(proteome_fasta) # Remove the proteome fasta we made, if indeed we made it


if __name__ == "__main__":
    main()
