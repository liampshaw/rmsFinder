#!/usr/bin/env python
# vim: set fileencoding=<utf-8> :
# Copyright 2021 Liam Shaw

import subprocess
import os
from Bio import SeqIO


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
    tmp_file = 'tmp.out'
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
                    if line[4] < df[line[0]][3]:
                        df[line[0]] = line[1:]
                    else: # otherwise no (keep top hit)
                        pass
                else: # if not in dict, add the hit
                    df[line[0]] = line[1:]
    os.remove(tmp_file)
    return df


def extractHitsFasta(input_fasta, stored_hmm_search_results, family='all'):
    '''Extracts sequences from REBASE fasta which have a hit to a specified
    HMM profile.
    Sped up by computing the hmmsearch once (with prepREBASE) beforehand.

    Args:
        rebase_fasta (str)
            Filename of input fasta
        stored_hmm_search_results (dict)
            Dict of hmmsearch results (from searchHMM)
        family (str)
            Optional string specifying key of family within HMM profile.
            Default is to process all families.

    Returns:
        proteome_seqs_subset (dict)
            Dict of subsetted results
    '''
    results = stored_hmm_search_results
    if family!='all':
        results_subset = dict()
        for r in stored_hmm_search_results.keys():
            if stored_hmm_search_results[r][1]==family:
                results_subset[r] = stored_hmm_search_results[r]
            else:
                pass
        results = results_subset
    proteome_seqs = SeqIO.to_dict(SeqIO.parse(input_fasta, 'fasta'))
    proteome_seqs_subset = {k:v for k,v in proteome_seqs.items()
                            if k in results.keys()} # subset the fasta
    return proteome_seqs_subset


def separateHMMs(hmm_profile):
    '''Fetches all profiles from a HMM profile.'''
    hmmsearch_command = ['hmmsearch', '--noali', '--tblout', tmp_file,
                        hmm_file, query_protein_file]
    hmmsearch_process = subprocess.Popen(hmmsearch_command,
                        stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    # Not done


def prepREBASE(input_fasta, hmm_profile, hmm, dir):
    '''Preps a fasta for a given HMM profile, matching each protein in the
    input to its top family hit and writing them to a separate file.

    Args:
        rebase_fasta (str)
            Filename of input fasta (presumed from REBASE)
        hmm_profile (str)
            Filename of HMM profile
        hmm (str)
            Name of HMM (for output)
        dir (str)
            Path of directory to write output

    Returns:
        None
    '''
    hmm_search_results = searchHMM(input_fasta, hmm_profile)
    families = list(set([hmm_search_results[k][1] for k in hmm_search_results.keys()]))
    for fam in families:
        # Fetch the family HMM out
        # Fetch the sequences out for that family
        rebase_subset = extractHitsFasta(input_fasta, hmm_search_results, family=fam)
        print(rebase_subset)
        with open(dir+'/'+hmm+'.'+fam+'.fa', 'w') as f:
            for id in rebase_subset:
                f.write('>%s\n%s\n' % (id, str(rebase_subset[id].seq)))
    return
