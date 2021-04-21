#!/usr/bin/env python
# vim: set fileencoding=<utf-8> :
# Copyright 2021 Liam Shaw

import subprocess
import os
from Bio import SeqIO
import re


def extractSequences(rebase_proteins, enzyme_type, output_file):
    '''Extracts sequences with known RS and enzyme type.

    Args:
        rebase_proteins (str)
            Filename of input fasta (downloaded from REBASE: protein_seqs.txt)
        enzyme_type (str)
            Name of desired enzyme type to extract e.g. Type I methytransferase
        output_file (str)
            Filename of output fasta file

    Returns:
        None
    '''
    RecSeqFlag = False # Bool to track if sequence has known sequence
    i = 0; seqs = 0 # Counters
    with open(output_file, 'w') as output:
        with open(rebase_proteins, 'r') as f:
            for line in f.readlines():
                line = line.strip('\n')
                line = re.sub('<>', '', line)
                if line.startswith('>'):
                    i += 1
                    if 'RecSeq' in line and enzyme_type in line:
                        #if want only non-putative: 'and 'putative' not in line
                        seqs += 1
                        RecSeqFlag = True
                        # Check if bracket in line and add underscore
                        # to avoid duplicate lines
                        if '(' in line:
                            line = re.sub(' \\(', '_(', line)
                        output.write('%s\n' % line)
                    else:
                        RecSeqFlag = False
                else:
                    if line.startswith('<>'):
                        pass
                    elif line=='':
                        pass
                    elif RecSeqFlag==True:
                        output.write('%s\n' % line)
        print(i, seqs)
    return


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
                    if float(line[4]) < float(df[line[0]][3]):
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


def hmmer2Align(input_fasta, hmm_profile, output_aln):
    '''Aligns a fasta and hmm profile (originally in hmmer3) with hmmer2.

    Args:
        input_fasta (str)
            Filename of input fasta
        hmm_profile (str)
            Filename of input hmm profile (hmmer3 format)
        output_aln (str)
            Filename of output file

    Returns:
        None
    '''
    # Convert the profile to hmmer2 format.
    hmm2_profile = hmm_profile+'.2'
    hmm2_profile_file = open(hmm2_profile, 'w')
    #hmmoutput = subprocess.run(['ls', '-l'], stdout=PIPE).stdout.splitlines()
    hmmconvert_command = ['hmmconvert', '-2', hmm_profile]
    hmmconvert_process = subprocess.Popen(hmmconvert_command,
                        stdout = hmm2_profile_file)
    hmmconvert_out, _ = hmmconvert_process.communicate()
    hmmconvert_process.wait()
    hmm2_profile_file.close()

    # Align with hmmer2.
    hmm2align_command = ['/Users/liam/Applications/hmmer-2.3.2/src/hmmalign',
                        '-o', output_aln,
                        hmm2_profile,
                        input_fasta]
    hmm2align_process = subprocess.Popen(hmm2align_command,
                        stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    hmm2align_out, _ = hmm2align_process.communicate()
    hmm2align_process.wait()
    return

def hmmer3Align(input_fasta, hmm_profile, output_aln):
    '''Aligns a fasta and hmm profile with hmmer3.

    Args:
        input_fasta (str)
            Filename of input fasta
        hmm_profile (str)
            Filename of input hmm profile (hmmer3 format)
        output_aln (str)
            Filename of output file

    Returns:
        None
    '''
    # Align with hmmer3.
    hmmalign_command = ['hmmalign',
                        '-o', output_aln,
                        hmm_profile,
                        input_fasta]
    hmmalign_process = subprocess.Popen(hmmalign_command,
                        stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    hmmalign_out, _ = hmmalign_process.communicate()
    hmmalign_process.wait()
    return

def fetchHMMs(hmm_file, out_dir):
    '''Fetches all HMMs out of a profile.'''
    fam_file_str = 'tmp.fams'
    with open(fam_file_str, 'w') as fam_file:
        hmmstat_command = ['hmmstat', hmm_file]
        hmmstat_process = subprocess.Popen(hmmstat_command,
                            stdout = fam_file)
        hmmstat_process.wait()
    families = []
    with open(fam_file_str, 'r') as fam_file:
        for line in fam_file.readlines():
            if line.startswith('#'):
                pass
            else:
                families.append(line.split()[1]) # append family
    for fam in families:
        fam_hmm_file = out_dir+'/'+re.sub('.*\/', '', re.sub('hmm', '', hmm_file))+fam+'.hmm'
        hmmfetch_command = ['hmmfetch',
                            '-o', fam_hmm_file,
                            hmm_file, fam]
        hmmfetch_process = subprocess.Popen(hmmfetch_command,
                            stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        hmmfetch_out, _ = hmmfetch_process.communicate()
        # Convert to hmm2
        hmm2_profile = fam_hmm_file+'.2'
        hmm2_profile_file = open(hmm2_profile, 'w')
        #hmmoutput = subprocess.run(['ls', '-l'], stdout=PIPE).stdout.splitlines()
        hmmconvert_command = ['hmmconvert', '-2', fam_hmm_file]
        hmmconvert_process = subprocess.Popen(hmmconvert_command,
                            stdout = hmm2_profile_file)
        hmmconvert_out, _ = hmmconvert_process.communicate()
        hmmconvert_process.wait()
        hmm2_profile_file.close()
        print('Fetched profile for', fam, 'and wrote to', fam_hmm_file, 'and', hmm2_profile)
    return

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
    if not os.path.exists(dir):
        os.makedirs(dir)

    # Fetch HMMs out of main profile
    fetchHMMs(hmm_profile, dir)
    # Search HMM
    hmm_search_results = searchHMM(input_fasta, hmm_profile)
    families = list(set([hmm_search_results[k][1]
                        for k in hmm_search_results.keys()]))
    for fam in families:
        print(fam)
        # Fetch the family HMM out
        # Fetch the sequences out for that family
        rebase_subset = extractHitsFasta(input_fasta, hmm_search_results,
                                        family=fam)
        subset_file_str = dir+'/'+hmm+'.'+fam+'.fa'
        with open(subset_file_str, 'w') as f:
            for id in rebase_subset:
                f.write('>%s\n%s\n' % (id, str(rebase_subset[id].seq)))
        # Then align these sequences
        family_hmm = dir+'/'+hmm+'.'+fam+'.hmm'
        print(family_hmm)
        hmmer2Align(subset_file_str, family_hmm, subset_file_str+'.aln')
    return
