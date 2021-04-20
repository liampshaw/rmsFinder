#!/usr/bin/env python
# vim: set fileencoding=<utf-8> :
# Copyright 2021 Liam Shaw

import rebase_processing as rp
from Bio import SeqIO
import os
import subprocess

# Function to create MSA for each hit with the REBASE Gold standards within that family
def createMSA(input_fasta, reference_fasta, hmm_file, family):
    '''Creates MSA for hit against REBASE sequences from the same family.

    Args:
        protein (dict item)
            A protein with its sequence
        '''
    # get the rebase gold hits with that family out of the fasta file
    #rebase_hits = extractHitsGold(rebase_gold_fasta, hmm, family)
    # Load the REBASE hits (prepREBASE)
    if os.path.isfile(reference_fasta):
        reference = SeqIO.to_dict(SeqIO.parse(reference_fasta, 'fasta'))
        input = SeqIO.to_dict(SeqIO.parse(input_fasta, 'fasta'))
        # write these to file
        tmp_fasta = 'tmp.fa'
        with open(tmp_fasta, 'w') as f:
            for r in reference:
                f.write('>%s\n%s\n' % (r, str(reference[r].seq)))
            for i in input:
                f.write('>%s\n%s\n' % (i, str(input[i].seq)))
        ##### Fetch just the family of interest and write to tmp.hmm
        hmm_fetch_command = ['hmmfetch', '-o', 'tmp.hmm', hmm_file, family]
        hmm_fetch_process = subprocess.Popen(hmm_fetch_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE) # Runs the process
        hmm_fetch_out, _ = hmm_fetch_process.communicate() # Read the output from stdout
        #### Compute the MSA
        print('creating MSA...')
        hmmalign_command = ['hmmalign',  '-o', 'tmp.aln', 'tmp.hmm', 'tmp.fa'] # align, output to tmp.aln
        hmmalign_process = subprocess.Popen(hmmalign_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE) # Runs the process
        hmmalign_out, _ = hmmalign_process.communicate() # Read the output from stdout
        #os.remove(tmp_fasta)
        return(1)
    else:
        print('That family has zero matches in REBASE Gold!')
        return(0)

def getPID(a, b):
    '''% identity for two protein strings a and b.'''
    n_identical = len([x for x in zip(a, b) if x[0]==x[1] and x[0]!='-'])
    min_length = min(len([x for x in a if x!='-']),
                    len([x for x in b if x!='-']))
    pid = float(n_identical)/float(min_length) * 100.0
    return(pid)

def getTopMatch(query, msa_file):
    '''Returns top match (by pairwise similarity) for a protein in a MSA.'''
    msa = SeqIO.to_dict(SeqIO.parse(msa_file, 'stockholm'))
    query_seq = str(msa[query].seq)
    hit_dict = dict()
    for seq in msa.keys():
        if seq!=query:
            hit_dict[seq] = getPID(str(msa[seq].seq), query_seq)
    return(max(hit_dict.values()))


def calculatePID(protein, msa_file):
    '''Returns the pids (% ids) for a protein in an MSA.'''
    print('calculating pid...')
    esl_alipid_command = ['esl-alipid', msa_file] # gets all pids
    esl_alipid_process = subprocess.Popen(esl_alipid_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE) # Runs the process
    output = subprocess.check_output(('grep', protein), stdin=esl_alipid_process.stdout) # Find just pairwise comparisons involving protein
    output_split = output.decode().split('\n') # split output
    output_df = pd.DataFrame(data=[row.split() for row in output_split]) # make dataframe
    output_df = output_df.dropna() # remove None/NA
    output_df[2] = pd.to_numeric(output_df[2]) # convert the pid (col 2) to numeric
    output_df = output_df.sort_values(2, ascending=False) # sort by decreasing pid
    return(output_df)

def getMSAs(proteome_fasta, hits, hmm, hmm_file, msa_dir):
    '''Gets all MSAs for a proteome dict.'''
    proteome = SeqIO.to_dict(SeqIO.parse(proteome_fasta, 'fasta'))
    # Take the unique families
    families = list(set([h[1] for h in hits.values()]))
    print(families)
    for fam in families:
        tmp_fasta = 'tmp_'+fam+'.fa' # Open tmp fasta and write sequences to it
        with open(tmp_fasta, 'w') as f:
            for prot, hit in hits.items():
                if hit[1]==fam:
                    f.write('>%s\n%s\n' % (prot, str(proteome[prot].seq)))
        createMSA(tmp_fasta, msa_dir+'/'+hmm+'.'+fam+'.fa',
                hmm_file, fam)


    for prot in hits.keys():
        fam = hits[prot][1] # family
        print('')
        print(prot)
        print(fam)
        msa = createMSA(proteome[prot], str(proteome[prot].seq),
                        msa_dir+'/'+hmm+'.'+fam+'.fa', hmm_file, fam)
        if msa==0:
            pass
        else:
            output_pid = calculatePID(prot, 'tmp.aln')
            print(output_pid.iloc[[0]])
            if 'stored_results' in locals():
                output_pid['hmm'] = hmm
                output_pid['family'] = fam
                stored_results = stored_results.append(output_pid.iloc[[0]])
                #stored_results['family'][nrow()]
            else:
                stored_results = output_pid.iloc[[0]]
                stored_results['hmm'] = hmm
                stored_results['family'] = fam
    return stored_results
