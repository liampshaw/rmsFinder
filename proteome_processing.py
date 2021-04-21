#!/usr/bin/env python
# vim: set fileencoding=<utf-8> :
# Copyright 2021 Liam Shaw

import rebase_processing as rp
from Bio import SeqIO
import os
import subprocess

# Function to create MSA for each hit with the REBASE Gold standards within that family
def createMSA(input_fasta, reference_fasta, hmm_file, family, msa_output):
    '''Creates MSA for hit against REBASE sequences from the same family.

    Args:
        input_fasta (str)
            Fasta file with protein sequences to align.
        reference_fasta (str)
            Fasta file with reference sequences.
        hmm_file (str)
            File with hmm profile.
        family (str)
            Name of family (within hmm profile).
        msa_output (str)
            Filename to write MSA to.

    Returns:
        1 (success) or 0 (failure)
    '''
    # Load the REBASE hits (typically prepared with prepREBASE)
    if os.path.isfile(reference_fasta):
        reference = SeqIO.to_dict(SeqIO.parse(reference_fasta, 'fasta'))
        input = SeqIO.to_dict(SeqIO.parse(input_fasta, 'fasta'))
        # write these to file
        tmp_fasta = 'tmp.with.reference.'+input_fasta
        tmp_hmm = 'tmp.hmm'
        with open(tmp_fasta, 'w') as f:
            for r in reference:
                f.write('>%s\n%s\n' % (r, str(reference[r].seq)))
            for i in input:
                f.write('>%s\n%s\n' % (i, str(input[i].seq)))
        ##### Fetch just the family of interest and write to tmp.hmm
        hmm_fetch_command = ['hmmfetch', '-o', tmp_hmm, hmm_file, family]
        hmm_fetch_process = subprocess.Popen(hmm_fetch_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE) # Runs the process
        hmm_fetch_out, _ = hmm_fetch_process.communicate() # Read the output from stdout
        #### Compute the MSA
        print('creating MSA...')
        hmmalign_command = ['hmmalign',  '-o', msa_output, tmp_hmm, tmp_fasta] # align, output to tmp.aln
        hmmalign_process = subprocess.Popen(hmmalign_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE) # Runs the process
        hmmalign_out, _ = hmmalign_process.communicate() # Read the output from stdout
        os.remove(tmp_fasta)
        os.remove(tmp_hmm)
        return(1)
    else:
        print('That family has zero matches in the reference fasta provided!')
        return(0)

def getPID(a, b):
    '''% identity for two protein strings a and b.

    Args:
        a, b (str)
            Both protein strings (although could be any character string).

    Returns:
        pid (float)
            Percentage identity of a,b.
    '''
    n_identical = len([x for x in zip(a, b) if x[0]==x[1] and x[0]!='-'])
    min_length = min(len([x for x in a if x!='-']),
                    len([x for x in b if x!='-']))
    pid = float(n_identical)/float(min_length) * 100.0
    return(pid)

def getTopMatch(query, msa_file, msa_format='stockholm', match='REBASE'):
    '''Returns top match (by pairwise similarity) for a protein in a MSA.

    Args:
        query (str)
            Name of protein.
        msa_file (str)
            Filename of MSA (which query must be within).
        (opt) msa_format (str)
            Format of MSA. Default is Stockholm since this is hmmer's default.
        (opt) match (str)
            String which starts legitimate matches.

    Returns:
        top_hits (dict)
            Dict of top hit(s).
            (could be tweaked)
    '''
    print('Finding match for', query, 'in', msa_file)
    msa = SeqIO.to_dict(SeqIO.parse(msa_file, format=msa_format))
    query_seq = str(msa[query].seq)
    hit_dict = dict()
    for seq in msa.keys():
        if seq!=query:
            hit_dict[seq] = getPID(str(msa[seq].seq), query_seq)
    hit_dict_2 = {k:v for k, v in hit_dict.items() if k.startswith(match)}
    max_value = max(hit_dict_2.values())
    top_hits = {k:v for k, v in hit_dict_2.items() if v==max_value}
    return(top_hits)


def hmmer2AlignToReference(input_fasta, reference_msa, hmm_profile, output_aln):
    '''Aligns an input fasta to a reference MSA with hmmer2.

    Args:
        input_fasta (str)
            Filename of fasta file to align.
        reference_msa (str)
            Filename of (pre-computed) alignment
        hmm_profile (str)
            Filename of hmm profile
        output_aln (str)
            Filename of output

    Returns:
        None
    '''
    # To be completed.
    hmm2align_command = ['/Users/liam/Applications/hmmer-2.3.2/src/hmmalign',
                        '--withali', reference_msa,
                        '-o', output_aln,
                        hmm_profile,
                        input_fasta]
    hmm2align_process = subprocess.Popen(hmm2align_command,
                        stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    hmm2align_out, _ = hmm2align_process.communicate()
    #print(hmm2align_out)
    return

def getHits(proteome_fasta, hits, hmm, hmm_file, msa_dir):
    '''Returns the top hits for each protein.

    Args:
        proteome_fasta (str)
            Filename of fasta file with proteins in.
        hits (dict)
            The best family match for each protein from hmmsearch.
        hmm (str)
            The name of the hmm profile.
        hmm_file (str)
            Filename of the hmm profile.
        msa_dir (str)
            Directory to save output to.

    Returns:
        results_dict (dict)
            Dict of hits (themselves dicts)
    '''
    proteome = SeqIO.to_dict(SeqIO.parse(proteome_fasta, 'fasta'))
    # Take the unique families
    families = list(set([h[1] for h in hits.values()]))
    print(families)
    results_dict = {}
    for fam in families:
        tmp_fasta = 'tmp_'+fam+'.fa' # Open tmp fasta and write sequences to it
        proteins_of_interest = []
        with open(tmp_fasta, 'w') as f:
            for prot, hit in hits.items():
                if hit[1]==fam:
                    proteins_of_interest.append(prot)
                    f.write('>%s\n%s\n' % (prot, str(proteome[prot].seq)))
        tmp_aln = 'tmp_'+fam+'.aln'
        reference_aln = msa_dir+'/'+hmm+'.'+fam+'.fa.aln'
        reference_hmm = msa_dir+'/'+hmm+'.'+fam+'.hmm.2'
        #print(reference_aln)
        if os.path.isfile(reference_aln):
            hmmer2AlignToReference(tmp_fasta, reference_aln,
                            reference_hmm, tmp_aln)
        else:
            print('No reference alignment available...!')
            break
        for prot in proteins_of_interest:
            print()
            print(prot)
            match = getTopMatch(prot, tmp_aln)
            print(match)
            results_dict[prot] = [fam, match]
    return results_dict

def getRS(query, fasta_file):
    '''Returns the RecSeq for a query in a (REBASE) fasta file.'''
    seqs = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    description_query = seqs[query].description.split('\t')
    rec_seq = [x for x in description_query if 'RecSeq' in x][0]
    rec_seq = rec_seq.split(':')[1]
    return(rec_seq)
