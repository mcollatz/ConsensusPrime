#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 09:20:45 2021

@author: Maximilian Collatz
"""


import argparse
import pandas as pd
import os
import subprocess
import re
import sys

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile", type=str, help="Input nucletide alignment .fna.")
parser.add_argument("-o", "--outdir", type=str, help="Name of the output directory.", default='')
# parser.add_argument("-d", "--delim", type=str, help="Column delimiter of the fasta header. Default '\t'",default='\t')
# parser.add_argument("-p", "--position", type=int, help="Position of geneID in fasta header. Default = 0", default=0)
parser.add_argument("-t", "--threads", type=str, help="Number of threads used by MAFFT. Default = -1 (all)", default=-1)
parser.add_argument("-k", "--keepduplicates", type=str2bool, help="Keep duplicate sequences. Default = False", default=False)
parser.add_argument("-c", "--consensusthreshold", type=float, help="Consensus threshold bitween 0 and 1 with 1 beeing a perfect consensus for regions to be considered for primer prediction in the final alignment. Default = 0.95", default=0.95)
#parser.add_argument("-g", "--gapthreshold", type=float, help="Percentage of gaps in sequences to be considered as partial sequence for removal. Default = 0.2", default=0.2)
parser.add_argument("-s", "--consensussimilarity", type=float, help="Minimum similarity threshold for sequences in the input alignment when comparing each sequence to the consensus sequence. Default = 0.8", default=0.8)
parser.add_argument("-x", "--primer3", type=str, help="Primer3 input parameter file.")
parser.add_argument("--primers", type=str, help="Known primers for visualization in primer plot in .fasta format.", default="")
parser.add_argument("--negativesequences", type=str, help="File with sequences that get their consensus sequence added to the final alignment .fna.", default='')
args = parser.parse_args()

# splitting the fasta header is not really needed
# delim = args.delim
delim = '\t'
# idpos = args.position
idpos = 0
alignment = args.infile
threshold = args.consensusthreshold
#gapthreshold = args.gapthreshold
consensussimilarity = args.consensussimilarity
keepduplicates = args.keepduplicates
primer3file = args.primer3
primers = args.primers
negativesequencefile = args.negativesequences
threads = str(args.threads)

if args.outdir:
	outdir = args.outdir
	if not os.path.isabs(outdir):
		outdir = f'{os.getcwd()}/{args.outdir}'
	if not os.path.exists(outdir):
		os.makedirs(outdir)
else:
	outdir = os.getcwd()

if not os.path.exists(outdir + '/results'):
	os.makedirs(outdir + '/results')
outdir = outdir + '/results'
print(f'\nOut dir set to : {outdir}\n')

if not os.path.exists(outdir + '/primer3_runfiles'):
	os.makedirs(outdir + '/primer3_runfiles')
primer3_runfiles_dir = outdir + '/primer3_runfiles'

if not os.path.exists(outdir + '/primer3_results'):
	os.makedirs(outdir + '/primer3_results')
primer3_result_dir = outdir + '/primer3_results'


def read_fasta(multifasta, delim, idpos):
    """reading input fasta file
    returns fasta dictionary with key=accessionNumber and value=Sequence
    returns fastaheader dictionary with key=accesstionNumber and value=originalFastaHeader"""
    fasta = {}
    fastaheader = {}
    with open(multifasta, 'r') as infile:
        acNumber = ''
        for line in infile:
            if line.startswith('>'):
                if delim:
                    acNumber = line.split(delim)[idpos].strip().strip('>')
                    fastaheader[acNumber] = line.strip()
                else:
                    acNumber = line.split()[idpos].strip().strip('>')
                    fastaheader[acNumber] = line.strip()
            else:
                if acNumber in fasta:
                    fasta[acNumber] += line.strip().upper()
                else:
                    fasta[acNumber] = line.strip().upper()
    return fasta, fastaheader

def avg(liste):
    return sum(liste)/len(liste)

# function to add line breaks to sequence output for nicely formatted alignment files
def insert_newlines(string, linelen=64): # linelen= 64 is recommended for a well readable file.
	return '\n'.join(string[i:i+linelen] for i in range(0, len(string), linelen))

# function to build mafft alignments
def align(infile,outfile):
    with open(outfile, 'w') as outfile:
        subprocess.call(['mafft','--auto' , '--adjustdirection', '--thread', threads, infile], stdout=outfile, stderr=subprocess.DEVNULL)

# reverse complement
def revcomp(seq):
    seq = seq.upper().replace("-","")
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return("".join(complement.get(base, base) for base in reversed(seq)))

# create html table from array of arrays [[row1,element2],[row2,element2]] with optional header row [colheader1,colheader2]
def to_html_table(table,header=[]):
    html_table = '<table border=1>\n'
    if header:
        html_table += '<tr><th>' + '</td><th>'.join([str(x) for x in header]) + '</th></tr>\n'
    for row in table:
        html_table += '<tr><td>' + '</td><td>'.join([str(x) for x in row]) + '</td></tr>\n'
    html_table += '</table>\n'
    return(html_table)

# returns the percentage similarity of two sequency as value betwenn 0 and 1
def compare_similarity(seq1, seq2):
    seqlen = len(seq1)
    mismatch_count = len([pos for pos in range(seqlen) if seq1[pos] != seq2[pos] ])
    return(1 - mismatch_count/seqlen)


# create table for Number of Sequences
alignment_statistic_table = []

# read input fasta/alignment
fasta, fastaheader = read_fasta(alignment, delim, idpos)
number_of_sequences = len(fasta)
alignment_statistic_table.append(['Input sequences', number_of_sequences])
print(f'Number of sequences: {number_of_sequences}')

# make normal fasta from alignment if alignment was submitted
for gene in fasta:
    fasta[gene] = fasta[gene].replace('-','')

# build base alignment
print(f'Building and writing base alignment to: {outdir}/base_alignment.fna')
outfile = f'{outdir}/base_fasta.fna'
with open(outfile, 'w') as outfile:
    for header in fasta:
        outfile.write(f'>{header}\n')
        outfile.write(f'{insert_newlines(fasta[header])}\n')

if len(fasta) > 1:
    with open(f'{outdir}/base_alignment.fna', 'w') as outfile:
            subprocess.call(['mafft','--auto' , '--adjustdirection', '--reorder', '--thread', threads, f'{outdir}/base_fasta.fna'], stdout=outfile, stderr=subprocess.DEVNULL)
else:
    with open(f'{outdir}/base_alignment.fna', 'w') as outfile:
            subprocess.call(['mafft','--auto' , '--adjustdirection', '--thread', threads, f'{outdir}/base_fasta.fna'], stdout=outfile, stderr=subprocess.DEVNULL)

# read base_alignment fasta
fasta, fastaheader = read_fasta(f'{outdir}/base_alignment.fna', delim, idpos)



# remove duplicates from fasta and build new alignment
unique = ''
if not keepduplicates:
    unique = 'unique_'
    outfile = f'{outdir}/unique_fasta.fna'
    sequences = set()
    number_of_uniques = 0
    with open(outfile,'w') as outfile:
        for header in fasta:
            seqeuence_gapless = fasta[header].replace("-","")
            if not seqeuence_gapless in sequences:
                outfile.write(f'>{header}\n')
                outfile.write(f'{insert_newlines(fasta[header])}\n')
                sequences.add(seqeuence_gapless)
                number_of_uniques += 1
    print(f'Building and writing unique alignment to: {outdir}/unique_alignment.fna\n')
    if number_of_uniques > 1:
        align(infile=f'{outdir}/unique_fasta.fna', outfile=f'{outdir}/unique_alignment.fna')
    else:
        command = f'cp {outdir}/unique_fasta.fna {outdir}/unique_alignment.fna'
        os.system(command)
    
    # read unique_alignment fasta
    fasta, fastaheader = read_fasta(f'{outdir}/unique_alignment.fna', delim, idpos)
    number_of_sequences = len(fasta)
    alignment_statistic_table.append(['Unique sequences', number_of_sequences])
    print(f'Number of unique sequences: {number_of_sequences}')


####### filter by consensussimilarity #######

# convert fasta to panda dataframe to calculate consensus scores
fasta_for_panda = {}
for gene in fasta:
    fasta_for_panda[gene] = list(fasta[gene])
fasta_pd = pd.DataFrame.from_dict(fasta_for_panda, orient='index')
# calculate nucleotide ratio for every position in the alighment
# NOTE: non-standard nucleotides are counted as '-'
consensus_dict = {}
for col in fasta_pd:
    nucl_score = {'A':0,'T':0,'C':0,'G':0,'-':0}
    for nt in fasta_pd[col]:
        if nt in nucl_score:
            nucl_score[nt] += 1
        else:
            nucl_score['-'] += 1
    for key in nucl_score:
        nucl_score[key] = nucl_score[key] / number_of_sequences
    consensus_dict[col] = nucl_score

# get consensus sequence
# The score is the proportion of the most common amino acid in a column of the alignment.
consensus_seq = ''
consensus_score = []
for pos in consensus_dict:
    consensus_nt = ''
    max_consensus_score = 0
    for nt in consensus_dict[pos]:
        score = consensus_dict[pos][nt]
        if score > max_consensus_score:
            consensus_nt = nt
            max_consensus_score = score
    consensus_seq += consensus_nt
    consensus_score.append(max_consensus_score)
    max_consensus_score = 0

# filter unsimilar sequences
similarity_filtered_dict = {}
for header, sequence in fasta.items():
    if compare_similarity(sequence, consensus_seq) >= consensussimilarity:
        similarity_filtered_dict[header] = sequence

outfile = f'{outdir}/{unique}similarity_filtered_fasta.fna'
with open(outfile,'w') as outfile:
    for header in similarity_filtered_dict:
        outfile.write(f'>{header}\n')
        outfile.write(f'{insert_newlines(similarity_filtered_dict[header])}\n')
# build new alignment without unsimilar sequences
print(f'Building and writing similarity filtered alignment to: {outdir}/{unique}similarity_filtered_alignment.fna\n')
if len(similarity_filtered_dict) > 1:
    align(infile=f'{outdir}/{unique}similarity_filtered_fasta.fna', outfile=f'{outdir}/{unique}similarity_filtered_alignment.fna')
else:
    os.system(f'cp {outdir}/{unique}similarity_filtered_fasta.fna {outdir}/{unique}similarity_filtered_alignment.fna')

# read similarity_filtered_alignment fasta
fasta, fastaheader = read_fasta(f'{outdir}/{unique}similarity_filtered_alignment.fna', delim, idpos)
number_of_sequences = len(fasta)
print(f'Number of {unique}similar sequences: {number_of_sequences}\n')
if keepduplicates:
    alignment_statistic_table.append(['Similar sequences', number_of_sequences])
else:
    alignment_statistic_table.append(['Unique similar sequences', number_of_sequences])


# convert fasta to panda dataframe to calculate consensus scores
fasta_for_panda = {}
for gene in fasta:
    fasta_for_panda[gene] = list(fasta[gene])
fasta_pd = pd.DataFrame.from_dict(fasta_for_panda, orient='index')

# calculate nucleotide ratio for every position in the alighment
# NOTE: non-standard nucleotides are counted as '-'
consensus_dict = {}
for col in fasta_pd:
    nucl_score = {'A':0,'T':0,'C':0,'G':0,'-':0}
    for nt in fasta_pd[col]:
        if nt in nucl_score:
            nucl_score[nt] += 1
        else:
            nucl_score['-'] += 1
    for key in nucl_score:
        nucl_score[key] = nucl_score[key] / number_of_sequences
    consensus_dict[col] = nucl_score

# get consensus sequence
# The score is the proportion of the most common amino acid in a column of the alignment.
consensus_seq = ''
consensus_score = []
for pos in consensus_dict:
    consensus_nt = ''
    max_consensus_score = 0
    for nt in consensus_dict[pos]:
        score = consensus_dict[pos][nt]
        if score > max_consensus_score:
            consensus_nt = nt
            max_consensus_score = score
    consensus_seq += consensus_nt
    consensus_score.append(max_consensus_score)
    max_consensus_score = 0



# remove gaps from the consensus sequence as primers can not be designed with gaps
gap_positions = []
for i in range(len(consensus_seq)-1,-1,-1):
    if '-' == consensus_seq[i]:
        gap_positions.append(i)

consensus_seq_gapless = list(consensus_seq)
for pos in gap_positions:
    del consensus_seq_gapless[pos]
    del consensus_score[pos]
consensus_seq_gapless = ''.join(consensus_seq_gapless)



# get regions above the consensus threshold (ratio of most common amino acid per posisiton)
# a consensus region must be at least of the length 20 to be considered for primer design
print(f'Consensus threshold set to: {threshold}\n')
minlength = 20
start = 0
stop = 0
new = True
old_start = -1
consensus_regions = []
for pos in range(len(consensus_score)):
    if consensus_score[pos] >= threshold:
        if new:
            start = pos
            stop = pos
            new = False
        else:
            stop = pos
    elif stop - start >= minlength and start != old_start:
        consensus_regions.append([start,stop])
        old_start = start
        new = True
    else:
        new = True
if stop - start >= minlength and start != old_start:
    consensus_regions.append([start,stop])

# write identified consensus regions to a csv table
outfile = f'{outdir}/consensus_regions.csv'
print(f'Writing consensus regions to: {outfile}\n')
with open(outfile, 'w') as outfile:
    outfile.write('Start\tStop\tLength\tConsensus_score\tSequence\n')
    outfile.write(f'1\t{len(consensus_seq_gapless)}\t{len(consensus_seq_gapless)}\t{"{:.3f}".format(avg(consensus_score))}\t{consensus_seq_gapless}\n')
    if not consensus_regions:
        print('#####################################################################\n### No Consensus regions detected for current parameter settings. ###\n#####################################################################')
        print('\nDone.')
        exit()
    for region in consensus_regions:
        scores = consensus_score[region[0]:region[1]+1]
        score_avg = avg(scores)
        scores = '\t'.join([str("{:.2f}".format(x)) for x in scores])
        print(f'Start: {region[0]+1} Stop: {region[1]+1} Length: {region[1]+1 - region[0]} Consensus_score: {"{:.3f}".format(score_avg)} Sequence: {consensus_seq_gapless[region[0]:region[1]+1]}')
        outfile.write(f'{region[0]+1}\t{region[1]+1}\t{region[1]+1 - region[0]}\t{"{:.3f}".format(score_avg)}\t{consensus_seq_gapless[region[0]:region[1]+1]}\n')
#    print(scores)

# get regions that have to be excluded from the primer design
# Necessary because primer3 can only automatically search for primers over multiple excluded areas in multiple regions.
consensus_flag = [0] * len(consensus_score)
for region in consensus_regions:
    for i in range(region[0],region[1]+1):
        consensus_flag[i] = 1
#consensus_flag_inverted = [1 if x == 0 else 0 for x in consensus_flag]
excluded_regions = []
start = 0
length = 0
new = True
for pos in range(len(consensus_flag)):
    if consensus_flag[pos] == 0:
        if new:
            start = pos
        length += 1
        new = False
    elif not new:
        excluded_regions.append([start,length])
        new = True
        length = 0
if excluded_regions:
    if not excluded_regions[-1] == [start,length] and length >=1:
        excluded_regions.append([start,length])

print('\nExcluded regions are (start,length):')
excluded_regions_out = ' '.join([','.join([str(y) for y in x]) for x in excluded_regions])
print(excluded_regions_out)




# read submitted primer3 parameter file
primer3_params = {}
if primer3file:
    print(f'Reading provided primer3 parameter file: {primer3file}')
    with open(primer3file) as infile:
        for line in infile:
            if not line.startswith('#'):
                line = line.strip().split('=')
                primer3_params[line[0]] = line[1]


# write adapted primer3 parameter files with consensus regions and excluded regions for consensus primer design
if primer3file:
    print(f'Writing primer3 run file to: {primer3_runfiles_dir}')
    primer3_runfile = f'{primer3_runfiles_dir}/consensus_threshold_{threshold}.txt'
    primer3_params['SEQUENCE_ID'] = f'consensus_sequence'
    primer3_params['SEQUENCE_TEMPLATE'] = consensus_seq_gapless
    primer3_params['SEQUENCE_EXCLUDED_REGION'] = excluded_regions_out
    primer3_params['SEQUENCE_INTERNAL_EXCLUDED_REGION'] = excluded_regions_out
    with open(primer3_runfile, 'w') as outfile:
        for param in primer3_params:
            outfile.write(f'{param}={primer3_params[param]}\n')



# use primer3 to calculate primers
print('\nUsing consensus regions to calculate primers with primer3.\n')
outfile = f'{primer3_result_dir}/consensus_threshold_{threshold}.txt'
with open(outfile, 'w') as outfile:
    subprocess.call(['primer3_core',primer3_runfile], stdout=outfile)


# read primers predicted by primer3 for html summary
print('Reading primer3 results.')
primer_dict = {}
primer3_dict = {}
with open(f'{outdir}/primer3_results/consensus_threshold_{threshold}.txt') as infile:
    for line in infile:
        # read primer sequences for primer alignment
        if re.match('^PRIMER_.*_SEQUENCE=',line):
            line1 = line.strip().split('=')
            if 'RIGHT' in line1[0]:
                primer_dict[line1[0]] = revcomp(line1[1])
            else:
                primer_dict[line1[0]] = line1[1]
        # read for command line output
        if re.match('^PRIMER_.*_EXPLAIN=',line):
            print('\t'.join(line.strip().split('=')))
        # read info for html output
        line = line.strip().split('=')
        param = line[0]
        value = line[1]
        if param:
            primer3_dict[param] = value 



# print message if no primers could be predicted for current parameters and input sequences
if not primer_dict:
    print('\n########################################################\n### No Primers found for current parameter settings. ###\n########################################################\n')
#    print('\nDone.')
#    exit()



# read user provided primers for visualization in primer alignment
if primers:
    known_primers, known_primer_headers = read_fasta(primers, delim, idpos)


# calculate negative consensus and add it to the final alignment (not considered during consensus primer design)
if negativesequencefile:
    negativesequences, negativesequence_header = read_fasta(negativesequencefile, delim, idpos)

    # convert fasta to panda dataframe to calculate consensus scores
    fasta_for_panda = {}
    for gene in negativesequences:
        fasta_for_panda[gene] = list(negativesequences[gene])
    fasta_pd_negative = pd.DataFrame.from_dict(fasta_for_panda, orient='index')

    # calculate nucleotide ratio for every position in the alighment
    # NOTE: non-standard nucleotides are counted as '-'
    consensus_dict_negative = {}
    for col in fasta_pd_negative:
        nucl_score = {'A':0,'T':0,'C':0,'G':0,'-':0}
        for nt in fasta_pd_negative[col]:
            if nt in nucl_score:
                nucl_score[nt] += 1
            else:
                nucl_score['-'] += 1
        for key in nucl_score:
            nucl_score[key] = nucl_score[key] / number_of_sequences
        consensus_dict_negative[col] = nucl_score

    #get consensus sequence
    consensus_seq_negative = ''
    consensus_score_negative = []
    for pos in consensus_dict_negative:
        consensus_nt = ''
        max_consensus_score = 0
        for nt in consensus_dict_negative[pos]:
            score = consensus_dict_negative[pos][nt]
            if score > max_consensus_score:
                consensus_nt = nt
                max_consensus_score = score
        consensus_seq_negative += consensus_nt
        consensus_score_negative.append(max_consensus_score)
        max_consensus_score = 0

    # remove gaps
    gap_positions = []
    for i in range(len(consensus_seq_negative)-1,-1,-1):
        if '-' == consensus_seq_negative[i]:
            gap_positions.append(i)
    consensus_seq_negative_gapless = list(consensus_seq_negative)
    for pos in gap_positions:
        del consensus_seq_negative_gapless[pos]
        del consensus_score_negative[pos]
    consensus_seq_negative_gapless = ''.join(consensus_seq_negative_gapless)


# the further alignment writing steps are necessary to guarantee proper alignments with short (primers) and longer sequences 

# write predicted and provided primers to fasta file for alignment with mafft
if primer_dict:
    outfile = f'{outdir}/primer_fasta.fna'
    with open(outfile, 'w') as outfile:
        for primer in primer_dict:
            outfile.write(f'>{primer}\n')
            outfile.write(f'{primer_dict[primer]}\n')
        if primers:
            for primer in known_primers:
                outfile.write(f'>{primer}\n')
                outfile.write(f'{known_primers[primer]}\n')

#write output fasta with consensus regions for the mafft alignment
outfile = f'{outdir}/consensusregions_fasta.fna'
with open(outfile, 'w') as outfile:
    for region in consensus_regions:
        scores = consensus_score[region[0]:region[1]+1]
        score_avg = avg(scores)
        fastahead = f'>Start:{region[0]+1}_Stop:{region[1]+1}_Length:{region[1]+1 - region[0]}_Consensusscore:{"{:.3f}".format(score_avg)}\n'
        outfile.write(fastahead)
        outfile.write(f'{consensus_seq_gapless[region[0]:region[1]+1]}\n')

# write consensus sequence and negative consensus sequence for alignment with mafft
outfile = f'{outdir}/consensussequences_fasta.fna'
with open(outfile, 'w') as outfile:
    outfile.write('>complete_gapless_consensus_sequence\n')
    outfile.write(f'{consensus_seq_gapless}\n')
    if negativesequencefile:
        outfile.write('>negative_consensus_sequence\n')
        outfile.write(f'{consensus_seq_negative_gapless}\n')

# write consensus regions for alignment with mafft
outfile = f'{outdir}/consensus_regions.fna'
with open(outfile, 'w') as outfile:
    for region in consensus_regions:
        scores = consensus_score[region[0]:region[1]+1]
        score_avg = avg(scores)
        fastahead = f'>Start:{region[0]+1}_Stop:{region[1]+1}_Length:{region[1]+1 - region[0]}_Consensusscore:{"{:.3f}".format(score_avg)}\n'
        outfile.write(fastahead)
        outfile.write(f'{consensus_seq_gapless[region[0]:region[1]+1]}\n')


# build final alignment in several steps with different parameters with all sequences
if primer_dict:
    print(f'\nBuilding and writing final_alignment to: {outdir}/final_alignment.fna\n')
# add consensus and negative consensus
with open(f'{outdir}/consensus_alignment.fna', 'w') as outfile:
    subprocess.call(['mafft','--auto', '--thread', threads,  '--add', f'{outdir}/consensussequences_fasta.fna', f'{outdir}/{unique}similarity_filtered_alignment.fna'], stdout=outfile, stderr=subprocess.DEVNULL)
# add consensusregions
with open(f'{outdir}/consensusregions_alignment.fna', 'w') as outfile:
    subprocess.call(['mafft','--6merpair', '--thread', threads,  '--addfragments', f'{outdir}/consensusregions_fasta.fna', f'{outdir}/consensus_alignment.fna'], stdout=outfile, stderr=subprocess.DEVNULL)
# add primers
if primer_dict:
    with open(f'{outdir}/final_alignment.fna', 'w') as outfile:
        subprocess.call(['mafft','--6merpair', '--thread', threads, '--adjustdirection',  '--addfragments', f'{outdir}/primer_fasta.fna', f'{outdir}/consensusregions_alignment.fna'], stdout=outfile, stderr=subprocess.DEVNULL)
    
    # read final_alignment fasta
    fasta, fastaheader = read_fasta(f'{outdir}/final_alignment.fna', delim, idpos)


## generate alignment html table

# Define color codes for each nucleotide
color_codes = {
    'A': 'lightgreen',
    'C': 'lightyellow',
    'G': 'lightsalmon',
    'T': 'lightblue',
    'U': 'lightcyan',
    '-': 'white'  # Gap symbol
}
if primer_dict:
    alignment_table = ''

    alignment_table += '<style>'
    alignment_table += 'table { border-collapse: collapse; }\n'
    alignment_table += 'table.outer { border: 1px solid black; }\n'
    alignment_table += 'table.outer td { font-family: monospace; padding: 1px; white-space:nowrap;}\n'
    alignment_table += 'table.outer .scrollable-col { max-height: 300px; overflow-y: auto; }\n'
    alignment_table += 'table.outer td.seq-name { font-weight: bold; font-size: 10px; padding: 0px; position: sticky; left: 0; z-index: 1; background-color: white; }\n'
    alignment_table += 'table.outer td.seq-data { position: relative; }\n'
    alignment_table += 'table.seq-table { border-collapse: collapse; }\n'
    alignment_table += 'table.seq-table td { border: none; font-size: 10px; padding: 0px; monospace;}\n'


    for nucleotide, color in color_codes.items():
        alignment_table += f'table.outer td.seq-data .{nucleotide} {{ background-color: {color}; }}\n'

    alignment_table += '</style>\n</head>\n<body>\n'
    alignment_table += '<div style="overflow-x:auto;">\n'  # Add a scrollable div container
    alignment_table += '<table class="outer">\n'

    for header, seq in fasta.items():
        seq_name = header
        alignment_table += '<tr>\n'
        alignment_table += f'<td class="seq-name">{seq_name}</td>\n'
        alignment_table += '<td class="seq-data scrollable-col">\n'
        alignment_table += '<table class="seq-table">\n'
        alignment_table += '<tr>\n'

        for nucleotide in seq:
            color = color_codes.get(nucleotide, 'white')
            alignment_table += f'<td style="background-color: {color};">{nucleotide}</td>'

        alignment_table += '\n</tr>\n'
        alignment_table += '</table>\n'
        alignment_table += '</div>\n'
        alignment_table += '</td>\n'
        alignment_table += '</tr>\n'


    alignment_table += '</table>\n'
    alignment_table += '</div>\n'  # Close the div container
    alignment_table += '</body>\n'

# make HTML output from primer3 result file
html_file = f'{outdir}/consensus_prime_summary.html'
with open(html_file, 'w') as outfile:
    primercount = 0
    summary_table = []
    primer_resulttables = {}
    primer_resulttable = []
    for param in primer3_dict:
        if 'STUCT' in param:
            primer3_dict[param] = '<p style="font-family:\'Lucida Console\', monospace"><pre>' + primer3_dict[param].replace("\\n","<br>").replace("U+2510", "&#x2510;").replace("U+2518","&#x2518;").replace("U+2502", "&#x2502;") + '</pre></p>'
        if str(primercount) in param or str(primercount+1) in param:
            if str(primercount+1) in param:
                primer_resulttables[primercount] = primer_resulttable
                primercount += 1
                primer_resulttable = []
            if re.match('PRIMER_LEFT_.+_SEQUENCE',param):
                primer_resulttable.append([param + "_FWD (5'-3')", primer3_dict[param]])
            elif re.match('PRIMER_RIGHT_.+_SEQUENCE',param):
                primer_resulttable.append([param + '_REV', revcomp(primer3_dict[param])])
                primer_resulttable.append([param + "_REVCOMP (5'-3')", primer3_dict[param]])
            else:
                primer_resulttable.append([param,primer3_dict[param]])
        else:
            if param == 'SEQUENCE_TEMPLATE':
                primer3_dict[param] = '<pre>' + insert_newlines(primer3_dict[param], linelen=100) + '</pre>'
            summary_table.append([param, primer3_dict[param]])
    primer_resulttables[primercount] = primer_resulttable
    # input parameters of ConsensusPrime run
    consensusparametertable = [
        ['Input fasta',alignment],
        ['Output directory',outdir],
        ['Keep Duplicate sequences', str(keepduplicates)],
        ['Consensus Similarity Threshold', str(consensussimilarity)],
        ['Consensus Threshold for Primerdesign', str(threshold)]
        ]
    if primer3file:
        consensusparametertable.append(['Primer3 parameter File',primer3file])
    if primers:
        consensusparametertable.append(['Provided Primers for alignment',primers])
    if negativesequencefile:
        consensusparametertable.append(['Sequences for negative consensus alignment', negativesequencefile])
    # write HTML document
    outfile.write('<html>\n<head>\n<style>\n')
    # define document font
    outfile.write('h1, h2, body { font-family: sans-serif;}\n')
    outfile.write('pre { font-size: 140%;}\n')
    # define table style
    outfile.write('table{ font-family: sans-serif; border-collapse: collapse;}\n')
    outfile.write('td { border: 1px solid #909090; text-align: left; padding: 8px; }\n')
    outfile.write('th { border: 1px solid #909090; text-align: center; padding: 8px; }\n')
    outfile.write('tr:nth-child(even) { background-color: #dddddd; }\n')
    outfile.write('</style>\n</head>\n<body>\n')
    outfile.write('<title>ConsensusPrime Report</title>\n')
    outfile.write('<h1>ConsensusPrime Report</h1>\n<br><hr><br>\n<h2>Content</h2>\n')
    # content overview
    outfile.write('<ul>\n')
    outfile.write('<li><a href="#userparams">User Defined Parameters</a></li>\n')
    outfile.write('<li><a href="#command">Run Command</a></li>\n')
    outfile.write('<li><a href="#alignment">Number of Sequences</a></li>\n')
    outfile.write('<li><a href="#finalalignment">Final Alignment</a></li>\n')
    outfile.write('<li><a href="#primer3param">Primer3 Parameter Summary</a></li>\n')
    outfile.write('<li><a href="#primers">Primers</a></li>\n')
    outfile.write('</ul>\n')
    outfile.write('<br><hr><br>\n')
    # content
    outfile.write('<h2 id="userparams">User Defined Parameters</h2>\n')
    outfile.write(to_html_table(consensusparametertable, header=['Parameter','Value']))
    outfile.write('<br><hr><br>\n')
    outfile.write('<h2 id="command">Run Command</h2>\n')
    outfile.write(' '.join(sys.argv) + '\n')
    outfile.write('<br><hr><br>\n')
    outfile.write('<h2 id="alignment">Number of Sequences</h2>\n')
#    outfile.write(HTML.table(alignment_statistic_table, header_row=['Fasta/Alignment', 'Number of Sequences']) + '\n')
    outfile.write(to_html_table(alignment_statistic_table, header=['Fasta/Alignment', 'Number of Sequences']))
    outfile.write('<br><hr><br>\n')
    # alignment table
    if primer_dict:
        outfile.write('<h2 id="finalalignment">Final Alignment</h2>\n')
        outfile.write(os.path.abspath(f'{outdir}/final_alignment.fna<p>\n'))
        outfile.write(alignment_table)
        outfile.write('<br><hr><br>\n')
    #outfile.write(HTML.link('Final alignment',f'{outdir}/final_alignment.fna'))
    outfile.write('<h2 id="primer3param">Primer3 Parameter Summary</h2>\n')
    outfile.write(to_html_table(summary_table, header=['Parameter','Value']))
    outfile.write('<br><hr><br>\n')
    outfile.write('<h2 id="primers">Primers</h2>\n')
    if primer_dict:
        for primerID in primer_resulttables:
            outfile.write(f'<h3>Primer {str(primerID)}</h3>\n')
            outfile.write(to_html_table(primer_resulttables[primerID], header=['Parameter','Value']))
    else:
        outfile.write('No Primers Detected.')
    outfile.write('\n</html>\n</body>\n')




print('(>°o°)> Done. <(°o°<)')

