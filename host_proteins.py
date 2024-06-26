import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

import warnings
warnings.filterwarnings("ignore")

import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
import time

################################################################################

def prepare_pointer_results(df, system):
    """
    df - dataframe of one host padloc/defensefinder merged results
    """
    df = df[df.sys_id == system]
    median = df.start_global.median()
    df['diff_med'] = abs(df.start_global -  median)
    name = df[df.diff_med == df.diff_med.min()].protein.tolist()[0]
    sequence = df[df.diff_med == df.diff_med.min()].sequence.tolist()[0]
    cs = df[df.diff_med == df.diff_med.min()].checksum.tolist()[0]

    return name, sequence, cs

################################################################################

def read_proteom(file):
    '''
    method to read a list of names and sequences from a fasta file
    '''
    lon = []
    los = []
    with open(file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            lon.append(record.id)
            los.append(str(record.seq))
            
    return lon, los

def sequence_filling_method(range_t, names, sequences):
    '''
    a method supplementing the incomplete protein context around the central protein in a situation where the central protein is located at a distance of less than +/- 9 proteins from the proteome boundary, i.e. the protein list
    '''
    ctx_proteins = []
    ctx_names = []
         
    for idx in range(range_t[0], range_t[1]):
        if idx < 0:
            idx = 'to_fill'
        elif idx > range_t[1]:
            idx = 'to_fill'
        else:
            pass

        try:
            ctx_proteins.append(sequences[idx])
            ctx_names.append(names[idx])
        except:
            ctx_proteins.append('')
            ctx_names.append('None')

    if len(ctx_proteins) != 19:
        raise ValueError('Sequences list has not allowed dimension.')

    return ctx_proteins, ctx_names

################################################################################

def encode_bacteria_proteom(f, csv_path, system):
    '''
    f - host proteom .faa file path
    csv_path - padloc/defensefider unified format result csv file path
    
    '''
    #import proteom and related csv file
    names, sequences = read_proteom(f)
    df = pd.read_csv(csv_path)
    
    #select protein candidate for selected system
    name_sys, sequence_sys, _ = prepare_pointer_results(df, system)
    index_sys = names.index(name_sys)
    
    range_sys = (index_sys - 9, index_sys + 10)
    _, names_sys = sequence_filling_method(range_sys, names, sequences)
    
    #prepare output dataframe
    host_output_dataframe = pd.DataFrame()
    host_output_dataframe[system] = names_sys
    
    return host_output_dataframe, index_sys


################################################################################

#Arguments
parser = argparse.ArgumentParser()
parser.add_argument("-path_host_proteom", help="host proteom fasta file path")
parser.add_argument("-path_host_csv", help="host padloc/defensefinder csv file path")
args = parser.parse_args()


#Encode Phage
print('Prepare Collections...')
enc_h_dataframe, i_s,  = encode_bacteria_proteom(args.path_host_proteom,
                              args.path_host_csv, 
                              'cas')

print('System index: {}'.format(i_s))

#Export Results
print('Export Results...')
timehost = time.strftime("%Y%m%d_%H%M%S_host.csv")
enc_h_dataframe.to_csv(os.path.join('./output/host/', timehost), index=False)