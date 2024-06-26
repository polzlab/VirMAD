import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy
from Bio.Seq import Seq
from Bio import SeqIO
import random
import collections
from collections import OrderedDict
import pylab
import math

class CGREncoder:
    def __init__(self, kmer):
        self.kmer = kmer
        
    def read_sequence(self, filepath, file_type = 'fasta'):
        if file_type == 'fasta':
            seq = ''
            name = []
            for seq_record in SeqIO.parse(filepath, "fasta"):
                name.append(seq_record.id)
                seq += str(seq_record.seq)
        elif file_type == 'genbank':
            for seq_record in SeqIO.parse(filepath, "genbank"):
                name.append(seq_record.id)
                seq += str(seq_record.seq)
        else:
            print('Wrong type of file...')
            seq = None
            name = None
        
        seq = seq.replace('N', random.choice(['A', 'C', 'T', 'G']))
        seq = seq.replace('W', random.choice(['A', 'T']))
        seq = seq.replace('U', random.choice(['A', 'C', 'T', 'G']))
        seq = seq.replace('R', random.choice(['A', 'G']))
        seq = seq.replace('Y', random.choice(['C', 'T']))
        seq = seq.replace('K', random.choice(['T', 'G']))
        seq = seq.replace('M', random.choice(['A', 'C']))
        seq = seq.replace('S', random.choice(['C', 'G']))
        seq = seq.replace('B', random.choice(['C', 'T', 'G']))
        seq = seq.replace('D', random.choice(['A', 'T', 'G']))
        seq = seq.replace('H', random.choice(['A', 'C', 'T']))
        seq = seq.replace('V', random.choice(['A', 'C', 'G']))
        return seq, name
    
    def count_kmers(self, sequence):
        d = collections.defaultdict(int)
        for i in range(len(sequence)-(self.kmer-1)):
            d[sequence[i:i+self.kmer]] +=1
        return d
    
    def probabilities(self, kmer_count, sequence):
        probabilities = collections.defaultdict(float)
        N = len(sequence)
        for key, value in kmer_count.items():
            probabilities[key] = float(value) / (N - self.kmer + 1)
        return probabilities
        
    def chaos_game_representation(self, probabilities):
        array_size = int(math.sqrt(4**self.kmer))
        chaos = []
        for i in range(array_size):
            chaos.append([0]*array_size)
 
        maxx = array_size
        maxy = array_size
        posx = 1
        posy = 1
        
        for key, value in probabilities.items():
            for char in key:
                if char == "T":
                    posx += maxx / 2
                elif char == "C":
                    posy += maxy / 2
                elif char == "G":
                    posx += maxx / 2
                    posy += maxy / 2
                maxx = maxx / 2
                maxy /= 2
                
            posy = int(posy)
            posx = int(posx)

            chaos[posy-1][posx-1] = value
            maxx = array_size
            maxy = array_size
            posx = 1
            posy = 1
 
        return chaos