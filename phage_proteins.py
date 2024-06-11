import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 

import warnings
warnings.filterwarnings("ignore")

import tensorflow as tf

if tf.test.is_gpu_available():
    physical_devices = tf.config.list_physical_devices('GPU') 
    tf.config.experimental.set_memory_growth(physical_devices[0], True)
    tf.compat.v1.disable_eager_execution()
else:
    pass
    
import argparse
from tensorflow.keras.models import load_model
import pickle
from Bio import SeqIO
import numpy as np
import pandas as pd
import time


################################################################################

def init_pointer_models():
    '''
    models import
    '''
    model_counterdef_xgb = pickle.load(open('./models/counterdef_XGB.pkl', "rb"))
    model_counterdef_cnn = load_model('./models/counterdef_CNN.h5', compile=False)
    model_orf_stage1 = load_model('./models/orf_stage1_pointer_model.h5', compile=False)
    model_orf_stage2 = load_model('./models/orf_stage2_pointer_model.h5', compile=False)
    return model_counterdef_xgb, model_counterdef_cnn, model_orf_stage1, model_orf_stage2

################################################################################

def read_proteom(file):
    '''
    fasta file reading function
    '''
    lon = []
    los = []
    with open(file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            lon.append(record.id)
            los.append(str(record.seq))
            
    return lon, los


#################################################################################

def encode_lengths_orf(names, sequences, encoded):
    '''
    orf protein lengths encoding
    '''
    lengths = []
    for s in sequences:
        lengths.append(len(s))

    lengths = np.array(lengths).reshape(len(lengths), 1) 
    enc = np.concatenate([encoded, lengths], axis = 1)
    return names, enc

def encode_lengths_counterdef(names, sequences, encoded):
    '''
    counterdef protein lenths encoding
    '''
    lengths = []
    for s in sequences:
        lengths.append(len(s))
        
    lengths = np.array(lengths).reshape(len(lengths), 1) / 100
    enc = np.concatenate([encoded, lengths], axis = 1)
    return names, enc

def sequence_filling_method(range_t, names, sequences):
    '''
    context filling method
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

#################################################################################

def do_prediction_orf(virus, names, sequences, encoded, model, smodel):
    '''
    orf candidate selecting function
    
    virus - virus name
    names - protein names list in proteom file
    sequences - protein sequence list in protein file
    encoded - encoded matrix with lengths of proteins
    model - orf protein pointer model
    smodel - structural protein filtering model
    
    '''
    
    #matrix transform
    encoded = np.array(encoded)
    encoded = encoded.reshape(encoded.shape[0], encoded.shape[1], 1)
    
    #orf prediction
    poss = model.predict(encoded)
    poss = [e[0] for e in poss]
    predictions = [round(value) for value in poss]
    positives = [idx for idx, element in enumerate(predictions) if element == 1]


    if 1 not in predictions:
        best_score = max(poss)
        index_best = poss.index(best_score)
        name = names[index_best]
        seq = sequences[index_best]
    
    else:
        #structural proteins filtering
        encodeds = encoded[positives]
        poss_f = [poss[i] for i in positives]
        names_f = [names[i] for i in positives]
        seq_f = [sequences[i] for i in positives]
        poss_st = smodel.predict(encodeds)
        poss_st = [e[0] for e in poss_st]
        
        positives_st = [idx for idx, element in enumerate(poss_st) if element < 0.5]
        poss_f = [poss_f[i] for i in positives_st]
        names_f = [names_f[i] for i in positives_st]
        seq_f = [seq_f[i] for i in positives_st]

        try:
            best_score = max(poss_f)
            index_best = poss_f.index(best_score)
            name = names_f[index_best]
            seq = seq_f[index_best]
        except:
            best_score = max(poss)
            index_best = poss.index(best_score)
            name = names[index_best]
            seq = sequences[index_best]
    
    return virus, name, seq, best_score



def do_prediction_counterdef(virus, names, sequences, encoded, model, cmodel):
    '''
    counterdef candidates selecting function
    model - XGB classifier
    cmodel - CNN classifier
    '''
    X = encoded
    X_cnn = X.reshape(X.shape[0], X.shape[1], 1)
    counterdef_r1 = model.predict_proba(X)[:,1]
    counterdef_r2 = cmodel.predict(X_cnn).reshape(X_cnn.shape[0],)
    counterdef_r = (counterdef_r1 * 0.4 + counterdef_r2 * 0.6)

    idx = np.argsort(np.array(counterdef_r))

    counterdef_r = np.array(counterdef_r)[idx]
    names = np.array(names)[idx]
    sequences = np.array(sequences)[idx]

    counterdef1 = names[-1]
    counterdef2 = names[-2]
    
    result1 = counterdef_r[-1]
    result2 = counterdef_r[-2]
    
    seq1 = sequences[-1]
    seq2 = sequences[-2]
    
    return virus, counterdef1, counterdef2, result1, result2, seq1, seq2

#################################################################################

def encode_virus_proteom(f, encoded, model_orf_stage1, model_orf_stage2, model_counterdef_xgb, model_counterdef_cnn):
    '''
    select phage protein regions for anti-defense candidates
    
    f - .faa file proteom path
    encoded - .faa encoded proteins in np.array format path
    models_orf_stage1 - stage 1 orf classifier
    models_orf_stage2 - stage 2 orf classifier
    model_counterdef_xgb - anti-defense XGB classifier
    model_counterdef_cnn - anti-defense CNN classifier
    '''
    
    virus = f.split(',')[-1][:-4]
    names, sequences = read_proteom(f) #read proteom file
    encoded_m = np.load(encoded)
    
    print(len(names), len(sequences), encoded_m.shape)
    
    names, enc_orf = encode_lengths_orf(names, sequences, encoded_m) #orf lengths encoding
    names, enc_counterdef = encode_lengths_counterdef(names, sequences, encoded_m) #counterdef lengths encoding
  
    
    #select orf candidate
    virus, orf_name, seq, best_score = do_prediction_orf(virus, names, sequences, enc_orf, model_orf_stage1, model_orf_stage2)
    #select counterdef candidates
    virus, counterdef1, counterdef2, result1, result2, seq1, seq2 = do_prediction_counterdef(virus, names, sequences, enc_counterdef, model_counterdef_xgb, model_counterdef_cnn) 

    #check proteom list indices for each 
    index_orf = names.index(orf_name) 
    index_counterdef1 = names.index(counterdef1)
    index_counterdef2 = names.index(counterdef2)
    
    
    #select range of interesting orf and counterdef proteins
    range_orf = (index_orf - 9, index_orf + 10)
    range_counterdef1 = (index_counterdef1 - 9, index_counterdef1 + 10)
    range_counterdef2 = (index_counterdef2 - 9, index_counterdef2 + 10)
    
    #select names of orf and counterdef proteins to encode
    _, names_orf = sequence_filling_method(range_orf, names, sequences)
    _, names_counterdef1 = sequence_filling_method(range_counterdef1, names, sequences)
    _, names_counterdef2 = sequence_filling_method(range_counterdef2, names, sequences)        
    
    
    #prepare dataframe with protein candidates to encode
    virus_output_dataframe = pd.DataFrame()
    virus_output_dataframe['orf'] = names_orf
    virus_output_dataframe['counterdef1'] = names_counterdef1
    virus_output_dataframe['counterdef2'] = names_counterdef2


    return virus_output_dataframe, index_orf, index_counterdef1, index_counterdef2



#################################################################################

#Arguments
parser = argparse.ArgumentParser()
parser.add_argument("-path_phage_proteom", help="phage proteom fasta file path")
parser.add_argument("-path_phage_encoded", help="phage encoded numpy file path")
args = parser.parse_args()

#Init Models
print('Init Models...')
model_counterdef_xgb, model_counterdef_cnn, model_orf_stage1, model_orf_stage2 = init_pointer_models()

#Encode Phage
print('Prepare Collections...')
enc_v_dataframe, i_s, i_a1, i_a2 = encode_virus_proteom(args.path_phage_proteom, 
                              args.path_phage_encoded, 
                              model_orf_stage1, 
                              model_orf_stage2, 
                              model_counterdef_xgb, 
                              model_counterdef_cnn)

print('orf index: {}, counterdef1 index: {}, counterdef2 index: {}'.format(i_s, i_a1, i_a2))

#Export Results
print('Export Results...')
timephage = time.strftime("%Y%m%d_%H%M%S_phage.csv")
enc_v_dataframe.to_csv(os.path.join('./output/phage/', timephage), index=False)