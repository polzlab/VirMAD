import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 

import warnings
warnings.filterwarnings("ignore")

import argparse
import time
from FCGR.enc import CGREncoder
import numpy as np
import matplotlib.image
import joblib
import tensorflow as tf
import cv2

if tf.test.is_gpu_available():
    physical_devices = tf.config.list_physical_devices('GPU') 
    tf.config.experimental.set_memory_growth(physical_devices[0], True)
    tf.compat.v1.disable_eager_execution()
else:
    pass

from tensorflow.keras.layers import Dense, Input, Conv2D, GlobalMaxPooling2D, Dropout
from tensorflow.keras.models import Model
from tensorflow.keras.models import load_model

#################################################################################

def init_pointer_models():
    '''
    taxonomy models init
    '''
    model_taxo_bact = load_model('./models/order_taxo_model_bac.h5', compile=False)
    model_taxo_vir = load_model('./models/order_taxo_model_vir.h5', compile=False)
    return model_taxo_bact, model_taxo_vir

def init_taxo_emb_models(model_bact, model_vir):
    '''
    emb models init
    '''
    
    model_tsne_b = get_model_taxo_b(tsne=True)
    model_tsne_b.set_weights(model_bact.get_weights()[:-2])

    model_tsne_v = get_model_taxo_v(tsne=True)
    model_tsne_v.set_weights(model_vir.get_weights()[:-2])
    
    return model_tsne_b, model_tsne_v

def init_scalers():
    '''
    taxonomy encoders init
    '''
    scaler_host = joblib.load('./scalers/taxonomy_host_image_scaler.save')
    scaler_phage = joblib.load('./scalers/taxonomy_phage_image_scaler.save')

    return scaler_host, scaler_phage

#################################################################################


def encode_image_by_FCGR(filename, encoder, output_type = 'list'):
    '''
    FCGR encoding method
    '''
    seq, name = encoder.read_sequence(filename, file_type = 'fasta')
    org_kmers = encoder.count_kmers(seq)
    org_prob = encoder.probabilities(org_kmers, seq)
    encoded_img = encoder.chaos_game_representation(org_prob)
    
    if output_type == 'array':
        encoded_img = np.array(encoded_img)
        
    return encoded_img, name

#################################################################################

def save_image(arr, path):
    '''
    save np array as image
    '''
    rescaled = cv2.resize(arr, dsize=(640, 640), interpolation=cv2.INTER_LINEAR)
    matplotlib.image.imsave(path, rescaled)
    
#################################################################################

def scale_data_taxo_img(X, scaler):
    '''
    FCGR matrix scaling method
    '''
    X = X.reshape(1, 128*128)
    X_scaled = scaler.transform(X)
    X_scaled = X_scaled.reshape(1, 128, 128)
    return X_scaled

#################################################################################

def get_model_taxo_b(tsne=False):
    '''
    host taxonomy model structure
    '''
    # create model
    x_input = Input(shape=(128, 128, 1))
    x = Conv2D(128, kernel_size=(3,3), activation='relu')(x_input)
    x = GlobalMaxPooling2D()(x)
    x = Dense(128, activation='relu')(x)
    x = Dropout(0.5)(x)
    if not tsne:
        x = Dense(79, activation='softmax')(x)    
    # Compile model
    model = Model(inputs=x_input, outputs=x, 
                  name='CNN_Model')
    return model

def get_model_taxo_v(tsne=False):
    '''
    phage taxonomy model structure
    '''
    # create model
    x_input = Input(shape=(128, 128, 1))
    x = Conv2D(48, kernel_size=(5,5), activation='relu')(x_input)
    x = Dropout(0.5)(x)
    x = Conv2D(64, kernel_size=(4,4), activation='relu')(x)
    x = GlobalMaxPooling2D()(x)
    x = Dense(128, activation='relu')(x)
    x = Dropout(0.5)(x)
    if not tsne:
        x = Dense(5, activation='softmax')(x)    
    # Compile model
    model = Model(inputs=x_input, outputs=x, 
                  name='CNN_Model')
    return model

#################################################################################

def encode_taxonomy(f_bact, f_vir, scaler_host, scaler_phage, encoder_taxo, model_taxo_b, model_taxo_v):
    '''
    taxonomy vectors encoding method
    
    f_bact - host fasta file path
    f_vir - phage fasta file path
    scaler_host - FCGR host scaler
    scaler_phage - FCGR phage scaler 
    encoder_taxo - (CGREncoder(kmer=7)) instance
    model_taxo_b - emb host taxonomy model
    model_taxo_v - emb phage taxonomy model

    '''

    #FCGR genome encoding
    encoded_img_host, host_name = encode_image_by_FCGR(f_bact, encoder_taxo, output_type='array')
    encoded_img_phage, phage_name = encode_image_by_FCGR(f_vir, encoder_taxo, output_type='array')
    
    #FCGR matrix scaling
    scaled_host = scale_data_taxo_img(encoded_img_host, scaler_host)
    scaled_phage = scale_data_taxo_img(encoded_img_phage, scaler_phage)
       
    #embedding generation
    X_emb_b = model_taxo_b.predict(scaled_host.reshape(1, 128, 128, 1))
    X_emb_v = model_taxo_v.predict(scaled_phage.reshape(1, 128, 128, 1))
    

    #result dictionary
    taxonomy_vector_dict = {'host' : X_emb_b,
                            'phage' : X_emb_v}

    return taxonomy_vector_dict, scaled_host.reshape(128, 128), scaled_phage.reshape(128, 128)



#################################################################################

#Arguments
parser = argparse.ArgumentParser()
parser.add_argument("-path_host", help="host fasta file path")
parser.add_argument("-path_phage", help="phage fasta file path")
args = parser.parse_args()


#Init Models
print('Init Models...')
model_taxo_bact, model_taxo_vir = init_pointer_models()

#Init Emb Models
model_taxo_emb_b, model_taxo_emb_v = init_taxo_emb_models(model_taxo_bact, model_taxo_vir)

#Init Scalers
print('Init Scalers...')
scaler_host, scaler_phage = init_scalers()

#Init Encoder
print('Init Encoder...')
encoder_taxo = CGREncoder(kmer=7)

#Encode Taxonomy
print('Prepare Embeddings...')
enc_t_dict, host_img, phage_img = encode_taxonomy(args.path_host, 
                         args.path_phage, 
                         scaler_host, scaler_phage,
                         encoder_taxo,
                         model_taxo_emb_b, model_taxo_emb_v)


#Export Results
print('Export Results...')
timepair = time.strftime("%Y%m%d_%H%M%S_pair.npy")
timehost = time.strftime("%Y%m%d_%H%M%S_host.png")
timephage = time.strftime("%Y%m%d_%H%M%S_phage.png")


np.save(os.path.join('./output/taxonomy/', timepair), enc_t_dict)
save_image(host_img, os.path.join('./output/taxonomy/', timehost))
save_image(phage_img, os.path.join('./output/taxonomy/', timephage))