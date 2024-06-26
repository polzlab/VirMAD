import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

import warnings
warnings.filterwarnings("ignore")

import argparse
import pandas as pd
import numpy as np
import shap
import tensorflow as tf
import plotly.express as px
import time
from shap.utils._legacy import kmeans

tf.compat.v1.disable_v2_behavior()
if tf.test.is_gpu_available():
    physical_devices = tf.config.list_physical_devices('GPU') 
    tf.config.experimental.set_memory_growth(physical_devices[0], True)
    tf.compat.v1.disable_eager_execution()
else:
    pass


#################################################################################

def init_model():
    '''
    init classification model
    '''
    model = tf.keras.models.load_model('./models/Base_Model.h5')
    return model

#################################################################################
def normalize_data(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def prepare_activation_vis(arr, savepath):
    '''
    html heatmap export function
    '''
    fig = px.imshow(arr, aspect='auto')
    fig.update_layout(width=1000, height=400, margin=dict(l=10, r=10, b=10, t=10))
    fig.update_xaxes(showgrid=True, color = 'black', gridcolor='red', gridwidth=2, linewidth=5, layer='above traces')
    fig.update_yaxes(showgrid=True, color = 'black', gridcolor='red', gridwidth=2, linewidth=5, layer='above traces')
    fig.write_html(savepath)


def shap_image(shap_values, pixel_values=None, labels=None, width=20, aspect=0.2, hspace=0.2, labelpad=None):
    # support passing an explanation object
    if str(type(shap_values)).endswith("Explanation'>"):
        shap_exp = shap_values
        feature_names = [shap_exp.feature_names]
        ind = 0
        if len(shap_exp.base_values.shape) == 2:
            shap_values = [shap_exp.values[..., i] for i in range(shap_exp.values.shape[-1])]
        else:
            raise Exception("Number of outputs needs to have support added!! (probably a simple fix)")
        if pixel_values is None:
            pixel_values = shap_exp.data
        if labels is None:
            labels = shap_exp.output_names

    multi_output = True
    if type(shap_values) != list:
        multi_output = False
        shap_values = [shap_values]

    # make sure labels
    if labels is not None:
        labels = np.array(labels)
        assert labels.shape[0] == shap_values[0].shape[0], "Labels must have same row count as shap_values arrays!"
        if multi_output:
            assert labels.shape[1] == len(shap_values), "Labels must have a column for each output in shap_values!"
        else:
            assert len(labels.shape) == 1, "Labels must be a vector for single output shap_values."

    label_kwargs = {} if labelpad is None else {'pad': labelpad}

    # plot our explanations
    x = pixel_values

    fig_size = np.array([20, shap_values[0].shape[0] * 3])

    # additional feature to return list of transformed images:
    transformed_images = []

    for row in range(x.shape[0]):
        x_curr = x[row].copy()

        # make sure we have a 2D array for grayscale
        if len(x_curr.shape) == 3 and x_curr.shape[2] == 1:
            x_curr = x_curr.reshape(x_curr.shape[:2])
        if x_curr.max() > 1:
            x_curr /= 255.

        # get a grayscale version of the image
        if len(x_curr.shape) == 3 and x_curr.shape[2] == 3:
            x_curr_gray = (0.2989 * x_curr[:, :, 0] + 0.5870 * x_curr[:, :, 1] + 0.1140 * x_curr[:, :, 2])
            x_curr_disp = x_curr
        elif len(x_curr.shape) == 3:
            x_curr_gray = x_curr.mean(2)
            flat_vals = x_curr.reshape([x_curr.shape[0] * x_curr.shape[1], x_curr.shape[2]]).T
            flat_vals = (flat_vals.T - flat_vals.mean(1)).T
            means = kmeans(flat_vals, 3, round_values=False).data.T.reshape([x_curr.shape[0], x_curr.shape[1], 3])
            x_curr_disp = (means - np.percentile(means, 0.5, (0, 1))) / (
                        np.percentile(means, 99.5, (0, 1)) - np.percentile(means, 1, (0, 1)))
            x_curr_disp[x_curr_disp > 1] = 1
            x_curr_disp[x_curr_disp < 0] = 0
        else:
            x_curr_gray = x_curr
            x_curr_disp = x_curr
        if len(shap_values[0][row].shape) == 2:
            abs_vals = np.stack([np.abs(shap_values[i]) for i in range(len(shap_values))], 0).flatten()
        else:
            abs_vals = np.stack([np.abs(shap_values[i].sum(-1)) for i in range(len(shap_values))], 0).flatten()
        max_val = np.nanpercentile(abs_vals, 99.9)
        for i in range(len(shap_values)):
            sv = shap_values[i][row] if len(shap_values[i][row].shape) == 2 else shap_values[i][row].sum(-1)
            transformed_images.append(sv)
    return normalize_data(np.array(transformed_images))


def find_extremum_values(X, mode='both'):
    """
    Prepare most significant protein positions from X matrix

    """

    highest_values_indices = []
    highest_values = []

    for i in range(X.shape[0]):
        if mode == 'both':
            max_indices = np.argmax(X[i], axis=1)
            max_values = [X[i][j][max_indices[j]] for j in range(6)]
            total_max = np.argmax(np.array(max_values))
            highest_values_indices.append(f'{total_max}-{max_indices[total_max]}')
            highest_values.append(f'{max_values[total_max]}')

        elif mode == 'bact':
            max_values = []
            max_sort = np.argsort(X[i], axis=1)
            max_indices = [max_sort[0][-1], max_sort[0][-2], max_sort[0][-3]]
            max_indices_transformed = []
            max_values = []
            for j in range(3):
                max_indices_transformed.append(f'{i}-{max_indices[j]}')
                max_values.append(X[i][i][max_indices[j]])
            highest_values_indices = max_indices_transformed
            highest_values = max_values

        elif mode == 'vir':
            max_values = []
            max_indices = np.argmax(X[i], axis=1)
            max_indices_transformed = []
            for j in range(3):
                max_indices_transformed.append(f'{j + 3}-{max_indices[j]}')
                max_values.append(X[i][j][max_indices[j]])
            total_max = np.argmax(np.array(max_values))
            highest_values_indices = max_indices_transformed
            highest_values = max_values

    return highest_values_indices, highest_values


def analyze_results(X, sys_id, mode='max'):
    """
    Transform shapley values 128 -> 1
    """
    X = np.swapaxes(X, 1, 0)

    vector_sums = []
    for i in range(6):
        Xf = X[i]
        protein_sums = []
        for j in range(19):
            if mode == 'median':
                protein_sums.append(np.median(Xf[:, j * 128:(j + 1) * 128], axis=1))
            elif mode == 'mean':
                protein_sums.append(np.mean(Xf[:, j * 128:(j + 1) * 128], axis=1))
            elif mode == 'sum':
                protein_sums.append(np.sum(Xf[:, j * 128:(j + 1) * 128], axis=1))
            elif mode == 'std':
                protein_sums.append(np.std(Xf[:, j * 128:(j + 1) * 128], axis=1))
            elif mode == 'max':
                protein_sums.append(Xf[:, j * 128:(j + 1) * 128].max(axis=1))
            elif mode == 'min':
                protein_sums.append(Xf[:, j * 128:(j + 1) * 128].min(axis=1))
        vector_sums.append(protein_sums)

    X = np.swapaxes(np.array(vector_sums), 2, 1)
    result = np.swapaxes(X, 1, 0)

    dict_system = {0: sys_id,
                   1: 'taxo bact',
                   2: 'taxo phage',
                   3: 'cds2',
                   4: 'cds1_1',
                   5: 'cds1_2'}

    for i in range(6):
        max_v = X[i].max()
        min_v = X[i].min()
        print(f'System: {dict_system[i]}, max: {max_v}, min: {min_v}')

    return np.array(result)


#################################################################################


def prepare_shap_analysis(X_path, background_path, sys_id):
    '''
    X_path - path to encoded X pair matrix
    background_path - path to model related background matrix
    sys_id - system related to model
    '''
    #init model, X and background
    model = init_model()
    X = np.load(X_path)
    background = np.load(background_path)
    
    # init shap mechanism and calculate results
    e = shap.DeepExplainer(model, background)
    shap_values = e.shap_values(X)
    img = shap_image(shap_values[0], -X)
    vector_max = analyze_results(img, sys_id, mode='max')

    
    # split matrix into host and vir section
    vector_host_max = vector_max[:, 0:1]
    vector_vir_max = vector_max[:, 3:]
    
    dict_matrix_max = {'bact': vector_host_max,
                        'vir': vector_vir_max}
    
    
    # ###################################### BACTERIA SHAP MODE ##########################################
    mode = 'bact'
    highest_host_indices, highest_host_values = find_extremum_values(dict_matrix_max[mode], mode=mode)

    # ##################################### VIRUS SHAP MODE ##############################################
    mode = 'vir'
    highest_virus_indices, highest_virus_values = find_extremum_values(dict_matrix_max[mode], mode=mode)

    
    return highest_host_indices, highest_virus_indices, highest_host_values, highest_virus_values, vector_max


def decode_results(highest_host_indices, highest_virus_indices, highest_host_values, highest_virus_values, sys_id, host_csv_path, phage_csv_path, output_path):
    '''
    decode activation positions to protein names
    '''
    
    dict_system = {0: sys_id,
               1: 'taxo bact',
               2: 'taxo phage',
               3: 'cds2',
               4: 'cds1_1',
               5: 'cds1_2'}
    
    dfh = pd.read_csv(host_csv_path)
    dfp = pd.read_csv(phage_csv_path)
    
    host_scores = []
    phage_scores = []
    
    host_candidates = []
    for i, r in enumerate(highest_host_indices):
        host_candidates.append(dfh[dict_system[int(r.split('-')[0])]].tolist()[int(r.split('-')[1])])
        host_scores.append(highest_host_values[i])
        
    phage_candidates = []
    for i, r in enumerate(highest_virus_indices):
        phage_candidates.append(dfp[dict_system[int(r.split('-')[0])]].tolist()[int(r.split('-')[1])])
        phage_scores.append(highest_virus_values[i])
    
    dfo_b = pd.DataFrame()
    dfo_v = pd.DataFrame()
    
    max_host_score = max(host_scores)
    max_phage_score = max(phage_scores)
    
    dfo_b['host_protein_candidate'] = host_candidates
    dfo_b['SCORE_bact'] = [(e / max_host_score) * 100 for e in host_scores]
    
    dfo_v['phage_protein_candidate'] = phage_candidates
    dfo_v['SCORE_virus'] = [(e / max_phage_score) * 100 for e in phage_scores]
    
    dfo = pd.merge(left=dfo_v, right=dfo_b, how='cross')
    dfo['SCORE'] = (dfo.SCORE_virus + dfo.SCORE_bact) / 2
    dfo = dfo.sort_values(by='SCORE', ascending=False).reset_index(drop=True)
    dfo = dfo[['SCORE', 'host_protein_candidate', 'phage_protein_candidate']]
    dfo.to_csv(output_path, index=False)
    


#################################################################################

#Arguments
parser = argparse.ArgumentParser()
parser.add_argument("-path_X", help="encoded Guest Host matrix")
parser.add_argument("-path_background", help="SHAP background related to model")
parser.add_argument("-path_phage_proteins", help="phage proteins result csv")
parser.add_argument("-path_host_proteins", help="host proteins result csv")
args = parser.parse_args()

#Init Models
print('Init Models...')
hhi, hvi, hhv, hvv, vm = prepare_shap_analysis(args.path_X, args.path_background, 'cas')

print('Export Results...')
timepairvis = time.strftime("%Y%m%d_%H%M%S_activation_map.html")
timepaircsv = time.strftime("%Y%m%d_%H%M%S_protein_candidates.csv")
prepare_activation_vis(vm[0], os.path.join('./output/postanalysis/', timepairvis))
decode_results(hhi, hvi, hhv, hvv, 'cas', args.path_host_proteins, args.path_phage_proteins, os.path.join('./output/postanalysis/', timepaircsv))