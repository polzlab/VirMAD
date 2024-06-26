import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 

import warnings
warnings.filterwarnings("ignore")

import argparse

import tensorflow as tf

if tf.test.is_gpu_available():
    physical_devices = tf.config.list_physical_devices('GPU') 
    tf.config.experimental.set_memory_growth(physical_devices[0], True)
    tf.compat.v1.disable_eager_execution()
else:
    pass

from tensorflow.keras.models import load_model
from tensorflow.keras.models import Sequential

import numpy as np
import pandas as pd
import time
import plotly
import plotly.graph_objs as go
import random
import umap

################################################################################

def plot_html(df:pd.DataFrame, family_list, to_save = False, add_custom_legend = False, file_out = 'default.html'):
    list_of_plots = []
    for i, family in enumerate(family_list):
        adapted_embedded_f = df[df.family == family]
        color = random.randint(0, 0xFFFFFF)
        
        x1 = adapted_embedded_f['V1'].values
        y1 = adapted_embedded_f['V2'].values
        z1 = adapted_embedded_f['V3'].values
        
        if family == 'Negative':
            plot_f = go.Scatter3d(
                    x=x1,
                    y=y1,
                    z=z1,
                    hovertext = adapted_embedded_f.to_view.tolist(),
                    mode='markers',
                    marker=dict(
                        size=4,
                        color='red',            
                        symbol='circle',
                        line=dict(
                                color=color,
                                width=1
                            )
                    ),
                    name=family,
                )
            list_of_plots.append(plot_f)
        
        elif family == 'Positive':
            plot_f = go.Scatter3d(
                    x=x1,
                    y=y1,
                    z=z1,
                    hovertext = adapted_embedded_f.to_view.tolist(),
                    mode='markers',
                    marker=dict(
                        size=4,
                        color='green',            
                        symbol='circle',
                        line=dict(
                                color=color,
                                width=1
                            )
                    ),
                    name=family,
                )
            list_of_plots.append(plot_f)

        else:
            plot_f = go.Scatter3d(
                    x=x1,
                    y=y1,
                    z=z1,
                    hovertext = adapted_embedded_f.to_view.tolist(),
                    mode='markers',
                    marker=dict(
                        size=6,
                        color='yellow',            
                        symbol='circle',
                        line=dict(
                                color='yellow',
                                width=1
                            )
                    ),
                    name=family,
                )
            list_of_plots.append(plot_f)
            
        
    
    fig = go.Figure(data=list_of_plots)
    
    if add_custom_legend:
        fig.update_layout(
        legend=dict(
                x=0,
                y=1,
                traceorder="reversed",
                title_font_family="Calibri",
                font=dict(
                    family="Courier",
                    size=12,
                    color="black"
                ),
                bgcolor="White",
                bordercolor="Black",
                borderwidth=1
                )
            )
        
    plotly.offline.iplot(fig)
    if to_save:
        plotly.offline.plot(fig, filename=file_out)

################################################################################

def init_model():
    '''
    init GH CAS model instance
    '''
    model_gh = load_model('./models/Base_Model.h5')

    return model_gh


def get_model_GH(x, y, tsne=False):
    model = Sequential()                                
    model.add(tf.keras.layers.Conv2D(16, activation='elu', kernel_size=(5, 196), input_shape=(x, y, 1)))  
    model.add(tf.keras.layers.MaxPooling2D(pool_size=(2, 2177)))
    model.add(tf.keras.layers.Flatten())
    model.add(tf.keras.layers.Dense(16, activation='elu'))
    model.add(tf.keras.layers.Dropout(0.6))
    if not tsne:
        model.add(tf.keras.layers.Dense(1, activation='sigmoid'))
    
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
    return model

################################################################################

def make_prediction_on_GH(X, model_GH):
    '''
    X - encoded matrix for GH prediction
    model_GH - classifier instance
    '''
    model_tsne = get_model_GH(X.shape[1], X.shape[2], tsne=True)
    model_tsne.set_weights(model_GH.get_weights()[:-2])      
    prediction = model_GH.predict(X)[0][0]
    embedding = model_tsne.predict(X)
    
    return prediction, embedding


def reduce_dimmension(X_emb, labels):
    '''
    X_emb - GH model embedding with prediction sample
    labels - transformed to str labels
    '''
    metric_umap = 'correlation'
    adapted_embedded = umap.UMAP(n_neighbors=25,
                      n_components=3,
                      min_dist=0.5,
                      metric=metric_umap).fit_transform(X_emb)

    adapted_embedded = pd.DataFrame(adapted_embedded, columns=['V1', 'V2', 'V3'])
    adapted_embedded['to_view'] = ['pair_{}'.format(e) for e in range(len(X_emb))]
    adapted_embedded['family'] = labels
    return adapted_embedded


#################################################################################

#Arguments
parser = argparse.ArgumentParser()
parser.add_argument("-encoded_pair", help="host phage encoded array")
parser.add_argument("-background_X", help="VIS background X")
parser.add_argument("-background_y", help="VIS background y")
args = parser.parse_args()

#Init Model
print('Init Models...')
model_GH = init_model()

#Load Pair Matrix
print('Load Data...')
X = np.load(args.encoded_pair)
background_X = np.load(args.background_X)
background_y = np.load(args.background_y)

#Make Prediction
print('Make Prediction...')
p, e = make_prediction_on_GH(X,
                             model_GH)


output_dict = {'prediction' : p,
               'embedding' : e}


#Prepare Backround Embedding
_, e_background = make_prediction_on_GH(background_X,
                             model_GH)


print('Prepare Vis...')
X_emb = np.concatenate((e_background, e))
y = background_y.tolist() + [2]

labels = []

for prediction in y:
    if prediction == 1:
        labels.append('Positive')
    elif prediction == 0:
        labels.append('Negative')
    else:
        labels.append('Predicted')

print('Transform Results...')
adapted_embedded = reduce_dimmension(X_emb, labels)
label_list = adapted_embedded.family.value_counts().index.tolist()

print('Export Results...')
timepair = time.strftime("%Y%m%d_%H%M%S_pair.npy")
np.save(os.path.join('./output/prediction/', timepair), output_dict)

timehtml = time.strftime("%Y%m%d_%H%M%S_pair.html")
plot_html(adapted_embedded, label_list, to_save=True, file_out = os.path.join('./output/prediction/', timehtml))
