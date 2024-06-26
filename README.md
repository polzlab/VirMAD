# VirMAD
## Interpretable Machine Learning for Bacterial Defense and Viral Antidefense Pair Inference
VirMAD works by:
1. Encoding viral and bacterial proteins using a BERT-architecture protein language model
2. Encoding bacterial and viral taxonomy using CNNs on Chaos-Game-Representations of nucleotide sequences
3. Using a Convolutional Neural Network (CNN) to predict infection outcome based on coupled bacterial-virus input tensors and experimental cross-infection matrix derived labels
4. Identify putative defense and antidefense pairs as co-activated viral proteins and bacterial defenses using a Shapley value analysis.

## Table of Contents

### Genomic Analysis
1. [Comparative Genomics](#comparative-genomics)

### VirMAD
2. [Introduction](#introduction)
3. [Environment Setup](#environment-setup)
4. [Running Scripts](#running-scripts)
5. [Jupyter Notebook Demo](#jupyter-notebook-demo)

## Comparative Genomics
The subfolder Comparative_Genomics contains julia and bash scripts to carry out a comparative analysis on a group of host and viral genomes, while utilizing experimental cross-infection matrix data. The scripts are tailored for the specific case of the Nahant cross infection matrix, but can be trivially adjusted to any other genome resolved cross-infectio matrix.

## VirMAD
## Introduction

This repository contains scripts and data for analyzing VirMAD project related data. The scripts are designed to process genomic and proteomic data related to VirMAD, perform various analyses, and generate visualizations.

## Environment Setup

To set up the environment for running the scripts, follow these steps:

1. Clone this repository to your local machine
 

2. Download the required data from [Zenodo](https://zenodo.org/records/12548715).
   unzip the file and place the "models", "data", and "output" folders 
   under the  VirMAD folder you just cloned

3. Navigate to the repository directory:
   ```
   cd VirMAD
   ```

4. Create a conda environment named `VirMAD-env` with Python 3.10:
   ```
   conda create --name VirMAD-env python=3.10
   ```

5. Activate the `VirMAD-env` environment:
   ```
   conda activate VirMAD-env
   ```

6. Install the required packages:
   ```
   pip install -r requirements.txt
   ```

## Running Scripts

To run the scripts for this project, follow these steps:

1. Ensure that the `VirMAD-env` conda environment is activated. If it's not already activated, you can activate it using:
    ```
    conda activate VirMAD-env
    ```

2. Navigate to the directory of the GitHub VirMAD Repository:
    ```
    cd VirMAD
    ```

3. Once in the repository directory, you can run the desired scripts using Python:

   - **`taxonomy.py`**:
     ```bash
     python taxonomy.py -path_host ./data/genome/host/GCA_900186385.1.fna -path_phage ./data/genome/phage/NC_031911.fna
     ```
     Here, `-path_host` and `-path_phage` are command-line arguments specifying the input fasta formatted files for the host genome and the phage genome respectively. 

     After running the script, you can expect the following output:
     - Two PNG files representing encoded by FCGR algorithm phage and host genomes.
     - Numpy file - dictionary with taxonomy embeddings for phage and host.

     The output directory for this script is `./output/taxonomy`.

   - **`phage_proteins.py`**:
     ```bash
     python phage_proteins.py -path_phage_proteom ./data/proteom/phage/NC_031911.faa -path_phage_encoded ./data/encoded/NC_031911.npy
     ```
     Here: 
     - `-path_phage_proteom` is the proteom fasta file containing a set of proteins in the organism, which is a MultiPhATE2 output of the genome fna file. 
     - `-path_phage_encoded` is a NumPy file containing the encoding of all proteins from the proteom fasta file. 

     After running the script, you can expect the following output:
     - A CSV formatted file with list of phage 3 counter defence proteins included in analysis (cds1_1, cds1_2, cds2).

     The output directory for this script is `./output/phage`.

   - **`host_proteins.py`**:
     ```bash
     python host_proteins.py -path_host_proteom ./data/proteom/host/GCA_900186385.1.faa -path_host_csv ./data/csv/GCA_900186385.1.csv
     ```
     Here: 
     - `-path_host_proteom` is the proteom fasta file containing a set of proteins for the host organism, which is the Prodigal software output. 
     - `-path_host_csv` is the input CSV file with output important defense protein hits from Padloc and defense-finder softwares.

     After running the script, you can expect the following output:
     - A CSV formatted file with list of host CAS proteins included in analysis.

     The output directory for this script is `./output/host`.

   - **`pair_prediction.py`**:
     ```bash
     python pair_prediction.py -encoded_pair ./data/samples/GCA_900186385.1_NC_031911.npy -background_X ./data/background/background_X.npy -background_y ./data/background/background_y.npy
     ```
     Here: 
     - `-encoded_pair` is a NumPy file with the encoded pair of phage and host with phage, host proteins, and their taxonomy.
     - `-background_X` and `-background_y` are NumPy files required to prepare classifier decision space visualization. 

     After running the script, you can expect the following output:
     - Decision space of the classifier presenting position of the predicted pair on the background manifold.
     - NumPy file containing dictionary with probability of interaction and model embedding of delivered pair.

     The output directory for this script is `./output/prediction`.

   - **`post_analysis.py`**:
     ```bash
     python post_analysis.py -path_X ./data/samples/GCA_900186385.1_NC_031911.npy -path_background ./data/background/background_X_SHAP.npy -path_phage_proteins ./output/phage/Xp -path_host_proteins ./output/host/Xh
     ```
     Here: 
     - `-path_X` is the encoded representation for the phage and host pair with proteins and taxonomy information. 
     - `-path_background` is a group of encoded phage and host pairs required for Shapley value analysis. 
     - `-path_phage_proteins` and `-path_host_proteins` are output files from the scripts `phage_proteins.py` and `host_proteins.py`. The user needs to edit the paths to these files after running both scripts previously!!!

     After running the script, you can expect the following output:
     - An activation matrix of encoded elements in the phage and host matrix, saved as an HTML visualization.
     - A CSV formatted file with protein pairs of phage and host which are important in the classification decision process.

     The output directory for this script is `./output/postanalysis`.

## Jupyter Notebook Demo

For a comprehensive demonstration of running all the scripts and presenting their outputs, you can use the `sections_demo.ipynb` Jupyter Notebook provided in this repository.

To use the notebook, follow these steps:

1. Ensure that you have Jupyter Notebook installed. If not, you can install it using pip:
   ```
   pip install notebook
   ```

2. Navigate to the directory of the GitHub VirMAD Repository:
    ```
    cd VirMAD
    ```

3. Start the Jupyter Notebook server by running the following command:
    ```
    jupyter notebook
    ```

4. This will open a new tab in your web browser showing the contents of the repository. Navigate to the `sections_demo.ipynb` file and open it.

5. Follow the instructions provided within the notebook to run the scripts and view their outputs.

The notebook provides an interactive way to explore the functionality of the scripts and visualize their outputs.
