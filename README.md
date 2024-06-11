# VirMAD
## Interpretable Machine Learning for Bacterial Defense and Viral Antidefense Pair Inference
Python scripts to identify viral antidefenses against specific bacterial defenses.

VirMAD works by:
1) encoding viral and bacterial proteins using a BERT-architecture protein language model
2) encoding bacterial and viral taxonomy using CNNs on CGR nucleotide inputs
3) using a CNN to predict infection outcome based on coupled bacterial-virus input tensors and experimental cross-infection matrix derived labels
4) identify putative defense and antidefense pairs as co-activated viral proteins and bacterial defenses using a Shapley value analysis.

Data can be downloaded from [google drive](https://drive.google.com/drive/folders/1F6MdBPKILQesGUqA8AS2x6j-p4sVctfu?usp=drive_link)
Downloaded data should be placed in the cloned folder in a "data" subfolder

## Comparative Genomics
The folder Comparative_Genomics contains julia and bash scripts to carry out a comparative analysis on a group of host and viral genomes, while utilizing experimental cross-infection matrix data. The scripts are tailored for the specific case of the Nahant cross infection matrix, but can be trivially adjusted to any other genome resolved cross-infectio matrix. For questions and clarifications about the comparative genomics analysis please email Shaul Pollak at shaul.pollak.pasternak@univie.ac.at
