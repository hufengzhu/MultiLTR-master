# MultiLTR-master

## Introduction
MultiLTR: An ensemble machine learning toolkit with combined features for LTR-RT recognition and classification, 

which takes advantage of six machine learning algorithms and three types of features.

## Directory
[Feature] Used for encoding, which contain four feature encoding schemes, Kmer, reverse complement kmer (RCKmer), and electron-ion interaction pseudopotentials of trinucleotide (PseEIIP). 

[Model] 14 different binary classification models for different classification levels, including distinguishing LTR-RTs and non-LTR-RTs, superfamilies and lineages. 

## Usage
load_stage.py: Used to distinguish LTR-RTs and non-LTR-RTs, superfamilies and lineages.
result_metric.py: Evaluate the prediction results of MultiLTR.

O.sativa_Inpactor2_library.fasta: O.sativa genome TEs annotations were obtained by combining Inpactor2 with Repeat Mask.
