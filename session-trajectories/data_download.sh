#! /usr/bin/env bash

# Jules GILET <jules.gilet@curie.fr>

# Download of the example datasets used in the trajectory inferrence lab 
# ELIXIR scRNAseq school


wget -nH -nd	http://blood.stemcells.cam.ac.uk/data/nestorowa_corrected_log2_transformed_counts.txt.gz \
	http://blood.stemcells.cam.ac.uk/data/nestorowa_corrected_population_annotation.txt.gz \
	ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81682/suppl/GSE81682%5FHTSeq%5Fcounts%2Etxt%2Egz

mkdir data
gunzip *
mv *.txt data/

exit 0
