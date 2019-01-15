# Deriving the onset and offset times of planning units from acoustic and articulatory measurements
Joe Rodd 1,2 
Hans Rutger Bosker 1,3
Louis ten Bosch 2,1,3
Mirjam Ernestus 2,1

1) Max Planck Institute for Psycholinguistics, P.O. Box 310, 6500 AH, Nijmegen, The Netherlands
2) Radboud University, Centre for Language Studies, P.O. Box 9103, 6500 HD, Nijmegen, The Netherlands
3) Radboud University, Donders Institute for Brain, Cognition and Be- haviour, P.O. Box 9104, 6500 HE, Nijmegen, The Netherlands

# This repository

This repository contains code pertaining to the above publication. This code was run on Mac OS and on the Linux cluster at the MPI, which runs SuSE 13.2.

There are some absolute paths in some of the scripts, search for "/data/clusterfs/pol/joerod/" and correct as appropriate for your environement.

# Prerequisites
To replicate, you'll need parts of the mgnu0 dataset (available from http://www.mngu0.org/ after registration):
- Day1 basic EMA data, head corrected and unnormalised
- Day1 audio data processed for experiments
- Day1 transcriptions, Festival utterances and ESPS label files

Required software:
- R, various libraries
- HTK

# Run order
1) `process_mngu0.R` is a preprocessing script that prepares the raw mngu0 data for use. This is optimised for our cluster, the multidplyr call should probably be modified to suit your environment.
2) `prepare_for_annotation.R` performs further preprocessing before hand annotation of stable periods from the principal components of the articulatory signals
3) `annotate/ui.R`, `annotate/server.R` and `annotate/global.R` are a shiny app that facilitates the hand annotation
4) `gather_acoustic_onsets_and_offsets_for_annotated.R` does what it says on the tin
5) `target_comparison_on_grid_multi_within.R` performs the analysis of the articulatory data
6) `correlation_analysis.R` performs the correlation analysis reported in the article.

