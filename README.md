# Deriving the onset and offset times of planning units from acoustic and articulatory measurements
[Joe Rodd <sup>1,2</sup>](https://www.mpi.nl/people/rodd-joe)
[Hans Rutger Bosker <sup>1,3</sup>](https://www.mpi.nl/people/bosker-hans-rutger)
[Louis ten Bosch <sup>2,1,3</sup>](https://www.ru.nl/english/people/bosch-l-ten/)
[Mirjam Ernestus <sup>2,1</sup>](http://www.mirjamernestus.nl)

1) [Max Planck Institute for Psycholinguistics](https://www.mpi.nl), P.O. Box 310, 6500 AH, Nijmegen, The Netherlands
2) Radboud University, [Centre for Language Studies](https://www.ru.nl/cls/), P.O. Box 9103, 6500 HD, Nijmegen, The Netherlands
3) Radboud University, [Donders Institute for Brain, Cognition and Behaviour](https://www.ru.nl/donders/), P.O. Box 9104, 6500 HE, Nijmegen, The Netherlands

> Many psycholinguistic models of speech sequence planning make claims about the onset and offset times of planning units, such as words, syllables, and phonemes. These predictions typically go untested, however, since psycholinguists have assumed that the temporal dynamics of the speech signal is a poor index of the temporal dynamics of the underlying speech planning process. It is argued that this problem is tractable, and two simple metrics that derive planning unit onset and offset times from the acoustic signal and articulatographic data are presented and validated.

This article is in press at the Journal of the Acoustical Society of America. [Here's our postprint](https://github.com/joerodd/deriving-onset-offset-times-planning-units/raw/master/manuscript/Manuscript.pdf).

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

