# iPSC_Rewind

This repository contains scripts needed to produce all graphs and images contained in figures for the paper:

Jain N, et al. Retrospective identification of intrinsic factors that mark pluripotency potential in rare somatic cells. (2023)
https://www.biorxiv.org/content/10.1101/2023.02.10.527870v1

For any questions, please contact NJ (naveen.jain.debate@gmail.com) and/or the corresponding author for this paper, AR (arjunrajlab@gmail.com).

In order to reproduce all graphs and images in this paper:

Download and install R v4.0.5.

Download all files in the structure provided for this repository from Dropbox: https://www.dropbox.com/sh/ulu6728tcp49dv2/AAAPwLYQiVLloH_JL38lvTj6a?dl=0 (original manuscript), https://www.dropbox.com/sh/zz958910t4fkj9w/AAAgTVwO5yAKZ1TpSQVfV6Qga?dl=0 (revised manuscript). Your project directory should have subdirectories extractedData, extractionScripts, figures, imageProcessing, plots, plotScripts, and rawData.

Edit directory variables (i.e. "homeDirectory") variables to reflect your local file paths for this project repository and the associated downloaded files. Set "plotDirectory" to reflect your desired location for finished plots.

Additional annotation/information available in individual R scripts. Scripts to produce extractedData from rawData can be found in extractionScripts. Scripts to produce plots from extractedData can be found in plotScripts.
