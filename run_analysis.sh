#!/bin/bash

####################################################################
# Modify this variables by adding the path to your singularity image
container='PATH TO THE CONTAINER'
####################################################################

mkdir -p procdata QC PCA PCA/notexcluded exp1
mkdir -p DE DE_controls DE_controls/DESeq2_reports exp1/exp1_DEreports/DESeq2_report DE/DESeq2_reports
mkdir -p DE/Max_Rad21_Meis1_allsgs/DESeq2_reports DE/Max_Rad21_Meis1_allsgs
mkdir -p NET modules CELLproj TFnet_DoTscore module_plot varfigures 
mkdir -p addfigures
mkdir -p QC_exp3 PCA_exp3 PCA_exp3/notexcluded
mkdir -p DE_exp3 DE_exp3/DESeq2_reports DE_controls_exp3 DE_controls_exp3/DESeq2_reports
mkdir -p NET_exp3 DEint_exp3

singularity exec ${container} R -e "knitr::spin('01_TFnet_QC.R')"
singularity exec ${container} R -e "knitr::spin('02_TFnet_PCA.R')"
singularity exec ${container} R -e "knitr::spin('03_TFnet_DE.R')"
singularity exec ${container} R -e "knitr::spin('04_TFnet_NET.R')"
singularity exec ${container} R -e "knitr::spin('05_TFnet_modules.R')"

singularity exec ${container} /usr/local/bin/jupytext --to notebook 06_TFnet_CELLproj.py
singularity exec ${container} /usr/local/bin/jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=0 06_TFnet_CELLproj.ipynb

singularity exec ${container} /usr/local/bin/jupytext --to notebook 07_TFnet_DoTscore.py
singularity exec ${container} /usr/local/bin/jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=0 07_TFnet_DoTscore.ipynb

singularity exec ${container} /usr/local/bin/jupytext --to notebook 08_TFnet_moduleplot.py
singularity exec ${container} /usr/local/bin/jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=0 08_TFnet_moduleplot.ipynb

singularity exec ${container} /usr/local/bin/jupytext  --to notebook 10_TFnet_varfigures.py
singularity exec ${container} /usr/local/bin/jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=0 10_TFnet_varfigures.ipynb

singularity exec ${container} R -e "knitr::spin('11_TFnet_exp1.R')"

singularity exec ${container} R -e "knitr::spin('12_sc5_addfigures.R')"

singularity exec ${container} R -e "knitr::spin('20_exp3_QC.R')"
singularity exec ${container} R -e "knitr::spin('21_exp3_PCA.R')"
singularity exec ${container} R -e "knitr::spin('22_exp3_DE.R')"
singularity exec ${container} R -e "knitr::spin('23_exp3_DEint.R')"
singularity exec ${container} R -e "knitr::spin('24_exp3_NET.R')"
