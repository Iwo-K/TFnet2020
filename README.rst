.. role:: pyth(code)
  :language: python

TFnet analysis - Kucinski et al. 2020
=====================================

This repository contains R and python scripts used to perform analysis published in the Kucinski et al. 2020 paper: "Cooperation of Transcription Factors Associated with Diverse Lineages Governs Haematopoietic Progenitor States".

Container
---------

To make the data analysis easy to reproduce we provide a pre-built `Singularity <https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps>`_ container with R and python installations and necessary packages. The container image is available here: http://128.232.224.252/TFnet2020/TFnet2020_container.simg


How to use
----------

To reproduce the analysis presented in the paper:

1. Clone and enter the repository:

.. code-block:: text

    git clone https://github.com/Iwo-K/TFnet2020
    cd TFnet2020

2. Download the container from: http://128.232.224.252/TFnet2020/TFnet2020_container.simg


3. Modify the container variable in the run_analysis.sh script to point to the singularity image


4. Download the data from: http://128.232.224.252/TFnet2020/TFnet_data.tar.gz, uncompress and move to the directory: data/TFnet_data


5. Run the analysis:

.. code-block:: text

   ./run_analysis.sh

This script will create several directories to hold output files from respective analysis script, e.g. figures or processed data. Additionally, the .R scripts will be converted into .html reports and .py scripts into Jupyter/IPython notebooks, to allow easy inspection.

