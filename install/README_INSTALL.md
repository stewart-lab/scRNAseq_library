To install new package:
* navigate to conda-forge https://conda-forge.org/feedstock-outputs/
* search for package
* click: https://anaconda.org/conda-forge/r-readr 
* run the command: `conda install -c conda-forge r-readr`
* find the version: `conda env export | grep readr`
* add that line to the environment_template.yml file