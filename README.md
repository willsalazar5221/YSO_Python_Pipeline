# sedfitting-analysis-ysos

Python scripts creating data files integrating multi-survey photometry (2MASS, Spitzer, Gaia) for spectral energy distribution modeling of YSOs in the SPICY catalog using Tom Robitaille's [sedfitter](https://github.com/astrofrog/sedfitter). In addition, we analyze the parameters returned from the SED fitting process from the parameter files (which we will refer to as pars files). We also set a basic tree to determine the "best fitted" model. 

Further procedures and details for going throught the SED fitting process can be found in Matt Povich's [sedfitting-ysos](https://github.com/mattpovich/sedfitting-ysos/tree/master) repository. 

In the 'python_scripts' directory, we find all the scripts needed for the entire analysis.

In the 'recipes' directory, we find one initial target recipe and two pipeline recipes for the SED fitter process and analysis. The recipe will detail how to get the target list of stars for your region as well as the csv file format needed for the following pipeline. The first pipeline is to get the data file in the correct format. More information on the [SED fitting website](https://sedfitter.readthedocs.io/en/stable/data.html). The second recipe describes how we will use the pars files for analysis and what outputs to expect. We will also use the tree for a final analaysis of our targets.

These pipelines have been tested on `python 3.8` and `3.9` with no adverse effects.

## Getting started 

* If you do not have csv file with the magnitudes from the SPICY catalog, start with the first recipe titled `retrieving_target_YSOs`
* If you have already have a csv file with the target YSO's, proceed to the first pipeline recipe titled `data_format_pipeline` 


