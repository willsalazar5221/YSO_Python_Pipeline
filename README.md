# sedfitting-analysis-ysos

Python scripts creating data files integrating multi-survey photometry (2MASS, Spitzer, Gaia) for spectral energy distribution modeling of YSOs in the SPICY catalog using Tom Robitaille's [sedfitter](https://github.com/astrofrog/sedfitter). In addition, we analyze the parameters returned from the SED fitting process from the parameter files (which we will refer to as pars files). We also set a basic tree to determine the "best fitted" model. 

Further procedures and details for going throught the SED fitting process can be found in Matt Povich's [sedfitting-ysos](https://github.com/mattpovich/sedfitting-ysos/tree/master) repository. 

In the `scripts_and_files` directory, we find all the scripts needed for the entire analysis.

In the 'recipes' directory, we find one initial target recipe and two pipeline recipes for the SED fitter process and analysis. The target recipe, `1_retrieving_target_YSOs` details how to get the target list of stars for your region as well as the csv file format needed for the following pipeline. This should be visited first. The second recipe, `2_data_format_pipeline` allows you to create the data file in the correct data file format. More information on the [SED fitting website](https://sedfitter.readthedocs.io/en/stable/data.html). The third recipe, `3_pars_analysis_pipeline`, describes how we will use the pars files for analysis and what outputs to expect. We will also use the tree for a final analaysis of our targets. A region file is outputted so you can visualize the YSO locations on a region file in DS9.

These pipelines have been tested on `python 3.8` and `3.9` with no adverse effects.

## Getting started 

* If you do not have csv file with the magnitudes from the SPICY catalog, start with the first recipe titled `retrieving_target_YSOs`. It is also a good idea to review it for the next pipeline.
* If you have already have a csv file with the target YSO's, proceed to the first pipeline recipe titled `data_format_pipeline`.

