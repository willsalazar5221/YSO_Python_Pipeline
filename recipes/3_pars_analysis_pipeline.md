# `3_pars_analysis_pipeline`


**Author: William B. Salazar**


## Description

This recipe guides you through parsing the "pars" files outputted by the SED Fitter. There are multiple functions embedded within the three scripts to help guide interpretation of the fitter's results. We note the three scripts for this pipeline are `extractor_pipeline_v6.py`, `analysis_functions_v2`, and `model_selection_v1.py`. There are multiple functions outputting the same type of diagram, however there are differences in what data is being outputted so it is important to read the docstrings of each function.

## Version History

- `extractor_pipeline_v6.py` (v6) - July 2025 by W. B. Salazar

This version has been streamlined to work directly from the plotting functions. Can also just use these functions directly to view dataframes.

- `analysis_functions_v2.py` (v2) - July 2025 by W. B. Salazar

Primarily HR diagrams to visually plot the YSOs. Added luminosity frequency distribution functions to visualize hypothesis tests on data.

- `model_selection_v1.py` (v1) - July 2025 by W. B. Salazar

First version of selecting the best model using a tree based on a the best chi-squared for the YSO models and the availability of Gaia data. 


Python blocks/lines will be preceded by >>> for clarity.


## INITIAL SETUP

For this procedure, we solely need to run Python scripts and functions. I primarily used Jupyter Notebooks from Anaconda to run this whole process. You will also need the pars files for the model types numbered 01, 02, 16, and 17 (sp_s_i, sp_h_i, spubsmi, and spubhmi, respectively). The functions for plotting will not work with the other models types.

You will also need the isochrones from Haemmerlé et al. (2019). We use the tracks for the ages 0.5 Myr, 1.0 Myr, 2 Myr, 5 Myr, 31.6 Myr. (these files can be found under `scripts_and_files` for your convenience). For reference:

Paper for Haemmerlé et al. (2019):  https://ui.adsabs.harvard.edu/abs/2019A%26A...624A.137H/abstract
Parent directory containing all the isochrone tracks from the paper: https://obswww.unige.ch/Research/evol/tables_PMS/isochrones/

You will also need a master list of stars in your region for the multiple plots to work. This master list is meant to be the csv file created in `1_retrieving_target_YSOs` recipe under `With Gaia Dataset`. We refer to it as `spicy_gaia_match_csv_name` in the second recipe.


## Process

### Background


In order to calculate the weighted mean, we following the following formula to obtain the weights (Povich, et al. 2013):

$$ P_i = P_n \cdot e^{-\chi_i^{2} / 2} $$

,$P_i$ is the weight of a singular value, $P_n$ is a value such that it normalizes the sum of all the $P_i$'s for the star, and $\chi_i^{2}$ is the chi-squared of that value.

We note 

$$ \sum_{n=1}^{i} P_i = 1 $$

,so we can say

$$ 1 = P_n \cdot \sum_{n=1}^{i} e^{-\chi_i^{2} / 2} $$

We clarify for the purposes of our code that $P_n$ applies to every value in the specific subarray, whereas $P_i$ is an entire subarray of values for that star. So, $P_n$ should equal the len(check_list) while $P_i$ should equal the len(line_list) - len(check_list) (check `extractor_pipeline_v6.py` for details on names of the variables).


The weighted mean is computed by the following formula, 

$$ A_v = \sum_{n=1}^{i} P_i \cdot A_{v,i} $$

,where we note $A_v$ is the weighted average.


### Functions

Within the three scripts, there are multiple helper functions to plot the graphs. The docstrings are very detailed and provide insight into the variables as well as outputs. For simplicity, I will outline the main functions I use outright in the example Notebooks. 

#### multi_single_hr_diagram_av_plots

**hr_diagram_and_dust_ext_single(model_type, star_index_or_name, star_names_list)**
<br>
&emsp; Generate and save plots for a selected star comparing IR-only and IR+Gaia model fits, including HR diagrams and dust extinction diagrams.

&emsp; **Parameters:&ensp; model_type &nbsp;: &nbsp;*int***
<br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; Which model family to use, This will 1,2,16, or 17.
<br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; **star_index_or_name &nbsp;: &nbsp;*int or str***
<br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; Either the index of the star in `star_names_list` or the star's name (string). 
<br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; **star_names_list &nbsp;: &nbsp;*pandas.Series***
<br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; A list or Series of star identifiers used for matching to model results.