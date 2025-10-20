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


## Background

In order to calculate the weighted mean, we following the following formula to obtain the weights (Povich, et al. 2013):

$$ P_i = P_n \cdot e^{-\chi_i^{2} / 2} $$

, $P_i$ is the weight of a singular value, $P_n$ is a value such that it normalizes the sum of all the $P_i$'s for the star, and $\chi_i^{2}$ is the chi-squared of that value.

We note 

$$ \sum_{n=1}^{i} P_i = 1 $$

,so we can say

$$ 1 = P_n \cdot \sum_{n=1}^{i} e^{-\chi_i^{2} / 2} $$

We clarify for the purposes of our code that $P_n$ applies to every value in the specific subarray, whereas $P_i$ is an entire subarray of values for that star. So, $P_n$ should equal the len(check_list) while $P_i$ should equal the len(line_list) - len(check_list) (check `extractor_pipeline_v6.py` for details on names of the variables).


The weighted mean is computed by the following formula, 

$$ A_v = \sum_{n=1}^{i} P_i \cdot A_{v,i} $$

,where we note $A_v$ is the weighted average.


## Functions

Within the three scripts, there are multiple helper functions to plot the graphs. The docstrings are very detailed and provide insight into the variables as well as outputs. For simplicity, I will outline the main functions I use outright in the example Notebooks. 


### ${\color{purple} Single }$

**multi_single_hr_diagram_av_plots(star_index_given, star_names_pd)**

The function iterates over four model types (`[1, 2, 16, 17]`) and calls `hr_diagram_and_dust_ext_single` for each, producing and saving plots for the given star. Generate and save plots for a selected star comparing IR-only and IR+Gaia model fits, including HR diagrams and dust extinction diagrams. Uses a master list to find stars that has both IR and Gaia data points.

&emsp; **Parameters:&ensp; star_index_given &nbsp;: &nbsp;*int or str***
<br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; Either the index of the star in `star_names_pd` or the star's name (string), depending on how the downstream function &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; handles it.
<br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; **star_names_pd &nbsp;: &nbsp;*pandas.Series***
<br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; A Series or DataFrame column containing star identifiers, used to match the given star index or name. 


### ${\color{purple}multi_region_hr_diagram_av_plots}$

**multi_region_hr_diagram_av_plots()**
<br>
Generate HR diagram and dust extinction plots for all model sets. Iterates over the four supported YSO model types and calls `hr_diagram_and_dust_ext_region()` for each. Produces one 2×2 figure per model set (IR-only HR diagram, Gaia HR diagram with isochrones, IR dust extinction, Gaia dust extinction).

&emsp; **Parameters:&ensp; None &nbsp;: &nbsp;*int***
<br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; Take in no parameters as it runs through `hr_diagram_and_dust_ext_region(model_type)` for each of the four model types we analyze. More details in the docstrings.


### ${\color{purple}multi_lum_freq_distribution_plot}$

**multi_lum_freq_distribution_plot()**
<br>
Compare cumulative luminosity distributions between disk-only and disk+envelope YSO model sets using IR-only and Gaia-constrained fits. Iterates over the two supported model combinations ((1,16) and (2,17)) and calls `lum_freq_distribution_plot()` for each, producing cumulative distribution comparisons between disk-only and disk+envelope YSO models. This function computes weighted average luminosities for the chosen pair of model types, then generates side-by-side cumulative frequency distributions. A Kolmogorov–Smirnov test is performed to quantify the statistical difference between the two model sets, and results (including sample sizes and p-values) are annotated on the plots.

&emsp; **Parameters:&ensp; None &nbsp;: &nbsp;*int***
<br>
Take in no parameters as it runs through `lum_freq_distribution_plot(model_combo_type)` for each of the four model types we analyze. More details in the docstrings.


### ${\color{purple}final_model_select}$

**final_model_select(master_list_IR, master_list_gaia, user_cdp)**
<br>
This function reads in the master SPICY catalog cutouts for IR and Gaia sources, along with model parameter files for four model types (1, 2, 16, and 17). It computes model likelihoods via `calc_p_dm_df()`, tags each model set with an identifier, and determines the best-fitting model for each source using `model_tree()`. If Gaia data are available for a star, its model selection is prioritized over IR-only fits, even when the chi-squared value is higher, due to the increased number of data points.


&emsp; **Parameters:&ensp; master_list_IR &nbsp;: &nbsp;*str***
<br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; File path to the IR master catalog (CSV) containing SPICY cutout data.
<br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; **master_list_gaia &nbsp;: &nbsp;*str***
<br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; File path to the Gaia master catalog (CSV) containing SPICY cutout data.
<br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; **user_cdp &nbsp;: &nbsp;*float***
<br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; Critical delta probability threshold used in the model likelihood calculation (Eq. 21 in Robitaille, T. P. 2017). Recorded in `2_data_format_pipeline` under `Next Steps`.


### ${\color{purple}multi_param_weight_avg_tree}$

**multi_param_weight_avg_tree(df, param_list)**
<br>
Compute weighted and unweighted means for multiple parameters. This function applies `weight_mean_tree` to several parameter columns (e.g., extinction, stellar radius, temperature, disk mass, luminosity) and combines the results into a single DataFrame. It preserves source metadata and appends mean statistics for each parameter side by side.

&emsp; **Parameters:&ensp; df &nbsp;: &nbsp;*pandas.DataFrame***
<br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; DataFrame containing model fit results for multiple sources. Must include columns required by `weight_mean_tree(df, col_name)`
<br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; **param_list &nbsp;: &nbsp;*list of str***
<br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; List of parameter column names for which to compute weighted and unweighted means.


### ${\color{purple}hr_diagram_and_dust_ext_region_tree}$

**hr_diagram_and_dust_ext_region_tree(df_pars)**
<br>
Plot HR diagram and dust extinction trends for combined model tree results. This function generates a two-panel figure summarizing the stellar properties from all model combinations (sp_s_i, sp_h_i, spubsmi, spubhmi).

&emsp; **Parameters:&ensp; df_pars &nbsp;: &nbsp;*pandas.DataFrame***
<br>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; DataFrame containing model-fitting results from multiple regions or model combinations.






