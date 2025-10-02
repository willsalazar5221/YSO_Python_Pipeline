# `retrieving_target_YSOs`

**Author: William B. Salazar**

## Description

This recipe guides you through creating the data file formatted for the SED fitting process. Refer to [Data format](https://sedfitter.readthedocs.io/en/stable/data.html) for further information on the data format. 

In addition, when introducing visible light magnitudes, we will use the text file produced by the `Perform aperture photometry on the MIPS 24 µm mosaic at each YSO position` step in Povich's [sedfitting_procedure_ysos](https://github.com/mattpovich/sedfitting-ysos/blob/master/recipes/sedfitting_procedure_ysos.md). 

## Version History

- Current (v3) - July 2025 by W. B. Salazar

This version has been streamlined such that there is only two lines to run for the process. Previously modified from a Jupyter notebook version. Python blocks/lines will be preceded by >>> for clarity.


## INITIAL SETUP

For this procedure, we solely need to run Python scripts and functions. I primarily used Jupyter Notebooks from Anaconda to run this whole process. You will need your initial csv file of the cutout of your region with the columns described in `retrieving_target_YSOs`. 

For the second part, you will need the text file outputted in the described Povich repository. The file name will most likely end in `_ugos_24um`. Add the `.txt` extension before running it through this pipeline. You can remove it again for the SED fitter. There should be no loss of information even if the txt extension is removed or added.

## Process

### Background

As background, we are calculating the necessary fluxes from the magnitudes. This performed in the following manner. In order to calculate the flux using the zero-point flux and magnitudes, we utilized the following equation

$$ Flux = ZP_v \cdot 1000 \cdot 10^{-\frac{m}{2.5}} $$

,where $ZP_v$ is the zero-point flux and $m$ is the magnitude. To calculate the uncertainties, we use the following equation 

$$ \Delta Flux = \sqrt{(-ZP_v \cdot 0.4 \cdot \ln{10} \cdot 10^{-\frac{m}{2.5}} \cdot \Delta m)^{2}} $$

Uncertainties for the Gaia magnitudes are given a 5% lower limit and calculated using the 'over error' columns. More information is detailed in lines 173 - 181 in the code.

### ${\color{red}IR-Only \space Datasets \space Dataset}$

First, we focus on using the infrared magnitudes first. For this run, we need the csv file as described in `retrieving_target_YSOs`. Once you have this, run the line

>>> 
`df_IR = flag_txt_file_IR('region_name.csv', 'data_file_name')`

where we define the following variables from the line:
- "df_IR" - You may rename it to any variable name. This is to display the outputted text file as a dataframe. Note it doesn't save this dataframe to a text file
- "region_name.csv" - Rename it to the csv file name, keeping the csv extension
- "data_file_name" - Rename it to something appropriate to remember it, such as "your_region_name_IR". It is a text file without the extension so it can pass through the SED fitter correctly.


### ${\color{purple}With \space Gaia \space Dataset}$

<i>Note</i>: This step should be done after running the SED Fitter once. Code for running that is provided in `NEXT STEPS` below.

Now, after running the SED fitter once, if there is available visible light data provided by Gaia, we will use it. Here, we will use the following function by running the line

>>>
`df_gaia = flag_txt_file_Gaia(ugos_24um_file_name, spicy_gaia_match_csv_name, data_file_name)`

where we define the following variables from the line:
- "ugos_24um_file_name" is the file created in the pre-SED fitting process. We will resuse it from the IR run
- "spicy_gaia_match_csv_name" is the second csv file we created in `retrieving_target_YSOs` recipe under `With Gaia Dataset`
- "data_file_name" is simply the name you chose for the text file to run in the SED Fitter. There is no need to add a .txt extension, however you may add it. It should not affect the SED fitting process. From here, you can use the file to run again in the SED Fitter.


## NEXT STEPS:

For both runs of the SED fitting process, you will use the created text files. This process can be found in Povich's repository with the recipe titled `sedfitting_procedure_ysos`, under the [`Fit the 1-24 µm SEDs of YSO candidates with R17 YSO model sets.` section] (https://github.com/mattpovich/sedfitting-ysos/blob/master/recipes/sedfitting_procedure_ysos.md). For me, the block of code is modified to the following for the IR run:

>>>
```
from astropy import units as u
from fit_multi import fit_multiyso
extinction24_file = 'ex_law_mc.par'
models_topdir = <path to your directory containing all the YSO model subdirectories> #Mine is called 'ysomodels'

data = 'data_glimpse+ysoc_ugos_24um'

#We only use IR in the first run
filters = ['2J', '2H', '2K', 'I1', 'I2', 'I3', 'I4', 'M1'] 
apertures = [3., 3., 3., 3., 3., 3., 3., 7.] * u.arcsec

fit_multiyso(data, filters, apertures, models_topdir,
	   extinction_file=extinction24_file,
	   n_data_min=3,  output_format=('F',6),
	   distance_range=[2.9, 3.3] * u.kpc,
	   av_range=[0.,50.], cpd=4.)
```

<b>IMPORTANT</b>: Change the string name in the variable `data` to the name of the text file you outputted above. For the `distance_range`, put in the distance range for your region. Note that `cpd` will be a number you must choose and remember for the next recipe.

For the Gaia run of the same region, we modify the following two line, where the rest are unchanged:

>>>

```
#We introduce Gaia data now
filters = ['GBP', 'GRP', '2J', '2H', '2K', 'I1', 'I2', 'I3', 'I4', 'M1'] 
apertures = [3., 3., 3., 3., 3., 3., 3., 3., 3. 7.] * u.arcsec
```

Once you run the SED fitting process, you now have pars files. Now proceed to the second pipeline recipe `pars_analysis_pipeline`.
