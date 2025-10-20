import numpy as np
import pandas as pd
from analysis_functions_v2 import *
from extractor_pipeline_v6 import *

log_lum_1,log_temp_1 = np.genfromtxt('Isochr_Z0.0140_Vini0.00_t05.700.dat',unpack=True,usecols=[4,5])

log_lum_2,log_temp_2 = np.genfromtxt('Isochr_Z0.0140_Vini0.00_t06.000.dat',unpack=True,usecols=[4,5])

log_lum_3,log_temp_3 = np.genfromtxt('Isochr_Z0.0140_Vini0.00_t06.300.dat',unpack=True,usecols=[4,5])

log_lum_4,log_temp_4 = np.genfromtxt('Isochr_Z0.0140_Vini0.00_t06.700.dat',unpack=True,usecols=[4,5])

log_lum_5,log_temp_5 = np.genfromtxt('Isochr_Z0.0140_Vini0.00_t07.500.dat',unpack=True,usecols=[4,5])

def final_model_select(master_list_IR, master_list_gaia, user_cdp):
    """
    Select the final best-fit YSO model for each source across IR-only and Gaia datasets.
    
    This function reads in the master SPICY catalog cutouts for IR and Gaia sources,
    along with model parameter files for four model types (1, 2, 16, and 17).  
    It computes model likelihoods via `calc_p_dm_df()`, tags each model set with an 
    identifier, and determines the best-fitting model for each source using `model_tree()`.  
    If Gaia data are available for a star, its model selection is prioritized over 
    IR-only fits, even when the chi-squared value is higher, due to the increased number of data points.
    
    Finally, a `.reg` region file is generated via `make_region_file()` for visualization in DS9,
    and a combined DataFrame of all selected models is returned.
    
    Args:
        master_list_IR : str
            File path to the IR master catalog (CSV) containing SPICY cutout data.
        master_list_gaia : str
            File path to the Gaia master catalog (CSV) containing SPICY cutout data.
        user_cdp : float
            Critical delta probability threshold used in the model likelihood calculation 
            (Eq. 21 in Robitaille, T. P. 2017).
    
    Returns
    -------
    pandas.DataFrame
        A DataFrame containing the final selected model parameters for each source, 
        including model flags and data type indicators.
    """

    #Read in master lists from SPICY catalog cutouts
    #### MUST STILL CONTAIN 'SPITZER' COLUMN ####
    master_df_IR = pd.read_csv(master_list_IR) #IR master list
    master_df_gaia = pd.read_csv(master_list_gaia) #gaia master list
    
    #Drops unnecessary parts of SPITZER name
    master_df_IR['Spitzer'] = master_df_IR['Spitzer'].str.replace(r'SSTGLMC ', '')
    master_df_gaia['Spitzer'] = master_df_gaia['Spitzer'].str.replace(r'SSTGLMC ', '')

    #Use read_extracted_file from analysis_functions_v2 to read in pars files for each region
    #Must follow the same naming convention to read the files in
    df_IR_01 = read_extracted_file('pars_01g04_IR.txt')
    df_IR_02 = read_extracted_file('pars_02g04_IR.txt')
    df_IR_16 = read_extracted_file('pars_16g04_IR.txt')
    df_IR_17 = read_extracted_file('pars_17g04_IR.txt')

    df_gaia_01 = read_extracted_file('pars_01g04_Gaia.txt')
    df_gaia_02 = read_extracted_file('pars_02g04_Gaia.txt')
    df_gaia_16 = read_extracted_file('pars_16g04_Gaia.txt')
    df_gaia_17 = read_extracted_file('pars_17g04_Gaia.txt') 

    #Calculate probability ffrom Eq. 21 using function below
    df_pdm_IR_01 = calc_p_dm_df(df_IR_01, user_cdp, 1)
    df_pdm_IR_02 = calc_p_dm_df(df_IR_02, user_cdp, 2)
    df_pdm_IR_16 = calc_p_dm_df(df_IR_16, user_cdp, 16)
    df_pdm_IR_17 = calc_p_dm_df(df_IR_17, user_cdp, 17)

    #Gaia data
    df_pdm_gaia_01 = calc_p_dm_df(df_gaia_01, user_cdp, 1)
    df_pdm_gaia_02 = calc_p_dm_df(df_gaia_02, user_cdp, 2)
    df_pdm_gaia_16 = calc_p_dm_df(df_gaia_16, user_cdp, 16)
    df_pdm_gaia_17 = calc_p_dm_df(df_gaia_17, user_cdp, 17)

    #Initialize an column for each df that has a simple column called 'Flag'. Then attached said flag to the whole dataframe
    df_pdm_IR_01['Model_Flag'] = 1
    df_pdm_IR_02['Model_Flag'] = 2
    df_pdm_IR_16['Model_Flag'] = 16
    df_pdm_IR_17['Model_Flag'] = 17

    #Gaia data
    df_pdm_gaia_01['Model_Flag'] = 1
    df_pdm_gaia_02['Model_Flag'] = 2
    df_pdm_gaia_16['Model_Flag'] = 16
    df_pdm_gaia_17['Model_Flag'] = 17

    #Now get the best model selection from two trees
    IR_tree_df = model_tree(df_pdm_IR_01, df_pdm_IR_02, df_pdm_IR_16, df_pdm_IR_17, master_df_IR) #IR
    gaia_tree_df = model_tree(df_pdm_gaia_01, df_pdm_gaia_02, df_pdm_gaia_16, df_pdm_gaia_17, master_df_gaia) #Gaia

    #Now we prioritze Gaia available stars. This is because even if they have a worse chi-squared, more data points are available
    #We always chose Gaia if available, otherwise, we pick the best IR model

    #Make a new flag for the star based on data points used
    IR_tree_df['Type_Flag'] = 'ir'
    gaia_tree_df['Type_Flag'] = 'gaia'

    #Since Gaia YSOs is a subset of SPICY, loop over IR cutout
    final_df = pd.DataFrame() #initialize df with all the stars

    #For every star in the master list, loop over its name
    for k in range(0,len(master_df_IR)):
        saved_row_2 = None
        name_check_2 = master_df_IR['Spitzer'][k]

        #Check if exists in gaia tree
        exists_gaia = name_check_2 in gaia_tree_df['MIR_NAME'].values
        #Check if it exists in IR tree
        exists_IR = name_check_2 in IR_tree_df['MIR_NAME'].values
        
        if exists_gaia == 1:
            #Get row and use that as saved_row_2
            saved_row_2 = gaia_tree_df[gaia_tree_df['MIR_NAME'] == name_check_2].iloc[[0]]
        elif exists_IR == 1:
            #Since it didn't exist in the Gaia, we choose the best one in IR
            saved_row_2 = IR_tree_df[IR_tree_df['MIR_NAME'] == name_check_2].iloc[[0]]
        final_df = pd.concat([final_df, saved_row_2], ignore_index=True)

    #Make a region file based on the 'Type_Flag' columns
    make_region_file(final_df)

    return final_df

def calc_p_dm_df(pars_file_df, cdp_choice, par_type):
    """
    Compute the model probability (P_DM) for each YSO based on chi-squared filtering criteria.
    
    This function applies Equation 20 from Robitaille 2017 to prune
    poor-fitting models from each star’s model array. It then calculates the number of
    “good” fits (n_good) remaining after this statistical cut and derives the normalized
    probability P_DM for each source by dividing n_good by the total number of models
    in the corresponding model grid.
    
    The function works on parameter files loaded via `read_extracted_file()` and is
    used as an intermediate step in selecting best-fit models.
    
    Args:
    pars_file_df : pandas.DataFrame
        DataFrame containing model parameters for each YSO, including arrays of
        chi-squared values and fitted parameters (e.g., luminosity, temperature, AV).
    cdp_choice : float
        The critical delta chi-squared multiplier (Δχ² threshold) used for filtering
        model fits, corresponding to the user-defined confidence level.
    par_type : int
        Identifier for the model set being processed:
            1 or 2 -> Disk-only models (10,000 total models)
            16     -> Disk + Envelope models (40,000 total models)
            17     -> Disk + Envelope models (80,000 total models)
    
    Returns:
    pandas.DataFrame
        The filtered DataFrame containing only accepted model fits, with additional
        columns for:
            - 'best_chi2' : minimum chi-squared value per source
            - 'n_good' : number of models passing the cut
            - 'P_DM' : normalized probability of detection or model likelihood
    """
    #To eliminate the fact that there would be more columns based on different models, I'm going to keep columns
    #I know I will want that are avaible in all models.

    #Note: This particular way creates a copy of the dataframe, which may raise warnings if you want to add a column
    #that's why we do .copy()
    columns_to_keep = ['MIR_NAME', 'n_data', 'n_fits', 'init_n_fits', 'chi_2_arr', 'av_arr',\
                       'star_rad_arr', 'star_temp_arr','d_mass_arr', 'lum_arr']
    pars_file_df = pars_file_df[columns_to_keep].copy()

    #Make a new column that saves the best chi2 for that particular model, which is always the first value in the array
    pars_file_df['best_chi2'] = np.nan
    for i in range(0, len(pars_file_df)):
        pars_file_df.loc[i, 'best_chi2'] = pars_file_df.loc[i, 'chi_2_arr'][0]

    #Use Eqn. 20 from Robitaille, T. P. 2017, A&A, 600, A11 to cut out models below this equation
    #Note logic is always switched around in python since we grab true values and leave false ones alone
    for j in range(0, len(pars_file_df)):
        index = np.where((pars_file_df['chi_2_arr'][j] - pars_file_df['best_chi2'][j]) > (cdp_choice * pars_file_df['n_data'][j]))
        pars_file_df.at[j, 'chi_2_arr'] = np.delete(pars_file_df['chi_2_arr'][j], index)
        pars_file_df.at[j, 'av_arr'] = np.delete(pars_file_df['av_arr'][j], index)
        pars_file_df.at[j, 'star_rad_arr'] = np.delete(pars_file_df['star_rad_arr'][j], index)
        pars_file_df.at[j, 'star_temp_arr'] = np.delete(pars_file_df['star_temp_arr'][j], index)
        pars_file_df.at[j, 'd_mass_arr'] = np.delete(pars_file_df['d_mass_arr'][j], index)
        pars_file_df.at[j, 'lum_arr'] = np.delete(pars_file_df['lum_arr'][j], index)

    #Make another column for the number of good fits now with the cuts
    pars_file_df['n_good'] = np.nan
    for k in range(0, len(pars_file_df)):
        pars_file_df.loc[k, 'n_good'] = float(np.size(pars_file_df['chi_2_arr'][k]))

    pars_file_df['P_DM'] = np.nan

    if par_type == 1 or par_type == 2: 
        N_tot_models = 10000
    elif par_type == 16:
        N_tot_models = 40000
    elif par_type == 17:
        N_tot_models = 80000
        
    for ii in range(0, len(pars_file_df)):
        pars_file_df.loc[ii, 'P_DM'] = pars_file_df['n_good'][ii] / N_tot_models
        
    return pars_file_df

def model_tree(model_01_df, model_02_df, model_16_df, model_17_df, master_list):
    """
    Select the best-fitting model for each source across all model grids.
    
    This function compares multiple model sets (01, 02, 16, 17) for each YSO
    and identifies the model with the lowest chi-squared value. The comparison
    is performed by matching source names from the provided master catalog
    against each model DataFrame and keeping only the best fit.
    
    Args:
    model_01_df : pandas.DataFrame
        Model parameter DataFrame for model type 01 (sp_s_i).
    model_02_df : pandas.DataFrame
        Model parameter DataFrame for model type 02 (sp_h_i).
    model_16_df : pandas.DataFrame
        Model parameter DataFrame for model type 16 (spubsmi).
    model_17_df : pandas.DataFrame
        Model parameter DataFrame for model type 17 (spubhmi).
    master_list : pandas.DataFrame
        Master catalog containing all sources, including the 'Spitzer' column
        used for matching against model DataFrames.
    
    Returns:
    pandas.DataFrame
        A DataFrame containing one row per source, representing the model
        with the minimum chi-squared value across all four model sets.
    """

    overall_df = pd.DataFrame() #initialize df with all the stars

    #For every star in the master list, loop over its name
    for i in range(0,len(master_list)):
        name_check = master_list['Spitzer'][i]
        saved_row = None
        for item in [model_01_df, model_02_df, model_16_df, model_17_df]:
            # Check if the name exists in df_stats
            exists = name_check in item['MIR_NAME'].values
            #If the star is actually in the pars df, check if its the lowest chi-squared
            if exists == 1:
                selected_row_df = item[item['MIR_NAME'] == name_check].iloc[[0]] 
                if saved_row is None:
                    saved_row = selected_row_df
                else:
                    if saved_row['best_chi2'].values[0] > selected_row_df['best_chi2'].values[0]:
                        saved_row = selected_row_df
        overall_df = pd.concat([overall_df, saved_row], ignore_index=True)
    return overall_df

def make_region_file(final_tree_df):
    """
    Generate a DS9 region file marking sources from the final model selection.
    
    This function creates a `.reg` file compatible with SAOImage DS9, marking 
    each source’s galactic coordinates with color-coded circular regions 
    depending on whether the source was classified as 'gaia' or 'ir'. The 
    regions are written using DS9 version 4.1 format.
    
    Args:
    final_tree_df : pandas.DataFrame
        DataFrame containing the final best-fit models for each source. Must 
        include columns:
            - 'MIR_NAME' : str
                Source identifier containing galactic longitude and latitude 
                encoded as part of the name.
            - 'Type_Flag' : str
                Classification flag, either 'gaia' or 'ir', used to set region 
                color and circle size.
    
    Outputs:
    whole_region_circles.reg : file
        A DS9 region file saved in the current working directory, where:
            - Gaia sources are shown as large yellow circles (10″ radius).
            - IR-only sources are shown as smaller red circles (6″ radius).
    
    Notes: 
        - The galactic coordinates are parsed directly from the `MIR_NAME` string
        by slicing its numeric substrings corresponding to longitude and latitude.
    """
    
    # Header for DS9 region file
    header = '''# Region file format: DS9 version 4.1
    global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
    galactic
    '''
    #Initialize txt string
    txt_string = ''
    
    for i in range(len(final_tree_df)):
        name_star_Spitzer = final_tree_df['MIR_NAME'][i]
        flag = final_tree_df['Type_Flag'][i] # either 'ir' or 'gaia'
    
        #Remove the 0th character (the 'G')
        name_star_trim = name_star_Spitzer[1:]
        
        #Slice characters 1–8 and 9–16 from the trimmed string
        l_str = name_star_trim[0:8]   # Characters 1 to 8 (from original string)
        b_str = name_star_trim[8:16]  # Characters 9 to 16 (from original string)
        
        l = float(l_str)
        b = float(b_str)
    
        # Format the values
        if flag == 'gaia':
            txt_string += f'circle({l:.4f},{b:.4f},10.000") # color=yellow width=4 tag={{Gaia}}\n'
        elif flag == 'ir':
            txt_string += f'circle({l:.4f},{b:.4f},6.000") # color=red width=3 tag={{IRonly}}\n'
    
    full_text = header + txt_string
    
    data_file = open('whole_region_circles' + '.reg', 'w') #write your file name in the open function
    n = data_file.write(full_text)
    data_file.close()

def weight_mean_tree(df, col_name):
    """
    Compute simple and weighted means of a specified parameter for each source.
    This function calculates both the arithmetic and weighted mean of a given
    parameter column (e.g., temperature, luminosity, or extinction) for each star
    in the input DataFrame. The weighting is based on the chi-squared values
    associated with each model fit, using the probability weights defined by 
    exp(-chi-squared/2). If only one model exists for a star, the simple mean is used 
    for both the regular and weighted mean.
    
    Args
    ----
    df : pandas.DataFrame
        DataFrame containing model fit results for multiple stars. Must include:
            - 'MIR_NAME' : str
                Identifier for each source.
            - 'chi_2_arr' : array-like
                Array of χ² values corresponding to model fits for each source.
            - 'n_data', 'n_good', 'Model_Flag', 'Type_Flag' : 
                Metadata columns preserved in the output.
            - `col_name` : array-like
                Parameter array (e.g., 'star_temp_arr', 'lum_arr', etc.) 
                for which means will be computed.
    
    col_name : str
        Name of the column in `df` containing the array of parameter values to 
        average.
    
    Returns
    -------
    pandas.DataFrame
        A new DataFrame containing the following columns:
            - 'MIR_NAME', 'n_data', 'n_good', 'Model_Flag', 'Type_Flag'
            - `<col>_avg`   : unweighted mean of the parameter
            - `<col>_w_avg` : weighted mean of the parameter based on χ² values
    """
    #Initialize empty dataframe

    rows = []
    
    for i in range(0, len(df)):
        star_name = df['MIR_NAME'][i]
        chi_sq_arr = df['chi_2_arr'][i]

        #Have to make sure the actual column has many model. If not, don't do a weighted mean,
        if np.size(df['chi_2_arr'][i]) == 1:
            new_row = {'MIR_NAME': star_name, 'n_data': df['n_data'][i], 'n_good': df['n_good'][i], \
                       'Model_Flag': df['Model_Flag'][i], 'Type_Flag': df['Type_Flag'][i], \
                       col_name.replace('arr', 'avg'): float(df[col_name][i]), col_name.replace('arr', 'w_avg'): float(df[col_name][i])}
        else:
            param_mean = df[col_name][i].mean() #parameter mean for that star
            #weighted means
            #Here we produce the normalization constants P_n for each star
            p_n = 1 / sum(np.exp((-1*(chi_sq_arr))/2))
            #print(p_n)
            #Here we compute the weights for each value
            p_i = p_n * np.exp((-1*(chi_sq_arr))/2)
            #print(p_i)
            
            param_weighted_mean = sum(p_i * df[col_name][i])
            
            # Create a new row and then append it to our rows list
            new_row = {'MIR_NAME': star_name, 'n_data': df['n_data'][i], 'n_good': df['n_good'][i], \
                       'Model_Flag': df['Model_Flag'][i], 'Type_Flag': df['Type_Flag'][i], \
                       col_name.replace('arr', 'avg'): param_mean, col_name.replace('arr', 'w_avg'): param_weighted_mean}
        rows.append(new_row)

    result_df = pd.DataFrame(rows)
    return result_df

def multi_param_weight_avg_tree(df, param_list):
    """
    Compute weighted and unweighted means for multiple parameters.
    This function applies `weight_mean_tree` to several parameter columns
    (e.g., extinction, stellar radius, temperature, disk mass, luminosity) 
    and combines the results into a single DataFrame. It preserves source
    metadata and appends mean statistics for each parameter side by side.
    
    Args:
    df : pandas.DataFrame
        DataFrame containing model fit results for multiple sources. Must include
        columns required by `weight_mean_tree`, such as:
            - 'MIR_NAME'
            - 'chi_2_arr'
            - 'n_data'
            - 'n_good'
            - 'Model_Flag'
            - 'Type_Flag'
            - Each parameter column listed in `param_list`.
    
    param_list : list of str
        List of parameter column names (e.g., 
        ['av_arr', 'star_rad_arr', 'star_temp_arr', 'd_mass_arr', 'lum_arr'])
        for which to compute weighted and unweighted means.
    
    Returns:
    pandas.DataFrame
        Combined DataFrame including:
            - Metadata columns from the first parameter processed.
            - Two new columns for each parameter:
                `<param>_avg`   : unweighted mean
                `<param>_w_avg` : weighted mean
    """
    for k in range(0, len(param_list)):
        if k == 0:
            big_df = weight_mean_tree(df, param_list[k])          
        else:
            new_df = weight_mean_tree(df, param_list[k])
            avg_cols = new_df.iloc[:, 5:7]
            big_df = pd.concat([big_df, avg_cols],axis=1)
            
    return big_df

def hr_diagram_and_dust_ext_region_tree(df_pars):
    """
    Plot HR diagram and dust extinction trends for combined model tree results.
    This function generates a two-panel figure summarizing the stellar properties
    from all model combinations (sp_s_i, sp_h_i, spubsmi, spubhmi). The top panel
    shows a Hertzsprung–Russell diagram using weighted-average stellar temperatures
    and luminosities, while the bottom panel shows dust extinction (A_V) as a
    function of effective temperature. Each model type is marked with a distinct
    symbol. Evolutionary isochrones are overplotted for age reference.
    
    Args:
    df_pars : pandas.DataFrame
        DataFrame containing model-fitting results from multiple regions or model
        combinations. Must include columns for:
            - 'star_temp_arr'
            - 'lum_arr'
            - 'av_arr'
            - 'Model_Flag' (1, 2, 16, 17)
            - plus all fields required by `multi_param_weight_avg_tree`.
    
    Returns:
    None
        Displays and saves the resulting two-panel figure ('full_region.png')
        showing:
            1. HR diagram (log T_eff vs. log L/L_⊙) with model type markers.
            2. Dust extinction A_V vs. log T_eff for the same sources.
    
    Notes:
        - Isochrones (e.g., 0.5, 1, 2, 5, 31.6 Myr) must be preloaded as arrays:
              `log_temp_1`, `log_lum_1`, ..., `log_temp_5`, `log_lum_5`.
        - The function automatically saves the figure as 'full_region.png'
          and displays it using matplotlib.
    """
    
    df_complete_tree = multi_param_weight_avg_tree(df_pars, ['star_temp_arr', 'lum_arr', 'av_arr'])

    fig, axs = plt.subplots(2,1,figsize=(14,18))

    axs[0].scatter(np.log10(df_complete_tree['star_temp_w_avg'][df_complete_tree['Model_Flag'] == 1]), \
                   np.log10(df_complete_tree['lum_w_avg'][df_complete_tree['Model_Flag'] == 1]), \
                   c='blue', marker='+', label='sp_s_i', alpha=0.5)

    axs[0].scatter(np.log10(df_complete_tree['star_temp_w_avg'][df_complete_tree['Model_Flag'] == 2]), \
                   np.log10(df_complete_tree['lum_w_avg'][df_complete_tree['Model_Flag'] == 2]), \
                   c='blue', marker='x', label='sp_h_i', alpha=0.5)

    axs[0].scatter(np.log10(df_complete_tree['star_temp_w_avg'][df_complete_tree['Model_Flag'] == 16]), \
                   np.log10(df_complete_tree['lum_w_avg'][df_complete_tree['Model_Flag'] == 16]), \
                   c='blue', marker='s', label='spubsmi', alpha=0.5)

    axs[0].scatter(np.log10(df_complete_tree['star_temp_w_avg'][df_complete_tree['Model_Flag'] == 17]), \
                   np.log10(df_complete_tree['lum_w_avg'][df_complete_tree['Model_Flag'] == 17]), \
                   c='blue', marker='D', label='spubhmi', alpha=0.5)
    
    #Plot isochrones
    axs[0].plot(log_temp_1, log_lum_1,'-',color='red',lw=2, label = '0.5 Myr')
    axs[0].plot(log_temp_2, log_lum_2,'-',color='red',lw=2, label = '1 Myr')
    axs[0].plot(log_temp_3, log_lum_3,'-',color='red',lw=2, label = '2 Myr')
    axs[0].plot(log_temp_4, log_lum_4,'-',color='red',lw=2, label = '5 Myr')
    axs[0].plot(log_temp_5, log_lum_5,'-',color='red',lw=2, label = '31.6 Myr')
    
    ###### Plot dust extinction plots now
    # axs[1].scatter(np.log10(df_complete_tree['star_temp_w_avg']), df_complete_tree['av_w_avg'], alpha=0.5)

    axs[1].scatter(np.log10(df_complete_tree['star_temp_w_avg'][df_complete_tree['Model_Flag'] == 1]), \
                   df_complete_tree['av_w_avg'][df_complete_tree['Model_Flag'] == 1], \
                   c='blue', marker='+', label='sp_s_i', alpha=0.5)
    
    axs[1].scatter(np.log10(df_complete_tree['star_temp_w_avg'][df_complete_tree['Model_Flag'] == 2]), \
                   df_complete_tree['av_w_avg'][df_complete_tree['Model_Flag'] == 2], \
                   c='blue', marker='x', label='sp_h_i', alpha=0.5)

    axs[1].scatter(np.log10(df_complete_tree['star_temp_w_avg'][df_complete_tree['Model_Flag'] == 16]), \
                   df_complete_tree['av_w_avg'][df_complete_tree['Model_Flag'] == 16], \
                   c='blue', marker='s', label='spubsmi', alpha=0.5)

    axs[1].scatter(np.log10(df_complete_tree['star_temp_w_avg'][df_complete_tree['Model_Flag'] == 17]), \
                   df_complete_tree['av_w_avg'][df_complete_tree['Model_Flag'] == 17], \
                   c='blue', marker='D', label='spubhmi', alpha=0.5)
    
    #Titles
    # fig.suptitle('Disk Only YSO Models (sp-s-i)')
        
    #PLot titles
    fig.supxlabel(r'log$(T_{eff})$ (K)')
    
    #HR plot graphics
    axs[0].set_ylabel(r'log$(L/L_{\odot})$', fontsize = 16) 
    
    ##Tick mark stuff
    axs[0].xaxis.set_minor_locator(AutoMinorLocator())
    axs[0].yaxis.set_minor_locator(AutoMinorLocator())
    axs[0].tick_params(axis="x", which = "both", direction="in")
    axs[0].tick_params(axis="y",which = "both", direction="in")
    
    axs[0].set_xlim(4.3,3.5)
    axs[0].set_ylim(-0.5, 4)
    axs[0].legend(loc="upper right")
    
    ##Accounting of YSOs in data sets
    num_stars = len(df_complete_tree)
    
    axs[0].text(4.25, -0.2,'YSOs in Sample = {}'.format(num_stars), \
                bbox = dict(facecolor = 'none', edgecolor='black'), fontsize = 12);
    
    #Dust extinction
    axs[1].set_ylabel('Dust Extinction', fontsize = 16)
    
    ##Tick mark stuff
    axs[1].xaxis.set_minor_locator(AutoMinorLocator())
    axs[1].yaxis.set_minor_locator(AutoMinorLocator())
    axs[1].tick_params(axis="x", which = "both", direction="in")
    axs[1].tick_params(axis="y",which = "both", direction="in")
    
    axs[1].set_xlim(4.3,3.5)
    axs[1].set_ylim(-0.5, 20)
    axs[1].legend(loc="upper left")
    
    fig.tight_layout()
    filename = 'full_region' + '.png'
    plt.savefig(filename, bbox_inches='tight', dpi=300)
    plt.show();