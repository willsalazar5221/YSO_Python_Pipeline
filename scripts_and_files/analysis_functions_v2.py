import numpy as np
import pandas as pd
from matplotlib.ticker import AutoMinorLocator
from scipy import stats
from extractor_pipeline_v6 import *

#Read in isochrone files from 

log_lum_1,log_temp_1 = np.genfromtxt('Isochr_Z0.0140_Vini0.00_t05.700.dat',unpack=True,usecols=[4,5])

log_lum_2,log_temp_2 = np.genfromtxt('Isochr_Z0.0140_Vini0.00_t06.000.dat',unpack=True,usecols=[4,5])

log_lum_3,log_temp_3 = np.genfromtxt('Isochr_Z0.0140_Vini0.00_t06.300.dat',unpack=True,usecols=[4,5])

log_lum_4,log_temp_4 = np.genfromtxt('Isochr_Z0.0140_Vini0.00_t06.700.dat',unpack=True,usecols=[4,5])

log_lum_5,log_temp_5 = np.genfromtxt('Isochr_Z0.0140_Vini0.00_t07.500.dat',unpack=True,usecols=[4,5])

def merge_funcs(IR_pd, gaia_pd):
    """
    Merge infrared (IR) and Gaia model-fit results into a single DataFrame, 
    including both raw data and weighted/unweighted parameter averages.

    The function:
      1. Merges the IR and Gaia DataFrames on `MIR_NAME`.
      2. Computes weighted and unweighted means (via `multi_param_weight_avg`) for 
         stellar temperature, luminosity, and extinction in each dataset.
      3. Merges these mean values back into the combined DataFrame with clear suffixes 
         (`_x` for IR, `_y` for Gaia).

    Args:
        IR_pd (pandas.DataFrame):  
            DataFrame of infrared model-fit results containing at least:
                - MIR_NAME (str): Star identifier.  
                - star_temp_arr, lum_arr, av_arr (array-like): Arrays of model-fit values.  
                - chi_2_arr (array-like): Chi-squared values for weighting.  
                - n_data (int), n_fits (int): Metadata about fit coverage.  

        gaia_pd (pandas.DataFrame):  
            DataFrame of Gaia model-fit results with the same structure as `IR_pd`.  

    Returns:
        pandas.DataFrame:  
            Merged DataFrame with columns including:
                - All original IR and Gaia columns (suffixes handled by `merge`).  
                - star_temp_avg_x, star_temp_w_avg_x (IR means).  
                - star_temp_avg_y, star_temp_w_avg_y (Gaia means).  
                - av_avg_x, av_w_avg_x (IR), av_avg_y, av_w_avg_y (Gaia).  
                - lum_avg_x, lum_w_avg_x (IR), lum_avg_y, lum_w_avg_y (Gaia).  
    """
    #Merge the two data sets
    df_IR_G = IR_pd.merge(gaia_pd, how='inner', on='MIR_NAME')
    
    #Get the weights and averages
    IR_pd_means = multi_param_weight_avg(IR_pd, ['star_temp_arr', 'lum_arr', 'av_arr'])
    gaia_pd_means = multi_param_weight_avg(gaia_pd, ['star_temp_arr', 'lum_arr', 'av_arr'])
    
    #Merge the two new pandas together
    df_IR_G_means = IR_pd_means.merge(gaia_pd_means, how='inner', on='MIR_NAME')
    
    #Final merged pd with new 
    merged_df = pd.merge(df_IR_G, df_IR_G_means[['MIR_NAME', 'star_temp_avg_x', 'star_temp_w_avg_x', 'star_temp_avg_y', 'star_temp_w_avg_y', \
                                                   'av_avg_x', 'av_w_avg_x', 'av_avg_y', 'av_w_avg_y', \
                                                   'lum_avg_x', 'lum_w_avg_x', 'lum_avg_y', 'lum_w_avg_y']], on='MIR_NAME', how='left')
    return merged_df
    
def hr_diagram_and_dust_ext_single(model_type, star_index_or_name, star_names_list):
    """
    Generate and save plots for a selected star comparing IR-only and IR+Gaia 
    model fits, including HR diagrams and dust extinction diagrams.

    The function:
      - Accepts either a star index or star name.  
      - Loads the appropriate model-fit parameter files (IR and Gaia) depending 
        on `model_type`.  
      - Merges IR and Gaia results with weighted and unweighted averages.  
      - Produces a 2×2 panel of plots:  
          (0,0): IR HR diagram  
          (0,1): Gaia HR diagram  
          (1,0): IR dust extinction  
          (1,1): Gaia dust extinction  
      - Saves the figure as a PNG named `<MIR_NAME>_<model_type>.png`.  

    Args:
        model_type (int):  
            Which model family to use:  
              - 1 -> "sp_s_i" (disk only, sp-s-i)  
              - 2 -> "sp_h_i" (disk only, sp-h-i)  
              - 16 -> "spubsmi" (disk + envelope)  
              - 17 -> "spubhmi" (disk + envelope)  
              
        star_index_or_name (int or str):  
            Either the index of the star in `star_names_list` or the star's name (string).  

        star_names_list (pandas.Series):  
            A list or Series of star identifiers used for matching to model results.   

    Returns:
        None  
            Displays the generated plots and saves them as PNG files in the 
            working directory.  
            Prints status messages about the processed star.  

    Notes:
        - Calls `read_extracted_file` to load model parameter files.  
        - Calls `merge_funcs` to combine IR and Gaia results.  
        - Isochrones (`log_temp_*`, `log_lum_*`) must be available in the 
          namespace for plotting.  
    """

    if type(star_index_or_name) == str:
        star_index_data = star_names_list[star_names_list.iloc[:] == star_index_or_name].index
        star_index = star_index_data[0]
    else:
        star_index = star_index_or_name

    if model_type == 1:
        model_name = 'sp_s_i'
        gaia_df = read_extracted_file('pars_01g04_Gaia.txt')
        ir_df = read_extracted_file('pars_01g04_IR.txt')
    elif model_type == 2:
        model_name = 'sp_h_i'
        gaia_df = read_extracted_file('pars_02g04_Gaia.txt')
        ir_df = read_extracted_file('pars_02g04_IR.txt')
    elif model_type == 16:
        model_name = 'spubsmi'
        gaia_df = read_extracted_file('pars_16g04_Gaia.txt')
        ir_df = read_extracted_file('pars_16g04_IR.txt')
    elif model_type == 17:
        model_name = 'spubhmi'
        gaia_df = read_extracted_file('pars_17g04_Gaia.txt')
        ir_df = read_extracted_file('pars_17g04_IR.txt')

    #Get merged of IR and Gaia parameters with weighted averages and averages
    big_data = merge_funcs(ir_df, gaia_df)
    name_to_check = star_names_list.loc[star_index]
    # Check if the name exists in df_stats
    exists = name_to_check in big_data['MIR_NAME'].values
    
    #true
    if exists == 1:

        #Get index of star in the merged big_data pandas df
        index_data_star_big_data = big_data[big_data['MIR_NAME'] == name_to_check].index
        index_star_big_data = index_data_star_big_data[0]
        
        print('Working on plots for model set {} for star '.format(model_name) + big_data['MIR_NAME'][index_star_big_data] + ':')
        fig, axs = plt.subplots(2,2,figsize=(20,14))
        
        size_avg_mark = 150
        weights_IR = big_data['chi_2_arr_x'][index_star_big_data]
        weights_gaia = big_data['chi_2_arr_y'][index_star_big_data]
        norm_IR = plt.Normalize(vmin=min(weights_IR), vmax=max(weights_IR))
        norm_gaia = plt.Normalize(vmin=min(weights_gaia), vmax=max(weights_gaia))
        
        #First plot (ax[0]) is IR only data HR diagram
        axs[0, 0].scatter(np.log10(big_data['star_temp_arr' + '_x'][index_star_big_data]),\
                          np.log10(big_data['lum_arr_x'][index_star_big_data]), \
                          c=weights_IR, cmap='rainbow', norm=norm_IR, alpha=0.5)
        
        axs[0,0].scatter(np.log10(big_data['star_temp_arr'[:-4] + '_w_avg' + '_x'][index_star_big_data]),\
                       np.log10(big_data['lum_arr'[:-4] + '_w_avg' + '_x'][index_star_big_data]), marker='^',facecolors='none',\
                       color = 'blue',s=size_avg_mark, linewidths=2, label = 'Weighted Average')
        
        axs[0,0].scatter(np.log10(big_data['star_temp_arr'[:-4] + '_avg' + '_x' ][index_star_big_data]),\
                       np.log10(big_data['lum_arr'[:-4] + '_avg' + '_x'][index_star_big_data]), marker='s',facecolors='none',\
                       color = 'black',s = size_avg_mark, linewidths=2, label = 'Average')
        
        #Plot isochrones on IR plot
        axs[0, 0].plot(log_temp_1, log_lum_1,'-',color='red',lw=3)
        axs[0, 0].plot(log_temp_2, log_lum_2,'-',color='red',lw=3)
        axs[0, 0].plot(log_temp_3, log_lum_3,'-',color='red',lw=3)
        axs[0, 0].plot(log_temp_4, log_lum_4,'-',color='red',lw=3)
        axs[0, 0].plot(log_temp_5, log_lum_5,'-',color='red',lw=3)
        
        #Second plot (ax[1]) is with Gaia data HR diagram
        axs[0,1].scatter(np.log10(big_data['star_temp_arr' + '_y'][index_star_big_data]),\
                         np.log10(big_data['lum_arr_y'][index_star_big_data]),\
                         c=weights_gaia, cmap='rainbow', norm=norm_gaia, alpha=0.5)
        
        axs[0,1].scatter(np.log10(big_data['star_temp_arr'[:-4] + '_w_avg' + '_y'][index_star_big_data]),\
                       np.log10(big_data['lum_arr'[:-4] + '_w_avg' + '_y'][index_star_big_data]), marker='^', facecolors='none',\
                       color = 'blue',s = size_avg_mark, linewidths=2, label = 'Weighted Average')
        
        axs[0,1].scatter(np.log10(big_data['star_temp_arr'[:-4] + '_avg' + '_y' ][index_star_big_data]),\
                       np.log10(big_data['lum_arr'[:-4] + '_avg' + '_y'][index_star_big_data]), marker='s',facecolors='none',\
                       color = 'black',s = size_avg_mark, linewidths=2, label = 'Average')
        
        #Plot isochrones on Gaia plot
        axs[0,1].plot(log_temp_1, log_lum_1,'-',color='red',lw=3, label = '0.5 Myr')
        axs[0,1].plot(log_temp_2, log_lum_2,'-',color='red',lw=3, label = '1 Myr')
        axs[0,1].plot(log_temp_3, log_lum_3,'-',color='red',lw=3, label = '2 Myr')
        axs[0,1].plot(log_temp_4, log_lum_4,'-',color='red',lw=3, label = '5 Myr')
        axs[0,1].plot(log_temp_5, log_lum_5,'-',color='red',lw=3, label = '31.6 Myr')
        
        ###### Plot dust extinction plots now
        
        #First plot (axs[1,0]) is IR only data
        
        caxs_1 = axs[1,0].scatter(np.log10(big_data['star_temp_arr' + '_x'][index_star_big_data]),\
                       big_data['av_arr' + '_x'][index_star_big_data], marker='o',\
                       c=weights_IR, cmap='rainbow', norm=norm_IR, alpha=0.5)
        
        axs[1,0].scatter(np.log10(big_data['star_temp_arr'[:-4] + '_w_avg' + '_x'][index_star_big_data]),\
                       big_data['av_arr'[:-4] + '_w_avg' + '_x'][index_star_big_data], marker='^',facecolors='none',\
                       color = 'blue',s=size_avg_mark, linewidths=2, label = 'Weighted Average')
        
        axs[1,0].scatter(np.log10(big_data['star_temp_arr'[:-4] + '_avg' + '_x' ][index_star_big_data]),\
                       big_data['av_arr'[:-4] + '_avg' + '_x'][index_star_big_data], marker='s',facecolors='none',\
                       color = 'black',s = size_avg_mark, linewidths=2, label = 'Average')
        
        #Second plot (axs[1,1]) is with Gaia data
        caxs_2 = axs[1,1].scatter(np.log10(big_data['star_temp_arr' + '_y'][index_star_big_data]), \
                       big_data['av_arr' + '_y'][index_star_big_data], marker='o',\
                       c=weights_gaia, cmap='rainbow', norm=norm_gaia, alpha=0.5)
        
        axs[1,1].scatter(np.log10(big_data['star_temp_arr'[:-4] + '_w_avg' + '_y'][index_star_big_data]),\
                       big_data['av_arr'[:-4] + '_w_avg' + '_y'][index_star_big_data], marker='^', facecolors='none',\
                       color = 'blue',s = size_avg_mark, linewidths=2, label = 'Weighted Average')
        
        axs[1,1].scatter(np.log10(big_data['star_temp_arr'[:-4] + '_avg' + '_y' ][index_star_big_data]),\
                       big_data['av_arr'[:-4] + '_avg' + '_y'][index_star_big_data], marker='s',facecolors='none',\
                       color = 'black',s = size_avg_mark, linewidths=2, label = 'Average')
        
        #Titles
        if model_type == 1:
            fig.suptitle('Disk Only YSO Models (sp-s-i) ' + '- ' + big_data['MIR_NAME'][index_star_big_data])
        elif model_type == 2:
            fig.suptitle('Disk Only YSO Models (sp-h-i)' + '- ' + big_data['MIR_NAME'][index_star_big_data])
        elif model_type == 16:
            fig.suptitle('Disk + Envelope YSO Models (spubsmi)' + '- ' + big_data['MIR_NAME'][index_star_big_data])
        elif model_type == 17:
            fig.suptitle('Disk + Envelope YSO Models (spubhmi)' + '- ' + big_data['MIR_NAME'][index_star_big_data])
            
        #PLot titles
        axs[0,0].set_title('IR Only')
        axs[0,1].set_title('With Gaia')
        fig.supxlabel(r'log$(T_{eff})$ (K)')
        
        #HR plot graphics
        axs[0,0].set_ylabel(r'log$(L/L_{\odot})$', fontsize = 12) 
        
        ##Tick mark stuff
        axs[0,0].xaxis.set_minor_locator(AutoMinorLocator())
        axs[0,0].yaxis.set_minor_locator(AutoMinorLocator())
        axs[0,0].tick_params(axis="x", which = "both", direction="in")
        axs[0,0].tick_params(axis="y",which = "both", direction="in")
        
        axs[0,1].xaxis.set_minor_locator(AutoMinorLocator())
        axs[0,1].yaxis.set_minor_locator(AutoMinorLocator())
        axs[0,1].tick_params(axis="x", which = "both", direction="in")
        axs[0,1].tick_params(axis="y",which = "both", direction="in")
        
        axs[0,0].set_xlim(4.3,3.5)
        axs[0,0].set_ylim(-0.5, 4)
        axs[0,1].set_xlim(4.3,3.5)
        axs[0,1].set_ylim(-0.5, 4)
        # axs[0,0].legend(loc="lower left")
        axs[0,1].legend(loc="upper right")
        
        ##Accounting of YSOs in data sets
        num_stars_ir = big_data['n_fits_x'][index_star_big_data]
        num_stars_gaia = big_data["n_fits_y"][index_star_big_data]
        axs[0,0].text(4.25, -0.2,'YSOs in Sample = {}'.format(num_stars_ir), \
                    bbox = dict(facecolor = 'none', edgecolor='black'), fontsize = 6.6);
        axs[0,1].text(4.25, -0.2,'YSOs in Sample = {}'.format(num_stars_gaia), \
                    bbox = dict(facecolor = 'none', edgecolor='black'), fontsize = 6.6);
        
        #Dust extinction
        axs[1,0].set_ylabel('Dust Extinction', fontsize = 12)
        
        ##Tick mark stuff
        axs[1,0].xaxis.set_minor_locator(AutoMinorLocator())
        axs[1,0].yaxis.set_minor_locator(AutoMinorLocator())
        axs[1,0].tick_params(axis="x", which = "both", direction="in")
        axs[1,0].tick_params(axis="y",which = "both", direction="in")
        
        axs[1,1].xaxis.set_minor_locator(AutoMinorLocator())
        axs[1,1].yaxis.set_minor_locator(AutoMinorLocator())
        axs[1,1].tick_params(axis="x", which = "both", direction="in")
        axs[1,1].tick_params(axis="y",which = "both", direction="in")
        
        axs[1,0].set_xlim(4.3,3.5)
        axs[1,0].set_ylim(-0.5, 20)
        
        axs[1,1].set_xlim(4.3,3.5)
        axs[1,1].set_ylim(-0.5, 20)
        
        axs[1,0].legend(loc="upper left")
        axs[1,1].legend(loc="upper left")
        
        fig.colorbar(caxs_1, ax=axs[1,0], orientation="horizontal")
        fig.colorbar(caxs_2, ax=axs[1,1], orientation="horizontal")
        
        #Accounting of YSOs in data sets
        axs[1,0].text(0.72 , 0.92, 'Number of Valid Models = {}'.format(num_stars_ir), \
                    bbox = dict(facecolor = 'white', edgecolor='black', alpha = 0.6), \
                    fontsize = 6.5, ha='left', va='top', transform=axs[1,0].transAxes);
        axs[1,1].text(0.72 , 0.92,'Number of Valid Models = {}'.format(num_stars_gaia), \
                    bbox = dict(facecolor = 'white', edgecolor='black', alpha = 0.6), \
                    fontsize = 6.5, ha='left', va='top', transform=axs[1,1].transAxes);
        
        fig.tight_layout()
        filename = big_data['MIR_NAME'][index_star_big_data] + '_' + model_name + '.png'
        plt.savefig(filename, bbox_inches='tight', dpi=300)
        plt.show();
    else:
        print(f"'{name_to_check}' does not have {model_name} model available")
        print('--------------------------------------------------------------')

def multi_single_hr_diagram_av_plots(star_index_given, star_names_pd):
    """
    Generate multiple HR diagram and dust extinction plots for a single star 
    using four model types.

    The function iterates over four model types 
    (`[1, 2, 16, 17]`) and calls `hr_diagram_and_dust_ext_single` for each, 
    producing and saving plots for the given star.

    Args:
        star_index_given (int or str):  
            Either the index of the star in `star_names_pd` or the star's name 
            (string), depending on how the downstream function handles it.  

        star_names_pd (pandas.Series):  
            A Series or DataFrame column containing star identifiers, used to 
            match the given star index or name.  

    Returns:
        None  
            Displays and saves plots for each model type.  
            Output files are named according to star identifier and model type.  

    Notes:
        - This function is a wrapper that simply loops over the four model 
          types and delegates plotting to `hr_diagram_and_dust_ext_single`.  
        - Model types are:  
            1 → disk-only (sp-s-i)  
            2 → disk-only (sp-h-i)  
            16 → disk+envelope (spubsmi)  
            17 → disk+envelope (spubhmi)  
    """
    model_types_arr = np.array([1,2,16,17])
    
    for i in range(0,4):
        model_type_i = model_types_arr[i]
        hr_diagram_and_dust_ext_single(model_type = model_type_i, star_index_or_name = star_index_given, star_names_list = star_names_pd)

def hr_diagram_and_dust_ext_region(model_type):
    """
    Create HR diagram and dust extinction plots for an entire YSO model set.

    Loads IR and Gaia parameter files for the specified `model_type`,
    computes weighted averages of stellar parameters, and produces a
    2×2 grid of plots:

    - Top-left: HR diagram (IR only, weighted averages).
    - Top-right: HR diagram (with Gaia, weighted averages).
    - Bottom-left: Dust extinction vs. temperature (IR only).
    - Bottom-right: Dust extinction vs. temperature (with Gaia).

    Each figure is saved to disk as `<model_name>_region.png`.

    Parameters
    ----------
    model_type : int
        Model set identifier:
        - 1  -> sp_s_i
        - 2  -> sp_h_i
        - 16 -> spubsmi
        - 17 -> spubhmi

    Returns
    -------
    None
        Displays the plots and saves them as a PNG file.

    Notes
    -----
    - Relies on `multi_param_weight_avg()` to compute weighted averages.
    - Isochrone arrays (`log_temp_*`, `log_lum_*`) must be available in scope.
    - HR diagram x-axes are reversed (hotter stars to the left).
    """

    if model_type == 1:
        model_name = 'sp_s_i'
        gaia_df = read_extracted_file('pars_01g04_Gaia.txt')
        ir_df = read_extracted_file('pars_01g04_IR.txt')
    elif model_type == 2:
        model_name = 'sp_h_i'
        gaia_df = read_extracted_file('pars_02g04_Gaia.txt')
        ir_df = read_extracted_file('pars_02g04_IR.txt')
    elif model_type == 16:
        model_name = 'spubsmi'
        gaia_df = read_extracted_file('pars_16g04_Gaia.txt')
        ir_df = read_extracted_file('pars_16g04_IR.txt')
    elif model_type == 17:
        model_name = 'spubhmi'
        gaia_df = read_extracted_file('pars_17g04_Gaia.txt')
        ir_df = read_extracted_file('pars_17g04_IR.txt')

    print('Working on HR diagram for model set {}:'.format(model_name))
    
    df_IR_pars = multi_param_weight_avg(ir_df, ['star_temp_arr', 'lum_arr', 'av_arr'])
    df_gaia_pars = multi_param_weight_avg(gaia_df, ['star_temp_arr', 'lum_arr', 'av_arr'])

    fig, axs = plt.subplots(2,2,figsize=(20,14))
    
    #First plot (ax[0]) is IR only data HR diagram
    axs[0, 0].scatter(np.log10(df_IR_pars['star_temp_w_avg']),\
                      np.log10(df_IR_pars['lum_w_avg']), alpha=0.5)
    
    #Plot isochrones on IR plot
    axs[0, 0].plot(log_temp_1, log_lum_1,'-',color='red',lw=2)
    axs[0, 0].plot(log_temp_2, log_lum_2,'-',color='red',lw=2)
    axs[0, 0].plot(log_temp_3, log_lum_3,'-',color='red',lw=2)
    axs[0, 0].plot(log_temp_4, log_lum_4,'-',color='red',lw=2)
    axs[0, 0].plot(log_temp_5, log_lum_5,'-',color='red',lw=2)
    
    #Second plot (ax[1]) is with Gaia data HR diagram
    axs[0,1].scatter(np.log10(df_gaia_pars['star_temp_w_avg']), \
                     np.log10(df_gaia_pars['lum_w_avg']), alpha=0.5)
    
    #Plot isochrones on Gaia plot
    axs[0,1].plot(log_temp_1, log_lum_1,'-',color='red',lw=2, label = '0.5 Myr')
    axs[0,1].plot(log_temp_2, log_lum_2,'-',color='red',lw=2, label = '1 Myr')
    axs[0,1].plot(log_temp_3, log_lum_3,'-',color='red',lw=2, label = '2 Myr')
    axs[0,1].plot(log_temp_4, log_lum_4,'-',color='red',lw=2, label = '5 Myr')
    axs[0,1].plot(log_temp_5, log_lum_5,'-',color='red',lw=2, label = '31.6 Myr')
    
    ###### Plot dust extinction plots now
    
    #First plot (axs[1,0]) is IR only data
    caxs_1 = axs[1,0].scatter(np.log10(df_IR_pars['star_temp_w_avg']), \
                              df_IR_pars['av_w_avg'], alpha=0.5)
    
    #Second plot (axs[1,1]) is with Gaia data
    caxs_2 = axs[1,1].scatter(np.log10(df_gaia_pars['star_temp_w_avg']), \
                              df_gaia_pars['av_w_avg'], alpha=0.5)
    
    #Titles
    if model_type == 1:
        fig.suptitle('Disk Only YSO Models (sp-s-i)')
    elif model_type == 2:
        fig.suptitle('Disk Only YSO Models (sp-h-i)')
    elif model_type == 16:
        fig.suptitle('Disk + Envelope YSO Models (spubsmi)')
    elif model_type == 17:
        fig.suptitle('Disk + Envelope YSO Models (spubhmi)')
        
    #PLot titles
    axs[0,0].set_title('IR Only')
    axs[0,1].set_title('With Gaia')
    fig.supxlabel(r'log$(T_{eff})$ (K)')
    
    #HR plot graphics
    axs[0,0].set_ylabel(r'log$(L/L_{\odot})$', fontsize = 12) 
    
    ##Tick mark stuff
    axs[0,0].xaxis.set_minor_locator(AutoMinorLocator())
    axs[0,0].yaxis.set_minor_locator(AutoMinorLocator())
    axs[0,0].tick_params(axis="x", which = "both", direction="in")
    axs[0,0].tick_params(axis="y",which = "both", direction="in")
    
    axs[0,1].xaxis.set_minor_locator(AutoMinorLocator())
    axs[0,1].yaxis.set_minor_locator(AutoMinorLocator())
    axs[0,1].tick_params(axis="x", which = "both", direction="in")
    axs[0,1].tick_params(axis="y",which = "both", direction="in")
    
    axs[0,0].set_xlim(4.3,3.5)
    axs[0,0].set_ylim(-0.5, 4)
    axs[0,1].set_xlim(4.3,3.5)
    axs[0,1].set_ylim(-0.5, 4)
    # axs[0,0].legend(loc="lower left")
    axs[0,1].legend(loc="upper right")
    
    ##Accounting of YSOs in data sets
    num_stars_ir = len(df_IR_pars)
    num_stars_gaia = len(df_gaia_pars)
    
    axs[0,0].text(4.25, -0.2,'YSOs in Sample = {}'.format(num_stars_ir), \
                bbox = dict(facecolor = 'none', edgecolor='black'), fontsize = 6.6);
    axs[0,1].text(4.25, -0.2,'YSOs in Sample = {}'.format(num_stars_gaia), \
                bbox = dict(facecolor = 'none', edgecolor='black'), fontsize = 6.6);
    
    #Dust extinction
    axs[1,0].set_ylabel('Dust Extinction', fontsize = 12)
    
    ##Tick mark stuff
    axs[1,0].xaxis.set_minor_locator(AutoMinorLocator())
    axs[1,0].yaxis.set_minor_locator(AutoMinorLocator())
    axs[1,0].tick_params(axis="x", which = "both", direction="in")
    axs[1,0].tick_params(axis="y",which = "both", direction="in")
    
    axs[1,1].xaxis.set_minor_locator(AutoMinorLocator())
    axs[1,1].yaxis.set_minor_locator(AutoMinorLocator())
    axs[1,1].tick_params(axis="x", which = "both", direction="in")
    axs[1,1].tick_params(axis="y",which = "both", direction="in")
    
    axs[1,0].set_xlim(4.3,3.5)
    axs[1,0].set_ylim(-0.5, 20)
    
    axs[1,1].set_xlim(4.3,3.5)
    axs[1,1].set_ylim(-0.5, 20)
    
    fig.tight_layout()
    filename = model_name + '_region' + '.png'
    plt.savefig(filename, bbox_inches='tight', dpi=300)
    plt.show();

    # return None

def multi_region_hr_diagram_av_plots():
    """
    Generate HR diagram and dust extinction plots for all model sets.

    Iterates over the four supported YSO model types and calls
    `hr_diagram_and_dust_ext_region()` for each. Produces one 2×2
    figure per model set (IR-only HR diagram, Gaia HR diagram with
    isochrones, IR dust extinction, Gaia dust extinction).

    Returns
    -------
    None
        Displays each figure and saves it to disk as `<model_name>_region.png`.
    """
    model_types_arr = np.array([1,2,16,17])
    
    for i in range(0,4):
        model_type_i = model_types_arr[i]
        hr_diagram_and_dust_ext_region(model_type = model_type_i)

def lum_freq_distribution_plot(model_combo_type):
    """
    Compare cumulative luminosity distributions between disk-only and 
    disk+envelope YSO model sets using IR-only and Gaia-constrained fits.
    
    This function computes weighted average luminosities for the chosen 
    pair of model types, then generates side-by-side cumulative frequency 
    distributions. A Kolmogorov–Smirnov test is performed to quantify 
    the statistical difference between the two model sets, and results 
    (including sample sizes and p-values) are annotated on the plots.
    
    Parameters
    ----------
    model_combo_type : tuple of (int, int)
        Model set pair to compare. Supported values are:
            (1, 16) → disk vs. disk+envelope (sp_s_i vs. spubsmi)
            (2, 17) → disk vs. disk+envelope (sp_h_i vs. spubhmi)
    
    Returns
    -------
    None
        Displays a matplotlib figure with two subplots:
        - Left: IR-only models
        - Right: Gaia-constrained models
    """

    if model_combo_type == 1_16:
        #First, only disk models
        label_1 = 'sp_s_i'
        gaia_df_disk = read_extracted_file('pars_01g04_Gaia.txt')
        ir_df_disk = read_extracted_file('pars_01g04_IR.txt')
        #Now read in with envelope models
        label_2 = 'spubsmi'
        gaia_df_env = read_extracted_file('pars_16g04_Gaia.txt')
        ir_df_env = read_extracted_file('pars_16g04_IR.txt')
    elif model_combo_type == 2_17:
        #First, only disk models
        label_1 = 'sp_h_i'
        gaia_df_disk = read_extracted_file('pars_02g04_Gaia.txt')
        ir_df_disk = read_extracted_file('pars_02g04_IR.txt')
        #Now read in with envelope models
        label_2 = 'spubhmi'
        gaia_df_env = read_extracted_file('pars_17g04_Gaia.txt')
        ir_df_env = read_extracted_file('pars_17g04_IR.txt')

    print('Working on Frequency of Luminosity Distributions plots for model set {} and {}:'.format(label_1, label_2))

    #Get weighted luminosities for each pars file
    df_IR_pars_disk = multi_param_weight_avg(ir_df_disk, ['lum_arr'])
    df_IR_pars_env = multi_param_weight_avg(ir_df_env, ['lum_arr'])

    df_gaia_pars_disk = multi_param_weight_avg(gaia_df_disk, ['lum_arr'])
    df_gaia_pars_env = multi_param_weight_avg(gaia_df_env, ['lum_arr'])
    
    fig, axs = plt.subplots(1,2,figsize=(14,6))

    axs[0].hist(np.log10(df_IR_pars_disk['lum_w_avg']), bins=len(df_IR_pars_disk), color='blue', histtype='step', alpha =0.5, label = label_1, \
                density = True, cumulative=True)
    
    axs[0].hist(np.log10(df_IR_pars_env['lum_w_avg']), bins=len(df_IR_pars_env), color = 'orange',histtype='step', alpha = 0.5, label= label_2, \
             density = True, cumulative=True)

    #p-value, two-sample Kolmogorov-Smirnov test for goodness of fit
    y_1 = stats.ks_2samp(np.log10(df_IR_pars_disk['lum_w_avg']), np.log10(df_IR_pars_env['lum_w_avg']))
    p_val_1 = y_1[1]
    
    axs[0].set_title('IR Only')
    
    axs[0].text(0.135, 0.75,'YSOs with {} models: {} \nYSOs with {} models: {} \np-Value = {:.2e}'.format(label_1, len(df_IR_pars_disk), label_2, len(df_IR_pars_env), p_val_1),\
             bbox = dict(facecolor = 'white', edgecolor='black', alpha = 0.6), fontsize = 7.5, \
             ha='left', va='top', transform=fig.transFigure);

    #Dashed lines to show percentiles at 10%, 50%, and 90%
    axs[0].axhline(y = 0.1, color = 'r', linestyle = 'dashed', linewidth = 0.5)
    axs[0].axhline(y = 0.5, color = 'r', linestyle = 'dashed', linewidth = 0.5)
    axs[0].axhline(y = 0.9, color = 'r', linestyle = 'dashed', linewidth = 0.5)
    
    axs[0].legend(loc="upper left")

    #plotting gaia plots
    axs[1].hist(np.log10(df_gaia_pars_disk['lum_w_avg']), bins=len(df_gaia_pars_disk), color='blue', histtype='step', alpha =0.5, label = label_1, \
                density = True, cumulative=True)
    
    axs[1].hist(np.log10(df_gaia_pars_env['lum_w_avg']), bins=len(df_gaia_pars_env), color = 'orange',histtype='step', alpha = 0.5, label= label_2, \
             density = True, cumulative=True)

    #p-value, two-sample Kolmogorov-Smirnov test for goodness of fit
    y_2 = stats.ks_2samp(np.log10(df_gaia_pars_disk['lum_w_avg']), np.log10(df_gaia_pars_env['lum_w_avg']))
    p_val_2 = y_2[1]

    #Dashed lines to show percentiles at 10%, 50%, and 90%
    axs[1].axhline(y = 0.1, color = 'r', linestyle = 'dashed', linewidth = 0.5)
    axs[1].axhline(y = 0.5, color = 'r', linestyle = 'dashed', linewidth = 0.5)
    axs[1].axhline(y = 0.9, color = 'r', linestyle = 'dashed', linewidth = 0.5)
    
    axs[1].set_title('With Gaia')
    
    axs[1].text(0.56, 0.75,'YSOs with {} models: {} \nYSOs with {} models: {} \np-Value = {:.2e}'.format(label_1, len(df_gaia_pars_disk), label_2, len(df_gaia_pars_env), p_val_2),\
             bbox = dict(facecolor = 'white', edgecolor='black', alpha = 0.6), fontsize = 7.5, \
             ha='left', va='top', transform=fig.transFigure);
    
    axs[1].legend(loc="upper left")
    
    fig.suptitle('Frequency of Luminosity Distribution of Models ' + label_1 + ' and ' + label_2)
    fig.supxlabel('Luminosity [log$(L/L_{\odot})$]')
    axs[0].set_ylabel('Frequency', fontsize = 12);
    
    # fig.savefig('lum_hist_dist_var' , bbox_inches='tight', dpi=300)
    plt.show();

def multi_lum_freq_distribution_plot():
    """
    Generate luminosity frequency distribution plots for all supported 
    model set pairs.
    
    Iterates over the two supported model combinations ((1,16) and (2,17)) 
    and calls `lum_freq_distribution_plot()` for each, producing cumulative 
    distribution comparisons between disk-only and disk+envelope YSO models.
    
    Returns:
        None
            Displays one figure per model pair, each with IR-only and 
            Gaia-constrained cumulative luminosity distributions.
    """
    model_types_arr = np.array([1_16,2_17])
    
    for i in range(0,2):
        model_type_i = model_types_arr[i]
        lum_freq_distribution_plot(model_combo_type = model_type_i)