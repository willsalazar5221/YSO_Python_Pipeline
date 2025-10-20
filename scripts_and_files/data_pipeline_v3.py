import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd

def read_file_and_flux_calc_IR(file_name):
    """
    Reads a CSV file containing infrared magnitudes and errors for stars, calculates the corresponding
    fluxes and flux uncertainties using standard zero-point values, and returns cleaned and structured data.

    Args:
        file_name (str): Path to the CSV file containing magnitude and error columns for IR bands.
                         Expected column names include:
                         - Magnitudes: 'mag3_6', 'mag4_5', 'mag5_8', 'mag8_0', 'j_synth', 'h_synth', 'k_synth'
                         - Errors: 'e_mag3_6', 'e_mag4_5', 'e_mag5_8', 'e_mag8_0'
                         - Identifier: 'Spitzer'

    Returns:
        tuple:
            df1 (pd.DataFrame): DataFrame containing magnitude values for the 7 IR bands (Spitzer and 2MASS).
            df2 (pd.DataFrame): DataFrame containing magnitude error values for the 4 Spitzer IRAC bands.
            df (pd.DataFrame): DataFrame with computed fluxes (in mJy, millijanskys) and flux uncertainties, with original
                               magnitude and error columns removed, formatted for downstream analysis.
            columns1 (list of str): Names of the magnitude columns used in the flux calculation.
            columns2 (list of str): Names of the magnitude error columns used for flux error propagation.
    """
    
    #First, we create columns of magnitude values for IR only values
    #Order of filters magnitudes: Spitzer/IRAC channel 1, channel 2, channel 3, channel 4
    #2MASS J, H, Ks
    
    columns1 = ['mag3_6','mag4_5','mag5_8','mag8_0','j_synth','h_synth','k_synth'] #Magnitude Values
    columns2 = ['e_mag3_6','e_mag4_5','e_mag5_8','e_mag8_0'] #Errors on magnitudes. Note 2MASS errors are unavailable

    #read into pandas dataframe
    df =  pd.read_csv(file_name)

    df.columns = df.columns.str.replace(' ', '') #Gets rid of unnecessary spaces in header names

    #Can rework code not to extract two dfs but this is legacy code kept from when I first was
    #learning Pandas and did this, so I never changed it
    df1 = pd.DataFrame(df, columns = columns1) #Extracts only the magnitudes
    df2 = pd.DataFrame(df, columns = columns2) #Extracts only part of the error on magnitudes

    #Drops unnecessary parts of SPITZER name
    df['Spitzer'] = df['Spitzer'].str.replace(r'SSTGLMC ', '')
    
    #This block creates the Flux and Flux error columns
    #Zero-point fluxes in Jansky from SVO Filter Profile Service (order of fluxes same as order of magnitudes)
    ZPv = [274.53, 177.66, 113.56, 63.70, 1594.00, 1024.00, 666.80]

    #Flux and Flux error columns
    Flux = ['Flux36','Flux45','Flux58','Flux80', 'FluxJ','FluxH','FluxK',]
    Fluxerr = ['Flux36err','Flux45err','Flux58err','Flux80err','FluxJerr','FluxHerr','FluxKerr']

    for i in range(len(ZPv)):
        df[Flux[i]] = 1000 * ZPv[i] * 10**(-df1[columns1[i]] / 2.5)

    #Here is where the error creation starts
    #Only does first 4 columns
    for j in range(0,4):
        df[Fluxerr[j]] = np.sqrt((-ZPv[j] * 0.4 * np.log(10) * 10**(-df1[columns1[j]] * 0.4) * df2[columns2[j]] )**2)*1000

    #Created 'conservative' errors on 2MASS bands' flux
    df[Fluxerr[4]] = 0.05 * df[Flux[4]]
    df[Fluxerr[5]] = 0.05 * df[Flux[5]]
    df[Fluxerr[6]] = 0.05 * df[Flux[6]]
    
    #Run this line for all data files. Run only after obtaining Flux and Errors
    df = df.drop(columns = ['mag3_6','mag4_5','mag5_8','mag8_0','j_synth','h_synth','k_synth',\
                            'e_mag3_6','e_mag4_5','e_mag5_8','e_mag8_0'])

    #This cell formats the data to be in scientific notation and replaces NaN values.
    ## This is only to display it here properly. In the last cell it formats it for the file
    pd.options.display.float_format = '{:,.4e}'.format
    df = df.replace(np.nan, 0.) #Replaces NaN values with 0.

    return df1, df2, df, columns1, columns2

def flag_txt_file_IR(file_name, data_file_name):
    """
    Reads a CSV file of infrared magnitudes, computes fluxes and flags indicating data availability, 
    and writes a structured text file containing positional data, flags, fluxes, and their uncertainties.

    Args:
        file_name (str): Path to the CSV file containing infrared magnitudes and errors for target stars.
                         Must be in a format compatible with `read_file_and_flux_calc_IR`, including 
                         'Spitzer', 'l', and 'b' columns.
        data_file_name (str): Desired name (without extension) of the output .txt file to be written.
                              The file will contain tabular data for use in Robitaille et al. (2017) 
                              SED fitting routine

    Returns:
        df (pd.DataFrame): DataFrame containing the fluxes, flux errors, and corresponding availability flags
                           for each target. Flags are binary indicators (1 if valid, 0 if missing).
        (saves a txt file to your computer)
    """
    df1, df2, df, columns1, columns2 = read_file_and_flux_calc_IR(file_name)
    
    #Flag values necessary to indicate if fluxes are available. See 'Data format' for more details

    #Flag1 - mag3_6, Flag 2 - mag4_5, Flag 3 - mag5_8, Flag 4 - mag8_0
    #Flag 5 - j_synth, Flag 6 - h_synth, Flag 7 - k_synth

    Flag = ['Flag1','Flag2','Flag3','Flag4','Flag5','Flag6','Flag7']

    for j in range(len(Flag)):
        df[Flag[j]] = np.isfinite(df1[columns1[j]]).astype(int)
        
    # Initialize the string to write to a text file. This format allows for great control over output
    #Flags are set this way since I wanted it to output a file in the following order:
    #2MASS JHK bands, Spitzer/IRAC (I1, I2, I3, I4) bands
    #Little note: notice for b, we have >8,.5f for the formatting. This is to accound for negative longitudes
    txt_string = ''
    
    for i in range(0,len(df)):
        txt_string = txt_string + str(df.Spitzer[i]) \
        + '  ' + str(f"{df.l[i]:.5f}") + '  ' + str(f"{df.b[i]:>8,.5f}") \
        + '  ' + str(df.Flag5[i]) + '  ' + str(df.Flag6[i]) + '  ' + str(df.Flag7[i]) \
        + '  ' + str(df.Flag1[i]) + '  ' + str(df.Flag2[i]) + '  ' + str(df.Flag3[i]) + '  ' + str(df.Flag4[i]) \
        + '  ' + str(f"{df.FluxJ[i]:.4e}") + '  ' + str(f"{df.FluxJerr[i]:.4e}") \
        + '  ' + str(f"{df.FluxH[i]:.4e}") + '  ' + str(f"{df.FluxHerr[i]:.4e}") \
        + '  ' + str(f"{df.FluxK[i]:.4e}") + '  ' + str(f"{df.FluxKerr[i]:.4e}") \
        + '  ' + str(f"{df.Flux36[i]:.4e}") + '  ' + str(f"{df.Flux36err[i]:.4e}") \
        + '  ' + str(f"{df.Flux45[i]:.4e}") + '  ' + str(f"{df.Flux45err[i]:.4e}") \
        + '  ' + str(f"{df.Flux58[i]:.4e}") + '  ' + str(f"{df.Flux58err[i]:.4e}") \
        + '  ' + str(f"{df.Flux80[i]:.4e}") + '  ' + str(f"{df.Flux80err[i]:.4e}") + '\n'
    
    data_file = open(str(data_file_name), 'w') #write your file name in the open function
    n = data_file.write(txt_string)
    data_file.close()

    return df

def flag_txt_file_Gaia(ugos_24um_file_name, spicy_gaia_match_csv_name, data_file_name):
    """
    Merges Spitzer-Gaia matched sources with supplemental 24μm data, computes Gaia fluxes and 
    uncertainties, sets flags indicating valid measurements, and writes a structured .txt file 
    combining positional, flag, and photometric data.

    Args:
        ugos_24um_file_name (str): Path to a fixed-width format (.txt) file containing 24μm data, 
                                   where the first column corresponds to Spitzer source IDs.
        spicy_gaia_match_csv_name (str): Path to a CSV file containing matched Gaia photometry 
                                         for Spitzer sources. Must include:
                                         - 'Spitzer'
                                         - 'phot_bp_mean_mag', 'phot_rp_mean_mag'
                                         - 'phot_bp_mean_flux_over_error', 'phot_rp_mean_flux_over_error'
        data_file_name (str): Name of the output text file (no `.txt` extension needed) to be written.

    Returns:
        merged_df (pd.DataFrame): DataFrame containing the merged dataset with Gaia fluxes, 
                                  errors, and data availability flags, along with all input 24μm data.
    """

    #Read in both files. Notice the difference in pandas functions due to ugos_24um file format
    df_csv =  pd.read_csv(spicy_gaia_match_csv_name)
    df_ugos_24um = pd.read_fwf(ugos_24um_file_name, header=None)

    #Must rename first column Spitzer for matching later
    df_ugos_24um.rename(columns={list(df_ugos_24um)[0]:'Spitzer'}, inplace=True)

    #Zero-point fluxes in Jansky from SVO Filter Profile Service (order of fluxes same as order of magnitudes)
    ZPv = [3552.01, 2554.95]
    
    #Flux and Flux error columns
    Flux = ['Flux_bp', 'Flux_rp']
    Fluxerr = ['Flux_bperr', 'Flux_rperr']
    
    df_csv[Flux[0]] = 1000 * ZPv[0] * 10**(-df_csv['phot_bp_mean_mag'] / 2.5)
    df_csv[Flux[1]] = 1000 * ZPv[1] * 10**(-df_csv['phot_rp_mean_mag'] / 2.5)
    
    #Now create errors for Gaia Fluxes
    #We set condition such that any value less than 20 in ...flux_over_error is set to 20 (5% uncertainty lower limit)
    df_csv.loc[df_csv['phot_bp_mean_flux_over_error'] > 20 , 'phot_bp_mean_flux_over_error'] = 20
    df_csv.loc[df_csv['phot_rp_mean_flux_over_error'] > 20 , 'phot_rp_mean_flux_over_error'] = 20
    
    #Next, we actually create errors on Gaia fluxes from ...flux_over_error
    #Have to explicitly call the columns instead of df2 since df2 keeps old errors
    df_csv[Fluxerr[0]] = (1 / df_csv['phot_bp_mean_flux_over_error']) * df_csv[Flux[0]]
    df_csv[Fluxerr[1]] = (1 / df_csv['phot_rp_mean_flux_over_error']) * df_csv[Flux[1]]

    #Flag8 - Gaia bp, Flag9 - Gaia rp
    #Chosen as 8 and 9 due to legacy flag numbering in initial versions
    Flag = ['Flag8','Flag9']

    df_csv[Flag[0]] = np.isfinite(df_csv['phot_bp_mean_mag']).astype(int)
    df_csv[Flag[1]] = np.isfinite(df_csv['phot_rp_mean_mag']).astype(int)

    #Run this line for all data files. Run only after obtaining Flux and Errors
    df_csv = df_csv.drop(columns=['phot_bp_mean_flux_over_error','phot_bp_mean_mag', 'phot_rp_mean_flux_over_error', 'phot_rp_mean_mag'])
    df_csv['Spitzer'] = df_csv['Spitzer'].str.replace(r'SSTGLMC ', '')

    merged_df = pd.merge(df_csv, df_ugos_24um, on='Spitzer')
    merged_df = merged_df.replace(np.nan, 0.) #Replaces NaN values with 0.

    #Drop rows where Gaia for some reason never actually has data
    filtered_df = merged_df[(merged_df['Flag8'] == 0) & (merged_df['Flag9'] == 0)]
    merged_df = merged_df.drop(filtered_df.index)
    merged_df = merged_df.reset_index(drop=True)

    txt_string = ''

    for i in range(0,len(merged_df)):
        txt_string = txt_string + str(merged_df.Spitzer[i]) \
        + '  ' + str(f"{merged_df[1][i]:.5f}") + '  ' + str(f"{merged_df[2][i]:.5f}") \
        + '  ' + str(merged_df.Flag8[i]) + '  ' + str(merged_df.Flag9[i]) + '  ' + str(merged_df[3][i]) \
        + '  ' + str(merged_df[4][i]) + '  ' + str(merged_df[5][i]) + '  ' + str(merged_df[6][i]) \
        + '  ' + str(merged_df[7][i]) + '  ' + str(merged_df[8][i]) + '  ' + str(merged_df[9][i]) + '  ' + str(merged_df[10][i]) \
        + '  ' + str(f"{merged_df.Flux_bp[i]:.4e}") + '  ' + str(f"{merged_df.Flux_bperr[i]:.4e}") \
        + '  ' + str(f"{merged_df.Flux_rp[i]:.4e}") + '  ' + str(f"{merged_df.Flux_rperr[i]:.4e}") \
        + '  ' + str(f"{merged_df[11][i]:.4e}") + '  ' + str(f"{merged_df[12][i]:.4e}") \
        + '  ' + str(f"{merged_df[13][i]:.4e}") + '  ' + str(f"{merged_df[14][i]:.4e}") \
        + '  ' + str(f"{merged_df[15][i]:.4e}") + '  ' + str(f"{merged_df[16][i]:.4e}") \
        + '  ' + str(f"{merged_df[17][i]:.4e}") + '  ' + str(f"{merged_df[18][i]:.4e}") \
        + '  ' + str(f"{merged_df[19][i]:.4e}") + '  ' + str(f"{merged_df[20][i]:.4e}") \
        + '  ' + str(f"{merged_df[21][i]:.4e}") + '  ' + str(f"{merged_df[22][i]:.4e}") \
        + '  ' + str(f"{merged_df[23][i]:.4e}") + '  ' + str(f"{merged_df[24][i]:.4e}") \
        + '  ' + str(f"{merged_df[25][i]:.4e}") + '  ' + str(f"{merged_df[26][i]:.4e}") + '\n'
    
    data_file = open(str(data_file_name), 'w') #write your file name in the open function
    n = data_file.write(txt_string)
    data_file.close()

    return merged_df