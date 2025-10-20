import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd


def read_extracted_file(e_file_name):
    """
    Reads, parses, and processes the pars files from Robitaille's SED fitter, returning
    a structured pandas DataFrame of stellar and circumstellar parameters. 
    This function only works for 01, 02, 16, and 17 model types.

    This function performs the following steps:
    1. Determines the model type from the file header.
    2. Extracts star names, metadata (number of data points, number of fits), 
       and parameter values.
    3. Organizes model parameters into arrays split by star.
    4. Applies physical cutoffs:
       - Removes models with stellar temperatures below the birthline threshold (10^3.6 K).
       - Removes models that fall below the luminosity–temperature relation from 
         Haemmerlé et al. (2019).
    5. Removes stars for which all models are eliminated by cutoffs.
    6. Returns results in a pandas DataFrame.

    Args:
        e_file_name (str): Path to the extracted model results file.

    Returns:
        pd.DataFrame: DataFrame containing merged and filtered model results, with columns:
            - 'model_type': Model type identifier.
            - 'star_names': Star identifiers.
            - 'n_datas': Number of data points per star.
            - 'num_fits': Number of model fits per star.
            - 'split_star_temp_arr': Stellar temperatures.
            - 'split_chi2_arr': χ² values.
            - 'split_of_avs': Extinction (A_V) values.
            - 'split_scale_arr': Scaling factors.
            - 'split_star_rad_arr': Stellar radii (in solar units).
            - 'split_d_mass_arr': Disk masses.
            - 'split_rmax_arr': Maximum disk radii.
            - 'split_d_beta_arr': Disk flaring exponents.
            - 'split_disk_p_arr': Disk density power-law indices.
            - 'split_disk_h100_arr': Disk scale heights at 100 AU.
            - 'split_scattering_arr': Scattering fractions.
            - 'split_inc_arr': Inclinations.
            - Additional columns depending on `model_type`:
                * "sp_h_i_02": 'split_disk_rmin_arr'
                * "spubsmi_16": envelope and cavity parameter arrays
                * "spubhmi_17": envelope, cavity, disk, and ambient parameter arrays
            - 'lum_array': Stellar luminosities computed from radius and temperature.

    Notes:
        - Stars with no valid models after filtering are removed entirely.
        - Stellar luminosities are calculated using the Stefan–Boltzmann law, 
          normalized to solar temperature (5772 K).
        - Luminosity cutoffs follow the empirical relation from 
          Haemmerlé et al. (2019).
    """
    
    #First we open to file
    init_file = open(e_file_name, 'r')
    
    #Following lines determine the model type in  the file
    test_list = init_file.readlines()[1];
    
    substr_1 = 'envelope'
    substr_2 = 'disk.rmin'
    
    if substr_1 in test_list and substr_2 in test_list:
        model_type = 'spubhmi_17'
    elif substr_1 in test_list and substr_2 not in test_list:
        model_type = 'spubsmi_16'
    elif substr_2 in test_list and substr_1 not in test_list:
        model_type = 'sp_h_i_02'
    else:
        model_type = 'sp_s_i_01'
    
    init_file.seek(0) #reset file reading for the rest to work properly
    
    #The line_list skips the first three lines in the file, as they are headers
    line_list = init_file.readlines()[3:];

    #OPTIONAL. Print the first 10 lines to see what line_list is (is a list of strings)
    #print(line_list[0:10])
    
    #Here we construct lists for the pandas dataframe later
    #The other thing is check_line, which records the line nuber at which the star name is located at
    #Counter is to check what line number we are on (a more rudimentary way of knowing)
    star_names = []
    n_datas = []
    num_fits = []
    count = 0
    check_line = []

    for line in line_list:
        count += 1
        if line.startswith('             '):
            #If Statement specifies that if the beginning of the line has 13 spaces, then it contains a star name
            check_line.append(count) 
            star_names.append(line[13:30])
            n_datas.append(float(line[39:41]))
            num_fits.append(float(line[48:52]))

    #OPTIONAL. Print the line numbers at which the star name occurs at and the star name list
    #printing the length is to find the number of stars in this fit
    #print(check_line)
    #print(len(check_line)) #This checks the number of stars in the list, theoretically
    
    #We create different arrays to contain the parameters. They are the length equal to the number of lines 
    #in the txt file we read in (minus the headers). 
    #This also erroneously includes the star_name lines. (solved later in code)

    chi2_arr = np.zeros(len(line_list)) #chi squared array
    av_arr = np.zeros(len(line_list)) #av array
    scale_arr = np.zeros(len(line_list)) #scale array
    star_rad_arr = np.zeros(len(line_list)) #star radius array
    star_temp_arr = np.zeros(len(line_list)) #star temperature array
    d_mass_arr = np.zeros(len(line_list)) #dust mass array
    d_rmax_arr = np.zeros(len(line_list))
    d_beta_arr = np.zeros(len(line_list))
    disk_p_arr = np.zeros(len(line_list))
    disk_h100_arr = np.zeros(len(line_list))
    
    scattering_arr = np.zeros(len(line_list))
    inc_arr = np.zeros(len(line_list)) #inclination array

    if model_type == 'sp_h_i_02':
        disk_rmin_arr = np.zeros(len(line_list))
    elif model_type == 'spubsmi_16':
        env_rho_0_arr = np.zeros(len(line_list)) #envelope density
        env_rc_arr = np.zeros(len(line_list))
        cav_power_arr = np.zeros(len(line_list))
        cav_theta_0_arr = np.zeros(len(line_list))
        cav_rho_0_arr = np.zeros(len(line_list))
        amb_den_arr = np.zeros(len(line_list))
        amb_temp_arr = np.zeros(len(line_list))
    elif model_type == 'spubhmi_17':
        env_rho_0_arr = np.zeros(len(line_list)) #envelope density
        env_rc_arr = np.zeros(len(line_list))
        cav_power_arr = np.zeros(len(line_list))
        cav_theta_0_arr = np.zeros(len(line_list))
        cav_rho_0_arr = np.zeros(len(line_list))
        disk_rmin_arr = np.zeros(len(line_list))
        env_rmin_arr = np.zeros(len(line_list))
        amb_den_arr = np.zeros(len(line_list))
        amb_temp_arr = np.zeros(len(line_list))
    
    #New counter
    counter_2 = 0

    #For loops runs through the line_list, grabbing out variable in its location. We note if it encounters a line
    #with the star name, nothing happens (actually a zero is added in the array as a consequence of nothing happening)

    ##Interesting comment for myself: Notice how counter_2 is put before any of the actual operations. That means
    #before anything happens, counter starts off at 1, and then passes to the next line. Since I'm not doing anything
    #at the passes, a zero is added to the array of arrays, acting as a marker or boundary between the stars.
    #However, more conventionally we put the counter after the operations, however here we do not. Notice 0 is not
    #in check_line, so its in the else clause. But this is bad since line 0 has nothing so an error occurs.

    for line in line_list:
        counter_2 += 1
        if counter_2 in check_line:
            pass
        else:
            chi_2 = float(line[45:52]) #potential to crossing into 100's so extra space in beginning is given
            av = float(line[56:63])
            scale = float(line[69:74])
            star_rad = float(line[76:85])
            star_temp = float(line[87:96])
            d_mass = float(line[98:108])
            d_rmax = float(line[109:118])
            d_beta = float(line[120:129])
            disk_p = float(line[130:140])
            disk_h100 = float(line[141:151])
            
            if model_type == 'sp_s_i_01':
                scattering = float(line[152:162])
                inc = float(line[163:173])
            elif model_type == 'sp_h_i_02':
                disk_rmin = float(line[152:162])
                scattering = float(line[163:173])
                inc = float(line[174:184])
            elif model_type == 'spubsmi_16':
                env_rho_0 = float(line[153:163])
                env_rc = float(line[164:174])
                cav_power = float(line[175:184])
                cav_theta_0 = float(line[185:195])
                cav_rho_0 = float(line[196:206])
                amb_den = float(line[207:217])
                amb_temp = float(line[218:229])
                scattering = float(line[230:240])
                inc = float(line[241:251])
            elif model_type == 'spubhmi_17':
                env_rho_0 = float(line[153:163])
                env_rc = float(line[164:174])
                cav_power = float(line[175:184])
                cav_theta_0 = float(line[185:195])
                cav_rho_0 = float(line[196:206])
                disk_rmin = float(line[207:217])
                env_rmin = float(line[218:229])
                amb_den = float(line[230:240])
                amb_temp = float(line[241:251])
                scattering = float(line[252:261])
                inc = float(line[262:273])
                
            #next we gather these parameters in the arrays
            chi2_arr[counter_2 - 1] = chi_2
            av_arr[counter_2 - 1] = av
            scale_arr[counter_2 - 1] = scale
            star_rad_arr[counter_2 - 1] = star_rad
            star_temp_arr[counter_2 - 1] = star_temp
            d_mass_arr[counter_2 - 1] = d_mass
            d_rmax_arr[counter_2 - 1] = d_rmax
            d_beta_arr[counter_2 - 1] = d_beta
            disk_p_arr[counter_2 - 1] = disk_p
            disk_h100_arr[counter_2 - 1] = disk_h100
            scattering_arr[counter_2 - 1] = scattering
            inc_arr[counter_2 - 1] = inc
            
            if model_type == 'sp_h_i_02':
                disk_rmin_arr[counter_2 - 1] = disk_rmin
            elif model_type == 'spubsmi_16':
                env_rho_0_arr[counter_2 - 1] = env_rho_0
                env_rc_arr[counter_2 - 1] = env_rc
                cav_power_arr[counter_2 - 1] = cav_power
                cav_theta_0_arr[counter_2 - 1] = cav_theta_0
                cav_rho_0_arr[counter_2 - 1] = cav_rho_0
                amb_den_arr[counter_2 - 1] = amb_den
                amb_temp_arr[counter_2 - 1] = amb_temp
            elif model_type == 'spubhmi_17':
                env_rho_0_arr[counter_2 - 1] = env_rho_0
                env_rc_arr[counter_2 - 1] = env_rc
                cav_power_arr[counter_2 - 1] = cav_power
                cav_theta_0_arr[counter_2 - 1] = cav_theta_0
                cav_rho_0_arr[counter_2 - 1] = cav_rho_0
                disk_rmin_arr[counter_2 - 1] = disk_rmin
                env_rmin_arr[counter_2 - 1] = env_rmin
                amb_den_arr[counter_2 - 1] = amb_den
                amb_temp_arr[counter_2 - 1] = amb_temp
            

    #We need to split the giant array into separate subarrays for each star. We do this by saying at the intervals
    #produced by check_line, we can split them.

    split_chi2_arr = np.split(chi2_arr, check_line)
    split_of_avs = np.split(av_arr, check_line)
    split_scale_arr = np.split(scale_arr, check_line)
    split_star_rad_arr = np.split(star_rad_arr, check_line)
    split_star_temp_arr = np.split(star_temp_arr, check_line)
    split_d_mass_arr = np.split(d_mass_arr, check_line)
    split_rmax_arr = np.split(d_rmax_arr, check_line)
    split_d_beta_arr = np.split(d_beta_arr, check_line)
    split_disk_p_arr = np.split(disk_p_arr, check_line)
    split_disk_h100_arr = np.split(disk_h100_arr, check_line)
    split_scattering_arr = np.split(scattering_arr, check_line)
    split_inc_arr = np.split(inc_arr, check_line)
    
    if model_type == 'sp_h_i_02':
        split_disk_rmin_arr = np.split(disk_rmin_arr, check_line)
    elif model_type == 'spubsmi_16':
        split_env_rho_0_arr = np.split(env_rho_0_arr, check_line)
        split_env_rc_arr = np.split(env_rc_arr, check_line)
        split_cav_power_arr = np.split(cav_power_arr, check_line)
        split_cav_theta_0_arr = np.split(cav_theta_0_arr, check_line)
        split_cav_rho_0_arr = np.split(cav_rho_0_arr, check_line)
        split_amb_den_arr = np.split(amb_den_arr, check_line)
        split_amb_temp_arr = np.split(amb_temp_arr, check_line)
    elif model_type == 'spubhmi_17':
        split_env_rho_0_arr = np.split(env_rho_0_arr, check_line)
        split_env_rc_arr = np.split(env_rc_arr, check_line)
        split_cav_power_arr = np.split(cav_power_arr, check_line)
        split_cav_theta_0_arr = np.split(cav_theta_0_arr, check_line)
        split_cav_rho_0_arr = np.split(cav_rho_0_arr, check_line)
        split_disk_rmin_arr = np.split(disk_rmin_arr, check_line)
        split_env_rmin_arr = np.split(env_rmin_arr, check_line)
        split_amb_den_arr = np.split(amb_den_arr, check_line)
        split_amb_temp_arr = np.split(amb_temp_arr, check_line)
    
    #Following line gets rid of that first 0 in the array
    split_chi2_arr = split_chi2_arr[1:]
    split_of_avs = split_of_avs[1:]
    split_scale_arr = split_scale_arr[1:] 
    split_star_rad_arr = split_star_rad_arr[1:]
    split_star_temp_arr = split_star_temp_arr[1:]
    split_d_mass_arr = split_d_mass_arr[1:]
    split_rmax_arr = split_rmax_arr[1:]
    split_d_beta_arr = split_d_beta_arr[1:]
    split_disk_p_arr = split_disk_p_arr[1:]
    split_disk_h100_arr = split_disk_h100_arr[1:]
    split_scattering_arr = split_scattering_arr[1:]
    split_inc_arr = split_inc_arr[1:]
    
    if model_type == 'sp_h_i_02':
        split_disk_rmin_arr = split_disk_rmin_arr[1:]
    elif model_type == 'spubsmi_16':
        split_env_rho_0_arr = split_env_rho_0_arr[1:]
        split_env_rc_arr = split_env_rc_arr[1:]
        split_cav_power_arr = split_cav_power_arr[1:]
        split_cav_theta_0_arr = split_cav_theta_0_arr[1:]
        split_cav_rho_0_arr = split_cav_rho_0_arr[1:]
        split_amb_den_arr = split_amb_den_arr[1:]
        split_amb_temp_arr = split_amb_temp_arr[1:]
    elif model_type == 'spubhmi_17':
        split_env_rho_0_arr = split_env_rho_0_arr[1:]
        split_env_rc_arr = split_env_rc_arr[1:]
        split_cav_power_arr = split_cav_power_arr[1:]
        split_cav_theta_0_arr = split_cav_theta_0_arr[1:]
        split_cav_rho_0_arr = split_cav_rho_0_arr[1:]
        split_disk_rmin_arr = split_disk_rmin_arr[1:]
        split_env_rmin_arr = split_env_rmin_arr[1:]
        split_amb_den_arr = split_amb_den_arr[1:]
        split_amb_temp_arr = split_amb_temp_arr[1:]
    
    #Following loop gets rid of that zero at the end of every array that is created by the counter
    #Note the range excludes last line since it isn't affected by it

    for i in range(0, len(check_line) - 1):
        split_chi2_arr[i] = split_chi2_arr[i][:-1]
        split_of_avs[i] = split_of_avs[i][:-1]
        split_scale_arr[i] = split_scale_arr[i][:-1]
        split_star_rad_arr[i] = split_star_rad_arr[i][:-1]
        split_star_temp_arr[i] = split_star_temp_arr[i][:-1]
        split_d_mass_arr[i] = split_d_mass_arr[i][:-1]
        split_rmax_arr[i] = split_rmax_arr[i][:-1]
        split_d_beta_arr[i] = split_d_beta_arr[i][:-1]
        split_disk_p_arr[i] = split_disk_p_arr[i][:-1]
        split_disk_h100_arr[i] = split_disk_h100_arr[i][:-1]
        split_scattering_arr[i] = split_scattering_arr[i][:-1]
        split_inc_arr[i] = split_inc_arr[i][:-1]
        
        if model_type == 'sp_h_i_02':
            split_disk_rmin_arr[i] = split_disk_rmin_arr[i][:-1]
        elif model_type == 'spubsmi_16':
            split_env_rho_0_arr[i] = split_env_rho_0_arr[i][:-1]
            split_env_rc_arr[i] = split_env_rc_arr[i][:-1]
            split_cav_power_arr[i] = split_cav_power_arr[i][:-1]
            split_cav_theta_0_arr[i] = split_cav_theta_0_arr[i][:-1]
            split_cav_rho_0_arr[i] = split_cav_rho_0_arr[i][:-1]
            split_amb_den_arr[i] = split_amb_den_arr[i][:-1]
            split_amb_temp_arr[i] = split_amb_temp_arr[i][:-1]
        elif model_type == 'spubhmi_17':
            split_env_rho_0_arr[i] = split_env_rho_0_arr[i][:-1]
            split_env_rc_arr[i] = split_env_rc_arr[i][:-1]
            split_cav_power_arr[i] = split_cav_power_arr[i][:-1]
            split_cav_theta_0_arr[i] = split_cav_theta_0_arr[i][:-1]
            split_cav_rho_0_arr[i] = split_cav_rho_0_arr[i][:-1]
            split_disk_rmin_arr[i] = split_disk_rmin_arr[i][:-1]
            split_env_rmin_arr[i] = split_env_rmin_arr[i][:-1]
            split_amb_den_arr[i] = split_amb_den_arr[i][:-1]
            split_amb_temp_arr[i] = split_amb_temp_arr[i][:-1]
        
    #Turn the list of arrays in a multidimensional array
    split_chi2_arr = np.array(split_chi2_arr, dtype=object)
    split_of_avs = np.array(split_of_avs, dtype=object)
    split_scale_arr = np.array(split_scale_arr, dtype=object)
    split_star_rad_arr = np.array(split_star_rad_arr, dtype=object)
    split_star_temp_arr = np.array(split_star_temp_arr, dtype=object)
    split_d_mass_arr = np.array(split_d_mass_arr, dtype=object)
    split_rmax_arr = np.array(split_rmax_arr, dtype=object)
    split_d_beta_arr = np.array(split_d_beta_arr, dtype=object)
    split_disk_p_arr = np.array(split_disk_p_arr, dtype=object)
    split_disk_h100_arr = np.array(split_disk_h100_arr, dtype=object)
    
    if model_type == 'sp_h_i_02':
        split_disk_rmin_arr = np.array(split_disk_rmin_arr, dtype=object)
    elif model_type == 'spubsmi_16':
        split_env_rho_0_arr = np.array(split_env_rho_0_arr, dtype=object)
        split_env_rc_arr = np.array(split_env_rc_arr, dtype=object)
        split_cav_power_arr = np.array(split_cav_power_arr, dtype=object)
        split_cav_theta_0_arr = np.array(split_cav_theta_0_arr, dtype=object)
        split_cav_rho_0_arr = np.array(split_cav_rho_0_arr, dtype=object)
        split_amb_den_arr = np.array(split_amb_den_arr, dtype=object)
        split_amb_temp_arr = np.array(split_amb_temp_arr, dtype=object)
    elif model_type == 'spubhmi_17':
        split_env_rho_0_arr = np.array(split_env_rho_0_arr, dtype=object)
        split_env_rc_arr = np.array(split_env_rc_arr, dtype=object)
        split_cav_power_arr = np.array(split_cav_power_arr, dtype=object)
        split_cav_theta_0_arr = np.array(split_cav_theta_0_arr, dtype=object)
        split_cav_rho_0_arr = np.array(split_cav_rho_0_arr, dtype=object)
        split_disk_rmin_arr = np.array(split_disk_rmin_arr, dtype=object)
        split_env_rmin_arr = np.array(split_env_rmin_arr, dtype=object)
        split_amb_den_arr = np.array(split_amb_den_arr, dtype=object)
        split_amb_temp_arr = np.array(split_amb_temp_arr, dtype=object)
    
    split_scattering_arr = np.array(split_scattering_arr, dtype=object)
    split_inc_arr = np.array(split_inc_arr, dtype=object)
        
    #Turn lists into arrays for numpy.delete not to throw an error
    check_line = np.array(check_line)
    star_names = np.array(star_names)
    n_datas = np.array(n_datas)
    num_fits = np.array(num_fits)
    
    #########################
    
    #Here we are removing all the temperatures along with their respective models that call below the cutoff
    for j in range(0, len(check_line)):
        index = np.where(split_star_temp_arr[j] < 10**3.6)
        split_star_temp_arr[j] = np.delete(split_star_temp_arr[j], index)
        split_chi2_arr[j] = np.delete(split_chi2_arr[j], index)
        split_of_avs[j] = np.delete(split_of_avs[j], index)
        split_scale_arr[j] = np.delete(split_scale_arr[j], index)
        split_star_rad_arr[j] = np.delete(split_star_rad_arr[j], index)
        split_d_mass_arr[j] = np.delete(split_d_mass_arr[j], index)   
        split_rmax_arr[j] = np.delete(split_rmax_arr[j], index)
        split_d_beta_arr[j] = np.delete(split_d_beta_arr[j], index)
        split_disk_p_arr[j] = np.delete(split_disk_p_arr[j], index)
        split_disk_h100_arr[j] = np.delete(split_disk_h100_arr[j], index)
        split_scattering_arr[j] = np.delete(split_scattering_arr[j], index)
        split_inc_arr[j] = np.delete(split_inc_arr[j], index)
        
        if model_type == 'sp_h_i_02':
            split_disk_rmin_arr[j] = np.delete(split_disk_rmin_arr[j], index)
        elif model_type == 'spubsmi_16':
            split_env_rho_0_arr[j] = np.delete(split_env_rho_0_arr[j], index)
            split_env_rc_arr[j] = np.delete(split_env_rc_arr[j], index)
            split_cav_power_arr[j] = np.delete(split_cav_power_arr[j], index)
            split_cav_theta_0_arr[j] = np.delete(split_cav_theta_0_arr[j], index)
            split_cav_rho_0_arr[j] = np.delete(split_cav_rho_0_arr[j], index)
            split_amb_den_arr[j] = np.delete(split_amb_den_arr[j], index)
            split_amb_temp_arr[j] = np.delete(split_amb_temp_arr[j], index)
        elif model_type == 'spubhmi_17':
            split_env_rho_0_arr[j] = np.delete(split_env_rho_0_arr[j], index)
            split_env_rc_arr[j] = np.delete(split_env_rc_arr[j], index)
            split_cav_power_arr[j] = np.delete(split_cav_power_arr[j], index)
            split_cav_theta_0_arr[j] = np.delete(split_cav_theta_0_arr[j], index)
            split_cav_rho_0_arr[j] = np.delete(split_cav_rho_0_arr[j], index)
            split_disk_rmin_arr[j] = np.delete(split_disk_rmin_arr[j], index)
            split_env_rmin_arr[j] = np.delete(split_env_rmin_arr[j], index)
            split_amb_den_arr[j] = np.delete(split_amb_den_arr[j], index)
            split_amb_temp_arr[j] = np.delete(split_amb_temp_arr[j], index)
            
    #First we get the index of empty array
    #IF no empty arrays, you can skip the next two blocks
    index_of_emp = []
    count_3 = -1
    for k in range(0,len(check_line)):
        count_3 += 1
        if split_star_temp_arr[k].size == 0:
            index_of_emp.append(count_3) 
            
    if index_of_emp != []:
        print("The following star(s) have been removed due to all models having temperatures falling below the birthline temperature:")
        for index_r in range(0, len(index_of_emp)):
            print("Star Name: " + star_names[index_of_emp[index_r]])
            
        index_of_emp = np.array(index_of_emp) #turning the index into an array
        #Now we delete the star since it has no models
        #First we delete the star attributes that were lists first
        check_line = np.delete(check_line, index_of_emp)
        star_names = np.delete(star_names, index_of_emp)
        n_datas = np.delete(n_datas, index_of_emp)
        num_fits = np.delete(num_fits, index_of_emp)

        #Now we delete the stellar parameter empty lists
        split_star_temp_arr = np.delete(split_star_temp_arr, index_of_emp)
        split_chi2_arr = np.delete(split_chi2_arr, index_of_emp)
        split_of_avs = np.delete(split_of_avs, index_of_emp)
        split_scale_arr = np.delete(split_scale_arr, index_of_emp)
        split_star_rad_arr = np.delete(split_star_rad_arr, index_of_emp)
        split_d_mass_arr = np.delete(split_d_mass_arr, index_of_emp)
        split_rmax_arr = np.delete(split_rmax_arr, index_of_emp)
        split_d_beta_arr = np.delete(split_d_beta_arr, index_of_emp)
        split_disk_p_arr = np.delete(split_disk_p_arr, index_of_emp)
        split_disk_h100_arr = np.delete(split_disk_h100_arr, index_of_emp)
        
        split_scattering_arr = np.delete(split_scattering_arr, index_of_emp)
        split_inc_arr = np.delete(split_inc_arr, index_of_emp)
        
        if model_type == 'sp_h_i_02':
            split_disk_rmin_arr = np.delete(split_disk_rmin_arr, index_of_emp)
        elif model_type == 'spubsmi_16':
            split_env_rho_0_arr = np.delete(split_env_rho_0_arr, index_of_emp)
            split_env_rc_arr = np.delete(split_env_rc_arr, index_of_emp)
            split_cav_power_arr = np.delete(split_cav_power_arr, index_of_emp)
            split_cav_theta_0_arr = np.delete(split_cav_theta_0_arr, index_of_emp)
            split_cav_rho_0_arr = np.delete(split_cav_rho_0_arr, index_of_emp)
            split_amb_den_arr = np.delete(split_amb_den_arr, index_of_emp)
            split_amb_temp_arr = np.delete(split_amb_temp_arr, index_of_emp)
        elif model_type == 'spubhmi_17':
            split_env_rho_0_arr = np.delete(split_env_rho_0_arr, index_of_emp)
            split_env_rc_arr = np.delete(split_env_rc_arr, index_of_emp)
            split_cav_power_arr = np.delete(split_cav_power_arr, index_of_emp)
            split_cav_theta_0_arr = np.delete(split_cav_theta_0_arr, index_of_emp)
            split_cav_rho_0_arr = np.delete(split_cav_rho_0_arr, index_of_emp)
            split_disk_rmin_arr = np.delete(split_disk_rmin_arr, index_of_emp)
            split_env_rmin_arr = np.delete(split_env_rmin_arr, index_of_emp)
            split_amb_den_arr = np.delete(split_amb_den_arr, index_of_emp)
            split_amb_temp_arr = np.delete(split_amb_temp_arr, index_of_emp)

    #Luminosity calculations
    #First initialize luminosity array. We need a luminosity for each individual model.
    lum_array = split_scale_arr.copy()

    #For loop to add in lum. values. We usee the Sun's temp. as 5772 K. Radii of stars are in solar units already
    for j in range(0, len(check_line)):
        lum_array[j] = (split_star_rad_arr[j])**2 * (split_star_temp_arr[j] / 5772)**4

    #OPTIONAL. Print lum_array. Should be another array of subarrays and the same length as check_line
    #print(len(lum_array))
    #print(lum_array)
    
    #Now we define an inequality for the luminosity
    #The line comes from the isochrone birthline in Haemmerlé et al. (2019).

    #Eqn. for cutoff: bx^k = y
    x = np.linspace(10, 10**2.5, 1000)
    y = 5.5 * np.log(x) - 20.7

    for k in range(0, len(check_line)):
        cutoff = np.log10(lum_array[k]) < 5.5 * np.log10(split_star_temp_arr[k]) - 20.7

        cut = np.where(cutoff == True)
        split_star_temp_arr[k] = np.delete(split_star_temp_arr[k], cut)
        lum_array[k] = np.delete(lum_array[k], cut)
        split_chi2_arr[k] = np.delete(split_chi2_arr[k], cut)
        split_of_avs[k] = np.delete(split_of_avs[k], cut)
        split_scale_arr[k] = np.delete(split_scale_arr[k], cut)
        split_star_rad_arr[k] = np.delete(split_star_rad_arr[k], cut)
        split_d_mass_arr[k] = np.delete(split_d_mass_arr[k], cut)
        split_rmax_arr[k] = np.delete(split_rmax_arr[k], cut)
        split_d_beta_arr[k] = np.delete(split_d_beta_arr[k], cut)
        split_disk_p_arr[k] = np.delete(split_disk_p_arr[k], cut)
        split_disk_h100_arr[k] = np.delete(split_disk_h100_arr[k], cut)
        split_scattering_arr[k] = np.delete(split_scattering_arr[k], cut)
        split_inc_arr[k] = np.delete(split_inc_arr[k], cut)
        
        if model_type == 'sp_h_i_02':
            split_disk_rmin_arr[k] = np.delete(split_disk_rmin_arr[k], cut)
        elif model_type == 'spubsmi_16':
            split_env_rho_0_arr[k] = np.delete(split_env_rho_0_arr[k], cut)
            split_env_rc_arr[k] = np.delete(split_env_rc_arr[k], cut)
            split_cav_power_arr[k] = np.delete(split_cav_power_arr[k], cut)
            split_cav_theta_0_arr[k] = np.delete(split_cav_theta_0_arr[k], cut)
            split_cav_rho_0_arr[k] = np.delete(split_cav_rho_0_arr[k], cut)
            split_amb_den_arr[k] = np.delete(split_amb_den_arr[k], cut)
            split_amb_temp_arr[k] = np.delete(split_amb_temp_arr[k], cut)
        elif model_type == 'spubhmi_17':
            split_env_rho_0_arr[k] = np.delete(split_env_rho_0_arr[k], cut)
            split_env_rc_arr[k] = np.delete(split_env_rc_arr[k], cut)
            split_cav_power_arr[k] = np.delete(split_cav_power_arr[k], cut)
            split_cav_theta_0_arr[k] = np.delete(split_cav_theta_0_arr[k], cut)
            split_cav_rho_0_arr[k] = np.delete(split_cav_rho_0_arr[k], cut)
            split_disk_rmin_arr[k] = np.delete(split_disk_rmin_arr[k], cut)
            split_env_rmin_arr[k] = np.delete(split_env_rmin_arr[k], cut)
            split_amb_den_arr[k] = np.delete(split_amb_den_arr[k], cut)
            split_amb_temp_arr[k] = np.delete(split_amb_temp_arr[k], cut)    
            
    index_of_emp_2 = []
    count_4 = -1
    for k in range(0,len(check_line)):
        count_4 += 1
        if split_star_temp_arr[k].size == 0:
            index_of_emp_2.append(count_4)

    if index_of_emp_2 != []:
        print("The following star(s) have been removed due to having no models with realistic luminosities:")
        for index_r_2 in range(0, len(index_of_emp_2)):
            print("Star Name: " + star_names[index_of_emp_2[index_r_2]])
        index_of_emp_2 = np.array(index_of_emp_2) #turning the index into an array

        #Now we delete the star since it has no models
        #First we delete the star attributes that were lists first
        check_line = np.delete(check_line, index_of_emp_2)
        star_names = np.delete(star_names, index_of_emp_2)
        n_datas = np.delete(n_datas, index_of_emp_2)
        num_fits = np.delete(num_fits, index_of_emp_2)

        #Now we delete the stellar parameter empty lists
        split_star_temp_arr = np.delete(split_star_temp_arr, index_of_emp_2)
        lum_array = np.delete(lum_array, index_of_emp_2)
        split_chi2_arr = np.delete(split_chi2_arr, index_of_emp_2)
        split_of_avs = np.delete(split_of_avs, index_of_emp_2)
        split_scale_arr = np.delete(split_scale_arr, index_of_emp_2)
        split_star_rad_arr = np.delete(split_star_rad_arr, index_of_emp_2)
        split_d_mass_arr = np.delete(split_d_mass_arr, index_of_emp_2)
        split_rmax_arr = np.delete(split_rmax_arr, index_of_emp_2)
        split_d_beta_arr = np.delete(split_d_beta_arr, index_of_emp_2)
        split_disk_p_arr = np.delete(split_disk_p_arr, index_of_emp_2)
        split_disk_h100_arr = np.delete(split_disk_h100_arr, index_of_emp_2)
        split_scattering_arr = np.delete(split_scattering_arr, index_of_emp_2)
        split_inc_arr = np.delete(split_inc_arr, index_of_emp_2)
        
        if model_type == 'sp_h_i_02':
            split_disk_rmin_arr = np.delete(split_disk_rmin_arr, index_of_emp_2)
        elif model_type == 'spubsmi_16':
            split_env_rho_0_arr = np.delete(split_env_rho_0_arr, index_of_emp_2)
            split_env_rc_arr = np.delete(split_env_rc_arr, index_of_emp_2)
            split_cav_power_arr = np.delete(split_cav_power_arr, index_of_emp_2)
            split_cav_theta_0_arr = np.delete(split_cav_theta_0_arr, index_of_emp_2)
            split_cav_rho_0_arr = np.delete(split_cav_rho_0_arr, index_of_emp_2)
            split_amb_den_arr = np.delete(split_amb_den_arr, index_of_emp_2)
            split_amb_temp_arr = np.delete(split_amb_temp_arr, index_of_emp_2)
        elif model_type == 'spubhmi_17':
            split_env_rho_0_arr = np.delete(split_env_rho_0_arr, index_of_emp_2)
            split_env_rc_arr = np.delete(split_env_rc_arr, index_of_emp_2)
            split_cav_power_arr = np.delete(split_cav_power_arr, index_of_emp_2)
            split_cav_theta_0_arr = np.delete(split_cav_theta_0_arr, index_of_emp_2)
            split_cav_rho_0_arr = np.delete(split_cav_rho_0_arr, index_of_emp_2)
            split_disk_rmin_arr = np.delete(split_disk_rmin_arr, index_of_emp_2)
            split_env_rmin_arr = np.delete(split_env_rmin_arr, index_of_emp_2)
            split_amb_den_arr = np.delete(split_amb_den_arr, index_of_emp_2)
            split_amb_temp_arr = np.delete(split_amb_temp_arr, index_of_emp_2)

    #Now to accurately update the number of good fits a star has (num_fits)
    #also keep the initial number of fits
    init_num_fits = num_fits.copy()
    for i_star in range(0,len(check_line)):
        num_fits[i_star] = len(split_chi2_arr[i_star])
        
    #Now we make the pandas data frame
    if model_type == 'sp_s_i_01':
        d = {"MIR_NAME": star_names, "n_data": n_datas, "n_fits": num_fits, "init_n_fits": init_num_fits, "chi_2_arr": split_chi2_arr,\
             "av_arr": split_of_avs, "scale_arr": split_scale_arr, "star_rad_arr": split_star_rad_arr,\
             "star_temp_arr": split_star_temp_arr, "d_mass_arr": split_d_mass_arr, "rmax_arr": split_rmax_arr,\
             "d_beta_arr": split_d_beta_arr, "disk_p_arr": split_disk_p_arr, "disk_h100_arr": split_disk_h100_arr,\
             "scattering_arr": split_scattering_arr, "inclination_arr": split_inc_arr, "lum_arr": lum_array}
    elif model_type == 'sp_h_i_02':
        d = {"MIR_NAME": star_names, "n_data": n_datas, "n_fits": num_fits, "init_n_fits": init_num_fits, "chi_2_arr": split_chi2_arr,\
             "av_arr": split_of_avs, "scale_arr": split_scale_arr, "star_rad_arr": split_star_rad_arr,\
             "star_temp_arr": split_star_temp_arr, "d_mass_arr": split_d_mass_arr, "rmax_arr": split_rmax_arr,\
             "d_beta_arr": split_d_beta_arr, "disk_p_arr": split_disk_p_arr, "disk_h100_arr": split_disk_h100_arr,\
             "disk_rmin_arr": split_disk_rmin_arr, "scattering_arr": split_scattering_arr,\
             "inclination_arr": split_inc_arr, "lum_arr": lum_array}
    elif model_type == 'spubsmi_16':
        d = {"MIR_NAME": star_names, "n_data": n_datas, "n_fits": num_fits, "init_n_fits": init_num_fits, "chi_2_arr": split_chi2_arr,\
             "av_arr": split_of_avs, "scale_arr": split_scale_arr, "star_rad_arr": split_star_rad_arr,\
             "star_temp_arr": split_star_temp_arr, "d_mass_arr": split_d_mass_arr, "rmax_arr": split_rmax_arr,\
             "d_beta_arr": split_d_beta_arr, "disk_p_arr": split_disk_p_arr, "disk_h100_arr": split_disk_h100_arr,\
             "env_rho_0_arr": split_env_rho_0_arr, "env_rc_arr": split_env_rc_arr, "cav_power_arr": split_cav_power_arr,\
             "cav_theta_0_arr": split_cav_theta_0_arr, "cav_rho_0_arr": split_cav_rho_0_arr,\
             "amb_den_arr": split_amb_den_arr, "amb_temp_arr": split_amb_temp_arr,\
             "scattering_arr": split_scattering_arr, "inclination_arr": split_inc_arr, "lum_arr": lum_array}
    elif model_type == 'spubhmi_17':
        d = {"MIR_NAME": star_names, "n_data": n_datas, "n_fits": num_fits, "init_n_fits": init_num_fits, "chi_2_arr": split_chi2_arr,\
             "av_arr": split_of_avs, "scale_arr": split_scale_arr, "star_rad_arr": split_star_rad_arr,\
             "star_temp_arr": split_star_temp_arr, "d_mass_arr": split_d_mass_arr, "rmax_arr": split_rmax_arr,\
             "d_beta_arr": split_d_beta_arr, "disk_p_arr": split_disk_p_arr, "disk_h100_arr": split_disk_h100_arr,\
             "env_rho_0_arr": split_env_rho_0_arr, "env_rc_arr": split_env_rc_arr, "cav_power_arr": split_cav_power_arr,\
             "cav_theta_0_arr": split_cav_theta_0_arr, "cav_rho_0_arr": split_cav_rho_0_arr, "disk_rmin_arr": split_disk_rmin_arr,\
             "env_rmin_arr": split_env_rmin_arr, "amb_den_arr": split_amb_den_arr, "amb_temp_arr": split_amb_temp_arr,\
             "scattering_arr": split_scattering_arr, "inclination_arr": split_inc_arr, "lum_arr": lum_array}

    df = pd.DataFrame(d)
    
    return df

def weight_mean(df, col_name):
    """
    Compute the simple and weighted mean of a given parameter for each star in the input DataFrame.
    Weights are based on chi-square values, with smaller chi-square values contributing more strongly.

    Args:
        df (pd.DataFrame): Input DataFrame from 'read_extracted_file' function. 
        DataFrame should contain, at minimum, the following columns:
            - 'MIR_NAME' (str): Identifier for each star.
            - 'n_data' (int): Number of data points used in the fit.
            - 'n_fits' (int): Number of models used in the fitting process.
            - 'chi_2_arr' (array-like of float): Array of chi-square values corresponding to model fits.
            - col_name (array-like of float): Column of parameter values to average (e.g., 'Teff_arr').
        col_name (str): Name of the column in `df` containing arrays of parameter values 
                        for which both simple and weighted means should be calculated.

    Returns:
        pd.DataFrame: DataFrame with one row per star and the following columns:
            - 'MIR_NAME' (str): Star identifier.
            - 'n_data' (int): Number of data points.
            - 'n_fits' (int): Number of model fits.
            - f"{col_name.replace('arr', 'avg')}" (float): Simple (unweighted) mean of the parameter values.
            - f"{col_name.replace('arr', 'w_avg')}" (float): Weighted mean of the parameter values, 
              using chi-square weighting.
    """
    #Initialize empty dataframe
    rows = []
    
    for i in range(0, len(df)):
        star_name = df['MIR_NAME'][i]
        #n_data = df['n_data'][i]
        #n_fit = df['n_fits'][i]
        chi_sq_arr = df['chi_2_arr'][i]
        
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
        new_row = {'MIR_NAME': df['MIR_NAME'][i], 'n_data': df['n_data'][i], 'n_fits': df['n_fits'][i], \
                               col_name.replace('arr', 'avg'): param_mean, col_name.replace('arr', 'w_avg'): param_weighted_mean}
        rows.append(new_row)

    result_df = pd.DataFrame(rows)
    return result_df

def multi_param_weight_avg(df, param_list):
    """
    Compute weighted and unweighted means for multiple parameters in a DataFrame.

    For each parameter in `param_list`, this function calls `weight_mean` to compute
    both the simple mean and the chi-squared–weighted mean across model fits for each star.
    The results are concatenated into a single DataFrame.

    Args:
        df (pandas.DataFrame):  
            Input DataFrame containing (at minimum) the following columns:
                - MIR_NAME (str): Unique star identifier.  
                - n_data (int): Number of data points used in the fit.  
                - n_fits (int): Number of models fitted.  
                - chi_2_arr (array-like): Array of chi-squared values for each model fit.  
                - <param> (array-like): Arrays of parameter values for each model fit
                  (e.g., "teff_arr", "av_arr"), where names correspond to `param_list`.  

        param_list (list of str):  
            List of column names (arrays of model-fit parameter values) for which
            weighted and unweighted averages will be calculated.

    Returns:
        pandas.DataFrame:  
            DataFrame containing one row per star, with the following columns:
                - MIR_NAME (str)  
                - n_data (int)  
                - n_fits (int)  
                - For each parameter in `param_list`:  
                    - `<param>_avg` (float): Unweighted mean value.  
                    - `<param>_w_avg` (float): Weighted mean value (using chi-squared weights).  
    """
    for k in range(0, len(param_list)):
        if k == 0:
            big_df = weight_mean(df, param_list[k])          
        else:
            new_df = weight_mean(df, param_list[k])
            avg_cols = new_df.iloc[:, 3:5]
            big_df = pd.concat([big_df, avg_cols],axis=1)
    return big_df