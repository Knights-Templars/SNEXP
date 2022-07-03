#!/usr/bin/env python
# coding: utf-8


## Light curves plotting and interpolating
## Anirban Dutta


# Import packages and module
#--------------------------------------------------------------------------------#
import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from sklearn.gaussian_process.kernels import Matern
from sklearn.gaussian_process.kernels import ConstantKernel
from sklearn.gaussian_process import GaussianProcessRegressor
#--------------------------------------------------------------------------------#

# A filter file containing all the required information on bandpasses
FILTER_directory = '/Users/anirbandutta/Dropbox/astromatic/'            # Filter information 
FILTER_data = FILTER_directory+'FILTERS.dat'

#--------------------------------------------------------------------------------#

# Read the filter file containing data on filters 

filter_df = pd.read_csv(FILTER_data, sep = '\s+')
filter_df = filter_df.replace('INDEF', np.nan).set_index(['FILTER', 'Marker', 'Color']).astype('float64')
filter_df = filter_df.reset_index().set_index('FILTER')
list_filters = filter_df.index.tolist()

# Cell that contains all the important functions of the code

class general:
    
    '''
    Remarks: A function to group similar kind of files in a directory.
    text_list: A text file used to store list of files
    common_text: A string (e.g. *.fits, *.list) used for grouping similar
    kinds of files
    exceptions: string of file name to exclude in grouping
    
    returns: list of grouped files
    '''
    
    def __init__(self, common_text):
        
        self.list_files = glob.glob(common_text)
        
    def group_similar_files(self, text_list, exceptions= ''):
        
        if exceptions != '':
            list_exceptions = exceptions.split(',')
            for text in list_exceptions:
                list_files = filter(lambda x: not re.search(text, x), self.list_files)
            
        self.list_files.sort()
        if len(text_list) != 0:
            with open(SN_directory + text_list, 'w') as f:
                for file_name in self.list_files:
                    f.write(file_name + '\n')
                    
        return self.list_files
        

class create_lc:
    
    '''
    The light curves should be in the form of eg. JD, U, UErr, B, BErr, V, VErr, R, RErr, I, IErr
    or otherwise JD, FILTER, MAG, MERR. 
    Other forms will raise error.
    '''
    
    def __init__(self, file_name):
        
        self.object_df = pd.read_csv(file_name, sep='\s+', engine='python', comment='#')
        
    def column_to_row(self, offset=0.0):
        
        '''
        Remarks: This function will reshape a dataframe from wide frame to long frame. 
        
        Output format: ['JD', 'FILTER', 'MAG', 'MERR']
        '''
        
        
        list_del = ['Date', 'Epoch', 'Seeing', 'EXPTIME', 'Telescope', 'Phase']
        
        if 'Epoch' in self.object_df.columns.values:
            self.object_df['JD'] = self.object_df['JD'] + self.object_df['Epoch']
        else:
            self.object_df['JD'] = self.object_df['JD'] + offset
        
        for column in list_del:
            if column in self.object_df.columns.values:
                self.object_df = self.object_df.drop(column, axis =1)
                
        list_filters = [x for x in self.object_df.columns.values if 'JD' not in x and 'Err' not in x]
        list_err = [x for x in self.object_df.columns.values if 'Err' in x]
        
        mag_df = pd.melt(self.object_df, id_vars=['JD'], value_vars=list_filters, 
                            var_name='FILTER', value_name='MAG')

        err_df = pd.melt(self.object_df, id_vars=['JD'], value_vars=list_err,
                        var_name='FILTER', value_name='MERR')
        
        output_df = mag_df.join(err_df['MERR'])
        
        return output_df
    
    
    def apparent_mag_df(self):
    
    
        '''
        Remarks: Creates an apparent magnitude dataframe. Assumes the light curves are in long frame.
        file_name: file containing the SN magnitudes
        (JD, Filter, FMAG, FERR)
        tempsub: Whether the magnitudes are template subtracted
        returns: apparent magnitude dataframe
        '''
    
    
        appmag_df = self.object_df.sort_values(by = ['FILTER', 'JD'], kind = 'mergesort')
        appmag_df = appmag_df[['JD','FILTER', 'FMAG', 'FERR']].reset_index(drop = True)
    
        return appmag_df

    def get_swift_lc(self):
        
        '''
        Get Swift light curve in a dataframe format
        '''
            
        self.object_df.rename(columns={'Filter': 'FILTER'}, inplace=True)
        total_df = self.object_df
        mag_df = self.object_df[['JD', 'FILTER', 'MAG', 'MERR']]
            
        return mag_df, total_df


    def jd_to_date(self, jd, input_fmt='jd', output_fmt='isot'):
        
        '''
        Remarks: Convert jd to date
        jd: Input jd
        input_fmt: format of input ('jd')
        output_fmt: format of output ('iso', 'isot')
        returns: date
        '''
    
        from astropy.time import Time
        times = Time(jd, format=input_fmt)
        date_time = times.to_value(output_fmt)
        date = date_time.split('T')[0]
    
        return date

    
    def restframe(self, t0, z, epoch, mode='row'):
        
        if mode == 'row':
            object_dataframe = self.apparent_mag_df()
            object_dataframe['Phase'] = object_dataframe['JD'].apply(lambda x: (x - t0) / (1 + z))
            object_dataframe['Date'] = object_dataframe['JD'].apply(lambda x: self.jd_to_date(x))
        
        elif mode == 'column':
            object_dataframe = self.column_to_row(offset=epoch)
            object_dataframe['Phase'] = object_dataframe['JD'].apply(lambda x: (x - t0) / (1 + z))
            object_dataframe['Date'] = object_dataframe['JD'].apply(lambda x: self.jd_to_date(x))
        
        return object_dataframe

    def ext_cor_mag(self, mag, band, ebv_Gal, ebv_host=0.0, host=False):
    
        if not host:
            ext_mag = mag - (filter_df.loc[band, 'RLambda'] * ebv_Gal)
        else:
            ext_mag = mag - (filter_df.loc[band, 'RLambda'] * ebv_Gal) - \
                        (filter_df.loc[band, 'RLambda'] * ebv_host)
        
        return ext_mag


    def ext_cor_err(self, merr, band, ebv_Gerr, ebv_herr=0.0, host=False):
    
        if not host:
            ext_err = np.sqrt(merr**2 + (filter_df.loc[band, 'RLambda'] * ebv_Gerr)**2)
        else:
            ext_err = np.sqrt(merr**2 + (filter_df.loc[band, 'RLambda'] * ebv_Gerr)**2 + \
                         (filter_df.loc[band, 'RLambda'] * ebv_herr)**2)
        
        return ext_err


    def color(self, df, color, ebv_Gal, ebv_Gerr, ebv_host=0.0, ebv_herr=0.0, ext_cor=True, host=False):
        
        colors = color.split('-')
        df = df.set_index('JD')
    
        dict_mag = {}
    
        for index, rows in df.iterrows():
            if index not in dict_mag.keys():
                dict_mag[index] = {}
            
            dict_mag[index][rows['FILTER']] = rows['MAG']
            dict_mag[index][rows['FILTER']+str('Err')] = rows['MERR']
    
        mag_col_df = pd.DataFrame(dict_mag).T
        mag_col_df.index.name = 'JD'
    
        mag_col_df = mag_col_df.reset_index()
    
        color_df = mag_col_df[['JD', colors[0], colors[0]+str('Err'), colors[1], colors[1]+str('Err')]]
        color_df = color_df.dropna()
    
        if not ext_cor:
            color_df[color] = color_df.apply(lambda x: x[colors[0]] -  x[colors[1]], axis=1)
        else:
        
            color_df[colors[0]+str('ex')] = color_df[colors[0]].apply(lambda x: self.ext_cor_mag(x, colors[0], ebv_Gal=ebv_Gal, ebv_host=ebv_host, host=host))
            color_df[colors[0]+str('ex')+'Err'] = color_df[colors[0]+'Err'].apply(lambda x: self.ext_cor_err(x, colors[0], ebv_Gerr=ebv_Gerr, ebv_herr=ebv_herr, host=host))
            color_df[colors[1]+str('ex')] = color_df[colors[1]].apply(lambda x: self.ext_cor_mag(x, colors[1], ebv_Gal=ebv_Gal, ebv_host=ebv_host, host=host))
            color_df[colors[1]+str('ex')+'Err'] = color_df[colors[1]+'Err'].apply(lambda x: self.ext_cor_err(x, colors[1], ebv_Gerr=ebv_Gerr, ebv_herr=ebv_herr, host=host))
        
            color_df[color] = color_df.apply(lambda x: x[colors[0]+str('ex')] - x[colors[1]+str('ex')], axis=1)
            color_df[color+'Err'] = color_df.apply(lambda x: np.sqrt(x[colors[0]+str('ex')+'Err']**2 + x[colors[1]+str('ex')+'Err']**2), axis=1)
        
        return color_df
    
    
    
    def plot_params(self, ax, ml_x, mil_x, ml_y, mil_y, invert=False):
    
        '''
        Remarks: Plotting parameters
        ax: axis object
        ml_x: major locator 'x'
        mil_x: minor locator 'x'
        ml_y: major locator 'y'
        mil_y: minor locator 'y'
        invert: Bool for the y-axis to be inverted 
    
        '''
    
        if invert:
            ax.invert_yaxis()
        
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.tick_params(axis = 'both', which = 'major', direction = 'in', length = 12, width = 2.0, labelsize = 12)
        ax.tick_params(axis = 'both', which = 'minor', direction = 'in', length = 6, width = 1.2, labelsize = 12)
        ax.xaxis.set_major_locator(MultipleLocator(ml_x))
        ax.xaxis.set_minor_locator(MultipleLocator(mil_x))
        ax.yaxis.set_major_locator(MultipleLocator(ml_y))
        ax.yaxis.set_minor_locator(MultipleLocator(mil_y))
    
    
    def plot_lc(self, df):
        
        '''
        Plot the light curves in each band in the rest frame.
        
        '''
        
        
        fig_app = plt.figure(figsize = (10, 10))
        ax = fig_app.add_subplot(111)
        self.plot_params(ax, ml_x = 20, mil_x = 5, ml_y = 2, mil_y = 0.5, invert=True)
        legend_properties = {'size':18, 'weight':'book'}
        plt.rcParams['font.family'] = 'DejaVu Sans'
        plt.rcParams['axes.linewidth'] = 3
        
        #U, B, V, R, I
        offset = [2.7, 1.2, 0.0, -1.0, -2.0, -3.0]
        str_offset = ['+2.7', '+1.2', ' ', '-1.0', '-2.0', '-3.0']
        i = 0
        line1 =[]
        for band, band_df in df.groupby('FILTER', sort=False):
            line1 += ax.plot(band_df['Phase'], band_df['MAG']+offset[i], mfc = filter_df.loc[band, 'Color'],
                             mec = filter_df.loc[band, 'Color'], markeredgewidth=3.5, marker = filter_df.loc[band, 'Marker'],  
                             markersize = 10, label = str(band) + str_offset[i], alpha = 0.7, ls = ' ')
                     
            ax.errorbar(band_df['Phase'], band_df['MAG']+offset[i], yerr = band_df['MERR'], fmt = '',
                        c = filter_df.loc[band, 'Color'], ls = '', lw = 0.7, capsize=2, capthick=1)
                    
    
            line1.append(line1)
    
            i = i+1
            
            
        ax.set_xlabel(r'Rest frame time since $B$-band maximum [days]', fontsize = 20)
        ax.set_ylabel('Apparent Magnitude [mag] + constant', fontsize = 20)
  

        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize = 20)
            tick.label1.set_fontweight('bold')
    
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize = 20)
            tick.label1.set_fontweight('bold')
    
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)   
        
        ax.legend(frameon= False, fancybox = False, shadow = False, framealpha=0.3, facecolor = 'white', 
                     prop = legend_properties, loc = 'upper right')
        
        
        plt.show()
    
    
class gp_interpolation:
    
    def __init__(self, df, band, supernova):
        
        self.object_df = df[df['FILTER']==band]
        self.x = np.ravel(self.object_df['JD'].values)           # x-value
        self.y = np.ravel(self.object_df['MAG'].values)          # y-value
        self.yerr = np.ravel(self.object_df['MERR'].values)      # error on y value
        self.band = band
        self.supernova = supernova
        
        
    def find_local_minima(self, arr, n):
    
        '''
        Remarks: Finds the minima
        arr: array to be used as input.
        n: length of the array.
    
        returns: index of minima
        '''
    
        idx_min = []
        
        for i in range(1, n-1):
            if (arr[i-1] > arr[i] < arr[i+1]):
                idx_min.append(i)
        if (arr[-1] < arr[-2]):
            idx_min.append(n-1)
        
        return idx_min
    
    
    def remove_outliers(self, list_, z=3):
    
        '''
        Remarks: Remove outliers from a list based on z-score
        list_: Input list
        z: z-score 
        returns: new_list removing outliers
        '''
        
        mean = np.mean(list_)
        std = np.std(list_)
        new_list = []
        [new_list.append(i) for i in list_ if abs(i-mean)/std < z]
    
        return new_list
    
    
    def iterate_outliers(self, list_, iterations=3):
    
        '''
        Remarks: Iterate over lists to remove outliers. Calls remove_outliers function
        list_: Input list
        iterations: No of iterations to performed, default=3
        returns: new_list removing outliers.
        '''
        
        new_list = list_[:]
        for i in range(iterations):
            new_list = self.remove_outliers(new_list)

        return new_list
    
    
    def interpolate_lc(self, amp=None, scale=None, diff_deg=None,
                       save_pred=False, find_two_peaks=False,
                       verbose=False, plot_interp=True, save_file=''):
    
        
        '''
        Remarks: Interpolate Light Curve with Gaussian Process
        Regression.
        df: SN LC dataframe
        band: Filter 
        amp: amplitude of kernel (typically set by the amount by which the function
        varies.default is set by the standard deviation of the magnitude.)
        scale: scale over which the function typically varies.(default is 10)
        diff_deg: degree of differentiability. Controls the smoothness of the function.
        save_pred (Bool): Bool to save the predictions as pickle file for later use.
        find_two_peaks (Bool): normal type Ia supernova has two peaks in the I-band.
        verbose (Bool): Print the results on the screen 
    
        returns: LC parameters after fit
        (jd, magnitude, magnitude_err) at maximum.
    
        '''
    
        mean = lambda x: x*0 + np.median(self.y)                 # construction of the mean function
            
        if amp is None:
            amp = np.std(self.y - mean(self.x))
        elif amp == 'LC':
            amp = np.std(self.y)
        if scale is None:
            scale = 10
        if diff_deg is None:
            diff_deg = 1
        
        # Define the kernel
        # A constant kernel plus the matern kernel
    
        kernel = ConstantKernel(amp, constant_value_bounds='fixed') *                  Matern(length_scale=scale, length_scale_bounds='fixed', nu=diff_deg+0.5)
    
        # Implement Gaussian Process Regressor
        gp = GaussianProcessRegressor(kernel=kernel, alpha=self.yerr**2)
        # Fit the gp model
        X = np.array([self.x]).T
        Y = self.y - mean(self.x)                                 # residual
        gauss_process = gp.fit(X, Y)
    
        ## Make predictions at points
        list_jd = np.atleast_2d(np.linspace(self.object_df['JD'].min()-8, self.object_df['JD'].max(), 1000)).T
    
        ## Make predictions
        x_pred = np.ravel(list_jd)
        res, sig_pred = gp.predict(list_jd, return_std=True)
        res += mean(x_pred)
    
        # Calculate the maxima of the LC and the jd for the maxima.
        n = len(res)
        idx_minima = self.find_local_minima(res, n)
        
        
    
        ## I-band will have a secondary maxima
        ## Need to modify this part of the code in an optimum way.
        if not find_two_peaks:
            
            jd_max = x_pred[idx_minima[0]]
            mag_max = res[idx_minima[0]]
            err_max = sig_pred[idx_minima[0]]
            jd_15 = jd_max + 15
            res_15, err_15 = gp.predict(np.atleast_2d(jd_15).T, return_std=True)
            mag_15 = res_15 + mean(jd_15)
            dm15 = mag_15 - mag_max
            dm15_err = np.sqrt(err_max**2 + err_15**2)
            
        else:
            jd_max = x_pred[idx_minima[0]]
            mag_max = res[idx_minima[0]]
            err_max = sig_pred[idx_minima[0]]
            jd_15 = jd_max + 15
            res_15, err_15 = gp.predict(np.atleast_2d(jd_15).T, return_std=True)
            mag_15 = res_15 + mean(jd_15)
            dm15 = mag_15 - mag_max
            dm15_err = np.sqrt(err_max**2 + err_15**2)
            if len(idx_minima) == 2:
                print ("There is a secondary maximum in I-band")
                jd_max2 = x_pred[idx_minima[1]]
                mag_max2 = res[idx_minima[1]]
                err_max2 = sig_pred[idx_minima[1]]  
            else:
                print ("There are no secondary maximum in I-band")
   
            
        # Save the predictions as pickle file
        if save_pred:
            interp_df = pd.DataFrame(list(zip(x_pred, res, sig_pred)), columns=['JD', 'FMAG', 'FERR'])
            output_filename = save_file+'interp_'+str(self.band)+'.pkl'
            if os.path.exists(output_filename):
                os.remove(output_filename)
            interp_df.to_pickle(output_filename)
            
        if verbose:
            print ("You are working on %s"%self.supernova)
            print ("The following are the interpolated values for the band %s "%self.band)
            print ("JD at maximum = %f "%jd_max)
            print ("Mag at maximum = %f +/- %f"%(mag_max, err_max))
            print ("dm15 = %f +/- %f"%(dm15, dm15_err))
            
        if plot_interp:
            
            fig = plt.figure(figsize=(6,6))
            ax = fig.add_subplot(111)
            ax.invert_yaxis()
        
            ax.plot(x_pred, res, 'k-', label='Prediction')
            ax.errorbar(self.x, self.y, self.yerr, fmt='r.', ms=10, 
                    label='Observation')
            
            if not find_two_peaks:
                ax.plot(jd_max, mag_max, '*', color='orange', ms=12)
            elif len(idx_minima)==2:
                ax.plot(jd_max, mag_max, '*', color='orange', ms=12)
                ax.plot(jd_max2, mag_max2, '*', color='orange', ms=12)
            else:
                ax.plot(jd_max, mag_max, '*', color='orange', ms=12)
            
            ax.fill_between(x_pred, res-sig_pred, res+sig_pred,
               alpha = 0.6, fc = '#a7c957', ec='None', label = '1-sigma')
        
            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize = 15)
                tick.label1.set_fontweight('bold')
    
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize = 15)
                tick.label1.set_fontweight('bold')
            ax.legend()
            
            
            
        if not find_two_peaks:
            return x_pred, res, sig_pred, jd_max, mag_max, err_max, jd_15, mag_15, err_15, dm15, dm15_err
        elif len(idx_minima)==2:
            return x_pred, res, sig_pred, jd_max, mag_max, err_max, jd_max2, mag_max2, err_max2, jd_15, mag_15, err_15, dm15, dm15_err
        else: 
            return x_pred, res, sig_pred, jd_max, mag_max, err_max, jd_15, mag_15, err_15, dm15, dm15_err
        
    
    
    def gp_fit(self, yerr, amp=None, scale=None, diff_deg=None,
               find_two_peaks=False):
    
        '''
        Remarks: Interpolate Light Curve with Gaussian Process Regression.
        x, y, yerr: the observed data array. 
        amp: amplitude of kernel (typically set by the amount by which the function
        varies.default is set by the standard deviation of the magnitude.)
        scale: scale over which the function typically varies.(default is 10)
        diff_deg: degree of differentiability. Controls the smoothness of the function.
    
        returns: jd at maximum.
    
        '''
    
        mean = lambda x: x*0 + np.median(self.y)
    
        if amp is None:
            amp = np.std(self.y - mean(self.x))
        elif amp == 'LC':
            amp = np.std(self.y)
        if scale is None:
            scale = 10
        if diff_deg is None:
            diff_deg = 1
    
        kernel = ConstantKernel(amp, constant_value_bounds='fixed') * Matern(length_scale=scale,
                                                                        length_scale_bounds='fixed',
                                                                        nu=diff_deg+0.5)
    
        gp = GaussianProcessRegressor(kernel=kernel, alpha=yerr**2)
        X = np.array([self.x]).T
        Y = self.y - mean(self.x)
        gauss_process = gp.fit(X, Y)
    
        ## Make predictions at points
        list_jd = np.atleast_2d(np.linspace(self.object_df['JD'].min()-8, self.object_df['JD'].max(), 1000)).T
    
        ## Make predictions
        x_pred = np.ravel(list_jd)
        res, sig_pred = gp.predict(list_jd, return_std=True)
        res += mean(x_pred)
    
        n = len(res)
        idx_minima = self.find_local_minima(res, n)
        if not find_two_peaks:
            jd_max = x_pred[idx_minima[0]]
        else:
            jd_max = x_pred[idx_minima[0]]
            jd_max2 = x_pred[idx_minima[1]]
        
        if not find_two_peaks:
            return jd_max
        else:
            return (jd_max, jd_max2)
        
        
    
    def monte_carlo(self, mc_trials, save_results,
                    store_params=False, kcorr=False, find_two_peaks=False, 
                   ):
    
        '''
        df: SN LC dataframe
        mc_trials: no of monte carlo realizations of the covariance of the LC
        band: Filter
        store_params: Bool to store parameters in a text file
        kcorr: Bool to store k-corr or not.
    
        returns: None
        
        '''
    
        
        param_value = self.interpolate_lc(verbose=False, amp=None, scale=30, diff_deg=2, plot_interp=False,
                              find_two_peaks=find_two_peaks)
    
        if not find_two_peaks:
            x_pred = param_value[0]
            res = param_value[1]
            sig_pred = param_value[2]
            jd_max = param_value[3]
            mag_max = param_value[4]
            err_max = param_value[5]
            jd_15 = param_value[6]
            mag_15 = param_value[7]
            err_15 = param_value[8]
            dm15 = param_value[9]
            dm15_err = param_value[10]
        else:
            x_pred = param_value[0]
            res = param_value[1]
            sig_pred = param_value[2]
            jd_max = param_value[3]
            mag_max = param_value[4]
            err_max = param_value[5]
            jd_max2 = param_value[6]
            mag_max2 = param_value[7]
            err_max2 = param_value[8]
            jd_15 = param_value[9]
            mag_15 = param_value[10]
            err_15 = param_value[11]
            dm15 = param_value[12]
            dm15_err = param_value[13]
       
        jd_max_list = []
        jd_max2_list = []
        np.random.seed(143)
        for trails in range(mc_trials):
            y_err = np.random.normal(scale=self.yerr, size=np.size(self.y))
            if not find_two_peaks:
                jd_max_list.append(self.gp_fit(y_err, scale=30, diff_deg=2))
            else:
                jd_max_list.append(self.gp_fit(y_err, scale=30, diff_deg=2)[0])
                jd_max2_list.append(self.gp_fit(y_err, scale=30, diff_deg=2)[1])
        
        if not find_two_peaks:
            jd_ = self.iterate_outliers(self.remove_outliers(jd_max_list))
            jd_max_err = np.std(jd_)
        else:
            jd_ = self.iterate_outliers(self.remove_outliers(jd_max_list))
            jd_max_err = np.std(jd_)
            jd2_ = self.iterate_outliers(self.remove_outliers(jd_max2_list))
            jd_max2_err = np.std(jd2_)
        
        if store_params:
            
            file_name = save_results+'Params_'+str(self.band)+'.txt'
            if os.path.exists(file_name):
                os.remove(file_name)
            f = open(file_name, 'w')
            
            f.write('This file contains the LC params of %s in filter-%s\n\n'%(self.supernova, self.band))
            f.write("#-----------------------------#\n\n")
            f.write("#-----------------------------#\n")
            f.write("Jd at maximum in %s band is: %f +/- %f\n" %(self.band, jd_max, jd_max_err))
            f.write("Magnitude at maximum in %s band is: %f +/- %f\n" %(self.band, mag_max, err_max))
            f.write("#-----------------------------#\n")
            f.write("#-----------------------------#\n")
            f.write("JD after 15-days from %s band maximum is: %f\n" %(self.band, jd_15))
            f.write("Magnitude after 15-days from %s band maximum is: %f +/- %f\n" %(self.band, mag_15, err_15))
            f.write("Estimated value of dm15 is: %f +/- %f\n" %(dm15, dm15_err))
            if not find_two_peaks:
                f.write("Secondary maximum is not applicable for this band\n")
            else:
                f.write("The values for the secondary maxima are\n")
                f.write("#-----------------------------#\n")
                f.write("JD for secondary maxima for %s band is: %f+/-%f\n" %(self.band, jd_max2, jd_max2_err))
                f.write("Magnitude for secondary maxima for %s band is: %f +/- %f\n" %(self.band, mag_max2, err_max2))
                f.write("#-----------------------------#\n")
            
            if not kcorr:
                f.write("The magnitudes are not k-corrected\n")
            
            f.write("Interpolated Light Curve through GP\n\n")
            f.write("#-----------------------------#\n\n")
            f.write("{0:>4s}{1:>20s}{2:>15s}\n\n".format('JD', 'MAG', 'MERR'))
            for l in zip(x_pred, res, sig_pred):
                f.write("{0:>4f}{1:>20f}{2:>15f}\n\n".format(*l))
        
        
            f.close()
        
    def save_file(self, file_name, z, EB_V_gal, EB_V_galerr, EB_V_host, EB_V_hosterr,
                 DM, DMerr):
        
        f = open(file_name, 'w')
        f.write("You are working on object %s\n"%Object)
        f.write("Redshit = %f\n"%z)
        f.write("E(B-V)G = %f +/- %f\n"%(EB_V_gal, EB_V_galerr))
        f.write("E(B-V)h = %f +/- %f\n"%(EB_V_host, EB_V_hosterr))
        f.write("Distance modulus = %f +/- %f"%(DM, DMerr))
        
        f.close()

#==========================================================================================================#