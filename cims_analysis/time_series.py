#!/usr/bin/python

import pandas as pd
import glob
import datetime as dt
import numpy as np
import scipy.stats
import operator
import re
import collections



class time_series_analysis(object):
    

    def __init__(self, reagent_ion):
       

        self.igor_time = 2082844800
        self.reagent_ion = reagent_ion
        self.time_series_raw = {}
        self.time_series_normalised = {}
        self.time_series_normalised_background = {}
        self.date = pd.Series()
        self.diurnal_time = []
        self.diurnal_time_series = {}
        self.correlation_dataframe = pd.DataFrame()
        self.time_series_conc = {}
        self.time_series_whitworth = {}
        self.background = {}
        self.background_conc = {}
        
    def linear_plot_params(self, x, y):

        """Calculate plot params."""

        mask = ~np.isnan(x.astype(float)) & ~np.isnan(y.astype(float))
        #mask = ~np.isnan(x) & ~np.isnan(y) 
	return scipy.stats.linregress(x[mask], y[mask])
        #return scipy.stats.linregress(x, y)


    def read_time_series_text_files(self, path):

        """
        Reads all files with *.txt in a dir assuming they are time series.
        Puts these time series into a dictionary with key of the ts name.
	Time series are a numpy array of numpy.float64s
        """

        print "reading in data to self.time_series_raw"

        allFiles = glob.glob(path + "/*.txt")
        for f in allFiles:
            ts = pd.read_csv(f)
            self.time_series_raw[ts.columns[0]] = np.array([x[0] for x in ts.values])

        
    def time_to_posix(self, time_column_name, igor_time=True):

        """
	Converts time in seconds to datetime object as self.'date' from time column name in self.time_series_raw.
	option to convert igor time offset.
	"""
        
        print "converting time"
        
        self.time_column_name = time_column_name
	if igor_time:
            new_time_column = pd.Series([dt.datetime.fromtimestamp(np.round(x) - self.igor_time) for x in self.time_series_raw[time_column_name]])
        else:
            new_time_column = pd.Series([dt.datetime.fromtimestamp(np.round(x)) for x in self.time_series_raw[time_column_name]])

        self.date = new_time_column


    def nan_zeros(self):

	"""NaN anything below 0 in traces."""

	for key in self.time_series_raw.keys():
   	    self.time_series_raw[key] = np.array([np.nan if i <= 0 else i for i in self.time_series_raw[key]])


    def nan_points_between_dates(self, start_time, end_time):

        """NaN all entries in self.time_series_raw between start and end time."""

        for time in [start_time, end_time]:
            if not isinstance(time, dt.datetime):
                raise TypeError ("%s is not a datetime object" % str(time))
                break

        for i, x in enumerate(self.date):
            if x < start_time:
	        sti = i+1
	    if x < end_time:
	        eti = i+1

	print "removing bad points", sti,":", eti

        for key in self.time_series_raw.keys():
            self.time_series_raw[key][sti:eti] = np.nan

    
    def remove_start(self, start_time):

        """Removes dataframe before start_time."""

	for i, x in enumerate(self.date):
 	    if x > start_time:
		print "removing start time ", str(start_time), " at index", i  
                break

	self.date = self.date[i:]  
        for key in self.time_series_raw.keys():
            self.time_series_raw[key] = self.time_series_raw[key][i:]

    
    def remove_end(self, end_time):

        """Removes dataframe after end_time."""

	for i, x in enumerate(self.date):
	    if x > end_time:
                print "removing end time ", str(end_time), " at index", i
	        break
	
	self.date = self.date[:i]        
        for key in self.time_series_raw.keys():
            self.time_series_raw[key] = self.time_series_raw[key][:i]    

 
    def normalise(self, reagent_ion=None, normalise_to=None):

        """
        Normalise all timeseries in self.time_series_raw by self.reagent_ion
        unless specificed in kwargs. Option to normalise counts to kwargs e.g. 1e6.
        """    
        
        print "normalising"

        if reagent_ion is None:
            reagent_ion = self.reagent_ion
    
        for y in self.time_series_raw.keys():
            if (y != self.time_column_name):          
            
                self.time_series_normalised[y] = np.nanmean(self.time_series_raw[reagent_ion]) *  (self.time_series_raw[y]  / self.time_series_raw[reagent_ion])


    def remove_background(self, start_time, end_time):

        """Removes background from time series."""

        print "remove background"
   
        start_index = self.date.tolist().index(start_time)
        end_index = self.date.tolist().index(end_time)

        for y in self.time_series_raw.keys():       
            if (y != self.time_column_name) and (y != self.reagent_ion):   
#                 mn = np.nanmin(self.time_series_normalised[y][start_index:end_index])
                self.background[y] = np.nanmin(self.time_series_normalised[y][start_index:end_index])
                if np.isnan(self.background[y]):
                    pass
                else:
                    self.time_series_normalised_background[y] = self.time_series_normalised[y] - self.background[y] 
    
#                 self.time_series_normalised_background[y] = self.time_series_normalised[y] - mn
        
#                 np.where(self.time_series_normalised_background[y] < 0, self.time_series_normalised_background[y], 0)

    
    def remove_bg(self):

        """Removes background from time series."""
    
        print "remove background"
   
        for y in self.time_series_raw.keys():       
            if (y != self.time_column_name) and (y != self.reagent_ion):                
                mn = 1.12 * np.nanmin(self.time_series_normalised[y])
                self.time_series_normalised_background[y] = self.time_series_normalised[y] - mn


                
    def get_reagent_ion_r2s(self, list_of_reagent_ions, y):
        
        "Returns a dict with reagent ion keys and R2 values to y."""
    
        #print "Reminder: These correlations are performed on the raw time series."
    
        trace_and_norm = {}       
        reagent_ions = {}
        for reagent_ion in list_of_reagent_ions:
            reagent_ions[reagent_ion] = round(self.linear_plot_params(self.time_series_raw[reagent_ion], self.time_series_raw[y]).rvalue**2, 4)
        return reagent_ions
    
    
    def best_reagent_ion(self, list_of_reagent_ions, y, R2_limit):  
        
        "Returns the key of the highest R2 value in the dict from self.norm_ion_r2s."""

        reagent_ions = self.get_reagent_ion_r2s(list_of_reagent_ions, y)  
        best = max(reagent_ions.iteritems(), key=operator.itemgetter(1))[0]

        if reagent_ions[best] < R2_limit:
            best = "best R2 of %f with %s is less than %f" % (reagent_ions[best], best, R2_limit)

        return best 
    
    
    def generate_correlation_dataframe(self):
        
        """Once data is prepared. Run this to get self.correlation_dataframe"""

        print "generating correlation dataframe on time_series_whitworth"
        
        result = {}
        for key in self.time_series_whitworth:
            result[key] = self.time_series_whitworth[key].astype(float)#.tolist()
        self.correlation_dataframe = pd.DataFrame(result).corr()
    
    
    def top_correlations(self, y, n):

        """Returns list of top n correlations for species y."""
      
        #list top 20 wave names that correlate
        if len(self.correlation_dataframe) == 0:
            self.generate_correlation_dataframe()
            
        names = list(self.correlation_dataframe[y].sort_values(ascending = False).head((n+1)).index)
        R2s = list(self.correlation_dataframe[y].sort_values(ascending = False).head((n+1)))
        top_hits = list(zip(names, R2s))
        
        return  top_hits[1:]
    
    
    def diurnalise(self, y):  

        """
        Returns descriptive stats on waves in a df.

        Must have 'date' set as datetime index i.e.
        df = df.set_index(pd.DatetimeIndex(df["date"]))
        """
        
        print "Generating Diurnal DataFrame for %s" % y
        
        ts_df = pd.DataFrame()
        ts_df['date'] = self.date.dt.strftime("%H:%M")
        ts_df = ts_df.set_index(pd.DatetimeIndex(self.date.dt.strftime("%H:%M")))

        ts_df[y] = self.time_series_normalised_background[y]

        df = ts_df.groupby('date').describe().unstack()
        df.index = pd.to_datetime(df.index.astype(str))

        returndf = df.copy()
        del df
        return returndf
    
    
    def calibrate(self, calibration_factors):
        
        
        """Applies calibration to self.time_series_normalised_background"""
        
        print "calibrating"

        ts_keys = self.time_series_normalised_background.keys()
        cf_keys = calibration_factors.keys()
        
        for ts_key in ts_keys:
            for cf_key in cf_keys:
                try:
                    self.time_series_conc[ts_key] = self.time_series_normalised_background[ts_key] / calibration_factors[ts_key]
                    self.background_conc[ts_key] = self.background[ts_key] / calibration_factors[ts_key]
               
                except KeyError:
                    
                    self.time_series_conc[ts_key] = self.time_series_normalised_background[ts_key]  / calibration_factors["C2H4IO2_"]
                    self.background_conc[ts_key] = self.background[ts_key] / calibration_factors["C2H4IO2_"]                    
                    
                    if "C" in self._Count_elements(ts_key).keys():
                        
                        nCarbons = self._Count_elements(ts_key)['C']  
                        
                        if nCarbons == 1:
                
                            self.time_series_conc[ts_key] = self.time_series_normalised_background[ts_key]  / calibration_factors["CH2IO2_"]
                            self.background_conc[ts_key] = self.background[ts_key] / calibration_factors["CH2IO2_"]

                    

                    
    def _Count_elements(self, moiety):

        """Returns dict with keys of elements and values of number of elements."""
         
        ret = []
        ans = re.findall(r'([A-Z][a-z]*)(\d*)', moiety)
        for i, value in enumerate(ans):
            val1 = value[1]
            if len(val1) == 0: val1 = 1
            ret.append((unicode(value[0]), int(val1)))

        return collections.OrderedDict(ret)
        

