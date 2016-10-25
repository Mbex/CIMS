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

        
        
    def linear_plot_params(self, x, y):

        """Calculate plot params."""

        mask = ~np.isnan(x.astype(float)) & ~np.isnan(y.astype(float))
        return scipy.stats.linregress(x[mask], y[mask])
        #return sc.stats.linregress(x, y)


    def read_time_series_text_to_dataframe(self, path):

        """
        Reads all files with *.txt in a dir assuming they are time series.
        Puts these time series into a dictionary with key of the ts name.
        """

        print "reading in data"

        allFiles = glob.glob(path + "/*.txt")
        for f in allFiles:
            ts = pd.read_csv(f)
            self.time_series_raw[ts.columns[0]] = ts.values
           
        
    def igor_time_to_posix(self, igor_time_column_name):

        """Converts igor time in seconds to datetime object as 'date' in self.time_series_raw."""
        
        print "converting time"
        
        self.igor_time_column_name = igor_time_column_name
        new_time_column = pd.Series([dt.datetime.fromtimestamp(np.round(x) - self.igor_time) for x in self.time_series_raw[igor_time_column_name]])
        self.date = new_time_column


    def nan_points_between_dates(self, start_time, end_time):

        """NaN all entries in self.time_series_raw between start and end time."""

        print "removing bad points"
        
        for time in [start_time, end_time]:
            if not isinstance(time, dt.datetime):
                raise TypeError ("%s is not a datetime object" % str(time))
                break

        mask = [i for i, x in enumerate(self.date) if (x > start_time) & (x < end_time)]
        for key in self.time_series_raw.keys():
            self.time_series_raw[key][mask] = None
        
    
    def remove_start(self, start_time):

        """Removes dataframe before start_time."""

        mask = [i for i, x in enumerate(self.date) if (x < start_time)]
        for key in self.time_series_raw.keys():
            self.time_series_raw[key][mask] = None

            
    def remove_end(self, end_time):

        """Removes dataframe after end_time."""

        mask = [i for i, x in enumerate(self.date) if (x > end_time)]
        for key in self.time_series_raw.keys():
            self.time_series_raw[key][mask] = None    
    
    
    def normalise(self, reagent_ion=None, normalise_to=None):

        """
        Normalise all timeseries in self.time_series_raw by self.reagent_ion
        unless specificed in kwargs. Option to normalise counts to kwargs e.g. 1e6.
        """    
        
        print "normalising"

        if not reagent_ion:
            reagent_ion = self.reagent_ion
    
        for y in self.time_series_raw.keys():
            if (y != self.igor_time_column_name) and (y != self.reagent_ion):          
                slope = self.linear_plot_params(self.time_series_raw[reagent_ion], self.time_series_raw[y]).slope            
                if normalise_to:
                    self.time_series_normalised[y] = (slope * self.time_series_raw[y] * normalise_to) / (self.time_series_raw[reagent_ion]  * slope)
                else: 
                    self.time_series_normalised[y] = (slope * self.time_series_raw[y] * np.nanmean(self.time_series_raw[reagent_ion])) / (self.time_series_raw[reagent_ion]  * slope)


    def remove_background(self, start_time, end_time):

        """Removes background from time series."""

        print "remove background"
        
        mask = [i for i, x in enumerate(self.date) if (x > start_time) & (x < end_time)]   
        for y in self.time_series_raw.keys():       
            if (y != self.igor_time_column_name) and (y != self.reagent_ion):                   
                mn = np.nanmin(self.time_series_raw[key][mask])
                self.time_series_normalised_background[y] = self.time_series_normalised[y] - mn
                [0 for x in self.time_series_normalised_background[y] if x < 0]

    
    def remove_bg(self):

        """Removes background from time series."""
    
        print "remove background"
   
        for y in self.time_series_raw.keys():       
            if (y != self.igor_time_column_name) and (y != self.reagent_ion):                
                mn = 1.12 * np.nanmin(self.time_series_normalised[y])
                self.time_series_normalised_background[y] = self.time_series_normalised[y] - mn
                [0 for x in self.time_series_normalised_background[y] if x < 0]
                
                
    def get_reagent_ion_r2s(self, list_of_reagent_ions, y):
        
        "Returns a dict with reagent ion keys and R2 values to y."""
    
        print "Reminder: These correlations are performed on the raw time series."
    
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
            best = "best R2 of %f with %s is less than %f" % (reagent_ions[best], best, round(R2_limit, 2))

        return best 
    
    def generate_correlation_dataframe(self):
        
        """Once data is prepared. Run this to get self.correlation_dataframe"""

        print "generating correlation dataframe on time_series_normalised_background"
        
        result = {}
        for key in bonfire.time_series_normalised_background:
            result[key] = bonfire.time_series_normalised_background[key].flatten().tolist()
        self.correlation_dataframe = pd.DataFrame(result).corr()
    
    
    
    def top_correlations(self, y, n):

        """Returns list of top n correlations for species y."""
      
        #list top 20 wave names that correlate
        if len(self.correlation_dataframe) == 0:
            self.generate_correlation_dataframe()
            
        names = list(bonfire.correlation_dataframe[y].sort_values(ascending = False).head((n+1)).index)
        R2s = list(self.correlation_dataframe[y].sort_values(ascending = False).head((n+1)))
        top_hits = zip(names, R2s)
        
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
    


