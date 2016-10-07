class time_series_analysis(object):
    
    def __init__(self):
        
        self.igor_time = 2082844800


    def time_series(self, df, y):

        """
        Plotly scatter plot of X and y.
        Requires dt.datetime 'date' column.
        """

        import plotly
        from plotly.graph_objs import Scatter, Layout

        return plotly.offline.iplot({
                "data" : [
                    Scatter(x = df["date"],
                            y = df[y]
                    )
                ],

                "layout" :
                    Layout(yaxis = dict( title = y),
                           xaxis = dict( title = "date")
                    ) 
        })


    def scatter(self, df, x, y):

        """Plotly scatter plot of x and y."""

        import plotly
        from plotly.graph_objs import Scatter, Layout

        return plotly.offline.iplot({
                "data" : [
                    Scatter(x = df[x],
                            y = df[y],
                            mode = "markers"
                    )
                ],

                "layout" :
                    Layout(title = "%s vs %s" % (x, y),
                           xaxis = dict( title = x),
                           yaxis = dict( title = y)
                    ) 
        })


    def linear_plot_params(self, x, y):

        """Calculate plot params."""

        mask = ~np.isnan(x.astype(float)) & ~np.isnan(y.astype(float))
        return sc.stats.linregress(x[mask], y[mask])


    def read_time_series_text_to_dataframe(self, path):

        """Joins all time series text files in a dir into one df."""

        allFiles = glob.glob(path + "/*.txt")
        frame = pd.DataFrame()
        list_ = []
        for file_ in allFiles:
            df = pd.read_csv(file_,index_col=None, header=0)
            list_.append(df)
        frame = pd.concat(list_, axis = 1)

        return frame



    def igor_time_to_posix(self, df, igor_time_column_name):

        "Converts igor time in seconds to datetime object."

        import datetime as dt
        import pandas as pd

        returndf = df.copy()
        
        new_time_column = pd.Series([dt.datetime.fromtimestamp(np.round(x) - self.igor_time) for x in df[igor_time_column_name]])
        returndf['date'] = new_time_column
        returndf.set_index('date')

        return returndf


    def remove_background(self, df, start_time, end_time):

        """
        Minus average value of time series found between bg_start and bg_end for
        the rest of the time series. Requires datetime "date" column to be present.
        """

        import datetime as dt

        
        returndf = df.copy()

        for time in [start_time, end_time]:
            if not isinstance(time, dt.datetime):
                raise TypeError ("%s is not a datetime object" % str(time))
                break

        start_index = df["date"].searchsorted(start_time)[0]
        end_index = df["date"].searchsorted(end_time)[0]

        for i, col in enumerate(df.columns):

            if col == "date":
                pass
            else:
                try:
                    bg_period = df.loc[start_index:end_index, col]
                    average_bg = bg_period.mean()
                    
                    returndf[col] = df[col].subtract(average_bg)
                except AttributeError:
                    pass
                except TypeError:
                    pass


        return returndf


    def nan_points_between_dates(self, df, start_time, end_time):

        """
        Nan all entries in dataframe between start and end time
        """

        import datetime as dt
        
        df_return = df.copy()

        for time in [start_time, end_time]:
            if not isinstance(time, dt.datetime):
                raise TypeError ("%s is not a datetime object" % str(time))
                break

        mask = (df['date'] > start_time) & (df['date'] <= end_time)
        df_return.loc[mask] = None

        return df_return


    def remove_start(self, df, start_time):

        """Removes dataframe before start_time."""

        df_return = df.copy()
        return df_return[df['date'] > start_time]


    def remove_end(self, df, end_time):

        """Removes dataframe after end_time."""

        df_return = df.copy()
        return df_return[df['date'] < end_time]

    

    def normalise(self, df, x, y):

        """Normalise column x in df by column y."""    

        df_return = df.copy()

        slope = self.linear_plot_params(df[x], df[y]).slope
        normalisation_factor = df[x] * slope

        df_return[y] = (df[y]  / normalisation_factor).multiply(np.mean(normalisation_factor))

        return df_return

