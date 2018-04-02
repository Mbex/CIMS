#!/usr/bin/python


#from bokeh.plotting import figure, output_notebook, show
#from bokeh.io import gridplot, output_file, show
#from bokeh.models import HoverTool, ColumnDataSource, OpenURL, DatetimeTickFormatter, FixedTicker, Range1d, LinearAxis
#from bokeh.models.tools import Tool
#from bokeh.palettes import Viridis256, Spectral11
#from bokeh.charts import BoxPlot
import glob
import pandas as pd
import datetime as dt
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt 
import matplotlib
import copy
import random
import operator
import timeit
import math
import re
import collections
import itertools
import matplotlib.dates as mdates




def quick_ts(x, y):




    tools="pan,wheel_zoom,box_zoom,reset,hover"

#    hover = HoverTool()
#    hover.tooltips = [
#        ("","@y"),
#    ]
#    TOOLS = [hover]

    # make figure
    p = figure(
        x_axis_type="datetime",
        plot_width=900,
        plot_height=500,
        tools=tools #TOOLS
    )


    p.yaxis.axis_label = "counts"
    p.xaxis.formatter=DatetimeTickFormatter(formats=dict(
            hours=["%d %b"],
             days=["%d %b"],
             months=["%d %b"],
             years=["%d %b"]
        ))

    p.title.text_font_size = '24pt'    
    p.yaxis.axis_label_text_font_size = '14pt'    
    p.yaxis.major_label_text_font_size = '14pt'    
    p.xaxis.major_label_text_font_size = '14pt'  
    p.line(x, y, line_color="blue")

    show(p)
    
    
def quick_ts_many(df, y):
    p = figure(
        x_axis_type="datetime",
        plot_width=900,
        plot_height=500,
    )

    p.yaxis.axis_label = "["+y[0]+"]"
    p.xaxis.formatter=DatetimeTickFormatter(formats=dict(
            hours=["%d %b"],
             days=["%d %b"],
             months=["%d %b"],
             years=["%d %b"]
        ))

    p.title.text_font_size = '24pt'    
    p.yaxis.axis_label_text_font_size = '14pt'    
    p.yaxis.major_label_text_font_size = '14pt'    
    p.xaxis.major_label_text_font_size = '14pt'  
    
    mypalette = Viridis256
    colour_step = range(0, len(mypalette), len(mypalette) / len(y) )
    for i, name in enumerate(y):

        if "date" not in df.columns:
            x_axis = df.index.values
            print "using datetime index"
        else:
            x_axis = df['date']
        
        p.line(x=x_axis,
               y=df[name],
               color = mypalette[colour_step[i]],
               legend = name)
    
    
    p.legend.location = "bottom_right"
    p.legend.orientation = "horizontal"

    show(p)



def plot_isotope_ratios(df, x35, x37, x37resolution, LOD35, LOD37, av_ratio=3.13, **kwargs):
  
    df_plot = df.copy()
    fig = plt.figure(figsize=(25,15))
    fig.suptitle(x35+" isotopes")
    gridsize = (3,3)
    
    c1 = "blue"
    c2 = "red"
    c3 = "m"
    c4 = "orange"
    # process data for plotting
    # throw away below LOD
    df_plot.loc[df_plot[x37] < LOD37, x37] = np.nan
    df_plot.loc[df_plot[x35] < LOD35, x35] = np.nan

    # throw away data where resolution is too low.
#     try:
#         # minimum allowed delta mass (float / series = series)
#         dm_min = kwargs["x37mass"] / df[x37resolution]
#         # actual delta mass (float - float = float)
#         dm = abs(kwargs["x37mass"] - kwargs["interference_mass"])
# #         df_plot.loc[dm < dm_min, x37] = np.nan
        
#         ax2 = plt.subplot2grid(gridsize, (2, 0), colspan=2)       
#         ax2.plot(df_plot.index, dm_min, label = "Min resolution for peak seperation")
#         ax2.set_ylabel("m/Q (Th)")
#         ax2.plot(df_plot.index, [dm if (x > dm_min.first_valid_index() and x < dm_min.last_valid_index()) else np.nan for x in df_plot.index],
#                  c ="k", ls= "--", label = "dm") # average Cl isotope ratio
#         ax2.legend()
        
#     except KeyError:
#         print "kwargs not given"

    mask = ~np.isnan(df_plot[x35]) & ~np.isnan(df_plot[x37])
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(df_plot[x37][mask], df_plot[x35][mask])
    df_plot["ratio"] = (df_plot[x35][mask] + (~int(intercept)+1)) / df_plot[x37][mask]

    # time series
    ax = plt.subplot2grid(gridsize, (0, 0), colspan=2)       
    [
      ax.axvspan(df_plot.index[i], df_plot.index[i+1], facecolor='0.2', alpha=0.2) 
      for i in range(len(df_plot.index)-1) if df_plot.index[i] > df_plot[x35].first_valid_index() and 
      df.index[i] < df_plot[x35].last_valid_index() and df_plot.nighttime[i]
    ]   
    ax.plot(df_plot[x35],c1)
    ax.plot(df_plot[x37],c2)
    ax.plot(df_plot.index, [LOD35 if (x > df_plot[x35].first_valid_index() and x < df_plot[x35].last_valid_index()) else np.nan for x in df_plot.index], c = c1, ls="--", label = " 35Cl LOD = "+str(int(LOD35))) 
    ax.plot(df_plot.index, [LOD37 if (x > df_plot[x37].first_valid_index() and x < df_plot[x37].last_valid_index()) else np.nan for x in df_plot.index], c = c2, ls="--", label = " 37Cl LOD = "+str(int(LOD37))) 
    ax.legend()
   
    #time series of ratio
    ax1 = plt.subplot2grid(gridsize, (1, 0), colspan=2)       
    [
      ax1.axvspan(df_plot.index[i], df_plot.index[i+1], facecolor='0.2', alpha=0.2) 
      for i in range(len(df_plot.index)-1) if df_plot.index[i] > df_plot["ratio"].first_valid_index() and 
      df_plot.index[i] < df_plot["ratio"].last_valid_index() and df_plot.nighttime[i]
    ]
    # measured ratio
    ax1.plot(df_plot["ratio"],"k", label="measured ratio")
    # measured average ratio
    ax1.plot(df_plot.index, [np.nanmean(df_plot["ratio"]) if (x > df_plot["ratio"].first_valid_index() and x < df_plot["ratio"].last_valid_index()) else np.nan for x in df_plot.index], c =c3, ls= "--", label = "average measured 35Cl/37Cl ratio") # average Cl isotope ratio
    # average ratio
    ax1.plot(df_plot.index, [av_ratio if (x > df_plot["ratio"].first_valid_index() and x < df_plot["ratio"].last_valid_index()) else np.nan for x in df_plot.index], c =c4, ls= "--", label = "average 35Cl/37Cl ratio") # average Cl isotope ratio
    ax1.set_ylabel(x35+"/"+x37) 
    ax1.legend()

    # 37 isotope vs ratio
    ax2 = plt.subplot2grid(gridsize, (2, 0), colspan=1)       
    ax2.scatter(df[x37], df_plot["ratio"], c = "purple")
    ax2.set_xlabel(x37)
    ax2.set_ylabel(x35+"/"+x37) 
    ax2.plot(df_plot[x37], [av_ratio if (x > df_plot["ratio"].first_valid_index() and x < df_plot["ratio"].last_valid_index()) else np.nan for x in df_plot.index], c =c4, ls= "--", label = "average 35Cl/37Cl ratio") # average Cl isotope ratio

    # MCP_monitor vs ratio
    ax2 = plt.subplot2grid(gridsize, (2, 1), colspan=1)       
    ax2.scatter(df["MCP_monitor"], df_plot["ratio"], c = c3)
    ax2.set_xlabel("MCP_monitor")
    ax2.set_ylabel(x35+"/"+x37) 
    ax2.plot(df_plot[x37], [av_ratio if (x > df_plot["ratio"].first_valid_index() and x < df_plot["ratio"].last_valid_index()) else np.nan for x in df_plot.index], c =c4, ls= "--", label = "average 35Cl/37Cl ratio") # average Cl isotope ratio    
    ax2.set_xlim((1825,1900))

    # resolution vs ratio
    ax3 = plt.subplot2grid(gridsize, (2, 2), colspan=1)       
    ax3.scatter(df_plot[x37resolution], df_plot["ratio"], c = c3)
    ax3.set_xlabel(x37resolution)
    ax3.set_xlim((2000, 4000))
    ax3.set_ylabel(x35+"/"+x37)
    ax3.plot(df_plot[x37resolution], [av_ratio if (x > df_plot["ratio"].first_valid_index() and x < df_plot["ratio"].last_valid_index()) else np.nan for x in df_plot.index], c =c4, ls= "--", label = "average 35Cl/37Cl ratio") # average Cl isotope ratio
    linear_trendline(ax3, df_plot[x37resolution], df_plot["ratio"], c='k')

    # 37 vs 35 scatter
    ax4 = plt.subplot2grid(gridsize, (0, 2), colspan=1, rowspan=2)       
    ax4.plot([0, np.nanmax(df_plot[x37])[0]],[0, np.nanmax(df_plot[x37])[0]*av_ratio], c ="lime", ls= "--", label="expected regression")
    ax4.scatter(df_plot[x37], df_plot[x35], c=c3, label="observed")
    linear_trendline(ax4, df_plot[x37], df_plot[x35], c='k')
    ax4.set_xlabel(x37)
    ax4.set_ylabel(x35)
    ax4.set_xlim((0,np.nanmax(df_plot[x37])[0]))
    ax4.set_ylim(0,np.nanmax(df_plot[x37])[0]*av_ratio)
    ax4.legend(loc=4)



def linear_trendline(ax, xd, yd, c='r', alpha=1):
    """Make a line of best fit"""

    #Calculate trendline
    mask = ~np.isnan(xd) & ~np.isnan(yd)
    
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(xd[mask], yd[mask])
    
    minxd = np.min(xd)
    maxxd = np.max(xd)

    xl = np.array([minxd, maxxd])
    xl = np.array([0, maxxd])
    yl = 0 * xl ** 2 + slope * xl + intercept

    #Plot trendline
    ax.plot(xl, yl, c, alpha=alpha, label = "observed trend")
    
    # Put on graph
#    ax.text((0.1) * maxxd + (0.90) * minxd, (0.90) * np.max(yd)+ (0.1) * np.min(yd),
#         '$R^2 = %0.2f$' % r_value**2, color = c)
#    ax.text((0.1) * maxxd + (0.85) * minxd, (0.85) * np.max(yd) + (0.1) * np.min(yd),
#          '$Y = %0.2f  X + %0.2f$' % (slope, intercept), color = c)
#    ax.text((0.1) * maxxd + (0.80) * minxd, (0.80) * np.max(yd) + (0.1) * np.min(yd),
#          '2sigma = %0.2f' % (2 * (std_err * np.sqrt(len(xd[mask])))), color = c)

#    ax.annotate( '$R^2 = %0.2f$' % r_value**2, xy=(0.2,0.9), xycoords='axes fraction',ha="right", va="bottom", color = c)



    return {"slope": slope, "intercept": intercept, "r_value": r_value, "p_value": p_value, "std_err": std_err}




def plot_diurnal_profile(df, ys, **kwargs):

    """"""
    
    # initialise plot
    fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
    
    try:
        fig.suptitle(kwargs["title"])
    except KeyError:
        fig.suptitle("Diurnal Profile")
    
    try:
        colors = kwargs["colors"]
    except KeyError:
        colors = itertools.cycle(('g','r', "b", "m","c","orange", "purple")) 

    try:
        markers = kwargs["markers"]
    except KeyError:
        markers = itertools.cycle(('s', '+', 'o', '*')) 

    # plot 
    try:
        ax1.plot(df.index, df[ys[0]]["mean"], c="k", marker="X", ls="--",label=ys[0])
    except KeyError:
        ax1.plot(df.index, df[ys[0]], c="k", marker="X", ls="--",label=ys[0])
        
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    if "J" in ys[0]:
        ax1.set_ylabel(ys[0] +" / s$^{-1}$")
    elif "global" in ys[0] or "direct" in ys[0] or "terrestrial" in ys[0] or "indirect" in ys[0]: 
        ax1.set_ylabel(ys[0] +" / W$^{-2}$")
    elif "I." in ys[0]:
        ax1.set_ylabel(ys[0].replace("I.","")+" / ppt")
    else:
        ax1.set_ylabel(ys[0])


    if len(ys) > 1:
        for i, y in enumerate(ys[1:]):
            ax = ax1.twinx()
	    try:
                c = colors[i]
	    except TypeError:
		c = colors.next()
            m = markers.next()
            try:
                ax.plot(df.index, df[y]["mean"], c=c, marker=m, ls="--",label=y)
            except KeyError:
                ax.plot(df.index, df[y], c=c, marker=m, ls="--",label=y)

            # y axis
            ax.spines['top'].set_visible(False)
            ax.spines["right"].set_position(('axes',1 + (0.1* i)))
            ax.set_ylabel(""+y+"".replace("I.",""),color=c)
            
            if "J" not in y:
                ax.set_ylabel(""+y+"".replace("I.","")+" / ppt",color=c)
            else:
                ax.set_ylabel(""+y+"".replace("I.","")+" / s-1",color=c)


    # x axis
#    max_xticks = 4
#    xloc = plt.MaxNLocator(max_xticks)
    ax1.xaxis.set_major_locator(mdates.HourLocator(interval=6))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#    ax1.xaxis.set_ticklabels(["00","01","02","03","04","05","06","07","08","09","10",
#   "11","12","13","14","15","16","17","18","19","20","21","22","23"])



def pairs(df, ys, **kwargs):

    """
    Plot to see relationships between species in DF, ys.
    df: dataframe.
    ys: array of column names from dataframe.
    """

    fig = plt.figure(figsize=(len(ys)*7, len(ys)*7))
    gridsize = (len(ys), len(ys))


    for i, gx in enumerate(range(gridsize[0])):
        for j, gy in enumerate(range(gridsize[1])):
            ax = plt.subplot2grid(gridsize, (gx, gy), colspan=1)       

            if i != j:
#            if i > j:

                x_c = ys[i]
                y_c = ys[j]

                x = df[x_c]
                y = df[y_c]
                
                # Set color
                try:
                    c = kwargs['c']
                except KeyError:
                    c = "blue"
                    
                # if dataframe is diurnal choose mean
                try:
                    x = x["mean"]
                    y = y["mean"]
                except KeyError:
                    pass
 
                ax.scatter(x, y, c=c)

                plot_params = linear_trendline(ax, x, y)

                ax.set_ylabel(y_c)
                ax.set_xlabel(x_c)
                ax.set_ylim((0, np.nanmax(y)[0]))
                ax.set_xlim((0, np.nanmax(x)[0]))
            else:
                ax.set_visible(False)
