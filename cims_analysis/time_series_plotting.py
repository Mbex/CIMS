from bokeh.plotting import figure, output_notebook, show
from bokeh.io import gridplot, output_file, show
from bokeh.models import HoverTool, ColumnDataSource, OpenURL, DatetimeTickFormatter, FixedTicker, Range1d, LinearAxis
from bokeh.models.tools import Tool
from bokeh.palettes import Viridis256, Spectral11
from bokeh.charts import BoxPlot



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


def diurnal_ts(df, y):
    
    """Quick time series plot of y in df."""

    p = figure(
        x_axis_type="datetime",
        plot_width=900,
        plot_height=500,
       title = y + " Diurnal profile",
    )

    p.yaxis.axis_label = "["+y+"]"
    
    p.xaxis.formatter=DatetimeTickFormatter(formats=dict(
            hours=["%H:%M"],
             days=["%H:%M"],
             months=["%H:%M"],
             years=["%H:%M"]
        ))
    p.title.text_font_size = '24pt'    
    p.yaxis.axis_label_text_font_size = '14pt'    
    p.yaxis.major_label_text_font_size = '14pt'    
    p.xaxis.major_label_text_font_size = '14pt'    


    #confidence intervals
    upperband = df[y]["75%"][:-1]
    lowerband = df[y]["25%"][:-1]
    band_y = np.append(lowerband, upperband[::-1])
    x_axis = df.index[:-1] 
    band_x = np.append(x_axis, x_axis[::-1])

    p.patch(x = band_x , y = band_y, fill_alpha = 0.2, line_alpha= 0.2, line_color = "blue")
    p.line(x=df.index, y=df[y]["mean"], line_color = "blue")

    show(p)


def diurnal_ts_many(df, y, mean=True):
    
    """Quick time series plot of y's in df."""

    p = figure(
        x_axis_type="datetime",
        plot_width=900,
        plot_height=700,
        title = " Diurnal profiles",
        y_axis_type="log",
   
    )

    
    p.title.text_font_size = '24pt'    
    p.yaxis.axis_label_text_font_size = '14pt'    
    p.yaxis.major_label_text_font_size = '14pt'    
    p.xaxis.major_label_text_font_size = '14pt'    
  
    xaxis = df.index[:-1].map(lambda x: int(x.strftime('%s')))
    
    mypalette = Viridis256
    colour_step = range(0, len(mypalette), len(mypalette) / len(y) )
    
    max_val = 1e-12
    min_val = 1e12
    for i, name in enumerate(y):
        
        if mean:
            y_trace = list(df[name]["mean"][:-1]) / np.mean(list(df[name]["mean"][:-1]))
            p.yaxis.axis_label = "normalised to mean value"  
        else:
            y_trace = df[name]["mean"][:-1]
            p.yaxis.axis_label = "cps"  
            
            
        p.line(df.index,
               y_trace,
               color = mypalette[colour_step[i]],
               legend = name)
       
        if np.max(y_trace) > max_val:
            max_val = np.max(y_trace)
            
        if np.min(y_trace) < min_val:
            min_val = np.min(y_trace)
        
    p.y_range= Range1d(min_val * 0.9, max_val * 1.1)

    p.xaxis.formatter=DatetimeTickFormatter(formats=dict(
        hours=["%H:%M"],
         days=["%H:%M"],
         months=["%H:%M"],
         years=["%H:%M"]
    ))
    
    p.legend.location = "bottom_right"
    p.legend.orientation = "horizontal"
    show(p)

    return p
