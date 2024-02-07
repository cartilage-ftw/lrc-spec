import matplotlib.pyplot as plt
import io_utils

import plotly.express as px
import plotly.graph_objects as go
#import plotly.plotly as py

from bokeh.plotting import figure


def make_spectrum_fig(data, file_name=''):
    tooltips = [("(x,y)", "($x{0.000},$y)")]
    fig = figure(x_axis_label='Wavenumber [cm-1]', y_axis_label='MS Fraction %',
                    height=600, width=800, tools='tap,hover,pan,reset,lasso_select,box_zoom,save',
                     tooltips=tooltips, active_inspect='hover')
    fig.line(data['Wavenumber'], data['MS Fraction'], color='gray', line_width=2.0)
    fig.scatter(data['Wavenumber'], data['MS Fraction'], size=10, #x_error=data['Wavenumber_err'],
                    color='deeppink')
    fig.xgrid.visible = False
    fig.ygrid.visible = False

    return fig