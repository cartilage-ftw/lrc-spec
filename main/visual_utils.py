import matplotlib.pyplot as plt
import io_utils

import plotly.express as px
import plotly.graph_objects as go
#import plotly.plotly as py

from bokeh.plotting import figure


"""
NOTE: deprecated
"""
def make_spectrum_fig_bokeh(data):
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


def make_spectrum_fig_plotly(data):
    fig = go.Figure()

    # Line
    fig.add_scatter(
        x=data["Wavenumber"],
        y=data["MS Fraction"],
        mode="lines",
        name="Spectrum (line)",
        line=dict(
            color="gray",
            width=2,
        ),
        hoverinfo="skip",  # hover handled by markers
    )

    # Markers
    fig.add_scatter(
        x=data["Wavenumber"],
        y=data["MS Fraction"],
        mode="markers",
        name="Spectrum (points)",
        marker=dict(
            color="deeppink",
            size=10,
        ),
        hovertemplate=(
            "(%{x:.3f}, %{y:.3f})"
            "<extra></extra>"
        ),
    )

    fig.update_layout(
        width=800,
        height=600,
        xaxis_title="Wavenumber [cm⁻¹]",
        yaxis_title="MS Fraction %",
        template="plotly_white",
        hovermode="closest",
        showlegend=False,
    )

    fig.update_xaxes(
        showline=True,
        linewidth=1,
        linecolor="black",
        mirror=False,
        ticks="outside",
        tickwidth=2,
        ticklen=6,
        showgrid=False,
    )

    fig.update_yaxes(
        showline=True,
        linewidth=1,
        linecolor="black",
        mirror=False,
        ticks="outside",
        tickwidth=2,
        ticklen=6,
        showgrid=False,
    )

    fig.update_layout(
        font=dict(
            family="Latin Modern Roman, Times New Roman, serif",
            #size=18,
            #color="black",
        )
    )
    
    return fig


def make_atd_fig_plotly(step_data, ms_cut_pos):
    arrival_us = step_data["cycle_time"] * 1e6

    fig = go.Figure()

    fig.add_histogram(
        x=arrival_us,
        xbins=dict(
            start=100, end=500, size=2.5
        ),
        name="ATD",
        hovertemplate="Arrival time: %{x:.1f} μs<extra></extra>",
    )

    # GS cutoff line
    fig.add_vline(
        x=ms_cut_pos * 1e3,
        line_width=2,
        line_dash="dash",
        line_color="black",
        annotation_text="GS cutoff",
        annotation_position="top",
    )

    fig.update_layout(
        xaxis_title="Arrival Time [μs]",
        yaxis_title="Counts",
        template="plotly_white",
        bargap=0.05,
    )

    return fig