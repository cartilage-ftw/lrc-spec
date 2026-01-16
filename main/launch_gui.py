import streamlit as st
import pandas as pd
import numpy as np
import time, os
import scipy.stats

import io_utils
import visual_utils
import line_fitting
import saturation_curve

from streamlit_bokeh import streamlit_bokeh
from bokeh.plotting import figure
from bokeh.models import Span
from pathlib import Path

st.set_page_config(page_title='LRC Analysis Pipeline', page_icon=':stars:')

# display logo
lrc_logo_url = 'https://www.lrc-project.eu/' + \
                r'Webpage%20of%20the%20European%20LRC%20project-Dateien/LRC-color-300x77.png'

st.write('![LRC Logo]({0})'.format(lrc_logo_url))

# Disclaimer
st.write('# LRC Data Analysis Interface')
st.write("""If anything breaks, blame Aayush. As simple as that""")


#spec = pd.read_csv('../data/spectra/2023-11-24-17-48-38_stacked_spec.csv')

"""
Start of GUI buttons, etc.
"""

file_type = st.radio("Content of File:", ['Raw Time Series', 'Extracted Spectrum'])

uploaded_file = st.file_uploader("Enter the file you'd like to read from")


def plot_spectrum_from_file(file):
    st.write('Reading file:', file.name)
    global spec_data
    spec_data = pd.read_csv(file)
    plot_spectrum_from_data(spec_data, file.name)


def plot_spectrum_from_data(spec_data):
    global spec_fig
    spec_fig = visual_utils.make_spectrum_fig(spec_data)
    st.write("Rendering spectrum plot.")
    streamlit_bokeh(spec_fig, key='spectrum_plot')

def g(x, A, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2*sigma**2))

def display_atd(atd, wavenum_req, wavenum_obs='#TODO'):
    step_data = atd[atd['wavenum_req'] == wavenum_req]
    median_wavenum = np.median(step_data['wavenum_obs'])
    weights, bin_edges = np.histogram(step_data['cycle_time']*1E6, bins='fd')

    global fig_atd
    fig_atd = figure(x_axis_label='Arrival Time [μs]',
                     y_axis_label='Frequency',
                     x_range=(100, 700),
                     y_range=(0, 1.05*np.max(weights)))
    fig_atd.quad(top=weights, bottom=0, left=bin_edges[:-1], right=bin_edges[1:],
                 fill_color='skyblue', line_color='white', legend_label=f'{median_wavenum:.2f} cm-1')
    gs_cutoff_line = Span(location=ms_cut_pos*1E3, dimension='height', line_color='black',
                          line_width=1)
     
    # for interactive fit/visualization
    # wanna_fit = st.radio("Display Manual Gaussian Fitter", options=['Yes', 'No'], index=1)
    # if wanna_fit == 'Yes':
    #     a1 = st.slider("G1 Amplitude", value=np.max(weights), max_value=5*np.max(weights))
    #     mu1 = st.slider("G1 Center", value=scipy.stats.mode(bin_edges)[0], max_value=np.max(bin_edges))
    #     s1 = st.slider("G1 Sigma", value=np.std(weights), max_value = 4*np.std(weights))

    #     a2 = st.slider("G2 Amplitude", value=0.1, max_value=5.*np.max(weights))
    #     mu2 = st.slider("G2 Center", value=0.1, max_value=np.max(bin_edges))
    #     s2 = st.slider("G2 Sigma", value=0.1, max_value = 4.*np.std(weights))
    #     x_grid = np.linspace(100, 700, 1000)
    #     y1_grid = g(x_grid, a1, mu1, s1)
    #     y2_grid = g(x_grid, a2, mu2, s2)
    #     y_grid = y1_grid + y2_grid
    #     fig_atd.line(x_grid, y1_grid, line_color='black')
    #     fig_atd.line(x_grid, y2_grid, line_color='black')
    #     fig_atd.line(x_grid, y_grid, line_color='green')

    fig_atd.renderers.extend([gs_cutoff_line])
    st.write("Rendering ATD plot")
    streamlit_bokeh(fig_atd, key='atd_plot', use_container_width=True)



st.write("""In some measurements you may have kept the laser off,
          in which case the wavemeter read out a garbage value such as -33333""")

garbage_wave_flag = st.radio("""Discard recorded hits with garbage wavelength readouts? (e.g. -33,333)""",
                            options=['Yes', 'No, keep'])


def load_atd():
    st.session_state.file_name=uploaded_file.name
    return io_utils.read_lrc_timeseries(uploaded_file, buncher_delay, 1/bunching_freq,
                         discard_garbage=(garbage_wave_flag == 'Yes'))


if 'display_spectrum' not in st.session_state:
    st.session_state.display_spectrum = True


#def toggle_spectrum_display():
#st.session_state.display_spectrum = True#not st.session_state.display_spectrum



if uploaded_file is not None:
    if 'spectrum' in file_type.lower():
        try:
            plot_spectrum_from_file(uploaded_file)
        except Exception as e:
            st.write(f"Error plotting spectrum from file {uploaded_file.name}! " + \
                  "\nPlease check if you're using the correct file/option")
            st.write(e)
    else:
        try:
            st.write("### Arrival Time Distribution")
            bunching_freq = st.number_input('Bunching Freq (in Hz)', value=100)
            buncher_delay = st.number_input('Buncher Delay (in seconds)')

            global atd
            atd = load_atd()
            display_wavenum = st.select_slider("Select wavenumber step to display",
                                        options=list(dict.fromkeys(atd['wavenum_req'])))
            st.write("Current setting at", display_wavenum*2, 'second harmonic. Please mind that ' + \
                                "there may be an offset between requested and actual. " + \
                                    "Trust the WS7 wavemeter readout instead.")
            ms_cut_pos = st.slider('GS Cutoff [ms]', value=0.295, step=1E-3, format='%.3f')
            display_atd(atd, wavenum_req=display_wavenum)
            
            filter1 = 0.#st.number_input('OD 1', value=0.)
            filter2 = 0.#st.number_input('OD 2', value=0.)
            #load_button = st.button("Hide/Show Spectrum!", on_click=toggle_spectrum_display())
            
        except Exception as e:
            st.write(f'Error plotting the ATD from {uploaded_file.name}!\n' +
                  'Please check if all wavenumbers are not garbage (e.g. -33333)')
            st.write(atd)
            st.write(e)


if uploaded_file is not None and st.session_state.display_spectrum == True:
    try:
        st.write("Creating spectrum")
        spectrum_data = io_utils.make_spectrum(atd, ms_cut=ms_cut_pos, filters=[filter1, filter2],
                            transm_percent=saturation_curve.get_transm(filter1, filter2),
                            save_file=False)
        st.write("NOTE: Currently I had to disable calculation of bootstrap uncertainty for the y-axis.\n" + \
                 " It was taking too long!")
        plot_spectrum_from_data(spectrum_data)
        #st.write(st.session_state)
        # also add a button to manually save this spectrum's data
        st.download_button(label='Save Spectrum to Device',
                           data=spectrum_data.to_csv(),
                           file_name=st.session_state.file_name[:-4] + '_spectrum.csv',
                           )
    except Exception as e:
        st.write('Something went wrong while plotting spectrum!\n', e)


st.write("## Line Fitting")

line1_pos = st.number_input('Enter mean value of first multiplet')
line2_pos = st.number_input('Second multiplet')
line3_pos = st.number_input('Third Multiplet')

st.write(""" The fitted profile will be a [Voigt](https://en.wikipedia.org/wiki/Voigt_profile). You can choose if you'd
          prefer the Gaussian (inhomogeneous) component to be fixed (e.g. if the laser bandwidth is known via measurement)
          or varied to obtain the best fit""")

gaussian_fixed_flag = st.radio("Gaussian component", ['Fixed', 'Variable'])

gaussian_fwhm = st.number_input(r'Enter $\Delta \nu_{\textrm{FWHM}}$ in GHz of the Gaussian component of the line profile',
                                value=4.6, disabled=gaussian_fixed_flag =='Variable')


fit_method_list = ['Simultaneously', 'Each line separately']
fit_method_option = st.radio("Fit all three lines", fit_method_list)

fit_button = st.button('Fit Lines!')


def print_fit_results(fit_results):
    if type(fit_results) == list: # if the three were fitted separately
        for i, result in zip(range(1,4), fit_results):
            gamma = result.params[f'v{i}_' + 'gamma']
            st.write(gamma, 'Standard Error:',gamma.stderr)
            #fit.params.
            st.write(result.params[f'v{i}_' + 'center'], 'Standard Error:',
                        result.params[f'v{i}_'+'center'].stderr)
            sigma = result.params[f'v{i}_' + 'sigma']
            fwhm = result.params[f'v{i}_fwhm']
            st.write(sigma, 'Standard Error:',sigma.stderr)
            st.write(result.params[f'v{i}_' + 'amplitude'], 'Standard Error:',
                      result.params[f'v{i}_'+'amplitude'].stderr)
            st.write(f"""Line {i} has $\Gamma=${float(gamma)*30:.2f} $\pm$ GHz""")
            st.write(f""" and FWHM {float(fwhm*30):.2f} GHz""")
            #st.write(result.fit_report())
            spec_fig.line(x_finer, result.eval(x=x_finer), color='royalblue', line_width=2.0)
            st.write("Rendering plot for fitted results display")
            streamlit_bokeh(spec_fig)
    else:
        spec_fig.line(x_finer, fit_results.eval(x=x_finer), color='royalblue', line_width=2.0)
        st.write("also rendering plot for fitted results display")
        streamlit_bokeh(spec_fig)
        for v in ['v1_', 'v2_', 'v3_']:
            gamma = fit_results.params[v + 'gamma']
            st.write(gamma, 'Standard Error:',gamma.stderr)
            #fit.params.
            st.write(fit_results.params[v + 'center'], 'Standard Error:',
                     fit_results.params[v+'center'].stderr)
            sigma = fit_results.params[v + 'sigma']
            fwhm = fit_results.params[v + 'fwhm']
            st.write(sigma, 'Standard Error:',sigma.stderr)
            st.write(fit_results.params[v + 'height'], 'Standard Error:', 
                     fit_results.params[v+'height'].stderr)
            st.write(f"""Line {v[1]} has $\Gamma=${float(gamma)*30:.2f} $\pm$ {float(gamma.stderr)*30:.2f} GHz""")
            st.write(f"Height, {fit_results.params[v+'height'].value:.2f} $\pm$ {fit_results.params[v+'height'].stderr:.2f}")
            st.write(f"""Gaussian $\sigma={float(sigma)*29.9:.2f}$ GHz (Gaussian FWHM: {2.355*sigma*30:.2f}) and total FWHM {float(fwhm*30):.2f} $\pm$ {float(fwhm.stderr*30):.2f} GHz""")
if fit_button:
    # make a finer grid to evaluate the fitted function over (instead of just the observed points)
    x_finer = np.linspace(np.min(spectrum_data['Wavenumber']), np.max(spectrum_data['Wavenumber']), 1000)

    if fit_method_option == 'Simultaneously':
        print_fit_results(line_fitting.fit_triple_voigt(spectrum_data, line1_pos, line2_pos, line3_pos,
                                         gaussian_fwhm/(29.99*2.355),
                                         vary_gauss_component=(gaussian_fixed_flag == 'Variable')))
    else:
        print_fit_results(line_fitting.fit_separate_voigts(spectrum_data,
                         line1_pos, line2_pos, line3_pos, gaussian_fwhm/(29.99*2.355),
                         vary_gauss_component=(gaussian_fixed_flag == 'Variable')))
    
    #st.write(f"""{fit.fit_report()}""")

"""
______
## Lab Notebook Summary
"""

show_meas_notes = st.radio("Do you want to look at the notes from your past measurements?",
                           options=['No', 'Yes, please'])


notes_file_path = '../notes/measurements/summary_key_measurements.md'

# NOTE: The working directory might be different when deployed on the web.
if 'notes' in os.listdir('./'):
    # get rid of .. and change that to ./
    notes_file_path = notes_file_path[1:]

if 'yes' in show_meas_notes.lower():
    md_file_contents = Path(notes_file_path).read_text()
    st.write(md_file_contents)