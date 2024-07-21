import streamlit as st
import pandas as pd
import numpy as np
import time, os

import io_utils
import visual_utils
import line_fitting
import saturation_curve

from bokeh.plotting import figure
from bokeh.models import Span

from pathlib import Path

# set title page and favicon displayed on the web app
st.set_page_config(page_title='LRC Analysis Pipeline', page_icon=':stars:')

# I want to show the LRC logo at the top
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


if 'file_modified' not in st.session_state:
    st.session_state.file_modified = True
# initialize a variable stored in the session state, which is accessible globally, in a safe way.
if 'file_name' not in st.session_state:
    st.session_state.file_name = None
if 'file_bytes_data' not in st.session_state:
    st.session_state.file_bytes_data = None
if 'buncher_delay' not in st.session_state:
    st.session_state.buncher_delay = None
if 'bunching_freq' not in st.session_state:
    st.session_state.bunching_freq = None

def handle_file_change():
    """
    This method is invoked each time the file uploader is touched, or the session
    is refreshed

    NOTE: After uploading a file, it will remain in memory for a while before being wiped out.
       => Make sure that the routine that reads data from this file isn't being called
       over and over, but only when the file is changed; This ensures you don't get "Broken Pipe"
        errors because the file is no longer there.
    """
    print("Checking file change", time.time())
    # make sure the variable at least exists.
    if 'uploaded_file' in globals():
        print("THIS TEST PASSED")
        # check if the file has likely changed
        if uploaded_file.name != st.session_state.file_name:
            print("FILE NAME CHANGED!")
            st.session_state.file_name = uploaded_file.name
        
        # To see if the uploaded file has *actually* changed, it's better to compare the bytes data
        # instead of just the file name as the file with the same name can get modified locally
        if True:#uploaded_file.getvalue() != st.session_state.file_bytes_data:
            # take action based on loaded file type (raw time series or spectrum)
            prepared = True
            if 'raw' in file_type.lower():
                print("CHECKING IF PREPARED".lower())
                if st.session_state.buncher_delay == None or st.session_state.bunching_freq == None:
                    prepared = False
                if prepared == True:
                    # update the global variable that contains the loaded series
                    st.session_state.atd = load_time_series(uploaded_file)
            elif 'spectrum' in file_type.lower():
                # TODO: call the method to load spectrum here.
                pass
            else:
                st.write("ERROR! Did you change the name of the specified file types?\n" + \
                         "The type could not be identified")
            # update the stored file bytes data
            st.session_state.file_bytes_data = uploaded_file.getvalue()


uploaded_file = st.file_uploader("Enter the file you'd like to read from", on_change=handle_file_change())

def plot_spectrum_from_file(file):
    """
    When the file type specified is an extracted "spectrum" which can be directly displayed
    """
    st.write('Reading file:', file.name)
    global spec_data
    spec_data = pd.read_csv(file)
    plot_spectrum_from_data(spec_data)


def plot_spectrum_from_data(spec_data):
    """
    A more generic method for displaying the spectrum, when the 
    """
    global spec_fig
    spec_fig = visual_utils.make_spectrum_fig(spec_data)
    st.bokeh_chart(spec_fig)


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
    gs_cutoff_line = Span(location=ms_cut_pos*1E3, dimension='height', line_color='red',
                          line_width=1)
    fig_atd.renderers.extend([gs_cutoff_line])
    st.bokeh_chart(fig_atd, use_container_width=True)


st.write(st.session_state)

# TODO: This info should be displayed only when loading raw time series.
st.write("""In some measurements you may have kept the laser off,
          in which case the wavemeter read out a garbage value such as -33333""")

garbage_wave_flag = st.radio("""Discard recorded hits with garbage wavelength readouts? (e.g. -33,333)""",
                            options=['Yes', 'No, keep'])


def load_time_series(uploaded_file):
    return io_utils.read_lrc_timeseries(uploaded_file,
                                         st.session_state.buncher_delay,
                                         1/st.session_state.bunching_freq,
                                         discard_garbage=(garbage_wave_flag == 'Yes'))


if 'display_spectrum' not in st.session_state:
    st.session_state.display_spectrum = True

if 'atd' not in st.session_state:
    st.session_state.atd = None
#def toggle_spectrum_display():
#st.session_state.display_spectrum = True#not st.session_state.display_spectrum


if uploaded_file is not None:
    handle_file_change()
    if 'spectrum' in file_type.lower():
        try:
            handle_file_change()
            plot_spectrum_from_file(uploaded_file)
        except Exception as e:
            st.write(f"Error plotting spectrum from file {uploaded_file.name}!\n" + \
                  "Please check if you're using the correct file/option")
            st.write(e)
    else:
        try:
            st.write("### Arrival Time Distribution")
            st.session_state.bunching_freq = st.number_input('Bunching Freq (in Hz)', value=100)
            st.session_state.buncher_delay = st.number_input('Buncher Delay (in seconds)')
            print("CHECKING FILE CHANGE!")
            handle_file_change()
            print("CHECKING IF ATD is DEFINED")
            if st.session_state.atd is not None:
                atd = st.session_state.atd
                print("WAS THIS REACHED?")
                display_wavenum = st.select_slider("Select wavenumber step to display",
                                        options=list(dict.fromkeys(st.session_state.atd['wavenum_req'])))
                st.write("Current setting at", display_wavenum*2, 'second harmonic. Please mind that ' + \
                                    "there may be an offset between requested and actual. " + \
                                        "Trust the WS7 wavemeter readout instead.")
                ms_cut_pos = st.slider('GS Cutoff [ms]', value=0.295, step=1E-3, format='%.3f')
                display_atd(st.session_state.atd, wavenum_req=display_wavenum)
            
            filter1 = 0.#st.number_input('OD 1', value=0.)
            filter2 = 0.#st.number_input('OD 2', value=0.)
            #load_button = st.button("Hide/Show Spectrum!", on_click=toggle_spectrum_display())
            
        except Exception as e:
            st.write(f'Error plotting the ATD from {uploaded_file.name}!\n' +
                  'Please check if all wavenumbers are not garbage (e.g. -33333)')
            st.write(e)
            if st.session_state.atd is not None:
                st.write(st.session_state.atd)
            else:
                st.write("WARNING: Instance of `atd` has disappeared")


if uploaded_file is not None and st.session_state.display_spectrum == True:
    try:
        st.write("Creating spectrum")
        if st.session_state.atd is not None:
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
                            file_name=st.session_state.file_name.replace('.csv', '') + '_spectrum.csv',
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
            st.write(f"""Line {i} has $\Gamma=${float(gamma)*30:.2f} $\pm$ {float(gamma.stderr)*30:.2f} GHz""")
            st.write(f""" and FWHM {float(fwhm*30):.2f} GHz""")
            #st.write(result.fit_report())
            spec_fig.line(x_finer, result.eval(x=x_finer), color='royalblue', line_width=2.0)
            st.bokeh_chart(spec_fig)
    else:
        spec_fig.line(x_finer, fit_results.eval(x=x_finer), color='royalblue', line_width=2.0)
        st.bokeh_chart(spec_fig)
        for v in ['v1_', 'v2_', 'v3_']:
            gamma = fit_results.params[v + 'gamma']
            st.write(gamma, 'Standard Error:',gamma.stderr)
            #fit.params.
            st.write(fit_results.params[v + 'center'], 'Standard Error:',
                     fit_results.params[v+'center'].stderr)
            sigma = fit_results.params[v + 'sigma']
            fwhm = fit_results.params[v + 'fwhm']
            st.write(sigma, 'Standard Error:',sigma.stderr)
            st.write(fit_results.params[v + 'amplitude'], 'Standard Error:', 
                     fit_results.params[v+'amplitude'].stderr)
            st.write(f"""Line {v[1]} has $\Gamma=${float(gamma)*30:.2f} $\pm$ {float(gamma.stderr)*30:.2f} GHz""")
            st.write(f""" and FWHM {float(fwhm*30):.2f} GHz""")
if fit_button:
    # make a finer grid to evaluate the fitted function over (instead of just the observed points)
    x_finer = np.linspace(np.min(spec_data['Wavenumber']), np.max(spec_data['Wavenumber']), 1000)

    if fit_method_option == 'Simultaneously':
        print_fit_results(line_fitting.fit_triple_voigt(spec_data, line1_pos, line2_pos, line3_pos,
                                         gaussian_fwhm/(29.99*2.355),
                                         vary_gauss_component=(gaussian_fixed_flag == 'Variable')))
    else:
        print_fit_results(line_fitting.fit_separate_voigts(spec_data,
                         line1_pos, line2_pos, line3_pos, gaussian_fwhm/(29.99*2.355),
                         vary_gauss_component=(gaussian_fixed_flag == 'Variable')))
    
    #st.write(f"""{fit.fit_report()}""")



"""
______
## Lab Notebook Summary
"""

show_meas_notes = st.radio("Do you want to look at the notes from your past measurements?",
                           options=['Yes, please', 'No'])


notes_file_path = '../notes/measurements/summary_key_measurements.md'

# NOTE: The working directory might be different when deployed on the web.
if 'notes' in os.listdir('./'):
    st.write("Found `notes` folder in the current working directory")
    # get rid of .. and change that to ./
    notes_file_path = notes_file_path[1:]

if 'yes' in show_meas_notes.lower():
    md_file_contents = Path(notes_file_path).read_text()
    st.write(md_file_contents)