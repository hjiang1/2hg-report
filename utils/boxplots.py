import re
import os

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob
from utils.common import unpack_config, generate_init, get_2hg_gln_glu, save_plot, print_delimiter, compile_inputs, print_scan_info

def _plot_against_normative(df_point, df_norm, name):
    # Plot normative ranges
    sns.set()
    g = sns.FacetGrid(df_norm, col="Metabolite", size=7.5, aspect=0.5, sharey=False)
    g.map(sns.boxplot, "Tissue Type", "value", order=["Normal", "Tumor"], hue_order=["Normal", "Tumor"])
    
    thg = round(df_point[' 2HG'], 3)
    thg_crlb = int(df_point[' 2HG %SD'])
    thg_cr = round(df_point[' 2HG/Cr+PCr'], 3)
    thg_gln_glu = round(get_2hg_gln_glu(df_point), 3)
    
    fig, ax = g.fig, g.axes
        
    # Add suptitle and ylabels
    fig.subplots_adjust(top=0.9, right=0.93, wspace=0.4)
    fig.suptitle(name, fontsize=14)
    
    ax[0, 0].set_title(f'2HG (CRLB = {thg_crlb}%)')
    ax[0, 0].set_ylabel('Concentration (mmol/Kg water)')
    ax[0, 1].set_title(f'2HG/Cr+PCr')
    ax[0, 1].set_ylabel('Ratio')
    ax[0, 2].set_title(f'(2HG+Gln)/Glu')
    ax[0, 2].set_ylabel('Ratio')
    
    # Plot pt reference
    dot_params = {
        'marker': 'D',
        'markersize': 8,
        'markerfacecolor': 'red',
        'markeredgewidth': 1,
        'markeredgecolor': 'white',
        'linewidth': 0,
        'zorder': 4
    }
    line_params = {
        'linestyle': 'dashed',
        'color': 'red',
        'zorder': 3
    }
    text_params = {
        'color': 'red',
        'verticalalignment': 'center'
    }
    
    ax[0, 0].plot([0, 1], [thg]*2, **dot_params)
    ax[0, 0].axhline(float(thg), **line_params)
    ax[0, 0].text(1.55, thg, float(thg), **text_params)

    ax[0, 1].plot([0, 1], [thg_cr]*2, **dot_params)
    ax[0, 1].axhline(float(thg_cr), **line_params)
    ax[0, 1].text(1.55, thg_cr, float(thg_cr), **text_params)

    ax[0, 2].plot([0, 1], [thg_gln_glu]*2, **dot_params)
    ax[0, 2].axhline(float(thg_gln_glu), **line_params)
    ax[0, 2].text(1.55, thg_gln_glu, float(thg_gln_glu), **text_params)

    return plt

def _load_normative_data():
    """Read and format normative ranges"""

    df_normative = pd.read_excel('normative_ranges.xlsx', sheet_name='SVS97')

    df_normative['Tissue Type'] = df_normative['Tissue Type'].map({'normal': 'Normal', 'tumor': 'Tumor'}) # Capitalize tissue types
    df_normative_filtered = df_normative.drop( # Drop unwanted columns in df
        df_normative.columns.difference(['Tissue Type', ' 2HG', ' 2HG/Cr+PCr', '(2HG+Gln)/Glu']),
        1
    )
    df_normative_long = pd.melt( # Melt into long format for seaborn
        df_normative_filtered,
        id_vars=['Tissue Type'],
        value_vars=[' 2HG', ' 2HG/Cr+PCr', '(2HG+Gln)/Glu'],
        var_name='Metabolite'
    )
    
    return (
        df_normative,
        df_normative_filtered,
        df_normative_long
    )
    
def _generate_plot(file, scan_id, scan_type, output_dir, verbose):
    df_normative, df_normative_filtered, df_normative_long = _load_normative_data()

    df_input = pd.read_csv(file)
    metabolites = df_input.iloc[0]

    # extract metadata
    pt_name = re.search('[a-zA-Z0-9_]*_[0-9]{8}', file).group(0)
    pipeline = file.split('-')[1].split('.')[0]
    roi = (re.search(f'_{scan_type}.*\.CSV', file.split('/')[-1])
        .group(0)
        .split('-')[0]
        .replace(f'_{scan_type}', '')
        .replace('_', ' ')
        .replace('CL', 'Contralateral')
        .capitalize()
    )
    roi = (' Lesion' if roi == '' else roi)[1:]

    # generate plot
    plot_title = f'{pt_name}; ROI = {roi.title()}; Pipeline = {pipeline}'
    plot = _plot_against_normative(df_input, df_normative_long, plot_title)
    output_fn = os.path.join(
        output_dir,
        scan_id,
        file.split(scan_id + '/')[1].replace('input', 'output')
    ).replace('.CSV', '.png')
    plot = save_plot(plot, output_fn)

    # display values
    if ((verbose > 0) and (verbose > 1 or 'AllPipelines' not in file)):        
        if ('AllPipelines' not in file):
            print('Best Pipeline')

        plot.show()
        
        print_delimiter()

    plot.close()
    
def generate_plots(config):
    scan_id, scan_type, output_dir, verbose = unpack_config(config)
    
    query_string = f'/external/SFA/SPECTRO_PROC/CCS_AutoProc/output/{scan_id}/**/*_{scan_type}*.CSV'
    query_results = generate_init(query_string, scan_id, scan_type)

    for file in query_results:
        try:
            _generate_plot(file, scan_id, scan_type, output_dir, verbose)
        except:
            print(f'File: {file}\n Processing failed\n')
            print_delimiter()
            
    print(f'Processing complete. Plots saved to {output_dir}{scan_id}')
    

