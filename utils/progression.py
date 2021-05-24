import os
import re
from glob import glob
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from .common import unpack_config, generate_init, get_2hg_gln_glu, compile_inputs, print_scan_info, print_delimiter, save_plot

def compile_progression_data(config_dict, axes):
    roi_filenames = {}
    
    for roi_name, dir_list in config_dict.items():
        file_list = []
        
        for directory in dir_list:
            filename = glob(directory + '/*.CSV')[0]
            file_list.append(filename)
            
        roi_filenames[roi_name] = file_list
        
    data = pd.DataFrame(columns=['date', 'roi', 'metabolite', 'value']) # metabolites correspond to axes
    
    for roi_name, file_list in roi_filenames.items(): 
        for filename in file_list:
            metabolites = pd.read_csv(filename).iloc[0] # read scan file
            date = re.search('[0-9]{8}', filename).group(0) # get scan date
            
            for met_i, met_name in enumerate(axes):
                val = get_2hg_gln_glu(metabolites) if met_name == '(2HG+Gln)/Glu' else metabolites[' ' + met_name]
                    
                data = data.append({
                    'date': date,
                    'roi': roi_name,
                    'metabolite': met_name,
                    'value': val
                }, ignore_index=True)
                
    return data

def plot_progression(data, axes, title):
    df_normative = pd.read_excel('normative_ranges.xlsx', sheet_name='SVS97')
    df_normative[' (2HG+Gln)/Glu'] = df_normative.apply(lambda row: get_2hg_gln_glu(row), axis=1)
    
    # sort and filter data
    roi_list = data['roi'].unique()
    data['sort_by_date'] = pd.to_datetime(data['date'], format='%Y%m%d').astype('int64')
    data['drop_by_date'] = pd.to_datetime(data['date']).dt.date
    data['display_date'] = pd.to_datetime(data['date'], format='%Y%m%d').dt.strftime('%m/%d/%Y')
    data.sort_values(by='sort_by_date')

    # plot data
    sns.set()
    plt.close()
    fig, ax = plt.subplots(1, len(axes), figsize=(6.4*len(axes), 4.8))
    fig.tight_layout(pad=4)
    
    for ax_i, ax_name in enumerate(axes):        
        for roi_name in roi_list:
            # Plot metabolite values
            vals = data[(data['metabolite'] == ax_name) & (data['roi'] == roi_name)]
            ax[ax_i].plot(
                'sort_by_date',
                'value',
                data=vals,
                marker='D',
                label=roi_name if ax_i == 0 else ''
            )

        # Plot normal normative range
        met_data = df_normative[df_normative['Tissue Type'] == 'normal'][' ' + ax_name]
        
        upper_quartile = np.percentile(met_data, 75)
        lower_quartile = np.percentile(met_data, 25)
        iqr = upper_quartile - lower_quartile
        mean_normal = np.mean(met_data)
        max_normal = met_data[met_data<=upper_quartile+1.5*iqr].max()
        min_normal = met_data[met_data>=lower_quartile-1.5*iqr].min()
        
        ax[ax_i].axhline(mean_normal, label='Normal Tissue Mean' if ax_i == 0 else '', linestyle='solid')
        ax[ax_i].axhline(max_normal, label='Normal Tissue Min/Max ' if ax_i == 0 else '', linestyle='dashed')
        ax[ax_i].axhline(min_normal, linestyle='dashed')

        date_labels = data[(data['metabolite'] == ax_name)].drop_duplicates(subset='drop_by_date')
        ax[ax_i].set_xticks(date_labels['sort_by_date'])
        ax[ax_i].set_xticklabels(date_labels['display_date'])

        # Set axis title
        ax[ax_i].set_title(ax_name)
        
    # Congfiguure figure
    fig.autofmt_xdate(rotation=75)
    fig.suptitle(title)
    fig.legend()

    return plt

def generate_progression(config, dates_to_exclude=[], scan_history=None):
    scan_id, scan_type, output_dir, verbose = unpack_config(config)
    
    query_string = f'/external/SFA/SPECTRO_PROC/CCS_AutoProc/output/{scan_id}/**/*_{scan_type}*.CSV'
    query_results = generate_init(query_string, scan_id, scan_type)
    
    best_pipelines = list(filter(lambda x: 'AllPipelines' not in x, query_results))

    metabolites = ['2HG', '2HG/Cr+PCr', '(2HG+Gln)/Glu']
    study_input_fn = glob(os.path.join(os.path.expanduser('~'), 'mnt/spectro_proc/CCS_Proc/Study_2HG/SVS-PRESS/Input*'))[0]
    study_output_fn = glob(os.path.join(os.path.expanduser('~'), 'mnt/spectro_proc/CCS_Proc/Study_2HG/SVS-PRESS/Output*'))[0]
    study_input_df = pd.read_excel(study_input_fn)
    study_output_df = pd.read_excel(study_output_fn, sheet_name='BestPipeline')

    if (scan_history):
        mode = 'manual'
        print('Manual Mode')

        data = pd.DataFrame(columns=['specn', 'id', 'date', 'roi', 'metabolite', 'value'])

        for roi, specn_list in scan_history.items():
            for specn in specn_list:
                if (type(specn) == str):
                    if (specn.endswith('.CSV')):
                        regex = specn
                    else:
                        regex = f'{specn}/*.CSV'
                    scan_file = glob(regex)[0]
                    scan_df = pd.read_csv(scan_file)

                    for met_name in metabolites:
                        data = data.append({
                            'date': re.search('[0-9]{8}', specn).group(0),
                            'roi': roi,
                            'metabolite': met_name,
                            'value': float(get_2hg_gln_glu(scan_df)) if met_name == '(2HG+Gln)/Glu' else scan_df.iloc[0][' ' + met_name]
                        }, ignore_index=True)
                else:
                    scan_df = study_output_df[study_output_df['SpecN'] == specn].iloc[0]

                    for met_name in metabolites:
                        data = data.append({
                            'specn': specn,
                            'id': scan_df['ID'],
                            'date': re.search('[0-9]{8}', scan_df['ID']).group(0),
                            'roi': roi,
                            'metabolite': met_name,
                            'value': scan_df['(2HG+Gln)/Glu' if met_name == '(2HG+Gln)/Glu' else ' ' + met_name]
                        }, ignore_index=True)

        data = data.sort_values(by='date')

    else:
        mode = 'auto'
        print('Automatic Mode')

        patient_name, scan_date = scan_id.split('_')
        pt_scans_df = study_input_df[study_input_df['folder'].str.contains(patient_name)]
        
        filenames = pd.Series(best_pipelines).str.split('/').str[-2]
        autoproc_uncon_rois = (
            filenames
                .where(
                    ~filenames.str.endswith(f'{scan_type}_CL')
                    & ~filenames.str.endswith(f'{scan_type}')
                )
                .dropna()
                .unique()
                .tolist()
        )
                
        db_uncon_rois = (
            pt_scans_df[
                ~pt_scans_df['name'].str.endswith(f'_{scan_type}.dat') &
                ~pt_scans_df['name'].str.endswith(f'_{scan_type}_CL.dat')
            ]['name']
                .str.split(f'_{scan_type}_').str[1]
                .str.replace('.dat', '')
                .unique()
                .tolist()
        )

        if (len(pt_scans_df) == 0):
            raise Exception(
                f'{patient_name} has no past scans.\n'
                f'If this is incorrect, double check that the folder name in'\
                f'input/ matches the past scans in the CCS_Database \'folder\' column.'
            )
        elif (len(autoproc_uncon_rois) > 0):
            raise Exception(
                f'{patient_name} has unconventional ROIs in AutoProc: {*autoproc_uncon_rois,}.\n'
                f'Please pass in scan_history to process this patient manually.'
            )
        elif (len(db_uncon_rois) > 0):
            raise Exception(
                f'{patient_name} has unconventional ROIs in CCS_Database: {*db_uncon_rois,}.\n'
                f'Please pass in scan_history to process this patient manually.'
            )
        else:
            data = pd.DataFrame(columns=['date', 'roi', 'metabolite', 'value'])
            
            for index, row in pt_scans_df.iterrows():
                for met_name in metabolites:
                    date = re.search('[0-9]{8}', row['folder']).group(0)
                    if (date not in dates_to_exclude):
                        data = data.append({
                            'date': date,
                            'roi': 'Lesion' if row['name'].endswith('_svs_se_97.dat') else 'Contralateral',
                            'metabolite': met_name,
                            'value': study_output_df[study_output_df['SpecN'] == index].iloc[0][
                                '(2HG+Gln)/Glu' if met_name == '(2HG+Gln)/Glu' else ' ' + met_name
                            ]
                        }, ignore_index=True)

            for file in best_pipelines:
                if (os.path.dirname(file).endswith('_svs_se_97') or os.path.dirname(file).endswith('_svs_se_97_CL')):
                    file_df = pd.read_csv(file)

                    for met_name in metabolites:
                        data = data.append({
                            'date': scan_date,
                            'roi': 'Lesion' if os.path.dirname(file).endswith('_svs_se_97') else 'Contralateral',
                            'metabolite': met_name,
                            'value': (float(get_2hg_gln_glu(file_df))
                                if met_name == '(2HG+Gln)/Glu'
                                else file_df.iloc[0][' ' + met_name])
                        }, ignore_index=True)

    display(data.head(len(data)))

    plot = plot_progression(data, metabolites, scan_id)
    plot = save_plot(plot, f'{output_dir}{scan_id}/{scan_id}_progression_{mode}.png')
    plot.show()

    print(f'Processing complete. Progression plot saved to {output_dir}{scan_id}')
