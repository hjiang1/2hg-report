import os

import pandas as pd
from utils.common import unpack_config, generate_init, parse_filename, get_2hg_gln_glu, print_delimiter, compile_inputs, print_scan_info

def generate_table(config, pipelines={}):
    scan_id, scan_type, output_dir, verbose = unpack_config(config)
    
    query_string = f'/external/SFA/SPECTRO_PROC/CCS_AutoProc/output/{scan_id}/**/*_{scan_type}*.CSV'
    query_results = generate_init(query_string, scan_id, scan_type)

    df = pd.DataFrame(columns=['ROI (Pipeline)', '2HG (CRLB)', '2HG/Cr+PCr', '(2HG+Gln)/Glu'])
    best_pipelines = list(filter(lambda x: 'AllPipelines' not in x, query_results))
    all_pipelines = list(filter(lambda x: 'AllPipelines' in x, query_results))
    best_pipelines.sort()

    rois = {}

    for file in best_pipelines:
        roi, pipeline = parse_filename(file, scan_type)
        rois[roi] = {
            'best_pipeline': pipeline
        }
        
    for file in all_pipelines:
        roi, pipeline = parse_filename(file, scan_type)
        rois[roi][pipeline] = file

    for roi, files in rois.items():
        if (roi in pipelines):
            pipeline = pipelines[roi]
        else:
            pipeline = files['best_pipeline']
            
        is_best = pipeline == files['best_pipeline']

        file = files[pipeline]
        file_df = pd.read_csv(file)
        file_row = file_df.iloc[0]
        df = df.append(
            {
                'ROI (Pipeline)': f'{roi} ({pipeline}*)' if is_best else f'{roi} ({pipeline})',
                '2HG (CRLB)': '{} ({}%)'.format(round(file_row[' 2HG'], 3), int(file_row[' 2HG %SD'])),
                '2HG/Cr+PCr': str(round(file_row[' 2HG/Cr+PCr'], 3)),
                '(2HG+Gln)/Glu': str(round(get_2hg_gln_glu(file_row), 3))
            },
            ignore_index=True
        )
    
    display(df)
    
    if (not os.path.exists(f'{output_dir}{scan_id}')):
        os.makedirs(f'{output_dir}{scan_id}')
    
    df.to_csv(f'{output_dir}{scan_id}/{scan_id}.csv', index=False)
          
    print(f'Processing complete. Table saved to {output_dir}{scan_id}')