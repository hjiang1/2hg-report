import os
import matplotlib.pyplot as plt
import pandas as pd
from glob import glob

def parse_filename(file, scan_type):
    roi, pipeline = file.split(f'{scan_type}')[-1].split('.')[0].split('-')
    
    if (roi == ''):
        roi = 'Lesion'
    else:
        roi = roi[1:]
        
    return (roi, pipeline)

def compile_inputs(query_string):
    query_results = glob(query_string, recursive=True)

    if (len(query_results) == 0):
        raise Exception(f'ERROR: No results found matching the following glob query: {query_string}')
    else:
        query_results.sort()
        
        return query_results
    
def print_scan_info(scan_id, query_results, scan_type):
    df = pd.DataFrame(columns=['ROI', 'Best Pipeline', 'File'])
    best_pipelines = list(filter(lambda x: 'AllPipelines' not in x, query_results))
    best_pipelines.sort()
    
    for file in best_pipelines:
        roi, pipeline = parse_filename(file, scan_type)
        
        df = df.append({
            'ROI': roi,
            'Best Pipeline': pipeline,
            'File': file
        }, ignore_index=True)
        
    print(f'Scan: {scan_id}')
    pd.set_option('display.max_colwidth', -1)
    display(df)
    
def generate_init(query_string, scan_id, scan_type):
    query_results = compile_inputs(query_string)
    print_scan_info(scan_id, query_results, scan_type)
    print_delimiter(char='=')
    
    return query_results

def get_2hg_gln_glu(ser, sigfig=None):
    if (float(ser[' Glu']) == 0):
        return -1
    
    val = (ser[' 2HG']+ser[' Gln'])/ser[' Glu']
    
    if (sigfig):
        return round(val, sigfig)
    else:
        return val

def save_plot(plot, file):
    directory = os.path.dirname(file)
    
    if (not os.path.exists(directory)):
        os.makedirs(directory)
        
    plot.savefig(file)
    
    return plot

def unpack_config(config):
    return tuple(config.values())

def print_delimiter(char='-', num=50, end='\n'):
    print(char * num, end=end)
