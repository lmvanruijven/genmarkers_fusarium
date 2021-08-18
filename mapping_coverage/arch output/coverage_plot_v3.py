"""
plot coverage
"""

from __future__ import print_function, division

import csv
import os
import re
import shlex
import sys
import glob
import os

try:
    import matplotlib
    import itertools
    #import cbook as cbook
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    import pandas as pd
except ImportError as e:
    print('[!] The required Python libraries could not be imported:', file=sys.stderr)
    print('\t{0}'.format(e))
    sys.exit(1)

def parse_list(fnames, sel_columns = 'all'):
    for fname in fnames:
        data = parse_item(fname, sel_columns=sel_columns)

    return data

def parse_item(fname, sel_columns='all'):
    """Remove all files from data folder"""
    loc = glob.glob('data/*')
    for f in loc:
        print(f)
        os.remove(f)

    """read xlsx file and convert to dataframe"""

    excel_file = fname
    all_sheets = pd.read_excel(excel_file, sheet_name = None)
    sheets = all_sheets.keys()

    for sheet_name in sheets:
        sheet = pd.read_excel(excel_file, sheet_name = sheet_name)
        sheet['Experiment'] = sheet_name
        sheet.to_csv("data/%s.csv" % sheet_name, index = False)

    all_files = glob.glob(os.path.join("data", "*.csv"))
    df_from_each_file = (pd.read_csv(f, sep = ',', ) for f in all_files)
    df_merged = pd.concat(df_from_each_file, ignore_index = True)
    #Add Gen marker en true case columns
    df_merged['Gen marker'] = df_merged['Reference sequence'].str.split('_').str[-1]
    df_merged['True case'] = df_merged['Experiment'].str.split('_').str[0]
    #save file when needed
    df_merged.to_csv("data/merged.csv")
    #

    list_agg = [ 'mean', 'min', 'max', 'std']
    data = df_merged.groupby( ['Gen marker','Experiment', 'True case'], as_index = False).agg({'Average coverage': list_agg})
    data.columns = data.columns.droplevel(level =0)
    data.columns = ['Gen marker', 'Experiment', 'True case', 'mean', 'min', 'max', 'std']


    return data

def plot_data_v1(gen_marker_group, dist_exp, data, outfile):
    f = plt.figure()
    ax = plt.gca()

    #gen_marker_group = ['ITS', '28S']
    print(gen_marker_group)

    #create marker figures
    markerP = itertools.cycle((',', '+', '.', 'o', '*',  'd', 'p', 'h', 'v', 'x'))
    markerN = itertools.cycle((',', '+', '.', 'o', '*',  'd', 'p', 'h', 'v', 'x'))

    #make subplot
    for exp in dist_exp:
        msft = data[(data['Experiment'] == exp)]
        msft = msft[(msft['Gen marker'].isin(gen_marker_group))]
        #print(msft)
        if msft['True case'].iloc[0] == 'P':
            color = 'g'
        else:
            color = 'r'
        #plot per experiment
        #sns.stripplot(msft['Gen marker'], msft['mean'], color = color, label = exp,
        #                marker = next(markerN), s = 11, jitter = 0.5)
        ax.scatter(msft['Gen marker'], msft['mean'], color = color, label = exp,
                        marker = next(markerN), s = 85 )


    # Formatting Labels & Appearance
    ax.set_xlabel("Gen marker")
    ax.set_ylabel("Average coverage")
    #ax.set_yscale('log')
    ax.grid('on')

    try:
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.92, box.height])
        #put legegd on the right side
        legend = ax.legend(loc = 'center left', bbox_to_anchor = (1 , 0.5))
        frame = legend.get_frame()
        #frame.set_facecolor(bg_color)
    except AttributeError as e:
        # No legend, likely because no labels
        pass

    if outfile:
        plt.savefig(outfile + str(gen_marker_group), dpi=300)

    plt.show()


def plot_coverage_version1(data, outfile, list_gen_markers_groups):
    #f, axes = plt.subplots(nrows = 1, ncols = 2)


    dist_exp = data['Experiment'].unique()

    #genereting plots
    for group in list_gen_markers_groups:
        plot_data_v1(group, dist_exp, data, outfile)


def plot_coverage_version2(data, outfile):
    print(data)
    n_series = len(data) - 1

    f = plt.figure()
    ax = plt.gca()

    #set list of colours for genetic markers
    #color_map = getattr(plt.cm, colormap)
    #dist_gen_marker = data['Gen marker'].unique()
    #print(dist_gen_marker)
    #color_list = color_map(np.linspace(0, 1,dist_gen_marker))
    #print(color_list)

    TC_list = ['N', 'P']
    markerP = itertools.cycle((',', '+', '.', 'o', '*',  'd', 'p', 'h', 'v', 'x'))
    markerN = itertools.cycle((',', '+', '.', 'o', '*',  'd', 'p', 'h', 'v', 'x'))
    color = itertools.cycle(('b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'))

    #make subselection of each Experiment
    dist_gen = data['Gen marker'].unique()

    for gen in dist_gen:
        msft = data[(data['Gen marker'] == 'ITS')]#gen)]
        #print(msft)
        if msft['True case'].iloc[0] == 'P':
            ax.scatter(msft['Experiment'], msft['mean'], color = next(color), label = msft['Gen marker'],
            marker = next(markerP))
        else:
            ax.scatter(msft['Experiment'], msft['mean'], color = next(color), label = msft['Gen marker'],
            marker = next(markerN))

    # Formatting Labels & Appearance
    ax.set_xlabel("Gen marker")
    ax.set_ylabel("Average coverage")
    ax.set_yscale('log')
    ax.grid('on')

    try:
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.92, box.height])
        #put legegd on the right side
        legend = ax.legend(loc = 'center left', bbox_to_anchor = (1 , 0.5))
        frame = legend.get_frame()
        #frame.set_facecolor(bg_color)
    except AttributeError as e:
        # No legend, likely because no labels
        pass

    if outfile:
        plt.savefig(outfile, dpi=300)

    plt.show()


def get_gen_markers_groups(data):
    gen_markers_group_1 = np.array(['18S', '28S', '5S', 'ITS'])
    index_group_1 = np.where(np.isin(data['Gen marker'].unique(), gen_markers_group_1))
    gen_markers_group_2 = np.delete(data['Gen marker'].unique(), index_group_1)

    gen_marker_groups_list = [gen_markers_group_1 , gen_markers_group_2 ]

    return gen_marker_groups_list


if __name__ == '__main__':

    import argparse
    from argparse import RawDescriptionHelpFormatter

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    ap.add_argument('-xls', type=str, help='input file(s)', metavar='input file', nargs='+')

    io_group = ap.add_mutually_exclusive_group(required=True)
    io_group.add_argument('-o', '--output', type=str, help='PDF output file')

    cmd = ap.parse_args()

    print(cmd.xls, cmd.output)
    metadata = parse_list(cmd.xls)

    list_gen_markers_groups = get_gen_markers_groups(metadata)

    plot_coverage_version1(metadata, cmd.output, list_gen_markers_groups)
    #plot_coverage_version2(metadata, cmd.output)


#    print("Hoi, parsed xvgs", "\n\n", data)
    #n_series = len(data[1:])
    #n_elements = sum(list(map(len, data[1:])))
    #print(n_elements)
    #print('[+] Read {0} series of data ({1} elements)'.format(n_series, n_elements))



#    plot_data(data, metadata, cmd.xvg_f,
#              window=cmd.window,
#              interactive=cmd.interactive, outfile=cmd.output,
#              colormap=cmd.colormap, bg_color=cmd.background_color)
