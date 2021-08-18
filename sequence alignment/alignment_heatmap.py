### Reading a single blast output file
### And getting information on the blast output
import sys
import numpy as np
import pandas as pd
import statistics
import glob
import os

import matplotlib.pyplot as plt
import seaborn as sns

def get_filepath():
    try:
        sys.argv[1]
    except:
        sys.exit('\nERROR: No file given\n')
    return(sys.argv[1])

# Load BLAST output file
def concatinate_files(file_path, output):
    #create format final DataFrame
    data = pd.DataFrame(columns = ['perc', 'gen_marker', 'true_case', 'num_rel', 'ref_seq'])


    all_files = glob.glob(file_path)
    #print(' all files')
    #print(all_files)
    for file in all_files:
        print(file)
        #Load file and split in lines
        try:
            blast_output = open(file,'r')
        except:
            sys.exit("\nERROR: Something went wrong when trying to load the file. Is the filename / path correct?\n")

        file_data = pd.read_csv(file, sep = ',')
        file_data = file_data[0:1]

        #add file_data to data(frame) per column of file_data
        #for col in file_data
        for key, value in file_data.iteritems():
            if len(key.split('_')) > 3:
                #print('key', key, '\n', 'value:', type(value), value.values[0])

                perc = value.values[0] # / 100
                true_case = key.split('_')[-1].replace(')',  "")
                gen_marker = key.split('_')[1].split(' ')[0] + ' ' + true_case
                num_rel = "+" + key.split('_')[-2]
                ref_seq = key.split(' ')[0]
                temp_data = pd.DataFrame([[perc, gen_marker, true_case, num_rel, ref_seq]], columns = ['perc', 'gen_marker', 'true_case', 'num_rel', 'ref_seq'])

                data = data.append(temp_data)

    return data


def results_plot(data, output):

    print(data.head())

    temp_data = data
    temp_data = temp_data.groupby(['gen_marker', 'num_rel'])['perc'].mean().reset_index()
    print(temp_data)
    temp_data = temp_data.pivot('gen_marker', 'num_rel', 'perc')
    temp_data = temp_data.reindex([ 'EF1 N', 'RPB2 N','RPB2 P','EF1 P', ])
    print("hoi")
    print(temp_data)

    sns.set(font_scale = 1.1)
    p = sns.heatmap(temp_data, annot = True, cmap="YlGnBu", cbar=False, fmt= '.2f')
    p.set(xlabel = 'relatives added to input data', ylabel = 'Genetic marker reference and consensus \n created with positive or negative beetle input data')
    p.set_title('avg. % identity')
    #plt.show()
    plt.savefig(output + '_heatmap', dpi=125)




if __name__ == '__main__':
    import argparse
    from argparse import RawDescriptionHelpFormatter

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    ap.add_argument('-i', '--input', type=str, help='input file(s)', metavar='input file', nargs='+')

    io_group = ap.add_mutually_exclusive_group(required=True)
    io_group.add_argument('-o', '--output', type=str, help='PDF output file')

    cmd = ap.parse_args()

    print(cmd.input, cmd.output)
    #filepath = get_filepath()

    #important parameters
    #pathogen = 'Fusarium circinatum\n'
    #write columsn of results table
    #f = open(cmd.output, 'w')
    #f.write('consensus;hits;% pathogen;avg evalue pathogen;avg coverage pathogen;gen marker;experiment\n')
    #f.close()

    table = concatinate_files(cmd.input[0], cmd.output)
    results_plot(table, cmd.output)
    #line_plot(table, cmd.output)

    #analyse results
    #results = analyse(cmd.output)

    #list_gen_markers_groups = get_gen_markers_groups(results)
    #Plot results table
    #for group in list_gen_markers_groups:
    #    table = results_plot(results, cmd.output, group)

    #print(table)
