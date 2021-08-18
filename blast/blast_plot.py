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
def concatinate_files(file_path, output, pathogen):
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
        blast_output = blast_output.readlines()

        blast_hits_dict = format_blast_output(blast_output)
        write_output(blast_hits_dict, output, pathogen, file)


def format_blast_output(blast_hits):
    blast_hits_dict = {}
    #print(blast_hits)
    while blast_hits:
        try:
            BLASTstart = blast_hits.index('# BLASTN 2.9.0+\n')
            BLASTend = blast_hits.index('# BLASTN 2.9.0+\n', BLASTstart+1)
            #print(" BLASTstart:", BLASTstart, " blastend:", BLASTend)

            if len(blast_hits[BLASTstart+3]) < 16:
                #print(len(blast_hits[BLASTstart+3]))
                #print("there are not hits")
                temp_key = blast_hits[BLASTstart+1].replace('# Query: ', "").strip('\n')
                blast_hits_dict[temp_key] = "no hits"
            else:
                #assign column names
                df_hits_columns = blast_hits[BLASTstart +3].replace('# Fields: ', '').replace('\n', "").split(", ")
                #print(df_hits_columns)

                #assign row names
                df_hits_rows = []
                for row in blast_hits[BLASTstart + 5: BLASTend]:
                    df_hits_rows.append(row.split("\t"))
                df_hits = pd.DataFrame(df_hits_rows, columns = df_hits_columns)
                #remove \n from query id value
                temp_key = blast_hits[BLASTstart+1].replace('# Query: ', "").strip('\n')
                blast_hits_dict[temp_key] = df_hits
                #print('temp key', temp_key, "df hits", df_hits['% identity'])

        except ValueError:
            print("BlastN 2.6.0 is not in list any more")
            del blast_hits[0]
        del blast_hits[BLASTstart:BLASTend]

    #print(blast_hits_dict)
    return (blast_hits_dict)

def write_output(blast_hits_dict, output, pathogen, file):
    r = open(output, 'a')
    for key, value in blast_hits_dict.items():
        if  isinstance(value, str): #there are no hits
            frac_pathogen = avg_e_val = avg_coverage = 'no hits'
            hits = 0
        else: #there are hits
            frac_pathogen, avg_e_val, avg_coverage, hits = do_analysis(value, pathogen)
            gen_marker = key.split('_')[10:]
            gen_marker = gen_marker[:-1]
            gen_marker = "_".join(gen_marker)
        rl = [key, ";", hits, ";",  frac_pathogen,";", avg_e_val,";", avg_coverage , ";", gen_marker, ";", file.replace('.txt', "") ]
        #convert all values to string
        rl = [str(x) for x in rl]
        #write values to
        for i in rl:
            r.write(i)
        r.write( "\n")

def do_analysis(value, pathogen):
    #convert blast_hits to dataframe
    data = pd.DataFrame(value)
    #select top 5 hits
    data = data.head(5)

    if data.shape[0] != 0:
        #count pathogen in top 5
        count_pathogen = list(data['subject sci names']).count(pathogen)
        if count_pathogen > 0 :
            #determine total number of hits
            hits = len(data.index)
            #get fraction pathogen of total number of hits
            frac_pathogen = (count_pathogen / hits ) * 100

            #select pathogen hits
            data = data.loc[data['subject sci names'] == pathogen]
            #calculate average coverge and evalue
            avg_e_val = statistics.mean(list(data['evalue'].apply(lambda x: float(x))))
            avg_coverage = statistics.mean(list(data['% identity'].apply(lambda x: float(x))))

        else:
            frac_pathogen = 0
            avg_e_val = 0
            avg_coverage = 0
            hits = len(data.index)
    else:
        frac_pathogen = 'not hits'
        avg_e_val = 'no hits'
        avg_coverage = 'no hits'
        hits = len(data.index)

    return frac_pathogen, avg_e_val, avg_coverage, hits

def analyse(input):
    data = pd.read_csv(input, sep = ';')
    print(data.columns)

    #voor testen
    #data = data.loc[data['experiment'] == 'N_6']

    #add columns
    data['input data'] = data['experiment'].str.split('_').str[1] + data['experiment'].str.split('_').str[-1]
    data['experiment'] =  "+" + data['experiment'].str.split('_').str[2]


    #####CALCULATE VARIABLE: % no hits
    data[['hits']] = data[['hits']].apply(pd.to_numeric)

    results = data.groupby(['gen marker','experiment', 'input data'])['hits'].mean().reset_index().rename(columns={"hits": 'average hits'})
    results = results.set_index(['gen marker','experiment', 'input data'])
    #results = data.join(data.groupby(['gen marker','experiment', 'input data'])['hits'].mean().reset_index().rename(columns={0: '%'}) ,
                #   on = ['gen marker', 'experiment'] , rsuffix = '_%')

    #####CALCULATE VARIABLE: from hits --> % pathogen, avg coverage pathogen and avg evalue
    #select only the 'hits' and make data numeric
    hits_data = data.loc[data['hits'] > 0 ]
    hits_data[['% pathogen', 'avg coverage pathogen', 'avg evalue pathogen']] = hits_data[['% pathogen', 'avg coverage pathogen', 'avg evalue pathogen']].apply(pd.to_numeric)

    #calculate variables
    ###### percentage pathogen
    var_perc_pathogen = hits_data.groupby(['gen marker', 'experiment', 'input data'])['% pathogen'].mean().reset_index(). rename(columns={ '% pathogen': 'avg. % pathogen hits' })
    var_perc_pathogen = var_perc_pathogen.set_index(['gen marker','experiment', 'input data'])
    #add to results
    results = results.merge(var_perc_pathogen, on = ['gen marker','experiment', 'input data'], how = 'outer')
    ####### average evalue
    var_avg_eval = hits_data.groupby(['gen marker', 'experiment', 'input data'])['avg evalue pathogen'].mean().reset_index(). rename(columns={ 0: 'avg evalue pathogen' })
    #print(var_avg_eval)
    var_avg_eval = var_avg_eval.set_index(['gen marker','experiment', 'input data'])
    #add to results
    results = results.merge(var_avg_eval, on = ['gen marker','experiment', 'input data'], how = 'outer')
    ####### average coverage
    #only select data with % pathogen > 0
    path_data = hits_data.loc[hits_data['% pathogen'] > 0 ]
    var_avg_coverage = path_data.groupby(['gen marker', 'experiment', 'input data'])['avg coverage pathogen'].mean().reset_index().rename(columns={ 0: 'avg coverage pathogen' })
    #var_avg_coverage = var_avg_coverage.loc[var_avg_coverage['avg. % pathogen hits'] > 0]
    var_perc_pathogen = var_avg_coverage.set_index(['gen marker','experiment', 'input data'])
    #add to results
    results = results.merge(var_avg_coverage, on = ['gen marker','experiment', 'input data'], how = 'outer')

    return results


def results_plot(data, output, gen_marker_group):
    hue_order_list = ['18S', '28S', '5S', 'ITS', 'CAL', 'EF1', 'RPB1', 'RPB2', 'TUB2_1-4', 'TUB2_4', 'TUB2_4-5']

    sns.set(font_scale = 1.7)
    plots = ['avg. % pathogen hits', 'avg coverage pathogen']
    for plot in plots:
        if plot == 'avg. % pathogen hits':
            p = sns.catplot(x = 'experiment', y = plot, data = data, hue = 'input data' ,
                kind = 'swarm', col = 'gen marker', col_wrap = 4, s = 10, legend_out = False ,
                col_order = hue_order_list )
            p.set_titles("{col_name}")
            p.set_axis_labels("relatives added to input data", 'avg. % pathogen hits')
            #plt.show()
            plt.savefig(output + str(gen_marker_group), dpi=125)



def get_gen_markers_groups(data):
    gen_markers_group_1 = np.array(['18S', '28S', '5S', 'ITS'])
    index_group_1 = np.where(np.isin(data['gen marker'].unique(), gen_markers_group_1))
    gen_markers_group_2 = np.delete(data['gen marker'].unique(), index_group_1)

    gen_marker_groups_list = [gen_markers_group_1 , gen_markers_group_2 ]

    return gen_marker_groups_list

if __name__ == '__main__':
    import argparse
    from argparse import RawDescriptionHelpFormatter

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    io_group = ap.add_mutually_exclusive_group(required=True)

    io_group.add_argument('-i', '--input', type=str, help='input file(s)', metavar='input file(s) in txt format. Command *.txt', nargs='+')
    ap.add_argument('-o', '--output', type=str, help='Folder path', metavar='Output (path and) filename')

    cmd = ap.parse_args()

    print(cmd.input, cmd.output)
    #filepath = get_filepath()

    #important parameters
    pathogen = 'Fusarium circinatum\n'
    #write columuns of results table
    f = open(cmd.output, 'w')
    f.write('consensus;hits;% pathogen;avg evalue pathogen;avg coverage pathogen;gen marker;experiment\n')
    f.close()

    concatinate_files(cmd.input[0], cmd.output, pathogen)

    #analyse results
    results = analyse(cmd.output)

    list_gen_markers_groups = get_gen_markers_groups(results)
    #Plot results table
    for group in list_gen_markers_groups:
        table = results_plot(results, cmd.output, group)
