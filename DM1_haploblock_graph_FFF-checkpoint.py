import argparse
# Initialize parser
parser = argparse.ArgumentParser()

# Adding optional argument
parser.add_argument("--hpt_file", help = "hpt file to use")
parser.add_argument("--out_dir", help = "path to output directory")
parser.add_argument("--patients_file", help = "patients file, 1 patients id per line")
parser.add_argument("--blocks_file", help = "path to blocks file")
parser.add_argument("--prefix", help = "prefix for output files")
parser.add_argument("--red_line", help = "position of the red line in bp")
parser.add_argument("--single_hpt", help = "will generate a graph for each hpt file as model by default, indicate the hpt to use if you want to generate a graph for a single hpt ex. : H001 ")
parser.add_argument("--main_color", help = "mai color for the graph")
parser.add_argument("--list_hpt" , help = "list of hpt to use for the graph, if not specified will use all hpt in the hpt file")
parser.add_argument("--percent", help = "tsv file with the pourcentage of each haplotype in the population with columns : haplo and pourcentage")


# Read arguments from command line 
args = parser.parse_args()


if args.hpt_file:
    print("hpt file : % s" % args.hpt_file)
else :
    raise Exception('need an hpt file use --help for help ')


if args.out_dir:
    print("output directory : % s" % args.out_dir)
else :
    raise Exception('need output directory use --help for help ')

if args.patients_file:
    print("patients_file file : % s" % args.patients_file)
else :
    raise Exception('need a patients_file use --help for help ')


if args.blocks_file:
    print("blocks_file : % s" % args.blocks_file)
else:
    raise Exception('need a patients_file use --help for help ')


if args.prefix:
    print("prefix : % s" % args.prefix)
else:
    args.prefix = ''


if args.red_line:
    print("red_line pos : % s" % args.red_line)
else:
    args.red_line = 0


if args.single_hpt:
    print("single_hpt : % s" % args.single_hpt)
else:
    args.single_hpt = False
    

if args.main_color:
    print("main color : % s" % args.main_color)
else :
    args.main_color = (0.91, 0.702, 0.02)


if args.percent:
    print("list_hpt file : % s" % args.list_hpt)
else :
    args.percent = False















import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colormaps
#from utile import pprint
import re
from collections import defaultdict


patients = pd.read_csv(args.patients_file,header=None)
#patients = pd.read_csv('/lustre03/project/6033529/DM1/results/genotypes/haplotypes/patients_ctrls_hts_possible_dm1_chr19_45000000-47000000.txt')
#patients = pd.read_csv('/home/jacobc/projects/rrg-girardsi/IBD_denovo/data/HPT_patients_clean_45-47mb.txt')
#patients = pd.read_csv('/lustre03/project/6033529/IBD_denovo/data/HPT_patients_clean_sureHPT_45-47mb.txt')

patients = patients[0].tolist()
#patients = patients['x'].tolist()
hpt_file = args.hpt_file
#hpt_file = '/lustre03/project/6033529/DM1/results/genotypes/haplotypes/haps_patients_ctrls_chr19_45000000-47000000.txt'
#hpt_file = '/lustre03/project/6033529/DM1/results/genotypes/haplotypes/haps_patients_ctrls_chr19_40116270-46929325.txt'
#hpt_file = '/lustre03/project/6033529/DM1/results/genotypes/haplotypes/haps_patients_ctrls_chr19_15274281-53809954.txt'
hpt = pd.read_csv(hpt_file,
                  sep='\t',
                  header=None,skiprows=5)

pos = pd.read_csv(hpt_file,
                  sep='\t',
                  header=None,skiprows=1,nrows=1)


hpt.columns = ['haplo'] + pos.iloc[0,1:-1].tolist() + ['ind']

hpt_all = hpt.copy()


hpt = hpt.iloc[:,:-1].drop_duplicates(ignore_index=True).iloc[:100,:]

hpt_unique = hpt_all.iloc[:,:-1].drop_duplicates(ignore_index=True)



def get_blocks(haplo_df,pos_list):
    
    select = np.zeros(haplo_df.shape[1])
    for i in pos_list:
        select[i] = 1
    select = select ==1
    haplo_df = haplo_df.loc[:,select]
        
    return haplo_df.apply(lambda x:(x.sum()),axis=1 )
    
blocksfile = args.blocks_file   
#blocksfile = '/lustre03/project/6033529/DM1/results/genotypes/haplotypes/dm1_CaG_sag_chr19_45000000-47000000_onlyhet_haps.haps.GABRIELblocks'
#blocksfile = '/lustre03/project/6033529/DM1/results/genotypes/haplotypes/dm1_CaG_sag_chr19_45000000-47000000_onlyhet_haps.haps.4GAMblocks'
#blocksfile = '/lustre03/project/6033529/DM1/results/genotypes/haplotypes/haps_patients_ctrls_chr19_40116270-46929325.haps.4GAMblocks'
#blocksfile = '/lustre03/project/6033529/DM1/results/genotypes/haplotypes/haps_patients_ctrls_chr19_15274281-53809954.haps.4GAMblocks'
def num2acgt(seq_list):
    final_list = []
    d = {'1':'A','2':'C','3':'G','4':'T'}
    for seq in seq_list:
        newseq = ''
        for i in range(len(seq)):
            newseq += d[seq[i]]
        final_list.append(newseq)
    return final_list
    
def read_blocks(blockfile):
    df = pd.DataFrame(columns=['blockname','markers','seq'])
    nline = 0
    nblock = 1
    markers = None
    with open(blockfile, 'r') as f:
        lines = f.readlines()
        length = len(lines)
        for line in lines:
            if line.startswith('BLOCK'):
                if markers is not None:
                    tempdf = pd.DataFrame({'blockname':'','markers':'','seq':''},index=[0])
                    tempdf.iloc[0,:] = [blockname,markers,num2acgt(inblock)]
                    
                    df = pd.concat([df,tempdf],ignore_index=True)
                
                var = line.split()
                markers = [int(x) for x in var[3:]]
                blockname = var[0]+str(nblock)
                
                nblock += 1
                inblock = []
            
            if re.match('^[0-4]',line) is not None:
                var = line.split()
                inblock.append(var[0])
                
            if nline == length-1:
                tempdf = pd.DataFrame({'blockname':'','markers':'','seq':''},index=[0])
                tempdf.iloc[0,:] = [blockname,markers,num2acgt(inblock)]

                df = pd.concat([df,tempdf],ignore_index=True)
                
            nline += 1
                
            
    return df
     
blocks = read_blocks(blocksfile)           


## match block color with haplotype
## input :
##        haplo : string of the haplotype name
##        haplo_df : dataframe of all haplotypes
##        blocks : dataframe of blocks from read_blocks function
## output :
##       blocks2 : dataframe of blocks with the right order depending on the haplotype
def match_block_color(haplo,haplo_df,blocks):
    blocks2 = blocks.copy()
    haplo_df = haplo_df.iloc[:,:-1].drop_duplicates(ignore_index=True)
    haplo_df.reset_index(drop=True,inplace=True)
    pos = np.where(haplo_df['haplo'] == haplo)[0][0]
   
    temp_hpt = haplo_df.iloc[pos:pos+2,:]
    
    for i in range(len(blocks2)):
        
        seq = get_blocks(temp_hpt,blocks2.iloc[i,1]).reset_index(drop=True)[0]
        
        if blocks2.iloc[i,2][0] != seq:
            try :
                blocks2.iloc[i,2].remove(seq)
                blocks2.iloc[i,2].insert(0,seq)
            except:
                print('exception when changing order of blocks')
    return blocks2
            
            


############################################################################################################
############################################################################################################


## graph only patients

haplo_list= hpt_all['haplo'].unique().tolist() ## liste des nom des haplotypes

##calculer le ratio de patients/crtls par haplotype ayant des patients
def ratio_patients(haplo,hpt_df):
    hpt_all = hpt_df
    n=0
    for i in hpt_all[ hpt_all['haplo'] == haplo]['ind'].tolist():
    
        if i in patients:
            n+=1
    
    ratio = n/len( hpt_all[ hpt_all['haplo'] == haplo]['ind'])
    if ratio > 0 :
        print(haplo,':',ratio)
        print(n)
'''
for i in haplo_list:
    ratio_patients(i,hpt_all)
   ''' 
    
## haplotypes ayant des patients avec les patients dans une liste 
haplo_ind = hpt_all.groupby(['haplo'],as_index=False)['ind'].apply(list)

### filtrer les haplotypes ayant des patients pour le graph
patients_set = set(patients)
n_patients = len(patients_set)

patients_haplo = []

for i in range(len(haplo_ind)):
    
    ratio = len(set(haplo_ind.iloc[i,1]).intersection(patients_set))/len(set(haplo_ind.iloc[i,1]))
    if ratio>0:
        patients_haplo.append(haplo_ind.iloc[i,0])




def make_haplo_graph(hpt_all,blocks,patients_haplo,MAIN_HAPLO,main_color = args.main_color):
    
    ## changer le main haplotype pour le mettre en premier
    blocks = match_block_color(MAIN_HAPLO,hpt_all,blocks)
    ## filtrer les le df de blocks d'haplotypes ayant des patients pour le graph
    hpt_patiens = hpt_unique[hpt_unique['haplo'].isin(patients_haplo)].reset_index(drop=True)


    plot_list = []
    test_df = hpt_patiens

    for j in range(len(test_df)):

        for i in range(len(blocks)):
            pos = blocks.iloc[i,1]
            start = (min(pos))
            end = max(pos)
            seq = get_blocks(test_df,pos)[j]
            

            d = {blocks.loc[i,'seq'][k] : str(k)  for k in range(len(blocks.loc[i,'seq']) )}
            try :
                plot_list.append(dict(Task = test_df.iloc[j,0],Start= start, Finish = end, Resource =d[seq] ))
            except:
                print('exception')
                plot_list.append(dict(Task = test_df.iloc[j,0],Start= start, Finish = end, Resource ='99' ))
            



        
    cmap = plt.get_cmap("Set3") 
    
    colors = { str(i):cmap(i)  for i in range(0,100)}
    colors['0'] = main_color
    
    #### CLAUDIA SI TU VEUX CHANGER LES COULEURS C'EST ICI #####
    #colors = {'0'  : (0.5,0,0), '1'  : (0,0,0), '2':(0,0,1),....} # mettre au moin autant de couleurs que de blocks différents,
    #                                                              # les couleurs sont en RGB de 0 à 1 
    
    colors['99'] =(1,1,1) ## À GARDER : couleur pour les blocks non trouvés dans le main haplotype 

    legends = []
    num = 0
    for i in colors:
        if num >10:
            break
        
        legends.append(mpatches.Patch(color=colors[i], label=i))
        num +=1
        


    d = defaultdict(lambda : 0)

    df = pd.DataFrame(plot_list)


    for i in plot_list:
        if i['Resource'] == '0':
            d[i['Task']] += i['Finish'] - i['Start']


    max(d, key=d.get)

    test = [(int(d[i]),i) for i in d]
    test.sort(key=lambda a: a[0],reverse=False)

    for i in test:
        plot_list.append(dict(Task = i[1],Start= 0, Finish = 0, Resource = '0' ))

    plot_list.reverse()
        
    red_line = int(args.red_line)
    pos_list = [ int(i) for i in hpt.columns.tolist()[1:]]
    pos_list = [pos_list[0]] + pos_list
    
    
    print('mordam')
    print(patients_haplo)    
   
   

    
    #fig.rcParams.update({'font.size': 6})
    #for test in plot_list:
        #if test['Task'] == 'H001' or test['Task'] == 'H119' or test['Task'] == 'H2190'  or test['Task'] == 'H3247'or test['Task'] == 'H2764' or test['Task'] == 'H2443' or test['Task'] == 'H2792':
        #print(test['Start'],test['Finish'])
        #ax1.hlines(test['Task'],pos_list[test['Start']]/1000000,pos_list[test['Finish']]/1000000,colors=colors[test['Resource']],lw = 4)
    #ax2 = ax1.twinx()
    
    #if args.percent != False:
    percent = pd.read_csv(args.percent,delim_whitespace=True,header=0)
    percent.set_index('haplo', inplace=True)
    print('mordam2')
    print(percent)
        #print("Ctrl column values:")
        #print(percent['ctrl'].values)
        # Convert the index to strings to ensure correct matching
        #percent.index = percent.index.astype(str)
        #haplo = percent.index
    # Print the unique haplo values from the percent file
        #percent_haplo_values = percent.index.unique()
        #print("Unique haplo values in percent file:", percent_haplo_values)
    
    # Print the unique haplo values from patients_haplo
        #patients_haplo_values = set(patients_haplo)
        #print("Unique haplo values in patients_haplo:", patients_haplo_values)

    # Check for matching values
        #matching_values = set(percent_haplo_values).intersection(patients_haplo_values)
        #print("Matching haplo values:", matching_values)
    
    # Check for non-matching values
        #non_matching_percent = set(percent_haplo_values) - patients_haplo_values
        #non_matching_patients = patients_haplo_values - set(percent_haplo_values)
        #print("Non-matching haplo values in percent file:", non_matching_percent)
        #print("Non-matching haplo values in patients_haplo:", non_matching_patients)

        #print(haplo)

    plt.clf()
    fig, ax1 = plt.subplots()
    
    #plt.figure(figsize=(8,3))
    fig.set_figheight(16)
    fig.set_figwidth(8)
    percent_is_plot = dict()


    for test in plot_list:
        ax1.hlines(test['Task'],pos_list[test['Start']]/1000000,pos_list[test['Finish']]/1000000,colors=colors[test['Resource']],lw = 4)
    ax2 = ax1.twinx()

    for haplo in percent.index:
        if haplo in (patients_haplo):
            value1 = percent.loc[haplo, 'ctrl']  # Value from the "ctrl" column
            print(value1)
            value2 = percent.loc[haplo, 'patient']   # Value from the "patient" column
            ax1.text(y=haplo,x= 29.15, s=f'{value1}\\{value2}', fontsize=4, verticalalignment='center', color='black')
                
                
            #plt.text(y=percent.iloc[i,0],x= 47_140_000 , s= f'{round(percent.iloc[i,1],3)} %', fontsize=6)
            #plt.axhline(y=percent.iloc[i,0], color='black', linestyle='--',lw=1)
    
    if red_line > 0:
        ax1.axvline(x=red_line/1000000, color='red', linestyle='--',lw=0.5)
        ax2.text(x=(red_line/1000000) -0.1,y = 1.05, s= 'ARSACS 8 Mb haploblocks',)
        
    ax1.set_ylabel('Haplotype')
    ax1.set_xlabel('Position on chromosome 13 (Mb)')
    
    ax2.set_ylabel('carriers\\patients Count', color='black')
    ax2.set_yticklabels([])
    ax2.set_yticks([])
    plt.subplots_adjust(left=0.5, right=1.5, top=0.90, bottom=0.15)
    fig.tight_layout()
    
    
    
    fig_name = f"{args.out_dir}/{args.prefix}plt_gant_"+MAIN_HAPLO+"8Mb.png"
    #plt.legend(handles= legends, bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0,prop={'size': 6})
    fig.savefig(fig_name,dpi = 500)
    '''
    plt.clf()
    
    
    plt.rcParams.update({'font.size': 5})
    for test in plot_list:
        if test['Task'] == 'H001' or test['Task'] == 'H025' or test['Task'] == MAIN_HAPLO :
        
            plt.hlines(test['Task'],test['Start'],test['Finish'],colors=colors[test['Resource']],lw = 10)
    
    if repeat_pos > 0:
        plt.axvline(x=repeat_pos, color='red', linestyle='--',lw=1)
    
    fig_name = "../results/patients_haploblocks_40mb/plt_gant_patients_"+MAIN_HAPLO+"_HGA_HGB.png"
    plt.legend(handles= legends, bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0,prop={'size': 6})
    plt.savefig(fig_name,dpi = 500)
    '''
    # Export haploblock segments for unified R plot
    df_export = pd.DataFrame(plot_list)
    df_export.rename(columns={'Task': 'patient', 'Start': 'start', 'Finish': 'end', 'Resource': 'color'}, inplace=True)
    df_export.to_csv(f"{args.out_dir}/{args.prefix}_haploblock_segments4MBhaplotypes.tsv", sep='\t', index=False)

    return None

    
    
if args.single_hpt != False:
    make_haplo_graph(hpt_all,blocks,patients_haplo,args.single_hpt)
    
else: 
    
    make_haplo_graph(hpt_all,blocks,patients_haplo,patients_haplo[0])



'''
###### get second most common haplogroup #### 
haplo2 = pd.read_csv('/lustre03/project/6033529/DM1/results/genotypes/haplotypes/second_haplogroupe_patients_DM1.txt',sep='\t',header=None)[0].tolist()

hpt_2 = hpt_all[hpt_all['ind'].isin(haplo2)].reset_index(drop=True)
hpt_2['haplo'].unique().tolist()

hpt_all[hpt_all['haplo'] == 'H3102']
hpt_all.iloc[:,200:204]

HGA = ['H001','H674','H2844','H2843','H038','H1891','H679','H2845','H2407','H1579','H2618','H406','H346','H2971']
HGB = ['H025','H061','H3102'] 

hpt_patiens = hpt_all[hpt_all['ind'].isin(patients)].reset_index(drop=True)

hpt_patiens[hpt_patiens['haplo'] == 'H2618']

with open('/lustre03/project/6033529/IBD_denovo/data/HGA_represent_45-47mb.txt','w') as f :
    for i in HGA:
        represent = hpt_patiens[hpt_patiens['haplo'] == i]['ind'].to_list()[0]
        f.write(represent+'\n')



with open('/lustre03/project/6033529/IBD_denovo/data/HGB_represent_45-47mb.txt','w') as f :
    for i in HGB:
        represent = hpt_patiens[hpt_patiens['haplo'] == i]['ind'].to_list()[-1]
        f.write(represent+'\n')





HGA_ind = hpt_all[hpt_all['haplo'].isin(HGA)]['ind'].to_list()
HGB_ind = hpt_all[hpt_all['haplo'].isin(HGB)]['ind'].to_list()
len(HGA_ind)

with open('/lustre03/project/6033529/IBD_denovo/data/HGA_patients_clean_sure_45-47mb.txt','w') as f :
    for i in HGA_ind:
        f.write(i+'\n')
    

with open('/lustre03/project/6033529/IBD_denovo/data/HGB_patients_clean_sure_45-47mb.txt','w') as f :
    for i in HGB_ind:
        f.write(i+'\n')
    
'''
