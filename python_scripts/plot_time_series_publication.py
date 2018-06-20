import pandas as pd
import re
from collections import Counter
from scipy import stats
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams
import os
from statsmodels.sandbox.stats.multicomp import multipletests 

# import metadata and data
meta = pd.read_table('group4/MScourse2018_Global_metadata_FINAL_JMG_6.11.18.txt', encoding = 'ISO-8859-1', index_col='filename')
meta.reset_index(inplace=True)
meta.rename(columns = {'filename':'#SampleID'}, inplace=True)

df = pd.read_csv('group4/FeatureTable_G4_yogurt_Clean.csv', encoding = 'ISO-8859-1', index_col='filename')
feat_start = np.where(df.columns.str.contains('^\d'))[0][0]
df = df[df.columns[feat_start:]]
df = df.apply(lambda a: a/sum(a))
df.reset_index(inplace=True)

# subselect metadata for studies of interest
meta = meta[['#SampleID', 'ATTRIBUTE_age', 'ATTRIBUTE_distributor_brand']]
matching_ids = pd.merge(meta, df, left_on='#SampleID', right_on='filename')

# select the target subset
kroger = matching_ids[matching_ids['ATTRIBUTE_distributor_brand']=='homemade (Horizon milk+Kroger yogurt)']
okroger = matching_ids[(matching_ids['ATTRIBUTE_distributor_brand']=='Kroger') & (matching_ids['ATTRIBUTE_age']==' t=0')]
okroger['ATTRIBUTE_age'] = 'Original' 
kroger = pd.concat([kroger, okroger])

# Calculate correlation
meanlist_kroger = [pd.DataFrame(kroger.groupby('ATTRIBUTE_age')[x].mean()) for x in kroger.columns[4:]]
meantab_kroger = pd.concat(meanlist_kroger, axis=1)
cortab_kroger = meantab_kroger.apply(lambda a: stats.spearmanr(a, range(7)))
pval_kroger = [x[1] for x in cortab_kroger]
#pval_kroger = multipletests(pval_kroger, method = 'fdr_bh')[1]

# Load associated GNPS note attributes
gnps = pd.read_table('group4_1f11fb76/clusterinfosummarygroup_attributes_withIDs/a736c603ddba4161852e254a0215580b..out')

# Match MS1 to MS2 in GNPS
def match2gnps(mz, rt, gnps, mzppm=20, rttol=20):
    idx = gnps.index[(((abs(gnps['parent mass']-mz)/mz)*10**6) < mzppm) & (abs(gnps['RTMean']-rt) < rttol)]
    idx = [str(x) for x in idx]
    if len(idx):
        return ','.join(idx)
    else:
        return ''

# Load associated network edges 
net = pd.read_table('group1_1520f711/networkedges/2adf392abab64b2db35ed2d15cafda80.pairsinfo', header=None)

# Add significance and score threshold
feat_list = []
for feat in cortab_kroger[((np.array(cor_kroger) > 0.7) | (np.array(cor_kroger) < -0.7)) & 
                         (np.array(pval_kroger)<0.05)
		         ].index:
    cidx, mz, rt = feat.split('_')
    mids = match2gnps(float(mz), float(rt)*60, gnps) 
    fn = 'yogurt/'+feat+'.png'
    ycor, ypval = cortab_kroger[feat] 
    feat_list.append([feat, mids, fn, ycor, ypval])

feat_tab = pd.DataFrame(feat_list)
feat_tab.columns = ['feat_id', 'gnps_pos', 'plot', 'kcor', 'kpval'] 

os.mkdir('yogurt')

# Plot boxplots for all selected features
xorder = meantab_kroger.index.tolist()
for feat in feat_tab['feat_id']: 
    ax = sns.boxplot(x='ATTRIBUTE_age', y=feat, 
                     data=kroger,
                    order=xorder)
    plt.xticks(rotation=30)
    #pdf.savefig()  # saves the current figure into a pdf page
    fig = ax.get_figure()
    fn = 'yogurt/'+feat+'.png'
    fig.savefig(fn) 
    plt.close()


def getLibraIDlink(i, tab):
    burl = 'https://gnps.ucsd.edu/ProteoSAFe/result.jsp?task='
    purl = '&view=view_all_annotations_DB#%7B%22main.Compound_Name_input%22%3A%22'
    eurl = '%22%7D'
    gid = tab.loc[0,'ProteoSAFeClusterLink'].split('&')[0].split('=')[1]
    if type(tab['LibraryID'][i]) == float:
        return '<a href=\"'+ tab.loc[i,'ProteoSAFeClusterLink'] + '\" target=\"_blank\">'+'N/A'+'</a>'
    else:       
        return '<a href=\"'+ burl+gid+purl+str(tab.loc[i,'LibraryID'])+eurl+ '\" target=\"_blank\">'+str(tab.loc[i,'LibraryID'])+'</a>'


# attach GNPS network component for each node
feat_tab['component'] = ''

for i in range(feat_tab.shape[0]):
    a = feat_tab.loc[i,'gnps_pos']
    if a!='': 
        tmp = [] 
        for x in list(a.split(',')):
            y = gnps.loc[int(x), 'cluster index']
            if y in netnode:
                tmp.append(str(net.loc[(net[0]==y) | (net[1]==y), 6].reset_index(drop=True)[0]))
        feat_tab.loc[i,'component'] = '#'.join(tmp)

def getComponentlink(comp):
    burl = 'https://gnps.ucsd.edu/ProteoSAFe/result.jsp?view=network_displayer&componentindex='
    purl = '&task=1520f7112d9f445384eb743ac4358c21#%7B%7D'
    return '<a href=\"'+ burl+comp+purl+ '\" target=\"_blank\">'+comp+'</a>'

# format and save table
feat_tab['gnps_pos'] = feat_tab['gnps_pos'].apply(lambda a: '#'.join([getLibraIDlink(int(x), gnps) for x in list(a.split(','))]) if a!='' else '')
feat_tab['component'] = feat_tab['component'].apply(lambda a: '#'.join([getComponentlink(x) for x in list(a.split(','))]) if a!='' else '')
feat_tab['plot'] = feat_tab['plot'].apply(lambda a: '<img src=\"'+a+'" />')
feat_tab.to_csv('yogurt.txt', sep='\t', index=None)
