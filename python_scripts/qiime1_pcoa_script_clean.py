import pandas as pd
import numpy as np
import re, os, time
from collections import Counter
import subprocess 
import shutil
from optparse import OptionParser

def qiime1_pcoa(feat_table, metadata, weight_norm=1, norm='PQN', index_col='filename', 
                qiime1_run='/home/msclass/.conda/envs/qiime1/bin/', dis_method='canberra', **kwargs):

    feat = pd.read_csv(feat_table, encoding = 'ISO-8859-1', index_col=index_col)
    meta = pd.read_table(metadata, encoding = 'ISO-8859-1', index_col=index_col)
    
    meta.reset_index(inplace=True)
    meta.rename(columns = {index_col:'#SampleID'}, inplace=True)
    if sum(meta.columns.str.contains('sample_name')):
        meta.rename(columns = {'sample_name':'sample_name_revised'}, inplace=True)
    meta['#SampleID'] = meta['#SampleID'].str.replace('\.mzxml','')
    mapping = pd.concat([meta['#SampleID'],
                         # create columns with required format
                         pd.DataFrame({'BarcodeSequence': meta.shape[0]*["GATACA"]}),   
                         pd.DataFrame({'LinkerPrimerSequence': meta.shape[0]*["GATACA"]}),
                         # add metadata from original matrix
                         meta[meta.columns[1:]],     
                         pd.DataFrame({'Description': meta.shape[0]*["Metabolome"]})
                        ], axis=1
                        )
    mapping['#SampleID'] = mapping['#SampleID'].str.replace('_','.')
    mapping['#SampleID'] = mapping['#SampleID'].str.replace('\.mzXML','')
    mapping = mapping[mapping['#SampleID']!='not collected']
    mapping.to_csv(metadata+"_mapping_qiime1.txt", index=None, sep='\t')
    
    # assume feature are column header and start with a number 
    feat_start = np.where(feat.columns.str.contains('^\d'))[0][0]
    feat = feat[feat.columns[feat_start:]]
    feat.reset_index(inplace=True)
    
    feat['filename'] = feat['filename'].apply(lambda x: re.sub('_','.',x))
    feat['filename'] = feat['filename'].apply(lambda x: re.sub('\.mzXML','',x))
    
    mapping = pd.merge(feat[feat.columns[:1]], mapping, left_on='filename', right_on='#SampleID') 
    
    print('Matching samples to metadata')
    print(int(mapping.shape[0]))
    
    feat.set_index('filename', inplace=True)
    feat = feat.loc[mapping['#SampleID']] 
    
    if weight_norm:
        weight = mapping['sample_in_tube_ethanol_weight']
        weight = weight.astype('float') 
        if sum(weight.isna()):
            print('Found nan in weights')
        if sum(weight==0):
            print('Found 0 in weights')
        feat = feat.apply(lambda a: a/np.array(weight))
    
    if norm=='TIC':
        feat = feat.apply(lambda a: a/sum(a))
        feat = feat.T
        feat.reset_index(inplace=True)
        feat.rename(columns = {'index':'#SampleID'}, inplace=True)
        #feat.replace(np.inf, 0, inplace=True)
        feat.fillna(0, inplace=True)
     
        feat.to_csv(feat_table+"_qiime1.txt", index=None, sep='\t')
    
    if norm=='PQN':
        feat.replace(0, 1, inplace=True)
        #variance = feat.apply(lambda a: np.std(a))
        #idxvec = np.argsort(variance)
        #sel = idxvec[:int(len(idxvec) * 0.3)]
        #ref = feat[feat.columns[sel]].apply(lambda a: np.median(a), axis=1)
        ref = feat.apply(lambda a: np.median(a), axis=1)
        quotient = feat.apply(lambda a: a/ref)
        
        quotient_median = quotient.apply(lambda a: np.median(a[(a.notna()) & (a!=np.inf)]))
        feat = feat.apply(lambda a: a/quotient_median, axis=1)
        
        feat = feat.T
        feat.reset_index(inplace=True)
        feat.rename(columns = {'index':'#SampleID'}, inplace=True)
        feat.replace(np.inf, 0, inplace=True)
        feat.fillna(0, inplace=True)
        feat.to_csv(feat_table+"_qiime1.txt", index=None, sep='\t')
    
    if norm=='None':
        feat = feat.T
        feat.reset_index(inplace=True)
        feat.rename(columns = {'index':'#SampleID'}, inplace=True)
        #feat.replace(np.inf, 0, inplace=True)
        feat.fillna(0, inplace=True)
     
        feat.to_csv(feat_table+"_qiime1.txt", index=None, sep='\t')
     
    print('Samples not present in the metadata')
    print(set(feat.columns[1:])-set(mapping['#SampleID']))
    
    
    subprocess.call([qiime1_run+'biom', 'convert', '-i', feat_table+"_qiime1.txt", 
                     '-o', feat_table+'.biom', '--to-hdf5', '--table-type', 'Metabolite table']) 
    os.remove(feat_table+"_qiime1.txt")

    subprocess.call([qiime1_run+'beta_diversity.py', '-i', feat_table+'.biom',
                    '-m', dis_method, '-o', feat_table+'_'+dis_method])
    os.remove(feat_table+'.biom')
    
    subprocess.call([qiime1_run+'principal_coordinates.py',  '-i', 
                    feat_table+'_'+dis_method+'/'+dis_method+'_'+feat_table+".txt",
                    '-o', feat_table+'_'+dis_method+'_div_coords.txt'])
    shutil.rmtree(feat_table+'_'+dis_method) 
    
    subprocess.call([qiime1_run+'make_emperor.py', '-i', feat_table+'_'+dis_method+'_div_coords.txt',
                    '-m', metadata+"_mapping_qiime1.txt", '-o', 'emperor_'+feat_table+'_'+dis_method])
    os.remove(feat_table+'_'+dis_method+'_div_coords.txt')
    os.remove(metadata+"_mapping_qiime1.txt")
       
def main():
    parser = OptionParser(usage="usage: %prog [options]",
                              version="%prog 1.0")
    parser.add_option("-i", "--input",
	                      help="ClusteApp style feature table")
    parser.add_option("-m", "--metadata",
	                      help="Matching metadata table")
    parser.add_option("-w", "--weight",
		              type="int",
                              default=1,
	                      help="weight normalization")
    parser.add_option("-n", "--norm",
                              default='PQN',
	                      help="weight normalization")
    parser.add_option("-I", "--indexcol",
                              default='filename',
	                      help="Index col, common to both tables")
    parser.add_option("-r", "--run",
                              default='/home/msclass/.conda/envs/qiime1/bin/',
	                      help="qiime1 run time")
    parser.add_option("-d", "--dissimilarity",
	                      help="dissimilarity measure")
    
    (options, args) = parser.parse_args()

    qiime1_pcoa(options.input, options.metadata, weight_norm=options.weight, 
                norm=options.norm, index_col=options.indexcol, 
                qiime1_run=options.run, dis_method=options.dissimilarity)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))
