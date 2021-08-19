"""
Written by Seyma Unsal Beyge
NETLAB
supervised by Assoc.Prof.Dr.Nurcan Tuncbag
Middle East Technical University - Informatics Institute

script to prepare interactome file that is necessary for network modeling
"""
import networkx as nx
import pandas as pd
from cmapPy.pandasGEXpress.parse import parse
import argparse, os, sys

def get_column_info(colname,docname):
    """
    takes column name and path to the L1000 inst file
    returns:
        [0]: cell type
        [1]: drug name
        [2]: drug dosage (concentration)
        [3]: drug introduction time
    """
#    import pandas as pd
    data = pd.read_csv(docname, sep='\t')
    cell = str(data[data['inst_id']==colname]['cell_id'].item())
    drug = str(data[data['inst_id']==colname]['pert_iname'].item())
    dosage = str(data[data['inst_id']==colname]['pert_dose'].item())
    time = str(data[data['inst_id']==colname]['pert_time'].item())
    time += str(data[data['inst_id']==colname]['pert_time_unit'].item())
    
    return cell,drug,dosage,time


def get_same_columns_l1000(cell_of_interest,drug_of_interest,col_list,docname):
    """
    takes drug name, cell name and the total list of columns
    returns the list of columns that dosage and time information is same
    """
    cols={}
    for col in col_list:
        info = get_column_info(col,docname)
        if info[0] == cell_of_interest and info[1] == drug_of_interest:
            otherinfo = info[2]+'_'+info[3]
            if otherinfo not in cols.keys():
                cols[otherinfo]=[col]
            else:
                cols[otherinfo]+=[col]
    return cols


def main():
    #Parsing arguments 
    # (run python 3.prepare_interactome.py -h to see all these decriptions)
    parser = argparse.ArgumentParser(
        description='Given the drugname, cell line name of interest, path to '\
        'the reference interactome paths to the raw L1000 data, gene-inst-pert '\
        'info files, path to the file including list of landmark genes  '\
        'processes the reference interactome to exclude hub nodes, low-expressed genes, '\
        'then perfroms link prediction and adds edges passing localization filters. '\
        'finally outputs the interactome file as "{cell_line}_{drugname}_iref_processed_interactome.txt".')
    #required arguments
    parser.add_argument("-c", "--cell", dest='cellline', help='Cell line of interest '\
                        'Choices: "A375","A549","MCF7","PC3","YAPC","NPC". ')
    parser.add_argument("-d", "--drug", dest='drugname', help='Drug of interest')
    
    #optional arguments
    parser.add_argument("-i", "--int", dest='interactome', help='path to the reference interactome ',
                        default='./raw_data/iref_mitab_miscore_2013_08_12_interactome_no_selfloop.txt')

    parser.add_argument("-lg", "--Lgene", dest='LgeneFile', help='Path to the L1000 gene info file '\
        'Default = "./raw_data/GSE101406_Broad_LINCS_L1000_gene_info.txt"', 
        default='./raw_data/GSE101406_Broad_LINCS_L1000_gene_info.txt')
        
    parser.add_argument("-la", "--Lanalyte", dest='LanalyteFile', help='Path to the L1000 analyte info file '\
        'Default = "./raw_data/GSE101406_Broad_LINCS_L1000_inst_info.txt"', 
        default='./raw_data/GSE101406_Broad_LINCS_L1000_inst_info.txt')

    parser.add_argument("-L", "--L1000", dest='l1000File', help='Path to the raw L1000 data file '\
        'Default = "./raw_data/GSE101406_Broad_LINCS_L1000_Level3_INF_mlr12k_n1667x12328.gctx"', 
        default='./raw_data/GSE101406_Broad_LINCS_L1000_Level3_INF_mlr12k_n1667x12328.gctx')
    
    parser.add_argument("--localization",  dest='localization', help='Path to the localization info file. '\
        'Default = "./raw_data/subcellular_location.txt"', 
        default='./raw_data/subcellular_location.tsv')
        
    parser.add_argument("--expthreshold", dest="expthreshold", default=2.0,
                       help='The level of expression value that will used for '\
                           'excluding the low expressed genes. Default is 2.0.')
        
    parser.add_argument("--outpath", dest = 'outputpath', help='Path to the directory which will '\
        'hold the output files. Default = "./raw_data/interactomes/{drugname}/"', default='./raw_data/interactomes/')
    
    parser.add_argument("--lpmethod", dest='lpmethod', default='Adamic/Adar',
                help='The link prediction method to use for interactome preparation. '\
                'Default is "Adamic/Adar". Options are "Adamic/Adar", '\
                '"Jaccards", "Preferential Attachment", "Resource Allocation" and None. '\
                'If None, any link prediction is not used and the preparation procedure is terminated.')
        
    options = parser.parse_args()
    #Check if outputpath exists
    outputpath=options.outputpath+options.drugname+'/'
    if not os.path.isdir(outputpath):
        # sys.exit('Outpath %s is not a directory'%options.outputpath)
        os.mkdir(outputpath)
    hubs=['UBC', 'APP', 'ELAVL1', 'SUMO2', 'CUL3']
    exp_threshold=float(options.expthreshold)
    
    ### CREATE DATAFRAMES for SOURCE DATA:
    l1000_metadata = parse(options.l1000File)
    l1000_data = l1000_metadata.data_df
    # lm_genes = pd.read_csv(options.lmGenesFile, sep='\t')
    genes = pd.read_csv(options.LgeneFile, sep='\t')
    # drug_data = pd.read_csv(options.pertFile, sep='\t')
    # druglist = drug_data['pert_iname'].tolist()
    
    # lmgenes=list(lm_genes['Entrez ID'])
    l1000_data.index=[int(x) for x in l1000_data.index]
    # if options.iflandmark:
    #     l1000_data=l1000_data.loc[lmgenes]

    drug=options.drugname
    cell=options.cellline
    print ('\nWorking on the interactome for the drug {0} in cell line {1}'.format(drug,cell))
    
    dmso = get_same_columns_l1000(cell,'DMSO',l1000_data.columns, options.LanalyteFile)
    dmso_key=list(dmso.keys())
    try:
        if len(dmso_key)==1:
           dmso_cols=dmso[dmso_key[0]]
    except:
        sys.exit('There is more than one type of DMSO treatment options')

    subdmso =l1000_data[dmso_cols]
    dmso_cols=subdmso.columns
    
    dmso_col = subdmso.loc[: , dmso_cols[0]:dmso_cols[-1]]
    subdmso['DMSO_ave'] = dmso_col.mean(axis=1)

    drg = get_same_columns_l1000(cell, drug, l1000_data.columns, options.LanalyteFile)
    drg_keys=list(drg.keys())
    try:
        if len(drg_keys)==1:
           drg_cols=drg[drg_keys[0]]
    except:
        sys.exit('There is more than one type of drug treatment options')

    sub_drg=l1000_data[drg_cols]
    drg_cols = sub_drg.loc[: , drg_cols[0]:drg_cols[-1]]
    sub_drg['DRUG_ave'] = drg_cols.mean(axis=1)
    
    print('\nNumber of drug sampling (',drug,'):',len(drg[drg_keys[0]]))
    print('Average of the expression of genes are calculated to use for filtering.\n')

    gene_ids=list(l1000_data.index)
    gene_names=[genes[genes['pr_gene_id']==x]['pr_gene_symbol'].item() for x in gene_ids]

    newframe=pd.concat([subdmso,sub_drg], axis=1, sort=False)
    newframe.index=gene_names
    
    print('...\n...Starting to process reference interactome.\n...')
    df=pd.read_csv(options.interactome, sep='\t', header=None)
    df.columns=['Node1', 'Node2', 'Weight']

    G=nx.Graph()
    edges=[list(df.iloc[x]) for x in df.index]
    #print edges[:5]
    G.add_weighted_edges_from(edges)
    print ('\nNumber of edges in the interactome:',len(G.edges()),'\n')
    
    all_nodes=G.nodes()
    new_interactome_path='{0}{1}_{2}_iref_processed_interactome_v1_without_artificial_edges.txt'.format(outputpath, cell,drug)
    filtered_edges=[]
    hub_edge_count=0
    low_expressed_count=0
    for edge in edges:
        n1,n2=edge[0],edge[1]
        if n1 in hubs or n2 in hubs:
            hub_edge_count+=1
        else:
            if n1 in gene_names:
                exp_n1=newframe['DMSO_ave'].loc[n1]
            else:
                exp_n1=exp_threshold
            if n2 in gene_names:
                exp_n2=newframe['DMSO_ave'].loc[n2]
            else:
                exp_n2=exp_threshold
            
            if exp_n1>=exp_threshold and exp_n2>=exp_threshold:
                filtered_edges.append(edge)
            else:
                low_expressed_count+=1
                G.remove_edge(n1,n2)
                
    print ('Number of edges after filtering hubs and the ones with expression values less than given threshold:\n',len(filtered_edges))
    print ('Number of edges with hubs:', hub_edge_count)
    print ('low_expressed_count:',low_expressed_count)
    
    new_interactome_df=pd.DataFrame(columns=['Protein1','Protein2','Weight'])
    data=open(new_interactome_path,'w')
    for edge in filtered_edges:
        new_interactome_df.append({'Protein1':edge[0],'Protein2':edge[1],'Weight':G.get_edge_data(*edge)['weight']}, ignore_index=True)
        word='%s\t%s\t%s\n'%(edge[0],edge[1],str(G.get_edge_data(*edge)['weight']))
        data.write(word)
    data.close()
    print ('\nProcessed interactome before application of link prediction is written to the file:{0}..\n'.format(new_interactome_path))
    
    docname='{0}{1}_{2}_log.txt'.format(outputpath, cell,drug)
    wrt=open(docname,'w')
    wrt.write(drug+':\n\n')
    landmark=set(gene_names)
    common=set(all_nodes).intersection(landmark)
    print ('Total gene number found in L1000 data:',len(landmark))
    wrt.write('Total gene number found in L1000 data: '+str(len(landmark))+'\n')
    print ('Total gene number found in the interactome:', len(all_nodes))
    wrt.write('Total gene number found in the interactome: '+str(len(all_nodes))+'\n\n')
    print ('Number of intersecting genes (L1000 genes and interactome genes):',len(common))
    wrt.write('Number of intersecting genes (L1000 genes and interactome genes): '+str(len(common))+'\n')
    wrt.write('\nNumber of edges after filtering hubs and those expression values are less than given threshold: %s\n'%len(filtered_edges))
    wrt.write('\nNumber of edges with hubs: %s\n'%hub_edge_count)
    wrt.write('\nNumber of edges with low expressed genes: %s\n'%low_expressed_count)

    l1 = 'DMSO - min: '+str(min(newframe['DMSO_ave']))+ ' - max: '+ str(max(newframe['DMSO_ave']))+'\n'
    l2 = drug+'- min: '+str(min(newframe['DRUG_ave']))+ ' - max: '+ str(max(newframe['DRUG_ave']))+'\n'
    print (l1,l2)
    print('Threshold used to eliminate nodes as low expressed:\t{0}'.format(options.expthreshold))
    wrt.write(l1)
    wrt.write(l2)
    wrt.close()
    
    ###########################################################################
    print('\n...\n...Starting to the application of selected link prediction algorithm to the interactome.\n...\n... ')
    ###########################################################################
    method=options.lpmethod
    if method == "Adamic/Adar":
        print('...Starting to calculate predictions based on the method {0}.\n...'.format(method))
        preds = nx.adamic_adar_index(G)
    elif method == "Jaccards":
        print('...Starting to calculate predictions based on the method {0}.\n...'.format(method))
        preds = nx.jaccard_coefficient(G)
    elif method == "Preferential Attachment":
        print('...Starting to calculate predictions based on the method {0}.\n...'.format(method))
        preds = nx.preferential_attachment(G)
    elif method == "Resource Allocation":
        print('...Starting to calculate predictions based on the method {0}.\n...'.format(method))
        preds = nx.resource_allocation_index(G)
    else:
        sys.exit('None or an undefined method selected. Terminating the process!')
    
    #Sorting the predictions by descending scores
    ilength=len(G.edges())
    print('Length of the current interactome:',ilength)
    sorted_pred=sorted(preds, key=lambda x: x[2], reverse=True)
    
    out='{0}{1}_{2}_edge_predictions.txt'.format(outputpath, cell, drug)
    outdata=open(out,'w')
    for u, v, p in sorted_pred:
        if p>0.00:
            word = '%s\t%s\t%.8f\n' % (u, v, p)
            outdata.write(word)
    outdata.close()
    print('\nPredicted edges with scores higher than 0.0 are written in the file {0}.'.format(out))

    #Create the dataFrame of subcellular localization info 
    locations=pd.read_csv(options.localization, sep='\t')
    locations=locations[['Gene name','Reliability','Main location','Additional location']]
    locations['Additional location'] = locations['Additional location'].fillna('')
    print
    newlocframe=pd.DataFrame(locations['Gene name'])
    newlocframe['all_locations']=locations['Main location']+';'+locations['Additional location']

    ############################################################
    ## sorted_pred: filter according to location information  ##
    ############################################################
    filtered_pred=pd.DataFrame(columns=['Protein1','Protein2','Score'])
    print('\nPredictions are sorted.\nThe same number of edges as the processed interactome are collected in descending order.\n')
    # for i in range(len(ilength)):
    # for pedge in sorted_pred:
    for p1, p2, w in sorted_pred[:ilength]:
        if p1 in hubs or p2 in hubs:
            continue
        p1loc=newlocframe['all_locations'][newlocframe['Gene name']==p1].values
        if len(p1loc)>0:
            p1loc=p1loc[0].split(';')
            p1loc=[x for x in p1loc if x!='']
        p2loc=newlocframe['all_locations'][newlocframe['Gene name']==p2].values
        if len(p2loc)>0:
            p2loc=p2loc[0].split(';')
            p2loc=[x for x in p2loc if x!='']
        p1loc=set(p1loc)
        p2loc=set(p2loc)
        if len(p1loc)==0 or len(p2loc)==0:
            row={'Protein1':p1,'Protein2':p2,'Score':w}
            filtered_pred=filtered_pred.append(row, ignore_index=True)
        else:
            common=p1loc.intersection(p2loc)
            if len(common)>0:
                row={'Protein1':p1,'Protein2':p2,'Score':w}
                filtered_pred=filtered_pred.append(row, ignore_index=True)
    
    print()
    print(filtered_pred.head())
    print(filtered_pred.shape)
    
    print('...\nFrom top scoring predicted edges, the ones with UBC and the ones with no common subcellular localizations are removed.\n...')
    
    predicted_edge_count=len(filtered_pred.index)
    print('Number of predicted edges left:\t{0}\n...'.format(predicted_edge_count))
    
    if len(filtered_pred.index)==0:
        print('There is no predicted edge that has passed the localization filters.\n...')

    max_score=float(filtered_pred['Score'].iloc[0])
    print ('max score of predictions after filtering:\t{0}\nScores are scaled to range between 0.0 and 0.5.'.format(max_score))
    filtered_pred['Weight']=[x/(max_score*2) for x in filtered_pred['Score']]
    
    new_interactome=new_interactome_df.append(filtered_pred, ignore_index=True)
    print('Shape of original interactome:{0}\nShape after predictions are added:{1}'.format(new_interactome_df.shape, new_interactome.shape))
    print ('Number of predicted edges after localization filters:\t{0}'.format(len(filtered_pred.index)))
    
    out='{0}{1}_{2}_LP_method_applied_interactome.txt'.format(outputpath, cell, drug)
    new_interactome[['Protein1','Protein2','Weight']].to_csv(out,sep='\t',header=False, index=False)
    print('Final interactome for {0} cell line and the drug {1} is written in the path:\n{2}'.format(cell,drug,out))
    
if __name__ == '__main__':
    main()
