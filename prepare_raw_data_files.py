"""
Written by Seyma Unsal Beyge
NETLAB
supervised by Assoc.Prof.Dr.Nurcan Tuncbag
Middle East Technical University - Informatics Institute

script to analyze raw files of P100 dataset,
to collect significantly perturbed cell line - drug cases,
and calculate log2FoldChanges 
"""

from cmapPy.pandasGEXpress.parse import parse
import scipy.stats as stats
import sys, argparse, os, math
import numpy as np
import pandas as pd

###############################################################################

def get_column_info(colname,docname):
    """
    takes column name and document name that includes inst info
    returns:
        [0]: cell type
        [1]: drug name
        [2]: drug dosage (concentration)
        [3]: drug introduction time
    """
    data = pd.read_csv(docname, sep='\t')
    cell = str(data[data['inst_id']==colname]['cell_id'].item())
    drug = str(data[data['inst_id']==colname]['pert_iname'].item())
    dosage = str(data[data['inst_id']==colname]['pert_dose'].item())
    time = str(data[data['inst_id']==colname]['pert_time'].item())
    time += str(data[data['inst_id']==colname]['pert_time_unit'].item())
    return cell,drug,dosage,time


def get_same_columns(cell_of_interest,drug_of_interest,col_list,docname):
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


def find_Ftest_results(first_dataframe,second_dataframe):
    """
    takes two dataframes consisting replicates of 
    cell line - drug condition of interest
    and runs columnwise F test between them
        [0]: F value
        [1]: p value
    returns: pvalues
    """
    F_dic = {}
    for col in first_dataframe.columns:
        first=first_dataframe[col].tolist()
        second=second_dataframe[col].tolist()
        res = stats.f_oneway(first,second)
        if col not in F_dic.keys():
            F_dic[col] = [res[1]]
        else:
            F_dic[col] += [res[1]]
    return F_dic


def find_sign_genes(anova_results_doc,gene_doc,pval):
    """
    takes a csv document that includes a dataframe with p values of gene ids
    and a txt document that includes gene ids vs gene names
    returns a list of significant genes
    """
    fdata = pd.read_csv(anova_results_doc)
    length=len(fdata[fdata.columns[0]])
    sign_geneids=[]
    for i in range(length):
        for col in fdata.columns[1:]:
            if fdata[col][i] <=pval:
                geneid = str(fdata[fdata.columns[0]][i])
                if geneid not in sign_geneids:
                    sign_geneids.append(geneid)
    return sign_geneids

def main():
    #Parsing arguments 
    # (run python 1.prepare_raw_data_files.py -h to see all these decriptions)
    parser = argparse.ArgumentParser(
        description='Find significantly perturbed genes and/or phosphoproteins '\
        'and create the files '\
        '(i) p values of each gene/phosphoprotein in each significantly perturbed condition '\
        '(ii) log2FoldChanges of those significantly perturbed gene/protein list ')
    #required arguments
    parser.add_argument("-c", "--cell", dest='cellline', help='Cell line of interest '\
                        'Choices: "A375","A549","MCF7","PC3","YAPC","NPC". ')
    
    #optional arguments
    parser.add_argument("-p", "--pval", dest='pvalThreshold', help='p-value threshold '\
        'that is desired to get the significant conditions (Default is 0.05)',
        default=0.05)
    parser.add_argument("-L", "--L1000", dest='l1000File', help='Path to the raw L1000 data file '\
        'Default = "./raw_data/GSE101406_Broad_LINCS_L1000_Level3_INF_mlr12k_n1667x12328.gctx"', 
        default='./raw_data/GSE101406_Broad_LINCS_L1000_Level3_INF_mlr12k_n1667x12328.gctx')
    parser.add_argument("-P", "--P100", dest='P100File', help='Path to the raw P100 data file '\
        'Default = "./raw_data/GSE101406_Broad_LINCS_L1000_Level3_INF_mlr12k_n1667x12328.gctx"', 
        default='./raw_data/GSE101406_Broad_LINCS_P100_Level3_QCNORM_n1684x96.gctx')
        
    parser.add_argument("-ci", "--cellinfo", dest='cellFile', help='Path to the cell info file '\
        'Default = "./raw_data/GSE101406_Broad_LINCS_cell_info.txt"', 
        default='./raw_data/GSE101406_Broad_LINCS_cell_info.txt')
    parser.add_argument("-pt", "--pert", dest='pertFile', help='Path to the perturbant info file '\
        'Default = "./raw_data/GSE101406_Broad_LINCS_pert_info.txt"', 
        default='./raw_data/GSE101406_Broad_LINCS_pert_info.txt')
    parser.add_argument("-lg", "--Lgene", dest='LgeneFile', help='Path to the L1000 gene info file '\
        'Default = "./raw_data/GSE101406_Broad_LINCS_L1000_gene_info.txt"', 
        default='./raw_data/GSE101406_Broad_LINCS_L1000_gene_info.txt')
    parser.add_argument("-la", "--Lanalyte", dest='LanalyteFile', help='Path to the L1000 analyte info file '\
        'Default = "./raw_data/GSE101406_Broad_LINCS_L1000_inst_info.txt"', 
        default='./raw_data/GSE101406_Broad_LINCS_L1000_inst_info.txt')
    parser.add_argument("-pi", "--Pinst", dest='PinstFile', help='Path to the P100 instrument info file '\
        'Default = "./raw_data/GSE101406_Broad_LINCS_P100_inst_info.txt"', 
        default='./raw_data/GSE101406_Broad_LINCS_P100_inst_info.txt')
    parser.add_argument("-pa", "--Panalyte", dest='PanalyteFile', help='Path to the P100 analyte info file '\
        'Default = "./raw_data/GSE101406_Broad_LINCS_P100_analyte_info.txt"', 
        default='./raw_data/GSE101406_Broad_LINCS_P100_analyte_info.txt')
    parser.add_argument("-lm", "--lmGenes", dest='lmGenesFile', help='Path to the L1000 Landmark genes file '\
        'Default = "./raw_data/landmark_genes_978.txt"', 
        default='./raw_data/landmark_genes_978.txt')
    parser.add_argument("--iflandmark", dest='iflandmark', default=False,
                        help='True or False to select which genes will be used in '\
                            'transcriptomic data (L1000). If True, only landmark genes ' \
                            'will be used in the analysis. Landmark genes are ' \
                            'listed in the file whose path can be given with --lmGenesoption.')
    
    parser.add_argument("--outpath", dest = 'outputpath', help='Path to the directory which will '\
        'hold the output files. Default = "./raw_data/"', default='./raw_data/')
    
    options = parser.parse_args()
    #Check if outputpath exists
    if not os.path.isdir(options.outputpath):
        sys.exit('Outpath %s is not a directory'%options.outputpath)

    ##########--------------------ACCESSING METADATA--------------------###########
    l1000_metadata = parse(options.l1000File)
    l1000_data = l1000_metadata.data_df
    
    p100_metadata = parse(options.P100File)
    p100_data = p100_metadata.data_df
    
    l1000_data.index=[int(x) for x in l1000_data.index]
    if options.iflandmark:
        lm_genes = pd.read_csv(options.lmGenesFile, sep='\t')
        lmgenes=list(lm_genes['Entrez ID'])
        l1000_data=l1000_data.loc[lmgenes]
        print('Landmark genes are selected for further analysis.\n'\
              'Therefore, the following analysis of L1000 data will be '\
              'conducted using only landmark genes using the file:\n',
              options.lmGenesFile)
    
    print ("...\n...accessed to L1000 & P100 metadata...\n...")

    ##########--------------PROCESS DATA & FIND F RESULTS---------------###########
    #############################--------L1000--------#############################
    drug_data = pd.read_csv(options.pertFile, sep='\t')
    druglist = drug_data['pert_iname'].tolist()
    dmso = get_same_columns(options.cellline,'DMSO',l1000_data.columns,options.LanalyteFile)
    dmso_key=list(dmso.keys())
    try:
        if len(dmso_key)==1:
            dmso_cols=dmso['nan_6h']
    except:
        sys.exit('There is more than one type of DMSO treatment options')
        
    L1000_subdmso =l1000_data[dmso_cols]
    subdmsoT=L1000_subdmso.transpose()
    df_list=[]
    print('...\n...Started processing the available drugs in L1000 data...\n')
    for dr in druglist:
        print('Now calculating transcriptomic one way Anova results for drug', dr)
        drug = get_same_columns(options.cellline,dr,l1000_data.columns,options.LanalyteFile)
        for k,v in drug.items():
            subdrug=l1000_data[v]
            subdrugT=subdrug.transpose()
     
            res = find_Ftest_results(subdrugT, subdmsoT)
            df = pd.DataFrame(res, index=[dr+'_'+k])
            dfT=df.transpose()
            df_list.append(dfT)
    
    Lresult = pd.concat(df_list, axis=1, sort=False)
    
    LFout='{0}{1}_L1000_Fresults_pvalues.csv'.format(options.outputpath,options.cellline)
    if options.iflandmark:
        LFout='{0}{1}_L1000_LM_Fresults_pvalues.csv'.format(options.outputpath,options.cellline)
    print ('\nWriting F results of L1000 data...')
    Lresult.to_csv(LFout)
    print ('\nwritten to the document:\n',LFout)
    
    ###############################################################################
    ##########--------------PROCESS DATA & FIND F RESULTS---------------###########
    #############################--------P100--------#############################
    dmso = get_same_columns(options.cellline,'DMSO',p100_data.columns,options.PinstFile)
    dmso_key=list(dmso.keys())
    try:
        if len(dmso_key)==1:
            dmso_cols=dmso[dmso_key[0]]
    except:
        sys.exit('There is more than one type of DMSO treatment options')
    
    subdmso =p100_data[dmso_cols]
    subdmsoT=subdmso.transpose()
    df_list=[]
    print('...\n...Started processing the available drugs in P100 data...\n')
    for dr in druglist:
        print('Now calculating phosphoproteomic one way Anova results for drug', dr)
        drug = get_same_columns(options.cellline,dr,p100_data.columns,options.PinstFile)
        
        for k,v in drug.items():
            subdrug=p100_data[v]
            subdrugT=subdrug.transpose()
            res = find_Ftest_results(subdrugT, subdmsoT)
            df = pd.DataFrame(res, index=[dr+'_'+k])
            dfT=df.transpose()
            df_list.append(dfT)
    
    Presult = pd.concat(df_list, axis=1, sort=False)
    
    PFout='{0}{1}_P100_Fresults_pvalues.csv'.format(options.outputpath,options.cellline)
    print ('\nWriting F results of P100 data...')
    Presult.to_csv(PFout)
    print ('\nwritten to the document:\n',PFout)
    
    ########------PROCESS F DATA & WRITE DOC WITH DIFFERENTIAL VALUES------########
    #############################--------L1000--------#############################
    print('\nfiltering the L1000 raw data for drugs with significant readouts in preferred cell line({0})...'.format(options.cellline))
    fdata=pd.read_csv(LFout,sep=',')
    sign_genes = find_sign_genes(LFout, options.LgeneFile, float(options.pvalThreshold))
    sign_genes = [int(x) for x in sign_genes]
    gene_data = pd.read_csv(options.LgeneFile, sep='\t')
    
    fcols=fdata.columns
    Lcols = l1000_data.columns

    Lcell_cols=[]
    Lcell_cols_revised=[]
    for c in Lcols:
        cell,drug,dosage,time = get_column_info(c, options.LanalyteFile)
        if cell == options.cellline:
            cc='%s_%s_%s'%(drug,dosage,time)
            if cc in fcols:
                Lcell_cols_revised.append(cc)
                Lcell_cols.append(c)
    
    print ('Number of drugs with significant L1000 readouts in cell line {0} is {1}'.format(options.cellline,len(Lcell_cols_revised)))
    l1000_newframe = l1000_data.loc[sign_genes][Lcell_cols]

    genelist = list(l1000_newframe.index)
    genenames=[gene_data[gene_data['pr_gene_id']==int(g)]['pr_gene_symbol'].item() for g in genelist]
    l1000_newframe.index = genenames
    l1000_newframe.columns = Lcell_cols_revised

    filtered_LFout = '{0}{1}_L1000_significant_genes_readouts_p_{2}.txt'.format(options.outputpath, options.cellline, options.pvalThreshold)
    if options.iflandmark:
        filtered_LFout = '{0}{1}_L1000_LM_significant_genes_readouts_p_{2}.txt'.format(options.outputpath, options.cellline, options.pvalThreshold)
    print ('\nWriting the readouts of significant genes... ')
    l1000_newframe.to_csv(filtered_LFout,sep='\t')
    print ('\nwritten to the document:\n',filtered_LFout)
    
    ########------PROCESS F DATA & WRITE DOC WITH DIFFERENTIAL VALUES------########
    #############################--------P100--------#############################
    print('\nfiltering the P100 raw data for drugs with significant readouts in preferred cell line({0})...'.format(options.cellline))
    sign_phosphosites = find_sign_genes(PFout, options.PanalyteFile, float(options.pvalThreshold))

    Pcols = p100_data.columns
    Pcell_cols=[]
    Pcell_cols_revised=[]
    for c in Pcols:
        cell,drug,dosage,time = get_column_info(c, options.PinstFile)
        if cell == options.cellline:
            Pcell_cols.append(c)
            Pcell_cols_revised.append('%s_%s_%s_%s'%(cell,drug,dosage,time))
    
    print ('Number of drugs with significant P100 readouts in cell line {0} is {1}'.format(options.cellline,len(Pcell_cols_revised)))
    
    p100_newframe = p100_data.loc[sign_phosphosites][Pcell_cols]
    p100_newframe.columns = Pcell_cols_revised


    filtered_PFout = '{0}{1}_P100_significant_genes_readouts_p_{2}.txt'.format(options.outputpath, options.cellline, options.pvalThreshold)
    print('Writing the values of union set of significant genes...')
    p100_newframe.to_csv(filtered_PFout,sep='\t')
    print ('\nwritten to the document:\n',filtered_PFout)
    
    ###########################################################################
    ##########--------------- CALCULATION OF Log2FC ---------------############
        #L1000 fold change calculation:
    print('\ncalculating the L1000 Log2FoldChanges of drugs with significant readouts in preferred cell line({0})...'.format(options.cellline))
    l1000_genes=l1000_newframe.index.tolist()
    l1000cols=l1000_newframe.columns
    l1000_col_dict={}
    drugs=[]
    for col in l1000cols:
        ind = col.rfind('.')
        if len(col)-ind == 2:
            drug = col[:ind]
        else:
            drug = col
        if 'DMSO' in drug:
            drug='DMSO'
        if drug not in drugs:
            drugs.append(drug)
        if drug not in l1000_col_dict.keys():
            l1000_col_dict[drug]=[col]
        else:
            l1000_col_dict[drug]+=[col]
    
    l1000_fold_change={}
    l1000_dmso_means = []
    for i in range(len(l1000_genes)):
        dmso_ave=np.mean(l1000_newframe[l1000_col_dict['DMSO']].values.tolist()[i])
        l1000_dmso_means.append(dmso_ave)
    for dr in drugs:
        print('Now processing the L1000 data of',dr)
        dr_fc=[]
        for i in range(len(l1000_genes)):
            dr_ave=np.mean(l1000_newframe[l1000_col_dict[dr]].values.tolist()[i])
            try:
                division = (float(dr_ave)/float(l1000_dmso_means[i]))
                if division<0:
                    fc=-1*(math.log(-1*division,2))
                else:
                    fc=math.log(division,2)
            except (ZeroDivisionError,RuntimeError,RuntimeWarning):
                # print ("can't calculate:",dr_ave,'/',l1000_dmso_means[i])
                fc=0.0
            except:
                # print ("unexpected error",dr_ave,'/',l1000_dmso_means[i])
                fc=0.0
            dr_fc.append(fc)
        l1000_fold_change[dr]=dr_fc

    df_fc = pd.DataFrame(l1000_fold_change,columns=drugs,index=l1000_genes)
    df_fc=df_fc.drop(['DMSO'], axis=1)

    LFCout='{0}{1}_L1000_fold_changes_p_{2}.txt'.format(options.outputpath, options.cellline, options.pvalThreshold)
    if options.iflandmark:
        LFCout='{0}{1}_L1000_LM_fold_changes_p_{2}.txt'.format(options.outputpath, options.cellline, options.pvalThreshold)
        print ('writing the fold change matrix...')
    df_fc.to_csv(LFCout,sep='\t')
    print ('written to the document:\n', LFCout)

        # P100 Fold change calculation:
    print('\ncalculating the P100 Log2FoldChanges of drugs with significant readouts in preferred cell line({0})...'.format(options.cellline))
    p100_genes=p100_newframe.index.tolist()
    p100cols=p100_newframe.columns
    p100_col_dict={}
    drugs=[]
    for col in p100cols:
        ind = col.rfind('.')
        if len(col)-ind == 2:
            drug = col[:ind]
        else:
            drug = col
        if 'DMSO' in drug:
            drug='DMSO'
        if drug not in drugs:
            drugs.append(drug)
        if drug not in p100_col_dict.keys():
            p100_col_dict[drug]=[col]
        else:
            p100_col_dict[drug]+=[col]
    
    p100_fold_change={}
    p100_dmso_means = []
    for i in range(len(p100_genes)):
        dmso_ave=np.mean(p100_newframe[p100_col_dict['DMSO']].values.tolist()[i])
        p100_dmso_means.append(dmso_ave)
    for dr in drugs:
        print('Now processing the P100 data of',dr)
        dr_fc=[]
        for i in range(len(p100_genes)):
            dr_ave=np.mean(p100_newframe[p100_col_dict[dr]].values.tolist()[i])
            try:
                division = (float(dr_ave)/float(p100_dmso_means[i]))
                if division<0:
                    fc=-1*(math.log(-1*division,2))
                else:
                    fc=math.log(division,2)
            except (ZeroDivisionError,RuntimeError,RuntimeWarning):
                # print ("error:",dr_ave,'/',p100_dmso_means[i])
                fc=0.0
            except:
                # print ("unexpected error",dr_ave,'/',p100_dmso_means[i])
                fc=0.0
            dr_fc.append(fc)
        p100_fold_change[dr]=dr_fc

    df_fc = pd.DataFrame(p100_fold_change,columns=drugs,index=p100_genes)
    df_fc=df_fc.drop(['DMSO'], axis=1)
        
    PFCout='{0}{1}_P100_fold_changes_p_{2}.txt'.format(options.outputpath, options.cellline, options.pvalThreshold)
    print ('writing the fold change matrix...')
    df_fc.to_csv(PFCout,sep='\t')
    print ('written to the document:\n', PFCout)

if __name__ == '__main__':
    main()




