"""
Written by Seyma Unsal Beyge
NETLAB
supervised by Assoc.Prof.Dr.Nurcan Tuncbag
Middle East Technical University - Informatics Institute

script to collect seed proteins for the cell line - drug condition of interest
"""


import pandas as pd
import numpy as np
import os, sys, argparse

def get_significant_genes(ftest_doc, gene_doc, pval):
    """
    Parameters
    ----------
    ftest_doc : TYPE str
        DESCRIPTION: path to file with F test pvalues 
    gene_doc : TYPE str
        DESCRIPTION: path to file including gene id & symbol information downloaded with original dataset
    pval : TYPE float
        DESCRIPTION. p-value threshold that will be used for significancy

    Returns
    -------
    sign_dict : TYPE dict
        DESCRIPTION: a dictionary that keys are drugs and values are genes
                     that are significantly perturbed based on the pvalue provided

    """
    gene_info=pd.read_csv(gene_doc, sep='\t')
    df=pd.read_csv(ftest_doc)
    df.index=df[df.columns[0]]
    df=df.drop(df.columns[0], axis=1)
    drugs=df.columns
    sign_dict={}
    for drug in drugs:
        dr=drug.split('_')[0]
        if 'DMSO' in drug:
            continue
        data=df[drug]
        data=data.loc[df[drug]<pval]
        x=list(data.index)
        x=[gene_info[gene_info['pr_gene_id']==int(num)]['pr_gene_symbol'].item() for num in x]
        sign_dict[dr]=x
    return sign_dict

def get_all_tf(tfdoc):
    """
    Parameters
    ----------
    tfdoc : TYPE str
        DESCRIPTION: path to the DBTF file

    Returns
    -------
    newdata : TYPE DataFrame
        DESCRIPTION: A dataframe of all Transcription Factors with their gene symbols and gene names

    """
    tfdata=pd.read_csv(tfdoc, sep='\t')
    newdata=tfdata[tfdata['DbTF']=='yes'][['Gene_symbol','Gene_Name']]
    return newdata

def get_my_tf(genelist,tfdoc):
    """
    Parameters
    ----------
    genelist : TYPE a list of genes
        DESCRIPTION: usually a list of significantly transcribed genes
    tfdoc : TYPE str
        DESCRIPTION. path to the DBTF file

    Returns
    -------
    mytf : TYPE: set
        DESCRIPTION: set of transcription factors intersecting the 
                    significantly transcribed genes and all transcription factors
                     collected from DBTF file

    """
    alltf=get_all_tf(tfdoc)['Gene_symbol'].tolist()
    all_tf=set(alltf)
    genelist=set(genelist)
    mytf=genelist.intersection(all_tf)
    return mytf

def find_related_tf(drugname,sign_gene_dict,regnet_doc,tfdoc):
    """
    Parameters
    ----------
    drugname : TYPE str
        DESCRIPTION: drug of the interest
    sign_gene_dict : TYPE dict
        DESCRIPTION: dictionary of drugs to significantly transcribed genes
    regnet_doc : TYPE str
        DESCRIPTION: path to the regulatory network file that includes
                    the list of interacting proteins where at least one of 
                    the pairs is transcriptioon factor 
    tfdoc : TYPE str
        DESCRIPTION: path to the DBTF file

    Returns
    -------
    tf_related_genelist: TYPE list
                DESCRIPTION: list of transcription factors related with significantly
                transcribed genes collected from transcriptomic data (L1000)
    interacting_genes: TYPE dict
                DESCRIPTION: a dictionary of TFs to the genes it regulates

    """
    regnet=pd.read_csv(regnet_doc, header=None, sep='\t')
    regnet.columns=['TF_symbol','TF_id','gene2_symbol','gene2_id']
    genelist=set(sign_gene_dict[drugname])
    print ('Number of significant genes for the drug {0}:\t{1}'.format(drugname,len(genelist)))
    tf_related_genelist=[]
    
    # find TFs that are already in the significantly transcribed gene list
    for g in genelist:
        related_tf=set(regnet[regnet['gene2_symbol']==g]['TF_symbol'].unique())
        if len(related_tf)>0:
            tf_related_genelist.append(related_tf)
    
    # update the list of TF sets as a list of TFs
    thelist=tf_related_genelist[0]
    for gset in tf_related_genelist:
        thelist=thelist.union(gset)
    finallist=get_my_tf(list(thelist),tfdoc)    
    tf_related_genelist=list(finallist)
    
    # Collect the regulated genes in significantly transcribed gene list
    # associated with the TF in order to use the mean of their log2FC 
    # later on prize designation
    interacting_genes={}
    for tf in tf_related_genelist:
        related_gene=set(regnet[regnet['TF_symbol']==tf]['gene2_symbol'].unique())
        related_gene=related_gene.intersection(genelist)
        interacting_genes[tf]=related_gene
        
    print('Number of TFs related with the significantly transcribed genes for the drug {0}:\t{1}'.format(drugname,len(tf_related_genelist)))
    return tf_related_genelist, interacting_genes

def find_sign_phosphosites(drugname, p_ftest_doc, p_fcdoc, analyte_doc, pval):
    """
    Parameters
    ----------
    drugname : TYPE str
        DESCRIPTION: drug of interest
    p_ftest_doc : TYPE str
        DESCRIPTION: path to the P100 Ftest pvalues file
    p_fcdoc : TYPE str
        DESCRIPTION: path to the P100 fold changes files
    analyte_doc : TYPE str
        DESCRIPTION: path to the P100 analyte info file
    pval : TYPE: float
        DESCRIPTION: pvalue threshold to evaluate significancy

    Returns
    -------
    fc_dic: TYPE dict
        DESCRIPTION: a dictionary of proteins to corresponding fold change values

    """
    ftest=pd.read_csv(p_ftest_doc)
    ftest.index=ftest[ftest.columns[0]]
    ftest=ftest.drop(ftest.columns[0], axis=1)
    columns=[c.split('_')[0] for c in ftest.columns]
    ftest.columns=columns

    f_drug=ftest[ftest[drugname]<pval][drugname]
    phosphosites=[x.strip() for x in list(f_drug.index)]
    
    pinfo=pd.read_csv(analyte_doc,sep='\t')
    pinfo.index=pinfo[pinfo.columns[0]]
    
    # get the protein gene symbol that has the significantly phosphorylated phosphosite
    sign_phospho_genes=[]
    for ph in phosphosites:
        xx = pinfo['pr_gene_symbol'][pinfo['pr_analyte_id']==ph].item()
        # if xx not in sign_phospho_genes:
        sign_phospho_genes.append(xx)
        
    print ('\nNumber of proteins that have significantly '\
           'phosphorylated sites for the drug {0}:\t{1}'.format(drugname,len(sign_phospho_genes)))

    # Collect the log2FC values of those significantly phosphorylated proteins
    fc_data=pd.read_csv(p_fcdoc,sep='\t')
    fc_data.index=fc_data[fc_data.columns[0]]
    fc_data=fc_data.drop(fc_data.columns[0], axis=1)
    fc_columns=[c.split('_')[1] for c in fc_data.columns]
    fc_data.columns=fc_columns

    fc_dic={}
    for i in range(len(phosphosites)):
        phospho=phosphosites[i]
        protein=sign_phospho_genes[i]
        fc_num=fc_data[drugname].loc[phospho]
        if protein not in fc_dic.keys():
            fc_dic[protein]=[fc_num]
        else:
            fc_dic[protein].append(fc_num)
    
    return fc_dic

"""
python 2.prepare_seed_protein_list.py -c NPC -d pazopanib -P 0.05 -L 0.05 -lf ./raw_data/NPC_L1000_Fresults_pvalues.csv -pf ./raw_data/NPC_P100_Fresults_pvalues.csv --l1000_fc ./raw_data/NPC_L1000_fold_changes_p_0.001.txt --p100_fc ./raw_data/NPC_P100_fold_changes_p_0.001.txt 
"""
def main():
    #Parsing arguments 
    # (run python 2.prepare_seed_protein_list.py -h to see all these decriptions)
    parser = argparse.ArgumentParser(
        description='Given the drugname, cell line name of interest and pvalues '\
        'for L1000 and P100 data, this script finds out seed protein list '\
        'including transcriptomic, phosphoproteomic and drug targetome data; '\
        'then calculates corresponding prize values and '\
        'outputs the prize file as "{cell_line}_{drugname}_prizefile.txt".')
    #required arguments
    parser.add_argument("-c", "--cell", dest='cellline', help='Cell line of interest '\
                        'Choices: "A375","A549","MCF7","PC3","YAPC","NPC". ')
    parser.add_argument("-d", "--drug", dest='drugname', help='Drug of interest')
    parser.add_argument("-P", "--pval_p100", dest='pval_p100', 
                        help='Pvalue necessary to collect significantly phophorylated '\
                            'proteins from P100 data.')
    parser.add_argument("-L", "--pval_l1000", dest='pval_l1000',
                        help='Pvalue necessary to collect TFs related with significantly '\
                            'transcribed genes from L1000 data.')
    
    parser.add_argument("-lf", "--LFfile", dest='LFtestfile', 
                        help='Path to the L1000 F test pvalues file.')
    parser.add_argument("-pf", "--PFfile", dest='PFtestfile', 
                        help='Path to the P100 F test pvalues file.')
    
    parser.add_argument("--l1000_fc", dest='l1000_fc', help='Path to the L1000 Fold Change file.')
    parser.add_argument("--p100_fc", dest='p100_fc', help='Path to the P100 Fold Change file.')
        
    #optional arguments
    parser.add_argument("-r", "--regnet", dest='regnet_file', help='Path to the Regulatory Network file '\
        'Default = "./raw_data/human_regulatory_network/human.source"', 
        default='./raw_data/human_regulatory_network/human.source')
        
    parser.add_argument("-t", "--tf_file", dest='tf_file', help='Path to the DBTF file. '\
        'Default = "./raw_data/DBTF_List.txt"', default='./raw_data/DBTF_List.txt')
    parser.add_argument("-dt", "--targetome_file", dest='targetome_file', help='Path to the drug targets file. '\
        'Default = "./raw_data/drug_targets.txt"', default='./raw_data/drug_targets.txt')
        
    parser.add_argument("-lg", "--Lgene", dest='LgeneFile', help='Path to the L1000 gene info file '\
        'Default = "./raw_data/GSE101406_Broad_LINCS_L1000_gene_info.txt"', 
        default='./raw_data/GSE101406_Broad_LINCS_L1000_gene_info.txt') 
    
    parser.add_argument("-pa", "--Panalyte", dest='PanalyteFile', help='Path to the P100 analyte info file '\
        'Default = "./raw_data/GSE101406_Broad_LINCS_P100_analyte_info.txt"', 
        default='./raw_data/GSE101406_Broad_LINCS_P100_analyte_info.txt')
        
    parser.add_argument("--outpath", dest = 'outputpath', help='Path to the directory which will '\
        'hold the output files. Default = "./raw_data/prize_files/"', default='./raw_data/prize_files/')
    
    options = parser.parse_args()
    #Check if outputpath exists
    if not os.path.isdir(options.outputpath):
        # sys.exit('Outpath %s is not a directory'%options.outputpath)
        os.mkdir(options.outputpath)

    pval_l1000=float(options.pval_l1000)
    pval_p100=float(options.pval_p100)
    drugname=options.drugname
    cell=options.cellline
    print('...\n...Started processing L1000 F test results and collecting significantly expressed genes '\
          'for the drug {0} in cell line {1}.\n'.format(drugname,cell))
    result = get_significant_genes(options.LFtestfile, options.LgeneFile, pval_l1000)
    print('Transcription Factors related with the significantly transcribed genes are being collected now...')
    tflist, gene_dic = find_related_tf(drugname, result, options.regnet_file, options.tf_file)

    #extract genes regulated by upstream transcription factors
    fcdata=pd.read_csv(options.l1000_fc, sep='\t')
    fcdata.index=fcdata[fcdata.columns[0]]
    fcdata=fcdata.drop(fcdata.columns[0], axis=1)
    columns= fcdata.columns
    columns=[c.split('_')[0] for c in columns]
    fcdata.columns=columns

    fc_dic={}
    for k,v in gene_dic.items():
        logfc=[]
        if len(v)>2:    ###--> in order to take TFs which regulates at least 3 genes in our list
            for gene in v:
                try:
                    fc=fcdata[drugname].loc[gene]
                    logfc.append(fc)
                except KeyError:
                    print('\nKeyError:\t{0}\nThe above gene cannot be found in fold change data.'.format(gene))
            fc_dic[k]=logfc

    print('...\n...Started processing P100 F test results and collecting significantly\nphosphorylated proteins '\
          'for the drug {0} in cell line {1}.'.format(drugname,cell))
    ##extract significantly phosphorylated proteins
    p100_dic = find_sign_phosphosites(drugname, options.PFtestfile, options.p100_fc, options.PanalyteFile, pval_p100)

    print('\nAll terminals are collected. Now prize values will be calculated.')
    prize_dic={}
    mxm=0.0
    for k,v in p100_dic.items():
        mean=np.mean(v)
        prize_dic[k]=mean
        if mean>mxm:
            mxm=mean
    for k,v in fc_dic.items():
        if k not in prize_dic.keys():
            mean=np.mean(v)
            prize_dic[k]=mean
            if mean>mxm:
                mxm=mean
        else:
            print('\n%s is in phosphorylated protein list,so the one coming from L1000 is ignored.'%k)

    
    drg_trgt=pd.read_csv(options.targetome_file,sep='\t')
    drug_targets=drg_trgt[drg_trgt['DRUGS']==drugname]['Target_Molecules'].tolist()[0]
    drug_targets=drug_targets.split(',')
    drug_targets=[drg.strip() for drg in drug_targets]
    if 'NA' in drug_targets:
        drug_targets=[]
    print ('\nTargets of drug {0}:\t{1}'.format(drugname, ','.join(drug_targets)))
    
    for t in drug_targets:
        prize_dic[t]=1.25000 #This is our selected constant drug target prize, can be changed optionally
        if mxm < 1.25000:
            mxm = 1.25000
    
    # Printing the summary of seed proteins: 
    print ('\n%s-p value<%s (L1000) and p-value<%s (P100):'%(drugname,options.pval_l1000,options.pval_p100))
    print ('L1000 length:', len(fc_dic))
    print ('P100 length:', len(p100_dic))
    print ('total length:', len(prize_dic))
    print ('Maximum prize weight is:',mxm)
    
    doc=options.outputpath+'%s_prize_file.txt'%drugname
    data=open(doc,'w')
    for k,v in prize_dic.items():
        word='%s\t%s\n'%(k, str(abs(v)))
        data.write(word)
    data.close()
    print ('\nPrize doc is written for {0} in cell line {1}!'.format(drugname, cell))
    
if __name__ == '__main__':
    main()
    
    
