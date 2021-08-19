
Integrative Network Modeling by use of Omics Integrator
--------------------------------------------------------

Omics Integrator is a package designed to integrate proteomic data, gene expression data and/or epigenetic data using a protein-protein interaction network. It is comprised of two modules, Garnet and Forest.

Contact: Seyma Unsal Beyge [unsalseyma@gmail.com]

All data included in P100 dataset is provided in "./raw_data/" folder.

The protocol of network modeling is as below:
---------------------------------------------
**1.** prepare_raw_data_files.py

**2.** prepare_seed_protein_list.py

**3.** prepare_interactome.py 

**4.** find_all_possible_drug_networks.py

**Running prepare_raw_data_files.py**

Find significantly perturbed genes and/or phosphoproteins and create the files:
- p values of each gene/phosphoprotein in each significantly perturbed condition 
- significant gene read outs
- log2FoldChanges of those significantly perturbed gene/protein list
	
If only landmark genes are chosen, all files will be prepared for 978 landmark genes.

**Options:**

-h, --help  show this help message and exit

-c, --cell			cellline
(Required) Cell line of interest. Choices: "A375","A549","MCF7","PC3","YAPC","NPC".

-p, --pval			pvalThreshold
p-value threshold that is desired to get the significant conditions (Default is 0.05)

-L, --L1000			l1000File

						Path to the raw L1000 data file
						Default = "./raw_data/GSE101406_Broad_LINCS_L1000_Level3_INF_mlr12k_n1667x12328.gctx"

-P, --P100			P100File

						Path to the raw P100 data file
						Default = "./raw_data/GSE101406_Broad_LINCS_L1000_Level3_INF_mlr12k_n1667x12328.gctx"

-ci, --cellinfo		cellFile

						Path to the cell info file
						Default = "./raw_data/GSE101406_Broad_LINCS_cell_info.txt"

-pt, --pert			pertFile

						Path to the perturbant info file
						Default = "./raw_data/GSE101406_Broad_LINCS_pert_info.txt"

-lg, --Lgene			LgeneFile

						Path to the L1000 gene info file
						Default = "./raw_data/GSE101406_Broad_LINCS_L1000_gene_info.txt"

-la, --Lanalyte		LanalyteFile

						Path to the L1000 analyte info file
						Default = "./raw_data/GSE101406_Broad_LINCS_L1000_inst_info.txt"

-pi, --PinstFile		PinstFile

						Path to the P100 instrument info file
						Default = "./raw_data/GSE101406_Broad_LINCS_P100_inst_info.txt"

-pa, --Panalyte		PanalyteFile

						Path to the P100 analyte info file
						Default = "./raw_data/GSE101406_Broad_LINCS_P100_analyte_info.txt"

-lm, --lmGenes		lmGenesFile

						Path to the L1000 Landmark genes file
						Default = "./raw_data/landmark_genes_978.txt"

--iflandmark			iflandmark

						True or False to select which genes will be used in transcriptomic data (L1000).
						If True, only landmark genes will be used in the analysis. Landmark genes are
                        listed in the file whose path can be given with --lmGenesoption.
						Default = False

--outpath   outpath

            Path to the directory which will hold the output files
						Default = "./raw_data/"

**Running prepare_seed_protein_list.py**
Given the drugname, cell line name of interest and pvalues for L1000 and P100 data, this script finds out seed protein list including transcriptomic, phosphoproteomic and drug targetome data; then calculates corresponding prize values and 
outputs the prize file as "{cell_line}_{drugname}_prizefile.txt".

Options:
  -h, --help            show this help message and exit
  -c, --cell			cellline
						(Required) Cell line of interest
                        Choices: "A375","A549","MCF7","PC3","YAPC","NPC".
  -d, --drug			drugname
						(Required) Drug of interest
  -P, --pval_p100		pval_p100
						(Required) Pvalue necessary to collect significantly phophorylated proteins from P100 data.
  -L, --pval_l1000		pval_l1000
						(Required) Pvalue necessary to collect TFs related with significantly transcribed genes from L1000 data.
  -lf, --LFfile			LFtestfile
						(Required) Path to the L1000 F test pvalues file.
  -pf, --PFfile			PFtesetfile
						(Required) Path to the P100 F test pvalues file.
  --l1000_fc			l1000_fc
						(Required) Path to the L1000 Fold Change file.
  --p100_fc				p100_fc
						(Required) Path to the P100 Fold Change file.
  -r, --regnet			regnet_file
						Path to the Regulatory Network file
						Default = "./raw_data/human_regulatory_network/human.source"
  -t, --tf_file			tf_file
						Path to the DBTF file.
						Default = "./raw_data/DBTF_List.txt"
  -dt, --targetome_file	targetome_file
						Path to the drug targets file.
						Default = "./raw_data/drug_targets.txt"
  -lg, --Lgene			LgeneFile
						Path to the L1000 gene info file
						Default = "./raw_data/GSE101406_Broad_LINCS_L1000_gene_info.txt"
  -ps, --Panalyte		PanalyteFile
						Path to the P100 analyte info file
						Default = "./raw_data/GSE101406_Broad_LINCS_P100_analyte_info.txt"
  --outpath				outpath
						Path to the directory which will hold the output files.
						Default = "./raw_data/prize_files/"

**Running prepare_interactome.py**
Given the drugname, cell line name of interest, path to the reference interactome paths to the raw L1000 data, gene-inst-pert info files, path to the file including list of landmark genes processes the reference interactome to exclude hub nodes, low-expressed genes, then perfroms link prediction and adds edges passing localization filters. Finally outputs the interactome file as "{cell_line}_{drugname}_iref_processed_interactome.txt".

Options:
  -c, --cell			cellline
						(Required) Cell line of interest
                        Choices: "A375","A549","MCF7","PC3","YAPC","NPC".
  -d, --drug			drugname
						(Required) Drug of interest
  -i, --int				interactome
						Path to the reference interactome
						default='./raw_data/iref_mitab_miscore_2013_08_12_interactome_no_selfloop.txt
  -lg, --Lgene			LgeneFile
						Path to the L1000 gene info file
						Default = "./raw_data/GSE101406_Broad_LINCS_L1000_gene_info.txt"
  -la, --Lanalyte		LanalyteFile
						Path to the L1000 analyte info file
						Default = "./raw_data/GSE101406_Broad_LINCS_L1000_inst_info.txt"
  -L, --L1000			l1000File
						Path to the raw L1000 data file
						Default = "./raw_data/GSE101406_Broad_LINCS_L1000_Level3_INF_mlr12k_n1667x12328.gctx"
  --localization		localization
						Path to the localization info file.
						Default = "./raw_data/subcellular_location.txt"
  --expthreshold		expthreshold
						he level of expression value that will used for excluding the low expressed genes.
						Default = 2.0.
  --outpath				outputpath
						Path to the directory which will hold the output files.
						Default = "./raw_data/interactomes/{drugname}/"
  --lpmethod			lpmethod
						The link prediction method to use for interactome preparation.
						Options are "Adamic/Adar", "Jaccards", "Preferential Attachment", "Resource Allocation" and None. 
						If None, any link prediction is not used and the preparation procedure is terminated.
						default='Adamic/Adar'
    
**Running find_all_possible_drug_networks.py**
This script designed as a case study which runs Forest module of Omics Integrator.
It includes several variables such as paths to the prize and edge files, parameters needed by Forest, outpaths etc.
These variables within the script should be updated accordingly.





