"""
script to run all parameter sets to select the best optimal 
parameters for network reconstruction with forest.py

Variables below should be updated accordingly.
"""

import os
import numpy as np

def conf_prep(mu,beta,D,w):
    file=open('conf.txt', 'w')
    file.writelines("w = %d\nb = %d\nD = %d\nmu = %f" % (w,beta,D,mu))
    file.close()

# Update this with the path to msgsteiner
msgpath = "/home/user/OmicsIntegrator/msgsteiner-1.3/msgsteiner"

# The parameter intervals, can be updated with selected ranges
mu_range = np.arange(0.000,0.101,0.005)
beta_range = np.arange(2,10.1,1)
w_range = np.arange(1,3.1,1)
D = 10

# cell line and drugname shoudl be specified!!
cell='cellline'
drug='drugname'

# Paths below should be updated !
prize_file = "./path_to_the_prize_file_of_selected_cell_line_drug_condition"
edge_file = "./path_to_the_processed_interactome_file_of_selected_cell_line_drug_condition"
conf_file = "conf.txt"
networks_path = "./results/{0}_{1}/".fomat(cell,drug)

# Create output directories if needed
if not os.path.exists(networks_path):
    os.makedirs(networks_path)

for mu in mu_range:
    for beta in beta_range:
        for w in w_range:
            conf_prep(mu,beta,D,w)
            out_label = "w%f_beta%d_D%d_mu%f" %(w,beta,D,mu)
            networks_path_new = networks_path+out_label+'/'
            if not os.path.exists(networks_path_new):
                os.makedirs(networks_path_new)
            
            # path to the forest.py should be updated !
            os.system("python /home/user/OmicsIntegrator/scripts/forest.py --prize %s --edge %s --conf conf.txt --msgpath %s --outpath %s --outlabel %s" %(prize_file,edge_file,msgpath,networks_path_new,out_label))
            