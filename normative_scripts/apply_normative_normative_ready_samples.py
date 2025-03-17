#!/usr/bin/env python
# coding: utf-8

import sys
# In[1]:
if len(sys.argv) != 5:
	print("Usage:", sys.argv[0], "<model_name>", "<site_names>", "<data_dir>", "<phenotypes>")
	exit(1)
else:
	model_name = sys.argv[1]	
	site_names = sys.argv[2]	
	datadir = sys.argv[3]
	phenotypes = sys.argv[4]

import os
import shutil

def add_folder_path(folder_path, file_names):
    return [os.path.join(folder_path, file_name) for file_name in file_names]

def remove_files_not_in_array(folder_path, file_array):
    file_array = add_folder_path(folder_path, file_array)
    for root, dirs, files in os.walk(folder_path):
        # Check files
        for file in files:
            file_path = os.path.join(folder_path, file)
            if file_path not in file_array:
                os.remove(file_path)
                print(f"Removed file: {file_path}")
        # Check directories
        for directory in dirs:
            dir_path = os.path.join(folder_path, directory)
            if dir_path not in file_array:
                shutil.rmtree(dir_path)
                print(f"Removed directory: {dir_path}")


allowed_files = [
    'resp_te.txt',
    'Z_predict.txt',
    'ys2_predict.txt',
    'yhat_predict.txt'
]

# In[2]:
#model_name = 'lifespan_57K_82sites'
#site_names = 'site_ids_ct_82sites.txt'
#phenotypes = 'phenotypes_sc_cth_test.txt'
#datadir = '/data_normative_long/df_mri/all/main-normsamples'


wd="/ess/p274/cluster/projects/p039_image_brain_change"
masterdir = wd + "/scripts/braincharts-master/"
normativedir = wd +  datadir

# where the analysis takes place
root_dir = masterdir

# where the data files live
data_dir = masterdir + '/docs'

# where the models live
out_dir = os.path.join(root_dir, 'models', model_name)
savedir = os.path.join(normativedir, model_name)
if os.path.exists(savedir):
    shutil.rmtree(savedir)

shutil.copytree(out_dir, savedir)
out_dir = savedir
# where the data will be saved


# scripts folder when we import the libraries in the code block below
#sys.path.insert(1, masterdir)
os.chdir(masterdir + 'scripts/') #this path is setup for running on Google Colab. Change it to match your local path
print(os.getcwd())
sys.path.insert(1, masterdir + 'scripts/')
# In[3]:


import numpy as np
import pandas as pd
import pickle
from matplotlib import pyplot as plt
import seaborn as sns

from pcntoolkit.normative import estimate, predict, evaluate
from pcntoolkit.util.utils import compute_MSLL, create_design_matrix
from nm_utils import remove_bad_subjects, load_2d


# In[4]:


# load a set of site ids from this model. This must match the training data
with open(os.path.join(root_dir,'docs', site_names)) as f:
    site_ids_tr = f.read().splitlines()


# In[7]:


test_data = os.path.join(normativedir, 'df_te.csv')
df_te = pd.read_csv(test_data)

# extract a list of unique site ids from the test set
site_ids_te =  sorted(set(df_te['site'].to_list()))


# In[8]:


# load the list of idps for left and right hemispheres, plus subcortical regions
with open(os.path.join(data_dir, phenotypes)) as f:
    idp_ids = f.read().splitlines()


# In[9]:


# which data columns do we wish to use as covariates? 
cols_cov = ['age','sex']

# limits for cubic B-spline basis 
xmin = 10 
xmax = 105

# Absolute Z treshold above which a sample is considered to be an outlier (without fitting any model)
outlier_thresh = 7


# In[10]:


for idp_num, idp in enumerate(idp_ids): 
    print('Running IDP', idp_num, idp, ':')
    idp_dir = os.path.join(out_dir, idp)
    os.chdir(idp_dir)
    
    # extract and save the response variables for the test set
    y_te = df_te[idp].to_numpy()
    
    # save the variables
    resp_file_te = os.path.join(idp_dir, 'resp_te.txt') 
    np.savetxt(resp_file_te, y_te)
        
    # configure and save the design matrix
    cov_file_te = os.path.join(idp_dir, 'cov_bspline_te.txt')
    X_te = create_design_matrix(df_te[cols_cov], 
                                site_ids = df_te['site'],
                                all_sites = site_ids_tr,
                                basis = 'bspline', 
                                xmin = xmin, 
                                xmax = xmax)
    np.savetxt(cov_file_te, X_te)
    
    # check whether all sites in the test set are represented in the training set
    if all(elem in site_ids_tr for elem in site_ids_te):
        print('All sites are present in the training data')
        
        # just make predictions
        yhat_te, s2_te, Z = predict(cov_file_te, 
                                    alg='blr', 
                                    respfile=resp_file_te, 
                                    model_path=os.path.join(idp_dir,'Models'))
    else:
        print('Some sites missing from the training data. Adapting model')
        
        # save the covariates for the adaptation data
        X_ad = create_design_matrix(df_ad[cols_cov], 
                                    site_ids = df_ad['site'],
                                    all_sites = site_ids_tr,
                                    basis = 'bspline', 
                                    xmin = xmin, 
                                    xmax = xmax)
        cov_file_ad = os.path.join(idp_dir, 'cov_bspline_ad.txt')          
        np.savetxt(cov_file_ad, X_ad)
        
        # save the responses for the adaptation data
        resp_file_ad = os.path.join(idp_dir, 'resp_ad.txt') 
        y_ad = df_ad[idp].to_numpy()
        np.savetxt(resp_file_ad, y_ad)
       
        # save the site ids for the adaptation data
        sitenum_file_ad = os.path.join(idp_dir, 'sitenum_ad.txt') 
        site_num_ad = df_ad['sitenum'].to_numpy(dtype=int)
        np.savetxt(sitenum_file_ad, site_num_ad)
        
        # save the site ids for the test data 
        sitenum_file_te = os.path.join(idp_dir, 'sitenum_te.txt')
        site_num_te = df_te['sitenum'].to_numpy(dtype=int)
        np.savetxt(sitenum_file_te, site_num_te)
         
        yhat_te, s2_te, Z = predict(cov_file_te, 
                                    alg = 'blr', 
                                    respfile = resp_file_te, 
                                    model_path = os.path.join(idp_dir,'Models'),
                                    adaptrespfile = resp_file_ad,
                                    adaptcovfile = cov_file_ad,
                                    adaptvargroupfile = sitenum_file_ad,
                                    testvargroupfile = sitenum_file_te)

    remove_files_not_in_array(idp_dir, allowed_files)
    
os.chdir(out_dir)

