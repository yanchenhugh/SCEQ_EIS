#!/usr/bin/env python
# coding: utf-8

# In[139]:


#!/usr/bin/env python
# coding: utf-8

import os
os.chdir(r'E:\my_model\IRR_code_v3')

import pandas as pd
import numpy as np

Ts = 61
Tstar = 60

df = pd.read_csv('model_sols_can_ghq7.csv')
#df = pd.read_csv('test_can.csv')
dfn = df.apply(pd.to_numeric,errors='coerce').dropna()
dt = dfn.to_numpy()
#print(dt[:,-3:])

err_cis = dt[:,-4:] # rgm succ succ_err err
#print(err_cis[0:10,])
dt_cis = dt[:,[1,3,4]] # y c inv

err2_cis = []
errv_cis = []
endo_cis = []
errcheck_cis = []      # Nsim x 2 matrix, [max mean]
errcheck_cis_bind = [] # Nsim x 3 matrix, [max sum #]
succ = 0 # count the number of successful simulations
for ji in range(1000):
    slic = err_cis[Ts*ji+1:Ts*ji+Tstar,]
    if min(slic[:,1]) == 1: # successful simulation and error checking
        succ = succ + 1
        #slic3 = slic[:,:,np.newaxis]
        #err2_cis = slic3 if ji < 1 else np.concatenate((err2_cis,slic3),axis=2)
        errv_cis = slic[:,-1] if ji < 1 else np.hstack((errv_cis,slic[:,-1]))
        
        slic_endo = dt_cis[Ts*ji+1:Ts*ji+Tstar,]
        slic3_endo = slic_endo[:,:,np.newaxis]
        endo_cis = slic3_endo if ji < 1 else np.concatenate((endo_cis,slic3_endo),axis=2)
        err_1 = np.ndarray((1,2))
        err_1[:,0] = np.max(slic[slic[:,2] == 1,-1])
        err_1[:,1] = np.mean(slic[slic[:,2] == 1,-1])
        errcheck_cis = err_1 if ji < 1 else np.concatenate((errcheck_cis,err_1),axis=0)
        
        err_2 = np.ndarray((1,3))
        err_2[:,2] = np.sum(slic[:,0]) 
        err_2[:,0] = np.max(slic[np.all(slic[:,[0,2]],axis=1) == 1,-1]) if err_2[:,2] > 0 else 0
        err_2[:,1] = np.sum(slic[np.all(slic[:,[0,2]],axis=1) == 1,-1]) if err_2[:,2] > 0 else 0
        errcheck_cis_bind = err_2 if ji < 1 else np.concatenate((errcheck_cis_bind,err_2),axis=0)

std_cis_all = np.std(endo_cis,axis=0)
std_cis = {'data': std_cis_all, 'mean': np.mean(std_cis_all,axis=1), 'std': np.std(std_cis_all,axis=1)}

print('Number of successful simulations: ' + str(succ))

errrep_cis = np.ndarray((1,3)) # [[mean of max mean of all],[max of bind mean of bind]]
errrep_cis[0,0] = np.max(errcheck_cis[:,0])
#errrep_cis[0,0] = np.max(errcheck_cis[:,1])
errrep_cis[0,1] = np.mean(errcheck_cis[:,1])
errrep_cis[0,2] = np.sum(errcheck_cis_bind[:,1]) / np.sum(errcheck_cis_bind[:,2])

print(errrep_cis)
#print(std_cis['mean'])

# Now deal with EIS

Ts = 60
Tstar = 60

#df = pd.read_csv('test_eis.csv')
df0 = pd.read_csv('model_sols_eis_ghq7_0.csv',skiprows=1) # 70
dfn0 = df0.apply(pd.to_numeric,errors='coerce').dropna()
dt0 = dfn0.to_numpy()

df1 = pd.read_csv('model_sols_eis_ghq7_1.csv',skiprows=1) # 461
dfn1 = df1.apply(pd.to_numeric,errors='coerce').dropna()
dt1 = dfn1.to_numpy()

df2 = pd.read_csv('model_sols_eis_ghq7_2.csv',skiprows=1) # 475
dfn2 = df2.apply(pd.to_numeric,errors='coerce').dropna()
dt2 = dfn2.to_numpy()

df3 = pd.read_csv('model_sols_eis_ghq7_3.csv',skiprows=1) # 475
dfn3 = df3.apply(pd.to_numeric,errors='coerce').dropna()
dt3 = dfn3.to_numpy()

dt = np.vstack((dt1,dt2,dt3))

err_eis = dt[:,-4:] # rgm succ succ_err err
dt_eis = dt[:,[1,3,4]] # y c inv

err2_eis = []
errv_eis = []
endo_eis = []
errcheck_eis = []      # Nsim x 2 matrix, [max mean]
errcheck_eis_bind = [] # Nsim x 3 matrix, [max sum #]
succ = 0 # count the number of successful simulations
for ji in range(1200):
    slic = err_eis[Ts*ji:Ts*ji+Tstar-1,]
    if min(slic[:,1]) == 1: # successful simulation and error checking
        succ = succ + 1
        #slic3 = slic[:,:,np.newaxis]
        #err2_eis = slic3 if ji < 1 else np.concatenate((err2_eis,slic3),axis=2)
        errv_eis = slic[:,-1] if ji < 1 else np.hstack((errv_eis,slic[:,-1]))
        
        slic_endo = dt_eis[Ts*ji:Ts*ji+Tstar-1,]
        slic3_endo = slic_endo[:,:,np.newaxis]
        endo_eis = slic3_endo if succ < 2 else np.concatenate((endo_eis,slic3_endo),axis=2)
        err_1 = np.ndarray((1,2))
        err_1[:,0] = np.max(slic[slic[:,2] == 1,-1])
        err_1[:,1] = np.mean(slic[slic[:,2] == 1,-1])
        errcheck_eis = err_1 if succ < 2 else np.concatenate((errcheck_eis,err_1),axis=0)
        
        err_2 = np.ndarray((1,3))
        err_2[:,2] = np.sum(slic[:,0]) 
        err_2[:,0] = np.max(slic[np.all(slic[:,[0,2]],axis=1) == 1,-1]) if err_2[:,2] > 0 else 0
        err_2[:,1] = np.sum(slic[np.all(slic[:,[0,2]],axis=1) == 1,-1]) if err_2[:,2] > 0 else 0
        errcheck_eis_bind = err_2 if succ < 2 else np.concatenate((errcheck_eis_bind,err_2),axis=0)

std_eis_all = np.std(endo_eis,axis=0)
std_eis = {'data': std_eis_all, 'mean': np.mean(std_eis_all,axis=1), 'std': np.std(std_eis_all,axis=1)}

#err_eis = 
print('Number of successful simulations: ' + str(succ))
#print(errcheck_eis)
#print(errcheck_eis_bind)
#print(std_eis['mean'])

errrep_eis = np.ndarray((1,3)) # [[mean of max mean of all],[max of bind mean of bind]]
errrep_eis[0,0] = np.max(errcheck_eis[0:1000,0])
errrep_eis[0,1] = np.mean(errcheck_eis[0:1000,1])
errrep_eis[0,2] = np.sum(errcheck_eis_bind[0:1000,1]) / np.sum(errcheck_eis_bind[0:1000,2])

# combine cis and eis into one table and print out
errrep = np.vstack((errrep_cis,errrep_eis))
infos = ["CIS","EIS"]
errmatr = np.vstack( (infos,np.transpose(errrep)) )
err_matriks = np.vstack( ( [" ","Linf","L1","L1 binding"], np.transpose(errmatr)) )
print(err_matriks)

# write to csv
import csv
err_file = open(r".\euler_err.txt",'w+')
with err_file:
    draft = csv.writer(err_file)
    draft.writerows(err_matriks)


# In[ ]:




