from scipy.io import loadmat
import numpy as np
import pickle as pkl
import pandas as pd
import time
import subprocess
import sys
#
#df = pd.read_csv('data_non_normalized_36575ev_887ch.csv', sep=' ', header=None)
df = pd.read_csv('data_non_normalized.csv', sep=' ', header=None)
U, S, Vh = np.linalg.svd(df.values, full_matrices=True)
print(df.shape)
print(U.shape)
print(S.shape)
print(Vh.shape)
#
#df_U=pd.DataFrame(U)
#df_U.to_csv('U.cvs',sep=' ',header=False,index=False)
#
df_S=pd.DataFrame(S)
df_S.to_csv('S.cvs',sep=' ',header=False,index=False)
#
df_Vh=pd.DataFrame(Vh)
df_Vh.to_csv('Vh.cvs',sep=' ',header=False,index=False)
#
n=45
reco = np.matrix(U[:, :n]) * np.diag(S[:n]) * np.matrix(Vh[:n, :])
df_reco=pd.DataFrame(reco)
df_reco[:100].to_csv('reco.cvs',sep=' ',header=False,index=False)
#
#start, end, step = 10, 800, 10
#for i in range(start, end, step):
#    reco = np.matrix(U[:, :i]) * np.diag(S[:i]) * np.matrix(Vh[:i, :])
#    df_reco=pd.DataFrame(reco)
#    df_reco[:100].to_csv(str('reco_' + str(i) + '.cvs'),sep=' ',header=False,index=False)
#    print(i)
#    #print(reco.shape)
