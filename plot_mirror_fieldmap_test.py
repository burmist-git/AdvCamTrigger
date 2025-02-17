import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd

df=pd.read_csv('mirror_CTA-LST-v20141201-198_short.csv')

#%matplotlib widget 

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i in range(0,len(df)):
    v_s = np.array([df['x'].values[i],df['y'].values[i],df['z'].values[i]])
    v_dir = np.array([df['n_x'].values[i],df['n_y'].values[i],df['n_z'].values[i]])
    vlength=np.linalg.norm(v_dir)*30000
    #print(vlength)
    ax.quiver(v_s[0],v_s[1],v_s[2],v_dir[1],v_dir[0],v_dir[2],
            pivot='tail',length=vlength,arrow_length_ratio=0.3)
ax.set_xlim([-1000,1000])
ax.set_ylim([-1000,1000])
ax.set_zlim([0,1000])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
