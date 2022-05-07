#!/usr/bin/env python
# coding: utf-8

# ### Load the data from Excel

# In[1]:


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pingouin as pg
from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox

f = pd.ExcelFile( "../data/Cumulus Expansion dataset_without_negative_CE_values.xlsx")                                                  
print(f.sheet_names)
dfs = {}
for sheet_name in f.sheet_names:
    dfs[sheet_name] = pd.read_excel(f, sheet_name, skiprows=[0,1], header=None,usecols=[1,2,3,4,5,6])


# In[3]:


dfs['Area method']


# In[4]:


def get_q(v):
    return np.log((v+100)/100)

def get_matrix(xls_df, judges=3, repeats=2, do_log=True):
    objects = len(xls_df.index)
    x = np.zeros((objects, judges, repeats))
    x[:,0,0] = xls_df[1]
    x[:,0,1] = xls_df[4]
    x[:,1,0] = xls_df[2]
    x[:,1,1] = xls_df[5]
    x[:,2,0] = xls_df[3]
    x[:,2,1] = xls_df[6]
    if do_log:
        x = get_q(x)
    return x


# In[215]:


def compute_kappas(method):
    x = get_matrix(dfs[method],do_log=False)
    N, J, R = x.shape
    judges = np.zeros(N*J)
    for i in range(J):
        judges[N*i:N*(i+1)]=i
    indexes = np.tile(np.arange(N),J)
    df_inter = pd.DataFrame({'oocyte':indexes,'measurer':judges})
    CI_low_inter = np.zeros(R)
    CI_high_inter = np.zeros(R)
    ICC_inter = np.zeros(R)
    for i in range(R):
        col_name = "r"+str(i)
        df_inter[col_name] = x[:,:,i].swapaxes(0,1).flatten()
        icc = pg.intraclass_corr(data=df_inter, targets='oocyte', raters='measurer', ratings=col_name, nan_policy="omit")
        icc = icc.set_index('Type')
        print("Inter expert ICCs for method " + method + " and measurement repetition "+str(i))
        print("ICC:",icc.loc["ICC2"]["ICC"]," CI95%:",icc.loc["ICC2"]["CI95%"])
        ICC_inter[i]=icc.loc["ICC2"]["ICC"]
        CI_low_inter[i] = float(icc.loc["ICC2"]["CI95%"][0])
        CI_high_inter[i] = float(icc.loc["ICC2"]["CI95%"][1])
    reps = np.zeros(N*R)
    for i in range(R):
        reps[N*i:N*(i+1)]=i
    indexesr = np.tile(np.arange(N),R)
    df_intra = pd.DataFrame({'oocyte':indexesr,'rep':reps})
    CI_low_intra = np.zeros(J)
    CI_high_intra = np.zeros(J)
    ICC_intra = np.zeros(J)
    for i in range(J):
        col_name = "a"+str(i)
        df_intra[col_name] = x[:,i,:].swapaxes(0,1).flatten()
        icc = pg.intraclass_corr(data=df_intra, targets='oocyte', raters='rep', ratings=col_name, nan_policy="omit")
        icc = icc.set_index('Type')
        print("Intra expert ICCs for method " + method + " and observer "+str(i))
        print("ICC:",icc.loc["ICC1"]["ICC"]," CI95%:",icc.loc["ICC1"]["CI95%"])
        ICC_intra[i]=icc.loc["ICC1"]["ICC"]
        CI_low_intra[i] = float(icc.loc["ICC1"]["CI95%"][0])
        CI_high_intra[i] = float(icc.loc["ICC1"]["CI95%"][1])
    return ICC_inter, CI_low_inter, CI_high_inter, ICC_intra, CI_low_intra, CI_high_intra
    


# In[216]:


iccs = {}
iccs["area"] = compute_kappas("Area method")
iccs["3 distance"] = compute_kappas('3distance method')
iccs["score"] = compute_kappas('Score method')


# In[228]:


x = ["Area", "3 distance", "Score"]

y = [np.mean(iccs["area"][0]), np.mean(iccs["3 distance"][0]), np.mean(iccs["score"][0])]
lowers = [y[0] - np.min(iccs["area"][1]), y[1] - np.min(iccs["3 distance"][1]), y[2]-np.min(iccs["score"][1])]
uppers = [np.max(iccs["area"][2])-y[0], np.max(iccs["3 distance"][2])-y[1], np.max(iccs["score"][2])-y[2]]
print(y)
print(lowers)
print(uppers)
yerr = [lowers,uppers]

plt.rcParams["figure.figsize"] = [16*0.4,9*0.4]
# plot:
fig, ax = plt.subplots()
#fig.figsize=(3,4)
ax.set_ylabel("ICC")
ax.errorbar(x, y, yerr, fmt='o', linewidth=3, capsize=6)
#ax.set_title("Inter expert agreement for the Area, 3 distance and Score methods")
ax.set(xlim=(-0.5, 2.5), xticks=np.arange(0, 3),
       ylim=(0, 1), yticks=np.arange(0,1,0.2 ))
r1=Rectangle((-0.5, 0.8), 3, 0.2, color="red", alpha=0.2)
r2=Rectangle((-0.5, 0.6), 3, 0.2, color="orange", alpha=0.2)
r3=Rectangle((-0.5, 0.4), 3, 0.2, color="yellow", alpha=0.2)
r4=Rectangle((-0.5, 0.2), 3, 0.2, color="cyan", alpha=0.2)
r5=Rectangle((-0.5, 0.), 3, 0.2, color="blue", alpha=0.2)
ax.add_patch(r1)
ax.add_patch(r2)
ax.add_patch(r3)
ax.add_patch(r4)
ax.add_patch(r5)
ax.legend([r1, r2, r3,r4,r5], ['Very good', 'Good', 'Moderate','Fair','Poor'], title="Levels of agreement", framealpha=1, bbox_to_anchor=(1., 1.))
plt.savefig("../figs/inter_expert_agreement.png",dpi=300.,bbox_inches=Bbox.from_extents(-0.5,-1,8,5))
#plt.show()


# In[234]:


x = [-0.1,0.,0.1,0.9,1.,1.1,1.9,2,2.1]
y = np.hstack((iccs["area"][3], iccs["3 distance"][3], iccs["score"][3]))
print(y)


lowers = np.hstack((iccs["area"][4], iccs["3 distance"][4], iccs["score"][4]))
lowers = y-lowers
uppers = np.hstack((iccs["area"][5], iccs["3 distance"][5], iccs["score"][5]))
uppers = uppers-y
print(lowers)
yerr = [lowers,uppers]

plt.rcParams["figure.figsize"] = [16*0.4,9*0.4]
# plot:
fig, ax = plt.subplots()
#fig.figsize=(3,4)
ax.set_ylabel("ICC")
ax.errorbar(x, y, yerr, fmt='o', linewidth=3, capsize=6)
#ax.set_title("Intra expert agreement for each of the three annotators\n when using the Area, 3 distance and Score methods")
ax.set(xlim=(-0.5, 2.5), ylim=(0, 1), yticks=np.arange(0,1,0.2 ))
ax.set_xticks(np.arange(0, 3),["Area","3 distance","Score"])
r1=Rectangle((-0.5, 0.8), 3, 0.2, color="red", alpha=0.2)
r2=Rectangle((-0.5, 0.6), 3, 0.2, color="orange", alpha=0.2)
r3=Rectangle((-0.5, 0.4), 3, 0.2, color="yellow", alpha=0.2)
r4=Rectangle((-0.5, 0.2), 3, 0.2, color="cyan", alpha=0.2)
r5=Rectangle((-0.5, 0.), 3, 0.2, color="blue", alpha=0.2)
ax.add_patch(r1)
ax.add_patch(r2)
ax.add_patch(r3)
ax.add_patch(r4)
ax.add_patch(r5)
ax.legend([r1, r2, r3,r4,r5], ['Very good', 'Good', 'Moderate','Fair','Poor'], title="Levels of agreement", framealpha=1, bbox_to_anchor=(1., 1.))
plt.savefig("../figs/intra_expert_agreement.png",dpi=300.,bbox_inches=Bbox.from_extents(-0.5,-1,8,5))
#plt.show()


# In[233]:


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox

#plt.style.use('_mpl-gallery')

# make data:
x = [-0.1,0.1,0.9,1.1,1.9,2.1]
#x = ["Area", "3 distance", "Score"]
y = np.hstack((iccs["area"][0], iccs["3 distance"][0], iccs["score"][0]))

lowers = np.hstack((iccs["area"][1], iccs["3 distance"][1], iccs["score"][1]))
lowers = y-lowers
uppers = np.hstack((iccs["area"][2], iccs["3 distance"][2], iccs["score"][2]))
uppers = uppers-y

yerr = [lowers,uppers]

plt.rcParams["figure.figsize"] = [16*0.4,9*0.4]
# plot:
fig, ax = plt.subplots()
#fig.figsize=(3,4)
ax.set_ylabel("ICC")
ax.errorbar(x, y, yerr, fmt='o', linewidth=3, capsize=6)
#ax.set_title("Inter expert agreement for the Area, 3 distance and Score methods")
r1=Rectangle((-0.5, 0.8), 3, 0.2, color="red", alpha=0.2)
r2=Rectangle((-0.5, 0.6), 3, 0.2, color="orange", alpha=0.2)
r3=Rectangle((-0.5, 0.4), 3, 0.2, color="yellow", alpha=0.2)
r4=Rectangle((-0.5, 0.2), 3, 0.2, color="cyan", alpha=0.2)
r5=Rectangle((-0.5, 0.), 3, 0.2, color="blue", alpha=0.2)
ax.set(xlim=(-0.5, 2.5), ylim=(0, 1), yticks=np.arange(0,1,0.2 ))
ax.set_xticks(np.arange(0, 3),["Area","3 distance","Score"])
ax.add_patch(r1)
ax.add_patch(r2)
ax.add_patch(r3)
ax.add_patch(r4)
ax.add_patch(r5)
ax.legend([r1, r2, r3,r4,r5], ['Very good', 'Good', 'Moderate','Fair','Poor'], title="Levels of agreement", framealpha=1, bbox_to_anchor=(1., 1.))
plt.savefig("../figs/inter_expert_agreement2.png",dpi=300.,bbox_inches=Bbox.from_extents(-0.5,-1,8,5))
#plt.show()


# In[ ]:




