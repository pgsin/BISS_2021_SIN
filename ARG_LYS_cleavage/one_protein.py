import Bio
from Bio import *
from Bio.PDB import *
from Bio import SeqIO
from Bio import Struct
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import nglview as nv
import pytraj as pt
from pytraj import view
import requests

import sys

NAME = str(sys.argv[1])

url = 'https://alphafold.ebi.ac.uk/files/AF-'+NAME+'-F1-model_v1.pdb'
print(NAME, url)
r = requests.get(url)

if len(r.content)<200:
	print('Bye')
	exit()
	
with open('AF-'+NAME+'-F1-model_v1.pdb', 'wb') as f:
    f.write(r.content)
    
for seq_record in SeqIO.parse("UP000005640_9606.fasta", "fasta"):
    if seq_record.id[:10] == "sp|"+NAME+"|":
        #print(seq_record.id)
        P = str(seq_record.seq)
        #print(P)
        #print(len(P))
        break
        
df = pd.read_csv('peptides.txt', sep="\t", lineterminator='\r')
peptides = list()
for i in range(len(df["Sequence"])):
    if df["Leading razor protein"].iloc[i] == NAME:
        peptides.append(df["Sequence"].iloc[i][1:])
            
ff = open("conv_"+NAME, "w")
ff.write('color grey, all \n')

cleavage = list()
KR = ['K', 'R']
for pep in peptides:
    i = P.find(pep)
    #print(P[i], P[i-1])
    if P[i] in KR:
        #print('qq')
        cleavage.append(i)
    elif P[i-1] in KR:
        #print('qq')
        cleavage.append(i-1)
    else:
        cleavage.append(i+1)
    ff.write('color red, resi ' + str(i) + "-" + str(i+len(pep)) + "\n")
    #print(i, cleavage[-1])
    #print(pep, P[i-1:i+len(pep)+1])
ff.close()



cleavage = sorted(cleavage)
#print([P[c] for c in cleavage])
cleavage = np.array(cleavage)+1
#print(cleavage)


def dist(r1, r2):
    return ( (r1[0]-r2[0])**2 + (r1[1]-r2[1])**2 + (r1[2]-r2[2])**2 )**(1/2)
    
cm_cleav = list()
cm_not = list()

f = open('forpym'+NAME, 'w')
f.write('color grey, all \n')

p=PDBParser(PERMISSIVE=1)
#s = p.get_structure("P02787", "AF-P02787-F1-model_v1.pdb")
s = Struct.read('AF-'+NAME+'-F1-model_v1.pdb')
cm = s.center_of_mass()
#print(cm)
arg_lys_cm_cl = list()
arg_lys_cm_not = list()
AL = ['ARG', 'LYS']
#prev = list(s.get_residues())[0]

for aa in s.get_residues():
    if aa.get_resname() in AL:
        if aa.get_full_id()[3][1] not in cleavage:
            arg_lys_cm_not.append(dist(cm, aa.center_of_mass()))
            f.write('color blue, resi '+str(aa.get_full_id()[3][1]) + "\n")
            cm_not.append(aa.center_of_mass())
    
    if aa.get_full_id()[3][1] in cleavage:
        #if aa.get_resname() not in AL:
            #print("!!!")
        #print(aa.get_full_id()[3][1], aa.get_resname(), aa.center_of_mass(), dist(cm, aa.center_of_mass()))
        arg_lys_cm_cl.append(dist(cm, aa.center_of_mass()))
        f.write('color green, resi '+str(aa.get_full_id()[3][1]) + "\n")
        cm_cleav.append(aa.center_of_mass())
    #prev = aa
f.close()


fig, ax = plt.subplots()
sns.histplot(arg_lys_cm_not, bins = 20, color = "blue", label = "not cleav", alpha = 0.5, kde = True, edgecolor=None)
sns.histplot(arg_lys_cm_cl,  bins = 20,  color = "green", label = "cleav", alpha = 0.5, kde=True, edgecolor=None)
ax.set_xlabel("r, A")
ax.set_ylabel("frequency")
ax.set_title("distance from proteinCM to ARG/LYS")
ax.legend()
plt.savefig("plot/"+NAME+"cm.png")

import os
os.system('pymol create_mesh.pml -- %s >nohup.out' % NAME)

import trimesh
mesh = trimesh.load(str(NAME+'_mesh.obj'))
#mesh.show()
(closest_points,distances_cleav,triangle_id) = mesh.nearest.on_surface(np.array(cm_cleav))
(closest_points,distances_not,triangle_id) = mesh.nearest.on_surface(np.array(cm_not))

#print(distances_cleav)
#print(distances_not)
fig, ax = plt.subplots()
sns.histplot(distances_not, bins = 30, color = "blue", label = "not cleav", alpha = 0.5, kde = True, edgecolor=None)
sns.histplot(distances_cleav,  bins = 30,  color = "green", label = "cleav", alpha = 0.5, kde=True, edgecolor=None)
ax.set_xlabel("d, A")
ax.set_ylabel("frequency")
ax.set_title("distance from surface to ARG/LYS")
ax.legend()
plt.savefig("plot/"+NAME+"surf.png")

count_not, bins_count_not = np.histogram(distances_not, bins=30)
count_cleav, bins_count_cleav = np.histogram(distances_cleav, bins=30)

pdf_not = count_not / sum(count_not)
pdf_cleav = count_cleav / sum(count_cleav)

cdf_not = np.cumsum(pdf_not)
cdf_not = cdf_not / cdf_not[-1]
cdf_cleav = np.cumsum(pdf_cleav)
cdf_cleav = cdf_cleav/cdf_cleav[-1]

fig, ax = plt.subplots()
plt.plot(bins_count_not[1:], pdf_not, color="blue", label="not cleav")
plt.plot(bins_count_cleav[1:], pdf_cleav, color="green", label="cleav")
plt.title('PDF')
plt.xlabel('d, A')
plt.legend()
plt.savefig("plot/"+NAME+"pdf.png")

fig, ax = plt.subplots()
plt.plot(bins_count_not[1:], cdf_not, color="blue", label="not cleav")
plt.plot(bins_count_cleav[1:], cdf_cleav, color="green", label="cleav")
plt.title('CDF')
plt.xlabel('d, A')
plt.legend()
plt.savefig("plot/"+NAME+"sdf.png")

from scipy.stats import kstest
from Bio.SeqUtils.ProtParam import ProteinAnalysis

test_stat = kstest(distances_not, distances_cleav)

mass = ProteinAnalysis(P).molecular_weight()/1000

with open('out.txt', "a") as fff:
    fff.write(str(NAME) + " " + str(mass) + " " + str(test_stat.pvalue) + "\n")

