
import csv
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage, cut_tree
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
import seaborn as sns
import os

pdbFileKeys=()
pdbFileindex=[]
CDRH3_dict={}
CDRL3_dict={}
with open("C:/Users/C00479477/PycharmProjects/pythonProject/sample_details_Corrected.csv") as f:
    reader=csv.reader(f)
    next(reader)
    for row in reader:
        path='C:/Users/C00479477/PycharmProjects/pythonProject/key_{}_CDRH3.keys_theta29_dist17'.format(row[3])
        isExist = os.path.exists(path)
        if isExist == True:
            CDRH3_dict[row[3]]=row[1]
            CDRL3_dict[row[3]]=row[2]
            file = open(path, 'r')
            lines = file.readlines()
            lines.remove(lines[0])
            DictKeys = {}
            for line in lines:
                keys = line.split()

                if len(keys)==2:

                    DictKeys[keys[0]] = int(keys[1])
            pdbFileKeys += (DictKeys,)
            pdbFileindex.append('{}_CDRH3'.format(row[3]))

df=pd.DataFrame(pdbFileKeys,index= pdbFileindex).fillna(0)

#Calculation of similarity Matrix
noOfPdbFiles=len(df)
similarity_Matrix=np.empty(shape=(noOfPdbFiles,noOfPdbFiles))
l=0
for index1,row1 in df.iterrows():
    k=0
    for index2, row2 in df.iterrows():
        #print(index1,row1)
        #print(index2,row2)
        sum1=0
        sum2=0
        for i in range(len(row1)):
            sum1+=min(row1[i],row2[i])
            sum2+=max(row1[i],row2[i])
        similarity=(sum1/sum2)*100
        similarity_Matrix[l][k]=int(similarity)
        k+=1
    l+=1

#Importing the similarity_Matrix to a csv file
df_2=pd.DataFrame(similarity_Matrix,columns=pdbFileindex,index=pdbFileindex)
df_2.to_csv("similarityMatrix.csv",index=True)


similarity_Matrix=abs(similarity_Matrix-100)
df_3=pd.DataFrame(similarity_Matrix,index=pdbFileindex,columns=pdbFileindex)
dists = squareform(df_3)
"""
#Separate dataframe to make specific file 'KeyFreq_theta29_dist17.csv'
first_column = list(df.iloc[:, 0])
pdbFileindex_2=[]
for i in range(len(df)):
    pdbFileindex_2.append(pdbFileindex[i]+";"+str(first_column[i]))
df.iloc[:, 0]=pdbFileindex_2
df.to_csv('KeyFreq_theta29_dist17.csv',header=False,index=False)
"""

# Plotting dendrogam (uses average distance)
linkage_matrix = linkage(dists, "average")
dendrogram(linkage_matrix, labels=pdbFileindex)
# plt.xlabel("PDB File(Protein+Drug )", fontsize=15)

plt.ylabel("Distance", fontsize=15)
plt.title("Dendogram", fontsize=15)
plt.savefig('Dendrogram.png')
plt.show()

# Plotting clustermap
sns.clustermap(df_3,row_linkage=linkage_matrix, col_linkage=linkage_matrix, annot=False)
plt.savefig('Heatmap.png')
plt.show()

n_clusters = range(1, len(pdbFileindex))
clusters=cut_tree(linkage_matrix,n_clusters )
clusters = np.insert(clusters,clusters.shape[1],range(clusters.shape[0]), axis=1)
clusters=clusters.T

Clustered_Data_File=open("Clustered_data_every_step.CDRH3",'w')
Clustered_Data_File.writelines('Number_of_Cluster----Clustered_PDB_ID\n')

amionAcidEachGroup=open('amionAcidEachGroup.CDRH3','w')
amionAcidEachGroup.writelines('Groups\tNo_aa\n')
for row in clusters[::-1]:
    groups={}
    for i,g in enumerate(row):
        if g not in groups:
            groups[g]=set([])
        groups[g].add(pdbFileindex[i][0:4])
    clustered_data=list(groups.values())
    Clustered_Data_File.writelines('{}----{}\n'.format(len(clustered_data),clustered_data))
    if len(clustered_data)==4:
        groupIndex = 1
        for i in clustered_data:
            amionAcidEachGroup.writelines('Group_{}\n'.format(groupIndex))
            for j in i:
                if j in CDRH3_dict:
                    amionAcidEachGroup.writelines('{}\t{}\n'.format(j, len(CDRH3_dict[j])))
                else:
                    print('{} not found'.format(j))


            groupIndex += 1

Clustered_Data_File.close()
amionAcidEachGroup.close()
