

```python
import pandas as pd
import numpy as np
```


```python
cd /your/path/../
```


```python
!ls
```

#### Creating 97% identity abundance tables from ASV tables


```python
#create a fasta file from ASV table
# first convert your ASV_abundance.txt file to an ASV_abundance.fasta with the 'make_fasta_from_ASV_table.sh' script
#see supplementary material for the 'make_fasta_from_ASV_table.sh' script. 
#just change names of the input and output files in the 'make_fasta_from_ASV_table.sh' script
#make sure the 'otu.py' file is in your directory
```


```python
# make your make_fasta_from_ASV_table.sh  and 'oty.py' files executable
chmod a+x ./make_fasta_from_zotu_table.sh 
chmod a+x ./otu.py
```


```python
#run the script
./make_fasta_from_zotu_table.sh 
```


```python
#load usearch
chmod a+x ./usearch11.0.667_i86linux32
```


```python
# use the cluster_fast algorithm using Usearch in cmd line to create clusters of ASVs at 97%. You need the .uc file later
# 
/your/path/to/usearch11.0.667_i86linux32 -cluster_fast ASV_abundance.fasta -id 0.97 -sort size -centroids Projectname_97otu.fasta -uc Projectname_97otuclusters.uc  
# xxxxx uniques (ASVs); yyyyy clusters (97%); max cluster size xxxxx; avg cluster size y.y; xxxxx singletons at 97%.
```


```python
#The uc record will have 'C'lusters, 'H'its to these clusters and 'S'centroids. 
#S and C are repeated info (see below same number of records for each)
!grep '^H' Projectname_97otuclusters.uc | wc -l
!grep '^S' Projectname_97otuclusters.uc | wc -l
!grep '^C' Projectname_97otuclusters.uc | wc -l
```


```python
# create text file of centroids
!grep '^S' Projectname_97otuclusters.uc > Projectname_97otuclusters_S.txt
```


```python
# cut only columns we want 
!cut -f2,9,10 Projectname_97otuclusters_S.txt > Projectname_97otuclusters_S_cut_col.txt
```


```python
#read in table and give headers to columns
B16S_S = pd.read_table('Derwent_99otuclusters_S_cut_col.txt',
                     names=['cluster', '#OTU ID', 'target'])
```


```python
#check
B16S_S.head()
```


```python
#duplicate the OTU seq into the 'target' column
B16S_S['target'] = B16S_S['#OTU ID']
```


```python
#check
B16S_S.head()
```


```python
#repeat above process for Hits
!grep '^H' Projectname_97otuclusters.uc > Projectname_97otuclusters_H.txt #| wc -l
```


```python
!cut -f2,9,10 Projectname_97otuclusters_H.txt > Projectname_97otuclusters_H_cut_col.txt
```


```python
B16S_H = pd.read_table('Projectname_97otuclusters_H_cut_col.txt',
                     names=['cluster', '#OTU ID', 'target'])
                   
```


```python
B16S_H.head()
```


```python
#define the centroid and Hits dataframes
dataframes = [B16S_H, B16S_S]
```


```python
#concatenate the centroids and Hits (this gets rid of replication of cluster data)
df_uc_full = pd.concat(dataframes)
```


```python
#head, tail and shape checks
df_uc_full.head()
```


```python
df_uc_full.tail()
```


```python
df_uc_full.shape #xxxxx number of ASVs, one record for each showing the 97% cluster it belongs to with target being centroid.
```


```python
df_split=df_uc_full.copy()
```


```python
from pandas import Series 
```


```python
#getting rid of the size annotation (after the ';') from #OTU ID column
df_split['#OTU ID'] = df_split['#OTU ID'].str.split(';').apply(Series, 1)#.stack()
```


```python
#set the index to #OTU ID
df_split = df_split.set_index('#OTU ID')
```


```python
df_split.head()
```


```python
#read in the original B16S ASV table
df_B16S_ASV=pd.read_table('ASV_abundance.txt',
                      sep='\t',
                     index_col='ASV',)
```


```python
#check
df_B16S_ASV.head()
```


```python
#merge 97% otu cluster/target table (df_split) with the original ASV table
B16S_merge = pd.merge(df_split, df_B16S_ASV, left_index=True, right_index=True, how='inner')
```


```python
B16S_merge.head()
```


```python
#define unique clusters in table
clusters = B16S_merge.cluster.unique()
```


```python
#define unique targets in table
targets = B16S_merge.target.unique()
```


```python
#check these are the same and agree with usearch cluster_fast output (number of clusters at 97%)
len(clusters), len(targets)
```


```python
#how many ASVs in original table 
!wc -l ASV_abundance.txt
```


```python
#check
B16S_merge.shape
```


```python
#group ASVs by cluster and target and sum reads
otu = B16S_merge.groupby(['cluster', 'target']).sum()
```


```python
#check shape, i.e. are there the correct number of columns (samples) and rows(clusters/otus)
otu.shape
```


```python
#check head and tail that all aligns
otu.head()
```


```python
#rename target to #OTU ID as is norm for otu tables ?not working
otu_rename = otu.rename(columns={'target': '#OTU ID'}, inplace=True)
```


```python
#print without renaming to csv file; opened in excel to fix column header
otu.to_csv('Projectname_97otu_table.csv', sep=',')
```


```python
#Sanity check
#open new 'Projectname_97otu_table.csv' and check sequences depth for a sample. 
# The sequence depth should be the same in 'Projectname_97otu_table.csv' and in the original B16S ASV table (ASV_abundance.txt file).
```
