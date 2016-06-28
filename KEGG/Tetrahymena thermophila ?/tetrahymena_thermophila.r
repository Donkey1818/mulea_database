#notes for creating up-to-date KEGG database

#########
# R
library(KEGGREST)
# all genes (ncbi-geneid) of all pathways
AllGenes2Path=keggLink("pathway", "tet")

# geneIDs of a single pathway:
names(AllGenes2Path[AllGenes2Path %in% 'path:tet00030'])

# create a list that contains the genes of each pathways
Paths=unique(AllGenes2Path)
Paths2geneIDs=list()
for(i in Paths) {
	Paths2geneIDs[[i]]=names(AllGenes2Path[AllGenes2Path %in% i])
}

#function to write a list to a file
fnlist=function(x, fil){
    z=deparse(substitute(x))
    cat(z, "\n", file=fil)
    nams=names(x)
    for (i in seq_along(x)){
        cat(nams[i], "\t", x[[i]], "\n", file=fil, append=T)
        }
    }
fnlist(Paths2geneIDs, "KEGG_Paths2geneIDs")

# shell
# delete tet: and path:
perl -pe 's/tet\:|path\:|Paths2geneIDs \n//g' KEGG_Paths2geneIDs | perl -pe 's/ \t /\t/g' | perl -pe 's/ \n/\n/g' > temp
mv temp KEGG_Paths2geneIDs

# names of pathways
wget http://rest.kegg.jp/list/pathway/tet
mv tet tet.pathways.list

perl -pe 's/path\:| - Tetrahymena thermophila \(\)//g' tet.pathways.list | perl -pe 's/\//-/g' | perl -pe 's/ /_/g' | awk '{print $1"_"$2}' > temp
mv temp tet.pathways.list



###############
# R

# create a list of DB IDs and associated genes from text file
x=strsplit(readLines("KEGG_Paths2geneIDs"), "[[:space:]]+")
KEGG_Paths2geneIDs=lapply(x, tail, n=-1)
names(KEGG_Paths2geneIDs)=lapply(x, head, n=1)
rm(x)

## conversion table from entrez to geneSymbol
library(mygene)

TryGenes=sort(KEGG_Paths2geneIDs[[1]])
Try=queryMany(TryGenes, scopes="entrezgene", fields=c("symbol", "uniprot"), species="")
TryD=as.data.frame(Try)
TryD
#    uniprot.TrEMBL uniprot.Swiss.Prot   X_id  symbol  query
# 1    Q5T621, ....             P14550  10327  AKR1A1  10327
# 2                             P07327    124   ADH1A    124
# 3    D6RHZ6, ....             P00325    125   ADH1B    125
# 4    A0A087WU....             P00326    126   ADH1C    126
# 5    A0A0D9SF....             P08319    127    ADH4    127
# 6    D6R9G2, ....             P11766    128    ADH5    128
# 7    A0A0A0MS....             P28332    130    ADH6    130
# 8    B8ZZ75, ....             Q96C23 130589    GALM 130589
# 9    A0A0C4DG....             P40394    131    ADH7    131
# 10     A0A087WUM2             Q6ZMR3 160287 LDHAL6A 160287

EntrezAll=sort(unique(unlist(KEGG_Paths2geneIDs)))
MyGeneList=queryMany(EntrezAll, scopes="entrezgene", fields=c("symbol", "uniprot"), species="")
ConvTab=as.data.frame(MyGeneList)

## change KEGG_Paths2geneIDs 2 KEGG_Paths2geneSymbols
library(plyr)

# mapvalues(x, from, to)
mapvalues(KEGG_Paths2geneIDs[[1]], ConvTab[, "query"], ConvTab[, "symbol"], warn_missing=F)


ConvFuncE2G=function(Vect) {
	mapvalues(Vect, ConvTab[, "query"], ConvTab[, "symbol"], warn_missing=F)
}

KEGG_Paths2geneSymbols=lapply(KEGG_Paths2geneIDs, ConvFuncE2G)



## give names of pathways
PWnames=readLines("tet.pathways.list")
# checking the order
names(PWnames)=names(KEGG_Paths2geneSymbols)

KEGG_Paths2geneSymbols2=KEGG_Paths2geneSymbols
names(KEGG_Paths2geneSymbols2)=PWnames

KEGG_Paths2geneIDs2=KEGG_Paths2geneIDs
names(KEGG_Paths2geneIDs2)=PWnames

fnlist(KEGG_Paths2geneSymbols2, "KEGG_Paths2geneSymbols")
fnlist(KEGG_Paths2geneIDs2, "KEGG_Paths2geneIDs")
############
# shell
perl -pe 's/KEGG_Paths2geneSymbols2 \n//g' KEGG_Paths2geneSymbols | perl -pe 's/ \t /\t/g' | perl -pe 's/ \n/\n/g' > temp
mv temp KEGG_Paths2geneSymbols

# delete tet: and path:
perl -pe 's/KEGG_Paths2geneIDs2 \n//g' KEGG_Paths2geneIDs | perl -pe 's/ \t /\t/g' | perl -pe 's/ \n/\n/g' > temp
mv temp KEGG_Paths2geneIDs























#######################
# creating conversion table

# Uniprot: http://www.uniprot.org/
# Advance search: Organism: Tetrahymena thermophila (it extends the search word automaticly)
# Columns: Search: geneID, delete everything else except Entry Gene names (primary ), save
# Download: all, Format: Tab separated, uncompressed -> All_uniprotID_2_genesymbol_entrezGeneID_

# make uni 2 other -> geneID 2 genesymbol and geneID 2 uni
perl -pe 's/; /;/g' All_uniprotID_2_genesymbol_entrezGeneID_ | perl -pe 's/;\n/\n/g' > temp
mv temp All_uniprotID_2_genesymbol_entrezGeneID_

awk -F"\t" '{print $1 "\t" $2}' All_uniprotID_2_genesymbol_entrezGeneID_ > Uni2geneName

#get list of kegg to uniprot conversion
wget http://rest.kegg.jp/conv/uniprot/tet
perl -pe 's/tet:|up://g' tet > Entrez2Uni
rm tet
## crate the same format as it is in Uni2 files, with ;

# multiple or NO GeneIDs can be connected to a UniprotID; run this on atlasz
# run IC_conv.sh on atlasz with Entrez2Uni KEGG_Paths2geneIDs KEGG_Paths2Uni (as outfile)
# run IC_conv.sh on atlasz with Uni2geneName KEGG_Paths2Uni KEGG_Paths2GeneName (as outfile)







#get list of pathways for Tetrahymena thermophila:
wget http://rest.kegg.jp/list/pathway/tet  # -> tet
mv tet tet.pathways.list
awk -F"\t" '{print $1}' tet.pathways.list > pathwayIDs

#get list of kegg to uniprot conversion
wget http://rest.kegg.jp/conv/uniprot/tet
mv tet tet.kegg_uniprot.list
#get list of kegg to geneID conversion
wget http://rest.kegg.jp/conv/ncbi-geneid/tet
mv tet tet.kegg_geneID.list

#get conversion table for geneIDs (CG...) to flybase IDs from http://flybase.org/static_pages/downloads/IDConv.html
#get list of genes in pathway:
#http://rest.kegg.jp//link/genes/path:dme00010

#download pathway list, run:
cat tet.pathways.list | python get_pathways.py > pathways_gene
#this gives a list of all gene IDs in these pathways and a gowinda file
#gene ids: 2535 
#create conversion table using  http://flybase.org/static_pages/downloads/IDConv.html

#run:
python change_ids_kegg_fb.py --in pathways_gene --ids FlyBase_IDs.txt --out kegg_pathways_flybase.gowinda
#to create a gowinda input file with flybase IDs

perl -pe 's/:/_/g' kegg_pathways_flybase.gowinda | perl -pe 's/\-Drosophilamelanogaster\(fruitfly\)\(|\)|:|\/|\-fly|\(fruitfly|\-Drosophilamelanogaster//g' | awk -F"\t" '{print $1 "_" $2 "\t" $3}' > kegg_dmel_database

#############
#Reactome:
#get reactome gene file from biomart:
http://central.biomart.org/
#get list of unique genes:
awk  -F"\t" ' /^REAC/ {print $NF}' reactome_pathways_ensemble_genes.txt | sort | uniq > reactome_genes_uniq
#create conversion table using  http://flybase.org/static_pages/downloads/IDConv.html
#run:
python create_gowinda_from_reactome.py --in reactome_pathways_ensemble_genes.txt --ids FlyBase_IDs_reactome.txt --out reactome_pathways_flybase.gowinda






