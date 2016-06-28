#notes for creating up-to-date KEGG database

#########
# R
library(KEGGREST)
# all genes (ncbi-geneid) of all pathways
AllGenes2Path=keggLink("pathway", "ptr")

# geneIDs of a single pathway:
names(AllGenes2Path[AllGenes2Path %in% 'path:ptr04726'])


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
# delete ptr: and path:
perl -pe 's/ptr\:|path\:|Paths2geneIDs \n//g' KEGG_Paths2geneIDs | perl -pe 's/ \t /\t/g' | perl -pe 's/ \n/\n/g' > temp
mv temp KEGG_Paths2geneIDs

# names of pathways
wget http://rest.kegg.jp/list/pathway/ptr
mv ptr ptr.pathways.list

perl -pe 's/path\:| - Pan troglodytes \(chimpanzee\)//g' ptr.pathways.list | perl -pe 's/\//-/g' | perl -pe 's/ /_/g' | awk '{print $1"_"$2}' > temp
mv temp ptr.pathways.list



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
Try=queryMany(TryGenes, scopes="entrezgene", fields=c("symbol", "uniprot"), species="chimpanzee")
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
MyGeneList=queryMany(EntrezAll, scopes="entrezgene", fields=c("symbol", "uniprot"), species="chimpanzee")
ConvTab=as.data.frame(MyGeneList)

## change KEGG_Paths2geneIDs 2 KEGG_Paths2geneSymbols
library(plyr)

# mapvalues(x, from, to)
mapvalues(KEGG_Paths2geneIDs[[1]], ConvTab[, "query"], ConvTab[, "symbol"], warn_missing=F)


ConvFuncE2G=function(Vect) {
	mapvalues(Vect, ConvTab[, "query"], ConvTab[, "symbol"], warn_missing=F)
}

KEGG_Paths2geneSymbols=lapply(KEGG_Paths2geneIDs, ConvFuncE2G) #hiba!!!!



## give names of pathways
PWnames=readLines("ptr.pathways.list")
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

# delete ptr: and path:
perl -pe 's/KEGG_Paths2geneIDs2 \n//g' KEGG_Paths2geneIDs | perl -pe 's/ \t /\t/g' | perl -pe 's/ \n/\n/g' > temp
mv temp KEGG_Paths2geneIDs
