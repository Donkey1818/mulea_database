http://flybase.org/cgi-bin/serveHTdata.cgi?dataset=modENCODE_mRNA-Seq_tissues&FBgn=FBgn0002521

wget "http://flybase.org/cgi-bin/serveHTdata.cgi?dataset=modENCODE_mRNA-Seq_tissues&FBgn=FBgn0265426" -O - >> out

No/Extremely low expression (0 - 0)
Very low expression (1 - 3)
Low expression (4 - 10)
Moderate expression (11 - 25)
Moderately high expression (26 - 50)
High expression (51 - 100)
Very high expression (101 - 1000)
Extremely high expression (>1000)

https://www.hgsc.bcm.edu/drosophila-modencode-project
http://www.modencode.org/

Enrichment: the expression relative to the average in whole flies (average of all tissues)
=G2/AVERAGE(G$2:G$30)


#Get all D.mel genes FBgn-s
cd /Volumes/Temp/Ari/Dmel_genome
wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/dmel-all-r5.55.gff.gz
gunzip dmel-all-r5.55.gff.gz
awk -F"\t" '$2=="FlyBase" && $3=="mRNA"' dmel-all-r5.55.gff | grep -o -e "FBgn......." | sort -u > dmel_FBgn_sorted
cp dmel_FBgn_sorted /Users/vetgrid09/Dropbox/Dsim_thermal/Results_hypergeom/09_Hypergeom/Analysis_of_organ_spec_modENCODE/.


#All FBgn-s in D.sim annotation:
/Volumes/Temp/Ari/Dsim_thermal/READS+GENOMES/gene_summary/dsim_FBgn_sorted
cp dsim_FBgn_sorted /Users/vetgrid09/Dropbox/Dsim_thermal/Results_hypergeom/09_Hypergeom/Analysis_of_organ_spec_modENCODE/.


nohup ./download_modENCODE_from_FB.sh dmel_FBgn_sorted &
mv modENCODE_raw.tsv dmel_modENCODE_raw.tsv
nohup ./download_modENCODE_from_FB.sh dsim_FBgn_sorted &
mv modENCODE_raw.tsv dsim_modENCODE_raw.tsv

#Did I get all the genes?
wc -l dmel_FBgn_sorted
#13939
grep -c "0K \.\." nohup.out
#13931
awk -F"\t" '{print $4}' dmel_modENCODE.tsv | sort -u | wc -l
#13932 (with the header)
awk -F"\t" '{print $4}' dmel_modENCODE.tsv | sort -u > dmel_FBgn_sorted_from_modENCODE
diff dmel_FBgn_sorted dmel_FBgn_sorted_from_modENCODE | grep "<"
#< FBgn0013673
#< FBgn0013680
#< FBgn0013681
#< FBgn0013683
#< FBgn0013684
#< FBgn0013685
#< FBgn0030765
#< FBgn0083077
#no modENCODE data for these genes

#####################
#calculate enrichment with R
tab=read.table("dmel_modENCODE.tsv", header=T)
#colnames(tab)=c("Library_Collection_ID", "Parent_Library_Name", "Gene_Symbol", "FBgn_ID", "Category_ID", "Category_Descriptive_Name", "RPKM_Expression_Value", "Expression_Value_BinID")

#check if all genes has all categories
categ=as.character(unique(tab$Category_Descriptive_Name))
categ_length=c()
for (i in 1:length(categ)) {
categ_length[i]=length(tab[tab$Category_Descriptive_Name==categ[i],1])
}
categ_length
# [1] 13931 13931 13931 13931 13931 13931 13931 13931 13931 13931 13931 13931
#[13] 13931 13931 13931 13931 13931 13931 13931 13931 13931 13931 13931 13931
#[25] 13931 13931 13931 13931 13931

#create a new col. for the average of RPKMs per gene
genes=as.character(unique(tab$FBgn_ID))
tab$Average=NA
for (i in 1:length(genes)) {
  tab[tab$FBgn_ID==genes[i],"Average"]=c(mean(tab[tab$FBgn_ID==genes[i],"RPKM_Expression_Value"]))
}

tab$Enrichment=tab$RPKM_Expression_Value/tab$Average

write.table(tab[,c(4,6,7,10)], file="dmel_modENCODE_with_enrichment.tsv", sep="\t", quote=F, col.names=T, row.names=F)

###########x
#shell

awk '{print $2}' dmel_modENCODE_with_enrichment.tsv | sort -u | grep -v "Category_Descriptive_Name" > categories

#awk separator
awk -F"\t" '{print $4 > $2}' dmel_modENCODE_with_enrichment2.tsv

perl -pe 's/\n/ /g' categories > All_modE_enrich_table

paste dmel_FBgn_sorted_from_modENCODE mE_mRNA_A_1d_carcass mE_mRNA_A_1d_dig_sys mE_mRNA_A_20d_carcass mE_mRNA_A_20d_dig_sys mE_mRNA_A_4d_carcass mE_mRNA_A_4d_dig_sys mE_mRNA_A_MateF_1d_head mE_mRNA_A_MateF_20d_head mE_mRNA_A_MateF_4d_head mE_mRNA_A_MateF_4d_ovary mE_mRNA_A_MateM_1d_head mE_mRNA_A_MateM_20d_head mE_mRNA_A_MateM_4d_acc_gland mE_mRNA_A_MateM_4d_head mE_mRNA_A_MateM_4d_testis mE_mRNA_A_VirF_1d_head mE_mRNA_A_VirF_20d_head mE_mRNA_A_VirF_4d_head mE_mRNA_A_VirF_4d_ovary mE_mRNA_L3_CNS mE_mRNA_L3_Wand_carcass mE_mRNA_L3_Wand_dig_sys mE_mRNA_L3_Wand_fat mE_mRNA_L3_Wand_imag_disc mE_mRNA_L3_Wand_saliv mE_mRNA_P8_CNS mE_mRNA_P8_fat mE_mRNA_WPP_fat mE_mRNA_WPP_saliv >> All_modE_enrich_table

perl -pe 's/\t/ /g' All_modE_enrich_table > All_modE_enrich_table2
mv All_modE_enrich_table2 All_modE_enrich_table

######
#R
AllEnrich_table=read.table(file="All_modE_enrich_table", header=T, row.names=1)

#delete genes with all NA-s
dim(AllEnrich_table)
#[1] 13931    29
dim(AllEnrich_table[complete.cases(AllEnrich_table),])
#[1] 13241    29
AllEnrich_table=AllEnrich_table[complete.cases(AllEnrich_table),]

min(AllEnrich_table)
#[1] 0
max(AllEnrich_table)
#[1] 29

summary(AllEnrich_table)


summary(AllEnrich_table<=1)
c(summary(AllEnrich_table<1))
# the result had been saved to number_of_genes_in_each_category.xls

#categories, colnames
categ=sort(categ)

#create lists for every category of genes that enrichment value is higher than 1
mE_list=list()
for (i in 1:length(categ)) {
  mE_list[[i]]=row.names(AllEnrich_table[AllEnrich_table[,i]>1,])
}
names(mE_list)=categ

#function for write a list
fnlist=function(x, fil){
    z=deparse(substitute(x))
    cat(z, "\n", file=fil)
    nams=names(x)
    for (i in seq_along(x)){
        cat(nams[i], "\t", x[[i]], "\n", file=fil, append=T)
        }
    }
fnlist(mE_list, "modENCODE_anatomy_database_dmel.txt")

######################
#shell
#delete the first line
grep -v "mE_list" modENCODE_anatomy_database_dmel.txt > modENCODE_anatomy_database_dmel_2.txt 

#delete all new line characters
sed ':a;N;$!ba;s/\n//g' modENCODE_anatomy_database_dmel_2.txt > modENCODE_anatomy_database_dmel_3.txt
#change mE_mRNA_ to new line
sed 's/mE_mRNA_/\n/g' modENCODE_anatomy_database_dmel_3.txt > modENCODE_anatomy_database_dmel_4.txt
#delete empty lines (contains only space)
sed -e '/^ *$/d' modENCODE_anatomy_database_dmel_4.txt > modENCODE_anatomy_database_dmel_5.txt

mv modENCODE_anatomy_database_dmel_5.txt modENCODE_anatomy_database_dmel.txt
rm modENCODE_anatomy_database_dmel_2.txt
rm modENCODE_anatomy_database_dmel_3.txt
rm modENCODE_anatomy_database_dmel_4.txt

#rename categories
cp modENCODE_anatomy_database_dmel.txt x.txt
sed -i 's/A_1d_carcass/1_day_adult_carcass/g' x.txt
sed -i 's/A_1d_dig_sys/1_day_adult_digestive_system/g' x.txt
sed -i 's/A_20d_dig_sys/20_day_adult_digestive_system/g' x.txt
sed -i 's/A_4d_carcass/4_day_adult_carcass/g' x.txt
sed -i 's/A_4d_dig_sys/4_day_adult_digestive_system/g' x.txt
sed -i 's/A_MateF_1d_head/1_day_adult_mated_female_head/g' x.txt
sed -i 's/A_MateF_20d_head/20_day_adult_mated_female_head/g' x.txt
sed -i 's/A_MateF_4d_head/4_day_adult_mated_female_head/g' x.txt
sed -i 's/A_MateF_4d_ovary/4_day_adult_mated_female_ovary/g' x.txt
sed -i 's/A_MateM_1d_head/1_day_adult_mated_male_head/g' x.txt
sed -i 's/A_MateM_20d_head/20_day_adult_mated_male_head/g' x.txt
sed -i 's/A_MateM_4d_acc_gland/4_day_adult_mated_male_accessory_gland/g' x.txt
sed -i 's/A_MateM_4d_head/4_day_adult_mated_male_head/g' x.txt
sed -i 's/A_MateM_4d_testis/4_day_adult_mated_male_testis/g' x.txt
sed -i 's/A_VirF_1d_head/1_day_adult_virgin_female_head/g' x.txt
sed -i 's/A_VirF_20d_head/20_day_adult_virgin_female_head/g' x.txt
sed -i 's/A_VirF_4d_head/4_day_adult_virgin_female_head/g' x.txt
sed -i 's/A_VirF_4d_ovary/4_day_adult_virgin_female_ovary/g' x.txt
sed -i 's/L3_CNS/larvae_L3_central_nervous_system/g' x.txt
sed -i 's/L3_Wand_carcass/wandering_larvae_L3_carcass/g' x.txt
sed -i 's/L3_Wand_dig_sys/wandering_larvae_L3_digestive_system/g' x.txt
sed -i 's/L3_Wand_fat/wandering_larvae_L3_fat_body/g' x.txt
sed -i 's/L3_Wand_imag_disc/wandering_larvae_L3_imaginal_disc/g' x.txt
sed -i 's/L3_Wand_saliv/wandering_larvae_L3_salivary_gland/g' x.txt
sed -i 's/P8_CNS/pupae_P8_central_nervous_system/g' x.txt
sed -i 's/P8_fat/pupae_P8_fat_body/g' x.txt
sed -i 's/WPP_fat/white_prepupae_fat_body/g' x.txt
sed -i 's/WPP_saliv/white_prepupae_salivary_gland/g' x.txt
mv x.txt modENCODE_anatomy_database_dmel.txt

#the database is ready to use!

#######################################x
#create lists for every category of genes that enrichment value is higher than 2
mE_list_2=list()
for (i in 1:length(categ)) {
  mE_list_2[[i]]=row.names(AllEnrich_table[AllEnrich_table[,i]>2,])
}
names(mE_list_2)=categ

#function for write a list
fnlist=function(x, fil){
    z=deparse(substitute(x))
    cat(z, "\n", file=fil)
    nams=names(x)
    for (i in seq_along(x)){
        cat(nams[i], "\t", x[[i]], "\n", file=fil, append=T)
        }
    }
fnlist(mE_list_2, "modENCODE_anatomy_database_enr2_dmel.txt")

######################
#shell
#delete the first line
grep -v "mE_list" modENCODE_anatomy_database_enr2_dmel.txt > modENCODE_anatomy_database_enr2_dmel_2.txt 

#delete all new line characters
sed ':a;N;$!ba;s/\n//g' modENCODE_anatomy_database_enr2_dmel_2.txt > modENCODE_anatomy_database_enr2_dmel_3.txt
#change mE_mRNA_ to new line
sed 's/mE_mRNA_/\n/g' modENCODE_anatomy_database_enr2_dmel_3.txt > modENCODE_anatomy_database_enr2_dmel_4.txt
#delete empty lines (contains only space)
sed -e '/^ *$/d' modENCODE_anatomy_database_enr2_dmel_4.txt > modENCODE_anatomy_database_enr2_dmel_5.txt

mv modENCODE_anatomy_database_enr2_dmel_5.txt modENCODE_anatomy_database_enr2_dmel.txt
rm modENCODE_anatomy_database_enr2_dmel_2.txt
rm modENCODE_anatomy_database_enr2_dmel_3.txt
rm modENCODE_anatomy_database_enr2_dmel_4.txt

#rename categories
cp modENCODE_anatomy_database_enr2_dmel.txt x.txt
sed -i 's/A_1d_carcass/1_day_adult_carcass/g' x.txt
sed -i 's/A_1d_dig_sys/1_day_adult_digestive_system/g' x.txt
sed -i 's/A_20d_dig_sys/20_day_adult_digestive_system/g' x.txt
sed -i 's/A_4d_carcass/4_day_adult_carcass/g' x.txt
sed -i 's/A_4d_dig_sys/4_day_adult_digestive_system/g' x.txt
sed -i 's/A_MateF_1d_head/1_day_adult_mated_female_head/g' x.txt
sed -i 's/A_MateF_20d_head/20_day_adult_mated_female_head/g' x.txt
sed -i 's/A_MateF_4d_head/4_day_adult_mated_female_head/g' x.txt
sed -i 's/A_MateF_4d_ovary/4_day_adult_mated_female_ovary/g' x.txt
sed -i 's/A_MateM_1d_head/1_day_adult_mated_male_head/g' x.txt
sed -i 's/A_MateM_20d_head/20_day_adult_mated_male_head/g' x.txt
sed -i 's/A_MateM_4d_acc_gland/4_day_adult_mated_male_accessory_gland/g' x.txt
sed -i 's/A_MateM_4d_head/4_day_adult_mated_male_head/g' x.txt
sed -i 's/A_MateM_4d_testis/4_day_adult_mated_male_testis/g' x.txt
sed -i 's/A_VirF_1d_head/1_day_adult_virgin_female_head/g' x.txt
sed -i 's/A_VirF_20d_head/20_day_adult_virgin_female_head/g' x.txt
sed -i 's/A_VirF_4d_head/4_day_adult_virgin_female_head/g' x.txt
sed -i 's/A_VirF_4d_ovary/4_day_adult_virgin_female_ovary/g' x.txt
sed -i 's/L3_CNS/larvae_L3_central_nervous_system/g' x.txt
sed -i 's/L3_Wand_carcass/wandering_larvae_L3_carcass/g' x.txt
sed -i 's/L3_Wand_dig_sys/wandering_larvae_L3_digestive_system/g' x.txt
sed -i 's/L3_Wand_fat/wandering_larvae_L3_fat_body/g' x.txt
sed -i 's/L3_Wand_imag_disc/wandering_larvae_L3_imaginal_disc/g' x.txt
sed -i 's/L3_Wand_saliv/wandering_larvae_L3_salivary_gland/g' x.txt
sed -i 's/P8_CNS/pupae_P8_central_nervous_system/g' x.txt
sed -i 's/P8_fat/pupae_P8_fat_body/g' x.txt
sed -i 's/WPP_fat/white_prepupae_fat_body/g' x.txt
sed -i 's/WPP_saliv/white_prepupae_salivary_gland/g' x.txt
mv x.txt modENCODE_anatomy_database_enr2_dmel.txt

#the database is ready to use!










