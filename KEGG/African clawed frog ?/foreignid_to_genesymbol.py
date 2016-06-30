with open("xenopus_kegg_genesymbol.csv", "r") as translate_file:
    foreign_new_id = {}
    translate_file.readline()
    for line in translate_file:
        line = line.strip().split("\t")
        if len(line) > 1:
            line[1] = line[1].split(";")
            for id in line[1]:
                if line[0] not in foreign_new_id:
                    foreign_new_id[line[0]] = []
                foreign_new_id[line[0]].append(id)

with open("KEGG_Paths2geneIDs", "r") as input_file:
    input_file.readline()
    pathway_gene = {}
    for line in input_file:
        line = line.strip().split("\t")
        kegg_ids = line[1].split(" ")
        if line[0] not in pathway_gene:
            pathway_gene[line[0]] = []
        pathway_gene[line[0]] = kegg_ids

pathway_genesymbol_list = {}
for pathway, gene_list in pathway_gene.items():
    for gene in gene_list:
        if gene in foreign_new_id:
            if pathway not in pathway_genesymbol_list:
                pathway_genesymbol_list[pathway] = []
            pathway_genesymbol_list[pathway].append(foreign_new_id[gene])
        else:
            if pathway not in pathway_genesymbol_list:
                pathway_genesymbol_list[pathway] = []
            pathway_genesymbol_list[pathway].extend(foreign_new_id[gene])
print (pathway_genesymbol_list)
#with open("xenopus_pathway_genesymbol.csv", "w") as output:
 #   for pathway, genesymbols in pathway_genesymbol_list.items():
  #      for gs in genesymbols:
   #         output.write(pathway + ";" + ";".join(gs) + "\n")

