with open("rhesus_monkey_keggid_genesmybol.txt", "r") as translate_file:
    foreign_new_id = {}
    translate_file.readline()
    for line in translate_file:
        line = line.strip().split("\t")
        if len(line) > 1:
            if "; " in line[1]:
                line[1] = line[1].split("; ")
                for id in line[1]:
                    if line[0] not in foreign_new_id:
                        foreign_new_id[line[0]] = []
                    foreign_new_id[line[0]].append(id)
            else:
                if line[0] not in foreign_new_id:
                    foreign_new_id[line[0]] = []
                foreign_new_id[line[0]].append(line[1])

with open("KEGG_Paths2geneIDs", "r") as input_file:
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
            try:
                if pathway not in pathway_genesymbol_list:
                    pathway_genesymbol_list[pathway] = []
                pathway_genesymbol_list[pathway].append(foreign_new_id[gene])
            except KeyError:
                print(gene)

with open("rhesus_monkey_pathway_genesymbol.csv", "w") as output:
    for pathway, genesymbols in pathway_genesymbol_list.items():
        for gs in genesymbols:
            output.write(pathway + ";" + ",".join(gs) + "\n")

def gmt_creator(input_file, sep, output_file):
    parameter_gene = {}
    with open(input_file, "r") as input_file:
        for line in input_file:
            line = line.strip().split(sep)
            if line[0] not in parameter_gene:
                parameter_gene[line[0]] = []
            parameter_gene[line[0]].append(line[1])
    with open(output_file, "w") as output:
        for parameter in parameter_gene:
            output.write(parameter + ";" + ",".join(parameter_gene[parameter]) + "\n")

gmt_creator("rhesus_monkey_pathway_genesymbol.csv",";", "KEGG_Macaca_mulatta_GeneSymbol.gmt")
