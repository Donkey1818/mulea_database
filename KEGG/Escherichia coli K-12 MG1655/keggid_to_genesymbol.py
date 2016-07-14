with open("KEGG_Paths2geneIDs", "r") as input_file:
    input_file.readline()
    id = []
    for line in input_file:
        line = line.strip().split("\t")
        kegg_ids = line[1].split(" ")
        for kegg_id in kegg_ids:
            id.append(kegg_id)
    id = set(id)

with open("Escherichia_coli_keggid_list.csv", "w") as output:
    for keggid in id:
        output.write(keggid + "\n")

