with open("Drosophila_melanogaster_modEncode.txt", "r") as raw_input_file:
    categori_expression = {}
    tissue_expression = {}
    raw_input_file.readline()
    for line in raw_input_file:
        line = line.strip().split("\t")
        if len(line) == 8 and line[0][0] != "#":
            expression_number = int(line[6])
            if expression_number == 0:
                expression_categori = "No/Extremely low expression"
            elif expression_number >= 1 and expression_number <= 3:
                expression_categori = "Very low expression"
            elif expression_number >= 4 and expression_number <= 10:
                expression_categori = "Low expression"
            elif expression_number >= 11 and expression_number <= 25:
                expression_categori = "Moderate expression"
            elif expression_number >= 26 and expression_number <= 50:
                expression_categori = "Moderately high expression"
            elif expression_number >= 51 and expression_number <= 100:
                line[6] = "High expression"
            elif expression_number >= 101 and expression_number <= 1000:
                expression_categori = "Very high expression"
            elif expression_number > 100:
                expression_categori = "Extremely high expression"
            if line[4] not in categori_expression:
                categori_expression[line[4]] = {}
            if expression_categori not in categori_expression[line[4]]:
                categori_expression[line[4]][expression_categori] = []
            categori_expression[line[4]][expression_categori].append(line[3])
with open ("Drosophila_melanogaster_modEncode.gmt", "w") as output:
    for tissue in categori_expression.keys():
        for expression_categori in categori_expression[tissue].keys():
            output.write(tissue+ "-" + expression_categori + ";"  + tissue+ "-" + expression_categori + ";" + ",".join(categori_expression[tissue][expression_categori]) + "\n")




