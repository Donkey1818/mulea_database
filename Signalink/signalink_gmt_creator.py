def gmt_creator(input_file, sep, cell_sep, output_file):
    parameter_gene = {}
    with open(input_file, "r") as input_file:
        input_file.readline()
        for line in input_file:
            line = line.strip().split(sep)
            if cell_sep in line[2]:
                line[2] = line[2].split(cell_sep)
                for parameter in line[2]:
                    if parameter not in parameter_gene:
                        parameter_gene[parameter] = []
                    parameter_gene[parameter].append(line[0])
    with open(output_file, "w") as output:
        for parameter in parameter_gene:
            output.write(parameter + ";" + parameter + ";" + ",".join(parameter_gene[parameter]) + "\n")

gmt_creator("human_0_1_signalink.csv","\t", ",", "human_signalink_pathway.gmt")

