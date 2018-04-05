in_file = open("all_proteins.fasta", "r")
out_file = open("all_proteins_names.txt", "w")

for line in in_file:
  line = line.rstrip()
  if len(line) > 1 and ">" in line:
    out_file.write(line[1:] + "\n")

in_file.close()
out_file.close()

  
