
import csv
geneidchrchoice = []
with open( 'all_genes.txt', 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter='\t')
    for r in reader:
        geneidchrchoice.append([r[1],r[2],r[0]])

geneidchrchoice =geneidchrchoice [1:]



print geneidchrchoice
