import csv
genenamechoice =[]
geneidchoice = []
with open( 'genelist.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    for r in reader:
        geneidchoice.append([r[1],r[0]])
        genenamechoice.append([r[1],r[2]])

geneidchoice =geneidchoice [1:]
genenamechoice = genenamechoice [1:]


with open( 'ccstrainmapid.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    strain_list = [[e for e in r] for r in reader]

strain_list = strain_list[1:]


with open( 'consequence.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    consequence_list = [[e for e in r] for r in reader]

consequence_list = consequence_list[1:]

print geneidchoice
print genenamechoice
print strain_list
print consequence_list
