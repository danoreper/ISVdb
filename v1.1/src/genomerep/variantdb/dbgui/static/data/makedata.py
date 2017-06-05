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

for i in range(8,81):
    strain_list[i][1]=strain_list[i][1].split('_')[0]


with open( 'consequence.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    consequence_list = [[e for e in r] for r in reader]

consequence_list = consequence_list[1:]

consequence_unique=[]
for x in consequence_list:
    for y in x[1].split('&'):
        consequence_unique.append(y)

consequence_unique = list(set(consequence_unique))

for x in consequence_unique:
    try: # replace the specific one with categorical one
        aa = [s for s in consequence_list if s[1]==x][0]
        x_index=consequence_list.index(aa)
        consequence_list[x_index][0] = ",".join([s[0] for s in consequence_list if x in s[1].split('&')])
    except: # add a new index of consequence
        aa1=",".join([s[0] for s in consequence_list if x in s[1].split('&')])
        aa2=x
        consequence_list.append([aa1,aa2])

consequence_list = sorted(consequence_list, key=lambda cq : len(cq[1]))

print geneidchoice
print genenamechoice
print strain_list
print consequence_list
