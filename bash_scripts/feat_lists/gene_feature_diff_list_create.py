#Script that outputs a list of the features unique to new_ncbi_feats and feat_list_from_GEO

# function for turning csv or tsv into list of lists
def sv_to_lines(file, sep):
 fileread = file.readlines()
 table = []
 for line in fileread:
         x = line.split(sep)
         table.append(x)
 return table

inpo = 'feat_list_from_GEO.txt'
inpn = 'feat_list_ncbi.txt'

inpold = open(inpo, 'r')
tableo = inpold.readlines()
listo = []
for i in tableo:
	listo.append(i.rstrip(' \n'))

inpnew = open(inpn, 'r')
tablen = inpnew.readlines()
listn = []
for i in tablen:
	listn.append(i.rstrip('\t\n'))

#Unique to Old ref
list_diff = []
#Listn is unique to new ref
for i in listo:
	if i in listn:
	    listn.remove(i)
	else:
		list_diff.append(i)

outn= open('new_ncbi_feats', 'a')
for i in listn:
	outn.write(i + '\n')

outu=open('outdated_feats', 'a')
for i in list_diff:
	outu.write(i + '\n')



