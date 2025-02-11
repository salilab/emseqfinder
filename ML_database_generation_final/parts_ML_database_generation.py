import sys
import os

mapfiletxt = sys.argv[1]

with open(mapfiletxt, 'r') as txtfile:
	expmappdb_list = txtfile.readlines()[1:]

for expmappdb in expmappdb_list:
	# EMDB, PDB, Chain Ids, resolution and thresholding values 
	expmappdb = expmappdb.rstrip()
	emdb = expmappdb.split('_')[0]
#	resolution = expmappdb.split('_')[3] # original format see file pdblist_part.txt
	resolution = expmappdb.split('_')[1]
#	threshold = expmappdb.split('_')[4]
	threshold = expmappdb.split('_')[2]
	# first normalize the map
	os.system('python3 ' + './normalize_map_for_parts_fiiting.py' + ' ' +  emdb + ' --thresh ' + threshold)  
	# now extract side chain densities
	datbase_path = emdb + '_ML_side.dat'
	os.system('python3 ' + './get_database_for_one_emdb_using_parts.py' + ' ' + emdb + ' ' + datbase_path + ' ' + resolution)  
