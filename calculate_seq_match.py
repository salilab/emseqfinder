import sys

result_file = sys.argv[1]
total_residue = 0
total_residue_matched = 0
total_residue_matched_abs = 0
with open(result_file, 'r') as infile:
	for lines in infile:
		line = lines.strip().split('|')
		emdb, resolution, seid = line[0].strip().split(' ')
		se_length = int(seid.split('_')[-1])
		if 'nan' in line[1]:
			print(line)
		else:
			total_residue += se_length
			rank = line[1].strip().split(' ')[0]
			if int(rank) <= se_length/3:
				total_residue_matched += se_length

			if int(rank) == 0:
				total_residue_matched_abs += se_length
print(result_file, "total percentage matched " , round(100*(total_residue_matched/total_residue), 3))
print(result_file, "total abs percentage matched " , round(100*(total_residue_matched_abs/total_residue), 3))
