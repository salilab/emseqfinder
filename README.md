\brief Assignment of sequence to backbone fragments traced in a cryo-EM map 

# Info

_Author(s)_: Dibyendu Mondal, Vipul Kumar

_Maintainer_: `benmwebb`

_License_: [LGPL](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Publications_:
 - D. Mondal, V. Kumar, T. Satler, R. Ramachandran, D. Saltzberg, I. Chemmama, K.B. Pilla, I. Echeverria, B.M. Webb, M. Gupta, K.A. Verba, A. Sali. Recognizing amino acid sidechains in a medium resolution cryo-electron density map. Prot Sci 34, e70217, 2025.
 - See [main IMP papers list](@ref publications).

# emseqfinder: command line tool to run emseqfinder protocol {#emseqfinder}

The protocol is typically run using the `emseqfinder` command line tool:

## How to run the program to get sequence assignments
## 1. We need a database folder. example "test"
## 2. create and folder with anything like "emdbid"
## 3. Make two folders '0system' and '1structure_elements' inside "emdbid" folder
## 4. Put sequence file (.fasta), EM map (.map) and PDB file (.pdb) associated to that EM map inside emdbid/0system
## 5. Get the PDB file and run stride to get the secodnary structure annotations.
## 6. Then run ./fragdb_generation/get_fraglib_from_native.py script to generate the fragements. Copy all the fragments to emdbid/1structure_elements

## 7. Run ./ML_database_generation_final/parts_ML_database_generation.py script on testrun.txt
## 8. Please check configuration of a testrun.txt file in the ./ML_database_generation_final folder.
## 9. once the ML database is generated from the EMDB voxel and some other informations, go head to run the ML prediction script in this folder.
## 10. CNN prediction
	Filenames: 1. CNN run script: final_ML_predict.py
		2. pickling the input before running CNN: convert_MLDB_topkl.py
##11. Now you can use evalute_output_database.py and calculate_seq_match.py for analysis.


###Note. Database genration for CNN prediction:
###Foldername: /wynton/home/sali/dibyendu/EMthreading_ML/gpu_sandbox/version2/database_generation_evaluation_scripts/EM_threading_Pipeline/ML_database_generation_final

	- We need a cryoEM map
	- We need backbone traces (please only consider helices and strands)	
	- We use the backbone traces to extract side chain densities from the cryoEM map. Before extracting the side chain denisties, we normalize the map.
	The map normalization script is normalize_map_for_parts_fiiting.py
