from subprocess import call
import sys

atoms = []
connects = []

residue_atom_names = {"ARG":['CA', 'C', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
                      "LYS":['CA', 'C', 'CB', 'CG', 'CD', 'CE', 'NZ'],
                      "ASP":['CA', 'C', 'CB', 'CG', 'OD1', 'OD2'],
                      "GLU":['CA', 'C', 'CB', 'CG', 'CD', 'OE1', 'OE2'],
                      "SER":['CA', 'C', 'CB', 'OG'],
                      "CYS":['CA', 'C', 'CB', 'SG'],
                      "GLY":['CA', 'C'],
                      "ALA":['CA', 'C', 'CB'],
                      "VAL":['CA', 'C', 'CB', 'CG1', 'CG2'],
                      "LEU":['CA', 'C', 'CB', 'CG', 'CD1', 'CD2'],
                      "MET":['CA', 'C', 'CB', 'CG', 'SD', 'CE'],
                      }


new_file_name = sys.argv[1] + ".temp" #"receptor_part.pdb"
new_connect_file_name = sys.argv[1] + ".connect" #"temp_conect"
call(["cp " + sys.argv[1] + " " + new_file_name], shell=True)
call(["grep \"CONECT\" " + new_file_name + " > " + new_connect_file_name], shell=True)

conect_lines = []
f = open(new_connect_file_name, 'r')
for line in f: conect_lines.append(line)
f.close()

prev_resi_num = -1
f = open(sys.argv[1], 'r')
residues = []
residues_conect = []
index_start = 0
residue = []
for line in f:

	line_arr = line.split()
	if line_arr[0] != "ATOM": continue
	resi_num = line_arr[5]
	
	if resi_num != prev_resi_num and len(residue) > 0: # new residue
		residues.append(residue)
		residues_conect.append(conect_lines[index_start:(index_start + len(residue))])
		index_start += len(residue)
		residue = []

	residue.append(line)
	prev_resi_num = resi_num


f.close()

for residue_index,residue in enumerate(residues):

	curr_residue_name = residue[0][17:20]
	curr_conect_line = residues_conect[residue_index]

	if curr_residue_name in residue_atom_names.keys():

		name_index = 0
		for line in residue:
			old_line = line[:11] + "  " + line[13:20]

			new_atom_name = residue_atom_names[curr_residue_name][name_index]
			space_to_residue_name = -1
			if len(new_atom_name) == 1:   space_to_residue_name = "   "
			elif len(new_atom_name) == 2: space_to_residue_name = "  "
			elif len(new_atom_name) == 3: space_to_residue_name = " "
			new_line = line[:11] + "  " + new_atom_name + space_to_residue_name + curr_residue_name

			call(["sed -i \"s/" + old_line + "/" + new_line + "/g\" " + new_file_name], shell=True)
			name_index += 1

	else:

		if curr_residue_name == "THR":

			atom_names = ['CA', 'C', 'CB', 'CG2']
			name_index = 0
			for line in residue:
				old_line = line[:11] + "  " + line[13:20]

				line_arr = line.split()
				if line_arr[2] == "O":
					new_atom_name = "OG1"
				else: 
					new_atom_name = atom_names[name_index]
					name_index += 1

				space_to_residue_name = -1
				if len(new_atom_name) == 1:   space_to_residue_name = "   "
				elif len(new_atom_name) == 2: space_to_residue_name = "  "
				elif len(new_atom_name) == 3: space_to_residue_name = " "
				new_line = line[:11] + "  " + new_atom_name + space_to_residue_name + curr_residue_name

				call(["sed -i \"s/" + old_line + "/" + new_line + "/g\" " + new_file_name], shell=True)

		elif curr_residue_name == "ASN":

			atom_names = ['CA', 'C', 'CB', 'CG']
			name_index = 0
			for line in residue:
				old_line = line[:11] + "  " + line[13:20]

				line_arr = line.split()
				if line_arr[2] == "O":
					new_atom_name = "OD1"
				elif line_arr[2] == "N":
					new_atom_name = "ND2"
				else: 
					new_atom_name = atom_names[name_index]
					name_index += 1

				space_to_residue_name = -1
				if len(new_atom_name) == 1:   space_to_residue_name = "   "
				elif len(new_atom_name) == 2: space_to_residue_name = "  "
				elif len(new_atom_name) == 3: space_to_residue_name = " "
				new_line = line[:11] + "  " + new_atom_name + space_to_residue_name + curr_residue_name

				call(["sed -i \"s/" + old_line + "/" + new_line + "/g\" " + new_file_name], shell=True)

		elif curr_residue_name == "GLN":

			atom_names = ['CA', 'C', 'CB', 'CG', 'CD']
			name_index = 0
			for line in residue:
				old_line = line[:11] + "  " + line[13:20]

				line_arr = line.split()
				if line_arr[2] == "O":
					new_atom_name = "OE1"
				elif line_arr[2] == "N":
					new_atom_name = "NE2"
				else: 
					new_atom_name = atom_names[name_index]
					name_index += 1

				space_to_residue_name = -1
				if len(new_atom_name) == 1:   space_to_residue_name = "   "
				elif len(new_atom_name) == 2: space_to_residue_name = "  "
				elif len(new_atom_name) == 3: space_to_residue_name = " "
				new_line = line[:11] + "  " + new_atom_name + space_to_residue_name + curr_residue_name

				call(["sed -i \"s/" + old_line + "/" + new_line + "/g\" " + new_file_name], shell=True)

		elif curr_residue_name == "ILE":

			CB_line = residue[2] # CA, C, then CB
			CB_line_arr = CB_line.split()
			CB_num = CB_line_arr[1] # first word is CONECT

			unknown_lines = curr_conect_line[3:]
			atom_mapping = {}

			# the atom that only has a single connection to CB_num is CG2
			# the atom that has two connections, one of them to CB_num, is CG1
			# the atom that has a single connection, not to CB_num, is CD1
			for line in unknown_lines:
				line_arr = line.split()
				curr_atom = line_arr[1]
				connected_atoms = line_arr[2:]
				if len(connected_atoms) == 1 and CB_num in connected_atoms: atom_mapping[curr_atom] = "CG2"
				elif len(connected_atoms) > 1: atom_mapping[curr_atom] = "CG1"
				elif len(connected_atoms) == 1 and CB_num not in connected_atoms: atom_mapping[curr_atom] = "CD1"


			atom_names = ['CA', 'C', 'CB']
			name_index = 0
			for line in residue:
				old_line = line[:11] + "  " + line[13:20]

				line_arr = line.split()
				curr_atom = line_arr[1]

				if curr_atom in atom_mapping.keys():
					new_atom_name = atom_mapping[curr_atom]
				else: 
					new_atom_name = atom_names[name_index]
					name_index += 1	

				space_to_residue_name = -1
				if len(new_atom_name) == 1:   space_to_residue_name = "   "
				elif len(new_atom_name) == 2: space_to_residue_name = "  "
				elif len(new_atom_name) == 3: space_to_residue_name = " "
				new_line = line[:11] + "  " + new_atom_name + space_to_residue_name + curr_residue_name

				call(["sed -i \"s/" + old_line + "/" + new_line + "/g\" " + new_file_name], shell=True)

		elif curr_residue_name == "HIS":

			CG_line = residue[3] # CA, C, CB, then CG
			CG_line_arr = CG_line.split()
			CG_num = CG_line_arr[1] # first word is CONECT

			unknown_lines_beg = 4
			atom_mapping = {}

			# finding the N atom that connects to CG_num must be ND1
			# finding the C atom that connects to CG_num must be CD2
			# finding the C atom remaining that connects to ND1 must be CE1
			# finding the N atom remaining that connects to CE1 must be NE2
			for mode in range(4):
				for line_index in range(unknown_lines_beg, len(residue)):
					c_line = curr_conect_line[line_index]
					a_line = residue[line_index]
					a_line_arr = a_line.split()
					curr_atom = a_line_arr[1]
					elem = a_line_arr[2]
					c_line_arr = c_line.split()
					connected_atoms = c_line_arr[2:]

					if mode == 0:
						if elem == 'N' and CG_num in connected_atoms: atom_mapping[curr_atom] = "ND1" 
					elif mode == 1:
						if elem == 'C' and CG_num in connected_atoms: atom_mapping[curr_atom] = "CD2"
					elif mode == 2:
						if elem == 'C' and curr_atom not in atom_mapping.keys():
							a0 = False
							if connected_atoms[0] in atom_mapping.keys():
								if atom_mapping[connected_atoms[0]] == 'ND1': a0 = True
							a1 = False
							if connected_atoms[1] in atom_mapping.keys():
								if atom_mapping[connected_atoms[1]] == 'ND1': a1 = True
							if a0 or a1: atom_mapping[curr_atom] = "CE1"
					elif mode == 3:
						if elem == 'N' and curr_atom not in atom_mapping.keys():
							a0 = False
							if connected_atoms[0] in atom_mapping.keys():
								if atom_mapping[connected_atoms[0]] == 'CE1': a0 = True
							a1 = False
							if connected_atoms[1] in atom_mapping.keys():
								if atom_mapping[connected_atoms[1]] == 'CE1': a1 = True
							if a0 or a1: atom_mapping[curr_atom] = "NE2"


			atom_names = ['CA', 'C', 'CB', 'CG']
			name_index = 0
			for line in residue:
				old_line = line[:11] + "  " + line[13:20]

				line_arr = line.split()
				curr_atom = line_arr[1]

				if curr_atom in atom_mapping.keys():
					new_atom_name = atom_mapping[curr_atom]
				else: 
					new_atom_name = atom_names[name_index]
					name_index += 1	

				space_to_residue_name = -1
				if len(new_atom_name) == 1:   space_to_residue_name = "   "
				elif len(new_atom_name) == 2: space_to_residue_name = "  "
				elif len(new_atom_name) == 3: space_to_residue_name = " "
				new_line = line[:11] + "  " + new_atom_name + space_to_residue_name + curr_residue_name

				call(["sed -i \"s/" + old_line + "/" + new_line + "/g\" " + new_file_name], shell=True)

		elif curr_residue_name == "PRO":

			atom_names = ['CA', 'C', 'CB', 'CG', 'CD']
			name_index = 0
			for line in residue:
				old_line = line[:11] + "  " + line[13:20]

				line_arr = line.split()
				if line_arr[2] == "N":
					new_atom_name = "N"
				else: 
					new_atom_name = atom_names[name_index]
					name_index += 1

				space_to_residue_name = -1
				if len(new_atom_name) == 1:   space_to_residue_name = "   "
				elif len(new_atom_name) == 2: space_to_residue_name = "  "
				elif len(new_atom_name) == 3: space_to_residue_name = " "
				new_line = line[:11] + "  " + new_atom_name + space_to_residue_name + curr_residue_name

				call(["sed -i \"s/" + old_line + "/" + new_line + "/g\" " + new_file_name], shell=True)

		elif curr_residue_name == "PHE":

			CD1_line = residue[4] # CA, C, CB, CG, then CD1
			CD1_line_arr = CD1_line.split()
			CD1_num = CD1_line_arr[1] # first word is CONECT

			unknown_lines_beg = 5
			atom_mapping = {}

			# finding the C atom that connects to CD1_num must be CE1
			# finding the C atom remaining that connects to CE1 must be CZ
			# finding the C atom remaining that connects to CZ must be CE2
			# finding the C atom remaining that connects to CE2 must be CD2
			for mode in range(4):
				for line_index in range(unknown_lines_beg, len(residue)):
					c_line = curr_conect_line[line_index]
					a_line = residue[line_index]
					a_line_arr = a_line.split()
					curr_atom = a_line_arr[1]
					elem = a_line_arr[2]
					c_line_arr = c_line.split()
					connected_atoms = c_line_arr[2:]

					if mode == 0:
						if elem == 'C' and CD1_num in connected_atoms: atom_mapping[curr_atom] = "CE1" 
					elif mode == 1:
						if elem == 'C' and curr_atom not in atom_mapping.keys():
							a0 = False
							if connected_atoms[0] in atom_mapping.keys():
								if atom_mapping[connected_atoms[0]] == 'CE1': a0 = True
							a1 = False
							if connected_atoms[1] in atom_mapping.keys():
								if atom_mapping[connected_atoms[1]] == 'CE1': a1 = True
							if a0 or a1: atom_mapping[curr_atom] = "CZ"
					elif mode == 2:
						if elem == 'C' and curr_atom not in atom_mapping.keys():
							a0 = False
							if connected_atoms[0] in atom_mapping.keys():
								if atom_mapping[connected_atoms[0]] == 'CZ': a0 = True
							a1 = False
							if connected_atoms[1] in atom_mapping.keys():
								if atom_mapping[connected_atoms[1]] == 'CZ': a1 = True
							if a0 or a1: atom_mapping[curr_atom] = "CE2"
					elif mode == 3:
						if elem == 'C' and curr_atom not in atom_mapping.keys():
							a0 = False
							if connected_atoms[0] in atom_mapping.keys():
								if atom_mapping[connected_atoms[0]] == 'CE2': a0 = True
							a1 = False
							if connected_atoms[1] in atom_mapping.keys():
								if atom_mapping[connected_atoms[1]] == 'CE2': a1 = True
							if a0 or a1: atom_mapping[curr_atom] = "CD2"


			atom_names = ['CA', 'C', 'CB', 'CG', 'CD1']
			name_index = 0
			for line in residue:
				old_line = line[:11] + "  " + line[13:20]

				line_arr = line.split()
				curr_atom = line_arr[1]

				if curr_atom in atom_mapping.keys():
					new_atom_name = atom_mapping[curr_atom]
				else: 
					new_atom_name = atom_names[name_index]
					name_index += 1	

				space_to_residue_name = -1
				if len(new_atom_name) == 1:   space_to_residue_name = "   "
				elif len(new_atom_name) == 2: space_to_residue_name = "  "
				elif len(new_atom_name) == 3: space_to_residue_name = " "
				new_line = line[:11] + "  " + new_atom_name + space_to_residue_name + curr_residue_name

				call(["sed -i \"s/" + old_line + "/" + new_line + "/g\" " + new_file_name], shell=True)
				
		elif curr_residue_name == "TYR":

			CD1_line = residue[4] # CA, C, CB, CG, then CD1
			CD1_line_arr = CD1_line.split()
			CD1_num = CD1_line_arr[1] # first word is CONECT

			unknown_lines_beg = 5
			atom_mapping = {}

			# finding the C atom that connects to CD1_num must be CE1
			# finding the C atom remaining that has three connections must be CZ
			# the O atom must be OH
			# finding the C atom remaining that connects to CZ must be CE2
			# finding the C atom remaining that connects to CE2 must be CD2
			for mode in range(5):
				for line_index in range(unknown_lines_beg, len(residue)):
					c_line = curr_conect_line[line_index]
					a_line = residue[line_index]
					a_line_arr = a_line.split()
					curr_atom = a_line_arr[1]
					elem = a_line_arr[2]
					c_line_arr = c_line.split()
					connected_atoms = c_line_arr[2:]

					if mode == 0:
						if elem == 'C' and CD1_num in connected_atoms: atom_mapping[curr_atom] = "CE1" 
					elif mode == 1:
						if elem == 'C' and len(connected_atoms) == 3: atom_mapping[curr_atom] = "CZ"
					elif mode == 2:
						if elem == 'O': atom_mapping[curr_atom] = "OH"
					elif mode == 3:
						if elem == 'C' and curr_atom not in atom_mapping.keys() and len(connected_atoms) == 2:
							a0 = False
							if connected_atoms[0] in atom_mapping.keys():
								if atom_mapping[connected_atoms[0]] == 'CZ': a0 = True
							a1 = False
							if connected_atoms[1] in atom_mapping.keys():
								if atom_mapping[connected_atoms[1]] == 'CZ': a1 = True
							if a0 or a1: atom_mapping[curr_atom] = "CE2"
					elif mode == 4:
						if elem == 'C' and curr_atom not in atom_mapping.keys() and len(connected_atoms) == 2:
							a0 = False
							if connected_atoms[0] in atom_mapping.keys():
								if atom_mapping[connected_atoms[0]] == 'CE2': a0 = True
							a1 = False
							if connected_atoms[1] in atom_mapping.keys():
								if atom_mapping[connected_atoms[1]] == 'CE2': a1 = True
							if a0 or a1: atom_mapping[curr_atom] = "CD2"


			atom_names = ['CA', 'C', 'CB', 'CG', 'CD1']
			name_index = 0
			for line in residue:
				old_line = line[:11] + "  " + line[13:20]

				line_arr = line.split()
				curr_atom = line_arr[1]

				if curr_atom in atom_mapping.keys():
					new_atom_name = atom_mapping[curr_atom]
				else: 
					new_atom_name = atom_names[name_index]
					name_index += 1	

				space_to_residue_name = -1
				if len(new_atom_name) == 1:   space_to_residue_name = "   "
				elif len(new_atom_name) == 2: space_to_residue_name = "  "
				elif len(new_atom_name) == 3: space_to_residue_name = " "
				new_line = line[:11] + "  " + new_atom_name + space_to_residue_name + curr_residue_name

				call(["sed -i \"s/" + old_line + "/" + new_line + "/g\" " + new_file_name], shell=True)


		elif curr_residue_name == "TRP":

			CG_line = residue[3] # CA, C, CB, then CG
			CG_line_arr = CG_line.split()
			CG_num = CG_line_arr[1] # first word is CONECT

			unknown_lines_beg = 4
			atom_mapping = {}

			# 0 the N atom is NE1
			# 1 the C atom that connects to CG_num and connects to NE1 must be CD1
			# 2 the C atom remaining that has three connections and connects to NE1 must be CE2
			# 3 the C atom remaining that has three connections and connects to CE2 must be CD2
			# 4 the C atom remaining that connects to CD2 must be CE3
			# 5 the C atom reminaing that connects to CE3 must be CZ3
			# 6 the C atom remaining that connects to CZ3 must be CH2
			# 7 the C atom remaining that connects to CH2 must be CZ2

			for mode in range(8):
				for line_index in range(unknown_lines_beg, len(residue)):
					c_line = curr_conect_line[line_index]
					a_line = residue[line_index]
					a_line_arr = a_line.split()
					curr_atom = a_line_arr[1]
					elem = a_line_arr[2]
					c_line_arr = c_line.split()
					connected_atoms = c_line_arr[2:]

					if mode == 0:
						if elem == 'N': atom_mapping[curr_atom] = "NE1" 
					elif mode == 1:
						if elem == 'C' and CG_num in connected_atoms and len(connected_atoms) == 2: atom_mapping[curr_atom] = "CD1"
					elif mode == 2:
						if elem == 'C' and curr_atom not in atom_mapping.keys() and len(connected_atoms) == 3:
							a0 = False
							if connected_atoms[0] in atom_mapping.keys():
								if atom_mapping[connected_atoms[0]] == 'NE1': a0 = True
							a1 = False
							if connected_atoms[1] in atom_mapping.keys():
								if atom_mapping[connected_atoms[1]] == 'NE1': a1 = True
							a2 = False
							if connected_atoms[2] in atom_mapping.keys():
								if atom_mapping[connected_atoms[2]] == 'NE1': a2 = True
							if a0 or a1 or a2: atom_mapping[curr_atom] = "CE2"
					elif mode == 3:
						if elem == 'C' and curr_atom not in atom_mapping.keys() and len(connected_atoms) == 3:
							a0 = False
							if connected_atoms[0] in atom_mapping.keys():
								if atom_mapping[connected_atoms[0]] == 'CE2': a0 = True
							a1 = False
							if connected_atoms[1] in atom_mapping.keys():
								if atom_mapping[connected_atoms[1]] == 'CE2': a1 = True
							a2 = False
							if connected_atoms[2] in atom_mapping.keys():
								if atom_mapping[connected_atoms[2]] == 'CE2': a2 = True
							if a0 or a1 or a2: atom_mapping[curr_atom] = "CD2"
					elif mode == 4:
						if elem == 'C' and curr_atom not in atom_mapping.keys() and len(connected_atoms) == 2:
							a0 = False
							if connected_atoms[0] in atom_mapping.keys():
								if atom_mapping[connected_atoms[0]] == 'CD2': a0 = True
							a1 = False
							if connected_atoms[1] in atom_mapping.keys():
								if atom_mapping[connected_atoms[1]] == 'CD2': a1 = True
							if a0 or a1: atom_mapping[curr_atom] = "CE3"
					elif mode == 5:
						if elem == 'C' and curr_atom not in atom_mapping.keys() and len(connected_atoms) == 2:
							a0 = False
							if connected_atoms[0] in atom_mapping.keys():
								if atom_mapping[connected_atoms[0]] == 'CE3': a0 = True
							a1 = False
							if connected_atoms[1] in atom_mapping.keys():
								if atom_mapping[connected_atoms[1]] == 'CE3': a1 = True
							if a0 or a1: atom_mapping[curr_atom] = "CZ3"
					elif mode == 6:
						if elem == 'C' and curr_atom not in atom_mapping.keys() and len(connected_atoms) == 2:
							a0 = False
							if connected_atoms[0] in atom_mapping.keys():
								if atom_mapping[connected_atoms[0]] == 'CZ3': a0 = True
							a1 = False
							if connected_atoms[1] in atom_mapping.keys():
								if atom_mapping[connected_atoms[1]] == 'CZ3': a1 = True
							if a0 or a1: atom_mapping[curr_atom] = "CH2"
					elif mode == 7:
						if elem == 'C' and curr_atom not in atom_mapping.keys() and len(connected_atoms) == 2:
							a0 = False
							if connected_atoms[0] in atom_mapping.keys():
								if atom_mapping[connected_atoms[0]] == 'CH2': a0 = True
							a1 = False
							if connected_atoms[1] in atom_mapping.keys():
								if atom_mapping[connected_atoms[1]] == 'CH2': a1 = True
							if a0 or a1: atom_mapping[curr_atom] = "CZ2"


			atom_names = ['CA', 'C', 'CB', 'CG']
			name_index = 0
			for line in residue:
				old_line = line[:11] + "  " + line[13:20]

				line_arr = line.split()
				curr_atom = line_arr[1]

				if curr_atom in atom_mapping.keys():
					new_atom_name = atom_mapping[curr_atom]
				else: 
					new_atom_name = atom_names[name_index]
					name_index += 1	

				space_to_residue_name = -1
				if len(new_atom_name) == 1:   space_to_residue_name = "   "
				elif len(new_atom_name) == 2: space_to_residue_name = "  "
				elif len(new_atom_name) == 3: space_to_residue_name = " "
				new_line = line[:11] + "  " + new_atom_name + space_to_residue_name + curr_residue_name

				call(["sed -i \"s/" + old_line + "/" + new_line + "/g\" " + new_file_name], shell=True)



call(["rm " + new_connect_file_name], shell=True)

import mdtraj as md

orig = md.load_pdb("receptor.pdb")

temp = md.load_pdb(new_file_name)
#orig_atoms = [a for a in orig.top.atoms]
temp_atoms = [a for a in temp.top.atoms]

for i, temp_atom in enumerate(temp_atoms):
    dash_index = str(temp_atom).find('-')
    residue_name = str(temp_atom)[:dash_index]
    element_type = str(temp_atom)[dash_index+1:]
    residue_index = int(residue_name[3:]) - 1
    
    orig_atoms_in_residue_indices = orig.top.select("resid == " + str(residue_index) + " and name == " + element_type)
    orig.xyz[0, orig_atoms_in_residue_indices[0], :] = temp.xyz[0, i, :]


orig.save_pdb(sys.argv[1] + ".complete")

sys.exit(0)





