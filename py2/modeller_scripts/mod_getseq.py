from modeller import *
# Get the sequence of the PDB file, and write to an alignment file
code = 'receptor_template'
e = environ()
m = model(e, file=code)
aln = alignment(e)
aln.append_model(m, align_codes=code)
aln.write(file=code+'.seq')
