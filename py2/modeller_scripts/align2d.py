from modeller import *

env = environ()
aln = alignment(env)
mdl = model(env, file='receptor_template') #, model_segment=('FIRST:A','LAST:C'))
aln.append_model(mdl, align_codes='receptor_template', atom_files='receptor_template.pdb')
aln.append(file='target_sequence.pir', align_codes='target_sequence')
aln.align2d()
aln.write(file='target_sequence-receptor_template.ali', alignment_format='PIR')
#aln.write(file='TvLDH-1bdmA.pap', alignment_format='PAP')