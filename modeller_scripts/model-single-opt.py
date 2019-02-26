# -*- coding: utf-8 -*-
# File: model-single-opt.py
# Reading the ali file and generating 1 models
# Example of changing the default optmisation schedule

from modeller import *
from modeller.automodel import *
import sys

log.verbose()
env = environ()

# Give less weight to all soft-sphere restraints:
env.schedule_scale = physical.values(default=1.0, soft_sphere=0.7)

#Considering heteroatoms and waters molecules
env.io.hetatm = env.io.water = True
# Directories with input atom files:
env.io.atom_files_directory = './:../atom_files'


# Modelling 'sequence' with file.ali
a = automodel(env, alnfile='target_sequence-receptor_template.ali',
              knowns='receptor_template',
              sequence='target_sequence',
              assess_methods=(assess.DOPE, assess.GA341)
              )
# Generating 3 models
a.starting_model = 1
a.ending_model = int(sys.argv[1])

# Very thorough Variable Target Function Method (VTFM) optimization:
a.library_schedule = autosched.slow
a.max_var_iterations = 300

# Thorough MD optimization:
a.md_level = refine.slow
 
# Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
a.repeat_optimization = 2
a.max_molpdf = 1e6

a.make()

# Get a list of all successfully built models from a.outputs
ok_models = filter(lambda x: x['failure'] is None, a.outputs)

# Rank the models by DOPE score
key = 'DOPE score'
ok_models.sort(lambda a,b: cmp(a[key], b[key]))

# Get top model
m = ok_models[0]
print "Top model: %s (DOPE score %.100f)" % (m['name'], m[key])
