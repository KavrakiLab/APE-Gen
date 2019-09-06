from subprocess import call
import time
import sys
import os
from modeller import *
from modeller.automodel import *
from Bio import pairwise2
from Bio import SeqIO
from Bio.PDB import *
import re
import argparse

defaults_location = os.path.dirname(os.path.abspath(__file__))

# python model_receptor.py < fasta file containing alpha chain> <template pdb>

def main(args):

    parser = argparse.ArgumentParser(description="Anchored Peptide-MHC Ensemble Generator", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('alpha_chain_seq_file', type=str, nargs=1, help='Fasta file containing the sequence of the alpha chain of HLA')
    parser.add_argument('template_pdb', type=str, nargs=1, help='PDB of the template HLA')

    parser.add_argument("-n", "--num_models", type=int, default=10, help='Number of models to sample with Modeller')

    args = parser.parse_args(args)

    alpha_chain_seq_file = args.alpha_chain_seq_file[0]
    template_pdb = args.template_pdb[0]

    num_models = args.num_models
    

    print("Removing signal and transmembrane portions of alpha chain seq")

    raw_seq = str(list(SeqIO.parse(alpha_chain_seq_file, "fasta"))[0].seq)
    
    parser = PDBParser()
    structure = parser.get_structure('asfd', template_pdb)
    ppb = PPBuilder()
    for pp in ppb.build_peptides(structure):
        alpha_chain = str(pp.get_sequence())
        break

    beg_index = [m.start() for m in re.finditer(alpha_chain[:5], raw_seq)][0]
    end_index = [m.start() for m in re.finditer(alpha_chain[-5:], raw_seq)][-1]

    #alignments = pairwise2.align.localxx(raw_seq[beg_index:end_index+1], alpha_chain)
    #print(alignments[0])

    target_sequence = raw_seq[beg_index:end_index+1]

    print("Length of original alpha chain seq:", len(raw_seq))
    print("Length of processed alpha chain seq:", len(target_sequence))
    print("Length of alpha chain in template:", len(alpha_chain))

    """
    f = open(args[0], 'r')
    f_str = ""
    for line in f:
        f_str += line
    f.close()
    target_sequence = "".join(f_str.split())
    """

    target_sequence += "/"
    # constant beta immunoglobin
    target_sequence += "MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM*" 
    receptor_template = template_pdb

    #log.verbose()
    #log.minimal()
    log.none()

    print("Preparing receptor template")
    #call(["cp " + receptor_template + " receptor_template.pdb"], shell=True)
    #call(["sed -i \"/[A-Z] C  /d\" receptor_template.pdb"], shell=True)
    #call(["python " + defaults_location + "/mod_getseq.py > log.txt"], shell=True) # creates receptor_template.seq
    mod_getseq(receptor_template)

    print("Preparing target sequence")
    f = open("target_sequence.pir", 'w')
    f.write(">P1;target_sequence\n")
    f.write("sequence::     : :     : :::-1.00:-1.00\n")
    f.write(target_sequence)
    f.close()

    print("Aligning target sequence with receptor template")
    starttime = time.time()
    #call(["python " + defaults_location + "/align2d.py >> log.txt"], shell=True) # creates target_sequence-receptor_template.ali
    align_2d()
    endtime = time.time()
    print("Aligning took " + str(endtime - starttime) + " seconds.")

    print("Creating model")
    starttime = time.time()
    #call(["python -u  " + defaults_location + "/model-single-opt.py 2 >> log.txt"], shell=True)
    #call(["cp $(grep \"Top model\" log.txt | awk '{ print $3 }') best_model.pdb"], shell=True)
    best_model, best_score = model_single_opt(num_models)
    call(["cp " + str(best_model) + " best_model.pdb"], shell=True)
    endtime = time.time()
    print("Homology modelling took " + str(endtime - starttime) + " seconds.")

def model_single_opt(num_models):

    #log.verbose()
    #log.minimal()
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
    a.ending_model = int(num_models)

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
    ok_models = [x for x in a.outputs if x['failure'] is None]

    # Rank the models by DOPE score
    key = 'DOPE score'
    #print(ok_models)
    ok_models.sort(key=lambda a: a[key])

    # Get top model
    m = ok_models[0]
    print("Top model: %s (DOPE score %.100f)" % (m['name'], m[key]))

    return m['name'], m[key]

def align_2d():

    env = environ()
    aln = alignment(env)
    mdl = model(env, file='receptor_template') #, model_segment=('FIRST:A','LAST:C'))
    aln.append_model(mdl, align_codes='receptor_template', atom_files='receptor_template.pdb')
    aln.append(file='target_sequence.pir', align_codes='target_sequence')
    aln.align2d()
    aln.write(file='target_sequence-receptor_template.ali', alignment_format='PIR')
    #aln.write(file='TvLDH-1bdmA.pap', alignment_format='PAP')

def mod_getseq(receptor_template):

    # Get the sequence of the PDB file, and write to an alignment file
    call(["cp " + receptor_template + " receptor_template.pdb"], shell=True)
    call(["sed -i \"/[A-Z] C  /d\" receptor_template.pdb"], shell=True)
    code = 'receptor_template'
    e = environ()
    m = model(e, file=code)
    aln = alignment(e)
    aln.append_model(m, align_codes=code)
    aln.write(file=code+'.seq')


if __name__ == "__main__":
    main(sys.argv[1:])
