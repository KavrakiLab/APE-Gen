from subprocess import call
import time
import sys
import os

from Bio import pairwise2
from Bio import SeqIO
from Bio.PDB import *
import re
import argparse

import bs4
import requests

try:
    from modeller import *
    from modeller.automodel import *
except:
    print("Error with importing Modeller: Make sure license key is correct.")
    sys.exit(0)

defaults_location = os.path.dirname(os.path.abspath(__file__))

# python model_receptor.py < fasta file containing alpha chain> <template pdb>

def main(args):

    parser = argparse.ArgumentParser(description="Homology Modeling of HLAs using Modeller", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('alpha_chain_seq', type=str, nargs=1, help='Fasta file (.fasta) containing the sequence of the alpha chain of HLA or name of HLA allele (ex. HLA-A*02:01). If allele name given, the program will try to download the sequence from EMBL-EBI.')
    parser.add_argument('template', type=str, nargs=1, help='PDB of the template HLA or name of HLA allele (ex. HLA-A*02:01). If allele name given, a template based on the allele\'s supertype (as defined in supertype_templates.csv) will be chosen.')

    parser.add_argument("-n", "--num_models", type=int, default=10, help='Number of models to sample with Modeller')

    args = parser.parse_args(args)

    alpha_chain_seq = args.alpha_chain_seq[0]
    template = args.template[0]

    num_models = args.num_models

    if template[-4:] == ".pdb": 
        print("Using user-inputted template PDB")
        template_pdb = template
    else:
        allele_name = template
        print("Allele to model:", allele_name)
        allele_name = allele_name.replace(":", "")

        supertypes = {}
        representatives = {}
        f = open(defaults_location + "/supertype_templates.csv", 'r')
        for line in f:
            allele, supertype, isRep, pdbid = line.split(",")
            pdbid = pdbid[:4] # remove newline
            if supertype not in supertypes.keys():
                supertypes[supertype] = []
            supertypes[supertype].append(allele)
            if isRep == "representative":
                representatives[supertype] = (allele, pdbid)
        f.close()

        #for k in supertypes.keys():
        #    if k != "Unclassified":
        #        print(k, len(supertypes[k]), representatives[k])

        foundAllele = False
        for s in supertypes.keys():
            if allele_name[4:] in supertypes[s]:
                foundAllele = True
                if s == "Unclassified":
                    print("Allele is Unclassified: Using a template from the same serotype")
                    if allele_name[4] == "A":
                        template_allele, template_pdbcode = representatives["A02"]
                    elif allele_name[4] == "B":
                        template_allele, template_pdbcode = representatives["B07"]
                    elif allele_name[4] == "C":
                        template_allele, template_pdbcode = representatives["C00"]
                    else:
                        print("Error: Unrecognized Serotype")
                        sys.exit(0)
                else:
                    print("Obtaining representative template from within supertype:")
                    template_allele, template_pdbcode = representatives[s]
                    print("Template allele, pdbcode:", template_allele, template_pdbcode)
                break
            elif allele_name[4] == "C":
                foundAllele = True
                print("This is an HLA-C allele: Using representative template")
                template_allele, template_pdbcode = representatives["C00"]
                break
                    
        if not foundAllele:
            print("Warning: Allele cannot be found in the internal database ... Defaulting to 2v2w.pdb")
            template_pdbcode = "2v2w"

        print("Downloading " + template_pdbcode)
        call(["python " + defaults_location + "/get_pMHC_pdb.py " + template_pdbcode], shell=True)
        template_pdb = template_pdbcode + ".pdb"    

    print("Removing signal and transmembrane portions of alpha chain seq")

    if alpha_chain_seq[-6:] != ".fasta": # attempt to download

        print("Attempting to Download sequence of allele", alpha_chain_seq)
        page = requests.get("https://www.ebi.ac.uk/cgi-bin/ipd/imgt/hla/get_allele.cgi?" + alpha_chain_seq[4:])
        soup = bs4.BeautifulSoup(page.content, 'lxml')

        list_of_sequences = soup.findAll('pre')
        if len(list_of_sequences) == 0:
            print("Error: Allele could not be found in database")
            sys.exit(0)

        raw_str = str(soup.findAll('pre')[0])
        raw_str = raw_str.replace("\n", "")
        raw_str = raw_str.replace("<pre>", "")
        raw_seq = raw_str.replace("</pre>", "")

    else: 
        alpha_chain_seq_file = alpha_chain_seq
        raw_seq = str(list(SeqIO.parse(alpha_chain_seq_file, "fasta"))[0].seq)
    
    parser = PDBParser()
    structure = parser.get_structure('asfd', template_pdb)
    ppb = PPBuilder()
    for pp in ppb.build_peptides(structure):
        alpha_chain = str(pp.get_sequence())
        break

    """
    alignments = pairwise2.align.localxx(raw_seq, alpha_chain, penalize_end_gaps=False)
    print(alignments[0])
    #sys.exit(0)

    print(alpha_chain[:10])
    print(alpha_chain[-5:])
    print(raw_seq)
    
    print(raw_seq[24:30])
    sys.exit(0)

    if alpha_chain[0] != "M":
        beg_index = [m.start() for m in re.finditer(alpha_chain[:5], raw_seq)][0]
    else:
        beg_index = [m.start() for m in re.finditer(alpha_chain[1:6], raw_seq)][0]
    """

    beg_index = 24 # The signaling portion in the beginning appear to all be the same length
    end_index = [m.start() for m in re.finditer(alpha_chain[-5:], raw_seq)][-1]

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
