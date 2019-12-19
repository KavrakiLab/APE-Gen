"""
* Software License Agreement (BSD License)
*
*  Copyright (c) 2019, Rice University.
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of the Rice University nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*********************************************************************/

Author: Jayvee Abella
"""

import sys
import os
from subprocess import call
from subprocess import check_output
from threading import Thread
import numpy as np
import mdtraj as md
import argparse
import glob
import uuid

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

one_letter_code = {'ARG':'R', 'HIS':'H', 'LYS':'K', 'ASP':'D', 'GLU':'E', \
                  'SER':'S', 'THR':'T', 'ASN':'N', 'GLN':'Q', 'CYS':'C', \
                  'GLY':'G', 'PRO':'P', 'ALA':'A', 'VAL':'V', 'ILE':'I', \
                  'LEU':'L', 'MET':'M', 'PHE':'F', 'TYR':'Y', 'TRP':'W'}

three_letter_code  = {v: k for k, v in one_letter_code.items()}

defaults_location = os.path.dirname(os.path.abspath(__file__)) #sys.path[0]
pymol_location = "pymol"
RCD_location = "rcd"
smina_location = "smina"
vina_location = "vina_split"
pdbfixer_location = "pdbfixer"

def main(args):

    parser = argparse.ArgumentParser(description="Anchored Peptide-MHC Ensemble Generator", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('peptide_input', type=str, nargs=1, help='Sequence of peptide to dock or pdbfile of crystal structure')
    parser.add_argument('receptor_class', type=str, nargs=1, help='Class descriptor of MHC receptor. Use REDOCK along with crystal input to perform redocking. Or pass a PDB file with receptor')

    parser.add_argument("-n", "--num_cores", type=int, default=8, help='Number of cores to use for RCD and smina computations.')
    parser.add_argument("-l", "--num_loops", type=int, default=100, help='Number of loops to generate with RCD. (Note that the final number of sampled conformations may be less due to steric clashes.')
    parser.add_argument("-t", "--RCD_dist_tol", type=float, default=1.0, help='RCD tolerance (in angstroms) of inner residues when performing IK')
    parser.add_argument("-r", "--rigid_receptor", action="store_true", help='Disable sampling of receptor degrees of freedom specified in flex_res.txt')
    parser.add_argument("-d", "--debug", action="store_true", help='Print extra information for debugging')
    parser.add_argument("-p", "--save_only_pep_confs", action="store_true", help='Disable saving full conformations (peptide and MHC)')
    #parser.add_argument("--n-mer-templates", default="", help='File with n-mer pdb templates.')
    #parser.add_argument("--receptor-class-templates", default="", help='File with pdb receptor class templates')
    #parser.add_argument("--flex_res", default="", help='File with flexible residues')
    parser.add_argument("-a", "--anchor_tol", type=float, default=2.0, help='Anchor tolerance (in angstroms) of first and last backbone atoms of peptide when filtering')
    parser.add_argument("-o", "--score_with_openmm", action="store_true", help='Rescore full conformations with openmm (AMBER)')
    parser.add_argument("-g", "--num_rounds", type=int, default=1, help='Number of rounds to perform.')
    parser.add_argument("-b", "--pass_type", type=str, default='receptor_only', choices=['receptor_only', 'pep_and_recept'], help="When using multiple rounds, pass best scoring conformation across different rounds (choose either 'receptor_only' or 'pep_and_recept')")
    parser.add_argument("-s", "--min_with_smina", action="store_true", help='Minimize with SMINA instead of the default Vinardo')
    parser.add_argument("--use_gpu", action="store_true", help='Use GPU for OpenMM Minimization step')
    parser.add_argument("--no_progress", action="store_true", help='Do not print progress bar')
    parser.add_argument("--clean_rcd", action="store_true", help='Remove RCD folder at the end of each round')

    args = parser.parse_args(args)

    peptide_input = args.peptide_input[0]
    receptor_class = args.receptor_class[0]

    num_cores = args.num_cores
    num_loops = args.num_loops
    RCD_dist_tol = args.RCD_dist_tol
    doReceptorMinimization = not args.rigid_receptor
    debug = args.debug
    saveFullConfs = not args.save_only_pep_confs
    anchor_tol = args.anchor_tol
    score_with_openmm = args.score_with_openmm
    num_rounds = args.num_rounds
    pass_type = args.pass_type
    min_with_smina = args.min_with_smina
    use_gpu = args.use_gpu
    if use_gpu: device = "OpenCL"
    else: device = "CPU"
    no_progress = args.no_progress
    printProgress = not no_progress
    cleanRCD = args.clean_rcd

    temp_location = "/dev/shm/" + uuid.uuid4().hex + "/"

    current_working_dir = os.getcwd()
    call(["mkdir -p " + temp_location], shell=True)
    os.chdir(temp_location)
    call(["cp " + current_working_dir + "/* ."], shell=True)

    print("Preparing peptide and MHC")

    native_loc = None
    performRedock = False
    # if input is PDB, then use it to define 'native' for later comparison
    if peptide_input[-4:] == ".pdb":

        native_pdb = peptide_input

        call(["grep \"[A-Z] C  \" " + native_pdb + " > native.pdb"], shell=True)
        native_loc = "native.pdb"
        f = md.load("native.pdb")
        seq_arr = [one_letter_code[r.name] for r in f.top.residues]
        seq = ''
        for s in seq_arr: seq += s
        peptide_sequence = seq
        print("Peptide sequence: " + peptide_sequence)
        print("Receptor class: " + receptor_class)

        if receptor_class == "REDOCK":
            performRedock = True
            receptor_template = native_pdb

    # otherwise, only sequence is needed
    else:
        peptide_sequence = peptide_input
        for resi in peptide_sequence:
            if resi not in list(three_letter_code.keys()):
                print("Unrecognized one letter residue code: " + resi)
                sys.exit(0)
        print("Peptide sequence: " + peptide_sequence)
        print("Receptor class: " + receptor_class)

        if receptor_class == "REDOCK":
            print("Must input native pdb as peptide input in order to perform redocking.")
            sys.exit(0)

    # retrieve peptide template
    n_mer_templates = {}
    f = open(defaults_location + "/n-mer-templates.txt")
    for line in f:
        line_arr = line.split()
        n_mer_templates[line_arr[0]] = line_arr[1]
    f.close()
    if str(len(peptide_sequence)) in list(n_mer_templates.keys()):
        peptide_template = n_mer_templates[str(len(peptide_sequence))]
    else:
        print(str(len(peptide_sequence)) + "-mers are not implemented yet. Check n-mer-templates.txt")
        sys.exit(0)

    call(["cp " + defaults_location + "/templates/" + peptide_template + " ."], shell=True) 

    # if pdb is provided, use MHC inside
    # elseif receptor allotype is specified, read in receptor templates
    if receptor_class[-4:] == ".pdb": receptor_template = receptor_class
    elif not performRedock:
        receptor_class_templates = {}
        f = open(defaults_location + "/receptor-class-templates.txt")
        for line in f:
            line_arr = line.split()
            receptor_class_templates[line_arr[0]] = line_arr[1]
        f.close()
        if receptor_class in list(receptor_class_templates.keys()):
            receptor_template = receptor_class_templates[receptor_class]
        else:
            print(receptor_class + " not yet implemented as a receptor class. Check receptor-class-templates.txt")
            sys.exit(0)
        call(["cp " + defaults_location + "/templates/" + receptor_template + " ."], shell=True)

    # read in flexible residues
    flexible_residues = ""
    f = open(defaults_location + "/flex_res.txt")
    for line in f:
        flexible_residues += line
    f.close()
    flexible_residues = flexible_residues.rstrip()

    for current_round in range(num_rounds): 

        if not os.path.exists(str(current_round)): call(["mkdir " + str(current_round)], shell=True)
        os.chdir(str(current_round))

        if current_round == 0:
            call(["cp ../*.pdb ."], shell=True)
        else:
            call(["cp ../*.pdb ."], shell=True)
            if pass_type == 'receptor_only':
                call(["cp ../" + str(current_round-1) + "/min_energy_system.pdb ."], shell=True)
                receptor_template = "min_energy_system.pdb"
            elif pass_type == 'pep_and_recept':
                call(["cp ../" + str(current_round-1) + "/min_energy_system.pdb ."], shell=True)
                receptor_template = "min_energy_system.pdb"
                peptide_template = "min_energy_system.pdb"
            call(["sed -i \"s/          /           /g\" min_energy_system.pdb"], shell=True) # fix annoying openbabel warning


        print("Aligning peptide anchors to MHC pockets")

        call([pymol_location + " -qc " + defaults_location + "/align.py " + peptide_template + " " + receptor_template + " >> align.log 2>&1"], shell=True)

        if native_loc != None:
            call([pymol_location + " -qc " + defaults_location + "/align.py " + native_pdb + " " + receptor_template + " >> align.log 2>&1"], shell=True)
            call(["grep \"[A-Z] C  \" aln-" + native_pdb + " > native.pdb"], shell=True)
            native_loc = "native.pdb"

        call(["cp " + receptor_template + " ./receptor.pdb; sed -i \"/[A-Z] C/d\" receptor.pdb"], shell=True)

        call(["grep \"[A-Z] C  \" aln-" + peptide_template + " > peptide_template.pdb"], shell=True)
        f = md.load("peptide_template.pdb")
        seq_arr = [r.name for r in f.top.residues]
        call(["sed -i \"s/" + seq_arr[0] + " C   1/" + three_letter_code[peptide_sequence[0]] + " C   1/g\" peptide_template.pdb"], shell=True)
        call(["sed -i \"s/" + seq_arr[1] + " C   2/" + three_letter_code[peptide_sequence[1]] + " C   2/g\" peptide_template.pdb"], shell=True)

        if len(peptide_sequence) < 10:
            last_anchor_resi_1 = " C   " + str(len(peptide_sequence) - 1)
            last_anchor_resi_2 = " C   " + str(len(peptide_sequence))
        elif len(peptide_sequence) == 10:
            last_anchor_resi_1 = " C   " + str(len(peptide_sequence) - 1)
            last_anchor_resi_2 = " C  " + str(len(peptide_sequence))
        elif len(peptide_sequence) > 10:
            last_anchor_resi_1 = " C  " + str(len(peptide_sequence) - 1)
            last_anchor_resi_2 = " C  " + str(len(peptide_sequence))

        call(["sed -i \"s/" + seq_arr[-2] + last_anchor_resi_1 + "/" + three_letter_code[peptide_sequence[-2]] + last_anchor_resi_1 + "/g\" peptide_template.pdb"], shell=True)
        call(["sed -i \"s/" + seq_arr[-1] + last_anchor_resi_2 + "/" + three_letter_code[peptide_sequence[-1]] + last_anchor_resi_2 + "/g\" peptide_template.pdb"], shell=True)
        if len(peptide_sequence) < 12:
            call(["sed -i \"/C   [3-" + str(len(peptide_sequence) - 2) + "]/d\" peptide_template.pdb"], shell=True)
        else:
            call(["sed -i \"/C   [3-9]/d\" peptide_template.pdb"], shell=True)
            for j in range(len(peptide_sequence)-1 - 10):
                call(["sed -i \"/C  1" + str(j) + "/d\" peptide_template.pdb"], shell=True)

        f = open("peptide_template.pdb", 'r')
        f_new = open("anchors.pdb", "w")
        for line in f:
            firstword = line.split()[0]
            elem = line.split()[2]
            resi = line.split()[3]
            if firstword == 'ATOM':
                if resi != 'GLY':
                    if elem in ['N', 'CA', 'C', 'O', 'CB']:
                        f_new.write(line)
                else:
                    if elem in ['N', 'CA', 'C', 'O']:
                        f_new.write(line)

        f.close()
        f_new.close()
        call(["cp anchors.pdb peptide.pdb"], shell=True)

        print("Sampling peptide backbone")
        if os.path.exists("RCD"):
            print("Found RCD folder. Skipping this step")
            ref = md.load("confs_top.pdb")
            ref_top = ref.top
            if native_loc != None:
                print("Native crystal found!")
                native = get_conf(native_loc, ref_top, "element != H", debug)
                print(native)
            else: native = None
            call(["cp confs_top.pdb peptide.pdb"], shell=True)
            call(["sed -i \"s/          /           /g\" peptide.pdb"], shell=True) # fix annoying openbabel warning

        else:
            call(["mkdir -p RCD/input"], shell=True)
            os.chdir("RCD/input")
            
            
            call(["cp ../../receptor.pdb anchored_pMHC.pdb"], shell=True)
            call(["cat ../../anchors.pdb >> anchored_pMHC.pdb"], shell=True)

            call(["cp " + defaults_location + "/loco.score " + defaults_location + "/dunbrack.bin ."], shell=True)

            last_loop_residue_1index = len(peptide_sequence) - 2
            call(["echo \"anchored_pMHC.pdb 3 " + str(last_loop_residue_1index) + " C " + peptide_sequence[2:last_loop_residue_1index] + "\" > loops.txt"], shell=True)

            #call([mpi_location + " " + str(num_cores) + " " + RCD_location + " -x dunbrack.bin --loco loco.score -o RCD -d " + str(RCD_dist_tol) + " -n " + str(num_loops) + " loops.txt"], shell=True)
            call([RCD_location + " -e 1 -x dunbrack.bin --energy_file loco.score -o RCD -d " + str(RCD_dist_tol) + " -n " + str(num_loops) + " loops.txt  >> rcd.log 2>&1"], shell=True)


            # ----------------------------------------
            print("Organizing RCD results")

            call(["mkdir models; cp RCD/anchored_pMHC_closed.pdb models/"], shell=True)        
            os.chdir("models")
            
            call([vina_location + " --input anchored_pMHC_closed.pdb --ligand partial >> vina.log 2>&1"], shell=True)

            array_splits = np.array_split(list(range(1, num_loops+1)), num_cores)
            folder_names = [str(s[0]) for s in array_splits]

            threads = []
            for loop_indices in array_splits:
                t = RefineThread(loop_indices, len(peptide_sequence), num_loops, doReceptorMinimization, flexible_residues, min_with_smina, debug)
                threads.append(t)
                t.start()
            
            #progress = int(check_output(["for i in " + " ".join(folder_names) + "; do grep \"MODEL\" $i/models_minimize.pdb | wc -l; done | paste -sd+ | bc"], shell=True))
            #print(progress)
            progress = 0
            while progress < num_loops:
                progress = 0
                for s in folder_names:
                    if os.path.exists(s + "/models_minimize.pdb"):
                        progress += int(check_output(["grep \"MODEL\" " + s + "/models_minimize.pdb | wc -l"], shell=True))
                #progress = int(check_output(["for i in " + " ".join(folder_names) + "; do grep \"MODEL\" $i/models_minimize.pdb | wc -l; done | paste -sd+ | bc"], shell=True))
                if progress == 0: continue
                if printProgress: printProgressBar(progress, 100, prefix = 'Progress:', suffix = 'Complete', length = 50)

            #for t in threads: t.join()

            call(["touch models_minimize.pdb"], shell=True)
            if doReceptorMinimization: call(["touch receptor_models_minimize.pdb"], shell=True)
            for folder_name in folder_names:
                call(["cat " + folder_name + "/models_minimize.pdb >> models_minimize.pdb"], shell=True)
                if doReceptorMinimization: call(["cat " + folder_name + "/receptor_new.pdb >> receptor_models_minimize.pdb"], shell=True)


            # find the minimum energy binding mode and use as reference
            f = open("models_minimize.pdb", "r")
            min_energy = 0.0
            model_index = 0
            hasHETATM = False
            min_model_index = -1
            for line in f:
                first_word = line.split()[0]
                if first_word == 'MODEL':
                    f_temp = open("temp.pdb", 'w')
                elif first_word == 'HETATM':
                    hasHETATM = True
                elif first_word == 'ATOM':
                    f_temp.write(line)
                elif first_word == 'ENDMDL':
                    model_index += 1
                    if debug: print("-------------------", model_index)
                    f_temp.close()
                    if hasHETATM:
                        hasHETATM = False
                        continue
                    else:
                        call(["sort -k6 -n temp.pdb > temp_sorted.pdb"], shell=True)
                        if energy < min_energy:
                            try:
                                potential_ref = md.load("temp_sorted.pdb")
                                resi_in_potential_ref = [r for r in potential_ref.top.residues]
                                if len(resi_in_potential_ref) != len(peptide_sequence): continue
                            except: continue
                            ref = potential_ref
                            min_energy = energy
                            min_model_index = model_index
                elif first_word == 'REMARK':
                    if line.split()[1] == 'minimizedAffinity': energy = float(line.split()[2])
            f.close()   

            print(min_energy)
            ref.save_pdb("../../../confs_top.pdb")
            ref = md.load("../../../confs_top.pdb", atom_indices=ref.top.select("element != H"))
            ref.save_pdb("../../../confs_top.pdb")
            ref_top = ref.top

            if native_loc != None:
                print("Native crystal found!")
                native = get_conf("../../../" + native_loc, ref_top, "element != H", debug)
                print(native)
            else:
                native = None

            call(["cp ../../../confs_top.pdb ../../../peptide.pdb"], shell=True)
            call(["sed -i \"s/          /           /g\" ../../../peptide.pdb"], shell=True) # fix annoying openbabel warning

            process_smina(ref, "../../../conf_data.npz", "../../../confs.dcd", native, min_model_index, debug)      

            os.chdir("..") # get out of models
            # ----------------------------------------

            os.chdir("../..") # get out of RCD/input 



       
        print("Loading sampled conformations")
        confs = md.load("confs.dcd", top="confs_top.pdb")
        energies = np.load("conf_data.npz")['energies']
        model_indices = np.load("conf_data.npz")['model_indices']
        min_model_index = np.load("conf_data.npz")['min_model_index']
        
        print("Num full confs:", len(confs))
        all_confs = confs
        all_energies = energies

        print("Saving filtered peptide confs")
        if os.path.exists("peptide_confs.pdb"):
            print("Found peptide_confs.pdb, Please move to recompute.")
        else:
            reference_bound = md.load("aln-" + peptide_template)
            ref_backbone = reference_bound.top.select("chainid == 2 and name == 'CA'")
            reference_bound = md.load("aln-" + peptide_template, atom_indices=ref_backbone)
            ref_a1 = reference_bound.xyz[0, 0, :]
            ref_a2 = reference_bound.xyz[0, -1, :]
            ref_a3 = reference_bound.xyz[0, 1, :]
            ref_a4 = reference_bound.xyz[0, -2, :]

            mdtraj_confs = []
            isFirst = True

            filtered_energies = []
            filtered_indices = []
            for i, conf in enumerate(all_confs):
                #print i
                conf.save_pdb("temp.pdb")
                sampled_conf = md.load("temp.pdb")
                sampled_backbone = sampled_conf.top.select("name == 'CA'")
                sampled_conf = md.load("temp.pdb", atom_indices=sampled_backbone)

                a1 = sampled_conf.xyz[0, 0, :]
                a2 = sampled_conf.xyz[0, -1, :]
                a3 = sampled_conf.xyz[0, 1, :]
                a4 = sampled_conf.xyz[0, -2, :]
                d_a1 = np.linalg.norm(a1-ref_a1)*10
                d_a2 = np.linalg.norm(a2-ref_a2)*10
                d_a3 = np.linalg.norm(a3-ref_a3)*10
                d_a4 = np.linalg.norm(a4-ref_a4)*10
                #print i, d_a1, d_a2, energies[i]

                if d_a1 < anchor_tol and d_a2 < anchor_tol and d_a3 < anchor_tol and d_a4 < anchor_tol:
                    #print(i)
                    filtered_energies.append(energies[i])
                    filtered_indices.append(i) # indices used for naming full system confs

                    if isFirst:
                        isFirst = False
                        mdtraj_confs = md.load("temp.pdb")
                    else:
                        mdtraj_confs += md.load("temp.pdb")

            
            print("Average filtered energy:", np.mean(filtered_energies))

            mdtraj_confs.save_pdb("peptide_confs.pdb")
            np.savez_compressed("filtered_energies.npz", filtered_energies=filtered_energies, filtered_indices=filtered_indices)


        filtered_indices = np.load("filtered_energies.npz")["filtered_indices"]
        filtered_energies = np.load("filtered_energies.npz")["filtered_energies"]
        print("Num filtered confs:", len(filtered_indices))
        print("Average filtered energy:", np.mean(filtered_energies))

        if saveFullConfs:
            print("Saving complete peptide-HLA complexes")
            if os.path.exists("full_system_confs"):
                print("Found full_system_confs/ folder. Please move to recompute.")
            else:
                if not doReceptorMinimization:
                    call(["mkdir full_system_confs"], shell=True)
                    for i, conf in enumerate(all_confs):
                        if i not in filtered_indices: continue
                        print(i)
                        conf.save_pdb("temp.pdb")
                        call(["sed -i \"s/ A  / C  /g\" temp.pdb"], shell=True)
                        call(["cat receptor.pdb temp.pdb | sed \"/MODEL/d\" | sed \"/ENDMDL/d\" | sed \"/END/d\" > system.pdb"], shell=True)
                        call(["cp system.pdb target.pdb"], shell=True)
                        #call(["echo \"14\r\n1\r\" | " + gromacs_location + " -f system.pdb -ignh -o target.pdb > /dev/null 2>&1"], shell=True)
                        call(["cp target.pdb full_system_confs/" + str(i) + ".pdb"], shell=True)
                        #if model_indices[i] == min_model_index: call(["cp target.pdb min_energy_system.pdb"], shell=True)
                    

                else:
                    """ # just use original receptor template as receptor conformation
                    call(["mkdir full_system_confs"], shell=True)
                    for i, conf in enumerate(all_confs):
                        print i
                        conf.save_pdb("temp.pdb")
                        call(["sed -i \"s/ A  / C  /g\" temp.pdb"], shell=True)
                        call(["cat receptor.pdb temp.pdb | sed \"/MODEL/d\" | sed \"/ENDMDL/d\" | sed \"/END/d\" > system.pdb"], shell=True)
                        call(["cp system.pdb target.pdb"], shell=True)
                        #call(["echo \"14\r\n1\r\" | " + gromacs_location + " -f system.pdb -ignh -o target.pdb > /dev/null 2>&1"], shell=True)
                        call(["cp target.pdb full_system_confs/" + str(i) + ".pdb"], shell=True)
                        if model_indices[i] == min_model_index: call(["cp target.pdb min_energy_system.pdb"], shell=True)
                    """
                    
                    os.chdir("RCD/input/models")
                    call([vina_location + " --input receptor_models_minimize.pdb --ligand receptor  >> vina.log 2>&1"], shell=True)
                    os.chdir("../../..")
                    call(["mkdir full_system_confs"], shell=True)

                    for j, conf in enumerate(all_confs):
                        #if printProgress: printProgressBar(j+1, len(all_confs), prefix = 'Progress:', suffix = 'Complete', length = 50)
                        
                        if j not in filtered_indices: continue

                        peptide_j = "temp" + str(model_indices[j]).zfill( len(str(num_loops)) ) + ".pdb"
                        conf.save_pdb(peptide_j)
                        call(["sed -i \"s/ A  / C  /g\" " + peptide_j], shell=True)


                    array_splits = np.array_split(list(range(len(all_confs))), num_cores)
                    folder_names = [str(s[0]) for s in array_splits]
                    num_confs = len(all_confs)

                    threads = []
                    for loop_indices in array_splits:
                        t = ReceptorThread(loop_indices, filtered_indices, model_indices, num_loops)
                        threads.append(t)
                        t.start()
                    for t in threads: t.join()

                    """
                    for j in range(num_confs):
                        
                        if j not in filtered_indices: continue                        

                        receptor_j = "receptor" + str(model_indices[j]).zfill( len(str(num_loops)) ) + ".pdbqt"
                        call(["cp RCD/input/models/" + receptor_j + " ."], shell=True)
                        
                        call(["python " + defaults_location + "/rename_atoms.py " + receptor_j], shell=True) # get receptor_j.temp and receptor_j.complete

                        complex_j = "target" + str(model_indices[j]).zfill( len(str(num_loops)) ) + ".pdb"
                        call(["cat " + receptor_j + ".complete " + peptide_j + " | sed \"/MODEL/d\" | sed \"/ENDMDL/d\" | sed \"/END/d\" > " + complex_j], shell=True)
                        call(["cp " + complex_j + " full_system_confs/" + str(j) + ".pdb"], shell=True)
                        call(["rm " + peptide_j + " " + receptor_j + " " + complex_j + " " + receptor_j + ".temp " + receptor_j + ".complete"], shell=True)
                    """


        # this comes last because calling md.rmsd centers the coordinates (messing up the alignment)
        all_confs = md.load("peptide_confs.pdb")
        all_energies = np.load("filtered_energies.npz")["filtered_energies"]
        alpha_carbon_atoms = all_confs.top.select("name == 'CA'")
        if native != None:
            rmsd_to_native = md.rmsd(all_confs, native, 0) * 10
            print("Min RMSD to native: ", np.min(rmsd_to_native), filtered_indices[np.argmin(rmsd_to_native)])
            if saveFullConfs: call(["cp full_system_confs/" + str(filtered_indices[np.argmin(rmsd_to_native)]) + ".pdb ./minRMSD.pdb"], shell=True)
            print("Energy of MinRMSD: ", all_energies[np.argmin(rmsd_to_native)])
            print("selected binding mode to native: ", md.rmsd(ref, native, 0)[0] * 10)
            rmsds = md.rmsd(all_confs, native, 0) * 10
            rmsds_sort = rmsds[np.argsort(all_energies)]
            pk_unnorm = np.exp(-np.arange(1, len(rmsds)+1))
            pk = pk_unnorm / pk_unnorm.sum()
            print("E-RMSD: ", np.dot(rmsds_sort, pk)) # try computing for just the filtered confs?

            carmsd_to_native = md.rmsd(all_confs, native, 0, atom_indices=alpha_carbon_atoms, ref_atom_indices=alpha_carbon_atoms) * 10
            print("alpha_carbon_MinRMSD: ", np.min(carmsd_to_native), filtered_indices[np.argmin(carmsd_to_native)])
            if saveFullConfs: call(["cp full_system_confs/" + str(filtered_indices[np.argmin(carmsd_to_native)]) + ".pdb ./minCaRMSD.pdb"], shell=True)
            print("Energy of alpha_carbon_MinRMSD: ", all_energies[np.argmin(carmsd_to_native)])
            print("selected binding mode to native alpha_carbon_rmsd: ", md.rmsd(ref, native, 0, atom_indices=alpha_carbon_atoms, ref_atom_indices=alpha_carbon_atoms)[0] * 10)
            rmsds = md.rmsd(all_confs, native, 0, atom_indices=alpha_carbon_atoms, ref_atom_indices=alpha_carbon_atoms) * 10
            rmsds_sort = rmsds[np.argsort(all_energies)]
            pk_unnorm = np.exp(-np.arange(1, len(rmsds)+1))
            pk = pk_unnorm / pk_unnorm.sum()
            print("E-alpha_carbon_RMSD: ", np.dot(rmsds_sort, pk))
        print("energy of selected binding mode:", np.min(all_energies), filtered_indices[np.argmin(all_energies)])
        if saveFullConfs: call(["cp full_system_confs/" + str(filtered_indices[np.argmin(all_energies)]) + ".pdb ./min_energy_system.pdb"], shell=True)

        if score_with_openmm:

            print("Scoring/Minimizing with OpenMM ...")
            if os.path.exists("full_system_confs/openmm-minimized"):
                print("Found full_system_confs/openmm-minimized folder. Please move to recompute.")
            else:

                os.chdir("full_system_confs/")

                filenames = [str(i)+".pdb" for i in filtered_indices]

                numTries = 10

                min_filenames = []
                energies = []
                for i in range(1, len(filenames)+1):

                    complex_model = "complex-" + str(i).zfill( len(str(len(filenames))) ) + ".pdb"

                    for j in range(numTries):
                        call(["pdbfixer " + filenames[i-1] + " --output=" + complex_model], shell=True)
                        
                        #call(["python " + defaults_location + "/minimize.py " + complex_model + " min-" + complex_model + " > temp.txt"], shell=True)

                        #f = open("temp.txt", 'r')
                        #for line in f:
                        #    line_arr = line.split()
                        #    if line_arr[0] == "total:": energy = float(line_arr[1]) # should go through here twice, second appearance is kept
                        #f.close()

                        energy = minimizeConf(complex_model, "min-" + complex_model, device)
                        #print(energy)

                        #print(i, filenames[i-1], j, energy)

                        if printProgress: printProgressBar((i-1)*10 + j, len(filenames)*10, prefix = 'Progress:', suffix = 'Complete', length = 50)

                        if energy < 0: 
                            min_filenames.append("min-" + complex_model)
                            energies.append(energy)
                            break
                        else:
                            if j == (numTries-1):
                                call(["rm " + complex_model + " min-" + complex_model], shell=True)


                if printProgress: printProgressBar(len(filenames)*10, len(filenames)*10, prefix = 'Progress:', suffix = 'Complete', length = 50)

                call(["mkdir openmm-minimized"], shell=True)
                call(["rm complex-*.pdb"], shell=True)
                call(["mv min-*.pdb openmm-minimized/"], shell=True)
                call(["cp openmm-minimized/" + min_filenames[np.argmin(energies)] + " ../openmm_min_energy_system.pdb"], shell=True)

                np.savez_compressed("openmm-minimized/openmm_energies.npz", min_filenames=min_filenames, energies=energies)

                os.chdir("..")
    
        if cleanRCD: call(["rm -r RCD"], shell=True)

        os.chdir("..") # get out of rounds index

        print("Writing output of round to disk ...")
        # make 0 folder in true location
        call(["mkdir -p " + current_working_dir + "/" + str(current_round)], shell=True)
        # copy all non-dir files 
        call(["cp " + str(current_round) + "/*.pdb " + current_working_dir + "/" + str(current_round) + "/"], shell=True)
        call(["cp " + str(current_round) + "/*.dcd " + current_working_dir + "/" + str(current_round) + "/"], shell=True)
        call(["cp " + str(current_round) + "/*.npz " + current_working_dir + "/" + str(current_round) + "/"], shell=True)
        call(["cp " + str(current_round) + "/*.log " + current_working_dir + "/" + str(current_round) + "/"], shell=True)
        # copy full_system_conf
        call(["cp -r " + str(current_round) + "/full_system_confs " + current_working_dir + "/" + str(current_round) + "/"], shell=True)


    call(["rm -r " + temp_location], shell=True)
    return #sys.exit(0)

def rescore_with_smina(models, receptor, output_loc, doReceptorMinimization, flexible_residues, useSMINA):

    if not useSMINA and doReceptorMinimization:
        call([smina_location + " -q --scoring vinardo --out_flex " + output_loc + "/receptor_new.pdb --ligand " + models + " --receptor " + receptor + " --autobox_ligand " + models + " --autobox_add 4 --local_only --minimize --flexres " + flexible_residues + " --energy_range 100 --out " + output_loc + "/models_minimize.pdb  > smina.log 2>&1"], shell=True)
    elif not useSMINA and not doReceptorMinimization:
        call([smina_location + " -q --scoring vinardo --ligand " + models + " --receptor " + receptor + " --autobox_ligand " + models + " --autobox_add 4 --local_only --minimize --energy_range 100 --out " + output_loc + "/models_minimize.pdb > smina.log 2>&1"], shell=True)  
    elif useSMINA and doReceptorMinimization:
        call([smina_location + " -q --out_flex " + output_loc + "/receptor_new.pdb --ligand " + models + " --receptor " + receptor + " --autobox_ligand " + models + " --autobox_add 4 --local_only --minimize --flexres " + flexible_residues + " --energy_range 100 --out " + output_loc + "/models_minimize.pdb  > smina.log 2>&1"], shell=True)
    elif useSMINA and not doReceptorMinimization:
        call([smina_location + " -q --ligand " + models + " --receptor " + receptor + " --autobox_ligand " + models + " --autobox_add 4 --local_only --minimize --energy_range 100 --out " + output_loc + "/models_minimize.pdb > smina.log 2>&1"], shell=True)       

def process_smina(ref, data_name, confs_name, native, min_model_index, debug):

    ref_top = ref.top

    f = open("models_minimize.pdb", "r")
    energies = []
    model_index = 0
    model_indices = []
    for line in f:
        first_word = line.split()[0]
        if first_word == 'MODEL':
            f_temp = open("temp.pdb", 'w')
        elif first_word == 'ATOM':
            f_temp.write(line)
        elif first_word == 'ENDMDL':
            model_index += 1
            if debug: print("-------------------", model_index)
            f_temp.close()
            call(["sort -k6 -n temp.pdb > temp_sorted.pdb"], shell=True)
            try:
                conf = get_conf("temp_sorted.pdb", ref_top, "element != H", debug)
                ref += conf
                energies.append(energy)
                model_indices.append(model_index)
            except:
                continue
        elif first_word == 'REMARK':
            if line.split()[1] == 'minimizedAffinity': energy = float(line.split()[2])
    f.close()   

    confs = ref[1:]
    ref = ref[0]

    print(confs)

    np.savez_compressed(data_name, energies=energies, model_indices=model_indices, min_model_index=min_model_index) #, rmsd_to_native=rmsd_to_native)
    confs.save_dcd(confs_name)

  
def get_conf(conf_loc, ref_top, selection, debug):

    ref_atoms = [a for a in ref_top.atoms]

    # load conf with no H
    x_top = md.load(conf_loc).top
    x = md.load(conf_loc, atom_indices=x_top.select(selection))

    if x.top.n_atoms != ref_top.n_atoms and debug: 
        print(x.top.n_atoms)
        print(ref_top.n_atoms)
        raise ValueError('x.top.n_atoms != ref_top.n_atoms')

    x_atoms = [a for a in x.top.atoms]

    # rearrange atom lines to link to reference
    new_xyz = np.zeros((1, x.xyz.shape[1], 3))
    for a, ref_atom in enumerate(ref_atoms):
        for x_index, x_atom in enumerate(x_atoms):
            if str(ref_atom) == str(x_atom):
                new_xyz[0, a, :] = x.xyz[0, x_index, :]
                break

    new_x = md.Trajectory(new_xyz, ref_top)
    
    return new_x  


class ReceptorThread(Thread):

    def __init__(self, loop_indices, filtered_indices, model_indices, num_loops):
        self.loop_indices = loop_indices
        self.filtered_indices = filtered_indices
        self.model_indices = model_indices
        self.num_loops = num_loops
        
        Thread.__init__(self)

    def run(self):
        for j in self.loop_indices:
            if j not in self.filtered_indices: continue 
            print(j)                       

            peptide_j = "temp" + str(self.model_indices[j]).zfill( len(str(self.num_loops)) ) + ".pdb"
            receptor_j = "receptor" + str(self.model_indices[j]).zfill( len(str(self.num_loops)) ) + ".pdbqt"
            call(["cp RCD/input/models/" + receptor_j + " ."], shell=True)
            
            call(["python " + defaults_location + "/rename_atoms.py " + receptor_j], shell=True) # get receptor_j.temp and receptor_j.complete

            complex_j = "target" + str(self.model_indices[j]).zfill( len(str(self.num_loops)) ) + ".pdb"
            call(["cat " + receptor_j + ".complete " + peptide_j + " | sed \"/MODEL/d\" | sed \"/ENDMDL/d\" | sed \"/END/d\" > " + complex_j], shell=True)
            call(["cp " + complex_j + " full_system_confs/" + str(j) + ".pdb"], shell=True)
            call(["rm " + peptide_j + " " + receptor_j + " " + complex_j + " " + receptor_j + ".temp " + receptor_j + ".complete"], shell=True)

            

class RefineThread(Thread):

    def __init__(self, loop_indices, pep_len, num_loops, doReceptorMinimization, flexible_residues, useSMINA, debug):
        self.loop_indices = loop_indices
        self.pep_len = pep_len
        self.num_loops = num_loops
        self.doReceptorMinimization = doReceptorMinimization
        self.flexible_residues = flexible_residues
        self.useSMINA = useSMINA
        self.debug = debug

        if pep_len < 10: self.last_anchor = "\"C   \"" + str(pep_len)
        else:            self.last_anchor = "\"C  \""  + str(pep_len)

        Thread.__init__(self)

    def run(self):
        for i in self.loop_indices:
            if self.debug: print(i)

            model_name_i = "model-" + str(i).zfill( len(str(self.num_loops)) ) + ".pdb"
            partial_name_i = "partial" + str(i).zfill( len(str(self.num_loops)) ) + ".pdbqt"
            fulltemp_model = "fulltemp-" + str(i).zfill( len(str(self.num_loops)) ) + ".pdb"
            complextemp_model = "complextemp-" + str(i).zfill( len(str(self.num_loops)) ) + ".pdb"
            complextemp_model_withH = "complextempwH-" + str(i).zfill( len(str(self.num_loops)) ) + ".pdb"
            full_model = "full-" + str(i).zfill( len(str(self.num_loops)) ) + ".pdb"

            # fill in anchors
            call(["touch " + model_name_i], shell=True)
            call(["grep \"C   1\" ../../../peptide.pdb > " + model_name_i], shell=True)
            call(["cat " + partial_name_i + " >> " + model_name_i], shell=True)
            call(["grep " + self.last_anchor + " ../../../peptide.pdb >> " + model_name_i], shell=True)

            # fill in sidechains
            call([pdbfixer_location + " " + model_name_i + " --output=" + full_model + " --add-atoms=heavy"], shell=True)

            #call(["rm " + model_name_i + " " + partial_name_i], shell=True)
            call(["rm " + partial_name_i], shell=True)

            """ # calling pdbfixer with peptide backbone in receptor, this ends up producing less valid conformations
            call(["cat ../../../receptor.pdb " + model_name_i + " > " + fulltemp_model], shell=True)
            call([pdbfixer_location + " " + fulltemp_model + " --output=" + complextemp_model + " --add-atoms=heavy"], shell=True)
            # create complexes with hydrogen
            #call([pdbfixer_location + " " + fulltemp_model + " --output=" + complextemp_model_withH], shell=True)
            call(["grep \"[A-Z] C \" " + complextemp_model + " > " + full_model], shell=True)

            call(["rm " + model_name_i + " " + partial_name_i + " " + fulltemp_model + " " + complextemp_model], shell=True)
            """

        folder_name = str(self.loop_indices[0])
        call(["mkdir " + folder_name], shell=True)
        call(["touch " + folder_name + "/all_models.pdb"], shell=True)
        for i in self.loop_indices:
            full_model = "full-" + str(i).zfill( len(str(self.num_loops)) ) + ".pdb"
            call(["echo \"MODEL        " + str(i) + "\" >> " + folder_name + "/all_models.pdb"], shell=True)
            call(["sed '/REMARK/d' " + full_model + " | sed '/TER/d' | sed '/END/d' >> " + folder_name + "/all_models.pdb"], shell=True)
            call(["echo \"ENDMDL\" >> " + folder_name + "/all_models.pdb"], shell=True)

        rescore_with_smina(folder_name + "/all_models.pdb", "../../../receptor.pdb", folder_name, self.doReceptorMinimization, self.flexible_residues, self.useSMINA)

# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '*'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()

def minimizeConf(filename, new_filename, device='CPU'):
    
    #print("Opening:", filename)

    pdb = PDBFile(filename)
    top = pdb.getTopology()
    positions = np.array(pdb.positions) #pdb.getPositions(asNumpy=True)
    numAtoms = len(positions)

    #print("Number of atoms:", numAtoms)
    #print("Number of residues:", top.getNumResidues())

    positions = np.reshape(positions, (3*numAtoms,1))

    # run file through pdb fixer first

    #forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
    forcefield = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')
    #forcefield = app.ForceField('amber03.xml', 'amber03_obc.xml')
    #forcefield = app.ForceField('amber10.xml', 'amber10_obc.xml')
    modeller = Modeller(pdb.topology, pdb.positions)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic, constraints=None)

    force_constant = 5000
    force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
    force.addGlobalParameter("k", force_constant)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    protein_particles = md.load(filename).top.select("backbone")

    particle_indices = []
    for protein_particle in protein_particles:
        particle_indices.append(force.addParticle(int(protein_particle), modeller.positions[protein_particle]) )
    system.addForce(force)


    forces = system.getForces()
    i = 0
    for f in forces:
        f.setForceGroup(i)
        i = i + 1


    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    #integrator = VerletIntegrator(0.002*picoseconds)
    platform = Platform.getPlatformByName(device)
    simulation = Simulation(modeller.topology, system, integrator, platform)


    simulation.context.setPositions(modeller.positions)

    #printForces(simulation)

    #print("Minimizing energy ... ")
    simulation.minimizeEnergy()

    #printForces(simulation)

    simulation.reporters.append(app.StateDataReporter(stdout, 100, step=True, 
    potentialEnergy=True, temperature=True, progress=False, remainingTime=True, 
    speed=True, totalSteps=250000, separator='\t'))

    #print "Equilibrating ..."
    #simulation.step(250000)

    #if len(sys.argv) == 2: r = PDBReporter('output.pdb', 1)
    #else: 
    r = PDBReporter(new_filename, 1)
    r.report(simulation, simulation.context.getState(getPositions=True, getEnergy=True))


    return simulation.context.getState(getEnergy=True).getPotentialEnergy() / kilojoule_per_mole

if __name__ == "__main__":
    main(sys.argv[1:])


