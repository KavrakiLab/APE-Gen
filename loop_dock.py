from subprocess import call
from subprocess import check_output
import sys, os

pdbid = sys.argv[1]

min_energy = 10000
for i in range(10):
	call(["mkdir " + str(i)], shell=True)
	os.chdir(str(i))

	if i == 0: call(["cp ../" + pdbid + ".pdb .; python ~/Dropbox/research/analyze-pMHC/dock/dock.py " + pdbid + ".pdb HLA-A*02:01-from3I6L -a 100 > log" + str(i) + ".txt"], shell=True)
	else: 
		#if os.path.exists("../" + str(i-1) + "/min_energy_system.pdb"):
		call(["cp ../" + str(i-1) + "/min_energy_system.pdb ."], shell=True)
		call(["cp ../" + str(i-1) + "/native.pdb ."], shell=True)
		#call(["cp ../1hhh.pdb .; python ~/Dropbox/research/analyze-pMHC/dock/dock.py 1hhh.pdb test -a 100000 -l 50 > log" + str(i) + ".txt"], shell=True)
		call(["python ~/Dropbox/research/analyze-pMHC/dock/dock.py min_energy_system.pdb REDOCK -a 100 > log" + str(i) + ".txt"], shell=True)
		#else:
		#	call(["cp ../4u6y.pdb .; python ~/Dropbox/research/analyze-pMHC/dock/dock.py 4u6y.pdb HLA-A*02:01-from3I6L -a 10 > log" + str(i) + ".txt"], shell=True)

	min_rmsd_i = float(check_output(["grep 'Min RMSD to native:' log" + str(i) + ".txt | awk '{ print $5 }'"], shell=True))
	print i, min_rmsd_i
	if min_rmsd_i < 2.5: print "Sampled native pose. Exiting."; break

	min_energy_i = float(check_output(["grep 'energy of selected binding mode:' log" + str(i) + ".txt | awk '{ print $6 }'"], shell=True))

	print min_energy_i
	#if min_energy_i < min_energy:
	#	min_energy = min_energy_i
	call(["sed -i \"s/          /           /g\" min_energy_system.pdb"], shell=True) # fix annoying openbabel warning
	#call(["cp min_energy_system.pdb ~/Dropbox/research/analyze-pMHC/dock/templates/best_model.pdb"], shell=True)

	os.chdir("..")

	if os.path.exists("STOP"): print "Detected STOP file."; break


