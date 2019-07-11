from subprocess import call
import time
import sys
import os

defaults_location = os.path.dirname(os.path.abspath(__file__))

# python model_receptor.py <file containing alpha chain sequence> <template pdb>

def main(args):

	f = open(args[0], 'r')
	f_str = ""
	for line in f:
		f_str += line
	f.close()
	target_sequence = "".join(f_str.split())
	target_sequence += "/"
	# constant beta immunoglobin
	target_sequence += "MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM*" 
	receptor_template = args[1]



	print("Preparing receptor template")
	call(["cp " + receptor_template + " receptor_template.pdb"], shell=True)
	call(["sed -i \"/[A-Z] C  /d\" receptor_template.pdb"], shell=True)
	call(["python " + defaults_location + "/mod_getseq.py > log.txt"], shell=True) # creates receptor_template.seq

	print("Preparing target sequence")
	f = open("target_sequence.pir", 'w')
	f.write(">P1;target_sequence\n")
	f.write("sequence::     : :     : :::-1.00:-1.00\n")
	f.write(target_sequence)
	f.close()

	print("Aligning target sequence with receptor template")
	starttime = time.time()
	call(["python " + defaults_location + "/align2d.py >> log.txt"], shell=True) # creates target_sequence-receptor_template.ali
	endtime = time.time()
	print("Aligning took " + str(endtime - starttime) + " seconds.")

	print("Creating model")
	starttime = time.time()
	call(["python -u  " + defaults_location + "/model-single-opt.py 2 >> log.txt"], shell=True)
	call(["cp $(grep \"Top model\" log.txt | awk '{ print $3 }') best_model.pdb"], shell=True)
	endtime = time.time()
	print("Homology modelling took " + str(endtime - starttime) + " seconds.")


if __name__ == "__main__":
    main(sys.argv[1:])