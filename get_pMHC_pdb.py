from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile
from subprocess import call
import sys

# usage: python get_pMHC_pdb.py <pdb code>
# assumes pdb code is of a peptide-MHC structure
# adds missing atoms/residues, removes all waters and ions, labels chains as A,B,C where chain C is the peptide
def main(args):
    
    pdbcode = args[0]

    if len(pdbcode) != 4:
        print("Please enter a correct 4 letter pdbid")
        sys.exit(0)

    call("pdbfixer --pdbid " + pdbcode + " --output=temp.pdb --add-atoms=heavy --keep-heterogens=none", shell=True)
    fixer = PDBFixer(filename="temp.pdb")
    found_A = False
    found_B = False
    found_C = False

    num_chains = len(list(fixer.topology.chains()))
    chains = fixer.topology.chains()
    chains_to_remove = []
    for i, c in enumerate(chains):
        num_residues = len(list(c.residues()))
        if num_residues > 250:
            if not found_A:
                found_A = True
                c.id = "A"
            else:
                c.id = "Z"
                chains_to_remove.append(i)
        elif num_residues > 50 and num_residues < 150:
            if not found_B:
                found_B = True
                c.id = "B"
            else:
                c.id = "Z"
                chains_to_remove.append(i)
        elif num_residues <= 15:
            if not found_C:
                found_C = True
                c.id = "C"
            else:
                c.id = "Z"
                chains_to_remove.append(i)
        else:
            c.id = "Z"
            chains_to_remove.append(i)
            #print "ERROR: Found chains with weird number of residues:", num_residues
            #sys.exit(0)

    fixer.removeChains(chains_to_remove)

    chains = fixer.topology.chains()
    chain_lengths = []
    for c in chains:
        num_residues = len(list(c.residues()))
        chain_lengths.append(num_residues)

    PDBFile.writeFile(fixer.topology, fixer.positions, open(pdbcode + ".pdb", 'w'))
    call(["rm temp.pdb"], shell=True)

    if chain_lengths[1] < chain_lengths[2]:
        call(["grep \"[A-Z] B  \" " + pdbcode + ".pdb > temp.pdb"], shell=True)
        call(["sed -i \"/[A-Z] B  /d\" " + pdbcode + ".pdb"], shell=True)
        call(["sed -i \"/END/d\" " + pdbcode + ".pdb"], shell=True)
        call(["sed -i \"/CONECT/d\" " + pdbcode + ".pdb"], shell=True)
        call(["less temp.pdb >> " + pdbcode + ".pdb"], shell=True)
        fixer = PDBFixer(filename=pdbcode + ".pdb") 
        PDBFile.writeFile(fixer.topology, fixer.positions, open(pdbcode + ".pdb", 'w'))
        call(["rm temp.pdb"], shell=True)

if __name__ == "__main__":
    main(sys.argv[1:])

    




