from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import sys
import numpy as np
import time
import mdtraj

# USAGE: python minimize.py <pdb file> <output name (optional)>
# minimize while fixing backbone atoms

def main():

    filename = sys.argv[1]

    print("Opening:", filename)

    pdb = PDBFile(filename)
    top = pdb.getTopology()
    positions = np.array(pdb.positions) #pdb.getPositions(asNumpy=True)
    numAtoms = len(positions)

    print("Number of atoms:", numAtoms)
    print("Number of residues:", top.getNumResidues())

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
    protein_particles = mdtraj.load(filename).top.select("backbone")

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
    platform = Platform.getPlatformByName('OpenCL')
    simulation = Simulation(modeller.topology, system, integrator, platform)


    simulation.context.setPositions(modeller.positions)

    printForces(simulation)

    print("Minimizing energy ... ")
    simulation.minimizeEnergy()

    printForces(simulation)

    simulation.reporters.append(app.StateDataReporter(stdout, 100, step=True, 
    potentialEnergy=True, temperature=True, progress=False, remainingTime=True, 
    speed=True, totalSteps=250000, separator='\t'))

    #print "Equilibrating ..."
    #simulation.step(250000)

    if len(sys.argv) == 2: r = PDBReporter('output.pdb', 1)
    else: r = PDBReporter(sys.argv[2], 1)
    r.report(simulation, simulation.context.getState(getPositions=True, getEnergy=True))

    printForces(simulation)

def printForces(simulation):

    for i in range(simulation.system.getNumForces()):
        f_name = simulation.system.getForce(i).__class__.__name__
        s0 = simulation.context.getState(getEnergy=True, groups=2**i)
        print(f_name + ": " + str(s0.getPotentialEnergy()))

    print("total:", simulation.context.getState(getEnergy=True).getPotentialEnergy())


if __name__ == "__main__":
    main()



