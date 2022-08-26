from openmm.app import PDBFile, ForceField
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator
from simtk.unit import *
from tqdm import tqdm
import mdtraj
import mdtools
from mdtools.utils import *
import pickle
import argparse
import subprocess
import os

# parameters
parser = argparse.ArgumentParser()
parser.add_argument("-E", "--E", type=str, help="Field strength of the THz pulse (MV/cm)", default=10e8)
parser.add_argument("-t0", "--t0", type=int, help="Initial duration of NVT equilibration (ns)", default=10)
parser.add_argument("-t1", "--t1", type=int, help="Number of pulse cycles for the production run")
parser.add_argument("-t2", "--n_pulses", type=int, help="Number of pulses per cycle", default=100)
parser.add_argument("-i", "--input", type=str, help="Input file for the crystal system", default="squeezed.pdb")
parser.add_argument("-o", "--output", type=str, help="Prefix for the output trajectory and state files")
parser.add_argument("-n", "--n_chains", type=int, help="Number of protein chains in the system", default=4)
args = parser.parse_args()

def get_file_path(file):
    root = '/n/hekstra_lab/people/ziyuan'
    for (root, dirs, files) in os.walk(root):
        if file in files:
            return os.path.join(root, file)

# initialize forcefield
# need to have openmmforcefields installed first
forcefield = ForceField('amber/ff14SB.xml',
                        'amber/tip3p_standard.xml',
                        'amber/tip3p_HFE_multivalent.xml',
                        get_file_path('cro.xml'))
efx = getFieldStrength(args.E * volts/meters)
print("Field strength is", efx)

crystal = PDBFile(get_file_path(args.input))
mdsystem = mdtools.LatticeMDSystem(crystal.topology,
                                   crystal.positions,
                                   forcefield, "P 21 21 21")
print("Started simulation.", flush=True)
mdsystem.buildSimulation(ensemble="NVT", posre=True,
                         saveTrajectory=False, saveStateData=False,
                         posre_sel="not water and not (element Na or element Ca) and not element H",
                         dt=0.002*picoseconds)
mdsystem.equilibrate(args.t0*nanoseconds, posre=True)
print("Finished NVT equilibration.", flush=True)
mdsystem.buildSimulation(ensemble="NVT",  filePrefix=args.output,
                         saveTrajectory=True, saveStateData=True,
                         trajInterval=50, stateDataInterval=50,
                         dt=0.002*picoseconds, efx=True) # record every 0.1ps

asu_ref = mdtraj.load(get_file_path('asu_ref.h5'))
unitcell_ref = mdtraj.load(get_file_path('unitcell_ref.h5'))
atom_selection = np.load(get_file_path('atoms_for_alignment.npy'))
print("Loaded auxiliary data files.")

for i in tqdm(range(args.t1)):
    for j in tqdm(range(args.n_pulses)):
        # Up
        mdsystem.simulation.context.setParameter('Ex', efx)
        mdsystem.simulate(1*picoseconds) # 1ps

        # Down
        mdsystem.simulation.context.setParameter('Ex', -1*efx)
        mdsystem.simulate(1*picoseconds) # 1ps

        # Off
        mdsystem.simulation.context.setParameter('Ex', 0.0)
        mdsystem.simulate(8*picoseconds) # 8ps

    # post-process this segment
    for reporter in mdsystem.simulation.reporters:
        try:
            reporter.close()
        except:
            continue

    # remove solvent, unwrap, align
    traj = mdtraj.load(f'{args.output}.h5')
    traj.remove_solvent()
    unwrap_time_axis(traj)
    for offset, phase in [(9, 'pos'), (19, 'neg'), (99, 'zero')]:
        fname=args.output+f'_epoch_{i}_chainwise_{phase}'
        align_and_split_by_chain(traj[offset::100], fname,
                                unitcell_ref=unitcell_ref, asu_ref=asu_ref,
                                sg=19, chainwise_alignment=True,
                                atom_selection=atom_selection)
    print("Aligned and split into subtrajs")

    # convert to snapshots, calculate structural factors, and average
    for phase in ['pos', 'neg', 'zero']:
        for k in range(args.n_chains):
            fname = args.output+f'_epoch_{i}_chainwise_{phase}_subtraj_{k}'
            save_snapshots_from_traj(mdtraj.load(f'{fname}.h5'), output_name=fname, frame_offset=0, d_frame=1) # should give 100 frames
            batch_annotate_spacegroup(fname, args.n_pulses, "P 21 21 21")
            batch_fmodel(fname, max_frame=args.n_pulses, resolution=1.5,
			phenix_command='source /n/hekstra_lab/people/ziyuan/egfp/phenix_env.sh; phenix.fmodel')
            average_structure_factors(fname, max_frame=args.n_pulses)
    print("Computed average reflection for pos, neg, and zero parts")

    # cleanup
    subprocess.run(["mv", "*avg*", "./data"])
    subprocess.run(["rm", "*snapshot*"])
    subprocess.run(["rm", f'{args.output}.h5'])
    print("Cleaned up intermediate files")

print("Finished terahertz pulse production run", flush=True)


