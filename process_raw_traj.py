import mdtraj
from mdtools.utils import *
import getopt
import numpy as np

try:
    opts, args = getopt.getopt(sys.argv[1:], "i:o:")
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

input_name = None
output_name = None
sg = 19

for o, a in opts:
	if o == "-i":
		input_name = a
		output_name = input_name
	elif o == "-o":
		output_name = a

def get_file_path(file):
    root = '/n/hekstra_lab/people/ziyuan'
    for (root, dirs, files) in os.walk(root):
        if file in files:
            return os.path.join(root, file)

traj = mdtraj.load(f'{input_name}.h5')
n_frames = len(traj)
print(n_frames)
asu_ref = mdtraj.load(get_file_path('asu_ref.h5'))
unitcell_ref = mdtraj.load(get_file_path('unitcell_ref.h5'))
atom_selection = np.load(get_file_path('atoms_for_alignment.npy'))
print('Loaded all files', flush=True) 

# remove water + ions
traj = traj.remove_solvent()
print('Removed solvent from traj', flush=True)

# wrap at periodic boundary
unwrap_time_axis(traj)
print('Unwrapped the proteins at the boundaries', flush=True)

# produce the first set of subtrajs by chain, by doing unitcell alignment
align_and_split_by_chain(traj, output_name+'_unitcell',
                         unitcell_ref=unitcell_ref, asu_ref=asu_ref, 
                         sg=19, chainwise_alignment=False, asu_reversion=False,                     
                         atom_selection=atom_selection)
print('First set of subtrajs produced', flush=True)

# the second set of subtrajs, chainwise alignment

for i in range(int(np.ceil(n_frames//10000))):
    align_and_split_by_chain(traj[i*10000:(i+1)*10000], output_name+'_chainwise',
                         	unitcell_ref=unitcell_ref, asu_ref=asu_ref, 
                                sg=19, chainwise_alignment=True, 
                         	atom_selection=atom_selection)
print('Done!', flush=True)
