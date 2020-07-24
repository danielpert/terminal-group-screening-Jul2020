import mdtraj as md
import numpy as np 
from scipy.integrate import simps
from scipy.stats import binned_statistic_2d

from util.helper.fileio import read_ndx

def calc_nematic_order(traj_filename, top_filename, output_filename,
                       ndx_filename, n_chains):
    """Calculate the nematic order of both monolayers in a two monolayer system
    Returns the nematic order of each monolayer at each frame of a trajectory
    in a dual monolayer system.
    Parameters
    ----------
    traj_filename : str
        Name of trajectory file
    top_filename : str
        Name of topology file
    output_filename : str
        Name of output file
    ndx_filename : str
        Name of Gromacs .ndx file which specifies atom groups
    n_chains : int
        Number of chains per monolayer
    """
    groups = read_ndx(ndx_filename)
    bottom_chains = [id-1 for id in groups['bottom_chains']]
    bottom_chains = [chain.tolist() for chain in
        np.array_split(bottom_chains, n_chains)]
    top_chains = [id-1 for id in groups['top_chains']]
    top_chains = [chain.tolist() for chain in np.array_split(top_chains, n_chains)]

    traj = md.load(traj_filename, top=top_filename)
    S2_bottom = md.compute_nematic_order(traj, indices=bottom_chains)
    S2_top = md.compute_nematic_order(traj, indices=top_chains)
    np.savetxt(output_filename, np.column_stack((traj.time, S2_bottom, S2_top)),
        header='Time\tBottom\tTop')
