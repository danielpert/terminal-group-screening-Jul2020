import flow
import numpy as np
import scipy
import os
import pathlib

from foyer import Forcefield
from mbuild.formats.lammpsdata import write_lammpsdata
from mbuild.lib.atoms import H
from util.helper.fileio import write_monolayer_ndx, read_ndx
from util.helper.recipes import Alkylsilane
from util.helper.recipes import DualSurface, SilicaInterface, SurfaceMonolayer
from util.helper.index_groups import generate_index_groups

def system_builder(seed, chainlength=17, backbone=Alkylsilane, terminal_group_1='methyl',
                   terminal_group_2='methyl', terminal_group_3='methyl', frac_1=0.5,
                   frac_2=0.5, num_chains=100):

    """ Define system variable"""
    chainlength = chainlength
    backbone = backbone
    seed = seed
    pattern_type = "random"
    terminal_group_1 = terminal_group_1
    terminal_group_2 = terminal_group_2
    terminal_group_3 = terminal_group_3
    frac_1 = frac_1
    frac_2 = frac_2
    num_chains = num_chains
    """
    -----------------------------------
    Generate amorphous silica interface
    -----------------------------------
    """
    surface_a = SilicaInterface(thickness=1.2, seed=seed)
    surface_b = SilicaInterface(thickness=1.2, seed=seed)

    """
    ------------------------------------------------------
    Generate prototype of functionalized alkylsilane chain
    ------------------------------------------------------
    """
    chain_prototype_A = backbone(
        chain_length=chainlength, terminal_group=terminal_group_1)
    chain_prototype_B = backbone(
        chain_length=chainlength, terminal_group=terminal_group_2)
    chain_prototype_C = backbone(
        chain_length=chainlength, terminal_group=terminal_group_3)
    """
    ----------------------------------------------------------
    Create monolayer on surface, backfilled with hydrogen caps
    ----------------------------------------------------------
    """
    # bottom monolayer is backfilled with the other terminal group
    # num_chains = num_chains * a_fraction
    monolayer_a = SurfaceMonolayer(
        surface=surface_a,
        chains=[chain_prototype_A, chain_prototype_B],
        n_chains=num_chains,
        seed=seed,
        fractions=[frac_1, frac_2],
        backfill=H(),
        rotate=False,
    )
    monolayer_a.name = "Bottom"
    monolayer_b = SurfaceMonolayer(
        surface=surface_b,
        chains=chain_prototype_B,
        n_chains=num_chains,
        seed=seed,
        backfill=H(),
        rotate=False,
    )
    monolayer_b.name = "Top"

    """
    ----------------------
    Create dual monolayers
    ----------------------
    """
    dual_monolayer = DualSurface(
        bottom=monolayer_a, top=monolayer_b, separation=2.0
    )

    """
    --------------------------------------------------------
    Make sure box is elongated in z to be pseudo-2D periodic
    --------------------------------------------------------
    """
    box = dual_monolayer.boundingbox
    dual_monolayer.periodicity += np.array([0, 0, 5.0 * box.lengths[2]])

    """
    -------------------------------------------------------------------
    - Save to .GRO, .TOP, and .LAMMPS formats
    - Atom-type the system using Foyer, with parameters from the OPLS
    force field obtained from GROMACS. Parameters are located in a
    Foyer XML file in  "../util/forcefield/oplsaa.xml".
    -------------------------------------------------------------------
    """

    if os.path.isfile("../util/forcefield/oplsaa.xml"):
        forcefield_filepath = "../util/forcefield/oplsaa.xml"
    elif os.path.isfile("../../../util/forcefield/oplsaa.xml"):
        forcefield_filepath = "../../../util/forcefield/oplsaa.xml"
    else:
        raise Exception('Forcefield file is not found')    


    dual_monolayer.save("init.gro", residues=["Top", "Bottom"], overwrite=True)

    structure = dual_monolayer.to_parmed(
        box=None, residues=["Top", "Bottom"]
        )
    ff = Forcefield(forcefield_files=forcefield_filepath)
    structure = ff.apply(structure)
    structure.combining_rule = "geometric"

    structure.save("init.top", overwrite=True)
    write_lammpsdata(filename="init.lammps", structure=structure)

    """
    --------------------------------------
    Specify index groups and write to file
    --------------------------------------
    """
    index_groups = generate_index_groups(
        system=dual_monolayer,
        terminal_group=terminal_group,
        freeze_thickness=0.5,
    )
    write_monolayer_ndx(rigid_groups=index_groups, filename="init.ndx")
    return dual_monolayer
