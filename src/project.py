import textwrap
import glob
import json
import logging
import os
import pathlib

import signac
import flow
import numpy as np
import scipy
from flow import FlowProject, environments
from foyer import Forcefield
from mbuild.formats.lammpsdata import write_lammpsdata
from mbuild.lib.atoms import H

from util.helper.functions import system_builder
from util.helper.recipes import Alkylsilane, DualSurface, SilicaInterface, SurfaceMonolayer
from util.helper.index_groups import generate_index_groups
from util.helper.fileio import write_monolayer_ndx, save_pattern, read_ndx


"""
test run with own version of lammps

set up for beginning em with gromacs using titan tonight

if time allows: then support nvt equil first
"""


def _setup_logger(logger_name, log_file, level=logging.INFO):
    mylog = logging.getLogger(logger_name)
    formatter = logging.Formatter("%(asctime)s : %(message)s")
    fileHandler = logging.FileHandler(log_file, mode="a")
    fileHandler.setFormatter(formatter)

    mylog.setLevel(level)
    mylog.addHandler(fileHandler)
    return mylog


def _mdrun_str(op_name):
    """return string formatted based on name of operation"""
    msg = (
        "gmx_mpi mdrun -v -deffnm {op} -s {op}.tpr -cpi {op}.cpt "
        "-ntomp 1 -gpu_id {}".format("0" * 16, op=op_name)
    )
    return msg


"""def _gmx_trjcat_shear(shear_force, extension):
    msg = ("gmx_mpi trjcat -f shear_{{shear}}nN*.{{ext}} -o shear_{{shear}}nN_combined.{{ext}}".format(
        shear_force, extension, shear_force, extension)
    )
    return msg
"""


def __trjcat_naive_shear(shear_force, file_list, extension):
    msg = "gmx_mpi trjcat -f {} -o shear_{}nN_combined.{}".format(
        " ".join(file_list), shear_force, extension
    )
    return msg


def _wraprun_add_tasks_helper(n_cpus, working_dir, command):
    """For serial jobs, build up a large string of executables, sep by :"""
    msg = "-n {} --w-cd {} {} ".format(n_cpus, working_dir, command)
    return msg


class SimpleProject(flow.environment.DefaultTorqueEnvironment):
    # Modify so that this work with hulk environment
    template = "hulk.sh"

    @classmethod
    def add_args(cls, parser):
        super(flow.environment.DefaultTorqueEnvironment, cls).add_args(parser)
        parser.add_argument(
            "--group", type=int, help="How many serial jobs at once"
        )
        parser.add_argument(
            "--walltime", type=float, help="walltime"
        )

        
class Project(FlowProject):
    pass


@Project.label
def is_initialized(job):
    result = bool(
        job.isfile("init.gro")
        and job.isfile("init.top")
        and job.isfile("init.lammps")
        and job.isfile("init.ndx")
    )
    return result


@Project.label
def overlaps_fixed(job):
    return job.isfile("minimize.xtc")


@Project.label
def converted_lmp(job):
    return job.isfile("minimize.gro")


@Project.label
def is_em_grompp(job):
    return job.isfile("em.tpr")


@Project.label
def em_completed(job):
    return job.isfile("em.gro")


@Project.label
def is_nvt_grompp(job):
    return job.isfile("nvt.tpr")


@Project.label
def nvt_completed(job):
    return job.isfile("nvt.gro")


@Project.label
def is_compress_grompp(job):
    return job.isfile("compress.tpr")


@Project.label
def compress_completed(job):
    return job.isfile("compress.gro")


@Project.label
def is_shear_5nN_grompp(job):
    return job.isfile("shear_5nN.tpr")


@Project.label
def shear_5nN_completed(job):
    return bool(glob.glob(job.fn("shear_5nN*.gro")))


@Project.label
def is_shear_15nN_grompp(job):
    return job.isfile("shear_15nN.tpr")


@Project.label
def shear_15nN_completed(job):
    return bool(glob.glob(job.fn("shear_15nN*.gro")))


@Project.label
def is_shear_25nN_grompp(job):
    return job.isfile("shear_25nN.tpr")


@Project.label
def shear_25nN_completed(job):
    return bool(glob.glob(job.fn("shear_25nN*.gro")))


@Project.label
def friction_calculated(job):
    is_calculated = list()
    for shear in [5, 15, 25]:
        is_calculated.append(job.isfile("friction_{}nN.txt".format(shear)))
    return all(is_calculated)


@Project.label
def cof_calculated(job):
    "if the job document has a listing for COF and intercept"
    is_calculated = (job.document.get("COF"), job.document.get("intercept"))
    return all(is_calculated)


@Project.operation
@Project.post.isfile("init.top")
@Project.post.isfile("init.gro")
@Project.post.isfile("init.lammps")
@Project.post.isfile("init.ndx")
def initialize_system(job):
    """ Generate the monolayer surfaces, parametrize, save LAMMPS, GRO, TOP.
    """
    
    """
    -------------------------------
    Declare the backbone dictionary
    -------------------------------
    """
    backbone_dict = {'Alkylsilane':Alkylsilane}
    """
    ---------------------------
    Read statepoint information
    ---------------------------
    """
    chainlength = job.statepoint()["chainlength"]
    backbone = backbone_dict[job.statepoint()["backbone"]]
    seed = job.statepoint()["seed"]
    pattern_type = job.statepoint()["pattern_type"]
    terminal_group_1 = job.statepoint()["terminal_group_1"]
    terminal_group_2 = job.statepoint()["terminal_group_2"]
    terminal_group_3 = job.statepoint()["terminal_group_3"]
    frac_1 = job.statepoint()["frac_1"]
    frac_2 = job.statepoint()["frac_2"]
    num_chains = job.statepoint()["n"]

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
    chain_prototype_1 = backbone(
        chain_length=chainlength, terminal_group=terminal_group_1)
    chain_prototype_2 = backbone(
        chain_length=chainlength, terminal_group=terminal_group_2)
    chain_prototype_3 = backbone(
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
        chains=[chain_prototype_1, chain_prototype_2],
        fractions=[frac_1, frac_2],
        n_chains=num_chains,
        seed=seed,
        backfill=H(),
        rotate=False,
    )
    monolayer_a.name = "Bottom"
    monolayer_b = SurfaceMonolayer(
        surface=surface_b,
        chains=chain_prototype_3,
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
    Foyer XML file in the `atools` git repo, with references provided
    as well as notes where parameters have been added or altered to
    reflect the literature.
    -------------------------------------------------------------------
    """
    # path for project root dir
    proj = signac.get_project()
    forcefield_filepath = pathlib.Path(
        proj.root_directory() + "/src/util/forcefield/oplsaa.xml"
    )
    # change into job directory
    _switch_dir(job)
    logging.info("at dir: {}".format(job.ws))
    dual_monolayer.save("init.gro", residues=["Top", "Bottom"], overwrite=True)

    if not (
        job.isfile("init.top")
        and job.isfile("init.lammps")
        and job.isfile("init.gro")
    ):

        structure = dual_monolayer.to_parmed(
            box=None, residues=["Top", "Bottom"]
        )
        ff = Forcefield(forcefield_files=forcefield_filepath.as_posix())
        structure = ff.apply(structure, verbose=True)
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
        terminal_groups=[terminal_group_1,
                         terminal_group_2,
                         terminal_group_3],
        freeze_thickness=0.5,
    )
    write_monolayer_ndx(rigid_groups=index_groups, filename="init.ndx")


def _switch_dir(job):
    p = pathlib.Path(job.workspace())
    os.chdir(str(p.absolute()))

'''
@Project.operation
@flow.directives(nranks=16)
@Project.pre.after(initialize_system)
@Project.post.isfile("system.ndx")
def create_system_ndx(job):
    ndx_file = pathlib.Path(job.workspace() + "/init.ndx")
    ndx_dict = read_ndx(str(ndx_file.resolve()))
    pathlib.Path(job.workspace() + "/system.ndx").touch(exist_ok=True)
    sys_ndx = pathlib.Path(job.workspace() + "/system.ndx")
    with sys_ndx.open(mode="wt") as outfile:
        outfile.write("[ System ]\n")
        atoms = "{}\n".format(" ".join(str(x+1) for x in ndx_dict["System"]))
        outfile.write(textwrap.fill(atoms, 80))
        outfile.write("\n")
#@Project.pre.after(create_system_ndx)
'''

@Project.operation
@Project.pre(is_initialized)
@Project.pre.after(initialize_system)
@Project.post.isfile("minimize.xtc")
@flow.cmd
def fix_overlaps(job):
    cmds = pathlib.Path(
        signac.get_project().root_directory() + "/src/util/mdp_files"
    )
    return "{} lmp_mpi -in {}/in.minimize -log {}/minimize.log".format(
        job.workspace(), str(cmds.absolute()), job.workspace()
    )


@Project.operation
@Project.pre(overlaps_fixed)
@Project.pre.after(fix_overlaps)
@Project.post.isfile("minimize.gro")
@flow.cmd
def lmp_to_gmx(job):
    msg = "cd {}; echo 0 | aprun -n 1 gmx_mpi trjconv -s init.gro -f minimize.xtc -o minimize.gro  -b 1.0 -e 1.0".format(job.workspace())
    return msg


@Project.operation
@Project.pre(converted_lmp)
@Project.pre.after(lmp_to_gmx)
@Project.post.isfile("em.tpr")
@flow.cmd
def em_grompp(job):
    em_mdp_path = pathlib.Path(
        signac.get_project().root_directory() + "/src/util/mdp_files/em.mdp"
    )
    #msg = "{}".format(job.workspace())
    msg = "cd {}; aprun -n 1 gmx_mpi grompp -f {} -c minimize.gro -p init.top \
           -n init.ndx -o em.tpr -maxwarn 1".format(job.workspace(), em_mdp_path)
    return msg


@Project.operation
@flow.directives(ngpu=1)
@Project.pre(is_em_grompp)
@Project.pre.after(em_grompp)
@Project.post.isfile("em.gro")
@flow.cmd
def mdrun_em(job):
    msg = _mdrun_str("em")
    return "{} -- {}".format(job.workspace(), msg)


@Project.operation
@Project.pre.after(mdrun_em)
@Project.post.isfile("nvt.tpr")
@flow.cmd
def nvt_equil_grompp(job):
    nvt_mdp_path = pathlib.Path(
        signac.get_project().root_directory() + "/src/util/mdp_files/nvt.mdp"
    )
    msg = "cd {}; aprun -n 1 gmx_mpi grompp -f {} -c {} -p {} -n {} -o {} -maxwarn 1".format(
        job.workspace(),
        nvt_mdp_path,
        "em.gro",
        "init.top",
        "init.ndx",
        "nvt.tpr",
    )
    return msg


@Project.operation
@flow.directives(ngpu=1)
@Project.pre(is_nvt_grompp)
@Project.pre.after(nvt_equil_grompp)
@Project.post.isfile("nvt.gro")
@flow.cmd
def mdrun_nvt(job):
    msg = _mdrun_str("nvt")
    return "{} -- {}".format(job.workspace(), msg)


@Project.operation
@Project.pre.after(mdrun_nvt)
@Project.post.isfile("compress.tpr")
@flow.cmd
def compress_grompp(job):
    compress_mdp_path = pathlib.Path(
        signac.get_project().root_directory()
        + "/src/util/mdp_files/compress.mdp"
    )
    msg = "cd {}; aprun -n 1 gmx_mpi grompp -f {} -c {} -p {} -n {} -o {} -maxwarn 1".format(
        job.workspace(),
        compress_mdp_path,
        "nvt.gro",
        "init.top",
        "init.ndx",
        "compress.tpr",
    )
    return msg


@Project.operation
@flow.directives(ngpu=1)
@Project.pre.after(compress_grompp)
@Project.post.isfile("compress.gro")
@flow.cmd
def mdrun_compress(job):
    msg = _mdrun_str("compress")
    return "{} -- {}".format(job.workspace(), msg)


@Project.operation
@Project.pre.after(mdrun_compress)
@Project.post.isfile("shear_5nN.tpr")
@flow.cmd
def shear_5nN_grompp(job):
    shear_5nN_mdp_path = pathlib.Path(
        signac.get_project().root_directory()
        + "/src/util/mdp_files/shear_5nN.mdp"
    )
    msg = "cd {}; aprun -n 1 gmx_mpi grompp -f {} -c {} -p {} -n {} -o {} -maxwarn 1".format(
        job.workspace(),
        shear_5nN_mdp_path,
        "compress.gro",
        "init.top",
        "init.ndx",
        "shear_5nN.tpr",
    )
    return msg


@Project.operation
@flow.directives(ngpu=1)
@Project.pre.after(shear_5nN_grompp)
@Project.post(shear_5nN_completed)
@flow.cmd
def mdrun_shear_5nN(job):
    msg = _mdrun_str("shear_5nN")
    return "{} -- gmx_mpi mdrun -s shear_5nN.tpr -deffnm shear_5nN -cpi shear_5nN.cpt -cpo shear_5nN.cpt -noappend -ntomp 1 -gpu_id 0000000000000000".format(
        job.workspace()
    )
    # return "cd {}; aprun -n 16 {} -append".format(job.workspace(), msg)


@Project.operation
@Project.pre.after(mdrun_compress)
@Project.post.isfile("shear_15nN.tpr")
@flow.cmd
def shear_15nN_grompp(job):
    shear_15nN_mdp_path = pathlib.Path(
        signac.get_project().root_directory()
        + "/src/util/mdp_files/shear_15nN.mdp"
    )
    msg = "cd {}; aprun -n 1 gmx_mpi grompp -f {} -c {} -p {} -n {} -o {} -maxwarn 1".format(
        job.workspace(),
        shear_15nN_mdp_path,
        "compress.gro",
        "init.top",
        "init.ndx",
        "shear_15nN.tpr",
    )
    return msg


@Project.operation
@flow.directives(ngpu=1)
@Project.pre(is_nvt_grompp)
@Project.pre.after(shear_15nN_grompp)
@Project.post(shear_15nN_completed)
@flow.cmd
def mdrun_shear_15nN(job):
    msg = _mdrun_str("shear_15nN")
    return "{} -- gmx_mpi mdrun -s shear_15nN.tpr -deffnm shear_15nN -cpi shear_15nN.cpt -cpo shear_15nN.cpt -noappend -ntomp 1 -gpu_id 0000000000000000".format(
        job.workspace()
    )
    # return "cd {}; aprun -n 16 {} -append".format(job.workspace(), msg)


@Project.operation
@Project.pre.after(mdrun_compress)
@Project.post.isfile("shear_25nN.tpr")
@flow.cmd
def shear_25nN_grompp(job):
    shear_25nN_mdp_path = pathlib.Path(
        signac.get_project().root_directory()
        + "/src/util/mdp_files/shear_25nN.mdp"
    )
    msg = "cd {}; aprun -n 1 gmx_mpi grompp -f {} -c {} -p {} -n {} -o {} -maxwarn 1".format(
        job.workspace(),
        shear_25nN_mdp_path,
        "compress.gro",
        "init.top",
        "init.ndx",
        "shear_25nN.tpr",
    )
    return msg


@Project.operation
@flow.directives(ngpu=1)
@Project.pre(is_nvt_grompp)
@Project.pre.after(shear_25nN_grompp)
@Project.post(shear_25nN_completed)
@flow.cmd
def mdrun_shear_25nN(job):
    # msg = _mdrun_str('shear_25nN')
    # return "cd {}; aprun -n 16 {} -append".format(job.workspace(), msg)
    return "{} -- gmx_mpi mdrun -s shear_25nN.tpr -deffnm shear_25nN -cpi shear_25nN.cpt -cpo shear_25nN.cpt -noappend -ntomp 1 -gpu_id 0000000000000000".format(
        job.workspace()
    )


@Project.operation
@Project.pre(shear_5nN_completed)
@Project.post(cof_calculated)
@flow.cmd
def shear_5nN_xtc_concat(job):
    file_list = glob.glob(job.workspace() + "/shear_5nN*.xtc")
    executable = __trjcat_naive_shear(5, file_list, "xtc")
    return _wraprun_add_tasks_helper(1, job.workspace(), executable)


@Project.operation
@Project.pre(shear_5nN_completed)
@Project.pre.after(mdrun_shear_5nN)
@Project.post(cof_calculated)
@flow.cmd
def shear_5nN_trr_concat(job):
    file_list = glob.glob(job.workspace() + "/shear_5nN*.trr")
    executable = __trjcat_naive_shear(5, file_list, "trr")
    return _wraprun_add_tasks_helper(1, job.workspace(), executable)


@Project.operation
@Project.pre(shear_15nN_completed)
@Project.pre.after(mdrun_shear_15nN)
@Project.post(cof_calculated)
@flow.cmd
def shear_15nN_xtc_concat(job):
    file_list = glob.glob(job.workspace() + "/shear_15nN*.xtc")
    executable = __trjcat_naive_shear(15, file_list, "xtc")
    return _wraprun_add_tasks_helper(1, job.workspace(), executable)


@Project.operation
@Project.pre(shear_15nN_completed)
@Project.pre.after(mdrun_shear_15nN)
@Project.post(cof_calculated)
@flow.cmd
def shear_15nN_trr_concat(job):
    file_list = glob.glob(job.workspace() + "/shear_15nN*.trr")
    executable = __trjcat_naive_shear(15, file_list, "trr")
    return _wraprun_add_tasks_helper(1, job.workspace(), executable)


@Project.operation
@Project.pre(shear_25nN_completed)
@Project.pre.after(mdrun_shear_25nN)
@Project.post(cof_calculated)
@flow.cmd
def shear_25nN_xtc_concat(job):
    file_list = glob.glob(job.workspace() + "/shear_25nN*.xtc")
    executable = __trjcat_naive_shear(25, file_list, "xtc")
    return _wraprun_add_tasks_helper(1, job.workspace(), executable)


@Project.operation
@Project.pre(shear_25nN_completed)
@Project.pre.after(mdrun_shear_15nN)
@Project.post(cof_calculated)
@flow.cmd
def shear_25nN_trr_concat(job):
    file_list = glob.glob(job.workspace() + "/shear_25nN*.trr")
    executable = __trjcat_naive_shear(25, file_list, "trr")
    return _wraprun_add_tasks_helper(1, job.workspace(), executable)


@Project.operation
@Project.pre(shear_5nN_completed)
@Project.pre(shear_15nN_completed)
@Project.pre(shear_25nN_completed)
@Project.post(cof_calculated)
@Project.post(friction_calculated)
def mdanalysis_concat_trr(job):
    import MDAnalysis
    print(job.workspace()) 
    for load in [5, 15, 25]:
        if len(glob.glob(job.workspace() + "shear_{}nN_combined.trr".format(load))) is 0:
            trr_list = sorted(glob.glob(job.workspace() + "/shear_{}nN.*.trr".format(load)))
            if len(trr_list) is 0:
                continue
            else:
                tpr_file = sorted(glob.glob(str(job.workspace() + "/shear_{}nN*tpr").format(load)), reverse=True)[0]
                u = MDAnalysis.Universe(tpr_file, trr_list)
                all_system = u.select_atoms("all")
                with MDAnalysis.Writer("{}/shear_{}nN_combined.trr".format(job.workspace(), load), all_system.n_atoms) as outfile:
                    for ts in u.trajectory:
                        outfile.write(all_system)
        else:
            continue
"""        
@Project.operation
@flow.directives(nranks=1)
@Project.pre.after(mdrun_shear_5nN)
@Project.pre.after(mdrun_shear_15nN)
@Project.pre.after(mdrun_shear_25nN)
@Project.post.isfile("shear_5nN_combined.trr")
@Project.post.isfile("shear_15nN_combined.trr")
@Project.post.isfile("shear_25nN_combined.trr")
@flow.cmd
def mdanalysis_concat_trr(job):
   return "python {}".format(str(signac.get_project().root_directory() + "/src/concat_wraprun.py"))
"""


@Project.operation
@Project.pre.after(mdanalysis_concat_trr)
@Project.post(cof_calculated)
@Project.post(friction_calculated)
def calc_friction_system(job):
    from atools.structure_mixed import calc_friction

    for load in [5, 15, 25]:
        trr_file = "{}/shear_{}nN_combined.trr".format(job.workspace(), load)
        out_file = "{}/friction_{}nN.txt".format(job.workspace(), load)
        ndx_file = "{}/init.ndx".format(job.workspace())
        calc_friction(
            trr_filename=trr_file,
            output_filename=out_file,
            ndx_filename=ndx_file,
        )


@Project.operation
@Project.pre(friction_calculated)
@Project.pre.after(calc_friction_system)
@Project.post(cof_calculated)
def calc_cof(job):
    from scipy import stats

    loads = [5, 15, 25]
    friction_forces = []
    for load in loads:
        friction_data = np.loadtxt(
            "{}/friction_{}nN.txt".format(job.workspace(), load)
        )
        friction_data = friction_data[int(len(friction_data) * 0.6) :, 1]
        friction_forces.append(np.mean(friction_data))
    cof, intercept, r, p, stderr = stats.linregress(loads, friction_forces)
    job.document["COF"] = cof
    job.document["intercept"] = intercept


if __name__ == "__main__":
    Project().main()
