import textwrap

import mbuild as mb


def write_monolayer_ndx(rigid_groups, filename):
    with open(filename, 'w') as f:
        for name, indices in rigid_groups.items():
            f.write('[ {name} ]\n'.format(name=name))
            atoms = '{}\n'.format(' '.join(str(x) for x in indices))
            f.write(textwrap.fill(atoms, 80))
            f.write('\n')


def save_pattern(filename, pattern, overwrite=False):
    lj_proto = mb.Compound(name='LJ')
    lj_box = mb.Compound()
    for pos in pattern:
        lj_particle = mb.clone(lj_proto)
        lj_particle.translate(pos)
        lj_box.add(lj_particle)
    lj_box.save(filename, overwrite=overwrite)


def read_ndx(filename):
    """Loads a Gromacs .ndx file into a dictionary. """
    ndx_dict = {}
    vals = []
    with open(filename, 'r') as ndx:
        for line in ndx:
            if '[' in line:
                if vals:
                    ndx_dict.update({group:vals})
                    vals = []
                group = line.strip()[2:-2]
            else:
                for atom_id in line.split():
                    vals.append(int(atom_id) - 1)
    ndx_dict.update({group:vals})
    return ndx_dict
