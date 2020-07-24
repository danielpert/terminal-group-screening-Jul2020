'''
Update the oplsaa.xml file with new atomtypes and bonded and nonbonded
potential energy parameters (from the OPLS all atom forcefield)
from an xml file generated by LigParGen.

http://traken.chem.yale.edu/ligpargen/index.html

Note that the SMARTS definitons and atom descriptions must be added
manually. This can be done by using the 'residues' section of the
input xml file, which is deleted in the process of adding to oplsaa.xml,
to create a drawing of the molecule with the atoms labeled by their atom
type number, in order to see which id is which atom, and then adding the
SMARTS definitions and descriptions.

Usage: python update_oplsaa.py -a <path to xml file to add> -t <path to xml file to add to>
       python update_oplsaa.py --add <path to xml file to add> --to <path to xml file to add to>

'''
# Still needs to be done:
#     automatically add overrides
#     automatically add correct indentation

# import relevent packages
import sys
import argparse
import numpy as np
import xml.etree.ElementTree as ET
from itertools import groupby
import signac
import pathlib


def fourier_coefs_to_RB_coefs(f1, f2, f3, f4):
    '''
    Convert fourier coefficients in periodic proper
    dihedral potential equation to Ryckaert-Bellemans
    parameters C0-C5
    
    Paramters
    ----------
    f1, f2, f3, f4 : float
        fourier coefficients 1 thru 4
    
    Returns
    ----------
    c0, c1, c2, c3, c4, c5 : float
        Ryckaert-Belleman parameters C1-C5
    '''
    c0 = f2 + (f1 + f3)/2.0
    c1 = (3*f3 - f1)/2.0
    c2 = 4*f4 - f2
    c3 = -2*f3
    c4 = -4*f4
    c5 = 0.0
    return (c0, c1, c2, c3, c4, c5)


def split_classname(classname):
    '''
    Splits a classname, i.e. 'H805'
    or 'Si1010' into its element
    and its type id
    
    Parameters
    ----------
    classname : str
        classname of atom type
        i.e. 'H805'
    Returns
    ----------
    element : str
        element symbol
        i.e. 'H'
    typeid : int
        number after element symbol
        i.e. 805
    '''
    element = ''
    for isalpha, grouper in groupby(classname, str.isalpha):
        if isalpha and not element:
            element = ''.join(grouper)
        elif not isalpha:
            typeid = int(''.join(grouper))
    return element, typeid


def update_classnums(root, parent_name, child_name, mapping):
    '''
    Given the root of an ElementTree, the name of a
    parent Element, and the name of the child Elements
    to modify:
    
    1) apply a mapping to remove child Elements containing
       the atoms whose IDs map to -1 in the 'mapping' dictionary.
    2) reset the 'class' attributes representing the atoms to
       reflect the new indexing based on the 'mapping' dictionary.
       
    Returns nothing, but rather directly modifies the xml document
    passed in.
    
    Parameters
    ----------
    root : xml.etree.ElementTree.Element
        root element of the ElementTree to be modified
    parent_name : str
        parent element of the elements
        to be altered, i.e. "HarmonicBondForce"
    child_name : str
        child elements to be updated, i.e. "Bond"
    mapping : dict
        mapping from old index to index after
        extra hydrogens are removed, or maps
        to -1 if old index is to be removed
    '''
    parent = root.find(parent_name)
    children = parent.findall(child_name)
    children_to_remove = []
    
    for child in children:
        attr2change = [attr for attr in child.attrib.keys() \
            if (('class' in attr) or ('type' in attr) or ('name' in attr))]
        
        for attr in attr2change:
            if (('type' in attr) or ('name' in attr)):
                index = int(child.get(attr).split('_')[-1])
            else:
                element, index = split_classname(child.get(attr))
            new_index = mapping[index]
            if new_index == -1:
                children_to_remove.append(child)
                break
            if (('type' in attr) or ('name' in attr)):
                child.set(attr, 'opls_{}'.format(new_index))
            else:
                child.set(attr, element + str(new_index))
    for ch in children_to_remove:
        parent.remove(ch)
    return


def add_atom(atom, atoms):
    '''
    Add info from Element from ElementTree
    with tag 'Atom' to a dict of atoms for
    which the keys are the atom IDs and the
    values are the elements
    '''
    element = split_classname(atom.get('name'))[0]
    ID = int(atom.get('type').split('_')[-1])-800
    atoms[ID] = element
    
    
def find_neighbors(atom_id, bonds):
    '''
    Given list of bonds, find ids of atoms
    that atom with `atom_id` is bonded to
    and return as list. List `bonds` should
    contain tuples
    of the form (atom_id_1, atom_id_2)
    '''
    neighbors = []
    for bond in bonds:
        if bond[0] == atom_id:
            neighbors.append(bond[1])
        elif bond[1] == atom_id:
            neighbors.append(bond[0])
    return neighbors


def find_atoms_and_bonds(residues):
    '''
    Given a 'Residues' Element in an ElementTree
    parsed from an oplsaa forcefield in xml
    format, create a dictionary with the atom info and
    a list with the bond info.
    
    Parameters
    ----------
    residues : Element
        contains data on the atoms
        and their bonds for the residue(s)
        
    Returns
    ----------
    atoms : dict
        mapping from atom index to element
    bonds : list
        contains tuples representing each
        bond of the form (atom_id_1, atom_id_2)
    '''
    atoms = dict()
    bonds = []
    for residue in residues:
        for elem in residue:
            if elem.tag == 'Atom':
                add_atom(elem, atoms)
            elif elem.tag == 'Bond':
                bonds.append((int(elem.get('to')), int(elem.get('from'))))
    return atoms, bonds


def find_H_to_remove(atoms, bonds):
    '''
    if more than one hydrogen is bonded to the
    same atom, all of the hydrogens except one
    should be deleted because they are equivilent
    and have the same atom type. This function
    finds all of the hydrogens that should be 
    deleted and returns their indices as keys
    in a dictionary, where the values are the
    indices of the equivilent hydrogens that
    are being kept.

    Parameters
    ----------
    atoms : dict
        mapping from atom index to element
    bonds : list
        contains tuples representing each
        bond of the form (atom_id_1, atom_id_2)

    Returns
    ----------
    to_remove : dict
        keys: index of hydrogen to remove
        value: index of equivilent hydrogen
            being kept 

    '''
    to_remove = dict()
    for atom_id, elem in atoms.items():
        if elem != 'H':
            h_neighbors = []
            for neighbor in find_neighbors(atom_id, bonds=bonds):
                if atoms[neighbor] == 'H':
                    h_neighbors.append(neighbor)
            if len(h_neighbors) > 1:
                kept = h_neighbors[0]
                for h_remove in h_neighbors[1:]:
                    to_remove[h_remove] = kept
    return to_remove


def map_atomtypes(atoms, to_remove, next_idx, starts_at=800):
    '''
    Create mapping from indices of atoms in forcefield
    file being added to new indices that account for the
    removal of duplicate hydrogens and that reflect the
    numbering in the file being added to

    Mapping of -1 means that atom should be removed
    
    Parameters
    ----------
    atoms : dict
        mapping from atom index to element
    to_remove : dict
        keys: index of hydrogen to remove
        value: index of equivilent hydrogen
            being kept 
    next_idx : int
        the next index to be added to oplsaa.xml
        (or whatever file is being appended to)
    starts_at : int, default=800
        The index that the class numbering starts at
        in the xml doc passed in, which is 800 from
        xml documents from the LigParGen server
        
    Returns
    ----------
    mapping : dict
        keys are old indexes, values are indexes that
        should replace the old index, or -1 if the atom
        with the index should be removed
    mapping4angles_dihedrals : dict
        similar to `mapping`, but hydrogens to be deleted
        map to the equivilent hydrogen that is being kept
        instead of mapping to -1. This is because some
        angles and dihedrals involve multiple equivilent
        hydrogens, and deleting all of the angles and dihedrals
        that contain a duplicate hydrogen will result in
        some angle and dihedral parameters being left out of
        the force field. This mapping is used for angles,
        propers, and impropers.
        
    '''
    def find_subtract_by(atom):
        '''
        Find the `subtract_by` variable
        for an atom given the atom id
        '''
        if mapping.get(atom + starts_at) is not None:
            return -mapping.get(atom + starts_at) + atom + next_idx
        else:
            subtract_by = 0
            for h_to_remove, equiv_h in to_remove.items():
                if h_to_remove < atom:
                    subtract_by += 1
            return subtract_by
    mapping = dict()
    mapping4angles_dihedrals = dict()
    for atom_id in list(atoms.keys()):
        subtract_by = 0
        for h_to_remove, equiv_h in to_remove.items():
            if h_to_remove < atom_id:
                subtract_by += 1
            elif h_to_remove == atom_id:
                mapping[atom_id + starts_at] = -1
                mapping4angles_dihedrals[atom_id + starts_at] = equiv_h - find_subtract_by(equiv_h) + next_idx
                break
        else:
            map_from = atom_id + starts_at
            map_to = atom_id - subtract_by + next_idx
            mapping[map_from] = map_to
            mapping4angles_dihedrals[map_from] = map_to
    return mapping, mapping4angles_dihedrals


def main():
    # use the oplsaa.xml by default as the file to add to
    proj = signac.get_project()
    forcefield_filepath = pathlib.Path(
        proj.root_directory() + "./src/util/forcefield/oplsaa.xml"
    )
    
    # parse arguments to get files and ensure correct usage
    parser = argparse.ArgumentParser(description='Process input files.')
    parser.add_argument('-a', '--add', type=str, required=True,
                        dest='file_to_add',
                        help='path to xml file to add')
    parser.add_argument('-t', '--to', type=str,
                        default=forcefield_filepath,
                        dest='file_to_add_to',
                        help='path to xml file to add to')

    args = parser.parse_args()
    file_to_add = vars(args)['file_to_add']
    file_to_add_to = vars(args)['file_to_add_to']
    
    # read and parse xml files
    with open(file_to_add, 'r') as f:
        new_oplsaa = ET.parse(f)
    with open(file_to_add_to, 'r') as f:
        oplsaa = ET.parse(f)
    
    root = oplsaa.getroot()
    new_root = new_oplsaa.getroot()
    
    # create dict of atoms and list of bonds in residue(s)
    residues = new_oplsaa.find("Residues")
    atoms, bonds = find_atoms_and_bonds(residues)
    new_root.remove(residues)
    
    # define next_idx as the number after the last index of the oplsaa file
    # being appended to
    existing_atomtypes = root.findall(".//Type")
    next_idx = int(existing_atomtypes[-1].get('name').split('_')[1]) + 1

    # find duplicate hydrogens to remove and create mapping
    to_remove = find_H_to_remove(atoms, bonds)
    mapping, mapping4angles_dihedrals = map_atomtypes(atoms, to_remove, next_idx)
    
    # remove duplicate hydrogens and update indexing
    forcefield_params = {'Type': 'AtomTypes',
                         'Bond': 'HarmonicBondForce',
                         'Angle': 'HarmonicAngleForce',
                         'Proper': 'PeriodicTorsionForce',
                         'Improper': 'PeriodicTorsionForce',
                         'Atom': 'NonbondedForce'}
    for child_name, parent_name in forcefield_params.items():
        if child_name in ['Type', 'Bond', 'Atom']:
            update_classnums(root=new_root, parent_name=parent_name,
                             child_name=child_name, mapping=mapping)
        else:
            update_classnums(root=new_root, parent_name=parent_name,
                             child_name=child_name, mapping=mapping4angles_dihedrals)
        
    # replace fourier parameters with RB parameters in torsions
    for proper in new_root.find("PeriodicTorsionForce").findall("Proper"):
        fi = list(proper.get('k'+str(i)) for i in range(1, 5))
        fi_converted = (float(f)*4.184/2.092 for f in fi)
        ci = fourier_coefs_to_RB_coefs(*fi_converted)
        for idx in range(6):
            proper.set('c'+str(idx), str(ci[idx]))
        for i in range(1, 5):
            for attr in ['periodicity', 'phase', 'k']:
                proper.attrib.pop(attr+str(i))
    
    # sort new atom types by the number in 'name' attribute
    # and append to the xml file to add to
    for elem_name, attr_name in [('Type', 'name'), ('Atom', 'type')]:
        parent = oplsaa.find(forcefield_params[elem_name])
        data = []
        for elem in new_oplsaa.findall('.//' + elem_name):
            attr = elem.get(attr_name)
            typenum = int(attr.split('_')[-1])
            data.append((typenum, elem))
        data.sort()
        parent.extend(ele[-1] for ele in data)
    
    # append bonded parameters to xml file
    for force in ['Bond', 'Angle', 'Proper', 'Improper']:
        parent_name = forcefield_params[force]
        if force == 'Proper':
            new_parent_name = 'RBTorsionForce'
        else:
            new_parent_name = parent_name
        parent = oplsaa.find(new_parent_name)
        parent.extend(elem for elem in new_oplsaa.findall('./{}/{}'.format(parent_name, force)))
        
    # write updated oplsaa forcefield to xml
    with open(file_to_add_to, 'wb') as f:
        oplsaa.write(f)

    
if __name__ == '__main__':
    main()
