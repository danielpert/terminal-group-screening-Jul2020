import mbuild as mb
from .alkane import Alkane
from mbuild.lib.moieties import Silane
from .one_port import Acetyl, Acetylene, Amino, Biphenyl, Carboxyl, Cyano, \
    Cyclopropyl, Ethylene, Fluorophenyl, Formyl, Hydroxyl, Isopropyl, \
    Methoxy, Methyl, Nitro, Nitrophenyl, Pentafluorophenyl, \
    Perfluoromethyl, Phenyl, Pyrrole, Triazole, Difluoromethyl, Phenol, \
    Toluene, Benzoicacid, Isopropylbenzene, Amide, Cyclohexyl, Sulfhydryl, \
    Nitroso


class Alkylsilane(mb.Compound):
    """A terminal-functionalized alkylsilane chain.

    An alkylsilane chain featuring a user-specified functional group at one
    terminus and a silane group (featuring an open port to attach to a surface)
    at the other terminus.

    Parameters
    ----------
    chain_length : int
        Length of the chain (number of carbons)
    terminal_group : str
       Functional group to attach to the chain terminus. Valid option for this
       repository is `methyl`, but more can be easily added by providing
       appropriate supplement structure files.
    """
    def __init__(self, chain_length, terminal_group):
        super(Alkylsilane, self).__init__()
        terminal_group_dict = \
            {'acetyl':Acetyl, 'acetylene':Acetylene, 'amino':Amino, 'biphenyl':Biphenyl,
             'carboxyl':Carboxyl, 'cyano':Cyano, 'cyclopropyl':Cyclopropyl,
             'ethylene':Ethylene, 'fluorophenyl':Fluorophenyl,
             'formyl':Formyl, 'hydroxyl':Hydroxyl, 'isopropyl':Isopropyl,
             'methoxy':Methoxy, 'methyl':Methyl, 'nitro':Nitro,
             'nitrophenyl':Nitrophenyl, 'pentafluorophenyl':
             Pentafluorophenyl, 'perfluoromethyl':Perfluoromethyl,
             'phenyl':Phenyl, 'pyrrole':Pyrrole, 'triazole':Triazole,
             'difluoromethyl':Difluoromethyl, 'phenol':Phenol,
             'toluene':Toluene, 'benzoicacid':Benzoicacid,
             'isopropylbenzene':Isopropylbenzene, 'amide':Amide,
             'cyclohexyl':Cyclohexyl, 'sulfhydryl':Sulfhydryl,
             'nitroso':Nitroso}
        tgroup = terminal_group_dict[terminal_group]()

        alkane = Alkane(chain_length, cap_front=False, cap_end=False)
        self.add(alkane, 'alkane')
        self.add(tgroup, 'terminal_group')
        mb.force_overlap(self['alkane'], self['alkane']['up'], 
                         self['terminal_group']['down'])
        silane = Silane()
        self.add(silane, 'silane')
        mb.force_overlap(self['silane'], self['silane']['up'], self['alkane']['down'])

        self.add(silane['down'], 'down', containment=False)

if __name__ == "__main__":
    chain = Alkylsilane(chain_length=8, terminal_group='toluene')
