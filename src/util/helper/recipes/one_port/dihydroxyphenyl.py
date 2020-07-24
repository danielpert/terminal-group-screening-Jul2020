from __future__ import division

import mbuild as mb
from mbuild.lib.atoms.h import H
from mbuild.lib.moieties.ch2 import CH2

from atools.lib.moieties.one_port.oh import OH
from atools.lib.moieties.one_port.phenyl_unfunctionalized import PhenylUnfunctionalized


class DihydroxyPhenyl(mb.Compound):
    """A phenyl group with hydroxyls at the 2 and 3 positions. """
    def __init__(self):
        super(DihydroxyPhenyl, self).__init__()

        phenyl_unfunctionalized = PhenylUnfunctionalized()
        hydroxyl = OH()
        ch2 = CH2()
        hydrogen = H()

        self.add(phenyl_unfunctionalized, 'ring')
        self.add(ch2, 'ch2')
        mb.force_overlap(self['ch2'], self['ch2']['up'], self['ring']['port0'])
        for i, port_number in enumerate([2, 3]):
            oh_clone = mb.clone(hydroxyl)
            self.add(oh_clone, 'oh_{}'.format(i))
            mb.force_overlap(self['oh_{}'.format(i)],
                             self['oh_{}'.format(i)]['up'],
                             self['ring']['port{}'.format(port_number)])
        for i, port_number in enumerate([1, 4, 5]):
            h_clone = mb.clone(hydrogen)
            self.add(h_clone, 'h_{}'.format(i))
            mb.force_overlap(self['h_{}'.format(i)],
                             self['h_{}'.format(i)]['up'],
                             self['ring']['port{}'.format(port_number)])

        self.add(ch2['down'], 'down', containment=False)

if __name__ == '__main__':
    molecule = DihydroxyPhenyl()
