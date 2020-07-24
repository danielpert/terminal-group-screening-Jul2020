import mbuild as mb
from . import Phenyl

class Biphenyl(mb.Compound):
    """ A biphenyl group. """
    def __init__(self):
        super(Biphenyl, self).__init__()

        phenyl1 = Phenyl()
        phenyl2 = Phenyl()

        phenyl1.remove(phenyl1[6])
        phenyl1.add(mb.Port(anchor=phenyl1[5], separation=0.07), 'up')

        self.add(phenyl1, 'phenyl1')
        self.add(phenyl2, 'phenyl2')
        mb.force_overlap(self['phenyl2'], self['phenyl2']['down'], 
                         self['phenyl1']['up'])

        self.add(phenyl1['down'], 'down', containment=False)

if __name__ == '__main__':
    biphenyl = Biphenyl()
    biphenyl.save('biphenyl-test.mol2', overwrite=True)
