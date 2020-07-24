import mbuild as mb
import numpy as np


class Methyl(mb.Compound):
    """ A methyl group. """
    def __init__(self):
        super(Methyl, self).__init__()

        mb.load('methyl.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)

        self.add(mb.Port(anchor=self[0], orientation=[0, -1, 0], separation=0.07),
            'down')

if __name__ == '__main__':
    methyl = Methyl()
    methyl.save('methyl-test.mol2', overwrite=True)
