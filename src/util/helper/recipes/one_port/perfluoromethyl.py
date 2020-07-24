import mbuild as mb
import numpy as np


class Perfluoromethyl(mb.Compound):
    """ A fluorinated methyl group. """
    def __init__(self):
        super(Perfluoromethyl, self).__init__()

        mb.load('perfluoromethyl.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)

        self.add(mb.Port(anchor=self[0], orientation=[0, -1, 0], separation=0.07),
            'down')

if __name__ == '__main__':
    pfmethyl = Perfluoromethyl()
    pfmethyl.save('perfluoromethyl-test.mol2', overwrite=True)
