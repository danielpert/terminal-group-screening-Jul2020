import mbuild as mb
import numpy as np


class Cyclohexyl(mb.Compound):
    """ A cyclohexyl group. """
    def __init__(self):
        super(Cyclohexyl, self).__init__()

        mb.load('cyclohexyl.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)

        self.add(mb.Port(anchor=self[0], orientation=[-0.5730, -0.6480, -0.6720], separation=0.076),
            'down')

if __name__ == '__main__':
    cyclohexyl = Cyclohexyl()
    cyclohexyl.save('cyclohexyl-test.mol2', overwrite=True)
