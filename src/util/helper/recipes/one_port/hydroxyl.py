import mbuild as mb
import numpy as np


class Hydroxyl(mb.Compound):
    """ A hydroxyl group. """
    def __init__(self):
        super(Hydroxyl, self).__init__()

        mb.load('hydroxyl.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)

        self.add(mb.Port(anchor=self[0], orientation=[0, -1, 0], separation=0.07),
            'down')

if __name__ == '__main__':
    hydroxyl = Hydroxyl()
    hydroxyl.save('hydroxyl-test.mol2', overwrite=True)
