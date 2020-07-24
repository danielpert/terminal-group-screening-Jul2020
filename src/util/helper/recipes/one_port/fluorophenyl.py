import mbuild as mb
import numpy as np


class Fluorophenyl(mb.Compound):
    """ A phenyl group with a fluorine para to the chain attachment. """
    def __init__(self):
        super(Fluorophenyl, self).__init__()

        mb.load('fluorophenyl.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)

        self.add(mb.Port(anchor=self[0], orientation=[0, -1, 0], separation=0.07),
            'down')

if __name__ == '__main__':
    fluorophenyl = Fluorophenyl()
    fluorophenyl.save('fluorophenyl-test.mol2', overwrite=True)
