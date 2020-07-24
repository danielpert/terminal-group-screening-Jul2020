import mbuild as mb
import numpy as np


class Amino(mb.Compound):
    """ A amino group. """
    def __init__(self):
        super(Amino, self).__init__()

        mb.load('amino.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)

        self.add(mb.Port(anchor=self[0], orientation=[0, -1, 0], separation=0.07),
            'down')

if __name__ == '__main__':
    amino = Amino()
    amino.save('amino-test.mol2', overwrite=True)
