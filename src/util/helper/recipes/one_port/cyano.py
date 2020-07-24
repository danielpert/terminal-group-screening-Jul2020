import mbuild as mb
import numpy as np


class Cyano(mb.Compound):
    """ A cyano group. """
    def __init__(self):
        super(Cyano, self).__init__()

        mb.load('cyano.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)

        self.add(mb.Port(anchor=self[0], orientation=[0, -1, 0], separation=0.07),
            'down')

if __name__ == '__main__':
    cyano = Cyano()
    cyano.save('cyano-test.mol2', overwrite=True)
