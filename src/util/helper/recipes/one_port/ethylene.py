import mbuild as mb
import numpy as np


class Ethylene(mb.Compound):
    """ A ethylene group. """
    def __init__(self):
        super(Ethylene, self).__init__()

        mb.load('ethylene.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)

        self.add(mb.Port(anchor=self[0], orientation=[0, -1, 0], separation=0.07),
            'down')

if __name__ == '__main__':
    ethylene = Ethylene()
    ethylene.save('ethylene-test.mol2', overwrite=True)
