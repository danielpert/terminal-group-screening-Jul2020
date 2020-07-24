import mbuild as mb
import numpy as np


class Acetylene(mb.Compound):
    """ A acetylene group. """
    def __init__(self):
        super(Acetylene, self).__init__()

        mb.load('acetylene.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)

        self.add(mb.Port(anchor=self[0], orientation=[0, -1, 0], separation=0.0735),
            'down')

if __name__ == '__main__':
    acetylene = Acetylene()
    acetylene.save('acetylene-test.mol2', overwrite=True)
