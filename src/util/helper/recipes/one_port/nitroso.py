import mbuild as mb
import numpy as np


class Nitroso(mb.Compound):
    """ A nitroso group. """
    def __init__(self):
        super(Nitroso, self).__init__()

        mb.load('nitroso.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)
        
        self.add(mb.Port(anchor=self[0], orientation=[0, -1, 0], separation=0.074),
            'down')

if __name__ == '__main__':
    nitroso = Nitroso()
    nitroso.save('nitroso-test.mol2', overwrite=True)