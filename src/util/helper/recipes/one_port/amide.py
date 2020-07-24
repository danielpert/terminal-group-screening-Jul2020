import mbuild as mb
import numpy as np


class Amide(mb.Compound):
    """ An amide group. """
    def __init__(self):
        super(Amide, self).__init__()

        mb.load('amide.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)

        self.add(mb.Port(anchor=self[0], orientation=[-0.5140, 0.4580, 0.7410], separation=0.075),
            'down')

if __name__ == '__main__':
    amide = Amide()
    amide.save('amide-test.mol2', overwrite=True)
