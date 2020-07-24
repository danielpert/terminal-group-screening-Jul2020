import mbuild as mb
import numpy as np


class Phenol(mb.Compound):
    """ A phenol group. """
    def __init__(self):
        super(Phenol, self).__init__()

        mb.load('phenol.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)
        # pop off bottom hydrogen on benzene ring
        direction = self[7].xyz - self[0].xyz
        self.remove(self[7])

        # add port anchored to newly hydrogen-less carbon in benzene ring
        self.add(
                mb.Port(anchor=self[0], orientation=direction.tolist()[0], separation=0.07), 'down')


if __name__ == '__main__':
    phenol = Phenol()
    phenol.save('phenol-test.mol2', overwrite=True)
