import mbuild as mb
import numpy as np


class Toluene(mb.Compound):
    """ A toluene group. """
    def __init__(self):
        super(Toluene, self).__init__()

        mb.load('toluene.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[4].pos)
        # pop off bottom hydrogen on benzene ring
        direction = self[12].xyz - self[4].xyz
        self.remove(self[12])

        # add port anchored to newly hydrogen-less carbon in benzene ring
        self.add(
                mb.Port(anchor=self[4], orientation=direction.tolist()[0], separation=0.07), 'down')


if __name__ == '__main__':
    toluene = Toluene()
    toluene.save('toluene-test.mol2', overwrite=True)
