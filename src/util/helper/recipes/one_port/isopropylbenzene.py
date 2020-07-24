import mbuild as mb
import numpy as np


class Isopropylbenzene(mb.Compound):
    """ A isopropylbenzene group. """
    def __init__(self):
        super(Isopropylbenzene, self).__init__()

        mb.load('isopropylbenzene.pdb', compound=self,
                relative_to_module=self.__module__)
        # pop off bottom hydrogen on benzene ring
        direction = self[18].xyz - self[6].xyz
        self.remove(self[18])

        # add port anchored to newly hydrogen-less carbon in benzene ring
        self.add(
                mb.Port(anchor=self[6], orientation=direction.tolist()[0], separation=0.07), 'down')


if __name__ == '__main__':
    isopropylbenzene = Isopropylbenzene()
    isopropylbenzene.save('isopropylbenzene-test.mol2', overwrite=True)
