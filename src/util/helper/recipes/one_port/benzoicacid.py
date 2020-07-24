import mbuild as mb
import numpy as np


class Benzoicacid(mb.Compound):
    """ A benzoicacid group. """
    def __init__(self):
        super(Benzoicacid, self).__init__()

        mb.load('benzoicacid.pdb', compound=self,
                relative_to_module=self.__module__)
        self.translate(self[12].pos)
        # pop off bottom hydrogen on benzene ring
        direction = self[12].xyz - self[6].xyz
        self.remove(self[12])

        # add port anchored to newly hydrogen-less carbon in benzene ring
        self.add(
                mb.Port(anchor=self[6], orientation=direction.tolist()[0], separation=0.07)
                , 'down')


if __name__ == '__main__':
    benzoicacid = Benzoicacid()
    benzoicacid.save('benzoicacid-test.mol2', overwrite=True)
