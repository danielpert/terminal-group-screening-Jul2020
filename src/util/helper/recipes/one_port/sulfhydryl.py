import mbuild as mb
import numpy as np


class Sulfhydryl(mb.Compound):
    """ A sulfhydryl group. """
    def __init__(self):
        super(Sulfhydryl, self).__init__()

        mb.load('sulfhydryl.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)
        
        self.add(mb.Port(anchor=self[0], orientation=[0, -1, 0], separation=0.09),
            'down')

if __name__ == '__main__':
    sulfhydryl = Sulfhydryl()
    sulfhydyl.save('sulfhydryl-test.mol2', overwrite=True)