import mbuild as mb


class OH(mb.Compound):
    """A hydroxyl group. """
    def __init__(self):
        super(OH, self).__init__()

        mb.load('oh.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos) # Move oxygen to origin

        self.add(mb.Port(anchor=self[0], separation=0.07), 'up')

if __name__ == '__main__':
    m = OH()
