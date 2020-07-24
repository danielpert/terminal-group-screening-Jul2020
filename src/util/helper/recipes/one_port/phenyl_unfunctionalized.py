from __future__ import division

import mbuild as mb


class PhenylUnfunctionalized(mb.Compound):
    """An unfunctionalized phenyl group. """
    def __init__(self):
        super(PhenylUnfunctionalized, self).__init__()

        mb.load('phenyl_unfunctionalized.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self.pos)

        for particle in self.particles():
            neighbors = self.bond_graph.neighbors(particle)
            midpoint = [(p1 + p0) / 2 for p1, p0 in zip(neighbors[1].pos,
                                                        neighbors[0].pos)]
            port_direction = particle.pos - midpoint
            self.add(mb.Port(anchor=particle, 
                             orientation=port_direction,
                             separation=0.075),
                     label='port{}'.format(len(self.referenced_ports())))

if __name__ == '__main__':
    molecule = PhenylUnfunctionalized()
