import os
from pkg_resources import resource_filename

import mbuild as mb

from .random_hemisphere_pattern import RandomHemispherePattern
from .monolayer import Monolayer

class SurfaceMonolayer(mb.Compound):

    def __init__(self, surface, chains, n_chains, seed, fractions=None,
                 backfill=None, **kwargs):
        super(SurfaceMonolayer, self).__init__()

        if surface.name == 'SilicaInterface' and n_chains > 0:
            pattern = mb.Random2DPattern(n_chains, seed=seed)
        elif surface.name == 'SilicaTip':
            shift = 0.25
            surface.translate_to([0, 0, 0])
            surface.translate([0, 0, -1 * min(surface.xyz[:,2]) - shift])
            radius = max(surface.xyz[:,2]) - min(surface.xyz[:,2])
            pattern = RandomHemispherePattern(n=n_chains, seed=seed)
            pattern.scale(radius)
        elif surface.name == 'SilicaAsperity':
            pass
        else:
            pattern = mb.Random2DPattern(n_chains, seed=seed)
        
        if chains and n_chains > 0:
            monolayer = Monolayer(surface=surface, chains=chains, pattern=pattern,
                                  fractions=fractions, backfill=backfill, **kwargs)
        else:
            monolayer = Monolayer(surface=surface, chains=backfill, guest_port_name='up', **kwargs)

        self.add(monolayer)

if __name__ == "__main__":
    from alkylsilane import Alkylsilane
    from silica_interface import SilicaInterface
    from mbuild.lib.atoms import H

    hydrogen = H()
    seed = 12345
    chain = Alkylsilane(chain_length=6, terminal_group='methyl')

    planar_surface = SilicaInterface(seed=seed)
    planar_monolayer = SurfaceMonolayer(surface=planar_surface, 
            chains=chain, n_chains=100, seed=seed, backfill=hydrogen)
