from __future__ import division

from math import cos, sin
import numpy as np

import mbuild as mb


class RandomHemispherePattern(mb.Pattern):
    def __init__(self, n, seed=12345, **kwargs):

        np.random.seed(seed)

        u = (np.random.random(n) * 2) - 1
        theta = np.random.random(n) * 2 * np.pi

        x = [(1 - a**2)**0.5 * cos(rot) for a,rot in zip(u,theta)]
        y = [(1 - a**2)**0.5 * sin(rot) for a,rot in zip(u,theta)]
        z = abs(u)

        points = np.stack((x, y, z), axis=1)

        super(RandomHemispherePattern, self).__init__(points=points, **kwargs)

if __name__ == "__main__":
    pattern = RandomHemispherePattern(n=1000, radius=1.0)

    lj_proto = mb.Compound(name='LJ')
    lj_box = mb.Compound()
    for pos in pattern:
        lj_particle = mb.clone(lj_proto)
        lj_particle.translate(pos)
        lj_box.add(lj_particle)
    lj_box.save('hemisphere.xyz', overwrite=True)
