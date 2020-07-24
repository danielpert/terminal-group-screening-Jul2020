import inspect
import os

import mbuild as mb
from .silica_interface_carve import SilicaInterfaceCarve


class SilicaInterface(SilicaInterfaceCarve):
    """Recipe for generating an amorphous silica interface.

    This is essentially a wrapper around the `SilicaInterfaceCarve` recipe,
    which is an extension of that included with
    the mBuild package (https://github.com/mosdef-hub/mbuild). Please refer to the
    mBuild documentation for further details into how the interface is carved.

    Parameters
    ----------
    tile_x : int, optional, default=1
        Number of times to replicate the surface in the x dimension. The default
        length in the x dimension is 5nm, so increasing tile_x will lead to higher
        multiples of this number.
    tile_y : int, optional, default=1
        Number of times to replicate the surface in the y dimension. The default
        length in the y dimension is 5nm, so increasing tile_y will lead to higher
        multiples of this number.
    thickness : float, optional, default=1.0
        Desired thickness of the surface (in nm; not including oxygen layers on the
        top and bottom of the surface)
    seed : int, optional, default=12345
        Seed for the random number generator used in bridging surface oxygens
    """
    def __init__(self, tile_x=1, tile_y=1, thickness=1.0, seed=12345):
        thickness = float(thickness)
        try:
            from mbuild.lib.bulk_materials import AmorphousSilica
            super(SilicaInterface, self).__init__(bulk_silica=AmorphousSilica(),
                tile_x=tile_x, tile_y=tile_y, thickness=thickness, seed=seed)
        except:
            from mbuild.lib.bulk_materials import AmorphousSilicaBulk
            super(SilicaInterface, self).__init__(bulk_silica=AmorphousSilicaBulk(),
                tile_x=tile_x, tile_y=tile_y, thickness=thickness, seed=seed)

if __name__ == "__main__":
    silica_interface = SilicaInterface(thickness=1.2, seed=10)
