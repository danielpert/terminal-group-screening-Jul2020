import numpy as np
import mbuild as mb

class DualSurface(mb.Compound):
    """ A recipe for creating a system with two opposing surfaces.
    """

    def __init__(self, top, bottom, separation=0.8):
        super(DualSurface, self).__init__()
        
        top.spin(np.pi, around=[0, 1, 0])
        top_box = top.boundingbox
        bot_box = bottom.boundingbox

        top_of_bot = bot_box.maxs[2]
        bot_of_top = top_box.mins[2]
        top.translate([0, 0, top_of_bot - bot_of_top + separation])

        self.add(bottom)
        self.add(top)
