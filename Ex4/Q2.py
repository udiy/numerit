"""

	Author: Udi Yosovzon
	
		
	Description: This python code evaluates Lagrange polynomial and plot it.
				 
	Execution: In cmd, navigate to this file directory, then type: "python Q2.py". IMPORTANT: this file has to be in the same folder as Q1

"""

import numpy as np
import matplotlib.pyplot as plt
from Q1 import lagrange


fx_mapping = { -1:0, 0:1, 2:9, 3:25, 4:67 }



def plot_polynomial_lagrange(fx_mapping, limits, png_name=None):
    
    """ calculate lagrange polynomial and plot them
    """
    llim, ulim = limits
    
    # get fx_mapping only the key:value pair within the limits
    fx_mapping_slice = {key:fx_mapping[key] for key in fx_mapping if key in range(llim,ulim+1)}
    x_values = fx_mapping_slice.keys()
    fx_values = fx_mapping_slice.values()
    
    x = np.arange(llim-1, ulim+1, 0.01)
    f = lagrange(fx_mapping_slice, x)

    plt.figure(figsize=(16,9))
    plt.plot(x, f)
    plt.plot(x_values,fx_values,"o",c="red", markersize=15)

    plt.grid()
    
    if png_name:
        plt.savefig(f"{png_name}.png", bbox_inches="tight")
    else:
        plt.show()


# plot_polynomial_lagrange(fx_mapping, (-1,4))
# plot_polynomial_lagrange(fx_mapping, (-1,2))
# plot_polynomial_lagrange(fx_mapping, (2,4))