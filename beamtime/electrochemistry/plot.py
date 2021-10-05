import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def plot_gc(cycles, which_cycles='all', chg=True, dchg=True, colours=None, x='C', y='U'):

	fig, ax = prepare_gc_plot()


	if which_cycles == 'all':
		which_cycles = [i for i, c in enumerate(cycles)]

	if not colours:
		chg_colour = (40/255, 70/255, 75/255) # Dark Slate Gray #28464B
		dchg_colour = (239/255, 160/255, 11/255) # Marigold #EFA00B
	
	

	for i, cycle in cycles:
		if i in which_cycles:
			if chg:
				cycle[0].plot(ax=ax)




	






def prepare_gc_plot(figsize=(14,7), dpi=None):

	fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
	

	return fig, ax