#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2022-05-07 12:41:20

import matplotlib.font_manager as font_manager
import matplotlib as mpl
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import seaborn as sns
import matplotlib.pyplot as plt
import proplot
import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path
from glob import glob
DATE = datetime.now().strftime('%Y-%m-%d')


sys.path.append(
    '/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/scripts/0_common_plotting_scripts')
from func_plot_bkg_and_trait_avg import plot_bkg_and_trait_avg, calc_p_adj, make_n_loci_and_label

###
# plotting imports
###

# get_ipython().run_line_magic('config', "InlineBackend.figure_format='retina'")
mpl.rcParams.update(mpl.rcParamsDefault)

%matplotlib inline
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
# %config InlineBackend.figure_format='retina'


fpath = '/dors/capra_lab/users/abraha1/conda/envs/py36_r_ml/lib/python3.6/site-packages/matplotlib/mpl-data/fonts/ttf/Arial.ttf'
font_dirs = ['/dors/capra_lab/users/abraha1/conda/envs/py36_r_ml/lib/python3.6/site-packages/matplotlib/mpl-data/fonts/ttf', ]
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
for font_file in font_files:
    font_manager.fontManager.addfont(font_file)
mpl.rcParams['font.family'] = 'Arial'



### PATHS
DATA_FREEZE_DIR = Path("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/scripts/1_manuscript/2022_10_29_dataFreeze")
PLOT_DATA = DATA_FREEZE_DIR.joinpath('2022-05-07_fig2a_height_bkg_and_trait_summary.tsv')

# %%
# -----------
# get data ready
# -----------
enrich_df = pd.read_csv( PLOT_DATA, sep="\t")
trait2num_loci_dict, trait2clean_label_dict = make_n_loci_and_label(enrich_df)


trait2clean_label_dict={ 'berg_paths':  'Berg 2019',
                         'neale_paths': 'Neale 2017',
                         'giant_paths': 'GIANT 2018',
                         'bolt_paths':  'Loh 2017'}

annos = ['argweave', 'betascore',
         'linsigh', 'phastCon100', 'phyloP100',
         'fst_eas_afr', 'fst_eur_afr', 'fst_eur_eas',
         'xpehh_afr2_eas', 'xpehh_afr2_eur', 'xpehh_eas_eur']

# %%
axs = plot_bkg_and_trait_avg(annos, enrich_df, trait2num_loci_dict, trait2clean_label_dict, fig_height=2.,
                             fig_length_mult=0.75, lg_bbox=(3,0.4), ncol=1, xlabel_fontsize=8, markersize=7, lwd=1.4)


[ax.legend().set_visible(False) for ax in axs]
[ax.tick_params(axis='x', which='major', length=5, width=0.75) for ax in axs]
[ax.spines['bottom'].set_linewidth(0.75) for ax in axs]
axs[0].set_ylabel('Height GWAS', fontsize=8)

fig = plt.gcf()
fig.text(0.5, 0.04, 'Evolutionary Measures', ha='center', fontsize=8)

### make custome legend
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

legend_elements = [Line2D([0], [0], color='gray', lw=1, label='Matched\nBackground'),
                   # Line2D([0], [0], marker='*', label='Height-associated regions',  markersize=0, lw=0),
                   Line2D([0], [0], marker='*', color='indianred', label='Trait (p<0.05)', markerfacecolor='indianred', markersize=8, lw=0),
                   Line2D([0], [0], marker='*', color='k', label='Trait (p≥0.05)', markerfacecolor='k', markersize=8, lw=0)]

# Create the figure
axs[0].legend(handles=legend_elements, bbox_to_anchor=(-1.9,-0.75), loc="lower left", fontsize=6, fancybox=True, frameon=True)




plt.subplots_adjust(top=0.99, bottom=0.45, left=0.12, right=0.99, wspace=0.25)



