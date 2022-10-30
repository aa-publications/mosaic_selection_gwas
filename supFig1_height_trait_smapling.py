#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2022-01-26 07:30:09


import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path

DATE = datetime.now().strftime('%Y-%m-%d')


import matplotlib.pyplot as plt
import seaborn as sns


## load matplotlib for manuscripts
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
# %config InlineBackend.figure_format='retina'


# fpath = '/dors/capra_lab/users/abraha1/conda/envs/py36_r_ml/lib/python3.6/site-packages/matplotlib/mpl-data/fonts/ttf/Arial.ttf'
# font_dirs = ['/dors/capra_lab/users/abraha1/conda/envs/py36_r_ml/lib/python3.6/site-packages/matplotlib/mpl-data/fonts/ttf', ]
# font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
# for font_file in font_files:
#     font_manager.fontManager.addfont(font_file)
# mpl.rcParams['font.family'] = 'Arial'
# mpl.rcParams.update(mpl.rcParamsDefault)
#

###
###    paths
###


### ORIGNIAL PATH TO FILES
# ROOT=Path("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/")
# DATA=ROOT.joinpath("data/consortia_traits/height/boltlmm_height/sample_lead_loci")
# with_replacement_files = list(DATA.glob("*_w_replacement.tsv"))
# no_replacement_files = list(DATA.glob("*_without_replacement.tsv"))
# OUTPUT_DIR=Path("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/scripts/1_manuscript/figs_height/")

DATA_FREEZE_DIR = Path("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/scripts/1_manuscript/2022_10_29_dataFreeze")


# -----------
# functions
# -----------
def load_all_n_sampled_experiments(files_to_load):

    combined_df = pd.DataFrame()
    for rfile in files_to_load:

        usecols = ['annotation', 'emp_pval', 'enrich_per_mean_diff_by_genomstd',
           'matched_5th',  'matched_95th',
            'mean_trait_loci', 'mean_matched_across_trait_loci',
           'median_matched_loci',  'n_lead_loci_final',
           'n_sample_lead_loci',
           'regions_summary', 'replacement',
           'trait_name', 'trait_summary']
        n_sampled_snps  = int(rfile.stem.split("_sample_")[1].split("_lead")[0])
        df = pd.read_csv( rfile, sep="\t", usecols=usecols)
        df['ideal_n_sampled'] = n_sampled_snps
        combined_df = combined_df.append(df)

    return combined_df

def svfig(label): return plt.savefig(
    OUTPUT_DIR.joinpath(f"{DATE}_{label}.pdf"))

# %%
# -----------
# main
# -----------

### SKIP IF TRYING TO REPRODUCE FIGURE
# no_replace_df = load_all_n_sampled_experiments(no_replacement_files)

### SKIP IF TRYING TO REPRODUCE FIGURE
# remove those with only 10 trait loci
# no_replace_df = no_replace_df[no_replace_df['n_lead_loci_final']!=10].copy()


# %%
# save to freeze directory
# no_replace_df.to_csv(DATA_FREEZE_DIR.joinpath(f"{DATE}_supFig1_height_trait_sampling.tsv"), sep="\t",index=False)

### LOAD IF ONLY TRYING TO REPRODUCE FIGURE

no_replace_df = pd.read_csv(DATA_FREEZE_DIR.joinpath(f"2022-05-07_supFig1_height_trait_sampling.tsv.gz"), sep="\t")
# %%
sns.set(style="ticks",  font_scale=1.0, rc={"figure.figsize": (18, 6)})
fig, axs = plt.subplots(nrows=2, ncols=int(np.ceil(no_replace_df['annotation'].nunique()/2)), sharex=True, sharey=False)

axs=axs.ravel()


for axi, anno in enumerate(no_replace_df['annotation'].unique()):
    ax = axs[axi]

    plot_df = no_replace_df.query("annotation==@anno").sort_values('n_lead_loci_final')
    ax.plot(plot_df['n_lead_loci_final'], plot_df['mean_trait_loci'], '.-')
    ax.plot(plot_df['n_lead_loci_final'], plot_df['mean_matched_across_trait_loci'], 'k.-')
    ax.plot(plot_df['n_lead_loci_final'], plot_df['matched_5th'], marker='.', linewidth=0.5, markersize=1, color='gray')
    ax.plot(plot_df['n_lead_loci_final'],plot_df['matched_95th'], marker='.', linewidth=0.5, markersize=1, color='gray')

    ax.fill_between(plot_df['n_lead_loci_final'],plot_df['matched_5th'], plot_df['matched_95th'], color='gray', alpha=0.5)

    ax.set_title(anno)

    ax.set_xticks(plot_df['n_lead_loci_final'].unique())
    ax.set_xticklabels(plot_df['n_lead_loci_final'].unique(), rotation=330)

plt.tight_layout()
# svfig('undersampling_without_replacement_height_boltlmm')


# %%



