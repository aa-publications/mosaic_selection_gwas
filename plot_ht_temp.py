#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2022-01-16 22:17:18



import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path
from glob import glob
DATE = datetime.now().strftime('%Y-%m-%d')


sys.path.append('/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/scripts/1_manuscript/2022_10_29_dataFreeze/0_common_plotting_scripts')
from func_plot_bkg_and_trait_avg import plot_bkg_and_trait_avg, calc_p_adj, make_n_loci_and_label
from func_plot_strat_by_bkg_and_trait_avg import plot_stratified_enrichments_clean, make_clean_labels

###
# plotting imports
###
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as font_manager
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import seaborn as sns
import proplot

%matplotlib inline

#
# fpath = '/dors/capra_lab/users/abraha1/conda/envs/py36_r_ml/lib/python3.6/site-packages/matplotlib/mpl-data/fonts/ttf/Arial.ttf'
# font_dirs = ['/dors/capra_lab/users/abraha1/conda/envs/py36_r_ml/lib/python3.6/site-packages/matplotlib/mpl-data/fonts/ttf', ]
# font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
# for font_file in font_files:
#     font_manager.fontManager.addfont(font_file)
# mpl.rcParams['font.family'] = 'Arial'
# mpl.rcParams.update(mpl.rcParamsDefault)
# mpl.rcParams['pdf.fonttype'] = 42
# mpl.rcParams['ps.fonttype'] = 42
### PATHS

### DO NOT USE THESE PATHS IF ONLY RECREATING THE FIGURES
# ROOT=Path("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/data/consortia_traits/height")
# ROOT1=Path("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/scripts/1_manuscript/2022_09_04_plos_genetics_revisions/height_w_tighter_gene_control/gsel_output")


DATA_FREEZE_DIR=Path("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/scripts/1_manuscript/2022_10_29_dataFreeze/supFig4_dir")
strat_beta_files = {'bolt':DATA_FREEZE_DIR.joinpath("bolt_ht_enrich_strat_by_beta_summary_exp1_indep_partitions.tsv"),
                    'gene_density10':DATA_FREEZE_DIR.joinpath("gene_dens_10_enrich_strat_by_beta_summary_exp1_indep_partitions.tsv"),
                    'gene_density25':DATA_FREEZE_DIR.joinpath("gene_dens_25_enrich_strat_by_beta_summary_exp1_indep_partitions.tsv"),
                    'gene_dist10':DATA_FREEZE_DIR.joinpath("gene_dist_10_per_enrich_strat_by_beta_summary_exp1_indep_partitions.tsv"),
                    'gene_dist25':DATA_FREEZE_DIR.joinpath("gene_dist_25_per_enrich_strat_by_beta_summary_exp1_indep_partitions.tsv")}

# labels and outputs
# OUTPUT_DIR = Path("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/scripts/1_manuscript/2022_09_04_plos_genetics_revisions/height_w_tighter_gene_control/figures")



def svfig(label): return plt.savefig(
    OUTPUT_DIR.joinpath(f"{DATE}_{label}.pdf"))


def wrapper_make_clean_partition_labels(temp_bolt_df):

    bolt_df = temp_bolt_df.copy()
    orderded_labels, label2_number_dict = make_clean_labels(bolt_df, 'indep')

    bar_df = bolt_df.loc[:, ['annotation', 'emp_pval', 'enrich_per_mean_diff_by_genomstd', 'mean_trait_loci', 'partition_id','trait']].copy()
    bar_df['partition_bin'] = bar_df['partition_id'].map(label2_number_dict)
    bar_df['partition_label'] = bar_df['partition_bin'].map(orderded_labels)


    bar_df['partition_label'].unique()
    [ "{:.2f}".format(float(x.split(",")[0][1:])) for x in bar_df['partition_label'].unique()]
    [ "{:.2f}".format(float(x.split(",")[1][:-1])) for x in bar_df['partition_label'].unique()]
    [f"{x} to {y}" for x,y in zip(['-0.17', '-0.02', '-0.01', '0.01', '0.02'], ['-0.02', '-0.01', '0.01', '0.02', '0.14'])]

    clean_partition_label_dict =  { '(-1.7E-01, -2.1E-02]': '-0.17 to -0.02',
                                     '(-2.1E-02, -1.2E-02]': '-0.02 to -0.01',
                                     '(-1.2E-02, 1.1E-02]': '-0.01 to 0.01',
                                     '(1.1E-02, 2.0E-02]': '0.01 to 0.02',
                                     '(2.0E-02, 1.4E-01]': '0.02 to 0.14'}

    bar_df['clean_partitition_label'] = bar_df['partition_label'].map(clean_partition_label_dict)

    return bar_df
# %%
# -----------
# main
# ----------

annos = ['argweave', 'betascore',
        'linsigh', 'phastCon100', 'phyloP100',
         'fst_eas_afr', 'fst_eur_afr', 'fst_eur_eas',
        'xpehh_afr2_eas', 'xpehh_afr2_eur', 'xpehh_eas_eur']


clean_anno_dict ={'argweave':'ARGweaver',
                    'betascore':'Beta Score',
                    'phastCon100':'PhastCons',
                    'phyloP100':'PhyloP',
                    'linsigh':'LINSIGHT',
                    'fst_eur_afr':'Fst\nafr-eur',
                    'fst_eas_afr':'Fst\nafr-eas',
                    'fst_eur_eas':'Fst\neas-eur',
                    'xpehh_eas_eur':'XP-EHH\neas-eur',
                    'xpehh_afr2_eas':'XP-EHH\nafr-eas',
                    'xpehh_afr2_eur':'XP-EHH\nafr-eur'
                    }

keepcols=['annotation', 'emp_pval', 'enrich_per_mean_diff_by_genomstd', 'error', 'mean_trait_loci', 'mean_matched_across_trait_loci',
            'matched_5th', 'matched_95th','median_matched_loci',  'n_lead_loci_final', 'n_lead_snp', 'partition_id',
            'percent_matched_and_ld_removed', 'percent_trait_and_ld_snps_removed', 'trait_name']

# load stratified file
label='gene_dist25'
bolt_df = pd.read_csv( strat_beta_files[label], sep="\t",  usecols=keepcols)
bolt_df.rename(columns={'trait_name':'trait'}, inplace=True)


# %%
# compare to ukbb bolt ht
all_datasets_df = pd.DataFrame()
all_mean_trait_df = pd.DataFrame()
for tr in [ "bolt","gene_density10","gene_density25","gene_dist10","gene_dist25"]:

    kcol= ['annotation', 'emp_pval','enrich_per_mean_diff_by_genomstd', 'mean_trait_loci', 'partition_id', 'trait_name']
    temp_df = pd.read_csv( strat_beta_files[tr], sep="\t",  usecols=kcol)
    temp_df.rename(columns={'trait_name':'trait'}, inplace=True)
    bar_df = wrapper_make_clean_partition_labels(temp_df)
    all_datasets_df = pd.concat([all_datasets_df, bar_df.pivot(index=['annotation', 'partition_bin'],  columns='trait', values='enrich_per_mean_diff_by_genomstd')], axis=1)


    all_mean_trait_df = pd.concat([all_mean_trait_df, bar_df.pivot(index=['annotation', 'partition_bin'],  columns='trait', values='mean_trait_loci')], axis=1)


    # all_datasets_df = pd.concat([all_datasets_df, bar_df.pivot(index=['annotation', 'partition_label', 'partition_bin'],  columns='trait', values='enrich_per_mean_diff_by_genomstd')], axis=1)


all_datasets_df.sort_values(['annotation', 'partition_bin'], inplace=True)
all_mean_trait_df.sort_values(['annotation', 'partition_bin'], inplace=True)



# convert to long df
all_datasets_df.reset_index(inplace=True)
long_all_df = pd.melt(all_datasets_df, id_vars=['annotation', 'partition_bin']).copy()
long_all_df.rename(columns={'value':'Enrichment'}, inplace=True)

all_mean_trait_df.reset_index(inplace=True)
long_trait_mean_df = pd.melt(all_mean_trait_df, id_vars=['annotation', 'partition_bin']).copy()
long_trait_mean_df.rename(columns={'value':'Evolutionary measure'}, inplace=True)


# %%
# plot all datasets for ENRICHMENT

sns.set(style="ticks",  font_scale=0.8, rc={"figure.figsize": (14, 5)})

fig, axs = plt.subplots(nrows=2, ncols=6, sharex=True, sharey=False)

axs = axs.ravel()
for id, anno in enumerate(annos):
    if id >= 5:
        id = id+1

    ax = axs[id]
    sns.stripplot(data=long_all_df.query('annotation==@anno').query("trait!='bolt_ht'"),  x='partition_bin', y='Enrichment', hue='trait', ax=axs[id], palette={'bolt_ht':'red', 'gene_dens_10':'lightgreen', 'gene_dens_25':'darkgreen', 'gene_dist_10_per':'lightblue','gene_dist_25_per':'darkblue'}, alpha=0.7, marker='o', edgecolor='k', linewidth=1)


    sns.stripplot(data=long_all_df.query('annotation==@anno').query("trait=='bolt_ht'"),  x='partition_bin', y='Enrichment', hue='trait', ax=axs[id], size=10,
                 palette={'bolt_ht':'red', 'gene_dens_10':'lightgreen', 'gene_dens_25':'darkgreen', 'gene_dist_10_per':'lightblue','gene_dist_25_per':'darkblue'}, alpha=0.5, marker='x', edgecolor='k', linewidth=2)

    ax.set_ylabel("")
    if id == 0 or id == 6:
        ax.set_ylabel("Enrichment")

    ax.set_xlabel("")
    if id >5:
        ax.set_xlabel("Effect Size Bins")






    ax.set_title(clean_anno_dict[anno], fontweight='bold')
    axs[id].legend().set_visible(False)
    if id == 11:
        axs[id].legend().set_visible(True)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
# svfig('enrich_across_matching_params')


# %%
