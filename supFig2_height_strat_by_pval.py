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
# -----------
# PICK ONE!
# -----------
### ORIGNIAL PATH TO FILES
# ROOT=Path("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/data/consortia_traits/height")
# strat_pval_files = {'berg':ROOT.joinpath("berg_et_al/doi_10.5061_dryad.mg1rr36__v1/by_chr/gsel_out/strat_enrich/berg_ht_enrich_strat_by_pval_summary_exp1_indep_partitions.tsv"),
#                     'neale': ROOT.joinpath("neale_UKBB/by_chr/strat_enrich/neale_ht_enrich_strat_by_pval_summary_exp1_indep_partitions.tsv"),
#                     'giant': ROOT.joinpath("giant_wood_2018/strat_enrich/giant_ht_enrich_strat_by_pval_summary_exp1_indep_partitions.tsv"),
#                     'bolt':ROOT.joinpath("boltlmm_height/strat_enrich/bolt_ht_enrich_strat_by_pval_summary_exp1_indep_partitions.tsv")}
#
#
# # labels and outputs
# OUTPUT_DIR = Path("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/scripts/1_manuscript/figs_height")


# location with file to remake figure
DATA_FREEZE_DIR = Path("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/scripts/1_manuscript/2022_10_29_dataFreeze")

def svfig(label): return plt.savefig(
    OUTPUT_DIR.joinpath(f"{DATE}_{label}.pdf"))


# %%
# -----------
# main
# ----------
keepcols=['annotation', 'emp_pval', 'enrich_per_mean_diff_by_genomstd', 'error', 'mean_trait_loci', 'mean_matched_across_trait_loci',
            'matched_5th', 'matched_95th','median_matched_loci',  'n_lead_loci_final', 'n_lead_snp', 'partition_id',
            'percent_matched_and_ld_removed', 'percent_trait_and_ld_snps_removed', 'trait_name']


### ORIGNIAL PATH TO FILES
### DO NOT LOAD THESE FILES IF ONLY TRYING TO REPRODUCE FIGURE
# # load stratified file
# bolt_df = pd.read_csv( strat_pval_files['bolt'], sep="\t",  usecols=keepcols)
# bolt_df.rename(columns={'trait_name':'trait'}, inplace=True)

annos = ['argweave', 'betascore',
        'linsigh', 'phastCon100', 'phyloP100',
         'fst_eas_afr', 'fst_eur_afr', 'fst_eur_eas',
        'xpehh_afr2_eas', 'xpehh_afr2_eur', 'xpehh_eas_eur']


# %% - plot only encirhment?
### ORIGNIAL PATH TO FILES
### DO NOT LOAD THESE FILES IF ONLY TRYING TO REPRODUCE FIGURE
# bar_df = bolt_df.loc[:, ['annotation', 'emp_pval', 'enrich_per_mean_diff_by_genomstd', 'mean_trait_loci', 'partition_id','trait']].copy()
# orderded_labels, label2_number_dict = make_clean_labels(bolt_df, 'indep')
# bar_df['partition_bin'] = bar_df['partition_id'].map(label2_number_dict)
# bar_df['partition_label'] = bar_df['partition_bin'].map(orderded_labels)
#
#
# bar_df['partition_label'].unique()
# output_path_r = str(OUTPUT_DIR)




bar_df = pd.read_csv(DATA_FREEZE_DIR.joinpath(f"2022-05-07_supFig2_height_strat_by_pval.tsv.gz"), sep="\t")

##
## grouped batplot per annotation with height = mean trait value, colored by enrichment
##

# %%
import rpy2.rinterface
%load_ext rpy2.ipython

# %%
%%R
library(ggplot2)

# %%
# %%R -i bar_df -i output_path_r -w 1000 -h 300 -u px
%%R -i bar_df  -w 1000 -h 300 -u px


unique(bar_df$partition_label)
lvls=c("(5.1E-09, 5.0E-08]",   "(1.1E-10, 5.1E-09]" , "(1.3E-13, 1.1E-10]", "(1.9E-20, 1.3E-13]",     "(0.0E+00, 1.9E-20]" )
bar_df$partition_label<- factor(bar_df$partition_label, ordered = TRUE, levels = lvls)

# %%
%%R
ggplot(data=bar_df, aes(x=partition_label, y=mean_trait_loci, fill=enrich_per_mean_diff_by_genomstd)) +
    geom_bar(stat="identity",   position=position_dodge()) + labs(fill="Enrichment") +
    facet_wrap(~annotation, scales="free_y", nrow=2) +
    scale_fill_gradient(low="gray", high="red") +
    theme(axis.text.x = element_text(angle = 330, hjust = 0, vjust = 0.5))

# svfile = file.path(output_path_r, sprintf("%s_bolt_height_strat_by_pvalue.pdf", Sys.Date()))
# ggsave(svfile, height=2, width=6.5)

# %%
### SAVE PLOTTING DATA TO DATA FREEZE
# bar_df.to_csv(DATA_FREEZE_DIR.joinpath(f"{DATE}_supFig2_height_strat_by_pval.tsv"), sep="\t", index=False)
