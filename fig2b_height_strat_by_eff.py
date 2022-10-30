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


sys.path.append(
    '/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/scripts/0_common_plotting_scripts')
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


fpath = '/dors/capra_lab/users/abraha1/conda/envs/py36_r_ml/lib/python3.6/site-packages/matplotlib/mpl-data/fonts/ttf/Arial.ttf'
font_dirs = ['/dors/capra_lab/users/abraha1/conda/envs/py36_r_ml/lib/python3.6/site-packages/matplotlib/mpl-data/fonts/ttf', ]
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
for font_file in font_files:
    font_manager.fontManager.addfont(font_file)
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
### PATHS
ROOT=Path("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/data/consortia_traits/height")
strat_beta_files = {'berg':ROOT.joinpath("berg_et_al/doi_10.5061_dryad.mg1rr36__v1/by_chr/gsel_out/strat_enrich berg_ht_enrich_strat_by_beta_summary_exp1_indep_partitions.tsv"),
                    'neale': ROOT.joinpath("neale_UKBB/by_chr/strat_enrich/neale_ht_enrich_strat_by_beta_summary_exp1_indep_partitions.tsv"),
                    'giant': ROOT.joinpath("giant_wood_2018/strat_enrich/giant_ht_enrich_strat_by_beta_summary_exp1_indep_partitions.tsv"),
                    'bolt':ROOT.joinpath("boltlmm_height/strat_enrich/bolt_ht_enrich_strat_by_beta_summary_exp1_indep_partitions.tsv")}

# labels and outputs
OUTPUT_DIR = Path("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/scripts/1_manuscript/figs_height")
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

# load stratified file
bolt_df = pd.read_csv( strat_beta_files['bolt'], sep="\t",  usecols=keepcols)
bolt_df.rename(columns={'trait_name':'trait'}, inplace=True)

annos = ['argweave', 'betascore',
        'linsigh', 'phastCon100', 'phyloP100',
         'fst_eas_afr', 'fst_eur_afr', 'fst_eur_eas',
        'xpehh_afr2_eas', 'xpehh_afr2_eur', 'xpehh_eas_eur']


# %% -- plot with bkg distribution
# _ = plot_stratified_enrichments_clean(annos, bolt_df,  plt_title=f"indep_bolt", savefile=False, w_anno=False, font_scale=0.8, xlabel_fontsize=10 , yticklabel_fontsize=10, fig_length_mult=1.5,  fig_height=3, lg_bbox=(2.6,0.4), ncol=1)
# plt.tight_layout()
# svfig("bolt_by_effect_size")

# %% - plot only encirhment?
bar_df = bolt_df.loc[:, ['annotation', 'emp_pval', 'enrich_per_mean_diff_by_genomstd', 'mean_trait_loci', 'partition_id','trait']].copy()

orderded_labels, label2_number_dict = make_clean_labels(bolt_df, 'indep')
bar_df['partition_bin'] = bar_df['partition_id'].map(label2_number_dict)
bar_df['partition_label'] = bar_df['partition_bin'].map(orderded_labels)


bar_df['partition_label'].unique()
# {:.8f}".format(float("8.99284722486562e-02"))
[ "{:.2f}".format(float(x.split(",")[0][1:])) for x in bar_df['partition_label'].unique()]
[ "{:.2f}".format(float(x.split(",")[1][:-1])) for x in bar_df['partition_label'].unique()]
[f"{x} to {y}" for x,y in zip(['-0.17', '-0.02', '-0.01', '0.01', '0.02'], ['-0.02', '-0.01', '0.01', '0.02', '0.14'])]

clean_partition_label_dict =  { '(-1.7E-01, -2.1E-02]': '-0.17 to -0.02',
                                 '(-2.1E-02, -1.2E-02]': '-0.02 to -0.01',
                                 '(-1.2E-02, 1.1E-02]': '-0.01 to 0.01',
                                 '(1.1E-02, 2.0E-02]': '0.01 to 0.02',
                                 '(2.0E-02, 1.4E-01]': '0.02 to 0.14'}

bar_df['clean_partitition_label'] = bar_df['partition_label'].map(clean_partition_label_dict)


# remove GERP
bar_df = bar_df.loc[~bar_df['annotation'].isin(['gerp'])].copy()
bar_df['annotation'].unique()

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

bar_df['clean_annotation'] = bar_df['annotation'].map(clean_anno_dict)
bar_df['clean_annotation'].unique()
output_path_r = str(OUTPUT_DIR)
# %%
# save to freeze
# bar_df.to_csv(DATA_FREEZE_DIR.joinpath(f"{DATE}_fig2b_height_strat_by_pval.tsv"), sep="\t", index=False)

# %%
import rpy2.rinterface
%load_ext rpy2.ipython

# %%
%%R
 # library("ggplot2")
install.packages("extrafont")
# %%
%%R -i bar_df -i output_path_r
unique(bar_df$partition_label)

# order x-tick-labels
lvls=c('-0.17 to -0.02', '-0.02 to -0.01', '-0.01 to 0.01','0.01 to 0.02', '0.02 to 0.14')
bar_df$clean_partitition_label<- factor(bar_df$clean_partitition_label, ordered = TRUE, levels = lvls)

# order facets
alvls=c( 'ARGweaver',  'Beta Score', 'LINSIGHT', 'PhastCons','PhyloP',
        'Fst\nafr-eas', 'Fst\nafr-eur', 'Fst\neas-eur',
        'XP-EHH\nafr-eas', 'XP-EHH\nafr-eur', 'XP-EHH\neas-eur' )
bar_df$clean_annotation<- factor(bar_df$clean_annotation, ordered = TRUE, levels = alvls)

# order
# arg, beta, LINSIGHT, phastcon, phylo, fst_eas-afr, Fst eur-afr, Fst eur-eas, XP-EHH, afr-eas, afr eur, eas-eur
# %%
%%R


# Stacked barplot with multiple groups
ggplot(data=bar_df, aes(x=clean_partitition_label, y=mean_trait_loci, fill=enrich_per_mean_diff_by_genomstd)) +
    geom_bar(stat="identity",   position=position_dodge()) + labs(fill="Enrichment") +
    facet_wrap(~clean_annotation, scales="free_y", nrow=2) +
    scale_fill_gradient(low="gray", high="red") + theme_classic() +
    xlab("GWAS Effect Size") + ylab("Evolutionary Measure") +
    theme(axis.text.x = element_text(angle = 300, hjust = 0, vjust = 0.5, size=7),
          axis.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, size=7),
          strip.background = element_rect(colour="white", fill="white", size=8),
          axis.title=element_text(size=8))



dev.new(width = 8.5, height = 3, unit="in", noRStudioGD = T);last_plot() #Set new windo size and replot whatever plot you just made.
# ggsave(svfile,width = dev.size()[1],height = dev.size()[2]);dev.off() #Save the plot and set the size using `dev.siz()` so you don't have to ever change that part and cannot possibly
# ggsave(svfile, height=3, width=8.5)

# %%