
#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2021-04-23 22:41:40

import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path

DATE = datetime.now().strftime('%Y-%m-%d')

import proplot
import matplotlib.pyplot as plt
import seaborn as sns
%matplotlib inline
%config InlineBackend.figure_format='retina'
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
### PATHS
DATA_FREEZE_DIR = Path("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/scripts/1_manuscript/2022_10_29_dataFreeze")
ANNO_FILE=DATA_FREEZE_DIR.joinpath("genome_wide_summary_of_annotations.tsv")


### path of original file
# ENRICH_DIR = Path("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/scripts/run_mosaic_round_3_2021_01_27/2021_03_07_updated_gsel_output/full_gsel_outputs")
# OUTPUT_DIR=Path("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/scripts/1_manuscript/figs_alzh")
# gwases = {"alz.IGAP_stage_1_2_combined":"enrich/set8/24162737.Alzheimer_disease.IGAP_stage_1_2_combined.format_extreme_regions_mean_enrichment_all_annotation.tsv",
#         "Proxy_and_clinically_alz":"enrich/setN/29777097.Proxy_and_clinically_diagnosed_Alzheimers_disease.4_UKB_IGAP_AD_meta_summary_output_June2019.format_extreme_regions_mean_enrichment_all_annotation.tsv",
#         "alz.Kunkle_etal_Stage1_results":"enrich/setH/30820047.Alzheimers.Kunkle_etal_Stage1_results.format_extreme_regions_mean_enrichment_all_annotation.tsv",
#         "Lambert-2022":"enrich/lambert_alzh_2022/gsel_output_extreme_regions_mean_enrichment_all_annotation.tsv",
#         "alz.GRACEStageI_dbGAP": "enrich/setC/31473137.Alzheimers_disease.GRACEStageI_dbGAP.format_extreme_regions_mean_enrichment_all_annotation.tsv"}

# address of original files
ENRICH_DIR = DATA_FREEZE_DIR.joinpath("fig4_enrich_dir")
gwases = {"alz.IGAP_stage_1_2_combined":"24162737.Alzheimer_disease.IGAP_stage_1_2_combined.format_extreme_regions_mean_enrichment_all_annotation.tsv",
        "Proxy_and_clinically_alz":"29777097.Proxy_and_clinically_diagnosed_Alzheimers_disease.4_UKB_IGAP_AD_meta_summary_output_June2019.format_extreme_regions_mean_enrichment_all_annotation.tsv",
        "alz.Kunkle_etal_Stage1_results":"30820047.Alzheimers.Kunkle_etal_Stage1_results.format_extreme_regions_mean_enrichment_all_annotation.tsv",
        "Lambert-2022":"gsel_output_extreme_regions_mean_enrichment_all_annotation.tsv",
        "alz.GRACEStageI_dbGAP": "31473137.Alzheimers_disease.GRACEStageI_dbGAP.format_extreme_regions_mean_enrichment_all_annotation.tsv"}


# sys.path.append("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/scripts/0_common_plotting_scripts")
sys.path.append(str(DATA_FREEZE_DIR.joinpath("0_common_plotting_scripts")))
from func_plot_bkg_and_trait_avg import plot_bkg_and_trait_avg

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

import matplotlib as mpl
# import matplotlib.font_manager as font_manager
# fpath='/dors/capra_lab/users/abraha1/conda/envs/py36_r_ml/lib/python3.6/site-packages/matplotlib/mpl-data/fonts/ttf/Arial.ttf'
# font_dirs = ['/dors/capra_lab/users/abraha1/conda/envs/py36_r_ml/lib/python3.6/site-packages/matplotlib/mpl-data/fonts/ttf', ]
# font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
# font_list = font_manager.createFontList(font_files)
# font_manager.fontManager.ttflist.extend(font_list)

mpl.rcParams['font.family'] = 'Arial'


# %%

anno_summary_df = pd.read_csv( ANNO_FILE, sep="\t")
mean_anno_dict = dict(zip(anno_summary_df['annotation'], anno_summary_df['mean']))


# %%

enrich_df = pd.DataFrame()
for key in gwases.keys():
    print(key)
    efile = ENRICH_DIR.joinpath(gwases[key])
    e_df = pd.read_csv( efile, sep="\t")
    e_df.rename(columns={'n_lead_loci_final':'n_lead_loci', 'emp_pval':'emp_pvalue'}, inplace=True)
    e_df['trait']  = key
    e_df['efile']  = efile.name
    enrich_df = enrich_df.append(e_df)

key='Lambert-2022'



# new
# Index(['annotation', 'emp_pval', 'enrich_per_mean_diff_by_genomstd',
#        'matched_25th', 'matched_5th', 'matched_75th', 'matched_95th',
#        'mean_matched_across_trait_loci', 'mean_trait_loci',
#        'median_matched_loci', 'median_trait_loci', 'n_lead_loci',
#        'n_matched_and_ld_snps_removed', 'n_matched_sets',
#        'n_no_na_matched_and_ld_snps', 'n_no_na_trait_and_ld_snps',
#        'n_og_matched_and_ld_snps', 'n_og_trait_and_ld_snps',
#        'n_trait_and_ld_snps_removed', 'percent_matched_and_ld_removed',
#        'percent_trait_and_ld_snps_removed', 'regions_summary',
#        'std_matched_across_trait_loci', 'trait_name', 'trait_summary', 'trait',
#        'efile'],
#       dtype='object')
#
#
# # og
# Index(['annotation', 'enrich_per_mean_diff_by_genomstd', 'mean_trait_loci',
#        'median_trait_loci', 'emp_pvalue', 'mean_matched_across_trait_loci',
#        'std_matched_across_trait_loci', 'genome_wide_anno_std', 'matched_5th',
#        'median_matched_loci', 'matched_95th', 'region_summary',
#        'trait_summary', 'n_lead_loci', 'n_matched_loci',
#        'matched_snps_at_start', 'matched_snps_removed',
#        'percent_matched_removed', 'trait_snps_removed', 'trait', 'efile'],
#       dtype='object')
# %%

# multiple testing correction
from statsmodels.stats.multitest import multipletests

reject, pvals_corrected, _ , _  = multipletests( enrich_df['emp_pvalue'], method='bonferroni')
reject, bh_pvals_corrected, _ , _  = multipletests( enrich_df['emp_pvalue'], method='fdr_bh')
enrich_df['pval.adj_bh'] = bh_pvals_corrected
enrich_df['pval_label'] = enrich_df['emp_pvalue'].apply(lambda x: '*' if x < 0.05 else '')

# count n_lead_snps per trait
count_per_trait_df = enrich_df.groupby(['trait'])['n_lead_loci'].mean().reset_index()



count_per_trait_df['abr_trait'] = count_per_trait_df['trait']
count_per_trait_df.sort_values('n_lead_loci', inplace=True)
count_per_trait_df['n_lead_loci_round'] = count_per_trait_df['n_lead_loci'].apply(lambda x: f"{np.int(x):,}")
count_per_trait_df['trait(n)']= count_per_trait_df['abr_trait']+ " (" + count_per_trait_df['n_lead_loci_round'] + ")"
n_label = dict(zip(count_per_trait_df['trait'],count_per_trait_df['trait(n)']))
n_loci = dict(zip(count_per_trait_df['trait'],count_per_trait_df['n_lead_loci']))


nicer_labels = {'alz.IGAP_stage_1_2_combined': 'IGAP (19)',
                 'alz.GRACEStageI_dbGAP': 'GRACE (24)',
                 'alz.Kunkle_etal_Stage1_results': 'Kunkle_etal (62)',
                 "Lambert-2022":"Bellenguez (132)",
                 'Proxy_and_clinically_alz': 'Marioni (71)'}

# %%

annos = ['argweave', 'betascore',
        'linsigh', 'phastCon100', 'phyloP100',
         'fst_eas_afr', 'fst_eur_afr', 'fst_eur_eas',
        'xpehh_afr2_eas', 'xpehh_afr2_eur', 'xpehh_eas_eur']

axs = plot_bkg_and_trait_avg(annos, enrich_df, n_loci,nicer_labels, fig_height=2, fig_length_mult=1.2, lg_bbox=(3,0.4), ncol=1, xlabel_fontsize=8, markersize=7, lwd=1.4)
plt.tight_layout()

# plt.savefig(OUTPUT_DIR.joinpath(f'{DATE}_trait_vs_background_alzh_for_manu.pdf'))

# %%
#
# plt.style.use('/dors/capra_lab/users/abraha1/bin/global_scripts/matplotlib_settings/publication.mplstyle')
#
# trait2num_loci_dict = {'berg_paths': 2470.3333333333335,
#  'neale_paths': 3393.2,
#  'giant_paths': 5127.266666666666,
#  'bolt_paths': 6697.133333333333}
#
# trait2clean_label_dict={ 'berg_paths':  'Berg 2019',
#                          'neale_paths': 'Neale 2017',
#                          'giant_paths': 'GIANT 2018',
#                          'bolt_paths':  'Loh 2017'}
#
#
# axs = plot_bkg_and_trait_avg(annos, enrich_df, trait2num_loci_dict, nicer_labels, fig_height=2.,
#                              fig_length_mult=0.75, lg_bbox=(3,0.4), ncol=1, xlabel_fontsize=8, markersize=7, lwd=1.4)
# # %%
# sns.set(style="ticks",  font_scale=0.8, rc={"figure.figsize": (1.5*len(annos), 4)})
# fig, axs =plt.subplots(ncols=len(annos), sharey=True)
#
#
# for axind, this_anno in enumerate(annos):
#     keepcols=['enrich_per_mean_diff_by_genomstd', 'mean_trait_loci', 'median_trait_loci', 'mean_matched_across_trait_loci', 'matched_5th', 'matched_95th','median_matched_loci','trait',  'emp_pvalue','pval.adj_bh']
#     bg_df= enrich_df.loc[enrich_df['annotation']==this_anno,keepcols ].copy()
#
#
#     bg_df['n_loci'] = bg_df['trait'].map(n_loci)
#     bg_df.sort_values('n_loci', inplace=True)
#     bg_df.reset_index(drop=True, inplace=True)
#
#     ax = axs[axind]
#     # background
#     ax.scatter(y=bg_df.index, x=bg_df['mean_matched_across_trait_loci'], marker='.', color='gray')
#     for ind, row in bg_df.iterrows():
#         if ind == 0:
#             label='matched-background'
#         else:
#             label =''
#         ax.plot([row['matched_5th'],row['matched_95th']], [ind,ind], color='gray', label=label, zorder=-1)
#
#     # trait plots
#     ax.scatter(y=bg_df.loc[bg_df['pval.adj_bh']>=0.05,:].index,
#                 x=bg_df.loc[bg_df['pval.adj_bh']>=0.05,'mean_trait_loci'],
#                 marker='x', color='black', label='p.adj>=0.05', zorder=2)
#
#     ax.scatter(y=bg_df.loc[bg_df['pval.adj_bh']<0.05,:].index,
#                 x=bg_df.loc[bg_df['pval.adj_bh']<0.05,'mean_trait_loci'],
#                 marker='x', color='royalblue', label='p.adj<0.05', zorder=2)
#
#     bg_ax = ax.scatter(y=bg_df.loc[bg_df['pval.adj_bh']==0,:].index,
#                 x=bg_df.loc[bg_df['pval.adj_bh']==0,'mean_trait_loci'],
#                 marker='x', color='indianred', label=f'p.adj==0', zorder=2)
#
#     ax.axvline(mean_anno_dict[this_anno], linestyle=":", color='black', label='genome-wide mean')
#
#     ax.xaxis.set_major_locator(plt.MaxNLocator(3))
#     ax.set_yticks(bg_df.index)
#     ax.set_yticklabels(bg_df['trait'].map(n_label))
#     ax.set_xlabel(this_anno)
#     ax.tick_params(axis='y', which='major', length=0)
#     sns.despine(ax=ax, top=True, right=True, bottom= True, left=True)
#     if axind == (len(annos)-1):
#         ax.legend(loc='center left', bbox_to_anchor=(1.2, 0.85))
#
#     ax.minorticks_off()
#     # plt.minorticks_off()
# # plt.suptitle("Early vs. Late onset Vitiligo")
# plt.tight_layout()
# # plt.savefig(OUTPUT_DIR.joinpath(f'{DATE}_trait_vs_background_alzh.pdf'))
