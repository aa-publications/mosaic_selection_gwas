#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2021-05-23 18:46:35


import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime


import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
DATE = datetime.now().strftime('%Y-%m-%d')


def plot_bkg_and_trait_avg(annos, enrich_df, trait2num_loci_dict,trait2clean_label_dict,
                            anno_fmt_dict={'argweave':"%.0e", 'betascore':"%.1f",
                                            'linsigh':"%.2f", 'phastCon100':"%.2f", 'phyloP100':"%.2f", 'gerp':"%.2f",
                                            'fst_eas_afr':"%.2f", 'fst_eur_afr':"%.2f", 'fst_eur_eas':"%.2f",
                                            'xpehh_afr2_eas':"%.2f", 'xpehh_afr2_eur':"%.2f", 'xpehh_eas_eur':"%.2f"},
                            anno_label_dict={'argweave': 'ARGweaver\n(TMRCA)','betascore': 'Beta Score',
                                            'linsigh': 'LINSIGHT',
                                            'phastCon100': 'PhastCons',
                                            'phyloP100': 'PhyloP',
                                            'gerp': 'GERP',
                                            'fst_eas_afr': 'F_ST\neas-afr',
                                            'fst_eur_afr': 'F_ST\neur-afr',
                                            'fst_eur_eas': 'F_ST\neur-eas',
                                            'xpehh_afr2_eas': 'XP-EHH\nafr-eas',
                                            'xpehh_afr2_eur': 'XP-EHH\nafr-eur',
                                            'xpehh_eas_eur': 'XP-EHH\neas-eur'},
                            fig_length_mult=1.5, fig_height=4, font_scale=0.8, xlabel_fontsize= 11, yticklabels_fontsize=9, legend_fontsize=12, lg_bbox=(-1.4, -0.35),ncol=5, markersize=5, lwd=2):

    """
        * enrich_df: columns
           - annotation: label of annotation
           - trait: name of trait
           - mean_matched_across_trait_loci: mean value across matched
           - mean_trait_loci: mean value across trait loci
           - matched_5th
           - matched_95th
           - pval.adj_bh
        * trait2num_loci_dict
        * trait2clean_label_dict
        * n_loci: 'trait' to 'n_loci dictionary
        * anno_fmt_dict: 'anno' to  fmt (how to represent values of this annotaiton)
        * anno_label_dict: anno to formmatted label

    """



    enrich_df.reset_index(drop=True, inplace=True)
    sns.set(style="ticks",  font_scale=font_scale, rc={"figure.figsize": (fig_length_mult*len(annos), fig_height)})
    fig, axs =plt.subplots(ncols=len(annos), sharey=True)


    for axind, this_anno in enumerate(annos):


        keepcols=['mean_trait_loci', 'mean_matched_across_trait_loci','matched_5th','matched_95th','trait','pval.adj_bh']
        bg_df= enrich_df.loc[enrich_df['annotation']==this_anno,keepcols ].copy()

        bg_df['n_loci'] = bg_df['trait'].map(trait2num_loci_dict)
        bg_df.sort_values('n_loci', inplace=True)
        bg_df.reset_index(drop=True, inplace=True)

        ax = axs[axind]

        # background: mean + 5th and 95th
        ax.scatter(y=bg_df.index, x=bg_df['mean_matched_across_trait_loci'], marker='.', color='gray', s=markersize)
        for ind, row in bg_df.iterrows():
            if ind == 0:
                label='matched-background'
            else:
                label =''
            ax.plot([row['matched_5th'],row['matched_95th']], [ind,ind],  linewidth=lwd ,color='gray', label=label, zorder=-1)

        # trait plots
        ax.scatter(y=bg_df.loc[bg_df['pval.adj_bh']>=0.05,:].index,
                    x=bg_df.loc[bg_df['pval.adj_bh']>=0.05,'mean_trait_loci'],
                    marker='*', color='black', label='p.adj>=0.05', zorder=2, alpha=0.6, clip_on=False, s=markersize*1.5)

        ax.scatter(y=bg_df.loc[bg_df['pval.adj_bh']<0.05,:].index,
                    x=bg_df.loc[bg_df['pval.adj_bh']<0.05,'mean_trait_loci'],
                    marker='*', color='indianred', label='p.adj<0.05', zorder=2, clip_on=False, s=markersize*1.5)

        # x-axis
        ax.xaxis.set_major_locator(plt.LinearLocator(numticks=3))
        ax.xaxis.set_major_formatter(plt.FormatStrFormatter(anno_fmt_dict[this_anno]))
        for tick in ax.get_xticklabels():
            tick.set_rotation(270)
        ax.set_xlabel(anno_label_dict[this_anno], fontsize=xlabel_fontsize)

        # yticks
        ax.set_yticks(bg_df.index)
        ax.set_yticklabels(bg_df['trait'].map(trait2clean_label_dict),fontsize=yticklabels_fontsize)
        ax.tick_params(axis='y', which='major', length=0)

        # spines
        sns.despine(ax=ax, top=True, right=True, bottom= False, left=True)
        if axind == (len(annos)-1):
            ax.legend(loc='center', bbox_to_anchor=lg_bbox,fancybox=False, shadow=False, ncol=ncol, fontsize=legend_fontsize)
        ax.minorticks_off()

    # plt.subplots_adjust(top=0.99, bottom=0.3, left=0.1, right=0.99, wspace=0.3)
    return axs

# load enrichment for all traits
def calc_p_adj(og_enrich_df):

    # columsn required:
    #       * emp_pvalue
    # will calculate fdr_bh, across all pvalue
    from statsmodels.stats.multitest import multipletests

    enrich_df = og_enrich_df.copy()
    _, bh_pvals_corrected, _ , _  = multipletests( enrich_df['emp_pvalue'], method='fdr_bh')
    enrich_df['pval.adj_bh'] = bh_pvals_corrected
    enrich_df['pval_nom_sig'] = enrich_df['emp_pvalue'].apply(lambda x: '*' if x < 0.05 else '')

    return enrich_df

def make_n_loci_and_label(og_enrich_df):
    # need column 'trait' and 'n_lead_loci'
    # will make dictionary with each trait (key) and value as label or n-lead-loci

    enrich_df = og_enrich_df.copy()

    # count n_lead_snps per trait
    count_per_trait_df = enrich_df.groupby(['trait']).apply(lambda x: x['n_lead_loci'].mean()).reset_index()
    count_per_trait_df['abr_trait'] = count_per_trait_df['trait']
    count_per_trait_df.rename({0:'n_lead_loci'}, axis=1, inplace=True)
    count_per_trait_df.sort_values('n_lead_loci', inplace=True)
    count_per_trait_df['n_lead_loci_round'] = count_per_trait_df['n_lead_loci'].apply(lambda x: f"{np.int(x):,}")
    count_per_trait_df['trait(n)']= count_per_trait_df['abr_trait']+ " (" + count_per_trait_df['n_lead_loci_round'] + ")"
    trait2clean_label_dict = dict(zip(count_per_trait_df['trait'],count_per_trait_df['trait(n)']))
    trait2num_loci_dict = dict(zip(count_per_trait_df['trait'],count_per_trait_df['n_lead_loci']))

    return trait2num_loci_dict, trait2clean_label_dict
