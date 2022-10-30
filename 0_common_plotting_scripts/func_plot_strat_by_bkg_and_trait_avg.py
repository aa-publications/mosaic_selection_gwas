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
from datetime import datetime

from pathlib import Path


from collections import OrderedDict
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
DATE = datetime.now().strftime('%Y-%m-%d')

# %config InlineBackend.figure_format='retina'

DATE = datetime.now().strftime('%Y-%m-%d')



def get_mean_anno_dict():

    ANNO_FILE=Path("/dors/capra_lab/projects/gwas_allele_age_evolution/scripts/pipeline/dev/gsel_vec/gsel_vec/create_annotations/genome_wide_summary_of_annotations.tsv")
    anno_summary_df = pd.read_csv( ANNO_FILE, sep="\t")
    mean_anno_dict = dict(zip(anno_summary_df['annotation'], anno_summary_df['mean']))

    return mean_anno_dict

def make_clean_labels(exp_df, type):


    if type == "indep":
        float_key_dict = {np.float(x.split(",")[0][1:]):'({:.1E}, {:.1E}]'.format(np.float(x.split(',')[0][1:]), np.float(x.split(',')[1][:-1])) for x in exp_df['partition_id'].unique()}

        orderded_labels = OrderedDict()
        for ii, key in enumerate(sorted(float_key_dict.keys())):
            orderded_labels[ii]= float_key_dict[key]



        temp_float_raw_dict = {np.float(x.split(",")[0][1:]):x for x in exp_df['partition_id'].unique()}
        temp_orderded_labels = OrderedDict()
        for ii, key in enumerate(sorted(temp_float_raw_dict.keys())):
            temp_orderded_labels[ii]= temp_float_raw_dict[key]

        label2_number_dict = {v:k for k,v in temp_orderded_labels.items()}

    elif type == "cumulative":
        # make number to label dictionary
        float_number = [exp_df['partition_id'].unique()[0].split(",")[0][1:]] + [x.split(",")[-1][:-1] for i,x in enumerate(exp_df['partition_id'].unique())]
        float_label = ["({:.1E} {:.1E}]".format( np.float(x.split(",")[0][1:]), np.float(x.split(",")[-1][:-1])) for x in exp_df['partition_id'].unique()]
        float_key_dict = dict(zip(float_number, float_label))
        orderded_labels = OrderedDict()
        for ii, key in enumerate(sorted(float_key_dict.keys())):
            orderded_labels[ii]= float_key_dict[key]


        # make parittion_id to ranked number
        temp_float_raw_dict = {np.float(x.split("_")[-1].split(",")[0][1:]):x for x in exp_df['partition_id'].unique()}
        temp_orderded_labels = OrderedDict()
        for ii, key in enumerate(sorted(temp_float_raw_dict.keys())):
            temp_orderded_labels[ii]= temp_float_raw_dict[key]

        label2_number_dict = {v:k for k,v in temp_orderded_labels.items()}

    elif type =="shuffled_indep":

        orderded_labels = dict(zip(shuf_indep_df['partition_id'], shuf_indep_df['partition_id']))

        label2_number_dict = dict(zip(shuf_indep_df['partition_id'], shuf_indep_df['partition_id']))

    elif type =="shuffled_cumulative":
        shuf_cumulat_df['partition_id'].unique()

        float_number = [int(x.split("_")[-1]) for x in shuf_cumulat_df['partition_id'].unique()]
        float_label=["{}to{}".format(x.split("_")[0],x.split("_")[-1])  for x in shuf_cumulat_df['partition_id'].unique()]
        float_key_dict = dict(zip(float_number, float_label))
        orderded_labels = OrderedDict()
        for ii, key in enumerate(sorted(float_key_dict.keys())):
            orderded_labels[ii]= float_key_dict[key]


        temp_float_raw_dict = {int(x.split("_")[-1]):x for x in shuf_cumulat_df['partition_id'].unique()}
        temp_orderded_labels = OrderedDict()
        for ii, key in enumerate(sorted(temp_float_raw_dict.keys())):
            temp_orderded_labels[ii]= temp_float_raw_dict[key]

        label2_number_dict = {v:k for k,v in temp_orderded_labels.items()}


    return orderded_labels, label2_number_dict

def plot_stratified_enrichments_clean(annotations, exp1_df, plt_title="Stratified Enrichments", savefile=None, w_anno=False,
                                      anno_fmt_dict={'argweave':"%.0e", 'betascore':"%.1f",
                                                    'linsigh':"%.2f", 'phastCon100':"%.2f", 'phyloP100':"%.2f", 'gerp':"%.2f",
                                                    'fst_eas_afr':"%.2f", 'fst_eur_afr':"%.2f", 'fst_eur_eas':"%.2f",
                                                    'xpehh_afr2_eas':"%.2f", 'xpehh_afr2_eur':"%.2f", 'xpehh_eas_eur':"%.2f"},
                                    anno_label_dict={'argweave': 'ARGweaver\n(TMRCA)','betascore': 'Beta Score',
                                                    'linsigh': 'LINSIGHT',
                                                    'phastCon100': 'PhastCons',
                                                    'phyloP100': 'PhyloP',
                                                    'gerp': 'GERP',
                                                    'fst_eas_afr': 'Fst\neas-afr',
                                                    'fst_eur_afr': 'Fst\neur-afr',
                                                    'fst_eur_eas': 'Fst\neur-eas',
                                                    'xpehh_afr2_eas': 'XP-EHH\nafr-eas',
                                                    'xpehh_afr2_eur': 'XP-EHH\nafr-eur',
                                                    'xpehh_eas_eur': 'XP-EHH\neas-eur'},
                                    xlabel_fontsize=10, legend_fontsize=10, fig_length_mult=1.5, fig_height=4, font_scale=0.8, yticklabel_fontsize=10, lg_bbox=(-1.4, -0.35),ncol=5):

    print("Assumes using non-overlapping independent p-value partitions....")
    ordered_labels, label2_number_dict = make_clean_labels(exp1_df, 'indep')
    exp1_df['partition_bin'] = exp1_df['partition_id'].map(label2_number_dict)


    nrows=1
    ncols=int(np.ceil(len(annotations)/nrows))
    sns.set(style="ticks",  font_scale=font_scale, rc={"figure.figsize": (fig_length_mult*len(annotations), fig_height)})
    fig, axs =plt.subplots(ncols=ncols, nrows=nrows, sharey=True)
    axs = axs.ravel()
    for axind, this_anno in enumerate(annotations):

        bg_df= exp1_df.loc[exp1_df['annotation']==this_anno,: ].copy()

        ax = axs[axind]

        # background
        ax.scatter(y=bg_df['partition_bin'], x=bg_df['mean_matched_across_trait_loci'], marker='.', color='gray')
        for ind, (_, row) in enumerate(bg_df.iterrows()):
            if ind ==( bg_df.shape[0]-1):
                label='matched-background'
            else:
                label =''
            ax.plot([row['matched_5th'],row['matched_95th']], [row['partition_bin'],row['partition_bin']], color='gray', label=label)
            ax.annotate(f'n_snps={row["n_lead_loci_final"]}', xy=(row['mean_matched_across_trait_loci'], row['partition_bin']+0.1), ha='center', color='gray', fontsize=8) if w_anno else None

        # trait plots
        ax.scatter(y=bg_df.loc[bg_df['emp_pval']>=0.05,'partition_bin'],
                    x=bg_df.loc[bg_df['emp_pval']>=0.05,'mean_trait_loci'],
                    marker='x', color='black', label='p.adj>=0.05', zorder=2, alpha=0.6)

        ax.scatter(y=bg_df.loc[bg_df['emp_pval']<0.05,'partition_bin'],
                    x=bg_df.loc[bg_df['emp_pval']<0.05,'mean_trait_loci'],
                    marker='*', color='indianred', label='p.adj<0.05',  zorder=2)


        # x-axis
        ax.xaxis.set_major_locator(plt.LinearLocator(numticks=3))
        ax.xaxis.set_major_formatter(plt.FormatStrFormatter(anno_fmt_dict[this_anno]))
        for tick in ax.get_xticklabels():
            tick.set_rotation(270)
        ax.set_xlabel(anno_label_dict[this_anno], fontsize=xlabel_fontsize)


        # y-axis
        ax.set_yticks(list(ordered_labels.keys()))
        ax.set_yticklabels([v for k,v in ordered_labels.items()], fontsize=yticklabel_fontsize)
        ax.set_ylabel("GWAS P-Value") if (axind == 0) else None


        if axind == (len(annotations)-1):
            ax.legend(loc='center', bbox_to_anchor=lg_bbox,fancybox=False, shadow=False, ncol=ncol, fontsize=legend_fontsize)

            # ax.legend(loc='center left', bbox_to_anchor=(1.2, 0.85))
            # ax.legend(loc='upper center', bbox_to_anchor=(-8, -0.15), fancybox=False, shadow=False, ncol=5)
        sns.despine(ax=ax, top=True,left=True, right=True)
        ax.tick_params(axis='y', which='major', length=0)
        ax.minorticks_off()

    # plt.suptitle(plt_title)
    # plt.tight_layout()
    if savefile:
        plt.savefig(savefile)

    return axs