#!/bin/python
# This script will ...
#
#
#
###
### TO REPRODUCE THE PLOT, load the file "2022-05-07_fig3a_bolt-lmm_bkg_trait_avg.tsv" in place of enrich_df on line 596; run all the other necessary varaible before
###
###
###
###
###
###
###

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
import textwrap
DATE = datetime.now().strftime('%Y-%m-%d')


# %load_ext autoreload
# %autoreload 2
# sys.path.append("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/scripts/0_common_plotting_scripts")
# from func_plot_bkg_and_trait_avg import plot_bkg_and_trait_avg
# from func_plot_radar import draw_radar

# plots
import proplot
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import seaborn as sns
%matplotlib inline

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

# import matplotlib.font_manager as font_manager
# fpath='/dors/capra_lab/users/abraha1/conda/envs/py36_r_ml/lib/python3.6/site-packages/matplotlib/mpl-data/fonts/ttf/Arial.ttf'
# font_dirs = ['/dors/capra_lab/users/abraha1/conda/envs/py36_r_ml/lib/python3.6/site-packages/matplotlib/mpl-data/fonts/ttf', ]
# font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
# font_list = font_manager.createFontList(font_files)
# font_manager.fontManager.ttflist.extend(font_list)
# mpl.rcParams['font.family'] = 'Arial'
# get_ipython().run_line_magic('config', "InlineBackend.figure_format='retina'")

from math import pi
import colorsys

### PATHS
DATA_FREEZE_DIR = Path("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/scripts/1_manuscript/2022_10_29_dataFreeze")

### ORIGINAL PATH TO FILES
### DO NOT USE IF ONLY TRYING TO RECREATE FIGURE
# ROOT= Path("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data")
# GSEL_DIR=ROOT.joinpath("scripts/run_mosaic_round_3_2021_01_27/2021_03_19_bolt-lmm_sstats/gsel_outputs")
# OUTPUT_DIR=ROOT.joinpath("scripts/1_manuscript/figs_blot_lmm")
# save_fig = lambda figname: plt.savefig(OUTPUT_DIR.joinpath(f"{DATE}_{figname}.pdf"))

# -----------
# FUNCTIONS
# -----------
import colorsys

def _get_colors(num_colors):
    colors=[]
    for i in np.arange(0., 360., 360. / num_colors):
        hue = i/360.
        lightness = (50 + np.random.rand() * 10)/100.
        saturation = (90 + np.random.rand() * 10)/100.
        colors.append(colorsys.hls_to_rgb(hue, lightness, saturation))
    return colors

def prep_radar_data(summary_df, pval_col):

    # required columsn
    #     enrich_per_mean_diff_by_genomstd
    #     annotation
    #     emp_pvalue -> This columsn shoudl be True or False

    plt_df = summary_df.copy()
    n_categories = plt_df.annotation.nunique()
    enrichs = summary_df['enrich_per_mean_diff_by_genomstd'].values.tolist()
    enrichs += enrichs[:1]

    pvals = list(plt_df[pval_col].values)
    angles = [n / float(n_categories) * 2 * pi for n in range(n_categories)]
    angles += angles[:1]

    sig_angles = [angles[i]  for i,x in enumerate(pvals) if x]
    sig_enrichs = [enrichs[i]  for i,x in enumerate(pvals) if x]

    return {'enrichs':enrichs, 'pvals':pvals, 'angles':angles, 'sig_angles':sig_angles, 'sig_enrichs':sig_enrichs}

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
                                            'fst_eas_afr': 'Fst\neas-afr',
                                            'fst_eur_afr': 'Fst\neur-afr',
                                            'fst_eur_eas': 'Fst\neur-eas',
                                            'xpehh_afr2_eas': 'XP-EHH\nafr-eas',
                                            'xpehh_afr2_eur': 'XP-EHH\nafr-eur',
                                            'xpehh_eas_eur': 'XP-EHH\neas-eur'}, fig_width=7.5,
                            fig_length_mult=1.5, fig_height=4, font_scale=0.8, xlabel_fontsize= 11, yticklabels_fontsize=9, legend_fontsize=12, lg_bbox=(-1.4, -0.35),ncol=5, markersize=5, lwd=2, ):

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
    # sns.set(style="ticks",  font_scale=font_scale, rc={"figure.figsize": (fig_length_mult*len(annos), fig_height)})
    sns.set(style="ticks",  font_scale=font_scale, rc={"figure.figsize": (fig_width, fig_height)})
    fig, axs =plt.subplots(ncols=len(annos), sharey=True)

    enrich_df.sort_values(['trait_category', 'n_lead_loci'], inplace=True)

    # # add a blank row in between traits
    # blank_row = pd.DataFrame.from_dict(enrich_df.head(1).to_dict(orient='index')[0], orient='index').T
    # blank_row['trait'] = 'blank'
    # trait2clean_label_dict['blank'] ='blank'
    #
    # new_df = pd.DataFrame()
    # for cat in enrich_df['trait_category'].unique():
    #
    #     this_df = enrich_df.query("trait_category == @cat").copy()
    #
    #     this_df = this_df.append(blank_row)
    #     new_df = new_df.append(this_df)
    #
    # enrich_df = new_df.copy()
    # print('hello')

    for axind, this_anno in enumerate(annos):


        keepcols=['mean_trait_loci', 'mean_matched_across_trait_loci','matched_5th','matched_95th','trait','pval.adj_bh']
        bg_df= enrich_df.loc[enrich_df['annotation']==this_anno,keepcols ].copy()

        # bg_df['n_loci'] = bg_df['trait'].map(trait2num_loci_dict)
        # bg_df.sort_values('n_loci', inplace=True)
        bg_df.reset_index(drop=True, inplace=True)

        bg_df['clean_label'] = bg_df['trait'].map(trait2clean_label_dict)
        bg_df.set_index('clean_label', inplace=True)
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
                    marker='*', color='gray', label='p.adj>=0.05', zorder=2, alpha=1, clip_on=False, s=markersize*1.1)

        ax.scatter(y=bg_df.loc[bg_df['pval.adj_bh']<0.05,:].index,
                    x=bg_df.loc[bg_df['pval.adj_bh']<0.05,'mean_trait_loci'],
                    marker='*', color='indianred', label='p.adj<0.05', zorder=2, clip_on=False, s=markersize*1.1)

        # x-axis
        ax.xaxis.set_major_locator(plt.LinearLocator(numticks=3))
        ax.xaxis.set_major_formatter(plt.FormatStrFormatter(anno_fmt_dict[this_anno]))
        for tick in ax.get_xticklabels():
            tick.set_rotation(270)
            tick.set_fontsize(yticklabels_fontsize)
        ax.set_xlabel(anno_label_dict[this_anno], fontsize=xlabel_fontsize)

        # yticks
        # ax.set_yticks(bg_df.index)
        # ax.set_yticklabels(bg_df['trait'].map(trait2clean_label_dict),fontsize=yticklabels_fontsize)
        # ax.tick_params(axis='x',
        ax.tick_params(axis='y', which='major', length=0, labelsize= yticklabels_fontsize)

        # spines
        sns.despine(ax=ax, top=True, right=True, bottom= False, left=True)
        if axind == (len(annos)-1):
            ax.legend(loc='center', bbox_to_anchor=lg_bbox,fancybox=False, shadow=False, ncol=ncol, fontsize=legend_fontsize)
        ax.minorticks_off()

    # plt.subplots_adjust(top=0.99, bottom=0.3, left=0.1, right=0.99, wspace=0.3)
    return axs

def plot_bkg_and_trait_avg1(annos, enrich_df, trait2num_loci_dict,trait2clean_label_dict,
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
                                            'xpehh_eas_eur': 'XP-EHH\neas-eur'}, fig_width=7.5,
                            fig_length_mult=1.5, fig_height=4, font_scale=0.8, xlabel_fontsize= 11, yticklabels_fontsize=9, legend_fontsize=12, lg_bbox=(-1.4, -0.35),ncol=5, markersize=5, lwd=2, ):

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
    # sns.set(style="ticks",  font_scale=font_scale, rc={"figure.figsize": (fig_length_mult*len(annos), fig_height)})
    sns.set(style="ticks",  font_scale=font_scale, rc={"figure.figsize": (fig_width, fig_height)})

    fig, axs =plt.subplots(ncols=len(annos),  sharey=True)
    # fig, axs =plt.subplots(ncols=len(annos), nrows=enrich_df['trait_category'].nunique(), sharey=True)
    enrich_df.sort_values(['trait_category', 'n_lead_loci'], inplace=True)

    new_df = pd.DataFrame()
    for cat in enrich_df['trait_category'].unique() :

        this_df = enrich_df.query('trait_category == @cat').copy()
        new_row= this_df.tail(1).copy()
        new_row['trait'] = 'blank+{}'.format(cat)
        trait2clean_label_dict['blank+{}'.format(cat)] = 'blank+{}'.format(cat)
        this_df = this_df.append(new_row)

        new_df = new_df.append(this_df)

    enrich_df = new_df.copy()

    for axind, this_anno in enumerate(annos):


        keepcols=['mean_trait_loci', 'mean_matched_across_trait_loci','matched_5th','matched_95th','trait','pval.adj_bh']
        bg_df= enrich_df.loc[enrich_df['annotation']==this_anno,keepcols ].copy()


        # add a blank ROW
        # new_row= bg_df.tail(1).copy()
        # new_row['trait'] = 'blank{}'.format()


        bg_df.reset_index(drop=True, inplace=True)

        bg_df['clean_label'] = bg_df['trait'].map(trait2clean_label_dict)
        bg_df.set_index('clean_label', inplace=True)
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
                    marker='*', color='gray', label='p.adj>=0.05', zorder=2, alpha=1, clip_on=False, s=markersize*1.1)

        ax.scatter(y=bg_df.loc[bg_df['pval.adj_bh']<0.05,:].index,
                    x=bg_df.loc[bg_df['pval.adj_bh']<0.05,'mean_trait_loci'],
                    marker='*', color='indianred', label='p.adj<0.05', zorder=2, clip_on=False, s=markersize*1.1)

        # x-axis
        ax.xaxis.set_major_locator(plt.LinearLocator(numticks=3))
        ax.xaxis.set_major_formatter(plt.FormatStrFormatter(anno_fmt_dict[this_anno]))
        for tick in ax.get_xticklabels():
            tick.set_rotation(270)
            tick.set_fontsize(yticklabels_fontsize)
        ax.set_xlabel(anno_label_dict[this_anno], fontsize=xlabel_fontsize)

        # yticks
        # ax.set_yticks(bg_df.index)
        # ax.set_yticklabels(bg_df['trait'].map(trait2clean_label_dict),fontsize=yticklabels_fontsize)
        # ax.tick_params(axis='x',
        ax.tick_params(axis='y', which='major', length=0, labelsize= yticklabels_fontsize)

        # spines
        sns.despine(ax=ax, top=True, right=True, bottom= False, left=True)
        if axind == (len(annos)-1):
            ax.legend(loc='center', bbox_to_anchor=lg_bbox,fancybox=False, shadow=False, ncol=ncol, fontsize=legend_fontsize)
        ax.minorticks_off()

    # plt.subplots_adjust(top=0.99, bottom=0.3, left=0.1, right=0.99, wspace=0.3)
    return axs

# %%
###
###    main
###
short_label_traits = {'disease_HI_CHOL_SELF_REP.sumstats':'high cholesterol',
        'disease_DERMATOLOGY.sumstats':'derm_disease',
        'blood_PLATELET_DISTRIB_WIDTH.sumstats':'platelet_dist_width',
        'body_BMIz.sumstats': 'bmi',
        'lung_FVCzSMOKE.sumstats': 'FVCzSmoke',
        'bp_DIASTOLICadjMEDz.sumstats': 'DBPadjMED',
        'pigment_HAIR_blonde.sumstats': 'hair_blond',
        'repro_NumberChildrenEverBorn_Pooled.sumstats': 'num_child_born',
        'bp_SYSTOLICadjMEDz.sumstats': 'SBP',
        'impedance_BASAL_METABOLIC_RATEz.sumstats': 'basal_metabolic_rate',
        'blood_PLATELET_COUNT.sumstats': 'platelet_counts',
        'disease_RESPIRATORY_ENT.sumstats': 'resp_disease',
        'blood_LYMPHOCYTE_COUNT.sumstats': 'lymphoctye_count',
        'blood_RED_COUNT.sumstats': 'red_cell_count',
        'blood_WHITE_COUNT.sumstats': 'white_cell_count',
        'body_WHRadjBMIz.sumstats': 'WHRadjBMI',
        'pigment_TANNING.sumstats': 'tanning',
        'pigment_SKIN.sumstats': 'skin_pigment',
        'other_MORNINGPERSON.sumstats': 'morning_person',
        'mental_NEUROTICISM.sumstats': 'neuroticism',
        'disease_HYPOTHYROIDISM_SELF_REP.sumstats': 'hypothyroidism_sr',
        'pigment_HAIR.sumstats': 'hair_pigment',
        'cov_EDU_YEARS.sumstats': 'edu-years',
        'blood_MEAN_PLATELET_VOL.sumstats': 'platelet_volume',
        'cov_EDU_COLLEGE.sumstats': 'edu-college',
        'disease_CARDIOVASCULAR.sumstats': 'cvd',
        'disease_AID_SURE.sumstats': 'AID_SURE',
        'disease_T2D.sumstats': 't2d',
        'body_BALDING4.sumstats': 'balding',
        'blood_MEAN_CORPUSCULAR_HEMOGLOBIN.sumstats': 'mean_corpus_hemoglobin',
        'disease_AID_ALL.sumstats': 'AID_all',
        'body_HEIGHTz.sumstats': 'height',
        'lung_FEV1FVCzSMOKE.sumstats': 'FEV1FVCzSMOKE',
        'body_BALDING1.sumstats': 'balding1',
        'repro_MENOPAUSE_AGE.sumstats': 'menopause_age',
        'cov_SMOKING_STATUS.sumstats': 'smoking_status',
        'bmd_HEEL_TSCOREz.sumstats': 'heel_tscorez',
        'repro_MENARCHE_AGE.sumstats': 'menarche_age',
        'disease_ALLERGY_ECZEMA_DIAGNOSED.sumstats': 'eczema',
        'disease_ASTHMA_DIAGNOSED.sumstats': 'asthma',
        'pigment_HAIR_darkbrown.sumstats': 'hair_darkbrown',
        'blood_RBC_DISTRIB_WIDTH.sumstats': 'ebc_dist_wdith',
        'pigment_SUNBURN.sumstats':'sunburn'}
nicer_label_traits =   {"blood_EOSINOPHIL_COUNT":'Eosinophil Count',
                        "blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT":'Reticulocyte Count',
                        "blood_LYMPHOCYTE_COUNT":'Lymphocyte Count',
                        "blood_MEAN_CORPUSCULAR_HEMOGLOBIN":'MCH',
                        "blood_MEAN_PLATELET_VOL":'Platelet Volume',
                        "blood_MEAN_SPHERED_CELL_VOL":'MSCV',
                        "blood_MONOCYTE_COUNT":'Monocyte Count',
                        "blood_PLATELET_COUNT":'Platelet Count',
                        "blood_PLATELET_DISTRIB_WIDTH":'Platelet Width',
                        "blood_RBC_DISTRIB_WIDTH":'RBC Width',
                        "blood_RED_COUNT":'RBC Count',
                        "blood_WHITE_COUNT":'WBC Count',
                        "bmd_HEEL_TSCOREz":'Heel T Score',
                        "body_BALDING1":'Balding(Type I)',
                        "body_BALDING4":'Balding(Type IV)',
                        "body_BMIz":'BMI',
                        "body_HEIGHTz":'Height',
                        "body_WHRadjBMIz":'Waist-hip Ratio',
                        "bp_DIASTOLICadjMEDz":'Diastolic BP',
                        "bp_SYSTOLICadjMEDz":'Systolic BP',
                        "cov_EDU_COLLEGE":'College Education',
                        "cov_EDU_YEARS":'Education(Yrs)',
                        "cov_SMOKING_STATUS":'Smoking',
                        "disease_AID_ALL":'Auto Immune',
                        "disease_AID_SURE":'Auto Immune (SURE)',
                        "disease_ALLERGY_ECZEMA_DIAGNOSED":'Eczema',
                        "disease_ASTHMA_DIAGNOSED":'Asthma',
                        "disease_CARDIOVASCULAR":'Heart Disease',
                        "disease_DERMATOLOGY":'Dermatologic',
                        "disease_HI_CHOL_SELF_REP":'High Cholesterol',
                        "disease_HYPOTHYROIDISM_SELF_REP":'Hypothyroidism',
                        "disease_RESPIRATORY_ENT":'Respiratory Diseases',
                        "disease_T2D":'Type 2 Diabetes',
                        "impedance_BASAL_METABOLIC_RATEz":'Basal Metabolic Rate',
                        "lung_FEV1FVCzSMOKE":'FEV1-FVC Ratio',
                        "lung_FVCzSMOKE":'FVC',
                        "mental_NEUROTICISM":'Neuroticism',
                        "other_MORNINGPERSON":'Morning Person',
                        "pigment_HAIR":'Hair Color',
                        "pigment_HAIR_blonde":'Blonde Hair',
                        "pigment_HAIR_darkbrown":'Dark Brown Hair',
                        "pigment_SKIN":'Skin Color',
                        "pigment_SUNBURN":'Sunburn',
                        "pigment_TANNING":'Tanning',
                        "repro_MENARCHE_AGE":'Menarche Age',
                        "repro_MENOPAUSE_AGE":'Menopause Age',
                        "repro_NumberChildrenEverBorn_Pooled":'Number Children'}
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
                'xpehh_eas_eur': 'XP-EHH\neas-eur'}
annos = ['argweave', 'betascore',
        'linsigh', 'phastCon100', 'phyloP100',
         'fst_eas_afr', 'fst_eur_afr', 'fst_eur_eas',
        'xpehh_afr2_eas', 'xpehh_afr2_eur', 'xpehh_eas_eur']

anno_col_list = _get_colors(5)
anno_color_dict = {'argweave':anno_col_list[0],
                     'betascore':anno_col_list[1],
                    'gerp':anno_col_list[2], 'linsigh':anno_col_list[2], 'phastCon100':anno_col_list[2], 'phyloP100':anno_col_list[2],
                    'fst_eas_afr':anno_col_list[3], 'fst_eur_afr':anno_col_list[3], 'fst_eur_eas':anno_col_list[3],
                    'xpehh_afr2_eas':anno_col_list[4], 'xpehh_afr2_eur':anno_col_list[4], 'xpehh_eas_eur':anno_col_list[4]}

# -----------
# check that all files ran
# -----------
#### DO NOT RUN THESE LINES IF ONLY TRYING TO RECREATE FIGURES
# enrich_files = list(GSEL_DIR.glob("*/*_extreme_regions_mean_enrichment_all_annotation.tsv"))
# bolt_traits = [x.name for x in list(GSEL_DIR.glob("*"))]

# len(enrich_files)
# -----------
# load enrichments
# -----------
# enrich_df = pd.DataFrame()
# for efile in enrich_files:
#     e_df = pd.read_csv( efile, sep="\t")
#     trait_name= efile.name.split('_extreme')[0]
#
#     if trait_name == "test_trait":
#         break
#     e_df['trait']  = trait_name
#     e_df['efile']  = efile.name
#     e_df.rename(columns={'n_lead_loci_final':'n_lead_loci','emp_pval':'emp_pvalue'},inplace=True)
#
#
#     keep_cols=['trait','efile','emp_pvalue', 'n_lead_loci', 'annotation', 'enrich_per_mean_diff_by_genomstd', "mean_matched_across_trait_loci","mean_trait_loci","matched_5th","matched_95th"]
#     enrich_df = enrich_df.append(e_df.loc[:, keep_cols])

# %%
#### DO NOT RUN THESE LINES IF ONLY TRYING TO RECREATE FIGURES
# trait_cat_df = pd.DataFrame.from_dict({x:x.split("_")[0] for x in enrich_df['efile'].unique()}, orient='index', columns=['category'])
# trait_cat_df['category'].replace('cov','social', inplace=True)
# trait_cat_df['category'].replace('other','social', inplace=True)

# %%
# -----------
# multiple testing correction
# -----------
#### DO NOT RUN THESE LINES IF ONLY TRYING TO RECREATE FIGURES
# from statsmodels.stats.multitest import multipletests
#
# ### multiptle testing correction over every trait and evo annotation
# reject, pvals_corrected, _ , _  = multipletests( enrich_df['emp_pvalue'], method='bonferroni')
# reject, bh_pvals_corrected, _ , _  = multipletests( enrich_df['emp_pvalue'], method='fdr_bh')
# enrich_df['pval.adj_bh'] = bh_pvals_corrected
# enrich_df['pval_label'] = enrich_df['emp_pvalue'].apply(lambda x: '*' if x < 0.05 else '')
#
# og_enrich_df = enrich_df.copy()
# # enrich_df=og_enrich_df.copy()
#
# # multiple test correct
# temp_df = pd.DataFrame()
# for anno in enrich_df['annotation'].unique():
#     anno_only_df = enrich_df.loc[enrich_df['annotation']==anno].copy()
#
#     reject, anno_pval_corr, _ , _ = multipletests(anno_only_df['emp_pvalue'], method='fdr_bh')
#     anno_only_df['pval.adj_bh_per_anno'] = anno_pval_corr
#     temp_df = temp_df.append(anno_only_df)
#
# enrich_df = temp_df.copy()
#
# # summarize
# count_per_trait_df = enrich_df.groupby(['trait'])['n_lead_loci'].mean().reset_index()
# count_per_trait_df['abr_trait'] = count_per_trait_df['trait'].apply(lambda x: nicer_label_traits[x.split('.sumstats')[0]])
# count_per_trait_df.sort_values('n_lead_loci', inplace=True)
# count_per_trait_df['n_lead_loci_round'] = count_per_trait_df['n_lead_loci'].apply(lambda x: f"{np.int(x):,}")
# count_per_trait_df['trait(n)']= count_per_trait_df['abr_trait']+ " (" + count_per_trait_df['n_lead_loci_round'] + ")"
# n_label = dict(zip(count_per_trait_df['trait'],count_per_trait_df['trait(n)']))
# n_loci = dict(zip(count_per_trait_df['trait'],count_per_trait_df['n_lead_loci']))


### SAVE n_label and n_loci
import pickle
# pickle.dump( n_label, open( DATA_FREEZE_DIR.joinpath("n_label.pickle"), "wb" ) )
# pickle.dump( n_loci, open( DATA_FREEZE_DIR.joinpath("n_loci.pickle"), "wb" ) )
n_label = pickle.load( open( DATA_FREEZE_DIR.joinpath("n_label.pickle"), "rb" ) )
n_loci = pickle.load( open( DATA_FREEZE_DIR.joinpath("n_loci.pickle"), "rb" ) )


# enrich_df = enrich_df.loc[~enrich_df['annotation'].isin(['B2','geva_allele_age','iES_Sabeti', 'gerp'])].copy()
# enrich_df['abr_trait'] = enrich_df['trait'].map(short_label_traits)
# enrich_df['nicer_trait_label'] = enrich_df['trait'].apply(lambda x: nicer_label_traits[x.split('.sumstats')[0]])
#
# # count of number of significant traits per annotation
# # n_traits n_enrich n_deplete n_tot prop
#
# sig_traits_df = pd.DataFrame()
# for anno_label, this_df in enrich_df.groupby('annotation'):
#
#     n_sig_bh_per_anno = sum(this_df['pval.adj_bh_per_anno'] < 0.05 )
#     n_sig_bh_all_corr = sum(this_df['pval.adj_bh'] < 0.05 )
#     n_tot_traits  = this_df.shape[0]
#     n_enrich = sum(this_df.loc[this_df['pval.adj_bh_per_anno'] < 0.05, 'enrich_per_mean_diff_by_genomstd'] >=0)
#     n_depelte = sum(this_df.loc[this_df['pval.adj_bh_per_anno'] < 0.05, 'enrich_per_mean_diff_by_genomstd'] < 0)
#
#     sig_traits_df = sig_traits_df.append(pd.DataFrame({'anno':[anno_label], "n_sig_bh_per_anno": [n_sig_bh_per_anno],"n_sig_bh_all_corr": [n_sig_bh_all_corr],"n_tot_traits": [n_tot_traits],"n_enrich": [n_enrich],"n_depelte": [n_depelte]}))
#
# enrich_df.query('annotation=="xpehh_afr2_eas"')
# enrich_df.loc[(enrich_df['annotation']=="xpehh_afr2_eas") & (enrich_df['pval.adj_bh'] < 0.05 )]
#
#
#
# new_col_dict = {'anno': 'annotation',
#  'n_sig_bh_per_anno': 'n_traits_sig_after_bh_over_annotations',
#  'n_sig_bh_all_corr': 'n_traits_sig_after_bh_over_annotations_and_traits',
#  'n_tot_traits': 'n_tot_traits',
#  'n_enrich':'n_traits_enrich',
#  'n_depelte':'n_traits_depelte'}
#
# sig_traits_df.rename(columns=new_col_dict, inplace=True)
# # sig_traits_df.to_csv(OUTPUT_DIR.joinpath(f"{DATE}_count_traits_after_bh_correction.tsv"), sep="\t", index=False)


# %%
# write
# nice_colums = ['nicer_trait_label','annotation', 'n_lead_loci', 'enrich_per_mean_diff_by_genomstd', 'emp_pvalue', 'mean_trait_loci', 'pval.adj_bh']
# enrich_df.loc[:, nice_colums].to_csv(OUTPUT_DIR.joinpath(f"{DATE}_enrichment_boltlmm.tsv"), sep="\t", index=False)
# %%
# enrich_df.groupby(['trait'])['pval_label'].apply(lambda x: sum(x=="*")).reset_index()
# enrich_df.groupby('annotation').apply(lambda x: sum(x['pval_label'] =="*"))

# %%
# -----------
# CLUSTER
# -----------
# wide_enrich_df = enrich_df.pivot(index='annotation', columns='trait', values='enrich_per_mean_diff_by_genomstd')
# wide_padj_df = enrich_df.pivot(index='annotation', columns='trait', values='pval_label')
# wide_enrich_df = wide_enrich_df.reindex(annos)
# wide_enrich_df.columns = [nicer_label_traits[x.split('.sumstats')[0]] for x in wide_enrich_df.columns]
# wide_enrich_df.index = [anno_label_dict[x] for x in wide_enrich_df.index]

# %%
# from scipy.cluster.hierarchy import fclusterdata
# clust_ = fclusterdata(wide_enrich_df.T, t=0.75, metric='euclidean',method='average', criterion='distance')
# trait_clust_df = pd.DataFrame({'trait': list(wide_enrich_df.T.index), 'cluster': clust_})
# cols =dict(zip(trait_clust_df['cluster'].unique(), _get_colors(trait_clust_df['cluster'].nunique())))
# trait_clust_df['color'] = trait_clust_df['cluster'].map(cols)
# trait_clust_df.set_index('trait',inplace=True)
# trait_clust_df['cluster'].nunique()
# trait_clust_df.reset_index(inplace=True)
#
# trait2clust = dict(zip(trait_clust_df['trait'],trait_clust_df['cluster']))



# %%
# -----------
# plot bkg adn trait evo value
# -----------


# add categories
# enrich_df['trait_category'] = enrich_df['efile'].map(trait_cat_df.to_dict()['category'])
#
# # add hclust ids
# annos = ['argweave', 'betascore',
#         'linsigh', 'phastCon100', 'phyloP100',
#          'fst_eas_afr', 'fst_eur_afr', 'fst_eur_eas',
#         'xpehh_afr2_eas', 'xpehh_afr2_eur', 'xpehh_eas_eur']
# enrich_df['clust'] = enrich_df['nicer_trait_label'].map(trait2clust)
# enrich_df['n_loci'] = enrich_df['trait'].map(n_loci)
# enrich_df.sort_values(['clust',"n_loci"], inplace=True)

# %%
# subset
# enrich_df['annotation'].unique()
keep_anno = ['argweave', 'betascore', 'fst_eur_afr', 'linsigh', 'phyloP100', 'xpehh_afr2_eur']

##### NOTE TO MAKE SUP FIGURE S3!!! ******
# to make Supplemental figure S3, repalce with all of the annotations
# keep_anno = ['argweave', 'betascore', 'fst_eur_afr', 'linsigh', 'phyloP100', 'xpehh_afr2_eur']
# %%
# save to freeze
# enrich_df.to_csv(DATA_FREEZE_DIR.joinpath(f"{DATE}_fig3a_bolt-lmm_bkg_trait_avg.tsv"), sep="\t", index=False)

enrich_df = pd.read_csv(DATA_FREEZE_DIR.joinpath(f"2022-05-07_fig3a_bolt-lmm_bkg_trait_avg.tsv.gz"), sep="\t",)
# %%

# axs = plot_bkg_and_trait_avg1(keep_anno, enrich_df, n_loci,n_label, fig_height=6.5, fig_width=6, font_scale=1.0, xlabel_fontsize=8, yticklabels_fontsize=8,
                             # legend_fontsize=8,lg_bbox=(-14, -.05), ncol=1,markersize=30, lwd=1)
axs = plot_bkg_and_trait_avg(keep_anno, enrich_df, n_loci,n_label, fig_height=6.5, fig_width=6, font_scale=1.0, xlabel_fontsize=8, yticklabels_fontsize=8,
                             legend_fontsize=8,lg_bbox=(-14, -.05), ncol=1,markersize=30, lwd=1)
plt.subplots_adjust(top=0.99, bottom=0.13, left=0.15, right=0.99, wspace=0.3)


# save_fig('bkg_bolt_lmm_sub_anno')
# %%
