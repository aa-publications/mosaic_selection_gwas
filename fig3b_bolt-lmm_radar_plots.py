#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2021-03-17 14:51:44

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
sys.path.append("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/scripts/0_common_plotting_scripts")
from func_plot_bkg_and_trait_avg import plot_bkg_and_trait_avg
from func_plot_radar import draw_radar

# plots
# import proplot
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import seaborn as sns
%matplotlib inline

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

import matplotlib.font_manager as font_manager
fpath='/dors/capra_lab/users/abraha1/conda/envs/py36_r_ml/lib/python3.6/site-packages/matplotlib/mpl-data/fonts/ttf/Arial.ttf'
font_dirs = ['/dors/capra_lab/users/abraha1/conda/envs/py36_r_ml/lib/python3.6/site-packages/matplotlib/mpl-data/fonts/ttf', ]
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
font_list = font_manager.createFontList(font_files)
font_manager.fontManager.ttflist.extend(font_list)
mpl.rcParams['font.family'] = 'Arial'
get_ipython().run_line_magic('config', "InlineBackend.figure_format='retina'")

from math import pi
import colorsys

### PATHS
ROOT= Path("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data")
GSEL_DIR=ROOT.joinpath("scripts/run_mosaic_round_3_2021_01_27/2021_03_19_bolt-lmm_sstats/gsel_outputs")
OUTPUT_DIR=ROOT.joinpath("scripts/1_manuscript/figs_blot_lmm")
OUTPUT_PER_TRAIT_DIR=ROOT.joinpath("scripts/1_manuscript/figs_blot_lmm/per_trait")
save_fig = lambda figname: plt.savefig(OUTPUT_DIR.joinpath(f"{DATE}_{figname}.pdf"))
save_per_trait_fig = lambda figname: plt.savefig(OUTPUT_PER_TRAIT_DIR.joinpath(f"{DATE}_{figname}.pdf"))
DATA_FREEZE_DIR = Path("/dors/capra_lab/projects/gwas_allele_age_evolution/mosaic_selection/abin_data/scripts/1_manuscript/2022_10_29_dataFreeze")
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


# %%
###
###    main
###

# -----------
# setting up variables for making these plots
# -----------
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
                        "blood_MEAN_CORPUSCULAR_HEMOGLOBIN":'Mean Corpular Hemoglobin',
                        "blood_MEAN_PLATELET_VOL":'Mean Platelet Volume',
                        "blood_MEAN_SPHERED_CELL_VOL":'Mean Sphered Cell Volume',
                        "blood_MONOCYTE_COUNT":'Monocyte Count',
                        "blood_PLATELET_COUNT":'Platelet Count',
                        "blood_PLATELET_DISTRIB_WIDTH":'Platelet Distribution Width',
                        "blood_RBC_DISTRIB_WIDTH":'Red Blood Cell Distribution Width',
                        "blood_RED_COUNT":'Red Blood Cell Count',
                        "blood_WHITE_COUNT":'White Blood Cell Count',
                        "bmd_HEEL_TSCOREz":'Heel T Score',
                        "body_BALDING1":'Balding Type I',
                        "body_BALDING4":'Balding Type IV',
                        "body_BMIz":'BMI',
                        "body_HEIGHTz":'Height',
                        "body_WHRadjBMIz":'Waist-hip Ratio',
                        "bp_DIASTOLICadjMEDz":'Diastolic Blood Pressure',
                        "bp_SYSTOLICadjMEDz":'Systolic Blood Pressure',
                        "cov_EDU_COLLEGE":'College Education',
                        "cov_EDU_YEARS":'Years of Education',
                        "cov_SMOKING_STATUS":'Smoking Status',
                        "disease_AID_ALL":'Auto Immune Traits',
                        "disease_AID_SURE":'Auto Immune Traits (Sure)',
                        "disease_ALLERGY_ECZEMA_DIAGNOSED":'Eczema',
                        "disease_ASTHMA_DIAGNOSED":'Asthma',
                        "disease_CARDIOVASCULAR":'Cardiovascular Diseases',
                        "disease_DERMATOLOGY":'Dermatologic Diseases',
                        "disease_HI_CHOL_SELF_REP":'High Cholesterol',
                        "disease_HYPOTHYROIDISM_SELF_REP":'Hypothyroidism',
                        "disease_RESPIRATORY_ENT":'Respiratory and Ear-nose-throat Diseases',
                        "disease_T2D":'Type 2 Diabetes',
                        "impedance_BASAL_METABOLIC_RATEz":'Basal Metabolic Rate',
                        "lung_FEV1FVCzSMOKE":'FEV1-FVC Ratio',
                        "lung_FVCzSMOKE":'Forced Vital Capacity (FVC)',
                        "mental_NEUROTICISM":'Neuroticism',
                        "other_MORNINGPERSON":'Morning Person',
                        "pigment_HAIR":'Hair Color',
                        "pigment_HAIR_blonde":'Blonde Hair',
                        "pigment_HAIR_darkbrown":'Dark Brown Hair',
                        "pigment_SKIN":'Skin Color',
                        "pigment_SUNBURN":'Sunburn Occasion',
                        "pigment_TANNING":'Tanning',
                        "repro_MENARCHE_AGE":'Age at Menarche',
                        "repro_MENOPAUSE_AGE":'Age at Menopause',
                        "repro_NumberChildrenEverBorn_Pooled":'Number Children (Pooled)'}
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
### SKIP THIS IF TRYING TO REPRODUCE FIGURES
### GSEL_DIR is too large of a directory to share easily
# enrich_files = list(GSEL_DIR.glob("*/*_extreme_regions_mean_enrichment_all_annotation.tsv"))
# bolt_traits = [x.name for x in list(GSEL_DIR.glob("*"))]


# -----------
# load enrichments
# -----------
### SKIP THIS IF TRYING TO REPRODUCE FIGURES
### GSEL_DIR is too large of a directory to share easily

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
#
#
# trait_cat_df = pd.DataFrame.from_dict({x:x.split("_")[0] for x in enrich_df['efile'].unique()}, orient='index', columns=['category'])
#
# trait_cat_df['category'].replace('cov','social', inplace=True)
# trait_cat_df['category'].replace('other','social', inplace=True)
#
#
#
#
# enrich_df = enrich_df.loc[enrich_df['annotation'].isin(annos), :].copy()
# enrich_df['annotation'].unique()





# %%
# -----------
# multiple testing correction
# -----------
### SKIP THIS IF TRYING TO REPRODUCE FIGURES
### GSEL_DIR is too large of a directory to share easily
#
# from statsmodels.stats.multitest import multipletests
#
# reject, pvals_corrected, _ , _  = multipletests( enrich_df['emp_pvalue'], method='bonferroni')
# reject, bh_pvals_corrected, _ , _  = multipletests( enrich_df['emp_pvalue'], method='fdr_bh')
# enrich_df['pval.adj_bh'] = bh_pvals_corrected
# enrich_df['pval_label'] = enrich_df['emp_pvalue'].apply(lambda x: '*' if x < 0.05 else '')
#
# count_per_trait_df = enrich_df.groupby(['trait'])['n_lead_loci'].mean().reset_index()
# count_per_trait_df['abr_trait'] = count_per_trait_df['trait'].apply(lambda x: nicer_label_traits[x.split('.sumstats')[0]])
# count_per_trait_df.sort_values('n_lead_loci', inplace=True)
# count_per_trait_df['n_lead_loci_round'] = count_per_trait_df['n_lead_loci'].apply(lambda x: f"{np.int(x):,}")
# count_per_trait_df['trait(n)']= count_per_trait_df['abr_trait']+ " (" + count_per_trait_df['n_lead_loci_round'] + ")"
# n_label = dict(zip(count_per_trait_df['trait'],count_per_trait_df['trait(n)']))
# n_loci = dict(zip(count_per_trait_df['trait'],count_per_trait_df['n_lead_loci']))
#
# enrich_df = enrich_df.loc[~enrich_df['annotation'].isin(['B2','geva_allele_age','iES_Sabeti', 'ger[]'])].copy()
# enrich_df['abr_trait'] = enrich_df['trait'].map(short_label_traits)
# enrich_df['nicer_trait_labe'] = enrich_df['trait'].apply(lambda x: nicer_label_traits[x.split('.sumstats')[0]])
# enrich_df['trait_category'] = enrich_df['efile'].map(trait_cat_df.to_dict()['category'])

# enrich_df.to_csv(DATA_FREEZE_DIR.joinpath(f"{DATE}_boltlmm-radar_plots.tsv"), sep="\t", index=False)
# %%
# LOAD HERE IF TRYING TO REPRODUCE FIGURES
erich_df = pd.read_csv( DATA_FREEZE_DIR.joinpath(f"2022-05-07_fig3b_boltlmm-radar_plots.tsv.gz"), sep="\t")
# %%
# -----------
# radar plots of specitic traits
# -----------
pval_col='psig'
enrich_df[pval_col]= enrich_df['pval.adj_bh'] < 0.05
anno_color_dict = {'argweave':'salmon',
                   'betascore':'olivedrab',
                   'fst_eas_afr':'plum', 'fst_eur_afr':'plum','fst_eur_eas':'plum',
                   'gerp':'seagreen', 'linsigh':'seagreen', 'phastCon100':'seagreen', 'phyloP100':'seagreen',
                   'xpehh_afr2_eas':'royalblue', 'xpehh_afr2_eur':'royalblue', 'xpehh_eas_eur':'royalblue'}
fig_ht=3

for this_trait in enrich_df['trait'].unique():
    this_trait_df = enrich_df.query('trait==@this_trait').copy()

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(fig_ht*1,fig_ht), subplot_kw=dict(polar=True))



    trait_color=dict(zip(this_trait_df['trait'].unique(), ['k']*this_trait_df['trait'].nunique()  ))
    TRAIT_LABEL=""

    max_y= 1
    min_y= -1
    # max_y= np.ceil(this_trait_df['enrich_per_mean_diff_by_genomstd'].max())+1
    # min_y= np.floor(this_trait_df['enrich_per_mean_diff_by_genomstd'].min())-1


    ax, labels = draw_radar(this_trait_df, pval_col, ax, TRAIT_LABEL,
                            anno_color_dict=anno_color_dict,
                            radar_fill='gray',
                            trait_color_dict=trait_color, lwd=1.5,
                            plt_params={'marker_size':2.5,'fill_alpha':0.2,'enrich_label_fontsize':7, 'enrich_label_size':7})

    for x in labels:
        x.set_text(anno_label_dict[x.get_text()])

    gxl = ax.get_xgridlines()
    for lbl, gx in zip(labels, gxl):
        this_color = lbl.get_color()
        gx.set_color('gray')
        gx.set_linewidth(0.5)
        gx.set_linestyle("-")

    gyl = ax.get_ygridlines()
    for lbl, gx in zip(labels, gyl):
        this_color = lbl.get_color()
        gx.set_color('gray')
        gx.set_linewidth(0.5)
        gx.set_linestyle("-")

    ax.spines['polar'].set_visible(False)
    L=ax.legend()
    all_new_traits = []
    for ll in L.get_texts():
        new_text = nicer_label_traits[ll.get_text().replace('.sumstats','')]
        ll.set_text(new_text)
        ll.set_fontsize(10)
        all_new_traits.append(new_text)

    ax.get_legend().remove()

    buffer = np.pi/12/10
    # buffer = 0
    start = 2*np.pi -1/12*np.pi
    end = 2*np.pi + 1/12*np.pi

    pieces = {'argweave':  [2*np.pi -1/12*np.pi + buffer,  2*np.pi + 1/12*np.pi - buffer],
             'betascore': [ 2*np.pi + 1/12*np.pi + buffer, 2*np.pi +3/12*np.pi - buffer ],
             'fst_eas_afr': [ 2*np.pi +3/12*np.pi + buffer, 2*np.pi +9/12*np.pi - buffer  ],
             'fst_eur_afr': None,
             'fst_eur_eas': None,
             'gerp': [ 2*np.pi +9/12*np.pi +buffer ,  2*np.pi +17/12*np.pi - buffer    ],
             'linsigh': None,
             'phastCon100': None,
             'phyloP100': None,
             'xpehh_afr2_eas':  [ 2*np.pi +17/12*np.pi + buffer  ,2*np.pi +23/12*np.pi -buffer ],
             'xpehh_afr2_eur':None,
             'xpehh_eas_eur': None,}

    new_order_ = ['argweave','betascore',
                    'fst_eas_afr',
                    'fst_eur_afr',
                    'fst_eur_eas',
                    'gerp',
                    'linsigh',
                    'phastCon100',
                    'phyloP100',
                    'xpehh_afr2_eas',
                    'xpehh_afr2_eur',
                    'xpehh_eas_eur']


    for an in new_order_:
        start_end = pieces[an]
        color = anno_color_dict[an]

        if start_end:
            ax.fill_between(np.linspace(start_end[0], start_end[1], 100)[1:-1], max_y-0.25, max_y, color=color, zorder=10, alpha=0.2)
    # modify enrichment
    txa = [a for a in ax.texts if a.get_text() == 'enrichment'][0]
    txa.set_x(-0.44)
    txa.set_y(0.4)
    # txa.set_rotation(0.4)

    ax.plot(np.linspace(0, 2*np.pi,100), [0]*100, linewidth=1, linestyle='--',color='r', marker='o', markerfacecolor='r',  markersize=0, markeredgecolor='r',  markeredgewidth=0)


    ax.set_ylim(min_y, max_y)
    ax.set_yticks(np.arange(min_y, max_y+1, 1))
    ax.set_yticklabels([f'{int(x)}' for x in np.arange(min_y, max_y+1, 1)], color='grey', size=6)
    ax.set_rlabel_position(360-45)
    ax.set_aspect('equal', adjustable='box')

    trait_name = this_trait_df['nicer_trait_labe'].unique()[0]
    ax.set_title(trait_name, fontsize=8)

    # save_per_trait_fig(trait_name)




