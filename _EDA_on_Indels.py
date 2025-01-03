# %% [markdown]
# Exploratory Data Analysis on Indels (Fiber-seq & Element AVITI)
# %%
import io
import os
import numpy as np
import pandas as pd
import gzip as gz
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib_venn import venn2, venn3
from supervenn import supervenn # type: ignore

sns.set_theme(font="Arial", font_scale=1.15, style='ticks') 
matplotlib.rcParams['figure.dpi'] = 150
plt.rc("axes.spines", top=True, right=True)

def read_vcf(path):
    if path[-3:] == ".gz": 
        with gz.open(path, 'rb') as f:
            lines = [l.decode('utf-8') for l in f if not l.startswith(b'##')]
            return pd.read_csv(
                io.StringIO(''.join(lines)),
                dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                       'QUAL': str, 'FILTER': str, 'INFO': str},
                       sep='\t'
                       ).rename(columns={'#CHROM': 'CHROM'})
    else:
        with open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]
            return pd.read_csv(
                io.StringIO(''.join(lines)),
                dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                       'QUAL': str, 'FILTER': str, 'INFO': str},
                       sep='\t'
                       ).rename(columns={'#CHROM': 'CHROM'})

def vcf_format_getter(df, field) -> pd.Series:
    """
    Parameters
    ----------
    
    df : pandas.core.frame.DataFrame
        vcf read through read_vcf()
    field : str
        GT, VAF, DP, AD
    
    return pd.Series

    Example vcf structure
    ----------
    FORMAT              COLO829T_PassageB_DSA  
    GT:GQ:DP:AD:VAF:PL  1/1:59:133:2,131:0.984962:65,60,0
    GT:GQ:DP:AD:VAF:PL  1/1:65:131:0,131:1:70,66,0
    GT:GQ:DP:AD:VAF:PL  1/1:53:49:0,49:1:58,54,0

    """
    sampleid = df.columns[9]
    format = list(set(df['FORMAT'].values))[0].split(':')

    if field == "GT":
        gtindex = format.index('GT')
        return df[sampleid].str.split(':').apply(lambda x: x[gtindex])
    
    elif field == "VAF":
        vafindex = format.index('VAF')
        return df[sampleid].str.split(':').apply(lambda x: float(x[vafindex]))
    
    elif field == "DP":
        dpindex = format.index('DP')
        return df[sampleid].str.split(':').apply(lambda x: int(x[dpindex]))
    
    elif field == "AD":
        adindex = format.index('AD')
        return df[sampleid].str.split(':').apply(lambda x: int(x[adindex].split(',')[1]))
    
    else:
        raise ValueError("field should be one of GT, VAF, DP and AD!")

def vcf_info_parser(info_string) -> dict:
    """
    Parameters
    ----------
    info_string : str
    SP=STK,RF,MT2,VN;RGN=Difficult;RGN_T=Tier2;VAF_Ill=0.965;VAF_PB=1

    return dict

    """
    return {i.split('=')[0]: i.split('=')[1] for i in info_string.split(';')}

def vcf_info_getter(df, field):
    """
    Parameters
    ----------
    
    df : pandas.core.frame.DataFrame
        vcf read through read_vcf()
    field : str
        VAF_Ill, VAF_PB
    
    return pd.Series

    Example vcf structure
    ----------
    INFO
    SP=STK,RF,MT2,VN;RGN=Difficult;RGN_T=Tier2;VAF_Ill=0.965;VAF_PB=1
    SP=STK,RF,MT2,VN;RGN=Easy;RGN_T=Tier0;VAF_Ill=0.981;VAF_PB=0.984
    """

    if field == "RGN":
        return df['INFO'].apply(lambda x: vcf_info_parser(x)[field])
    
    elif field == "RGN_T":
        return df['INFO'].apply(lambda x: vcf_info_parser(x)[field])
    
    elif field == "VAF_Ill" or field == "VAF_PB":
        return df['INFO'].apply(lambda x: float(vcf_info_parser(x)[field]) if vcf_info_parser(x)[field] != 'NA' else np.nan)

    else:
        raise ValueError("field should be one of RGN, RGN_T, VAF_Ill and VAF_PB!")

def make_site_list(df: pd.DataFrame, path: str, prefix: str) -> None:
    """
    Parameters
    ----------
    df : pandas.core.frame.DataFrame 
        vcf read through read_vcf()
    path : str

    prefix : str

    return target_list 
    """
    
    df = df[['CHROM', 'POS', 'POS']].drop_duplicates()
    
    df.to_csv(f"{os.path.join(path, prefix)}.sitelist", sep='\t', index=False, header=False)

dir_fiberseq="/mmfs1/gscratch/stergachislab/mhsohny/SMaHT/Improving_SomaticVariantCalling_through_DSA/Fiber-seq/VariantCalls_DeepVariant_1.6.1"
dir_element="/mmfs1/gscratch/stergachislab/mhsohny/SMaHT/Improving_SomaticVariantCalling_through_DSA/Element/VariantCalls_DeepVariant_1.6.1"

DSA="/mmfs1/gscratch/stergachislab/mhsohny/SMaHT/DSA/DSA_COLO829BL_v3.0.0.fasta"

# %%
colobl_fiberseq_vcf = read_vcf(f"{dir_fiberseq}/COLO829BL/deepvariant/COLO829BL.deepvariant.split.indel.vcf.gz")
colotb_fiberseq_vcf = read_vcf(f"{dir_fiberseq}/COLO829T_PassageB_DSA/deepvariant/COLO829T_PassageB_DSA.deepvariant.split.indel.vcf.gz")
colota_fiberseq_vcf = read_vcf(f"{dir_fiberseq}/COLO829T_PassageA_DSA/deepvariant/COLO829T_PassageA_DSA.deepvariant.split.indel.vcf.gz")
colobl_element_vcf = read_vcf(f"{dir_element}/COLO829BL_Element_DSA/deepvariant/COLO829BL_Element_DSA.deepvariant.split.indel.vcf.gz")

colobl_fiberseq_indels_pass = colobl_fiberseq_vcf[(colobl_fiberseq_vcf['FILTER'] == 'PASS')].reset_index(drop=True)
colotb_fiberseq_indels_pass = colotb_fiberseq_vcf[(colotb_fiberseq_vcf['FILTER'] == 'PASS')].reset_index(drop=True)
colota_fiberseq_indels_pass = colota_fiberseq_vcf[(colota_fiberseq_vcf['FILTER'] == 'PASS')].reset_index(drop=True)
colobl_element_indels_pass = colobl_element_vcf[(colobl_element_vcf['FILTER'] == 'PASS')].reset_index(drop=True)

colobl_fiberseq_indels_pass['INDELid'] = colobl_fiberseq_indels_pass[['CHROM', 'POS', 'REF', 'ALT']].astype(str).apply('_'.join, axis=1)
colotb_fiberseq_indels_pass['INDELid'] = colotb_fiberseq_indels_pass[['CHROM', 'POS', 'REF', 'ALT']].astype(str).apply('_'.join, axis=1)
colota_fiberseq_indels_pass['INDELid'] = colota_fiberseq_indels_pass[['CHROM', 'POS', 'REF', 'ALT']].astype(str).apply('_'.join, axis=1)
colobl_element_indels_pass['INDELid'] = colobl_element_indels_pass[['CHROM', 'POS', 'REF', 'ALT']].astype(str).apply('_'.join, axis=1)

# %%
colobl_fiberseq_flagger = pd.read_table(f"{dir_fiberseq}/COLO829BL/deepvariant/COLO829BL.deepvariant.split.indel.flagger", sep='\t', header=None)
colobl_fiberseq_flagger.columns = ['INDELid', 'FILTER', 'Flagger']
colotb_fiberseq_flagger = pd.read_table(f"{dir_fiberseq}/COLO829T_PassageB_DSA/deepvariant/COLO829T_PassageB_DSA.deepvariant.split.indel.flagger", sep='\t', header=None)
colotb_fiberseq_flagger.columns = ['INDELid', 'FILTER', 'Flagger']
colota_fiberseq_flagger = pd.read_table(f"{dir_fiberseq}/COLO829T_PassageA_DSA/deepvariant/COLO829T_PassageA_DSA.deepvariant.split.indel.flagger", sep='\t', header=None)
colota_fiberseq_flagger.columns = ['INDELid', 'FILTER', 'Flagger']

colobl_element_flagger = pd.read_table(f"{dir_element}/COLO829BL_Element_DSA/deepvariant/COLO829BL_Element_DSA.deepvariant.split.indel.flagger", sep='\t', header=None)
colobl_element_flagger.columns = ['INDELid', 'FILTER', 'Flagger']

colobl_fiberseq_indels_pass_annot = pd.merge(colobl_fiberseq_indels_pass, colobl_fiberseq_flagger[['INDELid', 'Flagger']], on='INDELid', how='left')
colotb_fiberseq_indels_pass_annot = pd.merge(colotb_fiberseq_indels_pass, colotb_fiberseq_flagger[['INDELid', 'Flagger']], on='INDELid', how='left')
colota_fiberseq_indels_pass_annot = pd.merge(colota_fiberseq_indels_pass, colota_fiberseq_flagger[['INDELid', 'Flagger']], on='INDELid', how='left')

colobl_element_indels_pass_annot = pd.merge(colobl_element_indels_pass, colobl_element_flagger[['INDELid', 'Flagger']], on='INDELid', how='left')

print(colobl_fiberseq_indels_pass_annot['Flagger'].value_counts())
print(colotb_fiberseq_indels_pass_annot['Flagger'].value_counts())
print(colota_fiberseq_indels_pass_annot['Flagger'].value_counts())

print(colobl_element_indels_pass_annot['Flagger'].value_counts())

# %%
colobl_fiberseq_indels_pass_annot_set: set[str] = set(colobl_fiberseq_indels_pass_annot['INDELid'].values)
colotb_fiberseq_indels_pass_annot_set: set[str] = set(colotb_fiberseq_indels_pass_annot['INDELid'].values)
colota_fiberseq_indels_pass_annot_set: set[str] = set(colota_fiberseq_indels_pass_annot['INDELid'].values)

colobl_element_indels_pass_annot_set: set[str] = set(colobl_element_indels_pass_annot['INDELid'].values)

fig, ax = plt.subplots(1,1, figsize=(10,5), constrained_layout=False)
sets: list[set] = [colobl_fiberseq_indels_pass_annot_set, colotb_fiberseq_indels_pass_annot_set, colota_fiberseq_indels_pass_annot_set, colobl_element_indels_pass_annot_set]
labels: list[str] = ['COLO829BL Fiber-seq', 'COLO829T Passage B Fiber-seq', 'COLO829T Passage A Fiber-seq', 'COLO829BL Element']
supervenn(sets, 
          labels, 
          color_cycle=['green', 'red', 'blue', 'darkolivegreen'],
          widths_minmax_ratio=0.2,
          sets_ordering=None,
          col_annotations_area_height=1.1,
          rotate_col_annotations=True,
          side_plots=True,
          fontsize=10, 
          ax=ax)

ax.set_xlabel("Overlapped counts")
ax.set_ylabel("Samples")
#plt.figure(figsize=(20, 10))
#supervenn(sets_list, species_names, widths_minmax_ratio=0.1,
#          sets_ordering='minimize gaps', rotate_col_annotations=True, col_annotations_area_height=1.2)

# %%
colobl_fiberseq_indels_pass_annot_hap_set: set[str] = set(colobl_fiberseq_indels_pass_annot[colobl_fiberseq_indels_pass_annot['Flagger'] == 'Hap']['INDELid'].values)
colotb_fiberseq_indels_pass_annot_hap_set: set[str] = set(colotb_fiberseq_indels_pass_annot[colotb_fiberseq_indels_pass_annot['Flagger'] == 'Hap']['INDELid'].values)
colota_fiberseq_indels_pass_annot_hap_set: set[str] = set(colota_fiberseq_indels_pass_annot[colota_fiberseq_indels_pass_annot['Flagger'] == 'Hap']['INDELid'].values)

colobl_element_indels_pass_annot_hap_set: set[str] = set(colobl_element_indels_pass_annot[colobl_element_indels_pass_annot['Flagger'] == 'Hap']['INDELid'].values)

fig, ax = plt.subplots(1,1, figsize=(10,5), constrained_layout=False)
sets: list[set] = [colobl_fiberseq_indels_pass_annot_hap_set, colotb_fiberseq_indels_pass_annot_hap_set, colota_fiberseq_indels_pass_annot_hap_set, colobl_element_indels_pass_annot_hap_set]
labels: list[str] = ['COLO829BL Fiber-seq', 'COLO829T Passage B Fiber-seq', 'COLO829T Passage A Fiber-seq', 'COLO829BL Element']
supervenn(sets, 
          labels, 
          color_cycle=['green', 'red', 'blue', 'darkolivegreen'],
          widths_minmax_ratio=0.2,
          sets_ordering=None,
          col_annotations_area_height=1.1,
          rotate_col_annotations=True,
          side_plots=True,
          fontsize=10, 
          ax=ax)

ax.set_xlabel("Overlapped counts")
ax.set_ylabel("Samples")