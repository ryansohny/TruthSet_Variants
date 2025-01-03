import io
import numpy as np
import pandas as pd
import gzip as gz
import seaborn as sns
import scipy.stats as stats
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

sns.set_theme(font="Arial", font_scale=1.15, style='ticks') 
matplotlib.rcParams['figure.dpi'] = 300
plt.rc("axes.spines", top=False, right=False)

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

# %%
dsa_hap_file="/mmfs1/gscratch/stergachislab/mhsohny/SMaHT/DSA/Flagger_v0.3.3/DSA_COLO829BL_v3.0.0.alt_removed.flagger_final_hap.bed.gz"
dsa_hap = pd.read_csv(dsa_hap_file, sep='\t', header=None)
dsa_hap.columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb']

# %%
repeatmasker_file="/mmfs1/gscratch/stergachislab/mhsohny/SMaHT/Improving_SomaticVariantCalling_through_DSA/Fiber-seq/VariantCalls_DeepVariant_1.6.1/Annotation/SNVs/Rhodonite/Intersect_RM_WITH_DSA_COLO829BL_v3.0.0.alt_removed.flagger_final_hap_sorted.bed.gz"

with gz.open(repeatmasker_file, 'rt') as dfh:
    repeatmasker = pd.read_csv(dfh, sep='\t', header=None)
    repeatmasker.columns = ['chrom', 'start', 'end', 'matchingrep', 'repsize', 'matchingstrand', 'repClass', 'repFamily', 'na', 'repID']

#haplotype1-0000001      6202    6636    L3      434     +       LINE    CR1     -1.0    3
#haplotype1-0000001      6761    6910    Plat_L3 149     +       LINE    CR1     -1.0    4
#haplotype1-0000001      8179    8575    MLT1K   396     +       LTR     ERVL-MaLR       -1.0    5
#haplotype1-0000001      9349    9570    MIR     221     -       SINE    MIR     -1.0    6
#haplotype1-0000001      10055   10242   L2b     187     +       LINE    L2      -1.0    7
#haplotype1-0000001      10318   10423   MIR     105     +       SINE    MIR     -1.0    8

# %%
# Size of the None-RE region
# bedtools subtract -a DSA_COLO829BL_v3.0.0.alt_removed.flagger_final_hap.bed.gz -b RM.bed.gz | awk 'BEGIN {sum=0} {sum+=$3-$2} END {print sum}'
# 2777444338
rec_df = pd.DataFrame(repeatmasker.groupby("repClass")["repsize"].sum())
rec_df.reset_index(inplace=True)
rec_df = pd.concat([pd.DataFrame([['None_RE', 2777444338]], columns=rec_df.columns), rec_df], axis=0).reset_index(drop=True)

rec_df['Exp_Proportion'] = rec_df['repsize'] / rec_df['repsize'].sum()

# %%
fig, ax = plt.subplots(1, 1, figsize=(10, 7), constrained_layout=True)
rec_df.plot(x='repClass', y='repsize', logy=True, kind="bar", color="grey", edgecolor="black", legend=None, ax=ax)
ax.set_title("Total Size of each RE class in DSAv3.0.0 (Flagger Hap)", fontsize=16)
ax.set_xlabel("Repeat Element Class", fontsize=14)
ax.set_ylabel("Total Size", fontsize=14)
ax.tick_params(axis='x', rotation=55, labelsize=12)

for p in ax.patches:
    ax.annotate(f'{p.get_height():,.0f}', (p.get_x() + p.get_width() / 2., p.get_height()),
                ha='center', 
                va='center', 
                xytext=(0, 10), 
                fontsize=5.5, 
                textcoords='offset points')

# %% [markdown]
## COLO829T (Passage B) SNVs
# %%
colotb_snvs = read_vcf("/mmfs1/gscratch/stergachislab/mhsohny/SMaHT/Improving_SomaticVariantCalling_through_DSA/Fiber-seq/VariantCalls_DeepVariant_1.6.1/Mutational_Spectrum/COLO829T_PassageB_DSA.deepvariant.split.snv.modified.final.tba.vcf.gz")

colotb_snvs_RM = pd.read_csv("/mmfs1/gscratch/stergachislab/mhsohny/SMaHT/Improving_SomaticVariantCalling_through_DSA/Fiber-seq/VariantCalls_DeepVariant_1.6.1/Annotation/SNVs/Rhodonite/Intersect-wao_COLO829T_PassageB_DSA.deepvariant.split.snv.modified.final.tba.vcf_WITH_RM.tab", sep="\t", header=None)

# %%
#colotb_snvs_RM[colotb_snvs_RM.iloc[:, 10] != "."]

rec_snv = pd.DataFrame(colotb_snvs_RM.groupby(colotb_snvs_RM.columns[16]).size(), columns=["SNV_count"])
rec_snv.index.name = "repClass"
rec_snv.reset_index(inplace=True)
rec_snv.at[0, 'repClass'] = 'None_RE'
rec_snv['Obs_Proportion'] = rec_snv['SNV_count'] / rec_snv['SNV_count'].sum()

fig, ax = plt.subplots(1,1, figsize=(10, 7), constrained_layout=True)
rec_snv.plot(kind="bar", x="repClass", y="SNV_count", logy=True, color=sns.color_palette(palette='terrain', n_colors=rec_snv.shape[0]), legend=None, ax=ax)
ax.set_title("SNV in different REs", fontsize=16)
ax.set_xlabel("Repeat Element Class", fontsize=14)
ax.set_ylabel("The number of SNV", fontsize=14)
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:,.0f}'.format(x)))
new_labels = [i.get_text() for i in ax.get_xticklabels()]
new_labels[0] = 'None RE'
ax.set_xticklabels(new_labels)
ax.tick_params(axis='x', rotation=55, labelsize=12)

for p in ax.patches:
    ax.annotate(f'{p.get_height():,.0f}', (p.get_x() + p.get_width() / 2., p.get_height()),
                ha='center', 
                va='center', 
                xytext=(0, 10), 
                fontsize=10, 
                textcoords='offset points')

# %%
rec_snv_obs_exp = pd.merge(rec_df, rec_snv, on="repClass", how="inner")
rec_snv_obs_exp['Obs-Exp_Ratio'] = rec_snv_obs_exp['Obs_Proportion'] / rec_snv_obs_exp['Exp_Proportion']

fig, ax = plt.subplots(1, 1, figsize=(10, 7), constrained_layout=True)
rec_snv_obs_exp.plot(kind="bar", x='repClass', y='Obs-Exp_Ratio', legend=None, color=sns.color_palette(palette='terrain', n_colors=rec_snv_obs_exp.shape[0]), ax=ax)
ax.tick_params(axis='x', rotation=55, labelsize=12)
ax.set_xlabel("Repeat Element Class", fontsize=14)
ax.set_ylabel("Observed/Expected Ratio", fontsize=14)
ax.hlines(y=1, xmin=-1, xmax=rec_snv_obs_exp.shape[0], color='red', linestyle='--')

for p in ax.patches:
    ax.annotate(f'{p.get_height():.2f}', (p.get_x() + p.get_width() / 2., p.get_height()),
                ha='center', 
                va='center', 
                xytext=(0, 10), 
                fontsize=10, 
                textcoords='offset points')
    


# %%
# for None_RE: contingency_table: list[list[int, int], list[int, int]] = [[19299, 2777425039], [32475, 3281919225]]
#contingency_table = [[8441, 422289847],[43333, 5637106191]]
contingency_table = [[19299, 2777425039], [32475, 3281919225]]
chi2, p, dof, expected = stats.chi2_contingency(contingency_table)
n = np.sum(contingency_table)
cramersv = np.sqrt(chi2 / (n * (min(len(contingency_table), len(contingency_table[0])) - 1)))
print(n)
print(cramersv)