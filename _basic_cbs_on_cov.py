import sys, os
import numpy as np
import pandas as pd
from linear_segment import segment
import ailist

try:
    inputdir = os.path.dirname(sys.argv[1])
    inputf = os.path.basename(sys.argv[1]) # mosdepth region-based file (gzipped) ex) COLO829BL_DSA_resetmapq.regions.bed.gz 
    contigid = sys.argv[2] # ex) haplotype2-0000195
    outputdir = os.path.join(inputdir, "CBS")
    os.system(f'mkdir -p {outputdir}')
    output_prefix = os.path.join(outputdir, inputf.rstrip(".bed.gz") + '.' + contigid) # bed file 
except IndexError:
    print("\nMissing Input Arguments!!\n")
    sys.exit()
    
def get_cbs_segment(df: pd.DataFrame, contigid: str) -> pd.DataFrame:
    """
    Applying Circular Binary Segmentation (CBS) to the Coverage Bed file (from mosdepth)
    Parameters
    ----------
    df : pandas.core.frame.DataFrame 
    gzipped regions.bed file from `mosdepth --by <size>` output
    df = pd.read_table(regions.bed.gz, sep="\t", names=['id', 'start', 'end', 'cov'], compression="gzip")
    
    contigid : str
    Specified contig name of the genome

    return pd.DataFrame of CBS result

    Example
    ----------
  
                        id  start   end cov
    0   haplotype1-0000001  0   1000    11.91
    1   haplotype1-0000001  1000    2000    20.17
    2   haplotype1-0000001  2000    3000    31.23
    3   haplotype1-0000001  3000    4000    39.86

    haplotype2-0000195      0       6000    2.01    0.69
    haplotype2-0000195      6000    14000   20.07625    19.49
    haplotype2-0000195      14000   29000   50.912  52.78
    haplotype2-0000195      29000   35000   39.54833    39.81
    haplotype2-0000195      35000   44735   18.506  18.255
    """
    df = df[df['id'] == contigid]
    cov_array = np.array(df['cov'])
    pseudo_labels = np.repeat('na', len(cov_array))

    segments = segment(cov_array, pseudo_labels, method="cbs")

    combined_df = pd.DataFrame()

    for i in range(len(segments)):
        seg = segments[i]
        seg_start, seg_end = seg.start, seg.end
        chrom = df.iloc[seg_start]['id']
        start = df.iloc[seg_start]['start']
        end = df.iloc[seg_end -1]['end']
        
        MeanCov = df.iloc[seg_start: seg_end]['cov'].mean()
        MedianCov = df.iloc[seg_start: seg_end]['cov'].median()
        
        seg_df = pd.DataFrame([[chrom, start, end, round(MeanCov, 5), round(MedianCov, 5)]])

        combined_df = pd.concat([combined_df, seg_df], ignore_index=True)

    return combined_df

df = pd.read_table(os.path.join(inputdir, inputf), sep="\t", names=['id', 'start', 'end', 'cov'], compression="gzip")

df_segment = get_cbs_segment(df, contigid)

df_segment.to_csv(f'{output_prefix}.cbs.bed', 
                  sep="\t",
                  header=False,
                  index=False)

bgzip="/mmfs1/gscratch/stergachislab/mhsohny/Miniconda/bin/bgzip"
tabix="/mmfs1/gscratch/stergachislab/mhsohny/Miniconda/bin/tabix"
os.system(f"{bgzip} -@ 5 {output_prefix}.cbs.bed && {tabix} -p bed {output_prefix}.cbs.bed.gz")