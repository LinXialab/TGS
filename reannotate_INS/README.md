## System requirements

`python 3.7.4`

python package `biopython` and `pandas`

`snakemake 6.2.1` (tested version)
`sdust 01(r2)`
`TRF 4.09.1`
`needle`, `EMBOSS 6.6.0.0`

## Installation guide

Install via `conda` by
`conda install -c bioconda snakemake emboss biopython pandas`

Install `sdust` as in
`https://github.com/lh3/sdust`

Install `TRF` as in
`https://github.com/Benson-Genomics-Lab/TRF`

## Demo

The demo contain 6 insertions with their `fasta` sequence and `RepeatMasker` results.

After download, change the column 2 value in `test_sample.txt` according to the path of the demo folder.

Takes 7 mins to run on 6 cores.

The final out file is `all.INS.sdust.trf.replaced.cor.type.final.tsv`.

## Instructions for use

1. Run `RepeatMasker` on all insertions. The `fasta` filename of every insertion should be `{insertion id}.fasta`. The input `fasta` sequence is the insertion sequence with flank `50 bp` reference genome sequence on both sides.

   The result folder structure should be as below.

   ```markdown
   ├── insertion1
   │	├── insertion1.fasta
   │   ├── RepeatMasker
   │   │   ├── insertion1.fasta.cat 
   │   │	├── insertion1.fasta.masked
   │   │	├── insertion1.fasta.out
   │   │	├── insertion1.fasta.tbl
   ```
   
2. Specify the path to the file recording information of the insertions to test in `Snakefile` `insert_df_path = xxx`. The file should has three columns with tab delimiter and no header. The content of each column is described below. 
| Column 1 | insertion id |
| :------------- | :------------- |
| **Column 2** | **path to the folder contain the `RepeatMasker` result of the insertion (a subfolder named `RepeatMasker`) and the `fasta` sequence of the insertion** |
| **Column 3** | **precise insert site on reference genome of the insertion** |

3. Specify the executable path of `TRF`, `sdust` and `needle` in `config.yaml` file `tool` section.

4. Specify the path to the `fasta` file of the reference genome in `Snakefile` `ref_file = xxx`.

5. Specify the result folder in `Snakefile` `workdir: xxx`.

6. Run the pipeline as `snakemake --cores xx`.

7. Output:

   The final output file is named `all.INS.sdust.trf.replaced.cor.type.final.tsv`. The explanation of some key columns is below.

   | trf_seq | repeat pattern found by TRF                                  |
   | :------------- | :------------- |
   | **replaced** | **`RepeatMasker` result is replaced or not** |
   | **segment_dup** | **the insertion is a one time duplication or not** |
   | **repeat_dup** | **the insertion is a multiple times duplication of a sequence pattern or not** |
   | **clean_decision** | **the final annotation of the insertion sequence with the coordinate information** |









































































