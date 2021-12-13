
# rDoAMP v0.1.1

<!-- badges: start -->
<!-- badges: end -->

興味ある配列データ（特定の分類群の配列など）が指定したプライマーで増幅しうるかどうかをチェックするためのパッケージです。`rentrez` (https://github.com/ropensci/rentrez) を使って興味ある配列を Entrez (https://www.ncbi.nlm.nih.gov/Web/Search/entrezfs.html) からダウンロードし、`seqkit` (https://bioinf.shenwei.me/seqkit/) を用いて指定したプライマーで増幅しうるかどうかチェックします。

## 事前にインストールする必要があるパッケージ
- `seqkit` (https://bioinf.shenwei.me/seqkit/)
- `rentrez` (https://github.com/ropensci/rentrez)

## インストール

Rコンソール内で以下のコードを入力するとインストールできます。

``` r
# install.packages("devtools")
devtools::install_github("ong8181/rDoAMP")
```

## クイックスタート
Rで以下を入力すると検索ワードでヒットしてダウンロードされた配列の中で、指定したプライマーで増幅しうる配列が出力フォルダに保存されます。

```r
# Load package
library(rDoAMP); packageVersion("rDoAMP") # 0.1.1, 2021.12.12

# Display Help
?doamp_auto
?doamp_custom
?expand_degenerate_primer

# Load primer set
data(primer_set)

# Download sequence data and check whether MiFish primer set can amplify downloaded sequences
doamp_auto("Trachurus AND mitochondrion AND 1000:20000[SLEN]",
           F_primer = primer_set$MiFish_U$forward,   # MiFish-U-F
           R_primer = primer_set$MiFish_U$reverse,   # MiFish-U-R
           n_mismatch = 3)
```

#### 出力フォルダに保存されるファイル
- `download.fa` (`doamp_auto()` 使用時) : 指定した検索ワードによりダウンロードされた配列データ
- `custom_db.fa` (`doamp_custom()` 使用時) : ユーザーが指定した増幅するかどうか調べてたい配列データ
- `amplified.fa`: プライマーとマッチした配列
- `parameter_list.txt`: 解析に使用したパラメータリスト
- `expanded_primer_list.tsv` (縮重プライマー使用時) : あり得るプライマーの組み合わせリスト 
- `stat.tsv`: ダウンロードした配列数、プライマーとマッチした配列数などの統計情報


# 詳しいマニュアル
## `doamp_auto()`
-  `rentrez` により検索ワードにヒットした配列をダウンロードします。
- 予め多めの配列の ID を検索し、その中から指定した数の配列をランダムに取得します。
- プライマー配列にマッチするものを `seqkit amplicon` により抽出します。
- 指定したプライマーに縮重塩基が含まれていても自動であり得る全てのプライマーの組み合わせをリストアップします。これにより`seqkit amplicon` の `--max-mismatch` も使用可能にします。

#### 引数
```r
doamp_auto(search_query,
           F_primer,
           R_primer,
           n_retmax = 20,
           n_mismatch = 0,
           output_dir = "rDoAMP_Out",
           random_sampling = TRUE,
           random_sampling_seed = 1234,
           n_retidmax = n_retmax * 10,
           save_parameter = TRUE,
           save_stat = TRUE,
           overwrite_output_dir = FALSE)
```

- `search_query` : データベースの検索ワード ("Trachurus AND mitochondrion AND 1000:20000[SLEN]" など)
- `F_primer` : フォワードプライマー配列
- `R_primer`: リバースプライマー配列
- `n_retmax`: 取得する配列の最大数 (v0.1.1 では大きすぎる配列数はエラーとなります)
- `n_mismatch` : 許容するプライマー - 鋳型間のミスマッチの数
- `output_dir`: 出力フォルダの名前。デフォルトは "rDoAMP_Out"。
- `random_sampling`: 多めの配列 ID を取得して、その中からランダムに配列を取得するかどうか
- `random_sampling_seed`: ランダムサンプリングのシード値
- `n_retidmax`: 取得する配列 ID の数。デフォルトでは取得予定の配列数の10倍の ID を取得し、その中からランダムに配列をダウンロード。
- `save_parameter`: 解析に使用したパラメータを保存するかどうか
- `save_stat`: 解析結果のサマリーを保存するかどうか
- `overwrite_output_dir`: 出力フォルダがすでに存在していた場合にフォルダを上書きするかどうか。上書きする場合は既存フォルダの中身は消去されるので注意。

#### 基本的な使用方法
```r
doamp_auto("Trachurus AND mitochondrion AND 1000:20000[SLEN]",
           F_primer = "GTCGGTAAAACTCGTGCCAGC",                    # MiFish-U-F
           R_primer = "CATAGTGGGGTATCTAATCCCAGTTTG",              # MiFish-U-R
           n_mismatch = 3)
```

## `doamp_custom()`
- ユーザーが指定した FASTA ファイルからプライマー配列にマッチするものを `seqkit amplicon` により抽出します。
- 指定したプライマーに縮重塩基が含まれていても自動であり得る全てのプライマーの組み合わせをリストアップします。これにより`seqkit amplicon` の `--max-mismatch` も使用可能にします。

#### 引数
```r
doamp_custom(target_fasta,
             F_primer,
             R_primer,
             n_mismatch = 0,
             output_dir = "rDoAMP_Out",
             save_parameter = TRUE,
             save_stat = TRUE,
             overwrite_output_dir = FALSE) 
```

- `farget_fasta` : ユーザ指定の配列データ
- `F_primer` : フォワードプライマー配列
- `R_primer`: リバースプライマー配列
- `n_mismatch` : 許容するプライマー - 鋳型間のミスマッチの数
- `output_dir`: 出力フォルダの名前。デフォルトは "rDoAMP_Out"。
- `save_parameter`: 解析に使用したパラメータを保存するかどうか
- `save_stat`: 解析結果のサマリーを保存するかどうか
- `overwrite_output_dir`: 出力フォルダがすでに存在していた場合にフォルダを上書きするかどうか。上書きする場合は既存フォルダの中身は消去されるので注意。

#### 基本的な使用方法
```r
doamp_custom("YOUR_FASTA.fasta",
             F_primer = "GTCGGTAAAACTCGTGCCAGC",                    # MiFish-U-F
             R_primer = "CATAGTGGGGTATCTAATCCCAGTTTG",              # MiFish-U-R
             n_mismatch = 3)
```

## データセット `primer_set`
- よく使用されるプライマーのリストです。`data(primer_set)` と入力してロードして下さい。


# Quick tutorial
R functions to extract amplicons from target sequences using a user-specified primer set.

(Convenient wrapper functions for `seqkit amplicon`)

## Prerequisites
- `seqkit` (https://bioinf.shenwei.me/seqkit/)
- `rentrez` (https://github.com/ropensci/rentrez)

## Installation

You can install the development version of rDoAMP from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ong8181/rDoAMP")
```

## Quick start
```r
# Load package
library(rDoAMP); packageVersion("rDoAMP") # 0.1.1, 2021.12.12

# Display Help
?doamp_auto
?doamp_custom
?expand_degenerate_primer

# Load primer set
data(primer_set)

# Download sequence data and check whether MiFish primer set can amplify downloaded sequences
doamp_auto("Trachurus AND mitochondrion AND 1000:20000[SLEN]",
           F_primer = primer_set$MiFish_U$forward,   # MiFish-U-F
           R_primer = primer_set$MiFish_U$reverse,   # MiFish-U-R
           n_mismatch = 3)
```

#### Output files
- `download.fa` (only for `doamp_auto()`) : Downloaded sequence data using a user-specified search query
- `custom_db.fa` (only for `doamp_custom()`) : User-specified sequence data to be checked by `seqkit amplicon`
- `amplified.fa`: Sequence data that could be amplified by a specified primer set
- `parameter_list.txt`: Parameter list used in the analysis
- `expanded_primer_list.tsv` (only for degenerated primers) : A list of primer combinations
- `stat.tsv`: Statistical information such as the numbers of downloaded and amplified sequences


# Detailed manual
## `doamp_auto()`
- Download sequences using `rentrez` package of R.
- Allow random sampling from searched sequence IDs.
- Extract amplicons using `seqkit amplicon`.
- Automatically expand degenerate primers and create a list of primer combinations to use `--max-mismatch` option of `seqkit amplicon` for degenerated primers.

#### Arguments
```r
doamp_auto(search_query,
           F_primer,
           R_primer,
           n_retmax = 20,
           n_mismatch = 0,
           output_dir = "rDoAMP_Out",
           random_sampling = TRUE,
           random_sampling_seed = 1234,
           n_retidmax = n_retmax * 10,
           save_parameter = TRUE,
           save_stat = TRUE,
           overwrite_output_dir = FALSE)
```

- `search_query` : Search query for Entrez
- `F_primer` : Forward primer sequence
- `R_primer`: Reverse primer sequence
- `n_retmax`: The maximum number of sequences retrieved from Entrez
- `n_mismatch` : The maximum number of primer-template mismatches allowed
- `output_dir`: Output directory name
- `random_sampling`: Logical. If TURE, n_retidmax IDs collected, and then n_retmax sequences are randomly collected
- `random_sampling_seed`: Random number seed for random_sampling
- `n_retidmax`: The maximum number of IDs collected from Entrez. Among the IDs, `n_retmax` sequences are retrieved 
- `save_parameter`: Logical. If TRUE, parameters used in the analysis saved
- `save_stat`: Logical. If TRUE, summary of the analysis saved
- `overwrite_output_dir`: Logical. If TRUE, overwrite the contents of output directory (the original contents of the output file will be deleted)


#### Basic usage
```r
doamp_auto("Trachurus AND mitochondrion AND 1000:20000[SLEN]",
           F_primer = "GTCGGTAAAACTCGTGCCAGC",                    # MiFish-U-F
           R_primer = "CATAGTGGGGTATCTAATCCCAGTTTG",              # MiFish-U-R
           n_mismatch = 3)
```

## `doamp_custom()`
- Extract amplicons from a user-specified, custom FASTA file using `seqkit amplicon`.
- Automatically expand degenerate primers and create a list of primer combinations to use `--max-mismatch` option of `seqkit amplicon` for degenerated primers.

#### Arguments
```r
doamp_custom(target_fasta,
             F_primer,
             R_primer,
             n_mismatch = 0,
             output_dir = "rDoAMP_Out",
             save_parameter = TRUE,
             save_stat = TRUE,
             overwrite_output_dir = FALSE) 
```

- `target_fasta` : FASTA file that contains target sequences
- `F_primer` : Forward primer sequence
- `R_primer`: Reverse primer sequence
- `n_mismatch` : The maximum number of primer-template mismatches allowed
- `output_dir`: Output directory name
- `save_parameter`: Logical. If TRUE, parameters used in the analysis saved
- `save_stat`: Logical. If TRUE, summary of the analysis saved
- `overwrite_output_dir`: Logical. If TRUE, overwrite the contents of output directory (the original contents of the output file will be deleted)


#### Basic usage
```r
doamp_custom("YOUR_FASTA.fasta",
             F_primer = "GTCGGTAAAACTCGTGCCAGC",                    # MiFish-U-F
             R_primer = "CATAGTGGGGTATCTAATCCCAGTTTG",              # MiFish-U-R
             n_mismatch = 3)
```

## `primer_set`
- A list of popular primer sets. Please load by typing `data(primer_set)`.
