
# rDoAMP v0.2.3
[![DOI](https://zenodo.org/badge/437441432.svg)](https://zenodo.org/badge/latestdoi/437441432)

<!-- badges: start -->
<!-- badges: end -->

興味ある配列データ（特定の分類群の配列など）が指定したプライマーで増幅しうるかどうかをチェックするためのパッケージです。`rentrez` を使って検索ワードにヒットした配列を Entrez (https://www.ncbi.nlm.nih.gov/Web/Search/entrezfs.html) からダウンロードし、`seqkit amplicon` で指定したプライマーで増幅しうるかどうかチェックします。

macOS と Windows で動作確認しています。

## 事前にインストールする必要があるパッケージ
- `seqkit` (https://bioinf.shenwei.me/seqkit/)
  - インストール: https://bioinf.shenwei.me/seqkit/download/ から対応する OS の *.tar.gz ファイルをダウンロード
    - macOS
      1. `tar -zxvf seqkit_darwin_amd64.tar.gz` でファイルを解凍
      2. `sudo cp seqkit /usr/local/bin` でファイルをコピー
      3. `seqkit version` でバージョン確認
    - Windows
      1. seqkit_windows_*.exe.tar.gz を解凍
      2. `seqkit.exe` を `C:\WINDOWS\system32` にコピー
      3. `seqkit version` でバージョン確認

## 依存パッケージ
- `rentrez` (https://github.com/ropensci/rentrez)

## インストール

ベータ版はRコンソール内で以下のコードを入力するとインストールできます。

``` r
# install.packages("devtools")
devtools::install_github("ong8181/rDoAMP")
```

## クイックスタート
Rで以下を入力すると検索ワードでヒットしてダウンロードされた配列の中で、指定したプライマーで増幅しうる配列が出力フォルダに保存されます。

```r
# Load package
library(rDoAMP); packageVersion("rDoAMP")

# Display Help
?doamp_auto
?doamp_custom
?expand_degenerate_primer

# Load primer set
data(primer_set)

# Download sequence data and check whether MiFish primer set can amplify downloaded sequences
doamp_auto("Trachurus AND mitochondrion AND 1000:20000[SLEN]",
           F_primer = primer_set$Animal12S_MiFish_U$forward,   # MiFish-U-F
           R_primer = primer_set$Animal12S_MiFish_U$reverse,   # MiFish-U-R
           n_mismatch = 3,
           overwrite_output_dir = FALSE)
```

#### 出力フォルダに保存されるファイル
- `download.fa` : 指定した検索ワードによりダウンロードされた配列データ
- `amplified.fa`: プライマーとマッチした配列
- `parameter_list.txt`: 解析に使用したパラメータリスト
- `expanded_primer_list.tsv` : あり得るプライマーの組み合わせリスト (縮重プライマー使用時) 
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
           n_retidmax = n_retmax * 10,
           n_mismatch = 0,
           random_sampling = TRUE,
           random_sampling_seed = 1234,
           output_dir = "rDoAMP_Out",
           save_parameter = TRUE,
           save_stat = TRUE,
           overwrite_output_dir = FALSE)
```

- `search_query` : データベースの検索ワード (例. "Trachurus AND mitochondrion AND 1000:20000[SLEN]" など)
- `F_primer` : フォワードプライマー配列
- `R_primer`: リバースプライマー配列
- `n_retmax`: 取得する配列の最大数。今の所、大きい値（だいたい > 500）を指定するとエラーとなります。
- `n_retidmax`: 取得する配列 ID の数。デフォルトでは取得予定の配列数の10倍の ID を取得し、その中からランダムに配列をダウンロード。ダウンロードする配列の多様性を担保するためのオプション。
- `n_mismatch` : 許容するプライマー - 鋳型間のミスマッチの数
- `random_sampling`: 多めの配列 ID を取得して、その中からランダムに配列を取得するかどうか
- `random_sampling_seed`: ランダムサンプリングのシード値
- `output_dir`: 出力フォルダの名前。デフォルトは "rDoAMP_Out"。
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
R functions to extract amplicons from target sequences using a user-specified primer set. - DOes my primer set AMPlify my targets? -

(Convenient wrapper functions for `seqkit amplicon`)

Executable in macOS and Windows.

## Prerequisites
- `seqkit` (https://bioinf.shenwei.me/seqkit/)
  - Installation: https://bioinf.shenwei.me/seqkit/download/ 


## Dependencies
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
library(rDoAMP); packageVersion("rDoAMP")

# Display Help
?doamp_auto
?doamp_custom
?expand_degenerate_primer

# Load primer set
data(primer_set)

# Download sequence data and check whether MiFish primer set can amplify downloaded sequences
doamp_auto("Trachurus AND mitochondrion AND 1000:20000[SLEN]",
           F_primer = primer_set$Animal12S_MiFish_U$forward,   # MiFish-U-F
           R_primer = primer_set$Animal12S_MiFish_U$reverse,   # MiFish-U-R
           n_mismatch = 3,
           overwrite_output_dir = FALSE)
```

#### Output files
- `download.fa` : Downloaded sequence data using a user-specified search query
- `amplified.fa`: Sequence data that could be amplified by a specified primer set
- `parameter_list.txt`: Parameter list used in the analysis
- `expanded_primer_list.tsv` : A list of primer combinations (only for degenerated primers)
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
           n_retidmax = n_retmax * 10,
           n_mismatch = 0,
           random_sampling = TRUE,
           random_sampling_seed = 1234,
           output_dir = "rDoAMP_Out",
           save_parameter = TRUE,
           save_stat = TRUE,
           overwrite_output_dir = FALSE)
```

- `search_query` : Search query for Entrez (e.g., "Trachurus AND mitochondrion AND 1000:20000[SLEN]")
- `F_primer` : Forward primer sequence
- `R_primer`: Reverse primer sequence
- `n_retmax`: The maximum number of sequences retrieved from Entrez. Current version does not accept a large number (e.g., > 500).
- `n_retidmax`: The maximum number of IDs collected from Entrez. Among the IDs, `n_retmax` sequences are retrieved. This option is set to increase the diversity of sequences retrieved.
- `n_mismatch` : The maximum number of primer-template mismatches allowed
- `random_sampling`: Logical. If TURE, n_retidmax IDs collected, and then n_retmax sequences are randomly collected
- `random_sampling_seed`: Random number seed for random_sampling
- `output_dir`: Output directory name
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
