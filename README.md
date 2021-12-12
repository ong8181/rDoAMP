
# rDoAMP v0.1.1

<!-- badges: start -->
<!-- badges: end -->

興味ある配列データ（特定の分類群の配列など）が指定したプライマーで増幅しうるかどうかをチェックするための関数です。`rentrez` (https://github.com/ropensci/rentrez) を使って興味ある配列を Entrez (https://www.ncbi.nlm.nih.gov/Web/Search/entrezfs.html) からダウンロードし、`seqkit` (https://bioinf.shenwei.me/seqkit/) を用いて指定したプライマーで増幅しうるかどうかチェックします。

(`seqkit` の機能の一つ、`seqkit amplicon` を使用しやすくしたラッパー関数です)

## 事前にインストールする必要があるパッケージ
- `seqkit` (https://bioinf.shenwei.me/seqkit/)
- `rentrez` (https://github.com/ropensci/rentrez)

## インストール

開発中のバージョンはRコンソール内で以下のコードを入力するとインストールできます。

``` r
# install.packages("devtools")
devtools::install_github("ong8181/rDoAMP")
```

## クイックスタート
Rで以下を入力するとダウンロードした配列･その中で指定したプライマーで増幅しうる配列が出力フォルダに保存されます。

```r
# Load package
library(rDoAMP); packageVersion("rDoAMP") # 0.1.1, 2021.12.12

# Load primer set
data(primer_set)

# Download sequence data and check whether MiFish primer set can amplify downloaded sequences
doamp_auto("Trachurus AND mitochondrion AND 1000:20000[SLEN]",
           F_primer = primer_set$MiFish_U$forward,   # MiFish-U-F
           R_primer = primer_set$MiFish_U$reverse,   # MiFish-U-R
           n_mismatch = 3)
```


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

# Load primer set
data(primer_set)

# Download sequence data and check whether MiFish primer set can amplify downloaded sequences
doamp_auto("Trachurus AND mitochondrion AND 1000:20000[SLEN]",
           F_primer = primer_set$MiFish_U$forward,   # MiFish-U-F
           R_primer = primer_set$MiFish_U$reverse,   # MiFish-U-R
           n_mismatch = 3)
```


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

#### Basic usage
```r
doamp_custom("YOUR_FASTA.fasta",
             F_primer = "GTCGGTAAAACTCGTGCCAGC",                    # MiFish-U-F
             R_primer = "CATAGTGGGGTATCTAATCCCAGTTTG",              # MiFish-U-R
             n_mismatch = 3)
```

## `primer_set`
- A list of popular primer sets. Please load by typing `data(primer_set)`.
