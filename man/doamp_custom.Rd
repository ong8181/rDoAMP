% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rDoAMP_funcs.R
\name{doamp_custom}
\alias{doamp_custom}
\title{Automatically extract potential amplicons from target sequences}
\usage{
doamp_custom(
  target_fasta,
  F_primer,
  R_primer,
  n_mismatch = 0,
  output_dir = "rDoAMP_Out",
  save_parameter = TRUE,
  save_stat = TRUE,
  overwrite_output_dir = FALSE
)
}
\arguments{
\item{target_fasta}{Character. Path to a FASTA file that contains target sequences}

\item{F_primer}{Character. Forward primer sequence. Degenerate base allowed.}

\item{R_primer}{Character. Reverse primer sequence. Degenerate base allowed.}

\item{n_mismatch}{Numeric. The maximum number of primer-template mismatches allowed}

\item{output_dir}{Character. Output directory name}

\item{save_parameter}{Logical. If TRUE, parameters used in the analysis saved}

\item{save_stat}{Logical. If TRUE, summary of the analysis saved}

\item{overwrite_output_dir}{Logical. If TRUE, overwrite the contents of output directory.}
}
\description{
\code{doamp_custom} Automatically extract potential amplicons from target sequences
}
