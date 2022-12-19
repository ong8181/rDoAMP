#' @title Visualize sequence length distribution
#' @description \code{seq_hist} Visualize sequence length distribution (a wrapper function of Bio).
#' @param path_to_fasta Path to FASTA file.
#' @return ghist Histogram as a ggplot object
#' @export
#' @examples
#' # seq_hist(path_to_fasta)

seq_hist <- function(path_to_fasta) {
  # Read FASTA file
  seq_d <- Biostrings::readDNAStringSet(path_to_fasta)
  s_len <- BiocGenerics::width(seq_d)

  # Draw figure
  ghist <- ggplot2::ggplot(data.frame(x = s_len), ggplot2::aes_string(x = "x")) +
    ggplot2::geom_histogram() +
    ggplot2::geom_vline(xintercept = mean(s_len), color = "red3", linetype = 1) +
    ggplot2::geom_vline(xintercept = stats::median(s_len), color = "red3", linetype = 2)

  # Return ggplot objects
  return(ghist)
}



#' A list of primer sets
#'
#' @name primer_set
#' @docType data
#' @keywords data
"primer_set"
