#' @title Automatically download sequences and extract potential amplicons
#' @description \code{doamp_auto} Automatically download sequences and extract potential amplicons
#' @importFrom utils write.table
#' @param search_query Character. Search query for Entrez
#' @param F_primer Character. Forward primer sequence. Degenerate base allowed.
#' @param R_primer Character. Reverse primer sequence. Degenerate base allowed.
#' @param n_retmax Numeric. The maximum number of sequences collected from Entrez
#' @param n_retidmax Numeric. The maximum number of IDs collected from Entrez. Default is 10. If random_sampling = TRUE, n_retidmax IDs will be collected and then n_retmax sequences will be downloaded.
#' @param n_mismatch Numeric. The maximum number of primer-template mismatches allowed
#' @param random_sampling Logical. If TURE, n_retidmax IDs collected, and then n_retmax sequences are randomly collected
#' @param random_sampling_seed Numeric. Random number seed for random_sampling
#' @param output_dir Character. Output directory name
#' @param save_parameter Logical. If TRUE, parameters used in the analysis saved
#' @param save_stat Logical. If TRUE, summary of the analysis saved
#' @param overwrite_output_dir Logical. If TRUE, overwrite the contents of output directory.
#' @export

doamp_auto <- function (search_query,
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
                        overwrite_output_dir = FALSE) {
  #-----------------------------------------------------------------#
  # Create output directory
  #-----------------------------------------------------------------#
  time_start <- proc.time() # Measure elapsed time
  # Check and create output directory
  if(overwrite_output_dir) {
    dir.create(output_dir, showWarnings = FALSE)
    if(length(list.files(sprintf("%s", output_dir))) > 0) {
      previous_files <- dir(path = output_dir, pattern = "*.*")
      file.remove(sprintf("%s/%s", output_dir, previous_files))
    }
  } else {
    if(dir.exists(output_dir)) {
      stop("Output directory already exists")
    } else {
      dir.create(output_dir)
    }
  }
  #-----------------------------------------------------------------#
  # Search data in Entrez
  #-----------------------------------------------------------------#
  if (!random_sampling) n_retidmax <- n_retmax
  rentrez_search <- rentrez::entrez_search(db = "nucleotide", term = search_query, retmax = n_retidmax, use_history = TRUE)
  # Random sampling from rentrez_search w/ n_retidmax IDs
  if (random_sampling & length(rentrez_search$ids) >= n_retmax) {
    set.seed(random_sampling_seed)
    retids <- sample(rentrez_search$ids, n_retmax, replace = FALSE)
  } else {
    retids <- rentrez_search$ids
  }

  #-----------------------------------------------------------------#
  # Retrieve using entrez_fetch() function (rettype = "fasta")
  # Save as FASTA file
  #-----------------------------------------------------------------#
  entrez_fasta <- rentrez::entrez_fetch(db = "nucleotide", id = retids,
                                        rettype="fasta", parsed = FALSE)
  write.table(entrez_fasta, file = sprintf("%s/entrez_fasta.fa", output_dir),
              quote = FALSE, col.names = FALSE, row.names = FALSE)

  #-----------------------------------------------------------------#
  # Compile the sequence file using seqkit in Shell
  #-----------------------------------------------------------------#
  shell_command1 <- sprintf("seqkit seq -w 0 %s/entrez_fasta.fa > %s/download.fa",
                            output_dir, output_dir)
  if (Sys.info()["sysname"] == "Windows") {
    shell(shell_command1)
  } else {
    system(shell_command1)
  }
  file.remove(sprintf("%s/entrez_fasta.fa", output_dir))

  #-----------------------------------------------------------------#
  # Generate degenerate primer list
  #-----------------------------------------------------------------#
  # Check degenerated bases
  deg_in_F <- !all(strsplit(F_primer, split = NULL)[[1]] %in% c("A", "T", "G", "C"))
  deg_in_R <- !all(strsplit(R_primer, split = NULL)[[1]] %in% c("A", "T", "G", "C"))

  if((deg_in_F | deg_in_R)) {
    # Expand degenerate primers
    if(deg_in_F) expanded_F_primer <- expand_degenerate_primer(F_primer) else expanded_F_primer <- F_primer
    if(deg_in_R) expanded_R_primer <- expand_degenerate_primer(R_primer) else expanded_R_primer <- R_primer
    # Generate primer list
    expanded_primer_list <- expand.grid(expanded_F_primer, expanded_R_primer)
    rownames(expanded_primer_list) <- sprintf("p%s", 1:nrow(expanded_primer_list))
    write.table(expanded_primer_list,
                file = sprintf("%s/expanded_primer_list.tsv", output_dir),
                col.names = FALSE, row.names = TRUE, sep = "\t", quote = FALSE)

    #-----------------------------------------------------------------#
    # Extract the barcoding region using seqkit amplicon
    #-----------------------------------------------------------------#
    shell_command2 <- sprintf("seqkit amplicon -p %s/expanded_primer_list.tsv -m %s %s/download.fa | seqkit rmdup -w 0 > %s/amplified.fa",
                              output_dir,
                              n_mismatch,
                              output_dir, output_dir)
    if (Sys.info()["sysname"] == "Windows") shell(shell_command2) else system(shell_command2)
  } else {
    #-----------------------------------------------------------------#
    # Extract the barcoding region using seqkit amplicon
    #-----------------------------------------------------------------#
    shell_command2 <- sprintf("seqkit amplicon -F %s -R %s -w 0 -m %s %s/download.fa > %s/amplified.fa",
                              F_primer, R_primer,
                              n_mismatch,
                              output_dir, output_dir)
    if (Sys.info()["sysname"] == "Windows") shell(shell_command2) else system(shell_command2)
  }

  #-----------------------------------------------------------------#
  # Save query and stats
  #-----------------------------------------------------------------#
  # Collect parameters
  parameter_list <- matrix(c(paste("NCBI search query:", search_query),
                             paste("Forward primer:", F_primer),
                             paste("Reverse primer:", R_primer),
                             paste("N of retrieved sequences:", length(retids)),
                             paste("Maximum N of retrieved sequences:", n_retmax),
                             paste("N of retrieved IDs:", length(rentrez_search$ids)),
                             paste("Maximum N of retrieved IDs:", n_retidmax),
                             paste("Maximum N of mismatches:", n_mismatch),
                             paste("Random sampling:", random_sampling),
                             paste("Random sampling seed:", random_sampling_seed),
                             paste("Output directory:", output_dir)),
                           ncol = 1)
  if (save_parameter) {
    write.table(parameter_list, sprintf("%s/parameter_list.txt", output_dir),
                quote = FALSE, col.names = FALSE, row.names = FALSE)
  }

  if (Sys.info()["sysname"] == "Windows") {
    if (save_stat) {
      setwd(sprintf("%s", output_dir))
      shell("seqkit stats -T amplified.fa download.fa > stat.tsv")
      shell("more stat.tsv"); setwd("..")
    } else {
      setwd(sprintf("%s", output_dir))
      shell("seqkit stats -T amplified.fa download.fa"); setwd("..")
    }
  } else {
    if (save_stat) {
      system(sprintf("seqkit stats -T %s/*.fa > %s/stat.tsv", output_dir, output_dir))
      system(sprintf("cat %s/stat.tsv", output_dir))
    } else {
      system(sprintf("seqkit stats -T %s/*.fa", output_dir))
    }
  }

  # Generate output message
  time_elapsed <- (proc.time() - time_start)[3]
  message(paste("\nNCBI search query:", search_query))
  message(paste("Primer set:", F_primer, "-", R_primer, "\n"))
  message(paste(round(time_elapsed, 2), "sec elapsed for downloading and extracting sequences\n"))
}


#' @title Automatically extract potential amplicons from target sequences
#' @description \code{doamp_custom} Automatically extract potential amplicons from target sequences
#' @importFrom utils write.table
#' @param target_fasta Character. Path to a FASTA file that contains target sequences
#' @param F_primer Character. Forward primer sequence. Degenerate base allowed.
#' @param R_primer Character. Reverse primer sequence. Degenerate base allowed.
#' @param n_mismatch Numeric. The maximum number of primer-template mismatches allowed
#' @param output_dir Character. Output directory name
#' @param save_parameter Logical. If TRUE, parameters used in the analysis saved
#' @param save_stat Logical. If TRUE, summary of the analysis saved
#' @param overwrite_output_dir Logical. If TRUE, overwrite the contents of output directory.
#' @export

doamp_custom <- function (target_fasta,
                          F_primer,
                          R_primer,
                          n_mismatch = 0,
                          output_dir = "rDoAMP_Out",
                          save_parameter = TRUE,
                          save_stat = TRUE,
                          overwrite_output_dir = FALSE) {
  #-----------------------------------------------------------------#
  # Creat output directory
  #-----------------------------------------------------------------#
  time_start <- proc.time() # Measure elapsed time

  # Check target.fasta file
  if(!file.exists(target_fasta)) stop("Your FASTA file does not exit!")

  # Check and create output directory
  if(overwrite_output_dir) {
    dir.create(output_dir, showWarnings = FALSE)
      if(length(list.files(sprintf("%s", output_dir))) > 0) {
        previous_files <- dir(path = output_dir, pattern = "*.*")
        if(target_fasta %in% sprintf("%s/%s", output_dir, previous_files)) {
          stop("Your target FASTA file is in a output directory that will be overwritten! Change the FASTA file location.")
        }
        file.remove(sprintf("%s/%s", output_dir, previous_files))
      }
  } else {
    if(dir.exists(output_dir)) {
      stop("Output directory already exists")
    } else {
      dir.create(output_dir)
    }
  }
  # Check target.fasta file
  if(!file.exists(target_fasta)) stop("Your FASTA file does not exit!")
  # Copy target fasta to output_dir
  file.copy(target_fasta, sprintf("%s/custom_db0.fa", output_dir))

  #-----------------------------------------------------------------#
  # Compile the sequence file using seqkit in Shell
  #-----------------------------------------------------------------#
  shell_command1 <- sprintf("seqkit seq -w 0 %s/custom_db0.fa > %s/custom_db.fa",
                            output_dir, output_dir)
  if (Sys.info()["sysname"] == "Windows") {
    shell(shell_command1)
  } else {
    system(shell_command1)
  }
  file.remove(sprintf("%s/custom_db0.fa", output_dir))

  #-----------------------------------------------------------------#
  # Generate degenerate primer list
  #-----------------------------------------------------------------#
  # Check degenerated bases
  deg_in_F <- !all(strsplit(F_primer, split = NULL)[[1]] %in% c("A", "T", "G", "C"))
  deg_in_R <- !all(strsplit(R_primer, split = NULL)[[1]] %in% c("A", "T", "G", "C"))

  if((deg_in_F | deg_in_R)){
    # Expand degenerate primers
    if(deg_in_F) expanded_F_primer <- expand_degenerate_primer(F_primer) else expanded_F_primer <- F_primer
    if(deg_in_R) expanded_R_primer <- expand_degenerate_primer(R_primer) else expanded_R_primer <- R_primer
    # Generate primer list
    expanded_primer_list <- expand.grid(expanded_F_primer, expanded_R_primer)
    rownames(expanded_primer_list) <- sprintf("p%s", 1:nrow(expanded_primer_list))
    write.table(expanded_primer_list,
                file = sprintf("%s/expanded_primer_list.tsv", output_dir),
                col.names = FALSE, row.names = TRUE, sep = "\t", quote = FALSE)

    #-----------------------------------------------------------------#
    # Extract the barcoding region using seqkit amplicon
    #-----------------------------------------------------------------#
    shell_command2 <- sprintf("seqkit amplicon -p %s/expanded_primer_list.tsv -m %s %s/custom_db.fa | seqkit rmdup -w 0 > %s/amplified.fa",
                              output_dir,
                              n_mismatch,
                              output_dir, output_dir)
    if (Sys.info()["sysname"] == "Windows") shell(shell_command2) else system(shell_command2)
  } else {
    #-----------------------------------------------------------------#
    # Extract the barcoding region using seqkit amplicon
    #-----------------------------------------------------------------#
    shell_command2 <- sprintf("seqkit amplicon -F %s -R %s -w 0 -m %s %s/custom_db.fa > %s/amplified.fa",
                              F_primer, R_primer,
                              n_mismatch,
                              output_dir, output_dir)
    if (Sys.info()["sysname"] == "Windows") shell(shell_command2) else system(shell_command2)
  }

  #-----------------------------------------------------------------#
  # Save query and stats
  #-----------------------------------------------------------------#
  # Collect parameters
  parameter_list <- matrix(c(paste("Target FASTA file:", target_fasta),
                             paste("Forward primer:", F_primer),
                             paste("Reverse primer:", R_primer),
                             paste("N of maximum mismatch:", n_mismatch),
                             paste("Output directory:", output_dir)),
                           ncol = 1)
  if (save_parameter) {
    write.table(parameter_list, sprintf("%s/parameter_list.txt", output_dir),
                quote = FALSE, col.names = FALSE, row.names = FALSE)
  }

  if (Sys.info()["sysname"] == "Windows") {
    if (save_stat) {
      setwd(sprintf("%s", output_dir))
      shell("seqkit stats -T amplified.fa custom_db.fa > stat.tsv")
      shell("more stat.tsv"); setwd("..")
    } else {
      setwd(sprintf("%s", output_dir))
      shell("seqkit stats -T amplified.fa custom_db.fa"); setwd("..")
    }
  } else {
    if (save_stat) {
      system(sprintf("seqkit stats -T %s/*.fa > %s/stat.tsv", output_dir, output_dir))
      system(sprintf("cat %s/stat.tsv", output_dir))
    } else {
      system(sprintf("seqkit stats -T %s/*.fa", output_dir))
    }
  }

  # Generate output message
  time_elapsed <- (proc.time() - time_start)[3]
  message(paste("\nTarget FASTA file:", target_fasta))
  message(paste("Primer set:", F_primer, "-", R_primer, "\n"))
  message(paste(round(time_elapsed, 2), "sec elapsed for extracting amplicons"))
}

