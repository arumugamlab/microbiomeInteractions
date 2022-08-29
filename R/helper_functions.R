################################
## Process two phyloseq objects
################################

#' Process two phyloseq objects for pairwise co-occurrence statistics
#'
#' @param p1             1st phyloseq object
#' @param p2             2nd phyloseq object
#' @param .parallel      if TRUE, run multithreading
#'
#' @return list(integer, dataframe) containing number of taxa and pairwise occurrence statistics
#' @export
#' 
#' @importFrom phyloseq tax_table otu_table prune_taxa
#'
#' @examples
#' \dontrun{
#' pairwise_res <- process_phyloseq_for_pairwise_interactions(ps_list[["ZellerG_2014"]])
#' }
process_two_phyloseqs_for_pairwise_interactions <- function(p1, p2, .parallel = FALSE) {

  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package \"phyloseq\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  # Sanity check
  
  if ( !identical(row.names(sample_data(p1)), row.names(sample_data(p2))) ) {
    stop("Function \"calculate_pairwise_2\" needs two data frames with same samples in columns. Please fix it.",
         call. = FALSE)
  }

  
  # Process the first phyloseq object
  
  # Prepare otu_table
  
  x <- as.data.frame(otu_table(p1))
  
  # CMD has species as rows and samples as columns.
  # calculate_pairwise() assumes this mode, so need to transpose() if that's not the case
  
  if (!taxa_are_rows(p1)) {
    x <- t(x)
  }
  
  # If this contains taxonomic features
  # Assign useful taxa names for x
  
  if (length(grep("tax_table", getslots.phyloseq(p1), fixed = TRUE)) > 0) {
    t <- tax_table(p1)
    rownames(x) <- paste0(t[,"phylum"], "::", t[,"species"])
  }
  
  # Process the 2nd phyloseq object
  
  # Prepare otu table
  
  y <- as.data.frame(otu_table(p2))

  # CMD has species as rows and samples as columns.
  # calculate_pairwise() assumes this mode, so need to transpose() if that's not the case
  
  if (!taxa_are_rows(p2)) {
    y <- t(y)
  }
  
  # If this contains taxonomic features
  # Assign useful taxa names for y
  
  if (length(grep("tax_table", getslots.phyloseq(p2), fixed = TRUE)) > 0) {
    t <- tax_table(p2)
    rownames(y) <- paste0(t[,"phylum"], "::", t[,"species"])
  }
  
  return(calculate_pairwise_2(x, y, .parallel = .parallel))

}


######## phyloseq object processing -------------------------------------------

################################
## Process one phyloseq object
################################

#' Process one phyloseq object for pairwise co-occurrence statistics
#'
#' @param p              phyloseq object
#' @param .parallel      if TRUE, run multithreading
#'
#' @return list(integer, dataframe) containing number of taxa and pairwise occurrence statistics
#' @export
#' 
#' @importFrom phyloseq tax_table otu_table prune_taxa
#'
#' @examples
#' \dontrun{
#' pairwise_res <- process_phyloseq_for_pairwise_interactions(ps_list[["ZellerG_2014"]])
#' }
process_phyloseq_for_pairwise_interactions <- function(p, .parallel = FALSE) {

  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package \"phyloseq\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  # Prepare otu_table
  
  t <- tax_table(p)
  x <- as.data.frame(otu_table(p))
  
  # CMD has species as rows and samples as columns.
  # calculate_pairwise() assumes this mode, so need to transpose() if that's not the case
  
  if (!taxa_are_rows(p)) {
    x <- t(x)
  }

  # Sort them so that comparison across datasets is easier

  x <- x[sort(rownames(x)), ]
  t <- t[sort(rownames(t)), ]
  
  # Assign useful taxa names for x
  
  rownames(x) <- paste0(t[,"phylum"], "::", t[,"species"])
  
  # Make sure that there isn't any invalid data. This happens in DN-metaHIT in CRC-meta
  
  valid_samples = which(!is.na(colSums(x)))
  x <- x[ , valid_samples]
  
  return(calculate_pairwise_1(x, .parallel = .parallel))
  #return(calculate_pairwise_2(x, x))
  
}


#############################
# Process a list of studies
#############################

#' Process a list of phyloseq objects for pairwise co-occurrence statistics
#'
#' @param study_list  list of study names
#' @param ps_list_1   1st list of phyloseq objects that must contain at least the 
#'                    studies in \code{study_list}
#' @param ps_list_2   Optional 2nd list of phyloseq objects that must contain at least the 
#'                    studies in \code{study_list}
#' @param min_log2DEP minimum log2DEP needed for an interaction to be included
#' @param .parallel   if TRUE, run multithreading
#'
#' @return list of 2 elements: containing list(exclusions) and list(dependencies)
#' @export
#'
#' @examples
#' \dontrun{
#' name_list <- c("ZellerG_2014", "YuJ_2015")
#' tmp_results <- process_phyloseq_study_list(name_list, ps_list)
#'
#' for (study in name_list) {
#'   exclusion[[study]]  <- tmp_results[[1]][[study]]
#'   equivalent[[study]] <- tmp_results[[2]][[study]]
#' }
#' }
process_phyloseq_study_list <- function(study_list, ps_list_1, ps_list_2 = NULL, min_log2DEP = 2, .parallel = FALSE) {
  
  my_exclusion  <- list()
  my_dependency <- list()
  
  for (name in study_list) {
    cat(paste("Processing Set:", name, "\n"))
    if (is.null(ps_list_2)) {
      pairwise_res <- process_phyloseq_for_pairwise_interactions(ps_list_1[[name]], .parallel = .parallel)
    } else {
      pairwise_res <- process_two_phyloseqs_for_pairwise_interactions(ps_list_1[[name]], ps_list_2[[name]], .parallel = .parallel)  
    }
    my_exclusion[[name]]  <- get_co_exclusions(pairwise_res, min_abs_log2DEP = abs(min_log2DEP))
    my_dependency[[name]] <- get_co_dependencies(pairwise_res)
  }
  return(list(exclusion=my_exclusion, dependency=my_dependency))
}


################################
## Get prevalent features
################################

#' Subset a phyloseq object with prevalent taxa 
#'
#' @param p phyloseq object
#' @param min_prevalence minimum prevalence below which taxa will be removed
#' @param min_occurrence minimum absolute occurrence below which feature will be removed (helps in smaller sample size)
#'
#' @return phyloseq object
#' @export
#'
#' @examples
#' \dontrun{
#' get_prevalent_taxa_from_phyloseq(p)
#' get_prevalent_taxa_from_phyloseq(p, 0.05, min_occurrence = 5)
#' }
get_prevalent_taxa_from_phyloseq <- function(p, min_prevalence = 0.01, min_occurrence = 1) {
  n_samples <- nsamples(p)
  min_occurrence <- max(min_prevalence * n_samples, min_occurrence)
  prev_taxa <- names(which((rowSums(as(otu_table(p), "matrix") != 0)) >= min_occurrence))
  return(prune_taxa(prev_taxa, p))
}


#' Subset a SummarizedExperiment object with prevalent features 
#'
#' @param se SummarizedExperiment object
#' @param min_prevalence minimum prevalence below which feature will be removed
#' @param min_occurrence minimum absolute occurrence below which feature will be removed (helps in smaller sample size)
#'
#' @return SummarizedExperiment object
#' @export
#'
#' @examples
#' \dontrun{
#' get_prevalent_feature_from_summarized_experiment(se)
#' get_prevalent_feature_from_summarized_experiment(se, 0.05, min_occurrence = 5, resource_name = "gene_families")
#' }
get_prevalent_features_from_summarized_experiment <- function(se, min_prevalence = 0.01, min_occurrence = 10, resource_name = "gene_families") {
  n_samples <- dim(se)[2]
  min_occurrence <- max(min_prevalence * n_samples, min_occurrence)
  prev_features <- which((rowSums(assays(se)[[resource_name]] != 0)) >= min_occurrence)
  return(se[prev_features, , drop = FALSE])
}

################################
## Remove invariant features
################################

#' Subset a phyloseq object to remove near-invariant taxa 
#'
#' @param p phyloseq object
#' @param max_prevalence maximum prevalence above which taxa will be removed
#'
#' @return phyloseq object
#' @export
#'
#' @examples
#' \dontrun{
#' remove_invariant_taxa_from_phyloseq(p, 0.99)
#' }
remove_invariant_taxa_from_phyloseq <- function(p, max_prevalence = 1.0) {
  n_samples <- nsamples(p)
  max_occurrence <- max_prevalence * n_samples
  prev_taxa <- names(which((rowSums(as(otu_table(p), "matrix") != 0)) < max_occurrence))
  return(prune_taxa(prev_taxa, p))
}

#' Subset a SummarizedExperiment object to remove near-invariant features 
#'
#' @param se SummarizedExperiment object
#' @param max_prevalence maximum prevalence above which features will be removed
#'
#' @return SummarizedExperiment object
#' @export
#'
#' @examples
#' \dontrun{
#' remove_invariant_features_from_summarized_experiment(se, 0.99)
#' }
remove_invariant_features_from_summarized_experiment <- function(se, max_prevalence = 1.0, resource_name = "gene_families") {
  n_samples <- dim(se)[2]
  max_occurrence <- max_prevalence * n_samples
  prev_features <- which((rowSums(assays(se)[[resource_name]] != 0)) < max_occurrence)
  return(se[prev_features, , drop = FALSE])
}

################################
## Get one rank
################################

#' Get abundance profiles of a single taxonomic rank from a MetaPhlan phyloseq object
#'
#' @param p        phyloseq object
#' @param tax_rank required taxonomic rank
#'
#' @return phyloseq object with subset of taxa belonging to the specified rank
#' @export
#'
#' @examples
#' \dontrun{
#' get_single_rank_metaphlan(p, "species")
#' }
get_single_rank_metaphlan <- function(p, tax_rank) {
  
  t <- tax_table(p)
  
  if (colnames(t)[8] == "Strain") {
    # This will be CMD or from MetaPhlan2
    taxa       <- rownames(t)
    tax_prefix <- paste0(substr(tax_rank, 1, 1), "__")
    tax_select <- taxa[startsWith(taxa, tax_prefix)]
    p <- prune_taxa(tax_select, p)
  }
  
  return(p)  
  
}

################################
# Plot interactions
################################

#' Plot interaction between two taxa as scatterplot
#'
#' @param ps_list list of phyloseq objects
#' @param t1 taxon 1
#' @param t2 taxon 2
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' 
#' t1="Blautia coccoides"
#' t2="Firmicutes bacterium CAG:95"
#'
#' plot_interaction_phyloseq(ps_list, t1, t2)
#' }
plot_interaction_phyloseq <- function(ps_list, t1, t2) {
  studies <- names(ps_list)
  df <- data.frame(double(), double(), character())
  for (study in studies) {
    p <- ps_list[[study]]
    t <- tax_table(p)
    if (length(grep(t1, row.names(t), fixed=TRUE)) == 0 |
        length(grep(t2, row.names(t), fixed=TRUE)) == 0) {
      next
    }
    ab <- otu_table(p)
    ab <- as.data.frame(t(ab[c(t1, t2), ]))
    ab <- ab[ab[ , 1]>0 | ab[ , 2]>0, ]
    ab[ab==0] <- 1e-4
    ab <- log10(ab)
    ab$study <- study
    df <- rbind(df, ab)
  }
  names(df) <- c(t1, t2, "study")
  g <- ggplot(df, aes(x=get(t1), y=get(t2))) + 
    geom_point() +
    xlab(t1) +
    ylab(t2) +
    facet_wrap(~study)
  return(g)
}
