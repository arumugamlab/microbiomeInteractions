######## phyloseq object processing -------------------------------------------

################################
## Process one phyloseq object
################################

#' Process one phyloseq object for pairwise co-occurrence statistiss
#'
#' @param p              phyloseq object
#' @param .parallel      if TRUE, run aaply function in parallel, using aaply/foreach
#'
#' @return list(integer, dataframe) containing number of taxapairwise occurrence statistics
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
  
  rownames(x) <- paste0(t[,"Phylum"], "::", t[,"Species"])
  
  # Make sure that there isnt any invalid data. THis happens in DN-metaHIT in CRC-meta
  
  valid_samples = which(!is.na(colSums(x)))
  x <- x[ , valid_samples]
  
  # Convert into 0-1 presence/absence matrix, since that's all we need
  
  x[x>0] <- 1
  
  return(calculate_pairwise_1(x, .parallel = .parallel))
  #return(calculate_pairwise_2(x, x))
  
}


#############################
# Process a list of studies
#############################

#' Process a list of phyloseq objects for pairwise co-occurrence statistics
#'
#' @param study_list  list of study names
#' @param ps_list     list of phyloseq objects that must contain at least the 
#'                    studies in \code{study_list}
#' @param min_log2DEP minimum log2DEP needed for an interaction to be included
#' @param .parallel   if TRUE, run aaply function in parallel, using aaply/foreach
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
process_phyloseq_study_list <- function(study_list, ps_list, min_log2DEP = 2, .parallel = FALSE) {
  
  my_exclusion  <- list()
  my_dependency <- list()
  
  for (name in study_list) {
    cat(paste("Processing Set:", name, "\n"))
    pairwise_res <- process_phyloseq_for_pairwise_interactions(ps_list[[name]], .parallel = .parallel)
    my_exclusion[[name]]  <- get_co_exclusions(pairwise_res, min_abs_log2DEP = abs(min_log2DEP))
    my_dependency[[name]] <- get_co_dependencies(pairwise_res)
  }
  return(list(exclusion=my_exclusion, dependency=my_dependency))
}


################################
## Get prevalent taxa
################################

#' Subset a phyloseq object with prevalent taxa 
#'
#' @param p phyloseq object
#' @param min_prevalence minimum prevalence below which taxa will be removed
#'
#' @return phyloseq object
#' @export
#'
#' @examples
#' \dontrun{
#' get_prevalent_taxa_from_phyloseq(p, 0.01)
#' }
get_prevalent_taxa_from_phyloseq <- function(p, min_prevalence) {
  n_samples <- nsamples(p)
  prev_taxa <- names(which((rowSums(as(otu_table(p), "matrix") != 0)/n_samples) > min_prevalence))
  return(prune_taxa(prev_taxa, p))
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

