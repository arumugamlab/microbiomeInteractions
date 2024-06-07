######## Co-occurrence calculations -------------------------------------------


################################
## Estimate cooccurrences from 1 dataframe
################################

#' Calculate pairwise co-occurrence relationship between entities in one dataframe.
#'
#' @param x          dataframe with entities in rows and sample in columns
#' @param .parallel  if TRUE, run analysis in parallel, using foreach and \%dopar\%
#'
#' @return dataframe containing pairwise occurrence statistics
#' @export
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom bit bitwhich
#'
#' @examples
#' \dontrun{
#' asvs <- as.data.frame(otu_table(p))
#' calculate_pairwise_1(asvs, asvs)
#' }
calculate_pairwise_1 <- function(x, shrink = FALSE, .parallel = FALSE) {

  # Estimate p_i

  n_features <- dim(x)[1]
  n_samples  <- dim(x)[2]
  p_i        <- rowSums(x>0) / n_samples

  # Start comparisons
  # NOTE / WARNING: I switched to bitwhich objects from "bit" package for speed. This affects how parallelization works.
  #                 If encoding happens in the main thread, worker threads do not interpret the encodings the same.
  #                 This leads to errors in computing overlaps in i_and_j. Therefore, I let the worker threads do their own
  #                 encoding, even though they all encode the full matrix repeatedly. But this keeps it accurate, so worth doing.
  #                 Even with this overhead, this is extremely fast, so double worth it!

  # Count the number of times two features co-occur in the samples
  if (.parallel) {
    i_and_j <- foreach (f1=1:(n_features-1), .combine = "c", .packages = c("bit")) %dopar%
      {
        # Create lists of bitwhich
        # For each feature in x, one bitwhich is created where values above 0 are TRUE and others are FALSE

        b1 <- bitwhich(n_samples, which(x[f1,]>0));
        y_bitwhich <- apply(X=x, MARGIN=1, FUN=function(row_data, n_items) { return(bitwhich(n_items, which(row_data > 0)))}, simplify=FALSE, n_items=n_samples);
        sapply((f1+1):n_features, function(f2) {sum(b1 & y_bitwhich[[f2]])})
      }
  } else {
    i_and_j <- foreach (f1=1:(n_features-1), .combine = "c", .packages = c("bit")) %do%
      {
        b1 <- bitwhich(n_samples, which(x[f1,]>0));
        y_bitwhich <- apply(X=x, MARGIN=1, FUN=function(row_data, n_items) { return(bitwhich(n_items, which(row_data > 0)))}, simplify=FALSE, n_items=n_samples);
        sapply((f1+1):n_features, function(f2) {sum(b1 & y_bitwhich[[f2]])})
      }
  }
  # Convert counts to probabilities
  i_and_j <- i_and_j / n_samples


  # Prepare results

  tmp_names <- rownames(x)

  # Make all combinations of species, and then estimate their co-occurrence probabilities
  combinations <- t(combn(n_features, m = 2))

  df <- data.frame(i = combinations[ ,2],
                   j = combinations[ ,1],
                   p_i_AND_j=i_and_j)
  df$s_i <- tmp_names[df$i]
  df$s_j <- tmp_names[df$j]
  df$p_i <- p_i[df$i]
  df$p_j <- p_i[df$j]
  df <- df %>%
    mutate(p_i_TIMES_p_j = p_i * p_j,
           p_i_OR_j      = p_i + p_j - p_i_AND_j,
           log2DEP       = log2(p_i_AND_j) - log2(p_i_TIMES_p_j)
    ) %>%
    select(i, s_i, j, s_j, p_i, p_j, p_i_TIMES_p_j, p_i_AND_j, p_i_OR_j, log2DEP) %>%
    arrange(log2DEP, s_i, s_j)

  return(df)
}



################################
## Estimate cooccurrences from 2 dataframes
################################

#' Calculate pairwise co-occurrence relationship between entities in two dataframes.
#'
#' @param x    1st dataframe with entities in rows and sample in columns
#' @param y    2nd dataframe with entities in rows and sample in columns
#' @param .parallel if TRUE, run analysis in parallel, using foreach and \%dopar\%
#'
#' @return dataframe containing pairwise occurrence statistics
#' @export
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom bit bitwhich
#'
#' @examples
#' \dontrun{
#' asvs <- as.data.frame(otu_table(p))
#' calculate_pairwise_2(asvs, asvs)
#'
#' # The above call is the same as running the single data frame function
#' # below, but will take twice as long.
#'
#' calculate_pairwise_1(asvs)
#'
#' # You can also run it against a function profile, let us say kegg_modules
#' calculate_pairwise_2(asvs, kegg_modules)
#' }
calculate_pairwise_2 <- function(x, y = NULL, shrink = FALSE, .parallel = FALSE) {

  self_comparison = FALSE
  if (is.null(y)) {
    self_comparison = TRUE
    y = x
  }

  # Make sure the 2 data frames have same samples

  if ( (dim(x)[2] != dim(y)[2]) || !identical(colnames(x), colnames(y)) ) {
    stop("Function \"calculate_pairwise_2\" needs two data frames with same samples in columns. Please fix it.",
         call. = FALSE)
  }

  # Estimate p_i

  n_feature1 <- dim(x)[1]
  n_feature2 <- dim(y)[1]
  n_samples  <- dim(x)[2]

  p_i        <- rowSums(x>0) / n_samples
  p_j        <- rowSums(y>0) / n_samples

  # Start comparisons
  # NOTE / WARNING: I switched to bitwhich objects from "bit" package for speed. This affects how parallelization works.
  #                 If encoding happens in the main thread, worker threads do not interpret the encodings the same.
  #                 This leads to errors in computing overlaps in i_and_j. Therefore, I let the worker threads do their own
  #                 encoding, even though they all encode the full matrix repeatedly. But this keeps it accurate, so worth doing.
  #                 Even with this overhead, this is extremely fast, so double worth it!


  # Count the number of times two features co-occur in the samples

  # For x vs y, we go i=1:n_feature1 and j=1:n_feature2
  loop1_max = n_feature1
  loop2_min = rep(1, loop1_max)
  # For x vs x, we go i=1:(n_feature1-1) and j=(i+1):n_feature1
  if (self_comparison == TRUE) {
    loop1_max = n_feature1 - 1
    loop2_min = seq(from = 2, to = n_feature1, by = 1)
  }
  if (.parallel) {
    i_and_j <- foreach (f1=1:loop1_max, .combine = "c", .packages = c("bit")) %dopar%
        {
          # Create lists of bitwhich
          # For each feature in x, one bitwhich is created where values above 0 are TRUE and others are FALSE

          b1 <- bitwhich(n_samples, which(x[f1,]>0));
          y_bitwhich <- apply(X=y, MARGIN=1, FUN=function(row_data, n_items) { return(bitwhich(n_items, which(row_data > 0)))}, simplify=FALSE, n_items=n_samples);
          sapply(loop2_min[f1]:n_feature2, function(f2) {sum(b1 & y_bitwhich[[f2]])})
        }
  } else {
    i_and_j <- foreach (f1=1:loop1_max, .combine = "c", .packages = c("bit")) %do%
        {
          b1 <- bitwhich(n_samples, which(x[f1,]>0));
          y_bitwhich <- apply(X=y, MARGIN=1, FUN=function(row_data, n_items) { return(bitwhich(n_items, which(row_data > 0)))}, simplify=FALSE, n_items=n_samples);
          sapply(loop2_min[f1]:n_feature2, function(f2) {sum(b1 & y_bitwhich[[f2]])})
        }
  }
  # Convert counts to probabilities
  i_and_j <- i_and_j / n_samples

  # Prepare results

  names_i <- rownames(x)
  names_j <- rownames(y)

  # Make all combinations of species, and then estimate their co-occurrence probabilities

  if (self_comparison == TRUE) {
    combinations <- t(combn(n_feature1, m = 2))
  } else {
    # expand.grid: first factor varies fastest.
    # calculation above: first factor is stable and 2nd factor varies 1:n
    # therefore, we do some switching around that might not sound logical.

    ilist <- seq_len(n_feature1)
    jlist <- seq_len(n_feature2)
    combinations <- as.matrix(expand.grid(jlist, ilist))
  }

  df <- data.frame(i = combinations[ ,2],
                   j = combinations[ ,1],
                   p_i_AND_j = i_and_j)
  df$p_i <- p_i[df$i]
  df$p_j <- p_j[df$j]
  df$s_i <- names_i[df$i]
  df$s_j <- names_j[df$j]
  df <- df %>%
        mutate(p_i_TIMES_p_j = p_i * p_j,
               p_i_OR_j      = p_i + p_j - p_i_AND_j,
               log2DEP       = log2(p_i_AND_j) - log2(p_i_TIMES_p_j)
               ) %>%
        select(i, s_i, j, s_j, p_i, p_j, p_i_TIMES_p_j, p_i_AND_j, p_i_OR_j, log2DEP) %>%
        arrange(log2DEP, s_i, s_j)

  if (isTRUE(shrink)) {
    df <- df %>%
            filter(log2DEP < -5 | log2DEP > 1)
  }

  return(df)
}


######## Co-dependences -------------------------------------------------------

#' Prioritize the top co-exclusions after rescoring absolute co-exclusions
#'
#' @param study_list        list of study names
#' @param exclusion_df_list list of dataframes with co-exclusion data
#' @param max_log2DEP       minimum value of log2DEP to keep a co-exclusion
#'                          in the returned dataframe. Note: valid
#'                          co-exclusions will be negative, so send a negative
#'                          number.
#'
#' @return data frame with co-exclusion data
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
#'
#' top_exc <- get_top_exclusions(name_list, exclusion)
#' }
get_top_dependencies <- function(study_list, pairwise_df_list, min_log2DEP = 2) {

  top_dependencies <- list()

  for (study in study_list) {
    df <- pairwise_df_list[[study]]

    top <- df %>%
      filter(log2DEP >= min_log2DEP) %>%
      mutate(p_i_OR_j = p_i + p_j - p_i_AND_j) %>%
      mutate(log2DEP = log2DEP + log10(p_i_OR_j)) %>%
      select(s_i, s_j, p_i, p_j, log2DEP, p_i_AND_j, p_i_OR_j, p_i_TIMES_p_j) %>%
      arrange(desc(log2DEP))

    top_dependencies[[study]] <- top
  }

  return(top_dependencies)
}


#' Obtain a list of co-exclusions that occur in multiple datasets
#'
#' @param study_list        list of study names
#' @param exclusion_df_list list of dataframes with exclusion info
#' @param min_exclusions    minimum number of studies where co-exclusion must occur
#' @param duplicate_pairs   whether to duplicate pairs by swapping *i* and *j*
#'                          labels of species. This helps in counting interactions
#'                          easily without checking both *s_i* and *s_j* columns
#'                          in the resulting dataframe.
#' @param exclusion_prob    background probability of exclusion, estimated from
#'                          datasets as number of co-exclusions among all pairwise
#'                          combinations. Default is 0.3, estimated from YuJ_2015,
#'                          FengQ_2015, ZellerG_2014.
#'
#' @return dataframe with overlapping co-exclusions.
#' @export
#'
#' @examples
#' \dontrun{
#' crc_studies = c("YuJ_2015", "ZellerG_2014", "FengQ_2015", "VogtmannE_2016")
#' crc_overlap = get_overlapping_exclusions(crc_studies, top_exc, min_exclusions = 3)
#' }
get_overlapping_dependencies <- function(study_list, dependencies_df_list, dependence_prob = 0.3, min_log2DEP = 1, min_dependencies = 2, duplicate_pairs = FALSE) {

  n <- length(study_list)
  df_list <- list()

  for (i in 1:n) {
    tmp_list1 <- dependencies_df_list[[study_list[i]]] %>%
      select(s_i, s_j, p_i, p_j, p_i_AND_j, log2DEP)
    if (duplicate_pairs == TRUE) {
      tmp_list2 <- tmp_list1 %>%
        select(s_j, s_i, p_j, p_i, p_i_AND_j, log2DEP)
      colnames(tmp_list2) <- colnames(tmp_list1)
      df_list[[i]] <- rbind(tmp_list1, tmp_list2)
    } else {
      df_list[[i]] <- tmp_list1
    }
  }

  # Collect overlapping entries of s_i and s_j
  overlap <- Reduce(function(x, y) dplyr::full_join(x, y, by=c("s_i", "s_j")), df_list, accumulate=FALSE)
  #overlap %>% dplyr::mutate(combined=apply(select(.,3,4),1,sum))

  # Set colnames to be more userfriendly
  # Note: This uses c(rbind(a,b)) construct, which cannot be used when length(a) <> length(b).
  colnames(overlap) <- c("species_i", "species_j", rbind(paste0(study_list, ".p_i"), paste0(study_list, ".p_j"), paste0(study_list, ".p_ij"), paste0(study_list, ".log2Dep")))

  offset <- 2 # represents s_i and s_j
  n_cols <- dim(df_list[[1]])[2] # columns coming from each result set
  n_cols <- n_cols - offset # since we do a join, common columns (offset) are only in the beginning

  # Calculate geo_means for s_i and s_j across studies

  # Version 1
  #    apply(., 1, function(x) {prod(x, na.rm=TRUE)}) %>%
  #    pracma::nthroot(n)
  # Version 2
  #    apply(., 1, function(x) {gm_mean(x, na.rm = TRUE, zero.propagate = TRUE)})

  indices_i <- offset+ ((1:n)-1)*n_cols + 1
  overlap$geo_mean_i <- overlap %>%
    select(all_of(indices_i)) %>%
    apply(., 1, function(x) {
      y = x[!is.na(x)];
      exp(sum(log(y))/length(y))
    })

  indices_j <- offset+ ((1:n)-1)*n_cols + 2
  overlap$geo_mean_j <- overlap %>%
    select(all_of(indices_j)) %>%
    apply(., 1, function(x) {
      y = x[!is.na(x)];
      exp(sum(log(y))/length(y))
    })

  # Calculate mean for p(s_i OR s_j) across studies

  indices_ij <- offset+ ((1:n)-1)*n_cols + 3
  overlap$geo_mean_ij <- overlap %>%
    select(all_of(indices_ij)) %>%
    apply(., 1, function(x) {
      y = x[!is.na(x)];
      exp(sum(log(y))/length(y))
    })

  # Calculate score for p(s_i & s_j) across studies

  indices_dep <- offset+ ((1:n)-1)*n_cols + 4

  x <- overlap %>%
    select(all_of(indices_dep)) %>%
    apply(., 1, function(x) {
      c(
        sum(!is.na(x)),
        sum(!is.na(x) & x > min_log2DEP),
        sum(x[x > min_log2DEP], na.rm=TRUE)
      )})
  overlap$n_observations <- x[1, ]
  overlap$n_dependencies <- x[2, ]
  overlap$mean_log2DEP   <- x[3, ]

  #########################################################################
  # Version 2
  # Calculate a score, which is the sum of log2DEP and log(dbinom)
  #########################################################################
  overlap <- overlap %>%
    mutate(
      log10_p_binomial = log10(pbinom(n_dependencies-1, n_observations, dependence_prob, lower.tail = FALSE)),
      mean_log2DEP     = mean_log2DEP/n_dependencies,
      score            = mean_log2DEP - log10_p_binomial
    )

  overlap <- overlap %>%
    filter(n_dependencies >= min_dependencies) %>%
    select(-ends_with("log2Dep", ignore.case = FALSE)) %>%
    select(species_i, species_j, score, n_dependencies, n_observations, log10_p_binomial, mean_log2DEP, geo_mean_ij, geo_mean_i, geo_mean_j, everything()) %>%
    arrange(desc(score), desc(n_dependencies), species_i, species_j)

  return(overlap)
}

######## Co-exclusions --------------------------------------------------------

#############################
# Summarize exclusions
#############################

#' Prioritize the top co-exclusions after rescoring absolute co-exclusions
#'
#' @param study_list        list of study names
#' @param exclusion_df_list list of dataframes with co-exclusion data
#' @param max_log2DEP       minimum value of log2DEP to keep a co-exclusion
#'                          in the returned dataframe. Note: valid
#'                          co-exclusions will be negative, so send a negative
#'                          number.
#'
#' @return data frame with co-exclusion data
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
#'
#' top_exc <- get_top_exclusions(name_list, exclusion)
#' }
get_top_exclusions <- function(study_list, exclusion_df_list, max_log2DEP = -2) {

  top_exclusions <- list()

  for (study in study_list) {
    df <- exclusion_df_list[[study]]

    # Make a pseudocount for p_ij when it is 0. I chose 2^-15 = (1/32768)
    # Then recalculate log2DEP
    p_ij <- df$p_i_AND_j
    p_ij[p_ij == 0] <- 2**(-15)
    df$log2DEP <- log2(p_ij)-log2(df$p_i_TIMES_p_j)

    top <- df %>%
            filter(log2DEP < max_log2DEP) %>%
            arrange(log2DEP, desc(p_i_TIMES_p_j)) %>%
            select(s_i, s_j, p_i, p_j, p_i_AND_j, log2DEP)

    top_exclusions[[study]] <- top
  }

  return(top_exclusions)
}


# Helper to get geometric mean
# Source: Paul McMurdie
# https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in

gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

#' Obtain a list of co-exclusions that occur in multiple datasets
#'
#' @param study_list        list of study names
#' @param exclusion_df_list list of dataframes with exclusion info
#' @param min_exclusions    minimum number of studies where co-exclusion must occur
#' @param duplicate_pairs   whether to duplicate pairs by swapping *i* and *j*
#'                          labels of species. This helps in counting interactions
#'                          easily without checking both *s_i* and *s_j* columns
#'                          in the resulting dataframe.
#' @param exclusion_prob    background probability of exclusion, estimated from
#'                          datasets as number of co-exclusions among all pairwise
#'                          combinations. Default is 0.3, estimated from YuJ_2015,
#'                          FengQ_2015, ZellerG_2014.
#'
#' @return dataframe with overlapping co-exclusions.
#' @export
#'
#' @examples
#' \dontrun{
#' crc_studies = c("YuJ_2015", "ZellerG_2014", "FengQ_2015", "VogtmannE_2016")
#' crc_overlap = get_overlapping_exclusions(crc_studies, top_exc, min_exclusions = 3)
#' }
get_overlapping_exclusions <- function(study_list, exclusion_df_list, exclusion_prob = 0.3, min_exclusions = 2, max_log2DEP = -5, duplicate_pairs = FALSE) {

  n <- length(study_list)
  df_list <- list()

  for (i in 1:n) {
    tmp_list1 <- exclusion_df_list[[study_list[i]]] %>%
                    select(s_i, s_j, p_i, p_j, p_i_AND_j, log2DEP)
    if (duplicate_pairs == TRUE) {
      tmp_list2 <- tmp_list1 %>%
        select(s_j, s_i, p_j, p_i, p_i_AND_j, log2DEP)
      colnames(tmp_list2) <- colnames(tmp_list1)
      df_list[[i]] <- rbind(tmp_list1, tmp_list2)
    } else {
      df_list[[i]] <- tmp_list1
    }
  }

  # Collect overlapping entries of s_i and s_j
  overlap <- Reduce(function(x, y) dplyr::full_join(x, y, by=c("s_i", "s_j")), df_list, accumulate=FALSE)
  #overlap %>% dplyr::mutate(combined=apply(select(.,3,4),1,sum))

  # Set colnames to be more userfriendly
  # Note: This uses c(rbind(a,b)) construct, which cannot be used when length(a) <> length(b).
  colnames(overlap) <- c("species_i", "species_j", rbind(paste0(study_list, ".p_i"), paste0(study_list, ".p_j"), paste0(study_list, ".p_ij"), paste0(study_list, ".log2Dep")))

  offset <- 2 # represents s_i and s_j
  n_cols <- dim(df_list[[1]])[2] # columns coming from each result set
  n_cols <- n_cols - offset # since we do a join, common columns (offset) are only in the beginning

  # Calculate geo_means for s_i and s_j across studies

  # Version 1
  #    apply(., 1, function(x) {prod(x, na.rm=TRUE)}) %>%
  #    pracma::nthroot(n)
  # Version 2
  #    apply(., 1, function(x) {gm_mean(x, na.rm = TRUE, zero.propagate = TRUE)})

  indices_i <- offset+ ((1:n)-1)*n_cols + 1
  overlap$geo_mean_i <- overlap %>%
                          select(all_of(indices_i)) %>%
                          apply(., 1, function(x) {
                            y = x[!is.na(x)];
                            exp(sum(log(y))/length(y))
                          })

  indices_j <- offset+ ((1:n)-1)*n_cols + 2
  overlap$geo_mean_j <- overlap %>%
                          select(all_of(indices_j)) %>%
                          apply(., 1, function(x) {
                            y = x[!is.na(x)];
                            exp(sum(log(y))/length(y))
                          })

  # Calculate how many studies does each co-exclusion occur
  # NOTE: Do NOT replace the sum(!is.na(x)) with sum(x, na.rm=TRUE).
  #       It is not summing all non-NA elements, it is counting all
  #       non-NA elements.

  # Calculate mean for p(s_i & s_j) across studies

  indices_ij <- offset+ ((1:n)-1)*n_cols + 3
  overlap$mean_ij <- overlap %>%
                      select(all_of(indices_ij)) %>%
                      apply(., 1, function(x) {mean(x, na.rm = TRUE)})

  # Calculate score for p(s_i & s_j) across studies

  indices_dep <- offset+ ((1:n)-1)*n_cols + 4

  #########################################################################
  # Version 1
  # Calculate a score, which is the geometric mean of all
  #########################################################################
  # overlap$score = overlap %>%
  #   select(indices_dep) %>%
  #   apply(., 1, function(x) {prod(x, na.rm=TRUE)})
  # overlap = overlap %>%
  #   mutate(score = pracma::nthroot(score, n_exclusions))

  #########################################################################
  # Version 2
  # Calculate a score, which is the sum of mean(log2DEP) and log(dbinom)
  #########################################################################
  #overlap$p_joint <- overlap %>%
  #                  select(all_of(indices_dep)) %>%
  #                  apply(., 1, function(x) {sum(log10(exclusion_prob_list[!is.na(x)])) + sum(log10(1-exclusion_prob_list[is.na(x)]))})

  x <- overlap %>%
    select(all_of(indices_dep)) %>%
    apply(., 1, function(x) {
      c(
        sum(!is.na(x)),
        sum(!is.na(x) & x < max_log2DEP),
        sum(x[x < max_log2DEP], na.rm=TRUE)
      )})
  overlap$n_observations <- x[1, ]
  overlap$n_exclusions   <- x[2, ]
  overlap$mean_log2DEP   <- x[3, ]

  overlap <- overlap %>%
              mutate(
                log10_p_binomial = log10(pbinom(n_exclusions-1, n_observations, exclusion_prob, lower.tail = FALSE)),
                mean_log2DEP     = mean_log2DEP/n_exclusions,
                score            = -(mean_log2DEP + log10_p_binomial),
                )
  overlap <- overlap %>%
              filter(n_exclusions >= min_exclusions) %>%
              select(-ends_with("log2Dep", ignore.case = FALSE)) %>%
              select(species_i, species_j, score, n_exclusions, n_observations, log10_p_binomial, mean_log2DEP, mean_ij, geo_mean_i, geo_mean_j, everything()) %>%
              arrange(desc(score), desc(n_exclusions), species_i, species_j)

  return(overlap)
}
