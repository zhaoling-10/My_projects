# ======= src/script/stability2.R =======
suppressPackageStartupMessages({
  library(Matrix)
  # install.packages(c("Hmisc","mclust"))  # run once if missing
})

# Focus: overall ploidy effect (4x vs 2x), independent of within-group genotype differences.
# This version averages replicates within each genotype (nesting_mode = "average")
# and runs LOO + 80% subsampling diagnostics. No permutation null here.

run_stability <- function(
    data, input_parameter, input_file,
    n_subsamples = 200, seed = 123,
    nesting_mode = c("average","residuals"),   # default "average" = genotype means
    verbose = FALSE, progress = TRUE,
    fixed_threshold = TRUE,                    # reuse reference threshold in all refits
    topM_edges = 0,                            # 0 => use all edges; >0 => top-M by |r|
    corr_type = c("pearson","spearman"),       # keep consistent with config.R (default pearson)
    threshold_shift = 0,                       # 0 => exact SWIM threshold; negative => more lenient
    threshold_min   = 0.45                     # do not go below this floor
) {
  nesting_mode <- match.arg(nesting_mode)
  corr_type    <- match.arg(corr_type)
  
  message("==> Stability diagnostics: start (mode=", nesting_mode,
          ", corr=", corr_type, ", fixed_threshold=", fixed_threshold, ")")
  
  # --- Ensure figure directories exist (SWIM may write PDFs here) ---
  if (!is.null(input_parameter$path)) {
    fig_dirs <- c(
      file.path(input_parameter$path, "switch",  "figure"),
      file.path(input_parameter$path, "network", "figure")
    )
    for (d in fig_dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  
  # ----- 0) Load expression and labels (use the same files as main pipeline) -----
  expr <- read.delim(input_parameter$filename_data, header = TRUE, check.names = FALSE)
  if (!("gene_id" %in% colnames(expr))) {
    stop("Expected 'gene_id' column in filtered matrix. Please check ", input_parameter$filename_data)
  }
  rownames(expr) <- expr$gene_id
  expr$gene_id <- NULL
  expr <- as.matrix(data.matrix(expr))  # coerce to numeric
  
  ctrl_ids <- readLines(input_parameter$filename_CTRL, warn = FALSE)
  case_ids <- readLines(input_parameter$filename_CASE, warn = FALSE)
  ctrl_ids <- trimws(ctrl_ids); ctrl_ids <- ctrl_ids[nzchar(ctrl_ids)]
  case_ids <- trimws(case_ids); case_ids <- case_ids[nzchar(case_ids)]
  
  stopifnot(all(ctrl_ids %in% colnames(expr)))
  stopifnot(all(case_ids %in% colnames(expr)))
  
  sample_ids <- colnames(expr)
  ploidy <- factor(ifelse(sample_ids %in% case_ids, "4x",
                          ifelse(sample_ids %in% ctrl_ids, "2x", NA)),
                   levels = c("2x","4x"))
  if (any(is.na(ploidy))) {
    bad <- sample_ids[is.na(ploidy)]
    stop("Some samples are not in CASE/CTRL lists: ", paste(bad, collapse=", "))
  }
  
  # ----- infer genotype IDs from naming scheme (matches your current dataset) -----
  rep_num <- as.integer(sub(".*\\.rep(\\d+)$", "\\1", sample_ids))
  genotype <- character(length(sample_ids))
  is_dip <- grepl("^Pn\\.diploid", sample_ids)
  is_tet <- grepl("^Pn\\.tetraploid", sample_ids)
  
  genotype[is_dip & rep_num %in% 1:3] <- "2xSex306"
  genotype[is_dip & rep_num %in% 4:6] <- "2xSex22"
  genotype[is_tet & rep_num %in%  1:3]  <- "4xSex216"
  genotype[is_tet & rep_num %in%  4:6]  <- "4xApo30"
  genotype[is_tet & rep_num %in%  7:9]  <- "4xApo34"
  genotype[is_tet & rep_num %in% 10:12] <- "4xApo115"
  
  if (any(!nzchar(genotype))) {
    stop("Could not infer genotype for samples: ",
         paste(sample_ids[!nzchar(genotype)], collapse=", "))
  }
  genotype <- factor(genotype)
  
  # ----- control for genotype nesting -----
  # average = compute a genotype-level matrix (captures overall 4x vs 2x pattern)
  # residuals = regress out genotype (keeps replicate-level matrix)
  if (nesting_mode == "residuals") {
    expr_for_corr <- t(apply(expr, 1, function(y) {
      res <- resid(lm(y ~ genotype))
      res[is.na(res)] <- 0
      res
    }))
    colnames(expr_for_corr) <- colnames(expr)
    rownames(expr_for_corr) <- rownames(expr)
  } else { # "average"
    grp <- split(seq_along(genotype), genotype)
    expr_for_corr <- sapply(grp, function(idx) rowMeans(expr[, idx, drop = FALSE]))
    sample_ids <- colnames(expr_for_corr)                    # now genotype IDs
    ploidy <- factor(ifelse(grepl("^2x", sample_ids), "2x", "4x"), levels = c("2x","4x"))
  }
  
  # remove genes with zero variance
  keep_genes <- apply(expr_for_corr, 1, function(v) var(v) > 0)
  if (!all(keep_genes)) expr_for_corr <- expr_for_corr[keep_genes, , drop = FALSE]
  
  # ----- override CASE/CTRL per run (so SWIM reads temporary lists) -----
  override_case_ctrl <- function(case_now, ctrl_now) {
    tf_case <- tempfile("CASE_", fileext = ".txt")
    tf_ctrl <- tempfile("CTRL_", fileext = ".txt")
    writeLines(case_now, tf_case)
    writeLines(ctrl_now, tf_ctrl)
    input_file$filename_CASE <<- tf_case
    input_file$filename_CTRL <<- tf_ctrl
    filename_CASE <<- tf_case
    filename_CTRL <<- tf_ctrl
    invisible(TRUE)
  }
  
  # ----- helpers -----
  if (!requireNamespace("Hmisc", quietly = TRUE))
    stop("Please install 'Hmisc' (needed for rcorr)")
  if (!requireNamespace("mclust", quietly = TRUE))
    stop("Please install 'mclust' for ARI.")
  
  adj_to_edges <- function(adj) {
    if (!inherits(adj, "dgCMatrix")) adj <- Matrix(adj, sparse = TRUE)
    U <- which(triu(adj, 1) != 0, arr.ind = TRUE)
    genes <- rownames(adj); if (is.null(genes)) genes <- seq_len(nrow(adj))
    data.frame(a = genes[U[,1]], b = genes[U[,2]], stringsAsFactors = FALSE)
  }
  edge_key <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(character(0))
    paste(pmin(df$a, df$b), pmax(df$a, df$b), sep = "_")
  }
  topM_edge_keys <- function(data_input, M, thr) {
    rc <- Hmisc::rcorr(t(data_input), type = corr_type); R <- rc$r
    diag(R) <- 0; R[is.na(R)] <- 0
    A <- abs(R) >= thr
    if (M <= 0) {
      A[lower.tri(A)] <- t(A)[lower.tri(A)]
      idx <- which(A, arr.ind = TRUE)
      if (nrow(idx) == 0) return(character(0))
      genes <- rownames(R)
      return(paste(pmin(genes[idx[,1]], genes[idx[,2]]),
                   pmax(genes[idx[,1]], genes[idx[,2]]), sep = "_"))
    } else {
      Rabs <- abs(R); Rabs[!A] <- 0
      Rabs[lower.tri(Rabs, diag = TRUE)] <- 0
      nz <- which(Rabs > 0, arr.ind = TRUE)
      if (nrow(nz) == 0) return(character(0))
      vals <- Rabs[cbind(nz[,1], nz[,2])]
      ord <- order(vals, decreasing = TRUE)
      take <- seq_len(min(M, length(ord)))
      nz <- nz[ord[take], , drop = FALSE]
      genes <- rownames(R)
      return(paste(genes[nz[,1]], genes[nz[,2]], sep = "_"))
    }
  }
  jaccard <- function(A, B) {
    A <- unique(A); B <- unique(B)
    u <- union(A, B); iu <- length(u)
    if (iu == 0) return(NA_real_)
    length(intersect(A, B)) / iu
  }
  module_ari <- function(mfull, msub) {
    common <- intersect(names(mfull), names(msub))
    if (length(common) < 2) return(0)
    mclust::adjustedRandIndex(mfull[common], msub[common])
  }
  switch_recovery <- function(Sfull, Ssub) {
    if (length(Sfull) == 0) return(0)
    sum(Sfull %in% Ssub) / length(Sfull)
  }
  
  # ----- run_once: rebuild adjacency from (shifted) SWIM threshold; catch plotting errors in SwitchMining -----
  .thr_log <- numeric(0)
  run_once <- function(subset_samples, case_now, ctrl_now, thr_override = NULL) {
    override_case_ctrl(case_now, ctrl_now)
    data_input <- expr_for_corr[, subset_samples, drop = FALSE]
    rn <- rownames(data_input)
    
    network_local <- NetworkAnalysis(data_input, checkNetIntegrity = FALSE, screePlot = FALSE)
    switch_local  <- tryCatch(
      SwitchMining(),
      error = function(e) {
        if (verbose) message("SwitchMining() error (continuing): ", e$message)
        list(switch_genes = character(0))
      }
    )
    
    thr <- if (!is.null(thr_override)) thr_override else network_local$threshold_corr
    if (!is.numeric(thr) || is.na(thr)) thr <- 0.7
    thr <- max(threshold_min, thr + threshold_shift)
    
    if (verbose) message("Using threshold_corr = ", signif(thr, 3), " (", corr_type, ")")
    .thr_log <<- c(.thr_log, thr)
    
    # drop zero-variance genes in this subset
    keep_g <- apply(data_input, 1, function(v) var(v) > 0)
    if (!all(keep_g)) {
      data_input <- data_input[keep_g, , drop = FALSE]
      rn <- rownames(data_input)
    }
    
    # Correlation + adjacency (for edge Jaccard, independent of SWIM internals)
    rc <- Hmisc::rcorr(t(data_input), type = corr_type)
    R  <- rc$r
    if (is.null(R)) stop("rcorr returned NULL correlation matrix (check input).")
    diag(R) <- 0
    R[is.na(R)] <- 0
    A <- abs(R) >= thr
    A[lower.tri(A)] <- t(A)[lower.tri(A)]
    adj <- Matrix::Matrix(A, sparse = TRUE, dimnames = list(rn, rn))
    
    # modules from SWIM object (fallback to single-module if missing)
    modules <- NULL
    if (!is.null(network_local$idx)) {
      modules <- network_local$idx
    } else {
      for (nm in c("modules","module","clusters","cluster","communities")) {
        if (!is.null(network_local[[nm]])) { modules <- network_local[[nm]]; break }
      }
    }
    if (is.null(modules)) { modules <- rep(1L, nrow(adj)); names(modules) <- rn }
    else if (is.null(names(modules)) && length(modules) == nrow(adj)) names(modules) <- rn
    
    # switches from SWIM object (robust to field name variants)
    switches <- NULL
    for (nm in c("switch_genes","SwitchGenes","switch","switches")) {
      if (!is.null(switch_local[[nm]])) { switches <- unique(switch_local[[nm]]); break }
    }
    if (is.null(switches) && !is.null(switch_local$genes)) switches <- unique(switch_local$genes)
    if (is.null(switches)) switches <- character(0)
    
    list(adj = adj, modules = modules, switches = switches, thr = thr, data_input = data_input)
  }
  
  set.seed(seed)
  
  # ----- 4) reference run (build the genotype-level 4x vs 2x contrast) -----
  if (nesting_mode == "residuals") {
    case_now <- sample_ids[ploidy == "4x"]
    ctrl_now <- sample_ids[ploidy == "2x"]
  } else {
    case_now <- sample_ids[grepl("^4x", sample_ids)]
    ctrl_now <- sample_ids[grepl("^2x", sample_ids)]
  }
  ref <- run_once(subset_samples = sample_ids, case_now = case_now, ctrl_now = ctrl_now)
  ref_thr <- ref$thr
  ref_switches <- unique(ref$switches)
  message("Reference switch genes (count): ", length(ref_switches))
  if (verbose && length(ref_switches) > 0) message("Example switches: ", paste(head(ref_switches, 10), collapse = ", "))
  write.table(ref_switches, file = "Reference_Switches.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Reference edge set for Jaccard
  if (topM_edges > 0) {
    E_full_keyed <- topM_edge_keys(ref$data_input, M = topM_edges, thr = ref_thr)
  } else {
    E_full <- adj_to_edges(ref$adj)
    E_full_keyed <- edge_key(E_full)
  }
  
  # ----- 5) LOO + 80% subsampling -----
  idx_4x <- which(ploidy == "4x")
  idx_2x <- which(ploidy == "2x")
  make_subsample <- function() {
    sub_4x <- sample(idx_4x, size = ceiling(0.8 * length(idx_4x)))
    sub_2x <- sample(idx_2x, size = ceiling(0.8 * length(idx_2x)))
    sort(c(sub_4x, sub_2x))
  }
  
  rob_list <- vector("list", length = length(sample_ids) + n_subsamples)
  switches_list <- vector("list", length = length(sample_ids) + n_subsamples + 1) # +1 for ref
  switches_list[[1]] <- ref_switches
  k <- 0
  
  if (progress) { pb_loo <- txtProgressBar(min = 0, max = length(sample_ids), style = 3); on.exit({ if (exists("pb_loo")) close(pb_loo) }, add = TRUE) }
  message("-> LOO phase...")
  for (i in seq_along(sample_ids)) {
    set.seed(seed + i)
    keep <- setdiff(seq_along(sample_ids), i)
    subsamples <- sample_ids[keep]
    case_now <- subsamples[ploidy[keep] == "4x"]
    ctrl_now <- subsamples[ploidy[keep] == "2x"]
    
    fi <- run_once(
      subset_samples = subsamples,
      case_now = case_now, ctrl_now = ctrl_now,
      thr_override = if (fixed_threshold) ref_thr else NULL
    )
    
    # Jaccard
    if (topM_edges > 0) {
      Ei_keyed <- topM_edge_keys(fi$data_input, M = topM_edges, thr = if (fixed_threshold) ref_thr else fi$thr)
      jacc <- jaccard(E_full_keyed, Ei_keyed)
    } else {
      Ei <- adj_to_edges(fi$adj)
      jacc <- jaccard(E_full_keyed, edge_key(Ei))
    }
    
    k <- k + 1
    rob_list[[k]] <- data.frame(
      type = "LOO", run = i,
      edge_jaccard = jacc,
      module_ari   = module_ari(ref$modules, fi$modules),
      switch_recovery = switch_recovery(ref$switches, fi$switches),
      thr_used = fi$thr,
      n_switches = length(unique(fi$switches))
    )
    switches_list[[k + 1]] <- unique(fi$switches)
    
    if (progress) setTxtProgressBar(pb_loo, i)
  }
  
  if (progress) { pb_sub <- txtProgressBar(min = 0, max = n_subsamples, style = 3); on.exit({ if (exists("pb_sub")) close(pb_sub) }, add = TRUE) }
  message("-> 80% subsampling phase...")
  for (b in seq_len(n_subsamples)) {
    set.seed(seed + 1000 + b)
    keep <- make_subsample()
    subsamples <- sample_ids[keep]
    case_now <- subsamples[ploidy[keep] == "4x"]
    ctrl_now <- subsamples[ploidy[keep] == "2x"]
    
    fb <- run_once(
      subset_samples = subsamples,
      case_now = case_now, ctrl_now = ctrl_now,
      thr_override = if (fixed_threshold) ref_thr else NULL
    )
    
    if (topM_edges > 0) {
      Eb_keyed <- topM_edge_keys(fb$data_input, M = topM_edges, thr = if (fixed_threshold) ref_thr else fb$thr)
      jacc <- jaccard(E_full_keyed, Eb_keyed)
    } else {
      Eb <- adj_to_edges(fb$adj)
      jacc <- jaccard(E_full_keyed, edge_key(Eb))
    }
    
    k <- k + 1
    rob_list[[k]] <- data.frame(
      type = "Subsample80", run = b,
      edge_jaccard = jacc,
      module_ari   = module_ari(ref$modules, fb$modules),
      switch_recovery = switch_recovery(ref$switches, fb$switches),
      thr_used = fb$thr,
      n_switches = length(unique(fb$switches))
    )
    switches_list[[k + 1]] <- unique(fb$switches)
    
    if (progress) setTxtProgressBar(pb_sub, b)
  }
  
  robust_df <- do.call(rbind, Filter(Negate(is.null), rob_list))
  if (is.null(robust_df) || nrow(robust_df) == 0) {
    stop("No robustness rows were collected. Check earlier warnings for failed refits.")
  }
  
  # Save per-run detailed results
  write.csv(robust_df, "Robustness_Runs_Detail.csv", row.names = FALSE)
  
  # ---------- Diagnostics before aggregation ----------
  diag_metric <- function(x) {
    c(N = length(x),
      N_NA = sum(is.na(x)),
      pct_NA = round(100 * mean(is.na(x)), 1),
      mean = if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE),
      sd   = if (all(is.na(x))) NA_real_ else sd(x,   na.rm = TRUE),
      q25  = if (all(is.na(x))) NA_real_ else as.numeric(quantile(x, 0.25, na.rm = TRUE)),
      med  = if (all(is.na(x))) NA_real_ else as.numeric(quantile(x, 0.50, na.rm = TRUE)),
      q75  = if (all(is.na(x))) NA_real_ else as.numeric(quantile(x, 0.75, na.rm = TRUE)))
  }
  message(sprintf("Rows collected: %d  |  LOO=%d  Subsample80=%d",
                  nrow(robust_df),
                  sum(robust_df$type == "LOO"),
                  sum(robust_df$type == "Subsample80")))
  message("Overall metric diagnostics:")
  print(rbind(edge_jaccard    = diag_metric(robust_df$edge_jaccard),
              module_ari      = diag_metric(robust_df$module_ari),
              switch_recovery = diag_metric(robust_df$switch_recovery)))
  message("Diagnostics by type:")
  print(aggregate(edge_jaccard    ~ type, data = robust_df, FUN = function(x) diag_metric(x)))
  print(aggregate(module_ari      ~ type, data = robust_df, FUN = function(x) diag_metric(x)))
  print(aggregate(switch_recovery ~ type, data = robust_df, FUN = function(x) diag_metric(x)))
  
  # Even if columns contain NA, compute mean/sd per type safely
  safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  safe_sd   <- function(x) if (all(is.na(x))) NA_real_ else sd(x,   na.rm = TRUE)
  
  robust_summary <- aggregate(cbind(edge_jaccard, module_ari, switch_recovery) ~ type,
                              data = robust_df,
                              FUN  = function(x) c(mean = safe_mean(x), sd = safe_sd(x)))
  
  robust_summary <- data.frame(
    type               = robust_summary$type,
    edge_jaccard_mean  = robust_summary$edge_jaccard[, "mean"],
    edge_jaccard_sd    = robust_summary$edge_jaccard[, "sd"],
    module_ari_mean    = robust_summary$module_ari[, "mean"],
    module_ari_sd      = robust_summary$module_ari[, "sd"],
    switch_rec_mean    = robust_summary$switch_recovery[, "mean"],
    switch_rec_sd      = robust_summary$switch_recovery[, "sd"]
  )
  write.csv(robust_summary, "Supplementary_Robustness_Summary.csv", row.names = FALSE)
  
  # Consensus switches across ref + robustness runs
  all_switches <- unlist(switches_list, use.names = FALSE)
  if (length(all_switches) > 0) {
    tab <- sort(table(all_switches), decreasing = TRUE)
    total_runs <- length(switches_list)
    consensus <- data.frame(
      gene = names(tab),
      times = as.integer(tab),
      prop  = as.integer(tab) / total_runs,
      in_reference = names(tab) %in% ref_switches,
      stringsAsFactors = FALSE
    )
    write.csv(consensus, "Supplementary_Consensus_Switches.csv", row.names = FALSE)
  } else {
    write.csv(data.frame(gene=character(0), times=integer(0), prop=numeric(0), in_reference=logical(0)),
              "Supplementary_Consensus_Switches.csv", row.names = FALSE)
  }
  
  # Robustness boxplots
  pdf("Supplementary_Fig_Robustness.pdf", width = 7, height = 6)
  par(mfrow = c(3,1), mar = c(4,4,2,1))
  boxplot(edge_jaccard ~ type, data = robust_df, main = "Edge Jaccard vs Full",
          ylab = "Jaccard", xlab = "", outline = FALSE)
  boxplot(module_ari ~ type, data = robust_df, main = "Module ARI vs Full",
          ylab = "ARI", xlab = "", outline = FALSE)
  boxplot(switch_recovery ~ type, data = robust_df, main = "Switch-gene Recovery",
          ylab = "Proportion", xlab = "", outline = FALSE)
  dev.off()
  
  message("==> Stability diagnostics: done")
}
