suppressPackageStartupMessages({
  library(Hmisc)       # rcorr
  library(igraph)
  library(data.table)
  library(mclust)      # adjustedRandIndex
  library(corpcor)     # shrinkage correlation
})

run_stability <- function(data,
                          input_parameter,
                          input_file,
                          n_subsamples = 100,
                          frac_keep = 0.8,
                          use_shrinkage = TRUE,
                          seed = 123) {
  
  set.seed(seed)
  
  # -----------------------------
  # Output paths
  # -----------------------------
  outdir <- file.path(input_parameter$path, "robustness")
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  file_detail   <- file.path(outdir, "Robustness_Runs_Detail.csv")
  file_summary  <- file.path(outdir, "Supplementary_Robustness_Summary.csv")
  file_fig      <- file.path(outdir, "Supplementary_Fig_Robustness.pdf")
  file_cons_sw  <- file.path(outdir, "Supplementary_Consensus_Switches.csv")
  
  # -----------------------------
  # Metadata + (opzionale) genotype centering
  # -----------------------------
  samples <- colnames(data)
  
  if (!is.null(input_file$metadata)) {
    meta <- input_file$metadata
    stopifnot(all(samples %in% meta$sample))
    meta <- meta[match(samples, meta$sample), ]
    stopifnot(all(c("sample", "ploidy", "genotype") %in% colnames(meta)))
  } else {
    ploidy <- ifelse(samples %in% input_file$control, "CTRL",
                     ifelse(samples %in% input_file$case, "CASE", NA))
    if (any(is.na(ploidy)))
      stop("Some samples are neither in control nor in case.")
    
    geno_guess <- sub("[._-].*$", "", samples)
    meta <- data.frame(
      sample   = samples,
      ploidy   = ploidy,
      genotype = geno_guess,
      stringsAsFactors = FALSE
    )
  }
  
  # Flag: usare o no il genotype-centering?
  use_gc <- TRUE
  if (!is.null(input_parameter$use_genotype_center)) {
    use_gc <- isTRUE(input_parameter$use_genotype_center)
  }
  
  genotype_center <- function(mat, meta) {
    Z <- mat
    for (g in unique(meta$genotype)) {
      idx <- which(meta$genotype == g)
      if (length(idx) > 1) {
        mu  <- rowMeans(Z[, idx, drop = FALSE])
        Z[, idx] <- sweep(Z[, idx, drop = FALSE], 1, mu, "-")
      }
    }
    Z
  }
  
  if (use_gc) {
    message(">> Using genotype-centered expression matrix")
    data_gc <- genotype_center(data, meta)
  } else {
    message(">> Using raw expression matrix (NO genotype centering)")
    data_gc <- data
  }
  
  # -----------------------------
  # Core helpers
  # -----------------------------
  corr_edge_list_from_matrix <- function(M,
                                         method_adj = input_parameter$correction_method,
                                         type = input_parameter$type_correlation,
                                         shrink = use_shrinkage) {
    if (!shrink) {
      rc <- Hmisc::rcorr(t(M), type = type)
      rho <- rc$r; pval <- rc$P
    } else {
      rho <- corpcor::cor.shrink(t(M))
      n   <- ncol(M)
      tval <- rho * sqrt((n - 2) / pmax(1e-8, 1 - rho^2))
      pval <- 2 * stats::pt(-abs(tval), df = n - 2)
      diag(pval) <- 0
    }
    
    pval_adj <- p.adjust(pval, method = method_adj)
    pval_adj <- matrix(
      pval_adj,
      ncol = ncol(pval),
      dimnames = list(rownames(rho), colnames(rho))
    )
    
    ut <- upper.tri(rho)
    data.frame(
      source      = rownames(rho)[row(rho)[ut]],
      target      = colnames(rho)[col(rho)[ut]],
      correlation = rho[ut],
      pval        = pval[ut],
      pval_adj    = pval_adj[ut],
      stringsAsFactors = FALSE
    )
  }
  
  build_network <- function(corr_edges,
                            thr_q = input_parameter$threshold_prc_corr,
                            thr_padj = input_parameter$threshold_pval_adj_corr) {
    thr_corr <- stats::quantile(corr_edges$correlation, thr_q, na.rm = TRUE)
    thr_corr <- round(as.numeric(thr_corr), 4)
    keep <- which((abs(corr_edges$correlation) >= thr_corr) &
                    (corr_edges$pval_adj <= thr_padj))
    corr_edges[keep, , drop = FALSE]
  }
  
  weighted_adj <- function(net) {
    if (nrow(net) == 0) {
      return(matrix(0, 0, 0))
    }
    src <- net$source
    tgt <- net$target
    w   <- net$correlation
    nodes <- unique(c(src, tgt))
    N <- length(nodes)
    W <- matrix(0, N, N, dimnames = list(nodes, nodes))
    for (i in seq_len(nrow(net))) {
      W[src[i], tgt[i]] <- w[i]
      W[tgt[i], src[i]] <- w[i]
    }
    W
  }
  
  kmeans_idx <- function(W,
                         k        = input_parameter$num_clusters,
                         iter_max = input_parameter$iter_max,
                         nstart   = input_parameter$num_repeats) {
    if (nrow(W) < k || k < 2) {
      cl <- rep(1, nrow(W))
      return(list(
        idx  = factor(cl, levels = 1),
        size = nrow(W),
        WSS  = NA,
        TWSS = NA
      ))
    }
    model <- stats::kmeans(W, centers = k, iter.max = iter_max, nstart = nstart)
    list(
      idx  = factor(model$cluster, levels = seq_len(k)),
      size = model$size,
      WSS  = model$withinss,
      TWSS = model$tot.withinss
    )
  }
  
  # -----------------------------
  # Cartography & switch helpers
  # -----------------------------
  apcc_fn <- function(W) {
    if (length(W) == 0) {
      return(setNames(numeric(0), character(0)))
    }
    Wz <- W
    Wz[Wz == 0] <- NA
    out <- rowMeans(Wz, na.rm = TRUE)
    names(out) <- rownames(W)
    out
  }
  
  binarize <- function(W) {
    out <- W
    out[out != 0] <- 1
    out
  }
  
  cluster_matrix <- function(idx, W) {
    if (length(idx) == 0) return(matrix(0, 0, 0))
    df <- data.frame(
      node = names(idx),
      cl   = as.integer(idx),
      stringsAsFactors = FALSE
    )
    l <- split(df$node, df$cl)
    pairs <- data.table::rbindlist(lapply(l, function(y) {
      if (length(y) == 1) {
        data.frame(source = y, target = y, stringsAsFactors = FALSE)
      } else {
        as.data.frame(expand.grid(y, y), stringsAsFactors = FALSE)
      }
    }))
    g <- igraph::graph_from_data_frame(pairs, directed = FALSE)
    g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
    C <- as.matrix(igraph::as_adj(g, type = "both"))
    C[match(rownames(W), rownames(C)), match(colnames(W), colnames(C))]
  }
  
  deg_lists <- function(W, idx) {
    if (length(idx) == 0) {
      return(list(
        deg  = setNames(numeric(0), character(0)),
        ideg = setNames(numeric(0), character(0))
      ))
    }
    C  <- cluster_matrix(idx, W)
    A  <- binarize(W)
    zc <- C * A
    d  <- rowSums(A)
    names(d) <- rownames(W)
    id <- rowSums(zc)
    names(id) <- rownames(W)
    list(deg = d, ideg = id)
  }
  
  compute_Pz <- function(nodes, deg, ideg) {
    m <- mean(deg[nodes])
    s <- stats::sd(deg[nodes])
    if (is.na(s) || s == 0) s <- 1
    z <- (ideg[nodes] - m) / s
    P <- 1 - (ideg[nodes] / pmax(1, deg[nodes]))^2
    list(
      P = setNames(P, nodes),
      z = setNames(z, nodes)
    )
  }
  
  node_role <- function(z, P) {
    hub    <- ifelse(z < 2.5, "non local hub", "local hub")
    region <- ifelse(z < 2.5 & P <= 0.04,  "R1",
                     ifelse(z < 2.5 & P <= 0.625, "R2",
                            ifelse(z < 2.5 & P <= 0.8, "R3",
                                   ifelse(z < 2.5,          "R4",
                                          ifelse(P <= 0.3,  "R5",
                                                 ifelse(P <= 0.75, "R6", "R7"))))))
    type   <- ifelse(z < 2.5 & P <= 0.04,  "Ultra-peripheral nodes",
                     ifelse(z < 2.5 & P <= 0.625, "Peripheral nodes",
                            ifelse(z < 2.5 & P <= 0.8, "Non-hub connectors",
                                   ifelse(z < 2.5,          "Non-hub kinless nodes",
                                          ifelse(P <= 0.3,  "Provincial hubs",
                                                 ifelse(P <= 0.75, "Connector hubs",
                                                        "Kinless hubs"))))))
    list(hub = hub, region = region, type = type)
  }
  
  hub_class_fn <- function(apcc_vec, deg_vec) {
    ifelse(apcc_vec > 0 & apcc_vec < 0.5 & deg_vec >= 5, "DATE",
           ifelse(apcc_vec >= 0.5 & deg_vec >= 5,        "PARTY",
                  ifelse(apcc_vec < 0 & deg_vec >= 5,    "FIGHT CLUB",
                         ifelse(deg_vec < 5,             "no hub", "not available"))))
  }
  
  # -----------------------------
  # Full reference run
  # -----------------------------
  corr_full <- corr_edge_list_from_matrix(data_gc)
  net_full  <- build_network(corr_full)
  W_full    <- weighted_adj(net_full)
  
  km_full   <- kmeans_idx(W_full)
  nodes_full <- rownames(W_full)
  
  apcc_full <- apcc_fn(W_full)
  dgl_full  <- deg_lists(W_full, km_full$idx)
  par_full  <- compute_Pz(nodes_full, dgl_full$deg, dgl_full$ideg)
  role_full <- node_role(par_full$z, par_full$P)
  hubc_full <- hub_class_fn(apcc_full, dgl_full$deg)
  
  attr_full <- data.frame(
    node                = nodes_full,
    Region              = role_full$region,
    Hub_classification  = hubc_full,
    stringsAsFactors    = FALSE
  )
  switches_full <- attr_full$node[
    attr_full$Region == "R4" & attr_full$Hub_classification == "FIGHT CLUB"
  ]
  
  # -----------------------------
  # Metrics
  # -----------------------------
  jaccard_edges <- function(netA, netB) {
    if (nrow(netA) == 0 && nrow(netB) == 0) return(1)
    if (nrow(netA) == 0 || nrow(netB) == 0) return(0)
    key <- function(df) paste(pmin(df$source, df$target),
                              pmax(df$source, df$target), sep = "||")
    a <- unique(key(netA))
    b <- unique(key(netB))
    length(intersect(a, b)) / length(unique(c(a, b)))
  }
  
  ari_modules <- function(WA, WB) {
    if (nrow(WA) < 3 || nrow(WB) < 3) return(NA_real_)
    ia <- intersect(rownames(WA), rownames(WB))
    if (length(ia) < 3) return(NA_real_)
    kma <- kmeans_idx(WA[ia, ia])$idx
    kmb <- kmeans_idx(WB[ia, ia])$idx
    mclust::adjustedRandIndex(as.integer(kma), as.integer(kmb))
  }
  
  switch_recovery <- function(sw_ref, sw_run) {
    if (length(sw_ref) == 0 && length(sw_run) == 0) return(1)
    if (length(sw_ref) == 0) return(NA_real_)
    length(intersect(sw_ref, sw_run)) / length(unique(sw_ref))
  }
  
  # -----------------------------
  # Single resampling run
  # -----------------------------
  single_run <- function(sel_samples, tag, run_id) {
    M <- data_gc[, sel_samples, drop = FALSE]
    corr <- corr_edge_list_from_matrix(M)
    net  <- build_network(corr)
    W    <- weighted_adj(net)
    
    jac  <- jaccard_edges(net_full, net)
    ari  <- tryCatch(ari_modules(W_full, W), error = function(e) NA_real_)
    
    if (nrow(W) >= 3) {
      km     <- kmeans_idx(W)
      apcc_v <- apcc_fn(W)
      dgl    <- deg_lists(W, km$idx)
      nodesR <- rownames(W)
      par    <- compute_Pz(nodesR, dgl$deg, dgl$ideg)
      role   <- node_role(par$z, par$P)
      hubc   <- hub_class_fn(apcc_v, dgl$deg)
      attr   <- data.frame(
        node               = nodesR,
        Region             = role$region,
        Hub_classification = hubc,
        stringsAsFactors   = FALSE
      )
      sw <- attr$node[
        attr$Region == "R4" & attr$Hub_classification == "FIGHT CLUB"
      ]
    } else {
      sw <- character(0)
    }
    
    sw_rec <- switch_recovery(switches_full, sw)
    
    row <- data.frame(
      run_type        = tag,
      run_id          = run_id,
      n_samples       = length(sel_samples),
      edge_jaccard    = jac,
      module_ARI      = ari,
      switch_recovery = sw_rec,
      stringsAsFactors = FALSE
    )
    
    list(row = row, switches = sw)
  }
  
  details_list <- list()
  all_switches <- list()
  
  # -----------------------------
  # LOO
  # -----------------------------
  for (s in seq_along(samples)) {
    sel <- samples[-s]
    res <- single_run(sel, "LOO", paste0("LOO_", samples[s]))
    details_list[[length(details_list) + 1]] <- res$row
    all_switches[[res$row$run_id]] <- res$switches
  }
  
  # -----------------------------
  # 80% subsampling stratificato per "genotype"
  # -----------------------------
  for (b in seq_len(n_subsamples)) {
    keep <- c()
    for (g in unique(meta$genotype)) {
      idx <- which(meta$genotype == g)
      k   <- max(1, floor(length(idx) * frac_keep))
      keep <- c(keep, sample(idx, k))
    }
    keep <- sort(unique(keep))
    sel_samples <- samples[keep]
    
    res <- single_run(sel_samples, "Subsample80", sprintf("SUB_%03d", b))
    details_list[[length(details_list) + 1]] <- res$row
    all_switches[[res$row$run_id]] <- res$switches
  }
  
  details_df <- data.table::rbindlist(details_list, fill = TRUE)
  data.table::fwrite(details_df, file_detail)
  
  # -----------------------------
  # Consensus switches (≥50% of runs)
  # -----------------------------
  all_sw_vec <- unlist(all_switches, use.names = FALSE)
  
  if (length(all_sw_vec) > 0) {
    tab <- sort(table(all_sw_vec), decreasing = TRUE)
    thr <- ceiling(0.5 * length(all_switches))
    cons <- names(tab)[tab >= thr]
    
    cons_df <- data.table::data.table(
      switch_gene = cons,
      count       = as.integer(tab[cons]),
      runs        = rep(length(all_switches), length(cons))
    )
  } else {
    # Nessuno switch ottenuto in nessun run:
    # creo una tabella vuota con 0 righe e 3 colonne
    cons_df <- data.table::data.table(
      switch_gene = character(),
      count       = integer(),
      runs        = integer()
    )
  }
  
  data.table::fwrite(cons_df, file_cons_sw)
  
  # -----------------------------
  # Summary + Figure
  # -----------------------------
  summary_df <- details_df[
    , .(
      n                     = .N,
      edge_jaccard_mean     = mean(edge_jaccard,    na.rm = TRUE),
      edge_jaccard_sd       = sd(edge_jaccard,      na.rm = TRUE),
      module_ARI_mean       = mean(module_ARI,      na.rm = TRUE),
      module_ARI_sd         = sd(module_ARI,        na.rm = TRUE),
      switch_recovery_mean  = mean(switch_recovery, na.rm = TRUE),
      switch_recovery_sd    = sd(switch_recovery,   na.rm = TRUE)
    ),
    by = run_type
  ]
  data.table::fwrite(summary_df, file_summary)
  
  pdf(file_fig, width = 7, height = 7)
  op <- par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
  boxplot(edge_jaccard ~ run_type, data = details_df,
          ylab = "Edge Jaccard vs Full", xlab = "", main = "Edge stability")
  abline(h = median(details_df$edge_jaccard[details_df$run_type=="Subsample80"], na.rm = TRUE),
         lty = 2)
  boxplot(module_ARI ~ run_type, data = details_df,
          ylab = "Adjusted Rand Index", xlab = "", main = "Module stability (k-means)")
  abline(h = median(details_df$module_ARI[details_df$run_type=="Subsample80"], na.rm = TRUE),
         lty = 2)
  boxplot(switch_recovery ~ run_type, data = details_df,
          ylab = "Switch recovery", xlab = "", main = "Switch robustness")
  abline(h = median(details_df$switch_recovery[details_df$run_type=="Subsample80"], na.rm = TRUE),
         lty = 2)
  par(op)
  dev.off()
  
  invisible(list(detail = details_df, summary = summary_df))
}

  
  # -----------------------------
  # Summary + Figure
  # -----------------------------
  summary_df <- details_df[
    , .(
      n                     = .N,
      edge_jaccard_mean     = mean(edge_jaccard,    na.rm = TRUE),
      edge_jaccard_sd       = stats::sd(edge_jaccard,    na.rm = TRUE),
      module_ARI_mean       = mean(module_ARI,      na.rm = TRUE),
      module_ARI_sd         = stats::sd(module_ARI,      na.rm = TRUE),
      switch_recovery_mean  = mean(switch_recovery, na.rm = TRUE),
      switch_recovery_sd    = stats::sd(switch_recovery, na.rm = TRUE)
    ),
    by = run_type
  ]
  data.table::fwrite(summary_df, file_summary)
  
  pdf(file_fig, width = 7, height = 7)
  op <- par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
  boxplot(edge_jaccard ~ run_type, data = details_df,
          ylab = "Edge Jaccard vs Full", xlab = "", main = "Edge stability")
  abline(h = stats::median(details_df$edge_jaccard[details_df$run_type == "Subsample80"], na.rm = TRUE),
         lty = 2)
  boxplot(module_ARI ~ run_type, data = details_df,
          ylab = "Adjusted Rand Index", xlab = "", main = "Module stability (k-means)")
  abline(h = stats::median(details_df$module_ARI[details_df$run_type == "Subsample80"], na.rm = TRUE),
         lty = 2)
  boxplot(switch_recovery ~ run_type, data = details_df,
          ylab = "Switch recovery", xlab = "", main = "Switch robustness")
  abline(h = stats::median(details_df$switch_recovery[details_df$run_type == "Subsample80"], na.rm = TRUE),
         lty = 2)
  par(op)
  dev.off()
  
  invisible(list(detail = details_df, summary = summary_df))
}
