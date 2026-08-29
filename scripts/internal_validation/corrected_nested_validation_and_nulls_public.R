#!/usr/bin/env Rscript

# Corrective, isolated re-analysis for the fertile internal validation and
# three negative controls. This script does not modify any historical asset.

suppressPackageStartupMessages({
  library(edgeR)
  library(dplyr)
})

OUT <- normalizePath(dirname(dirname(commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))] %>% sub("^--file=", "", .))), mustWork = FALSE)
if (!nzchar(OUT) || is.na(OUT)) OUT <- file.path(getwd(), "audit/corrected_nested_validation_20260824")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "01_inputs_and_provenance"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "02_corrected_primary_loocv"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "03_corrected_negative_controls"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "04_before_after_comparison"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "05_logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "06_scripts"), recursive = TRUE, showWarnings = FALSE)

# Public copy: local/raw inputs are intentionally not mirrored. Supply them at
# execution time through the documented environment variables below.
env_path <- function(name) Sys.getenv(name, unset = "")
counts_path <- env_path("GSE250130_COUNTS_CSV")
meta_path <- env_path("GSE250130_METADATA_CSV")
historical_script_path <- env_path("HISTORICAL_TIMING_SCRIPT")
historical_primary_path <- env_path("HISTORICAL_PRIMARY_METRICS_CSV")
historical_null_rows_path <- env_path("HISTORICAL_NULL_ROWS_CSV")
historical_null_summary_path <- env_path("HISTORICAL_NULL_SUMMARY_CSV")
historical_audit_path <- env_path("NEGATIVE_CONTROL_AUDIT_MD")
historical_source_table_path <- env_path("FROZEN_SIGNATURE_SOURCE_CSV")

stop_if_missing <- function(paths) {
  miss <- paths[!nzchar(paths) | !file.exists(paths)]
  if (length(miss)) stop("Missing locked input(s): ", paste(miss, collapse = " | "))
}
stop_if_missing(c(counts_path, meta_path, historical_script_path, historical_primary_path,
                  historical_null_rows_path, historical_null_summary_path,
                  historical_audit_path, historical_source_table_path))

safe_cor_test <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]
  if (length(x) < 4 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(c(rho = NA_real_, p = NA_real_))
  }
  z <- suppressWarnings(cor.test(x, y, method = "spearman"))
  c(rho = unname(z$estimate), p = z$p.value)
}

association_table <- function(expr, days) {
  ans <- t(vapply(seq_len(nrow(expr)), function(j) safe_cor_test(expr[j, ], days), numeric(2)))
  out <- data.frame(gene = rownames(expr), rho = ans[, "rho"], p = ans[, "p"], stringsAsFactors = FALSE)
  out$padj <- p.adjust(out$p, method = "BH")
  out
}

score_train_test <- function(train_expr, test_expr, late, early) {
  genes <- unique(c(intersect(late, rownames(train_expr)), intersect(early, rownames(train_expr))))
  genes <- intersect(genes, rownames(test_expr))
  if (length(genes) == 0) return(list(train = rep(NA_real_, ncol(train_expr)), test = rep(NA_real_, ncol(test_expr))))
  mu <- rowMeans(train_expr[genes, , drop = FALSE], na.rm = TRUE)
  sigma <- apply(train_expr[genes, , drop = FALSE], 1, sd, na.rm = TRUE)
  sigma[!is.finite(sigma) | sigma == 0] <- 1
  z_train <- sweep(sweep(train_expr[genes, , drop = FALSE], 1, mu, "-"), 1, sigma, "/")
  z_test <- sweep(sweep(test_expr[genes, , drop = FALSE], 1, mu, "-"), 1, sigma, "/")
  pos <- intersect(late, genes); neg <- intersect(early, genes)
  train_score <- colMeans(z_train[pos, , drop = FALSE], na.rm = TRUE) - colMeans(z_train[neg, , drop = FALSE], na.rm = TRUE)
  test_score <- colMeans(z_test[pos, , drop = FALSE], na.rm = TRUE) - colMeans(z_test[neg, , drop = FALSE], na.rm = TRUE)
  list(train = as.numeric(train_score), test = as.numeric(test_score), mu = mu, sigma = sigma)
}

nearest_stage <- function(x) {
  stages <- c(3, 5, 7, 9, 11)
  vapply(x, function(z) stages[which.min(abs(stages - z))], numeric(1))
}

pairwise_concordance <- function(df) {
  pairs <- expand.grid(a = seq_len(nrow(df)), b = seq_len(nrow(df))) %>%
    filter(a < b, df$true_day[a] != df$true_day[b])
  if (!nrow(pairs)) return(NA_real_)
  mean(sign(df$true_day[pairs$a] - df$true_day[pairs$b]) == sign(df$pred_day[pairs$a] - df$pred_day[pairs$b]))
}

metrics_from_predictions <- function(df) {
  keep <- is.finite(df$pred_day) & is.finite(df$true_day)
  if (!all(keep)) return(data.frame(rho = NA_real_, mae = NA_real_, pairwise_concordance = NA_real_, exact_nearest_stage_agreement = NA_real_, n_predictions = sum(keep), status = "non_estimable"))
  data.frame(
    rho = suppressWarnings(cor(df$pred_day, df$true_day, method = "spearman")),
    mae = mean(abs(df$pred_day - df$true_day)),
    pairwise_concordance = pairwise_concordance(df),
    exact_nearest_stage_agreement = mean(nearest_stage(df$pred_day) == df$true_day),
    n_predictions = nrow(df),
    status = "estimable",
    stringsAsFactors = FALSE
  )
}

read_counts <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(x)[1] <- "gene"
  x <- x %>% group_by(gene) %>% summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
  mat <- as.matrix(x[, -1, drop = FALSE])
  rownames(mat) <- x$gene
  storage.mode(mat) <- "numeric"
  mat
}

counts <- read_counts(counts_path)
meta <- read.csv(meta_path, check.names = FALSE, stringsAsFactors = FALSE)
meta$sample <- as.character(meta$sample)
meta$day <- dplyr::case_when(
  meta$group_simple == "LH3" ~ 3,
  meta$group_simple == "LH5" ~ 5,
  meta$group_simple == "Fertile_LH7" ~ 7,
  meta$group_simple == "LH9" ~ 9,
  meta$group_simple == "LH11" ~ 11,
  TRUE ~ NA_real_
)
fertile_meta <- meta %>% filter(!is.na(day)) %>% arrange(day, sample)
if (nrow(fertile_meta) != 18) stop("Locked fertile sample count mismatch: expected 18, observed ", nrow(fertile_meta))
if (!all(fertile_meta$sample %in% colnames(counts))) stop("Metadata/count sample mismatch")
if (!all(fertile_meta$sample %in% colnames(counts))) stop("Metadata/count sample mismatch")
counts <- counts[, fertile_meta$sample, drop = FALSE]

input_snapshot <- data.frame(
  item = c("counts", "metadata", "historical_script", "historical_primary", "historical_null_rows", "historical_null_summary", "negative_control_audit", "frozen_signature_source"),
  path = c(counts_path, meta_path, historical_script_path, historical_primary_path, historical_null_rows_path, historical_null_summary_path, historical_audit_path, historical_source_table_path),
  sha256 = NA_character_, bytes = NA_real_, stringsAsFactors = FALSE
)
for (j in seq_len(nrow(input_snapshot))) {
  input_snapshot$sha256[j] <- sub("  .*", "", system2("shasum", c("-a", "256", shQuote(input_snapshot$path[j])), stdout = TRUE))
  input_snapshot$bytes[j] <- file.info(input_snapshot$path[j])$size
}
write.csv(input_snapshot, file.path(OUT, "01_inputs_and_provenance", "INPUT_FILE_SNAPSHOT_MD5.csv"), row.names = FALSE)
write.csv(fertile_meta[, c("sample", "group_simple", "day")], file.path(OUT, "01_inputs_and_provenance", "GSE250130_FERTILE_SAMPLE_LOCK.csv"), row.names = FALSE)

# A per-fold object is built once. Filtering is training-only. To preserve the
# original scale definition as closely as possible, TMM factors and logCPM are
# then calculated on the union of training plus held-out counts after applying
# the training-derived keep mask. Thus the held-out sample can contribute to
# label-independent library normalization, but not to the gene universe.
build_fold <- function(i) {
  train_meta <- fertile_meta[-i, , drop = FALSE]
  test_meta <- fertile_meta[i, , drop = FALSE]
  train_samples <- train_meta$sample
  all_samples <- c(train_samples, test_meta$sample)
  y_train <- DGEList(counts = counts[, train_samples, drop = FALSE])
  keep <- filterByExpr(y_train, group = factor(train_meta$group_simple))
  kept_genes <- rownames(y_train)[keep]
  y_all <- DGEList(counts = counts[kept_genes, all_samples, drop = FALSE])
  y_all <- calcNormFactors(y_all)
  expr_all <- cpm(y_all, log = TRUE, prior.count = 1)
  train_expr <- expr_all[, train_samples, drop = FALSE]
  test_expr <- expr_all[, test_meta$sample, drop = FALSE]
  stats <- association_table(train_expr, train_meta$day)
  dynamic <- stats %>% filter(!is.na(padj), padj < 0.05)
  late <- dynamic %>% arrange(desc(rho), gene) %>% slice_head(n = 40) %>% pull(gene)
  early <- dynamic %>% arrange(rho, gene) %>% slice_head(n = 40) %>% pull(gene)
  strict_pool <- stats %>%
    mutate(expr_mean = rowMeans(train_expr[gene, , drop = FALSE], na.rm = TRUE)) %>%
    filter(!is.na(padj), gene %in% rownames(train_expr), padj > 0.50,
           abs(rho) <= quantile(abs(rho), 0.25, na.rm = TRUE),
           expr_mean >= quantile(expr_mean, 0.25, na.rm = TRUE)) %>% pull(gene) %>% unique()
  fallback <- FALSE
  fallback_pool <- character()
  pool <- strict_pool
  if (length(pool) < 80) {
    fallback <- TRUE
    fallback_pool <- stats %>% filter(gene %in% rownames(train_expr), !is.na(padj), padj > 0.25) %>% arrange(abs(rho), gene) %>% slice_head(n = 300) %>% pull(gene) %>% unique()
    pool <- fallback_pool
  }
  list(
    i = i, train_meta = train_meta, test_meta = test_meta,
    train_expr = train_expr, test_expr = test_expr, gene_ids = rownames(train_expr),
    keep_n = length(kept_genes), stats = stats, dynamic_n = nrow(dynamic),
    primary_late = late, primary_early = early,
    strict_pool = strict_pool, strict_pool_n = length(strict_pool),
    fallback = fallback, fallback_pool = fallback_pool, pool = pool,
    estimable_primary = length(late) == 40 && length(early) == 40,
    tmm_includes_test = TRUE
  )
}

folds <- lapply(seq_len(nrow(fertile_meta)), build_fold)

primary_preds <- vector("list", length(folds))
primary_audit <- vector("list", length(folds))
for (i in seq_along(folds)) {
  f <- folds[[i]]
  status <- if (f$estimable_primary) "estimable" else "non_estimable_insufficient_BH_dynamic_genes"
  pred <- NA_real_; intercept <- NA_real_; slope <- NA_real_
  if (f$estimable_primary) {
    sc <- score_train_test(f$train_expr, f$test_expr, f$primary_late, f$primary_early)
    train_score <- sc$train
    fit <- lm(f$train_meta$day ~ train_score)
    pred <- unname(predict(fit, newdata = data.frame(train_score = sc$test)))
    intercept <- unname(coef(fit)[1]); slope <- unname(coef(fit)[2])
  }
  primary_preds[[i]] <- data.frame(sample = f$test_meta$sample, observed_day = f$test_meta$day, predicted_day = pred,
                                    signed_error = pred - f$test_meta$day, absolute_error = abs(pred - f$test_meta$day),
                                    nearest_predicted_stage = if (is.finite(pred)) nearest_stage(pred) else NA_real_,
                                    fold = i, stringsAsFactors = FALSE)
  primary_audit[[i]] <- data.frame(
    fold = i, held_out_sample = f$test_meta$sample, observed_day = f$test_meta$day,
    training_n = nrow(f$train_meta), retained_gene_count = f$keep_n, dynamic_gene_count = f$dynamic_n,
    early_n = length(f$primary_early), late_n = length(f$primary_late),
    early_genes = paste(f$primary_early, collapse = ";"), late_genes = paste(f$primary_late, collapse = ";"),
    calibration_intercept = intercept, calibration_slope = slope,
    tmm_includes_heldout = f$tmm_includes_test, status = status, stringsAsFactors = FALSE
  )
}
primary_pred_df <- bind_rows(primary_preds)
primary_audit_df <- bind_rows(primary_audit)
primary_metrics <- metrics_from_predictions(data.frame(pred_day = primary_pred_df$predicted_day, true_day = primary_pred_df$observed_day))
primary_summary <- cbind(data.frame(analysis = "corrected_primary_nested_LOOCV", stringsAsFactors = FALSE), primary_metrics)
write.csv(primary_pred_df, file.path(OUT, "02_corrected_primary_loocv", "corrected_primary_loocv_predictions.csv"), row.names = FALSE)
write.csv(primary_audit_df, file.path(OUT, "02_corrected_primary_loocv", "corrected_primary_loocv_fold_audit.csv"), row.names = FALSE)
write.csv(primary_summary, file.path(OUT, "02_corrected_primary_loocv", "corrected_primary_loocv_summary.csv"), row.names = FALSE)

run_null <- function(control, base_seed, B = 250) {
  metric_rows <- vector("list", B)
  prov_rows <- list()
  for (iter in seq_len(B)) {
    pred_rows <- vector("list", length(folds))
    iter_status <- "estimable"
    failed <- character()
    for (i in seq_along(folds)) {
      f <- folds[[i]]
      seed_used <- base_seed + iter * 1000 + i
      set.seed(seed_used)
      train_day <- f$train_meta$day
      if (control == "permuted_LH") train_day <- sample(train_day, length(train_day), replace = FALSE)
      pool <- if (control == "data_derived_low_association") f$pool else f$gene_ids
      if (length(pool) < 80) {
        iter_status <- "non_estimable"; failed <- c(failed, as.character(i))
        pred_rows[[i]] <- data.frame(sample = f$test_meta$sample, true_day = f$test_meta$day, pred_day = NA_real_)
        next
      }
      sampled <- sample(pool, 80, replace = FALSE)
      sampled_stats <- association_table(f$train_expr[sampled, , drop = FALSE], train_day)
      sampled_stats$padj <- p.adjust(sampled_stats$p, method = "BH")
      late <- sampled_stats %>% arrange(desc(rho), gene) %>% slice_head(n = 40) %>% pull(gene)
      early <- sampled_stats %>% arrange(rho, gene) %>% slice_head(n = 40) %>% pull(gene)
      if (length(late) != 40 || length(early) != 40) {
        iter_status <- "non_estimable"; failed <- c(failed, as.character(i))
        pred_rows[[i]] <- data.frame(sample = f$test_meta$sample, true_day = f$test_meta$day, pred_day = NA_real_)
        next
      }
      sc <- score_train_test(f$train_expr, f$test_expr, late, early)
      train_score <- sc$train
      fit <- lm(train_day ~ train_score)
      pred <- unname(predict(fit, newdata = data.frame(train_score = sc$test)))
      pred_rows[[i]] <- data.frame(sample = f$test_meta$sample, true_day = f$test_meta$day, pred_day = pred)
      prov_rows[[length(prov_rows) + 1L]] <- data.frame(
        control = control, iteration = iter, fold = i, held_out_sample = f$test_meta$sample,
        seed_used = seed_used, training_gene_universe_n = length(f$gene_ids),
        strict_pool_n = f$strict_pool_n, fallback_used = f$fallback,
        fallback_pool_n = length(f$fallback_pool), sampled_n = length(sampled),
        sampled_80_genes = paste(sampled, collapse = ";"),
        late_40_genes = paste(late, collapse = ";"), early_40_genes = paste(early, collapse = ";"),
        status = "estimable", stringsAsFactors = FALSE
      )
    }
    pred_df <- bind_rows(pred_rows)
    m <- metrics_from_predictions(pred_df)
    metric_rows[[iter]] <- data.frame(control = control, iteration = iter, base_seed = base_seed,
                                      failed_folds = paste(failed, collapse = ","), status = iter_status,
                                      m, stringsAsFactors = FALSE)
  }
  list(metrics = bind_rows(metric_rows), provenance = bind_rows(prov_rows))
}

set.seed(20260824)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
random_res <- run_null("random_80_gene", 20260824, 250)
permuted_res <- run_null("permuted_LH", 20260825, 250)
low_res <- run_null("data_derived_low_association", 20260826, 250)

write.csv(random_res$metrics, file.path(OUT, "03_corrected_negative_controls", "corrected_random80_null_metrics.csv"), row.names = FALSE)
write.csv(permuted_res$metrics, file.path(OUT, "03_corrected_negative_controls", "corrected_permutedLH_null_metrics.csv"), row.names = FALSE)
write.csv(low_res$metrics, file.path(OUT, "03_corrected_negative_controls", "corrected_low_association_null_metrics.csv"), row.names = FALSE)
write.csv(bind_rows(random_res$provenance, permuted_res$provenance, low_res$provenance), file.path(OUT, "03_corrected_negative_controls", "corrected_negative_control_fold_provenance.csv"), row.names = FALSE)

obs_rho <- primary_metrics$rho
summarize_null <- function(df) {
  x <- df$rho[is.finite(df$rho)]
  b <- sum(x >= obs_rho)
  data.frame(
    control = unique(df$control)[1], B = nrow(df), b = b, raw_empirical_P = b / nrow(df), plus1_empirical_P = (b + 1) / (nrow(df) + 1),
    maximum_null_rho = max(x), median_null_rho = median(x), q025_null_rho = unname(quantile(x, 0.025)), q975_null_rho = unname(quantile(x, 0.975)),
    median_MAE = median(df$mae, na.rm = TRUE), median_pairwise_concordance = median(df$pairwise_concordance, na.rm = TRUE),
    observed_corrected_rho = obs_rho, observed_corrected_MAE = primary_metrics$mae, observed_corrected_pairwise = primary_metrics$pairwise_concordance,
    estimable_iterations = sum(df$status == "estimable"), non_estimable_iterations = sum(df$status != "estimable"), stringsAsFactors = FALSE
  )
}
null_summary <- bind_rows(summarize_null(random_res$metrics), summarize_null(permuted_res$metrics), summarize_null(low_res$metrics))
write.csv(null_summary, file.path(OUT, "03_corrected_negative_controls", "corrected_negative_control_summary.csv"), row.names = FALSE)

write.csv(data.frame(
  fold = seq_along(folds), held_out_sample = fertile_meta$sample, training_n = vapply(folds, function(x) nrow(x$train_meta), numeric(1)),
  training_gene_universe_n = vapply(folds, function(x) length(x$gene_ids), numeric(1)), strict_pool_n = vapply(folds, function(x) x$strict_pool_n, numeric(1)),
  fallback_used = vapply(folds, function(x) x$fallback, logical(1)), fallback_pool_n = vapply(folds, function(x) length(x$fallback_pool), numeric(1)),
  stringsAsFactors = FALSE
), file.path(OUT, "03_corrected_negative_controls", "corrected_fold_universe_and_low_association_audit.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(OUT, "05_logs", "R_sessionInfo.txt"))
write.csv(data.frame(
  execution_date = as.character(Sys.time()), working_directory = getwd(), R_version = R.version.string,
  edgeR_version = as.character(packageVersion("edgeR")), dplyr_version = as.character(packageVersion("dplyr")),
  RNGkind = paste(RNGkind(), collapse = ";"), master_seed = 20260824, random_seed = 20260824,
  permuted_seed = 20260825, low_association_seed = 20260826, B = 250, fertile_samples = nrow(fertile_meta),
  stage_distribution = paste(names(table(fertile_meta$day)), as.integer(table(fertile_meta$day)), collapse = ";"),
  tmm_normalization = "training-derived filter mask; TMM factors and logCPM calculated on training plus held-out counts; held-out contributes only label-independent library normalization",
  analysis_rerun_scope = "corrected primary LOOCV plus three negative controls only", formal_files_modified = "NO", stringsAsFactors = FALSE
), file.path(OUT, "01_inputs_and_provenance", "EXECUTION_PROVENANCE_CORE.csv"), row.names = FALSE)

cat("CORRECTED_RUN_COMPLETE\n")
cat("primary_rho=", primary_metrics$rho, " primary_mae=", primary_metrics$mae, " primary_pairwise=", primary_metrics$pairwise_concordance, " exact=", primary_metrics$exact_nearest_stage_agreement, "\n", sep = "")
print(null_summary)
