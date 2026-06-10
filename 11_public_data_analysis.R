#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

required_pkgs <- c("readxl", "edgeR", "ggplot2", "dplyr", "tidyr", "tibble")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(readxl)
  library(edgeR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

find_project_root <- function(start_dir) {
  current <- normalizePath(start_dir, mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(current, "0511")) && dir.exists(file.path(current, "github"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Could not locate project root containing both '0511' and 'github'.")
}

download_if_missing <- function(url, dest) {
  if (file.exists(dest) && file.size(dest) > 0) return(invisible(dest))
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  message("Downloading: ", basename(dest))
  ok <- TRUE
  tryCatch(
    download.file(url, dest, mode = "wb", quiet = FALSE, method = "libcurl"),
    error = function(e) ok <<- FALSE
  )
  if (!ok || !file.exists(dest) || file.size(dest) == 0) {
    if (file.exists(dest)) file.remove(dest)
    cmd <- sprintf("curl -L --retry 5 --retry-delay 2 '%s' -o '%s'", url, dest)
    status <- system(cmd)
    if (!identical(status, 0L) || !file.exists(dest) || file.size(dest) == 0) {
      stop("Failed to download: ", url)
    }
  }
  invisible(dest)
}

read_gse262053_table <- function(path) {
  dat <- read.csv(gzfile(path), check.names = FALSE)
  gene <- sub("^.*\\|", "", dat$gene_id)
  mat <- as.matrix(dat[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  agg <- rowsum(mat, group = gene, reorder = FALSE)
  data.frame(gene = rownames(agg), agg, check.names = FALSE, row.names = NULL)
}

merge_count_tables <- function(df1, df2) {
  merged <- merge(df1, df2, by = "gene", all = TRUE, sort = FALSE)
  rownames(merged) <- merged$gene
  merged <- merged[, -1, drop = FALSE]
  merged[is.na(merged)] <- 0
  mat <- as.matrix(merged)
  storage.mode(mat) <- "numeric"
  mat
}

calc_log2cpm <- function(count_mat) {
  dge <- DGEList(counts = count_mat)
  dge <- calcNormFactors(dge)
  cpm(dge, log = TRUE, prior.count = 1)
}

safe_t_test <- function(values, groups) {
  if (length(unique(groups)) != 2) return(NA_real_)
  split_vals <- split(values, groups)
  if (any(vapply(split_vals, length, integer(1)) < 2)) return(NA_real_)
  t.test(values ~ groups)$p.value
}

p_to_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

save_plot_both <- function(plot, path_stub, width, height) {
  ggsave(paste0(path_stub, ".pdf"), plot, width = width, height = height)
  ggsave(paste0(path_stub, ".png"), plot, width = width, height = height, dpi = 300)
}

make_group_point_plot <- function(df, label_df, title, subtitle, ncol = 6) {
  ggplot(df, aes(x = group, y = value, color = group)) +
    geom_jitter(width = 0.10, size = 2.2, alpha = 0.9) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.35, color = "black") +
    facet_wrap(~ gene, scales = "free_y", ncol = ncol) +
    geom_text(
      data = label_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 2.5,
      lineheight = 0.95,
      vjust = 0
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.32))) +
    coord_cartesian(clip = "off") +
    labs(title = title, subtitle = subtitle, x = NULL, y = "log2(CPM + 1)") +
    scale_color_brewer(palette = "Dark2") +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "#f0f0f0"),
      axis.text.x = element_text(angle = 20, hjust = 1),
      plot.margin = margin(10, 10, 10, 10)
    )
}

script_dir <- get_script_dir()
project_root <- find_project_root(script_dir)
analysis_dir <- file.path(script_dir, "output")
fig_dir <- file.path(analysis_dir, "figures")
result_dir <- file.path(analysis_dir, "results")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_bw(base_size = 11))
message("Project root: ", project_root)
message("Output directory: ", analysis_dir)

marker_path <- file.path(project_root, "0511", "TableS3_基因标志物表达差异.xlsx")
data_dir <- file.path(project_root, "0511", "data")
r1_path <- file.path(data_dir, "GSE262053_A549_R1_gene_count_matrix.csv.gz")
r2_path <- file.path(data_dir, "GSE262053_A549_R2_gene_count_matrix.csv.gz")

download_if_missing(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262053/suppl/GSE262053_A549_R1_gene_count_matrix.csv.gz",
  r1_path
)
download_if_missing(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262053/suppl/GSE262053_A549_R2_gene_count_matrix.csv.gz",
  r2_path
)

marker_tbl <- read_excel(marker_path, skip = 1) %>%
  transmute(gene = as.character(gene), avg_log2FC = as.numeric(avg_log2FC)) %>%
  filter(!is.na(gene), gene != "")

marker_genes <- unique(marker_tbl$gene)
rsi_genes <- c("AR", "JUN", "STAT1", "PRKCB", "RELA", "ABL1", "SUMO1", "PAK2", "HDAC1", "IRF1")
core_genes <- c("CDKN1A", "KDM4B", "CCNB1", "TOP2A", "GDF15", "POLR2A", "FDXR", "TP53I3", "BAX", "ZMAT3", "PHLDA3", "IFI16", "MDM2", "NINJ1")
candidate_genes <- unique(c(core_genes, marker_genes, rsi_genes))

gene_relevance <- c(
  CDKN1A = "p53 target / radiosensitization candidate",
  KDM4B = "DNA damage response / radiosensitization candidate",
  CCNB1 = "cell cycle / G2M axis",
  TOP2A = "cell cycle / G2M axis",
  GDF15 = "stress response / secreted signaling",
  POLR2A = "SCENIC regulon candidate",
  FDXR = "p53 target",
  TP53I3 = "p53 target",
  BAX = "apoptosis / p53 target",
  ZMAT3 = "p53 target",
  PHLDA3 = "p53 target",
  IFI16 = "DNA damage / innate sensing",
  MDM2 = "p53 feedback",
  NINJ1 = "stress response"
)

rsi_formula <- function(expr_named_vec) {
  -0.0098009 * expr_named_vec["AR"] +
    0.0128283 * expr_named_vec["JUN"] +
    0.0254552 * expr_named_vec["STAT1"] -
    0.0017589 * expr_named_vec["PRKCB"] -
    0.0038171 * expr_named_vec["RELA"] +
    0.1070213 * expr_named_vec["ABL1"] -
    0.0002509 * expr_named_vec["SUMO1"] -
    0.0092431 * expr_named_vec["PAK2"] -
    0.0204469 * expr_named_vec["HDAC1"] -
    0.0441683 * expr_named_vec["IRF1"]
}

sample_meta <- data.frame(
  sample = c(
    "A549_WT_DMSO_R1", "A549_WT_DMSO_R2",
    "A549_KO_DMSO_R1", "A549_KO_DMSO_R2",
    "A549_KO175_DMSO_R1", "A549_KO175_DMSO_R2",
    "A549_KO273_DMSO_R1", "A549_KO273_DMSO_R2"
  ),
  group = c("WT", "WT", "KO", "KO", "R175H", "R175H", "R273H", "R273H"),
  stringsAsFactors = FALSE
)

r1 <- read_gse262053_table(r1_path)
r2 <- read_gse262053_table(r2_path)
gse262053_counts <- merge_count_tables(r1, r2)

missing_samples <- setdiff(sample_meta$sample, colnames(gse262053_counts))
if (length(missing_samples) > 0) {
  stop("Missing expected samples: ", paste(missing_samples, collapse = ", "))
}

module_a_log2cpm <- calc_log2cpm(gse262053_counts[, sample_meta$sample, drop = FALSE])

module_a_genes <- intersect(candidate_genes, rownames(module_a_log2cpm))
module_a_long <- module_a_log2cpm[module_a_genes, , drop = FALSE] %>%
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "value") %>%
  left_join(sample_meta, by = "sample")

module_a_summary <- module_a_long %>%
  group_by(gene, group) %>%
  summarise(mean_log2cpm = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_log2cpm) %>%
  mutate(
    KO_vs_WT = KO - WT,
    R175H_vs_WT = R175H - WT,
    R273H_vs_WT = R273H - WT,
    relevance = unname(gene_relevance[gene])
  ) %>%
  arrange(gene)
write.csv(module_a_summary, file.path(result_dir, "moduleA_group_means_and_directions.csv"), row.names = FALSE)

module_a_pairwise <- bind_rows(lapply(c("KO", "R175H", "R273H"), function(g) {
  module_a_long %>%
    filter(group %in% c("WT", g)) %>%
    group_by(gene) %>%
    summarise(
      contrast = paste0(g, "_vs_WT"),
      mean_diff = mean(value[group == g]) - mean(value[group == "WT"]),
      p_value = safe_t_test(value, group),
      .groups = "drop"
    )
})) %>%
  mutate(relevance = unname(gene_relevance[gene])) %>%
  arrange(gene, contrast)
write.csv(module_a_pairwise, file.path(result_dir, "moduleA_pairwise_tests.csv"), row.names = FALSE)

module_a_pairwise_table <- module_a_pairwise %>%
  mutate(
    comparison = recode(
      contrast,
      KO_vs_WT = "KO vs WT",
      R175H_vs_WT = "R175H vs WT",
      R273H_vs_WT = "R273H vs WT"
    ),
    significance = dplyr::case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  select(gene, comparison, mean_diff, p_value, significance, relevance) %>%
  arrange(gene, comparison)
write.csv(module_a_pairwise_table, file.path(result_dir, "moduleA_pairwise_pvalue_table.csv"), row.names = FALSE)

module_a_pairwise_wide <- module_a_pairwise %>%
  mutate(
    significance = p_to_stars(p_value),
    p_text = paste0(signif(p_value, 3), significance)
  ) %>%
  select(gene, contrast, mean_diff, p_value, p_text, relevance) %>%
  pivot_wider(
    names_from = contrast,
    values_from = c(mean_diff, p_value, p_text),
    names_sep = "_"
  ) %>%
  arrange(gene)
write.csv(module_a_pairwise_wide, file.path(result_dir, "moduleA_pairwise_pvalue_table_wide.csv"), row.names = FALSE)

module_a_core_plot_genes <- intersect(c("CCNB1", "TOP2A", "POLR2A", "BAX", "ZMAT3", "PHLDA3"), module_a_genes)
module_a_plot_df <- module_a_long %>%
  filter(gene %in% module_a_core_plot_genes) %>%
  mutate(
    gene = factor(gene, levels = module_a_core_plot_genes),
    group = factor(group, levels = c("WT", "KO", "R175H", "R273H"))
  )

module_a_plot_p <- module_a_plot_df %>%
  group_by(gene) %>%
  summarise(
    x = 2.5,
    y = max(value) + 0.15,
    label = paste0(
      "KO vs WT: p=", signif(safe_t_test(value[group %in% c("WT", "KO")], droplevels(group[group %in% c("WT", "KO")])), 2),
      "\nR175H vs WT: p=", signif(safe_t_test(value[group %in% c("WT", "R175H")], droplevels(group[group %in% c("WT", "R175H")])), 2),
      "\nR273H vs WT: p=", signif(safe_t_test(value[group %in% c("WT", "R273H")], droplevels(group[group %in% c("WT", "R273H")])), 2)
    ),
    .groups = "drop"
  )

module_a_plot_star_p <- module_a_plot_df %>%
  group_by(gene, group) %>%
  summarise(y = max(value) + 0.15, .groups = "drop") %>%
  filter(group != "WT") %>%
  left_join(
    module_a_pairwise %>%
      filter(gene %in% module_a_core_plot_genes) %>%
      mutate(
        group = case_when(
          contrast == "KO_vs_WT" ~ "KO",
          contrast == "R175H_vs_WT" ~ "R175H",
          contrast == "R273H_vs_WT" ~ "R273H",
          TRUE ~ NA_character_
        ),
        label = case_when(
          is.na(p_value) ~ "",
          p_value < 0.001 ~ "***",
          p_value < 0.01 ~ "**",
          p_value < 0.05 ~ "*",
          TRUE ~ ""
        )
      ) %>%
      select(gene, group, label),
    by = c("gene", "group")
  ) %>%
  mutate(
    x = factor(group, levels = c("WT", "KO", "R175H", "R273H")),
    gene = factor(gene, levels = module_a_core_plot_genes)
  ) %>%
  filter(label != "")

module_a_points <- make_group_point_plot(
  module_a_plot_df,
  module_a_plot_p,
  "Module A: A549 TP53 WT / KO / mutant expression",
  "DMSO samples only; points are independent replicates; black bar is group mean",
  ncol = 6
)
save_plot_both(module_a_points, file.path(fig_dir, "moduleA_core_gene_points"), 16.8, 4.0)

module_a_points_star <- make_group_point_plot(
  module_a_plot_df,
  module_a_plot_star_p,
  "Module A: A549 TP53 WT / KO / mutant expression",
  "DMSO samples only; points are independent replicates; black bar is group mean",
  ncol = 6
)
save_plot_both(module_a_points_star, file.path(fig_dir, "moduleA_core_gene_points_star"), 16.8, 4.0)

module_a_heatmap_df <- module_a_summary %>%
  filter(gene %in% module_a_core_plot_genes) %>%
  select(gene, WT, KO, R175H, R273H) %>%
  pivot_longer(-gene, names_to = "group", values_to = "value") %>%
  mutate(
    gene = factor(gene, levels = rev(module_a_core_plot_genes)),
    group = factor(group, levels = c("WT", "KO", "R175H", "R273H"))
  )

module_a_heatmap_sig <- module_a_pairwise %>%
  filter(gene %in% module_a_core_plot_genes) %>%
  mutate(
    group = case_when(
      contrast == "KO_vs_WT" ~ "KO",
      contrast == "R175H_vs_WT" ~ "R175H",
      contrast == "R273H_vs_WT" ~ "R273H",
      TRUE ~ NA_character_
    ),
    significance = dplyr::case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    gene = factor(gene, levels = rev(module_a_core_plot_genes)),
    group = factor(group, levels = c("WT", "KO", "R175H", "R273H"))
  ) %>%
  filter(significance != "")

module_a_heatmap <- ggplot(module_a_heatmap_df, aes(x = group, y = gene, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 2.6) +
  geom_text(
    data = module_a_heatmap_sig,
    aes(x = group, y = gene, label = significance),
    inherit.aes = FALSE,
    nudge_y = 0.28,
    size = 3.0
  ) +
  scale_fill_gradient(low = "white", high = "#b2182b") +
  labs(title = "Module A: core candidate genes across TP53 states", x = NULL, y = NULL, fill = "log2(CPM + 1)") +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1))
save_plot_both(module_a_heatmap, file.path(fig_dir, "moduleA_core_gene_heatmap"), 6.8, 5.8)

module_a_heatmap_pvals <- module_a_pairwise %>%
  filter(gene %in% module_a_core_plot_genes) %>%
  mutate(
    group = case_when(
      contrast == "KO_vs_WT" ~ "KO",
      contrast == "R175H_vs_WT" ~ "R175H",
      contrast == "R273H_vs_WT" ~ "R273H",
      TRUE ~ NA_character_
    ),
    gene = factor(gene, levels = rev(module_a_core_plot_genes)),
    group = factor(group, levels = c("WT", "KO", "R175H", "R273H")),
    p_label = paste0("p=", signif(p_value, 2))
  ) %>%
  select(gene, group, p_label)

module_a_heatmap_with_p <- ggplot(module_a_heatmap_df, aes(x = group, y = gene, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 2.4, nudge_y = -0.12) +
  geom_text(
    data = module_a_heatmap_pvals,
    aes(x = group, y = gene, label = p_label),
    inherit.aes = FALSE,
    nudge_y = 0.18,
    size = 2.1
  ) +
  scale_fill_gradient(low = "white", high = "#b2182b") +
  labs(
    title = "Module A: core candidate genes across TP53 states",
    subtitle = "Exact P values are shown for KO, R175H, and R273H versus WT",
    x = NULL,
    y = NULL,
    fill = "log2(CPM + 1)"
  ) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1))
save_plot_both(module_a_heatmap_with_p, file.path(fig_dir, "moduleA_core_gene_heatmap_with_p"), 7.2, 6.2)

module_a_dotplot_df <- module_a_summary %>%
  filter(gene %in% module_a_core_plot_genes) %>%
  select(gene, WT, KO, R175H, R273H) %>%
  pivot_longer(-gene, names_to = "group", values_to = "value") %>%
  mutate(
    gene = factor(gene, levels = rev(module_a_core_plot_genes)),
    group = factor(group, levels = c("WT", "KO", "R175H", "R273H"))
  )

module_a_dotplot_p <- module_a_pairwise %>%
  filter(gene %in% module_a_core_plot_genes) %>%
  mutate(
    group = case_when(
      contrast == "KO_vs_WT" ~ "KO",
      contrast == "R175H_vs_WT" ~ "R175H",
      contrast == "R273H_vs_WT" ~ "R273H",
      TRUE ~ NA_character_
    ),
    neglog10_p = -log10(p_value)
  ) %>%
  select(gene, group, p_value, neglog10_p)

module_a_dotplot_df <- module_a_dotplot_df %>%
  left_join(module_a_dotplot_p, by = c("gene", "group")) %>%
  mutate(
    gene = factor(gene, levels = rev(module_a_core_plot_genes)),
    group = factor(group, levels = c("WT", "KO", "R175H", "R273H")),
    neglog10_p_plot = ifelse(as.character(group) == "WT", 0.2, neglog10_p),
    significance = dplyr::case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

module_a_dotplot <- ggplot(module_a_dotplot_df, aes(x = group, y = gene)) +
  geom_point(aes(size = neglog10_p_plot, fill = value), shape = 21, color = "#404040", stroke = 0.35) +
  geom_text(aes(label = sprintf("%.2f", value)), size = 2.4, color = "black") +
  geom_text(
    data = subset(module_a_dotplot_df, significance != ""),
    aes(label = significance),
    nudge_y = 0.24,
    size = 3.0,
    color = "black"
  ) +
  scale_fill_gradient(low = "#f7fbff", high = "#d7301f") +
  scale_size_continuous(range = c(5, 16), breaks = c(0.5, 1, 2), labels = c("WT", "0.1", "0.01"), name = "-log10(P)") +
  labs(
    title = "Module A: core candidate genes across TP53 states",
    subtitle = "Dot color indicates mean expression; dot size indicates significance versus WT",
    x = NULL,
    y = NULL,
    fill = "log2(CPM + 1)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_line(color = "#eeeeee", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "#f3f3f3", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1)
  )
save_plot_both(module_a_dotplot, file.path(fig_dir, "moduleA_core_gene_dotplot"), 7.4, 6.2)

missing_rsi_genes <- setdiff(rsi_genes, rownames(module_a_log2cpm))
if (length(missing_rsi_genes) > 0) {
  stop("Missing RSI genes in expression matrix: ", paste(missing_rsi_genes, collapse = ", "))
}

module_a_rsi_scores <- sapply(sample_meta$sample, function(s) rsi_formula(module_a_log2cpm[rsi_genes, s]))
module_a_rsi_df <- sample_meta %>%
  mutate(RSI = module_a_rsi_scores, group = factor(group, levels = c("WT", "KO", "R175H", "R273H")))
write.csv(module_a_rsi_df, file.path(result_dir, "moduleA_rsi_scores.csv"), row.names = FALSE)

module_a_rsi_plot <- ggplot(module_a_rsi_df, aes(x = group, y = RSI, color = group)) +
  geom_jitter(width = 0.10, size = 2.4, alpha = 0.9) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.35, color = "black") +
  labs(title = "Module A: RSI across A549 TP53 states", subtitle = "RSI formula applied to log2(CPM + 1) expression", x = NULL, y = "RSI") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")
save_plot_both(module_a_rsi_plot, file.path(fig_dir, "moduleA_rsi_points"), 6.5, 4.8)

summary_lines <- c(
  "# Module A Summary",
  "",
  paste0("- Project root: ", project_root),
  paste0("- Output figures: ", fig_dir),
  paste0("- Output tables: ", result_dir),
  paste0("- Samples used: ", paste(sample_meta$sample, collapse = ", ")),
  paste0("- Core genes visualized: ", paste(module_a_core_plot_genes, collapse = ", ")),
  paste0("- Total candidate genes tested: ", length(module_a_genes))
)
writeLines(summary_lines, file.path(analysis_dir, "README_output.md"))

message("Module A analysis completed.")
