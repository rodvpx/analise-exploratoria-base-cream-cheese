options(stringsAsFactors = FALSE)

# Cream Cheese - Exploracao de dados sem regressao/classificacao
# Foco: qualidade, estrutura sensorial, PCA, clusterizacao e MDS.

if (!requireNamespace("foreign", quietly = TRUE)) {
  stop("Pacote 'foreign' nao encontrado. Instale com install.packages('foreign').")
}

if (!requireNamespace("cluster", quietly = TRUE)) {
  stop("Pacote 'cluster' nao encontrado. Instale com install.packages('cluster').")
}

data_path <- "SPSS_CreamCheese/SPSS- for CreamCheese/CreamCheeseRawData2000.sav"
out_dir <- "outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cluster_palette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E")

add_cluster_legend <- function(labels, palette = cluster_palette, title = "Cluster") {
  labs <- sort(unique(as.integer(labels)))
  legend(
    "topright",
    legend = paste(title, labs),
    col = palette[labs],
    pch = 19,
    pt.cex = 1.2,
    bty = "n"
  )
}

plot_matrix_heatmap <- function(mat, main, low = "#2166AC", mid = "white", high = "#B2182B") {
  pal <- colorRampPalette(c(low, mid, high))(100)
  zlim <- range(mat, na.rm = TRUE)
  if (zlim[1] < 0 && zlim[2] > 0) {
    lim <- max(abs(zlim))
    zlim <- c(-lim, lim)
  }
  op <- par(mar = c(8, 8, 4, 5))
  image(
    x = seq_len(nrow(mat)),
    y = seq_len(ncol(mat)),
    z = mat,
    col = pal,
    zlim = zlim,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = main
  )
  axis(1, at = seq_len(nrow(mat)), labels = rownames(mat), las = 2, cex.axis = 0.75)
  axis(2, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.75)
  box()
  par(op)
}

cat("Lendo base SPSS em:", data_path, "\n")
df <- foreign::read.spss(
  file = data_path,
  to.data.frame = TRUE,
  use.value.labels = FALSE
)

required_cols <- c(
  "Productname", "Productnumber", "Panellist",
  "Replicate", "Session", "Servingorder"
)
missing_required <- setdiff(required_cols, names(df))
if (length(missing_required) > 0) {
  stop("Colunas obrigatorias ausentes: ", paste(missing_required, collapse = ", "))
}

id_vars <- c(
  "Productname", "Productnumber", "Panellist",
  "Replicate", "Session", "Servingorder"
)
sensory_vars <- setdiff(names(df), id_vars)

# =========================
# 1) Limpeza e preparo
# =========================
df_before <- df
n_rows_before <- nrow(df_before)

# Padroniza texto/codificacao em identificadores
df$Productname <- trimws(gsub("[[:space:]]+", " ", as.character(df$Productname)))
df$Productname[df$Productname == ""] <- NA_character_

# Corrige tipos numericos de identificadores
for (v in c("Productnumber", "Panellist", "Replicate", "Session", "Servingorder")) {
  df[[v]] <- suppressWarnings(as.integer(as.character(df[[v]])))
}

# Corrige tipo de atributos sensoriais
coercion_na <- setNames(rep(0L, length(sensory_vars)), sensory_vars)
for (v in sensory_vars) {
  old <- df[[v]]
  old_na <- is.na(old)
  new <- suppressWarnings(as.numeric(as.character(old)))
  coercion_na[v] <- sum(is.na(new) & !old_na)
  df[[v]] <- new
}

# Remove duplicatas exatas
dup_rows_before <- sum(duplicated(df))
if (dup_rows_before > 0) {
  df <- df[!duplicated(df), ]
}

# Faixas esperadas para IDs experimentais
invalid_id <- data.frame(
  variable = c("Panellist", "Replicate", "Session", "Servingorder"),
  invalid_n = c(
    sum(is.na(df$Panellist) | !(df$Panellist %in% 1:8)),
    sum(is.na(df$Replicate) | !(df$Replicate %in% 1:3)),
    sum(is.na(df$Session) | !(df$Session %in% 1:6)),
    sum(is.na(df$Servingorder) | !(df$Servingorder %in% 1:5))
  ),
  stringsAsFactors = FALSE
)

# Marca IDs fora de faixa como NA para imputacao segura posterior
df$Panellist[is.na(df$Panellist) | !(df$Panellist %in% 1:8)] <- NA_integer_
df$Replicate[is.na(df$Replicate) | !(df$Replicate %in% 1:3)] <- NA_integer_
df$Session[is.na(df$Session) | !(df$Session %in% 1:6)] <- NA_integer_
df$Servingorder[is.na(df$Servingorder) | !(df$Servingorder %in% 1:5)] <- NA_integer_

# Prepara vetor de painelista para imputacao por grupo
pan_for_imp <- df$Panellist
pan_for_imp[is.na(pan_for_imp)] <- 0L

# Diagnostico de missing antes
missing_before <- colSums(is.na(df))
miss_sens_before <- missing_before[sensory_vars]

# Imputacao robusta para atributos sensoriais: mediana por painelista -> mediana global
imputed_n <- setNames(rep(0L, length(sensory_vars)), sensory_vars)
for (v in sensory_vars) {
  x <- df[[v]]
  miss <- is.na(x)
  if (any(miss)) {
    med_by_pan <- tapply(x, pan_for_imp, median, na.rm = TRUE)
    fill <- med_by_pan[as.character(pan_for_imp[miss])]
    med_global <- median(x, na.rm = TRUE)
    fill[is.na(fill)] <- med_global
    x[miss] <- fill
    imputed_n[v] <- sum(miss)
    df[[v]] <- x
  }
}

# Imputacao para IDs: moda simples (mais frequente)
mode_int <- function(x) {
  xt <- table(x, useNA = "no")
  as.integer(names(xt)[which.max(xt)])
}
for (v in c("Panellist", "Replicate", "Session", "Servingorder")) {
  if (anyNA(df[[v]])) {
    df[[v]][is.na(df[[v]])] <- mode_int(df[[v]])
  }
}

# Winsorizacao para reduzir impacto de outliers extremos (1% - 99%)
winsor_low_n <- setNames(rep(0L, length(sensory_vars)), sensory_vars)
winsor_high_n <- setNames(rep(0L, length(sensory_vars)), sensory_vars)
for (v in sensory_vars) {
  x <- df[[v]]
  q <- quantile(x, probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE)
  lo <- q[1]
  hi <- q[2]
  winsor_low_n[v] <- sum(x < lo, na.rm = TRUE)
  winsor_high_n[v] <- sum(x > hi, na.rm = TRUE)
  x <- pmax(x, lo)
  x <- pmin(x, hi)
  df[[v]] <- x
}

# Revalida duplicacao no desenho experimental apos limpeza
dup_design_before <- sum(duplicated(df[id_vars]))
if (dup_design_before > 0) {
  # Consolida repeticoes exatas de desenho usando media dos atributos sensoriais
  df <- aggregate(
    df[sensory_vars],
    by = df[id_vars],
    FUN = mean
  )
}

# Tipagem final
product_order <- aggregate(
  Productnumber ~ Productname,
  data = df[!is.na(df$Productname) & !is.na(df$Productnumber), ],
  FUN = min
)
product_order <- product_order[order(product_order$Productnumber, product_order$Productname), ]
product_levels <- as.character(product_order$Productname)

df$Productname <- factor(df$Productname, levels = product_levels, ordered = TRUE)
df$Productnumber <- as.integer(df$Productnumber)
df$Panellist <- as.factor(df$Panellist)
df$Replicate <- factor(df$Replicate, levels = 1:3, ordered = TRUE)
df$Session <- factor(df$Session, levels = 1:6, ordered = TRUE)
df$Servingorder <- factor(df$Servingorder, levels = 1:5, ordered = TRUE)

# Relatorios de qualidade
quality_overview <- data.frame(
  metric = c(
    "rows_before_cleaning",
    "rows_after_cleaning",
    "exact_duplicates_removed",
    "design_duplicates_before_consolidation",
    "total_missing_before",
    "total_missing_after",
    "coercion_na_total",
    "winsorized_low_total",
    "winsorized_high_total"
  ),
  value = c(
    n_rows_before,
    nrow(df),
    dup_rows_before,
    dup_design_before,
    sum(missing_before),
    sum(is.na(df)),
    sum(coercion_na),
    sum(winsor_low_n),
    sum(winsor_high_n)
  ),
  stringsAsFactors = FALSE
)
write.csv(quality_overview, file.path(out_dir, "data_quality_overview.csv"), row.names = FALSE)

var_quality <- data.frame(
  variable = sensory_vars,
  coercion_na = as.integer(coercion_na[sensory_vars]),
  missing_before = as.integer(miss_sens_before[sensory_vars]),
  missing_after = as.integer(colSums(is.na(df[sensory_vars]))[sensory_vars]),
  imputed_n = as.integer(imputed_n[sensory_vars]),
  winsor_low_n = as.integer(winsor_low_n[sensory_vars]),
  winsor_high_n = as.integer(winsor_high_n[sensory_vars]),
  stringsAsFactors = FALSE
)
write.csv(var_quality, file.path(out_dir, "data_quality_by_variable.csv"), row.names = FALSE)
write.csv(invalid_id, file.path(out_dir, "data_quality_invalid_ids.csv"), row.names = FALSE)

cat("Dimensoes:", nrow(df), "linhas x", ncol(df), "colunas\n")
cat("Variaveis sensoriais:", length(sensory_vars), "\n\n")
cat("Nota: metricas ACC/AUC/F1/Specificity sao supervisionadas e exigem rotulo-alvo.\n")
cat("Neste script (sem classificacao/regressao), usamos metricas de clusterizacao e estrutura.\n\n")
cat("Limpeza/preparo concluido:\n")
print(quality_overview)
cat("\n")

# =========================
# 2) Qualidade e estrutura
# =========================
missing_total <- sum(is.na(df))
dup_rows <- sum(duplicated(df))
dup_design <- sum(duplicated(df[id_vars]))

cat("Missing total:", missing_total, "\n")
cat("Linhas duplicadas (dataset completo):", dup_rows, "\n")
cat("Duplicadas pelo desenho experimental:", dup_design, "\n\n")

cat("Contagem por produto:\n")
print(table(df$Productname))
cat("\nContagem por painelista:\n")
print(table(df$Panellist))
cat("\n")

prod_map <- table(df$Productname, df$Productnumber)
write.csv(as.data.frame(prod_map), file.path(out_dir, "productname_productnumber_map.csv"), row.names = FALSE)

multi_code <- rowSums(prod_map > 0) > 1
if (any(multi_code)) {
  cat("Aviso: produtos com mais de um Productnumber:\n")
  print(names(multi_code)[multi_code])
  cat("\n")
}

# =========================
# 3) Descritivas sensoriais
# =========================
sensory_mat <- as.matrix(df[sensory_vars])

sensory_summary <- data.frame(
  variable = sensory_vars,
  min = apply(sensory_mat, 2, min),
  q25 = apply(sensory_mat, 2, quantile, probs = 0.25),
  median = apply(sensory_mat, 2, median),
  mean = apply(sensory_mat, 2, mean),
  q75 = apply(sensory_mat, 2, quantile, probs = 0.75),
  max = apply(sensory_mat, 2, max),
  sd = apply(sensory_mat, 2, sd),
  cv = apply(sensory_mat, 2, function(x) sd(x) / mean(x))
)

sensory_summary <- sensory_summary[order(-sensory_summary$sd), ]
write.csv(sensory_summary, file.path(out_dir, "sensory_summary.csv"), row.names = FALSE)

cat("Top 8 atributos com maior dispersao (sd):\n")
print(head(sensory_summary[, c("variable", "sd", "cv")], 8))
cat("\n")

# =========================
# 4) Correcao por painelista
# =========================
# Z-score dentro de cada painelista para reduzir vies individual.
df_z <- df
for (v in sensory_vars) {
  means <- tapply(df[[v]], df$Panellist, mean)
  sds <- tapply(df[[v]], df$Panellist, sd)
  z <- (df[[v]] - means[df$Panellist]) / sds[df$Panellist]
  z[is.na(z)] <- 0
  df_z[[v]] <- as.numeric(z)
}

# Perfil medio por produto (bruto e corrigido por painelista)
product_profile_raw <- aggregate(df[sensory_vars], by = list(Productname = df$Productname), FUN = mean)
product_profile_z <- aggregate(df_z[sensory_vars], by = list(Productname = df_z$Productname), FUN = mean)

write.csv(product_profile_raw, file.path(out_dir, "product_profile_raw.csv"), row.names = FALSE)
write.csv(product_profile_z, file.path(out_dir, "product_profile_panelist_z.csv"), row.names = FALSE)

# =========================
# 5) PCA em perfil de produto
# =========================
pca_input <- as.matrix(product_profile_z[, sensory_vars])
rownames(pca_input) <- product_profile_z$Productname

pca <- prcomp(pca_input, center = TRUE, scale. = TRUE)
var_exp <- pca$sdev^2 / sum(pca$sdev^2)
cum_var <- cumsum(var_exp)

pca_scores <- data.frame(Productname = rownames(pca$x), pca$x, row.names = NULL)
pca_loadings <- data.frame(Attribute = rownames(pca$rotation), pca$rotation, row.names = NULL)
pca_variance <- data.frame(
  PC = paste0("PC", seq_along(var_exp)),
  variance_explained = var_exp,
  cumulative_variance = cum_var
)

write.csv(pca_scores, file.path(out_dir, "pca_scores.csv"), row.names = FALSE)
write.csv(pca_loadings, file.path(out_dir, "pca_loadings.csv"), row.names = FALSE)
write.csv(pca_variance, file.path(out_dir, "pca_variance.csv"), row.names = FALSE)

png(file.path(out_dir, "pca_scree.png"), width = 1000, height = 650)
op <- par(mar = c(5, 5, 4, 5))
barplot(
  100 * var_exp,
  names.arg = paste0("PC", seq_along(var_exp)),
  col = "#8DA0CB",
  border = "white",
  ylim = c(0, max(100 * var_exp) * 1.20),
  xlab = "Componente principal",
  ylab = "Variancia explicada (%)",
  main = "PCA - Variancia explicada por componente"
)
lines(seq_along(var_exp), 100 * cum_var, type = "b", pch = 19, col = "#D95F02", lwd = 2)
axis(4, at = seq(0, 100, by = 20), labels = seq(0, 100, by = 20), las = 1)
mtext("Variancia acumulada (%)", side = 4, line = 3)
legend(
  "topright",
  legend = c("Variancia por componente", "Variancia acumulada"),
  fill = c("#8DA0CB", NA),
  border = c("white", NA),
  lty = c(NA, 1),
  pch = c(NA, 19),
  col = c("#8DA0CB", "#D95F02"),
  bty = "n"
)
par(op)
dev.off()

png(file.path(out_dir, "pca_biplot_pc1_pc2.png"), width = 1100, height = 800)
scores12 <- pca$x[, 1:2, drop = FALSE]
load12 <- pca$rotation[, 1:2, drop = FALSE]
top_load <- order(rowSums(load12^2), decreasing = TRUE)[1:min(10, nrow(load12))]
score_lim <- range(scores12)
plot(
  scores12[, 1],
  scores12[, 2],
  pch = 19,
  cex = 1.5,
  col = "#1B9E77",
  xlab = paste0("PC1 (", round(100 * var_exp[1], 1), "%)"),
  ylab = paste0("PC2 (", round(100 * var_exp[2], 1), "%)"),
  main = "PCA - Produtos e principais atributos sensoriais",
  xlim = range(scores12[, 1]) * 1.25,
  ylim = range(scores12[, 2]) * 1.25
)
abline(h = 0, v = 0, col = "gray85", lty = 2)
text(scores12[, 1], scores12[, 2], labels = rownames(scores12), pos = 3, cex = 0.9)
arrow_scale <- 0.65 * max(abs(score_lim)) / max(abs(load12[top_load, ]))
arrows(
  0, 0,
  load12[top_load, 1] * arrow_scale,
  load12[top_load, 2] * arrow_scale,
  length = 0.08,
  col = "#D95F02",
  lwd = 1.4
)
text(
  load12[top_load, 1] * arrow_scale,
  load12[top_load, 2] * arrow_scale,
  labels = rownames(load12)[top_load],
  col = "#D95F02",
  cex = 0.85,
  pos = 4
)
legend(
  "topright",
  legend = c("Produtos", "Atributos com maior contribuicao"),
  col = c("#1B9E77", "#D95F02"),
  pch = c(19, NA),
  lty = c(NA, 1),
  bty = "n"
)
dev.off()

cat("Variancia explicada (primeiros 5 PCs):\n")
print(head(pca_variance, 5))
cat("\n")

# =========================
# 6) Clusterizacao de produtos
# =========================
# Funcoes auxiliares de metricas internas de clusterizacao
calc_davies_bouldin <- function(x, cl) {
  cl <- as.integer(cl)
  ids <- sort(unique(cl))
  centroids <- sapply(ids, function(k) colMeans(x[cl == k, , drop = FALSE]))
  if (is.vector(centroids)) {
    centroids <- matrix(centroids, ncol = 1)
  }
  centroids <- t(centroids)
  rownames(centroids) <- ids

  s <- sapply(ids, function(k) {
    pts <- x[cl == k, , drop = FALSE]
    ctd <- centroids[as.character(k), ]
    mean(sqrt(rowSums((pts - matrix(ctd, nrow(pts), ncol(pts), byrow = TRUE))^2)))
  })

  m <- as.matrix(dist(centroids))
  r <- matrix(0, nrow = length(ids), ncol = length(ids))
  for (i in seq_along(ids)) {
    for (j in seq_along(ids)) {
      if (i != j) {
        r[i, j] <- (s[i] + s[j]) / m[i, j]
      }
    }
  }
  mean(apply(r, 1, max))
}

calc_dunn <- function(x, cl) {
  cl <- as.integer(cl)
  ids <- sort(unique(cl))
  # max diametro intra-cluster
  intra_diam <- sapply(ids, function(k) {
    pts <- x[cl == k, , drop = FALSE]
    if (nrow(pts) <= 1) {
      return(0)
    }
    max(dist(pts))
  })
  max_intra <- max(intra_diam)
  if (max_intra == 0) {
    return(NA_real_)
  }

  # min distancia inter-clusters (single-link)
  min_inter <- Inf
  for (i in 1:(length(ids) - 1)) {
    for (j in (i + 1):length(ids)) {
      a <- x[cl == ids[i], , drop = FALSE]
      b <- x[cl == ids[j], , drop = FALSE]
      d <- as.matrix(dist(rbind(a, b)))
      n_a <- nrow(a)
      inter <- min(d[1:n_a, (n_a + 1):ncol(d), drop = FALSE])
      if (inter < min_inter) {
        min_inter <- inter
      }
    }
  }
  min_inter / max_intra
}

calc_ari <- function(labels_a, labels_b) {
  tab <- table(labels_a, labels_b)
  n <- sum(tab)
  if (n <= 1) {
    return(NA_real_)
  }
  comb2 <- function(v) v * (v - 1) / 2
  sum_ij <- sum(comb2(tab))
  sum_i <- sum(comb2(rowSums(tab)))
  sum_j <- sum(comb2(colSums(tab)))
  total <- comb2(n)
  expected <- (sum_i * sum_j) / total
  max_idx <- 0.5 * (sum_i + sum_j)
  (sum_ij - expected) / (max_idx - expected)
}

calc_withinss <- function(x, cl) {
  ids <- sort(unique(cl))
  wss <- 0
  for (k in ids) {
    pts <- x[cl == k, , drop = FALSE]
    ctr <- colMeans(pts)
    wss <- wss + sum((pts - matrix(ctr, nrow(pts), ncol(pts), byrow = TRUE))^2)
  }
  wss
}

calc_ch <- function(x, cl) {
  n <- nrow(x)
  k <- length(unique(cl))
  if (k <= 1 || n <= k) {
    return(NA_real_)
  }
  total_ss <- sum(scale(x, center = TRUE, scale = FALSE)^2)
  within_ss <- calc_withinss(x, cl)
  between_ss <- total_ss - within_ss
  (between_ss / (k - 1)) / (within_ss / (n - k))
}

fit_clusters <- function(method, k, x, dmat, hc_obj) {
  if (method == "kmeans") {
    return(as.integer(kmeans(x, centers = k, nstart = 100, iter.max = 200)$cluster))
  }
  if (method == "pam") {
    return(as.integer(cluster::pam(dmat, k = k, diss = TRUE)$clustering))
  }
  if (method == "hclust_ward") {
    return(as.integer(cutree(hc_obj, k = k)))
  }
  stop("Metodo nao suportado: ", method)
}

# Representacoes para comparar (estilo "varias abordagens" do script ILPD)
profile_raw <- as.matrix(product_profile_raw[, sensory_vars])
rownames(profile_raw) <- product_profile_raw$Productname
profile_raw <- scale(profile_raw)

profile_panel_z <- as.matrix(product_profile_z[, sensory_vars])
rownames(profile_panel_z) <- product_profile_z$Productname
profile_panel_z <- scale(profile_panel_z)

# Matriz produto x (painelista*atributo): aumenta sinal de interacao sensorial
panel_mean <- aggregate(df[sensory_vars], by = list(Productname = df$Productname, Panellist = df$Panellist), FUN = mean)
prod_levels <- levels(df$Productname)
pan_levels <- as.character(sort(unique(df$Panellist)))
panel_block <- matrix(NA_real_, nrow = length(prod_levels), ncol = length(sensory_vars) * length(pan_levels))
rownames(panel_block) <- prod_levels
col_idx <- 1
for (p in pan_levels) {
  subp <- panel_mean[panel_mean$Panellist == p, ]
  subp <- subp[match(prod_levels, as.character(subp$Productname)), ]
  block <- as.matrix(subp[, sensory_vars, drop = FALSE])
  colnames(block) <- paste0(colnames(block), "__P", p)
  panel_block[, col_idx:(col_idx + ncol(block) - 1)] <- block
  col_idx <- col_idx + ncol(block)
}
colnames(panel_block) <- as.vector(sapply(pan_levels, function(p) paste0(sensory_vars, "__P", p)))
panel_block <- scale(panel_block)

# Selecao exploratoria de atributos: variabilidade, redundancia e sinal de produto
attr_signal <- data.frame(
  variable = sensory_vars,
  sd = as.numeric(sapply(df[sensory_vars], sd)),
  iqr = as.numeric(sapply(df[sensory_vars], IQR)),
  unique_n = as.integer(sapply(df[sensory_vars], function(x) length(unique(x)))),
  eta_product = NA_real_,
  eta_panelist = NA_real_,
  eta_interaction = NA_real_,
  signal_score = NA_real_,
  low_information = FALSE,
  redundant_with = NA_character_,
  selected_uncorrelated = FALSE,
  stringsAsFactors = FALSE
)

for (i in seq_along(sensory_vars)) {
  v <- sensory_vars[i]
  fit <- aov(df[[v]] ~ df$Productname + df$Panellist + df$Productname:df$Panellist)
  tab <- summary(fit)[[1]]
  ss_total <- sum(tab[, "Sum Sq"], na.rm = TRUE)
  if (ss_total == 0 || is.na(ss_total)) {
    attr_signal$signal_score[i] <- 1e-6
    next
  }
  ss_prod <- tab[1, "Sum Sq"]
  ss_pan <- tab[2, "Sum Sq"]
  ss_int <- tab[3, "Sum Sq"]
  attr_signal$eta_product[i] <- ss_prod / ss_total
  attr_signal$eta_panelist[i] <- ss_pan / ss_total
  attr_signal$eta_interaction[i] <- ss_int / ss_total
  attr_signal$signal_score[i] <- max(attr_signal$eta_product[i] - attr_signal$eta_interaction[i], 1e-6)
}

attr_signal$low_information <- attr_signal$sd < 1e-8 | attr_signal$unique_n < 3
attr_signal <- attr_signal[order(-attr_signal$signal_score, -attr_signal$sd), ]

cor_for_selection <- cor(df[sensory_vars], method = "spearman")
selected_vars_uncor <- character()
corr_cutoff <- 0.90
for (v in attr_signal$variable) {
  if (attr_signal$low_information[attr_signal$variable == v]) {
    next
  }
  if (length(selected_vars_uncor) == 0) {
    selected_vars_uncor <- c(selected_vars_uncor, v)
    next
  }
  corr_to_selected <- abs(cor_for_selection[v, selected_vars_uncor])
  if (max(corr_to_selected, na.rm = TRUE) < corr_cutoff) {
    selected_vars_uncor <- c(selected_vars_uncor, v)
  } else {
    most_corr <- selected_vars_uncor[which.max(corr_to_selected)]
    attr_signal$redundant_with[attr_signal$variable == v] <- most_corr
  }
}
attr_signal$selected_uncorrelated <- attr_signal$variable %in% selected_vars_uncor

attr_weights <- attr_signal$signal_score
names(attr_weights) <- attr_signal$variable
attr_weights <- attr_weights[sensory_vars]
attr_weights <- attr_weights / sum(attr_weights)
attr_signal$normalized_weight <- as.numeric(attr_weights[attr_signal$variable])
write.csv(attr_signal, file.path(out_dir, "variable_selection_attribute_ranking.csv"), row.names = FALSE)

make_weighted_profile <- function(vars) {
  w <- attr_weights[vars]
  w <- w / sum(w)
  sweep(profile_panel_z[, vars, drop = FALSE], 2, sqrt(w), "*")
}

profile_weighted <- sweep(profile_panel_z, 2, sqrt(attr_weights[colnames(profile_panel_z)]), "*")

ranked_informative_vars <- attr_signal$variable[!attr_signal$low_information]
positive_signal_vars <- attr_signal$variable[
  !attr_signal$low_information &
    attr_signal$signal_score > 1e-6
]
candidate_top_n <- unique(pmin(c(20, 15, 10, 8), length(ranked_informative_vars)))
selection_scenarios <- list(
  selected_signal_positive = positive_signal_vars,
  selected_uncorrelated = selected_vars_uncor
)
for (n_vars in candidate_top_n) {
  selection_scenarios[[paste0("selected_top", n_vars)]] <- head(ranked_informative_vars, n_vars)
}

selection_summary <- data.frame()
for (scenario_name in names(selection_scenarios)) {
  vars <- selection_scenarios[[scenario_name]]
  selection_summary <- rbind(
    selection_summary,
    data.frame(
      scenario = scenario_name,
      n_variables = length(vars),
      variables = paste(vars, collapse = ";"),
      stringsAsFactors = FALSE
    )
  )
}
write.csv(selection_summary, file.path(out_dir, "variable_selection_scenarios.csv"), row.names = FALSE)

representations <- list(
  panelist_block = panel_block,
  panelist_weighted = profile_weighted,
  panelist_z = profile_panel_z,
  raw_scaled = profile_raw
)
for (scenario_name in names(selection_scenarios)) {
  vars <- selection_scenarios[[scenario_name]]
  representations[[scenario_name]] <- make_weighted_profile(vars)
}

# Representacoes adicionais para reduzir ruido mantendo interpretabilidade:
# PCA dos atributos selecionados e MDS a partir de distancia de correlacao.
positive_profile <- make_weighted_profile(positive_signal_vars)
max_reduced_dims <- min(4, nrow(positive_profile) - 1, ncol(positive_profile))
if (max_reduced_dims >= 2) {
  pca_positive <- prcomp(positive_profile, center = TRUE, scale. = FALSE)
  for (n_dim in 2:max_reduced_dims) {
    pca_rep <- pca_positive$x[, seq_len(n_dim), drop = FALSE]
    colnames(pca_rep) <- paste0("PC", seq_len(n_dim))
    representations[[paste0("selected_signal_positive_pca", n_dim)]] <- pca_rep
  }

  corr_products <- cor(t(positive_profile), method = "spearman")
  corr_products[is.na(corr_products)] <- 0
  corr_dist <- as.dist(1 - corr_products)
  mds_corr <- cmdscale(corr_dist, k = max_reduced_dims, eig = TRUE)
  for (n_dim in 2:max_reduced_dims) {
    corr_rep <- mds_corr$points[, seq_len(n_dim), drop = FALSE]
    colnames(corr_rep) <- paste0("CorrMDS", seq_len(n_dim))
    representations[[paste0("selected_signal_positive_corr_mds", n_dim)]] <- corr_rep
  }
}

methods <- c("kmeans", "pam", "hclust_ward")
cluster_results <- data.frame()
cluster_store <- list()
set.seed(123)

for (rep_name in names(representations)) {
  x_rep <- representations[[rep_name]]
  d_rep <- dist(x_rep)
  hc_rep <- hclust(d_rep, method = "ward.D2")
  max_k_rep <- min(5, nrow(x_rep) - 1)

  for (k in 2:max_k_rep) {
    for (m in methods) {
      cl <- fit_clusters(m, k, x_rep, d_rep, hc_rep)
      sil <- tryCatch(mean(cluster::silhouette(cl, d_rep)[, 3]), error = function(e) NA_real_)
      ch <- calc_ch(x_rep, cl)
      db <- calc_davies_bouldin(x_rep, cl)
      dunn <- calc_dunn(x_rep, cl)
      wss <- calc_withinss(x_rep, cl)
      cl_sizes <- as.numeric(table(cl))
      min_size <- min(cl_sizes)
      max_size <- max(cl_sizes)
      balance_ratio <- min_size / max_size

      key <- paste(rep_name, m, k, sep = "__")
      cluster_store[[key]] <- cl

      cluster_results <- rbind(
        cluster_results,
        data.frame(
          representation = rep_name,
          method = m,
          k = k,
          silhouette = sil,
          calinski_harabasz = ch,
          davies_bouldin = db,
          dunn_index = dunn,
          within_ss = wss,
          min_cluster_size = min_size,
          max_cluster_size = max_size,
          balance_ratio = balance_ratio,
          cophenetic_correlation = cor(as.numeric(d_rep), as.numeric(cophenetic(hc_rep))),
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

cluster_results$rank_silhouette <- rank(-cluster_results$silhouette, na.last = "keep", ties.method = "average")
cluster_results$rank_ch <- rank(-cluster_results$calinski_harabasz, na.last = "keep", ties.method = "average")
cluster_results$rank_db <- rank(cluster_results$davies_bouldin, na.last = "keep", ties.method = "average")
cluster_results$rank_dunn <- rank(-cluster_results$dunn_index, na.last = "keep", ties.method = "average")
cluster_results$rank_balance <- rank(-cluster_results$balance_ratio, na.last = "keep", ties.method = "average")
cluster_results$rank_total <- rowMeans(
  cbind(
    cluster_results$rank_silhouette,
    cluster_results$rank_ch,
    cluster_results$rank_db,
    cluster_results$rank_dunn,
    cluster_results$rank_balance
  ),
  na.rm = TRUE
)

cluster_results <- cluster_results[order(cluster_results$rank_total, -cluster_results$silhouette), ]
write.csv(cluster_results, file.path(out_dir, "cluster_model_comparison.csv"), row.names = FALSE)

# Regra de selecao para evitar clusterizacao trivial com grupos muito desbalanceados.
# PCA/MDS por correlacao entram como analise complementar; a solucao final preserva
# atributos sensoriais selecionados diretamente para manter interpretabilidade.
final_pool <- subset(
  cluster_results,
  k >= 3 &
    min_cluster_size >= 2 &
    balance_ratio >= 0.40 &
    !grepl("(_pca|_corr_mds)", representation)
)
sil_cut <- 0.85 * max(final_pool$silhouette, na.rm = TRUE)
candidate_cfg <- subset(final_pool, silhouette >= sil_cut)

if (nrow(candidate_cfg) > 0) {
  candidate_cfg <- candidate_cfg[order(candidate_cfg$rank_total, -candidate_cfg$silhouette), ]
  best_cfg <- candidate_cfg[1, ]
  selection_rule <- "segmentacao_balanceada"
} else if (nrow(final_pool) > 0) {
  final_pool <- final_pool[order(final_pool$rank_total, -final_pool$silhouette), ]
  best_cfg <- final_pool[1, ]
  selection_rule <- "melhor_rank_final_interpretavel"
} else {
  best_cfg <- cluster_results[1, ]
  selection_rule <- "melhor_rank_global"
}

best_key <- paste(best_cfg$representation, best_cfg$method, best_cfg$k, sep = "__")
best_labels <- cluster_store[[best_key]]
x_best <- representations[[best_cfg$representation]]
d_best <- dist(x_best)
hc_best <- hclust(d_best, method = "ward.D2")

cluster_assign <- data.frame(
  Productname = rownames(x_best),
  cluster = best_labels,
  representation = best_cfg$representation,
  method = best_cfg$method,
  k = best_cfg$k
)
cluster_assign <- cluster_assign[order(cluster_assign$cluster, cluster_assign$Productname), ]
write.csv(cluster_assign, file.path(out_dir, "product_clusters_best.csv"), row.names = FALSE)

cat("Top 10 configuracoes de clusterizacao (ranking composto):\n")
print(head(cluster_results, 10))
cat("\nMelhor configuracao:\n")
print(best_cfg)
cat("Regra de selecao:", selection_rule, "\n")
cat("\n")

# Gap statistic para cada representacao (kmeans)
gap_all <- data.frame()
for (rep_name in names(representations)) {
  x_rep <- representations[[rep_name]]
  max_k_rep <- min(5, nrow(x_rep) - 1)
  set.seed(123)
  gap <- cluster::clusGap(x_rep, FUN = kmeans, K.max = max_k_rep, B = 200, nstart = 50, iter.max = 100)
  gap_df <- data.frame(
    representation = rep_name,
    k = seq_len(max_k_rep),
    logW = gap$Tab[, "logW"],
    gap = gap$Tab[, "gap"],
    SE.sim = gap$Tab[, "SE.sim"],
    stringsAsFactors = FALSE
  )
  gap_all <- rbind(gap_all, gap_df)
}
write.csv(gap_all, file.path(out_dir, "gap_statistic_by_representation.csv"), row.names = FALSE)

# Estabilidade bootstrap da melhor configuracao (re-amostrando atributos)
set.seed(123)
n_boot <- 300
ari_boot <- numeric(n_boot)
coassoc <- matrix(0, nrow = nrow(x_best), ncol = nrow(x_best))
rownames(coassoc) <- rownames(x_best)
colnames(coassoc) <- rownames(x_best)

for (b in seq_len(n_boot)) {
  cols <- sample(seq_len(ncol(x_best)), size = ncol(x_best), replace = TRUE)
  xb <- x_best[, cols, drop = FALSE]
  db <- dist(xb)
  hb <- hclust(db, method = "ward.D2")
  cb <- fit_clusters(best_cfg$method, best_cfg$k, xb, db, hb)
  ari_boot[b] <- calc_ari(best_labels, cb)
  coassoc <- coassoc + outer(cb, cb, FUN = "==")
}

stability_summary <- data.frame(
  bootstrap_runs = n_boot,
  ari_mean = mean(ari_boot, na.rm = TRUE),
  ari_sd = sd(ari_boot, na.rm = TRUE),
  ari_q05 = unname(quantile(ari_boot, 0.05, na.rm = TRUE)),
  ari_median = median(ari_boot, na.rm = TRUE),
  ari_q95 = unname(quantile(ari_boot, 0.95, na.rm = TRUE)),
  stringsAsFactors = FALSE
)
write.csv(stability_summary, file.path(out_dir, "cluster_stability_bootstrap.csv"), row.names = FALSE)

coassoc <- coassoc / n_boot
write.csv(coassoc, file.path(out_dir, "cluster_consensus_matrix.csv"))

cluster_metrics <- data.frame(
  metric = c(
    "best_silhouette",
    "best_calinski_harabasz",
    "best_davies_bouldin",
    "best_dunn_index",
    "best_cophenetic_correlation",
    "bootstrap_ari_mean",
    "bootstrap_ari_q05",
    "bootstrap_ari_q95"
  ),
  value = c(
    best_cfg$silhouette,
    best_cfg$calinski_harabasz,
    best_cfg$davies_bouldin,
    best_cfg$dunn_index,
    best_cfg$cophenetic_correlation,
    stability_summary$ari_mean,
    stability_summary$ari_q05,
    stability_summary$ari_q95
  ),
  stringsAsFactors = FALSE
)
write.csv(cluster_metrics, file.path(out_dir, "cluster_quality_metrics.csv"), row.names = FALSE)

cat("Metricas finais da melhor configuracao:\n")
print(cluster_metrics)
cat("\nResumo de estabilidade bootstrap:\n")
print(stability_summary)
cat("\n")

# Visualizacao 1: dendrograma da melhor representacao
png(file.path(out_dir, "cluster_dendrogram_best.png"), width = 1100, height = 750)
plot(
  hc_best,
  main = paste("Dendrograma dos produtos -", best_cfg$representation),
  xlab = "Produtos",
  sub = paste("Metodo: Ward | k =", best_cfg$k),
  hang = -1,
  cex = 0.95
)
rect.hclust(hc_best, k = best_cfg$k, border = cluster_palette[seq_len(best_cfg$k)])
dev.off()

# Visualizacao 2: espalhamento no PCA da melhor representacao
pca_best <- prcomp(x_best, center = TRUE, scale. = TRUE)
pca_best_var <- pca_best$sdev^2 / sum(pca_best$sdev^2)
png(file.path(out_dir, "cluster_best_pca_scatter.png"), width = 1100, height = 750)
plot(
  pca_best$x[, 1], pca_best$x[, 2],
  col = cluster_palette[best_labels], pch = 19, cex = 1.8,
  xlab = paste0("PC1 (", round(100 * pca_best_var[1], 1), "%)"),
  ylab = paste0("PC2 (", round(100 * pca_best_var[2], 1), "%)"),
  main = paste("Produtos por cluster -", best_cfg$representation),
  sub = paste("Metodo:", best_cfg$method, "| k =", best_cfg$k)
)
abline(h = 0, v = 0, col = "gray85", lty = 2)
text(pca_best$x[, 1], pca_best$x[, 2], labels = rownames(x_best), pos = 3, cex = 0.9)
add_cluster_legend(best_labels)
dev.off()

# Visualizacao 3: silhouette da melhor configuracao
sil_best <- cluster::silhouette(best_labels, d_best)
png(file.path(out_dir, "cluster_best_silhouette.png"), width = 1000, height = 700)
plot(
  sil_best,
  main = paste("Silhouette por produto - media =", round(mean(sil_best[, 3]), 3)),
  col = cluster_palette[seq_len(best_cfg$k)],
  border = NA,
  cex.names = 0.85
)
abline(v = mean(sil_best[, 3]), col = "#D95F02", lwd = 2, lty = 2)
dev.off()

# Visualizacao 4: consenso de cluster (co-association)
png(file.path(out_dir, "cluster_consensus_heatmap.png"), width = 900, height = 800)
ord_cons <- order(best_labels, rownames(coassoc))
plot_matrix_heatmap(
  coassoc[ord_cons, ord_cons, drop = FALSE],
  main = "Consenso bootstrap - proporcao de vezes no mesmo cluster",
  low = "white",
  mid = "#B3CDE3",
  high = "#005B96"
)
dev.off()

# =========================
# 6) MDS (mapa de proximidade)
# =========================
mds <- cmdscale(d_best, k = 2, eig = TRUE)
mds_df <- data.frame(
  Productname = rownames(x_best),
  Dim1 = mds$points[, 1],
  Dim2 = mds$points[, 2],
  cluster = best_labels
)
write.csv(mds_df, file.path(out_dir, "mds_coordinates.csv"), row.names = FALSE)

png(file.path(out_dir, "mds_products.png"), width = 1100, height = 750)
plot(
  mds_df$Dim1, mds_df$Dim2,
  pch = 19, cex = 1.8, col = cluster_palette[mds_df$cluster],
  xlab = "Dim1", ylab = "Dim2",
  main = "MDS - proximidade sensorial entre produtos",
  sub = "Produtos proximos no grafico possuem perfis sensoriais mais semelhantes"
)
abline(h = 0, v = 0, col = "gray90", lty = 2)
text(mds_df$Dim1, mds_df$Dim2, labels = mds_df$Productname, pos = 3, cex = 0.9)
add_cluster_legend(mds_df$cluster)
dev.off()

# =========================
# 7) Interpretacao dos clusters
# =========================
cluster_factor <- factor(best_labels, levels = sort(unique(best_labels)))
profile_best <- product_profile_z
profile_best$cluster <- cluster_factor

cluster_profile <- aggregate(profile_best[, sensory_vars], by = list(cluster = profile_best$cluster), FUN = mean)
write.csv(cluster_profile, file.path(out_dir, "cluster_profile_panelist_z.csv"), row.names = FALSE)

cluster_profile_mat <- as.matrix(cluster_profile[, positive_signal_vars, drop = FALSE])
rownames(cluster_profile_mat) <- paste("Cluster", cluster_profile$cluster)
png(file.path(out_dir, "cluster_profile_heatmap.png"), width = 1100, height = 650)
plot_matrix_heatmap(
  cluster_profile_mat,
  main = "Perfil sensorial medio dos clusters (z-score por avaliador)"
)
dev.off()

anova_summary <- data.frame(
  attribute = sensory_vars,
  p_value_anova = NA_real_,
  p_value_kruskal = NA_real_,
  eta2 = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(sensory_vars)) {
  v <- sensory_vars[i]
  y <- profile_best[[v]]
  fit <- aov(y ~ cluster_factor)
  fit_sum <- summary(fit)[[1]]
  ss_between <- fit_sum[1, "Sum Sq"]
  ss_total <- sum(fit_sum[, "Sum Sq"])
  anova_summary$p_value_anova[i] <- fit_sum[1, "Pr(>F)"]
  anova_summary$p_value_kruskal[i] <- kruskal.test(y ~ cluster_factor)$p.value
  anova_summary$eta2[i] <- ss_between / ss_total
}

anova_summary <- anova_summary[order(-anova_summary$eta2, anova_summary$p_value_kruskal), ]
write.csv(anova_summary, file.path(out_dir, "cluster_attribute_significance.csv"), row.names = FALSE)

# =========================
# 8) Correlacoes entre atributos
# =========================
cor_mat <- cor(df[sensory_vars], method = "spearman")
write.csv(cor_mat, file.path(out_dir, "sensory_correlation_spearman.csv"))

pair_df <- data.frame(
  var1 = character(),
  var2 = character(),
  rho = numeric(),
  stringsAsFactors = FALSE
)
for (i in 1:(ncol(cor_mat) - 1)) {
  for (j in (i + 1):ncol(cor_mat)) {
    pair_df <- rbind(pair_df, data.frame(var1 = colnames(cor_mat)[i], var2 = colnames(cor_mat)[j], rho = cor_mat[i, j]))
  }
}
pair_df <- pair_df[order(-abs(pair_df$rho)), ]
write.csv(pair_df, file.path(out_dir, "top_correlations_ranked.csv"), row.names = FALSE)

png(file.path(out_dir, "correlation_heatmap.png"), width = 1000, height = 900)
cor_ord <- hclust(as.dist(1 - abs(cor_mat)))$order
plot_matrix_heatmap(
  cor_mat[cor_ord, cor_ord, drop = FALSE],
  main = "Correlacao de Spearman entre atributos sensoriais"
)
dev.off()

cat("=============================================================\n")
cat("RESULTADO FINAL - CLUSTERIZACAO NAO SUPERVISIONADA\n")
cat("=============================================================\n")
cat("Melhor configuracao:", best_cfg$representation, "|", best_cfg$method, "| k =", best_cfg$k, "\n")
cat("Silhouette:", round(best_cfg$silhouette, 4),
    "| CH:", round(best_cfg$calinski_harabasz, 4),
    "| DB:", round(best_cfg$davies_bouldin, 4),
    "| Dunn:", round(best_cfg$dunn_index, 4), "\n")
cat("Estabilidade bootstrap (ARI): media =", round(stability_summary$ari_mean, 4),
    "| q05 =", round(stability_summary$ari_q05, 4),
    "| q95 =", round(stability_summary$ari_q95, 4), "\n")
cat("=============================================================\n")
cat("Exploracao concluida.\n")
cat("Arquivos gerados em:", normalizePath(out_dir), "\n")
