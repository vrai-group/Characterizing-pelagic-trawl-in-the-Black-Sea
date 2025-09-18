library(data.table)
library(igraph)
library(ggplot2)

# --- 0) ROBUST PREP: rebuild dy_all with the correct columns ---

# If you already have dyadic_results1 (MMSI1, MMSI2, T, Proximity, DIh, DId, DItheta, ...)
stopifnot(exists("dyadic_results1"))

dy_all <- as.data.table(dyadic_results1)

# Normalize column names (to handle potential inconsistencies)
cn <- names(dy_all)
# Trim whitespace
setnames(dy_all, cn, trimws(cn))

# Check required columns; map alternative names here if you have any
req <- c("MMSI1","MMSI2","T","Proximity","DIh","DId","DItheta")
if (!all(req %in% names(dy_all))) {
  stop("Missing columns: ", paste(setdiff(req, names(dy_all)), collapse=", "))
}

# Coerce types to numeric/uniform
dy_all[, `:=`(
  MMSI1 = as.numeric(MMSI1),
  MMSI2 = as.numeric(MMSI2),
  T = as.numeric(T),
  Proximity = as.numeric(Proximity),
  DIh = as.numeric(DIh),
  DId = as.numeric(DId),
  DItheta = as.numeric(DItheta)
)]

# --- 1) HELPER FUNCTIONS (robust version) ---

edge_set <- function(g) {
  if (!inherits(g, "igraph") || igraph::gsize(g) == 0) return(character(0))
  el <- igraph::as_edgelist(g, names = TRUE)
  apply(el, 1, function(x) paste(sort(x), collapse = "-"))
}

deg_spearman_vs_base <- function(g, base_deg) {
  # base_deg: named numeric vector (degree), names = node ids (char)
  dg <- igraph::degree(g)
  all_nodes <- union(names(base_deg), names(dg))
  v0 <- base_deg[all_nodes]; v0[is.na(v0)] <- 0
  v1 <-  dg[all_nodes];      v1[is.na(v1)] <- 0
  if (length(unique(v0)) < 2 || length(unique(v1)) < 2) return(NA_real_)
  suppressWarnings(as.numeric(cor(v0, v1, method = "spearman")))
}

build_graph_from_thresholds <- function(dy, prox=0.5, dih=0.5, did=0.5, dith=0.5, min_T=0) {
  dy <- as.data.table(dy)
  req <- c("MMSI1","MMSI2","T","Proximity","DIh","DId","DItheta")
  if (!all(req %in% names(dy))) {
    stop("Required columns missing inside build_graph_from_thresholds(): ",
         paste(setdiff(req, names(dy)), collapse=", "))
  }
  keep <- dy[
    is.finite(T) & T >= min_T &
    is.finite(Proximity) & Proximity >= prox &
    is.finite(DIh)       & DIh       >= dih &
    is.finite(DId)       & DId       >= did &
    is.finite(DItheta)   & DItheta   >= dith
  ]
  if (nrow(keep) == 0L) {
    return(list(
      edges = keep,
      g     = igraph::make_empty_graph(),
      summary = data.table(
        n_edges=0L, n_nodes=0L, avg_deg=0, density=0, n_comp=0L, giant_frac=0
      )
    ))
  }
  g <- igraph::graph_from_data_frame(keep[, .(from = as.character(MMSI1),
                                              to   = as.character(MMSI2))],
                                     directed = FALSE)
  n_edges <- igraph::gsize(g)
  n_nodes <- igraph::gorder(g)
  avg_deg <- mean(igraph::degree(g))
  dens    <- igraph::edge_density(g)
  cmp     <- igraph::components(g)
  giant_frac <- if (n_nodes > 0) max(cmp$csize)/n_nodes else 0

  list(
    edges = keep,
    g     = g,
    summary = data.table(
      n_edges = n_edges,
      n_nodes = n_nodes,
      avg_deg = avg_deg,
      density = dens,
      n_comp  = cmp$no,
      giant_frac = giant_frac
    )
  )
}

# --- 2) BUILD BASELINE GRAPH ---
base_th <- list(prox=0.5, dih=0.5, did=0.5, dith=0.5, min_T=0)
base_res <- build_graph_from_thresholds(
  dy_all,
  prox = base_th$prox, dih = base_th$dih,
  did  = base_th$did,  dith= base_th$dith,
  min_T= base_th$min_T
)

base_edges <- edge_set(base_res$g)
base_deg   <- igraph::degree(base_res$g)  # named vector

# --- 3) SINGLE-METRIC SWEEP (0.40–0.60) ---
sweep_threshold_one <- function(dy, base_edges, base_deg,
                                metric = c("Proximity","DIh","DId","DItheta"),
                                grid   = seq(0.40, 0.60, by=0.05),
                                hold   = list(Proximity=0.50, DIh=0.50, DId=0.50, DItheta=0.50),
                                min_T  = 0) {
  metric <- match.arg(metric)
  out <- vector("list", length(grid))
  for (i in seq_along(grid)) {
    th <- hold; th[[metric]] <- grid[i]
    res <- build_graph_from_thresholds(dy,
              prox=th$Proximity, dih=th$DIh, did=th$DId, dith=th$DItheta, min_T=min_T)
    es <- edge_set(res$g)
    jacc <- if ((length(base_edges) + length(es)) == 0) NA_real_ else {
      inter <- length(intersect(base_edges, es))
      uni   <- length(union(base_edges, es))
      if (uni == 0) NA_real_ else inter/uni
    }
    rho  <- deg_spearman_vs_base(res$g, base_deg)
    out[[i]] <- data.table(
      metric = metric, thr = grid[i],
      n_edges = res$summary$n_edges,
      density = res$summary$density,
      giant_frac = res$summary$giant_frac,
      edge_jaccard_vs_base = jacc,
      deg_spearman_vs_base = rho
    )
  }
  rbindlist(out, use.names=TRUE, fill=TRUE)
}

sens_prox <- sweep_threshold_one(dy_all, base_edges, base_deg, metric="Proximity")
sens_dih  <- sweep_threshold_one(dy_all, base_edges, base_deg, metric="DIh")
sens_did  <- sweep_threshold_one(dy_all, base_edges, base_deg, metric="DId")
sens_dith <- sweep_threshold_one(dy_all, base_edges, base_deg, metric="DItheta")

# --- 4) SUMMARY TABLE ---
summarize_sweep <- function(dt, name){
  data.frame(
    metric      = name,
    edges_med   = median(dt$n_edges, na.rm=TRUE),
    edges_min   = min(dt$n_edges, na.rm=TRUE),
    edges_max   = max(dt$n_edges, na.rm=TRUE),
    dens_med    = median(dt$density, na.rm=TRUE),
    giant_med   = median(dt$giant_frac, na.rm=TRUE),
    jacc_med    = median(dt$edge_jaccard_vs_base, na.rm=TRUE),
    jacc_min    = min(dt$edge_jaccard_vs_base, na.rm=TRUE),
    jacc_max    = max(dt$edge_jaccard_vs_base, na.rm=TRUE),
    deg_rho_med = median(dt$deg_spearman_vs_base, na.rm=TRUE)
  )
}
robust_tbl <- rbind(
  summarize_sweep(sens_prox, "Proximity"),
  summarize_sweep(sens_dih,  "DIh"),
  summarize_sweep(sens_did,  "DId"),
  summarize_sweep(sens_dith, "DItheta")
)
print(robust_tbl)

# --- 5) EDGE TURNOVER (Proximity example) ---
baseline_edges_n <- nrow(as.data.frame(igraph::as_edgelist(base_res$g)))

sens_prox[, delta_edges_pct := 100 * (n_edges - baseline_edges_n) / baseline_edges_n]

edge_turnover <- function(dy, thr, hold, base_edges, min_T=0){
  th <- hold; th$Proximity <- thr
  res <- build_graph_from_thresholds(dy,
          prox=th$Proximity, dih=th$DIh, did=th$DId, dith=th$DItheta, min_T=min_T)
  es  <- edge_set(res$g)
  data.table(
    thr    = thr,
    gains  = length(setdiff(es, base_edges)),
    losses = length(setdiff(base_edges, es))
  )
}

turn_tbl <- rbindlist(
  lapply(sens_prox$thr, edge_turnover,
         dy = dy_all,
         hold = list(Proximity=0.5, DIh=0.5, DId=0.5, DItheta=0.5),
         base_edges = base_edges),
  use.names=TRUE, fill=TRUE
)
print(turn_tbl)

library(data.table)
library(igraph)
library(ggplot2)

# --- 1) Node retention in Proximity sweep ---
base_nodes <- V(base_res$g)$name

node_retention <- function(dy, thr, hold, min_T=0){
  th <- hold; th$Proximity <- thr
  res <- build_graph_from_thresholds(dy, prox=th$Proximity, dih=th$DIh, did=th$DId, dith=th$DItheta, min_T=min_T)
  nodes <- V(res$g)$name
  data.table(thr=thr,
             node_retention = ifelse(length(base_nodes)==0, NA_real_,
                                     length(intersect(nodes, base_nodes))/length(base_nodes)))
}
ret_tbl <- rbindlist(lapply(sens_prox$thr, node_retention,
                            dy=dy_all, hold=list(Proximity=0.5,DIh=0.5,DId=0.5,DItheta=0.5)))
print(ret_tbl)

# --- 2) Core edges in the Proximity sweep (edges appearing in ≥80% of thresholds) ---
edge_set <- function(g){
  if (!inherits(g,"igraph") || gsize(g)==0) return(character(0))
  apply(as_edgelist(g, names=TRUE), 1, function(x) paste(sort(x), collapse="-"))
}
# collect edge sets for all thresholds
edge_lists <- lapply(sens_prox$thr, function(th){
  res <- build_graph_from_thresholds(dy_all, prox=th, dih=0.5, did=0.5, dith=0.5, min_T=0)
  edge_set(res$g)
})
# count frequencies
edge_freq <- table(unlist(edge_lists))
thr_n <- length(edge_lists)
core_edges <- names(edge_freq)[edge_freq >= ceiling(0.8*thr_n)]
cat("Core edges (>=80% of thresholds):", length(core_edges), "\n")

# --- 3) Characterization of lost/gained edges (Proximity=0.55, 0.60) ---
edge_key <- function(a,b) paste(pmin(a,b), pmax(a,b), sep="-")
# collapse base metrics to edge level
dy_edge <- copy(dy_all)[, edge:=edge_key(MMSI1,MMSI2)]
dy_edge <- unique(dy_edge, by="edge")  # if dy_all is already summarized rows, this is sufficient

edge_set_at <- function(th){
  res <- build_graph_from_thresholds(dy_all, prox=th, dih=0.5, did=0.5, dith=0.5, min_T=0)
  edge_set(res$g)
}
es_base <- edge_set_at(0.50)
for (thr in c(0.55, 0.60)) {
  es_thr <- edge_set_at(thr)
  lost   <- setdiff(es_base, es_thr)
  gained <- setdiff(es_thr, es_base)
  cat(sprintf("\nProximity=%.2f  lost=%d gained=%d\n", thr, length(lost), length(gained)))
  if (length(lost)) {
    print(dy_edge[edge %in% lost, .(
      n=.N,
      prox_med = median(Proximity, na.rm=TRUE),
      dih_med  = median(DIh, na.rm=TRUE),
      did_med  = median(DId, na.rm=TRUE),
      dith_med = median(DItheta, na.rm=TRUE)
    )])
  }
}

# --- 4) Quick plots: threshold → Jaccard and threshold → density ---
sens_all <- rbindlist(list(sens_prox, sens_dih, sens_did, sens_dith), fill=TRUE)
sens_all[, metric := factor(metric, levels=c("Proximity","DIh","DId","DItheta"))]

ggplot(sens_all, aes(thr, edge_jaccard_vs_base)) +
  geom_line() + geom_point() +
  facet_wrap(~metric, ncol=2, scales="free_y") +
  labs(x="Threshold", y="Jaccard vs baseline", title="Threshold → Jaccard")

ggplot(sens_all, aes(thr, density)) +
  geom_line() + geom_point() +
  facet_wrap(~metric, ncol=2) +
  labs(x="Threshold", y="Network density", title="Threshold → Density")