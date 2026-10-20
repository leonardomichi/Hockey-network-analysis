# =============================================================================
#  Rete degli assist - Hockey Empoli "Flying Donkeys", Serie B 2025/26
#  Dati: tabellini ufficiali FISR (hockeyinline.fisr.it), 22 partite.
#  Grafo orientato: arco da chi fa l'assist a chi segna, peso = numero di gol.
#
#  File richiesti nella stessa cartella:
#    archi_assist_empoli_2526.csv   from, to, weight
#    nodi_empoli_2526.csv           attributi dei giocatori
#    gol_empoli_2526.csv            un gol per riga (marcatore, assist, partita, ...)
#    presenze_empoli_2526.csv       chi era a referto in ogni partita
#
#  Uso:  Rscript analisi_rete_empoli.R      (pacchetto: igraph)
# =============================================================================

suppressPackageStartupMessages(library(igraph))
set.seed(2526)
setwd("/Users/leonardomichi/Desktop/instagram hockey")
nodi     <- read.csv("nodi_empoli_2526.csv", stringsAsFactors = FALSE, encoding = "UTF-8")
archi    <- read.csv("archi_assist_empoli_2526.csv", stringsAsFactors = FALSE, encoding = "UTF-8")
gol      <- read.csv("gol_empoli_2526.csv", stringsAsFactors = FALSE, encoding = "UTF-8")
presenze <- read.csv("presenze_empoli_2526.csv", stringsAsFactors = FALSE, encoding = "UTF-8")

# hub_score()/authority_score() sono deprecate da igraph 2.1 a favore di hits_scores()
hits <- function(g, w) {
  if (exists("hits_scores", where = asNamespace("igraph"), mode = "function")) {
    h <- igraph::hits_scores(g, weights = w)
    list(hub = h$hub, authority = h$authority)
  } else {
    list(hub = hub_score(g, weights = w)$vector, authority = authority_score(g, weights = w)$vector)
  }
}

# -----------------------------------------------------------------------------
# 1. Il grafo
# -----------------------------------------------------------------------------
# Nodi: giocatori di movimento con almeno un legame (i portieri non hanno gol né assist)
legati <- unique(c(archi$from, archi$to))
v_att  <- subset(nodi, portiere == 0 & name %in% legati)
g  <- graph_from_data_frame(archi, directed = TRUE, vertices = v_att)
w  <- E(g)$weight
# proiezione non orientata: peso = gol in comune nelle due direzioni
gu <- as.undirected(g, mode = "collapse", edge.attr.comb = list(weight = "sum"))

cat("\n== Rete ==\n")
cat(sprintf("Nodi %d, archi %d, gol con assist %d\n", vcount(g), ecount(g), sum(w)))
cat(sprintf("Densità %.3f | reciprocità %.3f | transitività %.3f\n",
            edge_density(g), reciprocity(g), transitivity(gu, type = "global")))

# -----------------------------------------------------------------------------
# 2. Centralità
# -----------------------------------------------------------------------------
h <- hits(g, w)
centr <- data.frame(
  giocatore      = V(g)$name,
  ruolo          = V(g)$ruolo_rosa_2627,
  assist         = strength(g, mode = "out", weights = w),
  gol_su_assist  = strength(g, mode = "in", weights = w),
  serve_n        = degree(g, mode = "out"),
  servito_da_n   = degree(g, mode = "in"),
  connessioni    = degree(gu),
  coinvolgimenti = strength(gu, weights = E(gu)$weight),
  hub            = h$hub / max(h$hub),
  authority      = h$authority / max(h$authority),
  betweenness    = betweenness(g, weights = 1 / w, directed = TRUE, normalized = TRUE),
  autovettore    = eigen_centrality(gu, weights = E(gu)$weight)$vector,
  pagerank       = page_rank(g, weights = w)$vector,
  row.names = NULL
)
centr <- centr[order(-centr$coinvolgimenti), ]
num <- sapply(centr, is.numeric)
centr[num] <- lapply(centr[num], function(x) round(x, 3))
cat("\n== Centralità (hub, authority e autovettore: 1 = il primo) ==\n")
print(centr, row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Comunità (Louvain sulla rete non orientata pesata)
# -----------------------------------------------------------------------------
# Louvain dipende dall'ordine dei nodi: ripeto su permutazioni casuali e tengo la modularità massima
best <- NULL
for (k in 1:200) {
  perm <- sample(vcount(gu))
  gp <- permute(gu, perm)
  cl <- cluster_louvain(gp, weights = E(gp)$weight)
  memb <- membership(cl)[perm]          # riporta le etichette all'ordine originale
  q <- modularity(gu, memb, weights = E(gu)$weight)
  if (is.null(best) || q > best$q + 1e-12) best <- list(q = q, memb = memb)
}
cat(sprintf("\n== Comunità: modularità massima %.3f ==\n", best$q))
for (k in sort(unique(best$memb))) cat(sprintf("Gruppo %d: %s\n", k, paste(V(gu)$name[best$memb == k], collapse = ", ")))

# -----------------------------------------------------------------------------
# 4. Modello nullo: riassegnazione casuale degli assist
#    Si scambiano gli assistenti di due gol della stagione se
#    (a) ciascuno era a referto nella partita dell'altro gol e
#    (b) nessuno diventa assist di se stesso.
#    Ogni giocatore conserva i suoi assist e i suoi gol su assist.
# -----------------------------------------------------------------------------
ga <- subset(gol, assist != "")
s_ <- ga$marcatore; a_ <- ga$assist; partita <- as.character(ga$idp)
a_ref <- split(presenze$giocatore, presenze$idp)
n <- length(a_)
swap_chain <- function(a, steps) {
  for (k in seq_len(steps)) {
    i <- sample.int(n, 1); j <- sample.int(n, 1)
    ai <- a[i]; aj <- a[j]
    if (i == j || ai == aj || aj == s_[i] || ai == s_[j]) next
    if (!(aj %in% a_ref[[partita[i]]]) || !(ai %in% a_ref[[partita[j]]])) next
    a[i] <- aj; a[j] <- ai
  }
  a
}
B <- 1000
lev <- V(g)$name
obs <- table(factor(a_, lev), factor(s_, lev))
a_cur <- swap_chain(a_, 20000)
n_coppie <- hhi <- numeric(B)
W_null <- array(0L, c(B, length(lev), length(lev)), dimnames = list(NULL, lev, lev))
for (b in 1:B) {
  a_cur <- swap_chain(a_cur, 300)
  tb <- table(factor(a_cur, lev), factor(s_, lev))
  W_null[b, , ] <- tb
  n_coppie[b] <- sum(tb > 0)
  hhi[b] <- sum(tb^2) / sum(tb)^2
}
stopifnot(all(rowSums(W_null[B, , ]) == rowSums(obs)), all(colSums(W_null[B, , ]) == colSums(obs)))
cat("\n== Modello nullo (", B, " reti casuali) ==\n", sep = "")
cat(sprintf("Coppie diverse: reali %d, attese %.1f (sd %.1f), p = %.3f\n",
            sum(obs > 0), mean(n_coppie), sd(n_coppie), (1 + sum(n_coppie <= sum(obs > 0))) / (1 + B)))
cat(sprintf("Concentrazione (HHI): reale %.4f, attesa %.4f, p = %.4f\n",
            sum(obs^2) / sum(obs)^2, mean(hhi), (1 + sum(hhi >= sum(obs^2) / sum(obs)^2)) / (1 + B)))
coppia <- function(da, a) {
  x <- W_null[, da, a]
  cat(sprintf("%-11s -> %-11s reali %d, attesi %.1f, p(>=) = %.3f\n", da, a, obs[da, a], mean(x),
              (1 + sum(x >= obs[da, a])) / (1 + B)))
}
coppia("Carboncini", "Errico"); coppia("Carboncini", "Michi"); coppia("Cecchi", "Carboncini")
und <- W_null[, "Carboncini", "Simonelli"] + W_null[, "Simonelli", "Carboncini"]
cat(sprintf("Carboncini <-> Simonelli reali %d, attesi %.1f, p(<=) = %.3f\n",
            obs["Carboncini", "Simonelli"] + obs["Simonelli", "Carboncini"], mean(und),
            (1 + sum(und <= obs["Carboncini", "Simonelli"] + obs["Simonelli", "Carboncini"])) / (1 + B)))

# -----------------------------------------------------------------------------
# 5. Bootstrap sulle partite: quanto è solido il primo posto?
# -----------------------------------------------------------------------------
partite <- unique(gol$idp)
Bb <- 1000
primo <- setNames(numeric(length(lev)), lev)
for (b in 1:Bb) {
  pick <- sample(partite, length(partite), replace = TRUE)
  gb <- do.call(rbind, lapply(pick, function(p) ga[ga$idp == p, c("assist", "marcatore")]))
  colnames(gb) <- c("from", "to")
  gbb <- graph_from_data_frame(gb, directed = TRUE, vertices = data.frame(name = lev))
  E(gbb)$weight <- 1
  gub <- as.undirected(simplify(gbb, edge.attr.comb = list(weight = "sum")), mode = "collapse",
                       edge.attr.comb = list(weight = "sum"))
  sb <- strength(gub, weights = E(gub)$weight)
  top <- names(sb)[sb == max(sb)]
  primo[top] <- primo[top] + 1 / length(top)
}
cat("\n== Bootstrap (", Bb, " stagioni): quota di volte primo per coinvolgimenti ==\n", sep = "")
print(round(sort(primo[primo > 0] / Bb, decreasing = TRUE), 3))

# -----------------------------------------------------------------------------
# 6. Grafico
# -----------------------------------------------------------------------------
pal <- c("#B8893F", "#2643B8", "#4E8062", "#888888", "#AA5544")
png("rete_assist_empoli.png", width = 1600, height = 1300, res = 160)
par(mar = c(0, 0, 2, 0))
lay <- layout_with_fr(gu, weights = E(gu)$weight, niter = 2000)
plot(g, layout = lay,
     vertex.color = pal[best$memb], vertex.frame.color = "white",
     vertex.size = 6 + 2.2 * sqrt(strength(gu, weights = E(gu)$weight)),
     vertex.label.color = "#0E1642", vertex.label.cex = 0.8, vertex.label.dist = 1.6,
     vertex.label.family = "sans",
     edge.width = 0.6 + 0.7 * E(g)$weight, edge.arrow.size = 0.35,
     edge.color = adjustcolor("#0E1642", alpha.f = 0.45),
     edge.curved = 0.18)
title("Hockey Empoli 2025/26 - rete degli assist (assist -> gol)", cex.main = 0.95)
invisible(dev.off())
cat("\nGrafico salvato in rete_assist_empoli.png\n")

