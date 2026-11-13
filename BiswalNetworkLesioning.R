library(igraphdata)
library(igraph)


# pre-processing; threshold and clean-up
matrix <- read.delim('AveragedMatrix.txt', header = FALSE, sep = ' ')
short_labels <- read.delim('region_names_abbrev_file.txt', header = FALSE)
long_labels <- read.delim('region_names_full_file.txt', header = FALSE)
matrix <- abs(matrix)
quantile(c(as.matrix(matrix)), probs = seq(0, 1, 0.05))
matrix[matrix < 0.174] <- 0 # keeping the top 20% of all edges by weight
matrix[matrix > 0] <- 1 # create an unweighted, undirected matrix
matrix <- as.matrix(matrix)
colnames(matrix) <- short_labels$V1
rownames(matrix) <- short_labels$V1
merged <- data.frame(short_label = short_labels$V1, long_label = long_labels$V1)
head(merged)

# conversion to network object
brain_uu <- graph_from_adjacency_matrix(matrix, mode = 'undirected', weighted = NA) 
summary(brain_uu)
edge_density(brain_uu)

# exploring density of the brain network 
degree(brain_uu) |> mean() 
mean_distance(brain_uu)
transitivity(brain_uu, type = "global") 
edge_density(brain_uu)
diameter(brain_uu, directed = F) 
assortativity_degree(brain_uu, directed = F)

# erdos-renyi network generation and checking small-world network structure
set.seed(2)
er_graphs <- list()
for(i in 1:1000)
  er_graphs[[i]] <- sample_gnm(n = gorder(brain_uu), m = gsize(brain_uu), directed = F, loops = F)

sapply(er_graphs, mean_distance) |> mean() 
sapply(er_graphs, transitivity, type = 'global') |> mean() 
(0.4973652/0.2011473)/(2.00642/1.7992676) # small-world index value 

# detect louvain communities and obtain q of the network
# generating louvain clusters 1000 times
cluster_results <- list()

for(i in 1:1000) { 
  set.seed(i)
  
  cluster_results[[i]] <- cluster_louvain(brain_uu)
}

# compare Q to check consistency of results
sapply(cluster_results, modularity)

df1 <- lapply(
  cluster_results,
  function(x) membership(x)
)

df1_merged <- do.call(cbind, df1)

# save highest Q
brain_louvain <- cluster_results[[
  which.max(sapply(cluster_results, modularity))
]]

# build composite hub measure
braincomm <- membership(brain_louvain)

# compute participation coefficient of hub
PC <- sapply(seq_len(vcount(brain_uu)), function(i) {
  
  # identify neighbours of node 
  neighbours_i <- as.integer(neighbors(brain_uu, i))
  
  # total degree
  k <- length(neighbours_i)
  
  if (k == 0) return(0)
  
  # count connecting edges to each detected community
  k_is <- table(braincomm[neighbours_i])
  
  # calculate participation coefficient
  1 - sum((k_is / k)^2)
})

hub_measures <- data.frame(
  short_label = V(brain_uu)$name,
  degree = degree(brain_uu),
  
  betweenness = betweenness(brain_uu,
                            weights = NA,
                            directed = FALSE,
                            normalized = TRUE
  ),
  participation = PC
  )

# rank hubs by each measure
hub_measures$degree_rank <- rank(-hub_measures$degree)

hub_measures$betweenness_rank <- rank(-hub_measures$betweenness)

hub_measures$participation_rank <- rank(-hub_measures$participation)

# average rank across all three measures to obtain composite score
hub_measures$hub_score <- rowMeans(
  hub_measures[, c(
    "degree_rank",
    "betweenness_rank",
    "participation_rank"
  )]
)

# ranking of nodes by hubness score 
hub_ranked <- hub_measures[
  order(hub_measures$hub_score),
]

# select top 10% of all nodes to qualify hubs
top_hubs <- head(
  hub_ranked,
  ceiling(0.10 * nrow(hub_ranked))
)

top_hubs$long_label <- merged$long_label[
  match(top_hubs$short_label, merged$short_label)
]

top_hubs

# outcome metrics and removal
network_metrics <- function(g) {
  
  Q_results <- numeric(1000)

  for (i in 1:1000) {
    set.seed(i)
    Q_results[i] <- modularity(
      cluster_louvain(g, weights = NA)
    )
  }

  Q <- mean(Q_results)
  
  # average shortest path length (ASPL)
  L <- mean_distance(
    g,
    weights = NA,
    directed = FALSE,
    unconnected = FALSE
  )

  
  return(c(Q = Q, ASPL = L))
}


# baseline before removal
baseline <- network_metrics(brain_uu)

# independently remove each hub
for (i in seq_len(nrow(top_hubs))) {
  
  g <- delete_vertices(brain_uu, top_hubs$short_label[i])
  
  m <- network_metrics(g)
  
  top_hubs$delta_Q[i] <- m["Q"] - baseline["Q"]
  top_hubs$delta_ASPL[i] <- m["ASPL"] - baseline["ASPL"]
}

# dividing groups for comparison and combine measures
top_hubs$PC_group <- "Lower PC"

top_hubs$PC_group[
  top_hubs$participation > median(top_hubs$participation)
] <- "Higher PC"

top_hubs$abs_delta_Q <- abs(top_hubs$delta_Q)

aggregate(
  cbind(participation, degree, delta_ASPL, abs_delta_Q) ~ PC_group,
  data = top_hubs,
  FUN = mean
)

top_hubs$short_label

# rank by disruptions to ASPL and Q 
ranked_hubs <- top_hubs[
  order(-top_hubs$delta_ASPL),
]

ranked_hubs[, c(
  "short_label",
  "long_label",
  "degree",
  "participation",
  "delta_ASPL",
  "delta_Q"
)]

ranked_Q <- top_hubs[
  order(-top_hubs$delta_Q),
]

ranked_Q[, c(
  "short_label",
  "long_label",
  "degree",
  "participation",
  "delta_Q",
  "delta_ASPL"
)]


# export results

write.csv(top_hubs, "hub_lesion_results.csv", row.names = FALSE)

# inspect community membership of selected hubs
braincomm[c("LLOCsd2", "LPG10")]
V(brain_uu)$name[braincomm == braincomm["LPG10"]]
