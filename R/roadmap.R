## ---- echo=FALSE, eval=TRUE---------------------------------------------------
library(visNetwork)
nodes <- data.frame(id = c("W", "A", "Y"))
nodes$label <- nodes$id
edges <- data.frame(from = c("W", "W", "A"), to = c("A", "Y", "Y"))
network <- visNetwork(nodes, edges, height = "300px", width = "200px") %>%
  visEdges(arrows = list(to = TRUE)) %>%
  visLayout(randomSeed = 25)
network

