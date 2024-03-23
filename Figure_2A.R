###This code will make a plot similar to Figure 2A from the manuscript. For the actual figure, individual clusters were plotted and arranged. This code plots all clusters in a single heatmap in the same order as the manuscript.

###libraries needed
library(viridis)
library(superheat)

#####################################
genemeans <- plotclusters2 %>% dplyr::group_by(gene, time) %>% dplyr::summarise(mean_value = mean(value))

# Extract unique gene-cluster mappings
gene_cluster_unique <- plotclusters2 %>%
  dplyr::select(gene, CLUSTER) %>%
  dplyr::distinct()

# Define your desired cluster order
desired_cluster_order <- c(4,3,1,2,5)

genemeans_ordered <- genemeans %>%
  dplyr::inner_join(gene_cluster_unique, by = "gene") %>%
  dplyr::arrange(match(CLUSTER, desired_cluster_order), gene)

genemeans_wide <- genemeans_ordered %>%
  tidyr::pivot_wider(names_from = time, values_from = mean_value)

# Create a vector of colors, one for each cluster
cluster_colors <- RColorBrewer::brewer.pal(max(genemeans_wide$CLUSTER), "Set1")

# Create a named vector of colors for the clusters
cluster_color_map <- setNames(cluster_colors, unique(genemeans_wide$CLUSTER))

# Assign colors to each row based on their cluster
row_side_colors <- cluster_color_map[as.character(genemeans_wide$CLUSTER)]

#colors for heatmap
color_ramp <- viridis::cividis(20)


myclusters<-as.vector(genemeans_wide[,2])

superheat(genemeans_wide[c(-1,-2)], # Exclude the gene column for the heatmap
          heat.pal = color_ramp,
		   bottom.label.size = 0.05,
		   bottom.label.text.size = 2,
		   legend.height = 0.05,
          legend.width = 0.5,
          legend.text.size = 5,
		  membership.rows = myclusters$CLUSTER,
		   left.label.size = 0.06,
		  left.label.text.size = 4,
		  #grid.hline.size = 0.01,
		  grid.hline = FALSE)   # Use cviridis color palette


