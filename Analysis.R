#READING IN DATA 

setwd("C:/users/sgura/Desktop/DS_Projects/kaggle/METABRIC")
metabric <- read.csv("data/METABRIC_RNA_mutation.csv")

#FILTERING DATA 

library(dplyr)
names_metabric <- names(metabric)
#filters names_metabric for RNA expression columns
names_metabric <- names_metabric[seq(32,520)]
#filters metabric to exclude other cause deaths and selects expression columns
metabric_df <- na.omit(metabric %>% 
                        filter(death_from_cancer != "Died of Other Causes") %>% 
                        select(patient_id, 
                               death_from_cancer, 
                               all_of(names_metabric)))

#Reformats death_from_cancer category
metabric_df$death_from_cancer <- gsub("Living", 0, 
                                      metabric_df$death_from_cancer)
metabric_df$death_from_cancer <- gsub("Died of Disease", 1, 
                                      metabric_df$death_from_cancer)
metabric_df$death_from_cancer <- as.numeric(metabric_df$death_from_cancer)

#PCA

#drop nas and invert so genes are rows
metabric_mx <- as.matrix(na.omit(metabric_df[, 1:491]))
colnames_mt_mx <- metabric_mx[, 1]
metabric_mx <- t(metabric_mx[, 2:491])
#creating pca object
metabric_pca <- prcomp(metabric_mx[2:490, ], scale = TRUE)

#exploring PCA visually
library(ggplot2)
pc_scores <- metabric_pca$x
pc_scores_plot <- as.data.frame(pc_scores) %>%
                    ggplot(aes(x = PC1, y = PC2)) +
                    geom_point() +
                    labs(title = "PC Scores for RNA seq z scores data") +
                    theme_bw()
                  
#creating biplot
library(ggbiplot)
metabric_biplot <- ggbiplot(metabric_pca) +
                    theme_bw() +
                    scale_color_manual(labels = c("Died from Disease",
                                                  "Living"), 
                                       values = c("darkorange", "purple")) 
#drops text labels crowding graph                    
metabric_biplot$layers[[3]] <- NULL

#creating screeplot
screeplot(metabric_pca, 10, "barplot")

# K MEANS CLUSTERING 

# SSE estimations

#for loop calculating sse for k 1 to 20 clusters
metabric_sse <- vector()
for (i in 2:20) {
  metabric_sse[i] <- sum(kmeans(metabric_mx[2:490, ], centers = i)$withinss)
}

#creates df of metabric sse values
metabric_sse_df <- na.omit(data.frame(sse = metabric_sse,
                                       clusters = 1:20))
#Creates SSE plot 
sse_plot <- ggplot(metabric_sse_df, aes(x = clusters, y = sse)) +
              geom_point() +
              geom_smooth() +
              theme_bw() +
              labs(title = "SSEs for 2 to 20 clusters using kmeans",
                    y = "Sum of Squared Error") 

## Average Silhouette Width
library(cluster)

#for loop calculating silhouette width for 1 to 20 clusters 
metabric_sil <- vector()
for (i in 2:20) {
  k_one_to_twenty <- kmeans(metabric_mx[2:490, ], centers = i, nstart = 25, 
                            iter.max = 20)
  ss <- silhouette(k_one_to_twenty$cluster, dist(metabric_mx[2:490, ]))
  metabric_sil[i] <- mean(ss[, 3])
}
#creates silhoutte width df
metabric_sil_df <- na.omit(data.frame(AvgSilhouetteWidth = metabric_sil,
                                      cluster = 1:20))
#plots silhouette width
sil_plot <- ggplot(metabric_sil_df, aes(x = cluster, y = AvgSilhouetteWidth)) +
              geom_point() +
              geom_smooth() +
              labs(title = "Average Silhouette Width for 2 to 20 clusters",
                   y = "Average Silhouette Width") +
              theme_bw()

## Calinksy Criterion 
library(vegan)
#calculates calinsky index for 1 to 20 clusters 
metabric_calinksi <- cascadeKM(metabric_mx[2:490, ], 1, 20, iter = 100)
#Creates heatmap and index plot
plot(metabric_calinksi, sortg = TRUE, grpmts.plot = TRUE)


## Gap Statistic
set.seed(13)
#calculates Gap statistics 
metabric_gap <- clusGap(metabric_mx[2:490, ], kmeans, 20, B = 100, 
                        verbose = interactive())
#plot for gap statistics
plot(metabric_gap, main = "Gap statistic across 20 clusters")

# Clustering
set.seed(20)
kClust <- kmeans(metabric_mx[2:490, ], centers = 4,
                 nstart = 1000, iter.max = 20)
metabric_clust <- kClust$cluster

## Finding centroids; takes col means
cluster_centroid <- function(num, data, clusters) {
  x <- (clusters == num)
  colMeans(data[x, ])
}

#calculates centroids 
metabric_clust_centroids <- sapply(levels(factor(metabric_clust)), 
                                   cluster_centroid,
                                   metabric_mx[2:490, ],
                                   metabric_clust)
#assigns rownames to centroids 
rownames_centroids <- rownames(metabric_mx[2:490, ])
rownames(metabric_clust_centroids) <- rownames_centroids
#calculates correlation between centroids 
cor(metabric_clust_centroids)

## Visualizing cluster expression for 6 patients

#adding gene name as column and cluster# to metabric_mx_2 
metabric_mx_2 <- as.data.frame(metabric_mx)
rownames(metabric_mx_3) <- seq(1:489)
metabric_df_2 <- metabric_mx_3 %>% 
                  mutate(gene = names_mx_2) %>% 
                  mutate(cluster = as.vector(metabric_clust))
#Extracting expression data for 6 patients 
patient_col_names <- c("Expression", "gene", "cluster", "Patient")
patient1_df <- metabric_df_2 %>% 
  select(1418, gene, cluster) %>% 
  mutate(patient = "Patient 1") 
colnames(patient1_df) <- patient_col_names
patient2_df <- metabric_df_2 %>% 
  select(1423, gene, cluster) %>% 
  mutate(patient = "Patient 2")
colnames(patient2_df) <- patient_col_names
patient3_df <- metabric_df_2 %>% 
  select(1422, gene, cluster) %>% 
  mutate(patient = "Patient 3")
colnames(patient3_df) <- patient_col_names
patient4_df <- metabric_df_2 %>% 
  select(1421, gene, cluster) %>% 
  mutate(patient = "Patient 4")
colnames(patient4_df) <- patient_col_names
patient5_df <- metabric_df_2 %>% 
  select(1420, gene, cluster) %>% 
  mutate(patient = "Patient 5")
colnames(patient5_df) <- patient_col_names
patient6_df <- metabric_df_2 %>% 
  select(1419, gene, cluster) %>% 
  mutate(patient = "Patient 6")
colnames(patient6_df) <- patient_col_names
#combining patient dfs 
patientsdf <- rbind(patient1_df, patient2_df, patient3_df,
                     patient4_df, patient5_df, patient6_df)
#making plot
kmeans_patients <- ggplot(data = patientsdf, aes(x = cluster, y = Expression)) +
                    geom_bar(stat = "identity", fill = "#FF6666") +
                    facet_wrap(~Patient) +
                    labs(title = "Expression by kmeans cluster for 6 Patients",
                         x = "Cluster",
                         y = "Normalized Expression") +
                    theme_bw()

## Isolating core genes from cluster 1

#melt data to long form
metabric_clust_melt <- melt(metabric_clust_centroids)
#assign colnames 
colnames(metabric_clust_melt) <- c("sample", "cluster", "value")
#subset melted data frame
metabric_core_1 <- metabric_clust_melt %>% 
                    filter(cluster == "1")
#assign metabric_mx with only genes
metabric_mx_2 <- metabric_mx[2:490, ]

#pull out cluster 1
metabric_clust_1 <- metabric_mx_2[metabric_clust == 1, ]
#function that calculates correlation score for cluster 1
correlation_score_1 <- function(n) {
  cor(n, metabric_core_1$value)
}
#calculates cluster1 scores 
clust1_scores <- apply(metabric_clust_1, 1, correlation_score_1)
#identifies top 5 highest correlated genes 
core_clust_1_genes <- tail(sort(clust1_scores))

## Isolating core genes from rest of clusters
metabric_core_2 <- metabric_clust_melt %>% 
                    filter(cluster == "2")
metabric_clust_2 <- metabric_mx_2[metabric_clust == 2, ]
correlation_score_2 <- function(n) {cor(n, metabric_core_2$value)}
clust_2_scores <- apply(metabric_clust_2, 1, correlation_score_2)
core_clust_2_genes <- tail(sort(clust_2_scores))

metabric_core_3 <- metabric_clust_melt %>% 
                    filter(cluster == "3")
metabric_clust_3 <- metabric_mx_2[metabric_clust == 3, ]
correlation_score_3 <- function(n) {cor(n, metabric_core_3$value)}
clust_3_scores <- apply(metabric_clust_3, 1, correlation_score_3)
core_clust_3_genes <- tail(sort(clust_3_scores))

metabric_core_4 <- metabric_clust_melt %>% 
                    filter(cluster == "4")
metabric_clust_4 <- metabric_mx_2[metabric_clust == 4, ]
correlation_score_4 <- function(n) {cor(n, metabric_core_4$value)}    
clust_4_scores <- apply(metabric_clust_4, 1, correlation_score_4)
core_clust_4_genes <- tail(sort(clust_4_scores))

core_genes_df <- data.frame("cluster 1 genes" = names(core_clust_1_genes),
                            "cluster 1 scores" = core_clust_1_genes,
                            "cluster 2 genes" = names(core_clust_2_genes),
                            "cluster 2 scores" = core_clust_2_genes,
                            "cluster 3 genes" = names(core_clust_3_genes),
                            "cluster 3 scores"= core_clust_3_genes,
                            "cluster 4 genes" = names(core_clust_4_genes),
                            "cluster 4 scores" = core_clust_4_genes)
rownames(core_genes_df) <- c(1:6)

#HIERARCHICAL CLUSTERING 

#calculating gene distance for each row
metabric_gene_dist <- dist(metabric_mx_2)
#using hclust function on gene distance
metabric_hclust <- hclust(metabric_gene_dist, method = "complete")
#dendrogram of genes clustered
plot(metabric_hclust, labels = FALSE, main = "Gene Clustering Dendogram")
abline(h = 58, lwd = 2)

#seperating data by cluster and turning it into a tibble
metabric_hc_clusters <- as.vector(cutree(metabric_hclust, k = 6))
#adding gene name as column and cluster# to metabric_mx_2 
metabric_df_3 <- metabric_mx_3 %>% 
                  mutate(gene = names_mx_2) %>% 
                  mutate(cluster = metabric_hc_clusters) 
#Extracting expression data for 6 patients 
patient_1_df <- metabric_df_3 %>% 
                  select(1418, gene, cluster) %>% 
                  mutate(patient = "Patient 1") 
colnames(patient_1_df) <- patient_col_names
patient_2_df <- metabric_df_3 %>% 
                  select(1423, gene, cluster) %>% 
                  mutate(patient = "Patient 2")
colnames(patient_2_df) <- patient_col_names
patient_3_df <- metabric_df_3 %>% 
                  select(1422, gene, cluster) %>% 
                  mutate(patient = "Patient 3")
colnames(patient_3_df) <- patient_col_names
patient_4_df <- metabric_df_3 %>% 
                  select(1421, gene, cluster) %>% 
                  mutate(patient = "Patient 4")
colnames(patient_4_df) <- patient_col_names
patient_5_df <- metabric_df_3 %>% 
                  select(1420, gene, cluster) %>% 
                  mutate(patient = "Patient 5")
colnames(patient_5_df) <- patient_col_names
patient_6_df <- metabric_df_3 %>% 
                  select(1419, gene, cluster) %>% 
                  mutate(patient = "Patient 6")
colnames(patient_6_df) <- patient_col_names
#combining patient dfs 
patients_df <- rbind(patient_1_df, patient_2_df, patient_3_df,
                     patient_4_df, patient_5_df, patient_6_df)
#making plot
patient_plot <- ggplot(data = patients_df, aes(x = cluster, y = Expression)) +
                  geom_bar(stat = "identity", fill = "#FF6666") +
                  facet_wrap(~Patient) +
                  labs(title = "Expression by hierarchical cluster for 6 Patients",
                       x = "Cluster",
                       y = "Normalized Expression") +
                  theme_bw()
#Making heatmap
library(gplots)
hc_dist <- cor(t(metabric_mx_2), method = "pearson")
hc_hr <- hclust(as.dist(1 - hc_dist), method = "complete")

#heatmap.2 of hc_dist 
heatmap.2(hc_dist,
          Rowv = as.dendrogram(hc_hr),
          Colv = as.dendrogram(hc_hr),
          scale = "row",
          margins = c(2,2),
          cexCol = 0.7,
          labRow = FALSE,
          labCol = FALSE,
          main = "Hierarchically clustered Heatmap",
          trace = "none")
