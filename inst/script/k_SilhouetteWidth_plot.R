# Plot the mean of silhouette width for varying k
library(Biobase)
library(tidyverse)
library(cluster)
library(GenomicSuperSignature)

ov.pc = system.file("extdata", "PCs_10ovc_training.csv", package = "GenomicSuperSignature")
ov = read.csv(ov.pc, row.names = 1)
dat = ov

K_max = ncol(dat)/5 - 1
K_list = c((1:K_max)*5)
set.seed(123)

l_clst_kmeans = lapply(
    K_list,
    function(k) kmeans(x = t(dat), centers = k)$cluster
)

sw_kmeans = lapply(
    l_clst_kmeans,
    function(clst_tmp) silhouette(clst_tmp, daisy(t(dat)))
)

sw_kmeans_sil_width = sapply(sw_kmeans, function(sw) sw[,"sil_width"])

df_toreturnSWkmeans = data.frame(statistics = apply(sw_kmeans_sil_width, 2, mean),
                                 se = apply(sw_kmeans_sil_width, 2, sd),
                                 K = K_list,
                                 metric = "Wilhouette Width",
                                 method = "Kmeans")

# save the plot in pdf
pdf("k_sw_plot.pdf") 
plot(df_toreturnSWkmeans$statistics ~ df_toreturnSWkmeans$K)
dev.off() 
