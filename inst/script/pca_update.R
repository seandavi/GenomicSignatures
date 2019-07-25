rm( list=ls() )
library(siyuan)
library( ggplot2 )
library( grid )
library(tidyr)

results.dir <- dir.create2('results/PCA_update/')


# loading data and cleaning -----------------------------------------------

load('data/clean_new/setNames.RData')
load('data/clean_new/trainingSetNames.RData')
for (set in setNames) {
  load(set %>% paste0('data/clean_new/CRCSC_filled/', ., '.RData'))
}

expr.all.mapped <- lapply(trainingSetNames, function(set) set %>% get %>% exprs)
names(expr.all.mapped) <- trainingSetNames
genes.common <- lapply(expr.all.mapped, row.names) %>% Reduce('intersect', .)
expr.all.mapped <- lapply(expr.all.mapped, '[', genes.common, )

##removing MEOX2
expr.all.mapped<-lapply(expr.all.mapped,function(expr.dat)
{
  expr.dat<-(expr.dat[-(which(rownames(expr.dat)=="MEOX2")),])
  return(expr.dat)
})

data.sizes <- unlist( lapply( expr.all.mapped, function(expr.dat){ncol(expr.dat)} ) )
data.sizes <- data.frame( size=data.sizes, name=names(data.sizes), alias=paste0("ds", 1:8) )
rownames( data.sizes ) <- paste0("ds", 1:8)


# PCA ---------------------------------------------------------------------

#calculate PC for the training datasets
pca.expr.all.mapped<-lapply(expr.all.mapped,function(pc)
{
  dat.pca<-prcomp(t(pc))
  return(dat.pca)
})

# check PC variance proportions
sapply( pca.expr.all.mapped, function( dat.pca ) {
  return( c( sum( (dat.pca$sdev^2)[1:8] ), sum( (dat.pca$sdev^2)[1:20] ) ) /
            sum( dat.pca$sdev^2 )  )
} )

#First 20 PC loadings for each training dataset

data.loadings<-lapply(pca.expr.all.mapped,function(x)
{
  loadings<-x$rotation[,1:20]
  return(loadings)
})

loadings.dir <- dir.create2("results/PCA_update/loadings/")

# csv files for loadings
for (i in 1:length(data.loadings))
{
  write.csv(data.loadings[[i]],file=paste(loadings.dir,names(data.loadings)[i],".csv",sep=""))
}

# compare loadings
mat.cor.compare <- sapply(trainingSetNames, function(set) {
  df.loading1 <- set %>% paste0('results/PCA/loadings/', ., '.csv') %>% 
    read.csv(row.names = 1)
  df.loading2 <- set %>% paste0('results/PCA_update/loadings/', ., '.csv') %>% 
    read.csv(row.names = 1)
  genes.tmp <- intersect(rownames(df.loading1), rownames(df.loading2))
  cor.tmp <- sapply(1:20, function(i) {
    cor(df.loading1[genes.tmp, i], df.loading2[genes.tmp, i]) %>% return
  })
  return(cor.tmp)
})
dir.create2('results/PCA_update/checking/') %>% 
  paste0('correlation_plot.pdf') %>% 
  pdf(width = 10, height = 10)
mat.cor.compare %>% data.frame(n_loading = 1:20) %>% 
  gather(key = dataset, value = correlation, GSE13294_eset:NHS.HPFS_eset) %>% 
  ggplot(aes(x = n_loading, y = abs(correlation))) + geom_point() +
  geom_line() + facet_wrap(~dataset, nrow = 3) + theme_bw()
dev.off()

#PC scores
data.pca.scores<-lapply(pca.expr.all.mapped,function(x)
{
  score<-x$x[,1:8]
  return(score)
})
scores.dir <- dir.create2("results/PCA_update/scores/")
for (i in 1:length(data.pca.scores))
{
  write.csv(data.pca.scores[i],file=paste(scores.dir,names(data.pca.scores)[i],".csv",sep=""))
}

for(i in 1:length(data.loadings)){
  pca.data<-data.loadings[[i]]
  colnames(pca.data)<-paste("ds",i,".pc",1:20,sep="")
  data.loadings[[i]]<-pca.data
}
save.image(paste(results.dir,"pca_loadings_scores.RData",sep="")) #for average gsea

cor.list<-vector("list",length(data.loadings))#list of matrix with correlated PCs for each dataset
names(cor.list)<-names(expr.all.mapped)
cor.dir<-"results/PCA_update/correlations/"
dir.create(cor.dir)

for (i in 1:length(data.loadings)){
  pc.cor<-matrix(nrow=20)
  rownames(pc.cor)<-colnames(data.loadings[[i]])
  for(j in 1:length(data.loadings)){
    #labels[j]<-paste(names(pc.train[i]),"-",names(expr.all.test[j]), sep="")
    if(j==i)next
    #     expr.test<-expr.all.test[[j]]
    #     pc.test<-prcomp(t(expr.test))
    #     pc.test.loadings<-pc.test$rotation[,1:8]
    #     write.csv(pc.test.loadings,file=paste(output.dir[i],names(expr.all.test)[j],".csv",sep=""))
    #     colnames(pc.test.loadings)<-paste("test.","ds",j,".pc",1:8,sep="")
    cor.matrix<-cor(data.loadings[[i]],data.loadings[[j]], method="pearson")
    #     reference.loadings[[j]]<-pc.test.loadings
    pc.cor<-cbind(pc.cor,cor.matrix)#this will be saved as 8 20x20 in a list
  }
  #names(pc.cor)<-labels
  pc.cor<-pc.cor[,-1]
  cor.list[[i]]<-pc.cor
}

#CSV file for correlation matrices
for(i in 1:length(cor.list))
{
  fname<-paste(cor.dir,names(cor.list)[i],".csv",sep="")
  write.csv(cor.list[[i]],file=fname)
}

#Correlation Co-eff. >0.5 
#gives the number of correlation co-eff. >0.5
# pass<-numeric()
# 
# for(i in 1:length(cor.list))
# {
#   pass[i]<-nrow(which(cor.list[[i]]>0.5,arr.in=T))
# }
# names(pass)<-names(cor.list)
save.image(paste(cor.dir,"cor.list.RData",sep=""))

#Extract the highly correlated DSPC names
high.cor.dir<-"results/PCA_update/highly_correlated/"
dir.create(high.cor.dir)
fil.cor.list<-lapply(cor.list,function(x)
{
  x<-abs(x)
  high.cor<-which(x>0.5,arr.in=T)
  high.cor.list<-matrix(nrow=nrow(high.cor),ncol=3)
  #extract the rows
  for(i in 1:nrow(high.cor))
  {
    high.cor.list[i,1]<-rownames(x)[high.cor[i,1]]
  }
  #extract columns
  for(j in 1:nrow(high.cor))
  {
    high.cor.list[j,2]<-colnames(x)[high.cor[j,2]]
  }
  
  #value
  for(k in 1:nrow(high.cor))
  {
    high.cor.list[k,3]<-x[high.cor[k,1],high.cor[k,2]]
  }
  
  return(high.cor.list)
}
)

for(p in 1:length(fil.cor.list))
{
  fname<-paste(high.cor.dir,names(fil.cor.list)[p],".csv",sep="")
  write.csv(fil.cor.list[p],file=fname)
}
net.dir<-"results/PCA_update/network/"
dir.create(net.dir)
network.mat<-matrix(ncol=3)
for(i in 1:length(fil.cor.list))
{
  network.mat<-rbind(network.mat,fil.cor.list[[i]])
}
network.mat<-network.mat[-1,]
write.csv( network.mat,file=paste(net.dir,"network.matrix.csv",sep=""), quote=F )
nodes <- unique( network.mat[, 1] )
nodes.set <- unlist( lapply( strsplit( nodes, "\\." ), function(x){x[1]} ) )
nodes.attribute <- data.frame( node=nodes, size=data.sizes[nodes.set, 1], setname=data.sizes[nodes.set, 2], alias=data.sizes[nodes.set, 3])
write.csv( nodes.attribute, file=paste(net.dir, "nodes.csv", sep=""), quote=F )

#####################################################################
##Validation & Assignment of PC scores
#####################################################################
library( Biobase )
data.dir <- "data/clean_new/"

valiSetNames <- setdiff( setNames, trainingSetNames )

dat1 <- data.loadings[["GSE13294_eset"]]
dat2 <- data.loadings[["GSE14095_eset"]]
dat3 <- data.loadings[["GSE14333_eset"]]
dat4 <- data.loadings[["GSE17536_eset"]]
dat5 <- data.loadings[["GSE21815_eset"]]
dat6 <- data.loadings[["GSE26682.GPL570_eset"]]
dat7 <- data.loadings[["GSE26682.GPL96_eset"]]
dat8 <- data.loadings[["NHS.HPFS_eset"]]

merge.dat <- list()
# cor.list[[1]][1,][c('ds2.pc2', 'ds3.pc1', 'ds4.pc1', 'ds5.pc2', 'ds5.pc3', 'ds6.pc1', 'ds6.pc2', 'ds7.pc1', 'ds8.pc7')]
merge.dat[[1]] <- data.frame(dat1[, 1], dat2[, 2], -dat3[, 1], dat4[, 1],
                            dat5[, 2], -dat5[, 3], -dat6[, 1], dat6[, 2], 
                            dat7[, 1], dat8[, 7],row.names=rownames(dat1))
# cor.list[[4]][2,][c('ds1.pc2', 'ds2.pc3', 'ds3.pc2', 'ds3.pc3', 'ds5.pc4', 'ds7.pc3', 'ds8.pc8')]
merge.dat[[2]] <- data.frame(dat1[, 2], -dat2[, 3], dat3[, 2], dat3[, 3], 
                             dat4[, 2], -dat5[, 4], dat7[, 3], -dat8[, 8],
                             row.names=rownames(dat1))
# cor.list[[2]][5,][c('ds1.pc3', 'ds3.pc4', 'ds4.pc4', 'ds6.pc6', 'ds7.pc4', 'ds7.pc5')]
merge.dat[[3]] <- data.frame(dat1[, 3], dat2[, 5], dat3[, 4], dat4[, 4],
                             dat6[, 6], -dat7[, 4], -dat7[, 5], 
                             row.names=rownames(dat1))
# cor.list[[2]][4,][c('ds1.pc4', 'ds4.pc3', 'ds5.pc1', 'ds6.pc3', 'ds7.pc2')]
merge.dat[[4]] <- data.frame( dat1[, 4], dat2[, 4], dat4[, 3], dat5[ ,1], 
                              dat6[, 3], -dat7[, 2], row.names=rownames(dat1))

avg.loadings <- data.frame( lapply( merge.dat, function( dat ){ apply( dat, 1, mean, na.rm=F ) } ), row.names=rownames( merge.dat[[1]] ) )
colnames( avg.loadings ) <- paste0( 'cl', 1:4 )

write.csv( avg.loadings, paste0( loadings.dir, 'avg_loadings.csv' ) )

# check correlation between average loadings
avg.loadings1 <- read.csv('results/PCA/loadings/avg_loadings.csv', row.names = 1)
avg.loadings2 <- read.csv('results/PCA_update/loadings/avg_loadings.csv', row.names = 1)
genes.tmp <- intersect(rownames(avg.loadings1), rownames(avg.loadings2))
sapply(1:4, 
       function(i) cor(avg.loadings1[genes.tmp, i], avg.loadings2[genes.tmp, i])) %>% 
  data.frame(loading = 1:4, .) %>% 
  write.csv('results/PCA_update/checking/correlation_avg_loadings.csv')


# Assignment of score & validation at the same time
dir.create('results/PCA_update/new_scores/', showWarnings = F, recursive = T)
valiSetNames <- setdiff(setNames, trainingSetNames)
corrs <- list()
for ( i in 1:4 ) {
  corrs[[i]] <- matrix( NA, 8, length( valiSetNames ) )
  colnames( corrs[[i]] ) <- valiSetNames
}

for (set in setNames ) {
  
  print( set )
  eSet <- get(set)
  exprs <- exprs(eSet)[-(which(rownames(eSet)=="MEOX2")),]
  exprs <- exprs[apply(exprs, 1, function(x){!any(is.na(x)|(x==Inf)|(x==-Inf))}) ,]
  exprs <- apply(exprs, 1, function(x) x - mean(x)) %>% t
  genesCommon <- intersect( rownames(exprs), rownames(avg.loadings) )
  lcss <- t( exprs[genesCommon, ] ) %*% apply( avg.loadings[genesCommon, ], 2, function( x ) x / sqrt( sum( x^2, na.rm=T ) ) )
  colnames( lcss ) <- paste0( 'PCSS', 1:4 ) 
  
  if ( set %in% valiSetNames ) {
    #calculate PC loadings for the validation datasets
    prcomRes <- prcomp( t(exprs) )
    loadings <- prcomRes$rotation[,1:8]
    
    for ( i in 1:4 ) {
      for ( j in 1:8 ){
        corrs[[i]][j, set] <- cor( x=loadings[genesCommon, j], y=avg.loadings[genesCommon, i], use='pairwise.complete.obs' )
      }		
    }
  }
  write.csv(lcss, file = paste0('results/PCA_update/new_scores/', set, '.csv'))
}

# checking to see if the scores are agreeing
mat.cor.score <- setNames %>% sapply(function(set) {
  df.score1 <- get(set) %>% pData %>% '['( , paste0('LCSS', 1:4))
  df.score2 <- read.csv(paste0('results/PCA_update/new_scores/', set, '.csv'),
                        row.names = 1)
  samples.common <- intersect(rownames(df.score1), rownames(df.score2))
  sapply(1:4, function(i) {
    cor(df.score1[samples.common, i], df.score2[samples.common, i])
  })
}) %>% t
colnames(mat.cor.score) <- paste0('LCSS', 1:4)
pdf('results/PCA_update/checking/correlation_scores_new_centered.pdf', 
    width = 10, height = 8)
mat.cor.score %>% data.frame(dataset = rownames(mat.cor.score)) %>% 
  gather(key = 'LCSS', value = 'correlation', LCSS1:LCSS4) %>% 
  ggplot(aes(x = dataset, y = correlation, color = LCSS)) + 
  geom_point(alpha = 0.5) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

results <- matrix( NA, 4, 8 )
dimnames( results ) <- list( paste0("Loading for PCSS", 1:4), paste0("First ", 1:8, " PC Loadings") )
for ( i in 1:4 ) {
  for ( j in 1:8  ) results[i, j] <- sum( apply( matrix(t(abs(corrs[[i]][1:j, ])), nrow=j, byrow=T) > 0.5, 2, sum, na.rm=T  ) > 0, na.rm=T  )
}

vali.dir <- 'results/PCA_update/validation/'
dir.create( vali.dir, recursive=T )
write.csv( results, file=paste0(results.dir, "validation/validation.csv") )


dir.create('results/PCA_update/new_scores_centered/', showWarnings = F, recursive = T)

for (set in setNames ) {
  
  print( set )
  eSet <- get(set)
  exprs <- exprs(eSet)[-(which(rownames(eSet)=="MEOX2")),]
  exprs <- exprs[apply(exprs, 1, function(x){!any(is.na(x)|(x==Inf)|(x==-Inf))}) ,]
  exprs <- apply(exprs, 1, function(x) x - mean(x)) %>% t
  genesCommon <- intersect( rownames(exprs), rownames(avg.loadings1) )
  lcss <- t( exprs[genesCommon, ] ) %*% apply( avg.loadings1[genesCommon, ], 2, function( x ) x / sqrt( sum( x^2, na.rm=T ) ) )
  colnames( lcss ) <- paste0( 'PCSS', 1:4 ) 
  
  write.csv(lcss, file = paste0('results/PCA_update/new_scores_centered/', set, '.csv'))
}

# checking to see if the scores are agreeing
mat.cor.score <- setNames %>% sapply(function(set) {
  df.score1 <- get(set) %>% pData %>% '['( , paste0('LCSS', 1:4))
  df.score2 <- read.csv(paste0('results/PCA_update/new_scores_centered/', set, '.csv'),
                        row.names = 1)
  samples.common <- intersect(rownames(df.score1), rownames(df.score2))
  sapply(1:4, function(i) {
    cor(df.score1[samples.common, i], df.score2[samples.common, i])
  })
}) %>% t
colnames(mat.cor.score) <- paste0('LCSS', 1:4)
pdf('results/PCA_update/checking/correlation_scores_centered.pdf', 
    width = 10, height = 8)
mat.cor.score %>% data.frame(dataset = rownames(mat.cor.score)) %>% 
  gather(key = 'LCSS', value = 'correlation', LCSS1:LCSS4) %>% 
  ggplot(aes(x = dataset, y = correlation, color = LCSS)) + 
  geom_point(alpha = 0.5) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
