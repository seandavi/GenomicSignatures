# GenomicSignatures

<a href="https://github.com/shbrief/Diagrams/blob/master/GenomicSignatures.png"><img src="https://raw.githubusercontent.com/shbrief/Diagrams/master/GenomicSignatures.png?token=ADX67SWALLVD7DRT4IMU4N26HOT2Y"/></a> 

## Background
Dimension reduction such as PCA helps to interpret RNA sequencing and other genomic data by identifying a few components that explain most of the variability in the dataset at hand. Although more than 200k RNA sequencing profiles have been deposited in public archives, exploratory analysis of new studies is often performed without comparison to these archives. To help identify robust and replicable axes of variation in any given RNA-seq dataset, we performed dimension reduction on 146 recount2 studies with >= 50 samples in each (the total number of samples is 22,765) to identify recurrent loadings vectors. These recurrent vectors are interpreted through design of the originating studies as well as gene set analysis. The GenomicSignatures method matches the major axes of variation in new datasets to those of previously identified studies, helping to identify replicable axes of variation in a way that is robust to batch and study effects. 

