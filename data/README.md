# GSEA datasets
This is from the newer version (v.7.0) than PLIER's. Pre-processing is described
in `/data2/Genomic_Super_Signature/GSEA/priorKnowledge.Rmd`.

- `CanonicalPathways.rda`
- `ChemicalGeneticPerturvation.rda`
- `OncogenicSignatures.rda`
- `ImmunologicSignatures.rda`
- `MSigDB_c2c6c7.rds`: combine all 4 gene sets from MSigDB used in PLIER package.  

# `expr.all.mapped.rda`

# trainingSets_8ovc.rda
## 8 OVC training sets from curatedOvarianData (by Affymetrix)
GSE20565_eset   
GSE2109_eset   
GSE26193_eset   
GSE26712_eset   
GSE6008_eset   
GSE9891_eset   
PMID17290060_eset   
TCGA_eset   

# trainingSets_10ovc.rda
## 2 OVC training sets from curatedOvarianData (by Agilent)
GSE17260_eset   
GSE32062.GPL6480_eset   

# trainingSets_8crc.rda
## 8 CRC training sets used in Ma et al. Genome Biology (2018)
GSE13294_eset   
GSE14095_eset   
GSE14333_eset   
GSE17536_eset   
GSE21815_eset   
GSE26682.GPL570_eset   
GSE26682.GPL96_eset   
NHS.HPFS_eset   

# testSets_tcga.rda
## 4 RNAseq datasets from curatedTCGAData
(log2 transformation is required)    
TCGA_BRCA_rseq   
TCGA_COAD_rseq   
TCGA_LUAD_rseq   
TCGA_READ_rseq   

## 1 OVC RNAseq dataset from curatedOvarianData
(log2 transformation isn't requried (already done))     
TCGA.RNASeqV2_eset   

# TCGA_trainingSets.rda
5 TCGA RNAseq raw counts were log2-transformed and idMap to "SYMBOL"

# TCGA.RNASeqV2_eset.RData
195 TCGA colon cancer samples