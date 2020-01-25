These are gene sets from MSigDB (http://software.broadinstitute.org/gsea/msigdb/collections.jsp).   

"c2.cgp.v7.0.symbols.gmt" : chemical and genetic perturbations (3,302 gene sets)   
"c2.cp.v7.0.symbols.gmt"  : Canonical pathways (2,199 gene sets)   
"c6.all.v7.0.symbols.gmt" : oncogenic signatures (189 gene sets)   
"c7.all.v7.0.symbols.gmt" : immunologic signatures (4,872 gene sets)   

Prior knowledge included in the PLIER package is out of data, so the number of genes
and pathways are different from the version I download (v.7.0). Here is the brief summary
on the number of genes/pathways of MSigDB datasets.

For this package:   
> dim(c2_cgp)   
[1] 20181  3302   
> dim(c2_cp)   
[1] 11763  2199   
> dim(c6)   
[1] 10962   189   
> dim(c7)   
[1] 20437  4872   

PLIER package:   
> data(chemgenPathways)   
> dim(chemgenPathways)   
[1] 20554  3395   
> data(canonicalPathways)   
> dim(canonicalPathways)   
[1] 6023  545   
> data(oncogenicPathways)   
> dim(oncogenicPathways)   
[1] 11250   189   
> data(immunePathways)   
> dim(immunePathways)   
[1] 19841  1910   
