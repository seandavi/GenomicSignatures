#' Names of ten OVC training datasets
#'
#' A list of the dataset names from curatedOvarianData package, used as a training
#' dataset to build a model. Eight datasets were processed by Affymetrix, and two
#' datasets were processed by Agilent (GSE17260_eset, GSE32062.GPL6480_eset)
#'
#' @format a list object with 10 character strings
#' @references Ganzfried et al. (2013) Database
#' (\url{https://www.ncbi.nlm.nih.gov/pubmed/?term=curatedOvarianData%3A+clinically+annotated+data+for+the+ovarian+cancer+transcriptome}{PubMed})
#' @source \url{https://github.com/waldronlab/curatedOvarianData/tree/master/data}{curatedOvarianData}
#' @return
"trainingSets_10ovc"

#' Names of eight CRC training datasets
#'
#' @format a list object with 8 character strings
#' @references Ma et al. (2018) Genome Biol.
#' (\url{https://www.ncbi.nlm.nih.gov/pubmed/?term=Continuity+of+transcriptomes+among+colorectal+cancer+subtypes+based+on+meta-analysis}{PubMed})
#' @source \url{https://bitbucket.org/biobakery/crc-subtyping-paper/src/master/data/eSets/}{curatedCRCData}
#' @return
"trainingSets_8crc"

#' Eight CRC expression datasets
#'
#' List of length 8, where expressions are centered and scaled on a per-gene level
#' first, per usual PCA standard. "MEOX2" gene is removed.
#'
#' @format a list object with 8 matrices.
#' @references Ma et al. (2018) Genome Biol.
#' (\url{https://www.ncbi.nlm.nih.gov/pubmed/?term=Continuity+of+transcriptomes+among+colorectal+cancer+subtypes+based+on+meta-analysis}{PubMed})
#' @source \url{https://bitbucket.org/biobakery/crc-subtyping-paper/src/master/data/eSets/}{curatedCRCData}
#' @return
"expr.all.mapped"







#' Colon cancer (CRC) datasets from curatedCRCData for training
#'
#' @format A list object with 8 character strings
#'
#' @references Ma et al. (2018) Genome Biol.
#' (\url{https://www.ncbi.nlm.nih.gov/pubmed/?term=Continuity+of+transcriptomes+among+colorectal+cancer+subtypes+based+on+meta-analysis}{PubMed})
#' @source \url{https://bitbucket.org/biobakery/crc-subtyping-paper/src/master/data/eSets/}{curatedCRCData}
#' @return
"PCs_10ovc_training"
