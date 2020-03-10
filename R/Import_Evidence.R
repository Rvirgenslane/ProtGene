#' Prepares metadata file object for Proteus
#' @export
#' @inheritParams proteus::makePeptideTable
#' @return data.frame with sample metadata
#' @examples
#' devtools::install_github("bartongroup/proteusLabelFree")
#' library(proteusLabelFree)
#' evidenceFile <- system.file("extdata", "evidence.txt.gz", package="proteusLabelFree")
#' evi <- readEvidenceFile(evidenceFile)
#' Meta<-Proteus_Prep(evi=evi, measure.cols="intensity")
#' Meta

# devtools::create("ProGene")
Proteus_Prep<-function(evi=NULL, measure.cols="Intensity"){
# requires "condition", "sample", "experiment"
MS_metadata<-data.frame(experiment=unique(evi$experiment),
stringsAsFactors = FALSE)

# building required variables
MS_metadata$measure <- measure.cols
MS_metadata$sample <- MS_metadata$experiment

MS_metadata<-MS_metadata[order(MS_metadata$experiment),]
rownames(MS_metadata)<-NULL
return(MS_metadata)

}
