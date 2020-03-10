#' Prepares metadata file object for Proteus
#' @export
#' @inheritParams proteus::makePeptideTable
#' @return data.frame with sample metadata
#' @examples
#' library(proteusLabelFree)
#' evi<-proteusLabelFree::evi
#' evi$protein<-gsub(".*[sp|]([^.]+)[|].*", "\\1", evi$protein)
#' evi_symbols<-ProtGene::Evidence_Symbols(evi=evi, Species="Yeast")
#' head(evi_symbols)
#'
#' Meta<-Proteus_Prep(evi=evi_symbols, measure.cols="intensity")
#' Meta$condition<-gsub("-[1-9]","", Meta$experiment)
#' Meta$measure<-"Intensity"

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
