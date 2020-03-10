#' Prepares metadata file object for Proteus
#' @export
#' @inheritParams proteus::makePeptideTable
#' @return data.frame with sample metadata
#' @examples
#' # Getting evidence file
#' library(ProtGene)
#' library(proteusLabelFree)
#' evidenceFile <- system.file("extdata", "evidence.txt.gz", package="proteusLabelFree")
#' 
#' # The directory of the evidence file is required
#' Evidence_input <- proteus::readEvidenceFile(evidenceFile)
#' 
#' # extracting uniprot ids: "sp|P00175|CYB2_YEAST" to "P00175
#' # Your data may not need this step
#' Evidence_input$protein<-gsub(".*[sp|]([^.]+)[|].*", "\\1", Evidence_input$protein)
#' head(Evidence_input)
#' 
#' # Converting to gene name
#' evi_symbols<-ProtGene::Evidence_Symbols(evi=Evidence_input, Species="Yeast")
#' head(evi_symbols)
#' 
#' # Generating metadata from evidence file
#' Meta<-Proteus_Prep(evi=evi_symbols)
#' head(Meta)

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
