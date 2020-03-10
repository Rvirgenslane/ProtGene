#' Prepares adds gene symbols to evidence file object
#' @export
#' @inheritParams proteus::makePeptideTable
#' @param Species indicating the species. "Human" (default),
#' "Mouse", or "Yeast" are the only options.
#' @return data.frame with uniprot ids replaced with gene "SYMBOL".
#' For yeast, "GENENAME" is returned.
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Sc.sgd.db
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



Evidence_Symbols<-function(evi=NULL, Species="Human"){

# Converting to human symbols
if(Species == "Human"){
evi$SYMBOL<-AnnotationDbi::mapIds(org.Hs.eg.db,
keys=evi$protein, column='SYMBOL',
keytype='UNIPROT')
}else if(Species == "Mouse"){
evi$SYMBOL<-AnnotationDbi::mapIds(org.Mm.eg.db,
keys=evi$protein, column='SYMBOL',
keytype='UNIPROT')
}else if(Species == "Yeast"){
evi$SYMBOL<-AnnotationDbi::mapIds(org.Sc.sgd.db,
keys=evi$protein, column='GENENAME',
keytype='UNIPROT')
}

evi$protein<-evi$SYMBOL
evi$SYMBOL<-NULL
return(evi)

}
