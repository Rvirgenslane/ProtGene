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
#' library(proteusLabelFree)
#' evi<-proteusLabelFree::evi
#'
#' evi$protein<-gsub(".*[sp|]([^.]+)[|].*", "\\1", evi$protein)
#' evi_symbols<-ProtGene::Evidence_Symbols(evi=evi, Species="Yeast")
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
