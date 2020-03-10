#' Converts proteins table and metadata into eset_proteus
#' @export
#' @param Proteins_table Output of makeProteinTable function
#' @param MS_metadata Metadata object prepared by Proteus_Prep function
#' @return Expressset
#' @examples
#' library(proteusLabelFree)
#' data(proteusLabelFree)
#'
#' evi$protein<-gsub(".*[sp|]([^.]+)[|].*", "\\1", evi$protein)
#' evi_symbols<-ProtGene::Evidence_Symbols(evi=evi, Species="Yeast")
#' head(evi_symbols)
#'
#' Meta<-Proteus_Prep(evi=evi_symbols, measure.cols="intensity")
#' Meta$condition<-gsub("-[1-9]","", Meta$experiment)
#' Meta$measure<-"Intensity"
#'
#' pepdat <- proteus::makePeptideTable(evi_symbols, Meta, ncores=4)
#' Proteinsdat <- proteus::makeProteinTable(pepdat, ncores = 4)
#'
#' data_se<-eset_proteus(Proteins_table=Proteinsdat, MS_metadata=Meta)
#'
#' DEP::plot_numbers(data_se)


eset_proteus<-function(Proteins_table=NULL,
MS_metadata=NULL){

da<-data.frame(Proteins_table$tab, check.names = FALSE)

# Gets columns
Intensity_columns<-c(1:length(colnames(da)))

da$Gene.names<-rownames(da)

# merging with expression -------------------------------------------------
data_unique <- DEP::make_unique(proteins=da,
ids= "Gene.names",
names = "Gene.names", delim = ";")

# Updating meta data: requires label, condition, replicate
MS_metadata$label<-MS_metadata$sample

# adding replicate numbers
MS_metadata$c_c<-as.character(MS_metadata$condition)

condition_ids<-unique(MS_metadata$c_c)

Output_conditions <- setNames(vector(length(condition_ids),
mode="list"), condition_ids)
Output_conditions

MS_metadata$replicate<-c(NA)

# Getting output
for (i in condition_ids) {
max_length<-length(MS_metadata[MS_metadata$c_c == i,]$c_c)
MS_metadata[MS_metadata$c_c == i,]$replicate<-seq(1, max_length)

}

experimental_design <- MS_metadata[c("label", "condition", "replicate")]

# making expressionset ----------------------------------------------------
data_se <- DEP::make_se(proteins_unique = data_unique,
Intensity_columns,
experimental_design)

return(data_se)
}
