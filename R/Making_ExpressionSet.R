#' Converts proteins table and metadata into eset_proteus
#' @export
#' @param Proteins_table Output of makeProteinTable function
#' @param MS_metadata Metadata object prepared by Proteus_Prep function.
#' @param replicate_input a numeric vector indicate replicate number per 
#' condition.
#' @return ExpressSet
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
#' 
#' # Converting sample name into conditions
#' unique(Meta$experiment) # experiment names will be used to set conditions
#' Meta$condition<-gsub("-[1-9]","", Meta$experiment)
#' unique(Meta$condition) # experiment names will be used to set conditions
#' 
#' # measure columns between evidence and meta should match "Intensity"
#' # "Intensity" is case-sensitive
#' colnames(evi_symbols) # "intensity"  
#' unique(Meta$measure) # "Intensity"
#' 
#' measures_var<-c("intensity"="Intensity")
#' 
#' pepdat <- proteus::makePeptideTable(evi_symbols, Meta,
#' measure.cols=measures_var, ncores=2) 
#' 
#' Proteinsdat <- proteus::makeProteinTable(pepdat, ncores=2)
#' 
#' data_se<-eset_proteus(Proteins_table=Proteinsdat, MS_metadata=Meta)
#' 
#' DEP::plot_numbers(data_se)


eset_proteus<-function(Proteins_table=NULL,
MS_metadata=NULL, replicate_input=NULL){

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

if(!is.null(replicate_input)){
    MS_metadata$replicate<-replicate_input
}else if(is.null(replicate_input)){
    
# adding replicate numbers
MS_metadata$c_c<-as.character(MS_metadata$condition)
condition_ids<-unique(MS_metadata$c_c)

Output_conditions <- stats::setNames(vector(length(condition_ids),
mode="list"), condition_ids)
Output_conditions

MS_metadata$replicate<-c(NA)

# Getting output
for (i in condition_ids) {
max_length<-length(MS_metadata[MS_metadata$c_c == i,]$c_c)
MS_metadata[MS_metadata$c_c == i,]$replicate<-seq(1, max_length)
}
}

experimental_design <- MS_metadata[c("label", "condition", "replicate")]

# making expressionset ----------------------------------------------------
data_se <- DEP::make_se(proteins_unique = data_unique,
Intensity_columns,
experimental_design)

return(data_se)
}
