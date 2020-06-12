# ProtGene
This package reads the evidence file from mass spec data and formats it for differential expression analysis. This includes replacing uniprot names with gene symbols for both human and mouse data. Additional option for data imputation is also provided. 

# First install proteus
    install.packages("devtools")
    devtools::install_github("bartongroup/proteusLabelFree")
    devtools::install_github("bartongroup/proteusTMT")
    devtools::install_github("bartongroup/proteusSILAC")
    devtools::install_github("bartongroup/Proteus", 
    build_opts= c("--no-resave-data", "--no-manual"), build_vignettes=FALSE)
    
# Install additional packages
    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

    BiocManager::install(c("org.Hs.eg.db","org.Mm.eg.db", "org.Sc.sgd.db","vsn", "DEP"))

# Install ProtGene
    devtools::install_github("Rvirgenslane/ProtGene")
    
# Run example
    # Getting evidence file
    library(ProtGene)
    library(proteusLabelFree)
    evidenceFile <- system.file("extdata", "evidence.txt.gz", package="proteusLabelFree")
    
    # The directory of the evidence file is required
    Evidence_input <- proteus::readEvidenceFile(evidenceFile)
    
    # extracting uniprot ids: "sp|P00175|CYB2_YEAST" to "P00175
    # Your data may not need this step
    Evidence_input$protein<-gsub(".*[sp|]([^.]+)[|].*", "\\1", Evidence_input$protein)
    head(Evidence_input)
    
    # Converting to gene name
    evi_symbols<-ProtGene::Evidence_Symbols(evi=Evidence_input, Species="Yeast")
    head(evi_symbols)
    
    # Generating metadata from evidence file
    Meta<-Proteus_Prep(evi=evi_symbols)
    head(Meta)
    
    # Converting sample name into conditions
    unique(Meta$experiment) # experiment names will be used to set conditions
    Meta$condition<-gsub("-[1-9]","", Meta$experiment)
    unique(Meta$condition) # experiment names will be used to set conditions
    
    # measure columns between evidence and meta should match "Intensity"
    # "Intensity" is case-sensitive
    colnames(evi_symbols) # "intensity"  
    unique(Meta$measure) # "Intensity"
    
    measures_var<-c("intensity"="Intensity")
    
    pepdat <- proteus::makePeptideTable(evi_symbols, Meta,
        measure.cols=measures_var, ncores=2) 
    
    Proteinsdat <- proteus::makeProteinTable(pepdat, ncores=2)
    
    data_se<-eset_proteus(Proteins_table=Proteinsdat, MS_metadata=Meta)
    
    DEP::plot_numbers(data_se)
    
        
    # filtering data 
    
    data_filt <- filter_missval(data_se, thr = 0)
    plot_frequency(data_filt)
    plot_numbers(data_filt)
    
    # normalization 
    
    Data_norm <- normalize_vsn(data_filt)
    meanSdPlot(assay(Data_norm))
    plot_normalization(data_filt, Data_norm)
    
    # Missing values
    plot_missval(Data_norm)
    
    # imputing data -----------------------------------------------------------
    # Impute missing data using random draws from a 
    # Gaussian distribution centered around a minimal value (for MNAR)
    Data_MinProb_imputation <- DEP::impute(Data_norm, 
    fun = "MinProb", q = 0.1)
    
    plot_imputation(Data_norm, Data_MinProb_imputation)
    
    boxplot(assay(Data_norm))
    boxplot(assay(Data_MinProb_imputation))
    



