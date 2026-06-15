
############ Notes from my errors ########################
#if it crash due to fatal error when you run the mofa
#just uninstall the basilisk package and re-install it
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager::install("basilisk")

############ IF THE ABOVE DOESNT WORK #####################
##Go to terminal and type: conda env create -f mofa_arm_environment.yml
library(reticulate)
use_condaenv("mofa_arm", required = TRUE) ###OR 
use_python("/Users/ilmaqonaah/opt/miniconda3/envs/mofa_arm/bin/python", required = TRUE)
###########verify installation
py_run_string("
deps = {
    'numpy': '1.21.6',
    'scipy': '1.7.1',
    'pandas': '1.3.5',
    'h5py': '3.7.0',
    'sklearn': '1.0.2',
    'mofapy2': '0.7.0'
}

for pkg, ver in deps.items():
    try:
        m = __import__(pkg)
        assert m.__version__ == ver, f'{pkg} version mismatch'
        print(f'✓ {pkg} {ver}')
    except Exception as e:
        print(f'✗ {pkg}: {e}')
")
# Apply MOFA patch
py_run_string("
from mofapy2.run import entry_point
import sys
sys.modules['mofapy2'].run = type('run', (), {'entry_point': entry_point})
")

py_config() #SET USE_BASILISK to FALSE under runmofa
reticulate::py_run_string("
import numpy as np
import scipy
import mofapy2
print(f'numpy: {np.__version__}, scipy: {scipy.__version__}, mofapy2: {mofapy2.__version__}')
") #### check the versions
###########################################################
###########################################################
# FUNCTIONS
###########################################################
###########################################################

setup_project_environment <- function(use_renv = TRUE,
                                      renv_path = "renv") {
  print("Initialising project environment...")
  
  if (use_renv) {
    if (!requireNamespace("renv", quietly = TRUE)) {
      install.packages("renv")
      print("Installed 'renv' for virtual environment support.")
    } else {
      print("'renv' already installed.")
    }
    
    if (!file.exists(file.path(renv_path, "renv.lock"))) {
      print("Initialising new renv environment...")
      renv::init(bare = TRUE)
    } else {
      print("Restoring existing renv environment...")
      renv::restore(prompt = FALSE)
    }
  }
  
  # Ensure BiocManager is available
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
    print("Installed 'BiocManager'.")
  }
  
  # CRAN packages
  cran_packages <- c(
    "data.table",
    "ggplot2",
    "psych",
    "dplyr",
    "ggpubr",
    "GGally",
    "readxl",
    "readr",
    "tidyr",
    "stringi",
    "openxlsx",
    "reshape2",
    "writexl",
    "textclean",
    "stringr",
    "viridis",
    "grid",
    "gtable",
    "gridExtra",
    "gridGraphics",
    "stats",
    "clusterProfiler",
    "pathview",
    "magick",
    "tibble",
    "tidygraph",
    "ggraph",
    "igraph",
    "scales",
    "ggtext",
    "purrr",
    "XML",
    "png",
    "httr",
    "jsonlite",
    "VIM",
    "preprocessCore",
    "ggforce"
  )
  
  # Bioconductor packages
  bioc_packages <- c(
    "tximport",
    "DESeq2",
    "MOFA2",
    "KEGGREST",
    "MSnbase",
    "clusterProfiler",
    "pathview",
    "rgoslin",
    "fgsea",
    "sva",
    "limma",
    "FELLA",
    "biomaRt"
  )
  
  # Install missing CRAN packages
  install_cran <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      print(paste("Installing CRAN package:", pkg))
      install.packages(pkg, dependencies = TRUE)
    } else {
      print(paste("CRAN package already installed:", pkg))
    }
  }
  
  # Install missing Bioconductor packages
  install_bioc <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      print(paste("Installing Bioconductor package:", pkg))
      BiocManager::install(pkg, ask = FALSE, update = TRUE)
    } else {
      print(paste("Bioconductor package already installed:", pkg))
    }
  }
  
  print("Checking and installing CRAN packages...")
  sapply(unique(cran_packages), install_cran)
  
  print("Checking and installing Bioconductor packages...")
  sapply(bioc_packages, install_bioc)
  
  # Load all packages
  print("Loading all required packages into session...")
  invisible(lapply(c(cran_packages, bioc_packages), library, character.only = TRUE))
  
  print("All packages successfully installed and loaded.")
}

################### if it doesnt load the packages ####################
library(data.table)
library(ggplot2)
library(psych)
library(dplyr)
library(ggpubr)
library(GGally)
library(readr)
library(tidyr)
library(stringi)
library(openxlsx)
library(reshape2)
library(writexl)
library(textclean)
library(stringr)
library(viridis)
library(grid)
library(gtable)
library(gridExtra)
library(gridGraphics)
library(stats)
library(clusterProfiler)
library(pathview)
library(magick)
library(tibble)
library(tidygraph)
library(ggraph)
library(igraph)
library(scales)
library(ggtext)
library(purrr)
library(XML)
library(png)
library(httr)
library(jsonlite)
library(VIM)
library(preprocessCore)
library(ggforce)
library(tximport)
library(DESeq2)
library(MOFA2)
library(KEGGREST)
library(MSnbase)
library(clusterProfiler)
library(pathview)
library(rgoslin)
library(fgsea)
library(sva)
library(limma)
library(FELLA)
library(biomaRt)
#library(basilisk)

################### Define your file and sheet names here ##########
input_file <- "Input/All_Data.xlsx"
sheet_names <- c("Transcripts", "Metabolites")
new_names <- c("RNAseq", "LCMS")  # Internal names for MOFA views

###################### Load and Process Data #######################
load_and_clean_input_data <- function(input_file, sheet_names, new_names) {
  print("Loading and processing input data...")
  
  # Read Excel sheets and create named list, creates a data list where each element is a data frame for one omics layer
  print("Reading input sheets...")
  data_list <- setNames(lapply(sheet_names, function(sheet) { #lapply loops through each sheet name
    as.data.frame(readxl::read_excel(input_file, sheet = sheet)) #reads each sheet
  }), sheet_names) 
  
  # Rename the sheets
  print("Renaming sheets...")
  names(data_list) <- new_names
  
  # Define Greek letter replacement function
  replace_greek_letters <- function(text) {
    greek_map <- c(
      "α" = "alpha",
      "β" = "beta",
      "γ" = "gamma",
      "δ" = "delta",
      "ε" = "epsilon",
      "ζ" = "zeta",
      "η" = "eta",
      "θ" = "theta",
      "ι" = "iota",
      "κ" = "kappa",
      "λ" = "lambda",
      "μ" = "mu",
      "ν" = "nu",
      "ξ" = "xi",
      "ο" = "omicron",
      "π" = "pi",
      "ρ" = "rho",
      "σ" = "sigma",
      "τ" = "tau",
      "υ" = "upsilon",
      "φ" = "phi",
      "χ" = "chi",
      "ψ" = "psi",
      "ω" = "omega",
      "ς" = "sigma",
      "%" = "",
      "," = ""
    )
    for (char in names(greek_map)) {
      text <- gsub(char, greek_map[[char]], text, fixed = TRUE)
    }
    return(text)
  }
  
  # Clean 'Name' column
  clean_name_column <- function(df) {
    df$Name <- sapply(df$Name, function(x) {
      x <- replace_greek_letters(x)
      x <- iconv(x,
                 from = "UTF-8",
                 to = "ASCII//TRANSLIT",
                 sub = "")
      return(x)
    })
    return(df)
  }
  data_list <- lapply(data_list, clean_name_column)
  
  # Format dataframes: set row names and drop first column
  adjust_dataframe <- function(df) {
    df <- df %>% `row.names<-`(df$Name) %>% dplyr::select(-1)
    return(df)
  }
  data_list <- lapply(data_list, adjust_dataframe)
  
  # Get the reference column order from the first view
  ref_order <- colnames(data_list[[1]])
  
  # Reorder all views to match the reference
  data_list <- lapply(data_list, function(df)
    df[, ref_order])
  
  print("Data loading and initial cleaning complete.")
  return(data_list)
}

# Load the data
data_list <- load_and_clean_input_data(input_file, sheet_names, new_names)

# Optional: Preview data
str(data_list)
head(data_list$RNAseq)
head(data_list$LCMS)

#I dont need the lipid bit coz i'm doing RNA and metabolite
######### PROCESS LIPID annotations ###########

generate_lipid_annotations <- function(lipidomic_sheets,
                                       data_list,
                                       lipid_annotation_file = "Results/Lipid.annotations.csv") {
  if (length(lipidomic_sheets) == 0) {
    message("No lipidomic sheets provided. Skipping lipid annotation.")
    return(NULL)
  }
  
  print("Processing lipid annotations...")
  
  if (file.exists(lipid_annotation_file)) {
    message("Loading existing lipid annotations from file...")
    
    if (file.info(lipid_annotation_file)$size < 5) {
      print("Existing lipid annotation file is empty. Proceeding without annotations.")
      return(data.frame())
    } else {
      print("Lipid annotations loaded.")
      return(read.csv(lipid_annotation_file))
    }
  }
  
  print("No annotation file found. Building lipid annotations from scratch...")
  
  lipid_annotations <- data.frame(
    name = character(),
    CleanNames = character(),
    ShortAnnotation = character(),
    Description = character(),
    stringsAsFactors = FALSE
  )
  
  for (lipid_sheet in lipidomic_sheets) {
    print(paste("Processing sheet:", lipid_sheet))
    lipid_names <- rownames(data_list[[lipid_sheet]])
    
    for (lipid_name in lipid_names) {
      clean_lipid_name <- gsub("'", "", lipid_name)
      clean_lipid_names <- unlist(strsplit(clean_lipid_name, "/"))
      GoslinDF <- parseLipidNames(clean_lipid_names)
      GoslinDF$RefMet_ID <- NULL
      
      required_cols <- c(
        "Normalized.Name",
        "Message",
        "Grammar",
        "Mass",
        "Sum.Formula",
        "Extended.Species.Name",
        "Lipid.Maps.Main.Class",
        "Lipid.Maps.Category",
        "Functional.Class.Abbr",
        "RefMet_ID"
      )
      
      for (col in required_cols) {
        if (!(col %in% colnames(GoslinDF))) {
          GoslinDF[[col]] <- NA
        }
      }
      
      for (i in 1:nrow(GoslinDF)) {
        original.Name <- GoslinDF[i, "Original.Name"]
        query_url <- paste0(
          "https://www.metabolomicsworkbench.org/rest/refmet/name/",
          URLencode(original.Name, reserved = TRUE),
          "/all"
        )
        response <- httr::GET(query_url)
        
        if (httr::status_code(response) == 200) {
          json_data <- httr::content(response, as = "text", encoding = "UTF-8")
          parsed_data <- jsonlite::fromJSON(json_data, flatten = TRUE)
          
          if (length(parsed_data) > 1) {
            GoslinDF[i, "Normalized.Name"] <- parsed_data$name
            GoslinDF[i, "Message"] <- "NA"
            GoslinDF[i, "Grammar"] <- "refmet"
            GoslinDF[i, "Mass"] <- parsed_data$exactmass
            GoslinDF[i, "Sum.Formula"] <- parsed_data$formula
            GoslinDF[i, "Extended.Species.Name"] <- "ST"
            GoslinDF[i, "Lipid.Maps.Main.Class"] <- parsed_data$main_class
            GoslinDF[i, "Functional.Class.Abbr"] <- parsed_data$sub_class
            GoslinDF[i, "RefMet_ID"] <- parsed_data$refmet_id
          }
        }
      }
      
      GoslinDF$Functional.Class.Abbr <- gsub("\\[|\\]", "", GoslinDF$Functional.Class.Abbr)
      
      GoslinDF <- GoslinDF %>%
        dplyr::summarise(across(everything(), ~ paste(unique(.), collapse = " | ")))
      
      new_annotation <- data.frame(
        name = lipid_name,
        CleanNames = GoslinDF[, "Normalized.Name"],
        ShortAnnotation = GoslinDF[, "Lipid.Maps.Main.Class"],
        Description = paste0(
          "MainClass=[",
          GoslinDF[, "Lipid.Maps.Main.Class"],
          "],ExtendedClass=[",
          GoslinDF[, "Extended.Species.Name"],
          "],FunctionalClass=[",
          GoslinDF[, "Functional.Class.Abbr"],
          "],Mass=[",
          GoslinDF[, "Mass"],
          "],Formula=[",
          GoslinDF[, "Sum.Formula"],
          "]"
        ),
        stringsAsFactors = FALSE
      )
      
      lipid_annotations <- rbind(lipid_annotations, new_annotation)
    }
  }
  
  print("Finalising lipid annotations...")
  lipid_annotations <- unique(lipid_annotations)
  write.csv(lipid_annotations, file = lipid_annotation_file, row.names = FALSE)
  print("Lipid annotations written to file.")
  
  return(lipid_annotations)
}



###################### Set up the species KEGG code ##################
kegg_species <- "taes" ##Triticum aestivum
biomart="plants_mart"
biomart_url="https://plants.ensembl.org"
biomart_dataset="tarefseqv2_eg_gene"

metabolomic_sheets <- c("LCMS")
transcriptomic_sheets <- c("RNAseq")


######################  Get the KEGG metabolite pathway data #######################

generate_kegg_annotation_lists <- function(kegg_species, cache_dir = "Results/kegg_annotation_cache") { #download and organise KEGG pathway annotations 
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE) #create directory prevents crash if no directory exists
  taes_pathways_file <- file.path(cache_dir, paste0(kegg_species, "_pathways.rds"))
  
  if (file.exists(taes_pathways_file)) {
    message("Loading existing KEGG pathway list from file...")
    species_pathways <- readRDS(taes_pathways_file)
    message(paste("Loaded", length(species_pathways), "pathways."))
  } else {
    message("No cached pathway file found. Downloading from KEGG...")
    species_pathways <- KEGGREST::keggList("pathway", kegg_species)
    saveRDS(species_pathways, taes_pathways_file)
    message(paste(
      "Retrieved and saved",
      length(species_pathways),
      "pathways."
    ))
  }
  
  species_ids <- sub("path:", "", names(species_pathways))
  map_ids <- sub(paste0("^", kegg_species), "map", species_ids)
  
  TERM2NAME <- data.frame(
    ID = species_ids,
    Name = as.character(species_pathways),
    stringsAsFactors = FALSE
  )
  
  translator <- data.frame(species = species_ids,
                           map = map_ids,
                           stringsAsFactors = FALSE)
  
  compound2map <- KEGGREST::keggLink("pathway", "compound")
  compound_ids <- sub("^cpd:", "", names(compound2map)) #clean up IDs, remove cpd
  compound_pathways <- sub("^path:", "", compound2map) #clean up IDs, remove path
  
  compound_links <- data.frame(map = compound_pathways,
                               KEGG_ID = compound_ids,
                               stringsAsFactors = FALSE)
  
  compound_links_merged <- merge(translator, compound_links, by = "map")
  TERM2GENE_metabolite <- compound_links_merged[, c("species", "KEGG_ID")]
  colnames(TERM2GENE_metabolite) <- c("Pathway", "KEGG_ID")
  
  gene2pathway_raw <- KEGGREST::keggLink("pathway", kegg_species)
  gene_ids <- sub(paste0("^", kegg_species, ":"), "", names(gene2pathway_raw))
  pathway_ids <- sub("^path:", "", gene2pathway_raw)
  
  TERM2GENE_gene <- data.frame(Pathway = pathway_ids,
                               KEGG_ID = gene_ids,
                               stringsAsFactors = FALSE)
  
  TERM2GENE_list <- list(gene = TERM2GENE_gene, metabolite = TERM2GENE_metabolite)
  TERM2NAME_list <- list(gene = TERM2NAME, metabolite = TERM2NAME)
  
  return(list(TERM2GENE = TERM2GENE_list, TERM2NAME = TERM2NAME_list))
}

#kegg_annotations <- generate_kegg_annotation_lists(kegg_species) #call the function

###################### Define ensembl_to_entrez_mapping ###########################
ensembl <- useMart(
  biomart = "plants_mart",
  dataset = "tarefseqv2_eg_gene",
  host = "https://plants.ensembl.org"
)

# Retrieve mapping: Ensembl Gene ID -> EntrezGene ID
ensembl_to_entrez_mapping <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  mart = ensembl
)

# Preview the result
head(ensembl_to_entrez_mapping)

annotate_transcriptomic_data_with_kegg <- function(transcriptomic_sheets,
                                                   data_list,
                                                   TERM2GENE_list,
                                                   ensembl_to_entrez_mapping,
                                                   kegg_species,
                                                   output_file = "Results/KEGG_gene_annotations.csv",
                                                   overwrite = FALSE) {
  if (length(transcriptomic_sheets) == 0)
    return(NULL)
  print("Annotating transcriptomic data with KEGG IDs using Ensembl to Entrez mapping...")
  
  #if (file.exists(output_file)) {
    #message("Loading existing gene annotations from file.")
    #gene_annotation <- read.csv(output_file)
    #print(paste("Loaded", nrow(gene_annotation), "annotated genes."))
    #return(gene_annotation)
  #}
  ##replacement from chargpt####
  ##dumb bitch tip: if there is preexisting file just delete it manually##
  if (file.exists(output_file) && !overwrite) {
    message("Loading existing gene annotations from file.")
    gene_annotation <- read.csv(output_file)
    print(paste("Loaded", nrow(gene_annotation), "annotated genes."))
    return(gene_annotation)
  }
  
  
  if (is.null(kegg_species)) {
    warning("KEGG species code is NULL. Cannot perform annotation.")
    return(NULL)
  }
  
  # Extract Ensembl gene IDs from all transcriptomic sheets
  gene_list <- unlist(lapply(transcriptomic_sheets, function(sheet) {
    genes <- rownames(data_list[[sheet]])
    gsub(paste0("_", sheet), "", genes)
  }))
  gene_list <- unique(gene_list)
  print(paste(
    "Total unique Ensembl gene identifiers to annotate:",
    length(gene_list)
  ))
  
  # Step 1: Filter mapping for input genes
  mapping_filtered <- ensembl_to_entrez_mapping[
    ensembl_to_entrez_mapping$ensembl_gene_id %in% gene_list &
                                                  !is.na(ensembl_to_entrez_mapping$entrezgene_id), ]
  
  # Step 2: Prepare KEGG TERM2GENE table
  kegg_df <- TERM2GENE_list$gene
  kegg_df$EntrezID <- sub(paste0("^", kegg_species, ":"), "", kegg_df$KEGG_ID)
  
  # Step 3: Merge on Entrez ID to get KEGG mappings
  gene_annotation <- merge(mapping_filtered, kegg_df, by.x = "entrezgene_id", by.y = "EntrezID")
  
  # Step 4: Format and save
  gene_annotation <- gene_annotation[, c("ensembl_gene_id", "KEGG_ID")]
  colnames(gene_annotation) <- c("gene", "KEGG_ID")
  
  print(paste("Matched KEGG IDs for", nrow(gene_annotation), "genes."))
  write.csv(gene_annotation, output_file, row.names = FALSE)
  print("Gene annotations saved to file.")
  
  return(gene_annotation)
}

################# ILMA call function ##############################
transcriptomic_annot <- annotate_transcriptomic_data_with_kegg(
  transcriptomic_sheets = transcriptomic_sheets,
  data_list = data_list,
  TERM2GENE_list = kegg_annotations$TERM2GENE,
  ensembl_to_entrez_mapping = ensembl_to_entrez_mapping,
  kegg_species = kegg_species,
  output_file = "Results/KEGG_gene_annotations.csv",
  overwrite = TRUE
)

annotate_metabolomic_data_with_kegg <- function(metabolomic_sheets,
                                                data_list,
                                                output_file = "Results/KEGG_metabolite_annotations.csv") {
  if (length(metabolomic_sheets) == 0)
    return(NULL)
  print("Annotating metabolite data with KEGG IDs via direct KEGG queries...")
  
  if (file.exists(output_file)) {
    message("Loading existing KEGG metabolite annotations from file.")
    metabolite_annotation <- read.csv(output_file)
    print(paste(
      "Loaded",
      nrow(metabolite_annotation),
      "annotated metabolites."
    ))
    return(metabolite_annotation)
  }
  
  # Step 1: Clean metabolite names
  compound_list <- unlist(lapply(metabolomic_sheets, function(sheet) {
    compounds <- rownames(data_list[[sheet]])
    compounds <- gsub(paste0("_", sheet), "", compounds)
    compounds <- trimws(sub(".*\\|", "", compounds))  # Keep part after pipe
    return(compounds)
  }))
  compound_list <- unique(compound_list)
  print(paste("Total unique cleaned compound names:", length(compound_list)))
  
  # Step 2: Search KEGG for each compound
  print("Querying KEGG for compound identifiers...")
  
  get_kegg_id <- function(metabolite_name, max_retries = 2) {
    attempt <- 1
    while (attempt <= max_retries) {
      result <- tryCatch({
        KEGGREST::keggFind("compound", metabolite_name)
      }, error = function(e) {
        return(NULL)
      })
      if (!is.null(result) && length(result) > 0) {
        return(names(result)[1])  # First match
      }
      attempt <- attempt + 1
    }
    return(NA)
  }
  
  # Add progress bar
  n <- length(compound_list)
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  
  kegg_ids <- vector("character", n)
  for (i in seq_along(compound_list)) {
    kegg_ids[i] <- get_kegg_id(compound_list[i])
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Combine results
  metabolite_annotation <- data.frame(
    Metabolite = compound_list,
    KEGG_ID = kegg_ids,
    stringsAsFactors = FALSE
  )
  
  # Filter and continue as before
  
  # Filter out compounds without a KEGG match
  before_filter <- nrow(metabolite_annotation)
  metabolite_annotation <- metabolite_annotation[!is.na(metabolite_annotation$KEGG_ID), ]
  after_filter <- nrow(metabolite_annotation)
  
  print(paste(
    "Annotated",
    after_filter,
    "of",
    before_filter,
    "metabolites."
  ))
  
  # Clean KEGG ID formatting
  metabolite_annotation$KEGG_ID <- sub("^cpd:", "", metabolite_annotation$KEGG_ID)
  
  # Save and return
  write.csv(metabolite_annotation, output_file, row.names = FALSE)
  print("KEGG metabolite annotations written to file.")
  
  return(metabolite_annotation)
}

################# ILMA call function #############################
metabolite_annot <- annotate_metabolomic_data_with_kegg(metabolomic_sheets = metabolomic_sheets,
                                                        data_list = data_list,
                                                        output_file = "Results/KEGG_metabolite_annotations.csv")

##############################
# Get kegg graph data for FELLA analysis
##############################
load_or_build_fella_database <- function(kegg_species, fella_dir = "Results/fella_db") {
  if (is.null(kegg_species)) {
    warning("KEGG species code is NULL. Cannot proceed with FELLA initialisation.")
    return(NULL)
  }
  
  # Define paths
  fella_db_dir <- file.path(fella_dir, kegg_species)
  graph_file   <- file.path(fella_dir, "kegg_graph.rds")
  
  # Ensure base directory exists
  if (!dir.exists(fella_dir)) {
    dir.create(fella_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Load or create KEGG graph
  if (file.exists(graph_file)) {
    message("Loading KEGG graph from disk...")
    kegg_graph <- readRDS(graph_file)
  } else {
    message("Building KEGG graph from KEGG REST...")
    kegg_graph <- FELLA::buildGraphFromKEGGREST(organism = kegg_species)
    saveRDS(kegg_graph, graph_file)
  }
  
  # Load or build FELLA database
  if (dir.exists(fella_db_dir) &&
      file.exists(file.path(fella_db_dir, "diffusion.matrix.RData"))) {
    message("Loading FELLA database from disk...")
    fella.data <- FELLA::loadKEGGdata(
      databaseDir = fella_db_dir,
      internalDir = FALSE,
      loadMatrix  = "diffusion"
    )
  } else {
    message("Building FELLA database...")
    FELLA::buildDataFromGraph(
      keggdata.graph = kegg_graph,
      databaseDir    = fella_db_dir,
      internalDir    = FALSE,
      matrices       = c("diffusion"),
      # Add "pagerank" if needed
      normality      = c("diffusion"),
      niter          = 100
    )
    
    fella.data <- FELLA::loadKEGGdata(
      databaseDir = fella_db_dir,
      internalDir = FALSE,
      loadMatrix  = "diffusion"
    )
  }
  return(fella.data)
}



##############################
# Normalise data
##############################

normalise_data_mofa <- function(data_list,
                                log_base = 2,
                                na_limit = 1e-7) {
  print("Normalising data...")
  
  view_configs <- list(
    "RNAseq"   = list(
      log = TRUE,
      treat_zeros_as_NA = TRUE,
      library_scaled = TRUE
    ),
    "LCMS"     = list(
      log = TRUE,
      treat_zeros_as_NA = TRUE,
      library_scaled = TRUE
    )
  )
  
  match_config <- function(view) {
    matched <- names(view_configs)[startsWith(view, names(view_configs))]
    if (length(matched) > 0)
      return(view_configs[[matched[1]]])
    return(list(
      log = TRUE,
      treat_zeros_as_NA = TRUE,
      library_scaled = TRUE
    ))
  }
  
  log_transform <- function(mat, config, log_base, na_limit) {
    mat <- as.matrix(mat)
    mode(mat) <- "numeric"
    
    if (config$log) {
      if (config$treat_zeros_as_NA) {
        mat[abs(mat) <= na_limit] <- NA
      } else {
        mat[mat < 0] <- NA
        if (any(mat > 0, na.rm = TRUE)) {
          min_nonzero <- min(mat[mat > 0], na.rm = TRUE)
          pseudo_count <- 0.1 * min_nonzero
          mat <- mat + pseudo_count
        }
      }
      log_fn <- switch(
        as.character(log_base),
        "2" = log2,
        "10" = log10,
        "exp(1)" = log,
        stop("Unsupported log base.")
      )
      mat <- suppressWarnings(log_fn(mat))
      mat[!is.finite(mat)] <- NA
    }
    return(mat)
  }
  
  scale_to_zscore <- function(mat) {
    scaled <- t(apply(mat, 1, function(x) {
      if (sum(!is.na(x)) > 1)
        scale(x)
      else
        rep(NA, length(x))
    }))
    rownames(scaled) <- rownames(mat)
    colnames(scaled) <- colnames(mat)
    scaled
  }
  
  remove_duplicated_rows_by_name <- function(mat, verbose = TRUE) {
    if (!is.matrix(mat) && !is.data.frame(mat)) {
      stop("Input must be a matrix or data.frame")
    }
    
    rown <- rownames(mat)
    if (is.null(rown)) {
      stop("Matrix has no row names")
    }
    
    dup_idx <- duplicated(rown)
    n_dup <- sum(dup_idx)
    
    if (n_dup > 0 && verbose) {
      message("Removing ",
              n_dup,
              " duplicated rows: ",
              paste(unique(rown[dup_idx]), collapse = ", "))
    }
    
    mat[!dup_idx, , drop = FALSE]
  }
  
  result <- list()
  for (view_name in names(data_list)) {
    config <- match_config(view_name)
    
    mat <- data_list[[view_name]]
    mat_log <- log_transform(mat, config, log_base, na_limit)
    mat_scaled <- scale_to_zscore(mat_log)
    rownames(mat_scaled) <- paste0(rownames(mat_scaled), "_", view_name)
    mat_clean <- remove_duplicated_rows_by_name(mat_scaled)
    
    result[[view_name]] <- mat_clean
    message(
      sprintf(
        "Processed view '%s': %d → %d features after filtering.",
        view_name,
        nrow(mat),
        nrow(mat_clean)
      )
    )
  }
  return(result)
}
###################### ILMA call function 
#normalised_data <- normalise_data_mofa(data_list)

######################  Load and run MOFA #######################

initialise_and_run_mofa <- function(mofa_list, output_file = "Results/Full_MOFA_model.hdf5") {
  print("Initialising and running MOFA...")
  
  message("Creating MOFA object...")
  MOFAobject <- create_mofa(mofa_list)
  
  message("Setting MOFA options...")
  data_opts  <- get_default_data_options(MOFAobject)
  model_opts <- get_default_model_options(MOFAobject)
  
  if (model_opts$num_factors >= 6) {
    model_opts$num_factors <- 6
  }
  
  train_opts <- get_default_training_options(MOFAobject)
  
  train_opts$convergence_mode <- "slow"
  train_opts$maxiter <- 10000
  
  message("Preparing MOFA model...")
  MOFAobject <- prepare_mofa(
    MOFAobject,
    data_options     = data_opts,
    model_options    = model_opts,
    training_options = train_opts
  )
  
  if (file.exists(output_file)) {
    message("Loading existing MOFA model from file...")
    MOFAobject.trained <- load_model(output_file)
  } else {
    message("No existing MOFA model found. Training new model...")
    MOFAobject.trained <- run_mofa(MOFAobject, outfile = output_file, use_basilisk = FALSE)
    message("MOFA training complete.")
  }
  
  message("Annotating metadata with sample categories...")
  MOFAobject.trained@samples_metadata$category <- as.factor(MOFAobject.trained@samples_metadata$group)
  
  return(MOFAobject.trained)
}

######################  Post-process #######################

postprocess_mofa_model <- function(MOFAobject.trained,
                                   kegg_species,
                                   metabolite_annotation = NULL,
                                   gene_annotation = NULL,
                                   transcriptomic_sheets = character(),
                                   metabolomic_sheets = character()) {
  message("Extracting and scaling MOFA factor weights...")
  MOFAobject.trained@expectations$adjusted_W <- get_weights(MOFAobject.trained, scale = TRUE)
  print("MOFA initialisation and training complete.")
  
  message("Extracting KEGG-only weights from MOFA object...")
  weights <- MOFAobject.trained@expectations$adjusted_W
  annotation_maps <- list()
  
  if (!is.null(metabolite_annotation)) {
    message("Using metabolite annotations.")
    metabolite_annotation$Metabolite <- toupper(metabolite_annotation$Metabolite)
    map <- metabolite_annotation[!is.na(metabolite_annotation$KEGG_ID), ]
    map <- map %>% dplyr::distinct(Metabolite, .keep_all = TRUE)
    annotation_maps$metabolite <- map
  }
  
  if (!is.null(gene_annotation)) {
    message("Using gene annotations.")
    gene_annotation$Gene <- toupper(gene_annotation$gene)
    map <- gene_annotation[!is.na(gene_annotation$KEGG_ID), ]
    map <- map %>% dplyr::distinct(Gene, .keep_all = TRUE)
    annotation_maps$gene <- map
  }
  
  filtered_weights <- list()
  
  for (df_name in names(weights)) {
    df <- weights[[df_name]]
    clean_names <- toupper(trimws(gsub(
      "_.*$", "", sub("^.*\\|", "", rownames(df))
    )))
    
    row_map <- data.frame(
      original = rownames(df),
      cleaned = clean_names,
      stringsAsFactors = FALSE
    )
    
    matched_rows <- NULL
    
    if (df_name %in% metabolomic_sheets &&
        "metabolite" %in% names(annotation_maps)) {
      ann <- annotation_maps$metabolite
      by_col <- "Metabolite"
    } else if (df_name %in% transcriptomic_sheets &&
               "gene" %in% names(annotation_maps)) {
      ann <- annotation_maps$gene
      by_col <- "Gene"
    } else {
      ann <- NULL
    }
    
    if (!is.null(ann)) {
      matched_rows <- merge(row_map,
                            ann,
                            by.x = "cleaned",
                            by.y = by_col,
                            all.x = FALSE)
      message(sprintf(
        "View: %s | Matched: %d / %d",
        df_name,
        nrow(matched_rows),
        nrow(df)
      ))
    }
    
    if (!is.null(matched_rows) && nrow(matched_rows) > 0) {
      df_filtered <- df[matched_rows$original, , drop = FALSE]
      rownames(df_filtered) <- gsub(paste0(kegg_species, ":"), "", matched_rows$KEGG_ID)
      filtered_weights[[df_name]] <- df_filtered
    } else {
      message(sprintf("No KEGG mapping found for view: %s", df_name))
    }
  }
  
  MOFAobject.trained@expectations$KEGG_weights <- filtered_weights
  print("Finished extracting and storing KEGG-only weights.")
  
  message("Generating pairwise p-values across all views...")
  view_names <- names(MOFAobject.trained@data)
  
  MOFAobject.trained@expectations$Pvalues <- list()
  MOFAobject.trained@expectations$FDRPvalues <- list()
  MOFAobject.trained@expectations$logFC <- list()
  
  for (view_name in view_names) {
    message(sprintf("Processing view: %s", view_name))
    raw_data <- as.data.frame(MOFAobject.trained@data[[view_name]]$group1)
    colnames(raw_data) <- sub("_\\d+$", "", colnames(raw_data))
    
    group_labels <- factor(colnames(raw_data))
    design <- model.matrix( ~ 0 + group_labels)
    colnames(design) <- levels(group_labels)
    
    if (ncol(design) < 2) {
      message(sprintf("Skipping view '%s' — fewer than 2 sample groups.", view_name))
      MOFAobject.trained@expectations$Pvalues[[view_name]] <- NULL
      next
    }
    
    message(sprintf("Fitting linear model for view: %s", view_name))
    fit <- limma::lmFit(raw_data, design)
    comparisons <- combn(colnames(design), 2)
    contrast_names <- apply(comparisons, 2, function(x)
      paste(x, collapse = " vs "))
    contrast_defs <- apply(comparisons, 2, function(x)
      paste0(x[1], " - ", x[2]))
    contrast_matrix <- limma::makeContrasts(contrasts = contrast_defs, levels = design)
    
    fit2 <- limma::contrasts.fit(fit, contrast_matrix)
    fit2 <- limma::eBayes(fit2)
    
    pvals_df <- data.frame(row.names = rownames(raw_data))
    fdr_df <- data.frame(row.names = rownames(raw_data))
    logfc_df <- data.frame(row.names = rownames(raw_data))
    
    for (i in seq_along(contrast_names)) {
      tt <- limma::topTable(fit2,
                            coef = i,
                            number = Inf,
                            sort.by = "none")
      pvals_df[[paste0(contrast_names[i], " Pval")]] <- tt$P.Value
      fdr_df[[paste0(contrast_names[i], " FDR Pval")]] <- tt$adj.P.Val
      logfc_df[[paste0(contrast_names[i], " logFC")]] <- tt$logFC
    }
    
    MOFAobject.trained@expectations$Pvalues[[view_name]] <- pvals_df
    MOFAobject.trained@expectations$FDRPvalues[[view_name]] <- fdr_df
    MOFAobject.trained@expectations$logFC[[view_name]] <- logfc_df
  }
  
  print("Finished generating pairwise p-values for all views.")
  return(MOFAobject.trained)
}

######################  Add metadata to the experiment #######################
generate_sample_metadata <- function(view_names) {
  message("Generating metadata from sample names...")
  metadata <- data.frame(matrix(ncol = 0, nrow = length(view_names)))
  rownames(metadata) <- view_names
  metadata$sample <- rownames(metadata)
  metadata$category <- substr(metadata$sample, 1, nchar(metadata$sample) - 3)
  metadata$category <- factor(metadata$category, levels = unique(metadata$category))
  return(metadata)
}

assign_sample_metadata_to_mofa <- function(MOFAobject.trained, metadata) {
  message("Assigning group and binary indicators to MOFA object...")
  MOFAobject.trained@samples_metadata$group <- metadata$category
  for (cat in levels(metadata$category)) {
    metadata[[cat]] <- as.integer(metadata$category == cat)
  }
  samples_metadata(MOFAobject.trained) <- metadata
  return(MOFAobject.trained)
}

define_group_colours <- function(metadata) {
  category.colors <- c(
    #"#bcf60c",
    #"#3cb44b",
    #"#ffe119",
    #"#f58231",
    #"#99c0f0",
    #"#46f0f0",
    #"#911eb4",
    #"#3a13eb",
    #"#e6194b",
    #"#f032e6",
    #"#fabebe",
    #"#008080",
    #"#e6beff",
    #"#9a6324",
    #"#fffac8",
    "#ea9855",
    "#ffce71",
    "#343488",
    "#828be6",
    #"#5861bd",
    "#902f69",
    "#c26fad"
  )
  n <- length(levels(metadata$category))
  group.colors <- category.colors[1:n]
  names(group.colors) <- levels(metadata$category)
  return(group.colors)
}

define_view_colours <- function(view_names) {
  category.colors <- c(
    #"#bcf60c",
    #"#3cb44b",
    #"#ffe119",
    #"#f58231",
    #"#99c0f0",
    #"#46f0f0",
    #"#911eb4",
    #"#3a13eb",
    #"#e6194b",
    #"#f032e6",
    #"#fabebe",
    #"#008080",
    #"#e6beff",
    #"#9a6324",
    #"#fffac8",
    ##ea9855",
    "#ffce71",
    #"#343488",
    "#828be6",
    #"#5861bd",
    #"#902f69",
    "#c26fad"
    )
  n <- length(view_names)
  view.colors <- category.colors[1:n]
  names(view.colors) <- view_names
  return(view.colors)
}

standardise_feature_names <- function(MOFAobject.trained) {
  message("Standardising feature names...")
  MOFAobject.trained@features_metadata$cleanNames <- mapply(
    function(f, v) {
      sub(paste0("_", v, "$"), "", f)
    },
    MOFAobject.trained@features_metadata$feature,
    MOFAobject.trained@features_metadata$view
  )
  
  MOFAobject.trained.cleanNames <- MOFAobject.trained
  MOFAobject.trained.cleanNames@expectations$W <- lapply(names(MOFAobject.trained.cleanNames@expectations$W), function(view) {
    W <- MOFAobject.trained.cleanNames@expectations$W[[view]]
    rownames(W) <- sub(paste0("_", view, "$"), "", rownames(W))
    W
  })
  names(MOFAobject.trained.cleanNames@expectations$W) <- names(MOFAobject.trained@expectations$W)
  MOFAobject.trained.cleanNames@features_metadata$feature <- MOFAobject.trained@features_metadata$cleanNames
  return(MOFAobject.trained.cleanNames)
}



########################################
# VARIANCE PLOTS
########################################

generate_mofa_variance_plots <- function(MOFAobject.trained,
                                         view.colors,
                                         output_path = "Results/Factor_and_view_variance.png") {
  if (file.exists(output_path)) {
    return(invisible(NULL))
  }
  
  print("Generating variance explained plots from MOFA object...")
  
  # Get variance data
  message("Retrieving variance explained from MOFA model...")
  variance_data <- get_variance_explained(MOFAobject.trained)
  
  # Create group vs view variance data frame
  message("Preparing factor-wise variance table...")
  group_vs_view_variance_df <- as.data.frame(variance_data$r2_per_factor$group1)
  number_of_factors <- nrow(group_vs_view_variance_df)
  group_vs_view_variance_df$Factor <- rownames(group_vs_view_variance_df)
  variance_long <- reshape2::melt(
    as.data.frame(group_vs_view_variance_df),
    id.vars = "Factor",
    variable.name = "View",
    value.name = "Variance"
  )
  view_ratio <- (1 / (number_of_factors + 1) * 100)
  variance_long$Variance <- variance_long$Variance / number_of_factors
  
  # Prepare view-wise total variance explained
  message("Preparing view-level total variance table...")
  view_variance_df <- as.data.frame(variance_data$r2_total$group1)
  view_variance_df$View <- rownames(view_variance_df)
  view_variance_long <- reshape2::melt(
    view_variance_df,
    id.vars = "View",
    variable.name = "Factor",
    value.name = "VarianceExplained"
  )
  view_variance_long$Variance <- view_variance_long$VarianceExplained / number_of_factors
  
  # Sort view names alphabetically
  ordered_views <- sort(unique(variance_long$View))
  variance_long$View <- factor(variance_long$View, levels = ordered_views)
  view_variance_long$View <- factor(view_variance_long$View, levels = ordered_views)
  
  # Plot 1: Heatmap of factor-wise variance
  message("Creating factor-wise variance heatmap...")
  Totalvariance <- ggplot(variance_long, aes(x = View, y = Factor, fill = Variance)) +
    geom_tile() +
    scale_fill_gradient(low = "white",
                        high = "purple",
                        name = "Variance (R²)") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 22 #previously 12
      ),
      axis.text.y = element_text(size = 22), #previously 12
      plot.margin = margin(
        t = -10,
        r = 45,
        b = 0,
        l = 0
      ),
      text = element_text(size = 24), #previously 12
      panel.grid = element_blank()
    ) +
    labs(x = "MOFA Views", y = "MOFA Factors")
  
  # Plot 2: Barplot of total variance per view
  message("Creating total view-wise variance bar plot...")
  wheatvarExp <- ggplot(view_variance_long, aes(x = View, y = Variance, fill = View)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = view.colors) +
    theme_minimal() +
    theme(
      text = element_text(size = 24), #previously 15
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 22), #previously 12
      plot.margin = margin(
        t = 10,
        r = 10,
        b = 0,
        l = 25
      ),
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 22 #previously 12
      )
    ) +
    labs(
      title = "Variance Explained by MOFA Views and Factors",
      x = "",
      y = "Variance (R²)",
      fill = "View"
    )
  
  # Combine and save the plots
  message("Saving combined variance plot to disk...")
  combined_plot <- gridExtra::grid.arrange(wheatvarExp,
                                           Totalvariance,
                                           ncol = 1,
                                           nrow = 2)
  
  ggsave(
    filename = output_path,
    plot = combined_plot,
    width = 14, #previously 10 
    height = 15
  )
  
  # Ensure graphics devices are closed
  while (dev.cur() > 1)
    dev.off()
  
  print("Variance explained plots saved.")
}


########################################
# PCA PLOTS
########################################

generate_mofa_factor_matrix_plot <- function(MOFAobject.trained,
                                             view.colors,
                                             group.colors,
                                             output_path = "Results/Factor_Pair_PCA_Matrix.png") {
  if (file.exists(output_path)) {
    message("PCA plot matrix already exists. Skipping plot generation.")
    return(invisible(NULL))
  }
  
  print("Generating PCA matrix of MOFA factors...")
  
  # Extract factor numbers
  message("Extracting factor numbers...")
  factor_numbers <- 1:ncol(get_weights(MOFAobject.trained)[[1]])
  
  # Load feature weights
  message("Loading weights from MOFA object...")
  weights <- get_weights(MOFAobject.trained)
  
  # Convert weights into a long format for easier processing
  message("Melting weights into long format...")
  combined_weights <- do.call(rbind, lapply(names(weights), function(view_name) {
    df <- as.data.frame(weights[[view_name]])
    df$Feature <- rownames(df)
    df$View <- view_name
    return(df)
  }))
  melted_weights <- reshape2::melt(
    combined_weights,
    id.vars = c("Feature", "View"),
    variable.name = "Factor",
    value.name = "Weight"
  )
  
  # Reshape into wide format for feature scatter plotting
  message("Creating wide-format feature matrix...")
  feature_matrix <- reshape2::dcast(melted_weights, Feature + View ~ Factor, value.var = "Weight")
  
  # Initialise matrix to store plots
  plot_matrix <- matrix(list(),
                        nrow = length(factor_numbers),
                        ncol = length(factor_numbers))
  
  FactorZ <- MOFAobject.trained@expectations$Z$group1
  FactorZ <- as.data.frame(FactorZ)
  FactorZ$sample <- rownames(FactorZ)  # <- set rownames as a column
  FactorZ <- inner_join(metadata, FactorZ, by = c("sample" = "sample"))  # <- now join
  
  
  # Populate matrix with appropriate plots
  message("Populating plot matrix...")
  for (i in factor_numbers) {
    for (j in factor_numbers) {
      if (i < j) {
        plot_df <- FactorZ[, c("sample",
                               "category",
                               paste0("Factor", i),
                               paste0("Factor", j))]
        
        # Ensure 'category' is a factor and matches group.colors
        plot_df$category <- factor(plot_df$category, levels = names(group.colors))
        
        plot_matrix[[i, j]] <- ggplot(
          plot_df,
          aes_string(
            x = paste0("Factor", i),
            y = paste0("Factor", j),
            color = "category",
            fill = "category"
          )
        ) +
          theme_classic()+
          theme(
            axis.text = element_text(size = 18), #previously 8
            axis.title = element_text(size = 20), #previously 10
            legend.text = element_text(size = 18), #previously 8
            legend.title = element_text(size = 20), #previously 10
            plot.title = element_text(size = 20) 
          )+
          geom_point(size = 4) +
          geom_mark_ellipse(aes(group = category, label = NULL),
                            alpha = 0.25,
                            show.legend = FALSE) +
          scale_color_manual(values = group.colors) +
          scale_fill_manual(values = group.colors) +
          labs(title = paste("Sample PCA: Factor", i, "vs", j))
      } else if (i == j) {
        # Diagonal: density plots
        factor_values <- as.data.frame(get_factors(MOFAobject.trained)$group1[, i])
        colnames(factor_values) <- "FactorValue"
        factor_values$Sample <- rownames(factor_values)
        factor_values$SampleType <- stringr::str_replace(factor_values$Sample, "_\\d+$", "")
        
        plot_matrix[[i, j]] <- ggplot(factor_values, aes(x = FactorValue, fill = SampleType)) +
          geom_density(alpha = 0.6) +
          scale_fill_manual(values = group.colors) +
          theme_classic() +
          labs(
            title = paste("Density Plot: Factor", i),
            x = "Factor Value",
            y = "Density"
          ) +
          theme(
            axis.text = element_text(size = 18), #previously 8
            axis.title = element_text(size = 20), #previously 10
            legend.text = element_text(size = 18), #previously 8
            legend.title = element_text(size = 20),#previously 10
            plot.title = element_text(size = 20)
          )
      } else {
        # Lower triangle: feature weight scatter plots
        weight_data <- feature_matrix[, c("Feature",
                                          "View",
                                          paste0("Factor", i),
                                          paste0("Factor", j))]
        colnames(weight_data) <- c("Feature", "View", "FactorX", "FactorY")
        
        plot_matrix[[i, j]] <- ggplot(weight_data, aes(
          x = FactorX,
          y = FactorY,
          color = View
        )) +
          geom_point(size = 2, alpha = 0.7) +
          scale_color_manual(values = view.colors) +
          stat_ellipse(
            aes(color = View),
            geom = "polygon",
            alpha = 0.25,
            show.legend = FALSE
          ) +
          theme_classic() +
          labs(
            title = paste("Feature Weights: Factor", i, "vs", j),
            x = paste("Factor", i),
            y = paste("Factor", j)
          ) +
          theme(
            axis.text = element_text(size = 18), #previously 8
            axis.title = element_text(size = 20), #previously 10
            legend.text = element_text(size = 18), #previously 8
            legend.title = element_text(size = 20), #previously 10
            plot.title = element_text(size = 20)
          )
      }
    }
  }
  
  # Flatten plot matrix to list of grobs
  plot_list <- as.vector(plot_matrix)
  
  # Save as high-resolution PNG
  message("Saving PCA plot matrix to disk...")
  png(
    output_path,
    width = 35,
    height = 25,
    units = "in",
    res = 300
  )
  par(mar = c(5, 5, 5, 5))
  gridExtra::grid.arrange(
    grobs = plot_list,
    nrow = length(factor_numbers),
    ncol = length(factor_numbers)
  )
  
  # Close devices
  while (dev.cur() > 1)
    dev.off()
  
  print(paste("PCA plot matrix saved to", output_path))
}


########################################
### FACTOR COMPARISONS
########################################
generate_mofa_factor_comparison_plot <- function(MOFAobject.trained,
                                                 group.colors,
                                                 view.colors,
                                                 output_path = "Results/Factor_loading.png") {
  if (file.exists(output_path)) {
    message("Factor comparison plot already exists. Skipping plot generation.")
    return(invisible(NULL))
  }
  
  print("Creating MOFA factor comparison plots...")
  
  factor_numbers <- 1:ncol(get_weights(MOFAobject.trained)[[1]])
  plot_matrix <- matrix(list(), nrow = 2, ncol = length(factor_numbers))
  
  for (j in factor_numbers) {
    print(paste0("Processing factor ", j, "..."))
    
    # Extract factor values
    factor_values <- as.data.frame(get_factors(MOFAobject.trained)$group1[, j])
    colnames(factor_values) <- "FactorValue"
    factor_values$Sample <- rownames(factor_values)
    factor_values$SampleType <- stringr::str_replace(factor_values$Sample, "_\\d+$", "")
    
    # Row 1: Violin plot
    message(sprintf("Generating violin plot for Factor %d...", j))
    plot_matrix[[1, j]] <- ggplot(factor_values,
                                  aes(x = SampleType, y = FactorValue, fill = SampleType)) +
      geom_violin(trim = FALSE, alpha = 0.6) +
      geom_jitter(
        shape = 16,
        position = position_jitter(0.2),
        alpha = 0.5,
        size = 1
      ) +
      scale_fill_manual(values = group.colors) +
      theme_classic() +
      labs(title = paste("Factor", j),
           x = NULL,
           y = "Factor Value") +
      theme(
        axis.text = element_text(size = 12), #previously 8
        axis.title = element_text(size = 14), #previously 10
        legend.text = element_text(size = 12), #previously 8
        legend.title = element_text(size = 14) #previously 10
      )
    
    # Row 2: Feature mean vs weight
    message(sprintf("Generating feature-weight scatter plot for Factor %d...", j))
    views <- views_names(MOFAobject.trained)
    scatter_list <- vector("list", length(views))
    
    for (i in seq_along(views)) {
      view <- views[i]
      raw_data <- MOFAobject.trained@data[[view]]$group1
      weights  <- MOFAobject.trained@expectations$W[[view]][, j]
      feats    <- intersect(rownames(raw_data), names(weights))
      submat   <- raw_data[feats, , drop = FALSE]
      
      mean_mat <- sapply(feats, function(feat) {
        vals <- as.numeric(submat[feat, ])
        tapply(vals, factor_values$SampleType, mean, na.rm = TRUE)
      })
      
      mean_df <- as.data.frame(mean_mat, stringsAsFactors = FALSE)
      mean_df$SampleType <- rownames(mean_df)
      
      stacked <- stack(mean_df[feats])
      
      df <- data.frame(
        SampleType  = rep(mean_df$SampleType, times = length(feats)),
        Feature     = stacked$ind,
        FeatureMean = stacked$values,
        Weight      = abs(weights[stacked$ind]),
        view        = view,
        stringsAsFactors = FALSE
      )
      
      scatter_list[[i]] <- df
    }
    
    scatter_data <- do.call(rbind, scatter_list)
    
    plot_matrix[[2, j]] <- ggplot(scatter_data, aes(x = FeatureMean, y = Weight, color = view)) +
      geom_point(alpha = 0.8, size = 2) +
      theme_classic() +
      scale_color_manual(values = view.colors) +
      facet_wrap( ~ SampleType, scales = "free_x", nrow = 1) +
      labs(
        title = paste("Mean Feature Values vs Factor", j, "Weights"),
        y = paste("Factor", j, "Weight")
      ) +
      theme(
        axis.text = element_text(size = 12), #previously 8
        axis.title = element_text(size = 14), #previously 10
        legend.text = element_text(size = 12), #previously 8
        legend.title = element_text(size = 14) #previously 10
      )
  }
  
  # Assemble and save
  print("Assembling plot grid...")
  plot_list <- as.vector(plot_matrix)
  grid_plot <- gridExtra::arrangeGrob(
    grobs = plot_list,
    nrow = length(factor_numbers),
    ncol = 2,
    top = grid::textGrob("FACTOR vs Sample", gp = grid::gpar(
      fontsize = 14, fontface = "bold"
    )),
    padding = unit(30, "pt")
  )
  
  print("Saving factor comparison figure...")
  dynamic_width <- length(group.colors) * 2 + 2
  png(
    output_path,
    width = dynamic_width,
    height = 25, #previously 20
    units = "in",
    res = 300
  )
  grid::grid.draw(grid_plot)
  
  while (dev.cur() > 1)
    dev.off()
  
  print(paste("Factor comparison plot saved to", output_path))
}


########################################
# FACTOR CORRELATION PLOTS
########################################
generate_mofa_factor_correlation_plot <- function(MOFAobject.trained,
                                                  view.colors,
                                                  output_path = "Results/Correlation_analysis.png") {
  if (file.exists(output_path)) {
    message("Factor correlation plot already exists. Skipping plot generation.")
    return(invisible(NULL))
  }
  
  print("Generating MOFA factor correlation plots...")
  
  # Open PNG device
  png(
    filename = output_path,
    width = 9,
    height = 5,
    units = "in",
    res = 300
  )
  
  # Adjust layout for two plots with outer margin space
  par(
    mfrow = c(1, 2),
    mar = c(5, 4, 3, 2),
    oma = c(1, 0, 2, 0)
  )  # Outer margin added
  
  # Plot 1: Factor-to-Factor correlation
  message("Rendering factor-to-factor correlation heatmap...")
  plot_factor_cor(
    MOFAobject.trained,
    title = "Pearson's Correlation Coefficient",
    mar = c(1, 0, 1, 0),
    tl.cex = 1.2
  )
  mtext("Factor vs Factor", side = 1, line = 4)
  
  # Plot 2: Factor vs Covariate correlation
  message("Rendering factor-to-sample group correlation heatmap...")
  treatments <- levels(MOFAobject.trained@samples_metadata$category)
  correlate_factors_with_covariates(
    MOFAobject.trained,
    covariates = treatments,
    plot = "r",
    mar = c(1, 0, 1, 0),
    tl.cex = 1.2
  ) +
    scale_fill_manual(values = view.colors)
  
  mtext("Factor vs Sample Groups", side = 1, line = 4)
  
  # Close PNG device
  while (dev.cur() > 1)
    dev.off()
  
  print(paste("Saved combined factor correlation plot to", output_path))
}



########################################
# PER VIEW PLOTS
########################################

# Helper to ensure output directories exist
ensure_output_dirs <- function(factor_name, view_name) {
  output_dir <- file.path("Results", factor_name, view_name)
  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)
  return(output_dir)
}

# Save weight table per view
save_weight_table <- function(factor_name) {
  views <- names(MOFAobject.trained@expectations$W)
  for (view_name in views) {
    message(paste(
      "Saving weights for view:",
      view_name,
      "and factor:",
      factor_name
    ))
    output_dir <- ensure_output_dirs(factor_name, view_name)
    weight_df <- data.frame(
      Feature = gsub(
        paste0("_", view_name),
        "",
        rownames(MOFAobject.trained@expectations$W[[view_name]])
      ),
      Weight = MOFAobject.trained@expectations$W[[view_name]][, factor_name],
      AdjWeight = MOFAobject.trained@expectations$adjusted_W[[view_name]][, factor_name],
      MOFAobject.trained@expectations$Pvalues[[view_name]],
      MOFAobject.trained@expectations$FDRPvalues[[view_name]],
      MOFAobject.trained@data[[view_name]][["group1"]]
    )
    write.csv(weight_df,
              file.path(output_dir, "Weight_TPM.csv"),
              row.names = FALSE)
  }
}

# Volcano plots per view
plot_volcano_per_pval <- function(factor_name) {
  views <- names(MOFAobject.trained@expectations$W)
  for (view_name in views) {
    message(paste(
      "Generating volcano plot for view:",
      view_name,
      "and factor:",
      factor_name
    ))
    output_dir <- ensure_output_dirs(factor_name, view_name)
    pval_df <- data.frame(Weight = MOFAobject.trained@expectations$W[[view_name]][, factor_name],
                          MOFAobject.trained@expectations$Pvalues[[view_name]])
    pval_cols <- setdiff(colnames(pval_df), "Weight")
    plot_list <- list()
    for (pval_col in pval_cols) {
      pval_df$logPval <- -log10(pval_df[[pval_col]])
      pval_df$Color <- ifelse(
        pval_df[[pval_col]] < 0.05 & pval_df$Weight > 0,
        "red",
        ifelse(pval_df[[pval_col]] < 0.05 &
                 pval_df$Weight < 0, "blue", "grey")
      )
      p <- ggplot(pval_df, aes(x = Weight, y = logPval, color = Color)) +
        geom_point(size = 3) +
        scale_color_identity() +
        labs(title = pval_col, x = "Weight", y = "-log10(P-value)") +
        theme_minimal()
      plot_list[[pval_col]] <- p
    }
    ncol <- min(8, length(plot_list))
    nrow <- ceiling(length(plot_list) / ncol)
    png(
      file.path(output_dir, "Volcano_plots.png"),
      width = 14,
      height = 1.75 * nrow,
      units = "in",
      res = 300
    )
    grid.arrange(grobs = plot_list, ncol = ncol)
    while (dev.cur() > 1)
      dev.off()
  }
}

# Top weights bar plot per view
plot_top_weights_view <- function(factor_name) {
  views <- names(MOFAobject.trained@expectations$W)
  for (view_name in views) {
    message(
      paste(
        "Generating top weights plot for view:",
        view_name,
        "and factor:",
        factor_name
      )
    )
    output_dir <- ensure_output_dirs(factor_name, view_name)
    top_weights_plot <- plot_top_weights(
      MOFAobject.trained.cleanNames,
      view = view_name,
      factors = factor_name,
      nfeatures = 25,
      abs = TRUE,
      scale = TRUE,
      sign = "all"
    )
    ggsave(
      plot = top_weights_plot,
      filename = file.path(output_dir, "Top_Weights.png"),
      width = 15,
      height = 15
    )
  }
}

# Scatter plot per view
plot_scatter_per_view <- function(factor_name) {
  views <- names(MOFAobject.trained@expectations$W)
  for (view_name in views) {
    message(paste(
      "Generating scatter plot for view:",
      view_name,
      "and factor:",
      factor_name
    ))
    output_dir <- ensure_output_dirs(factor_name, view_name)
    scatter_plot <- plot_data_scatter(
      MOFAobject.trained,
      view = view_name,
      factor = factor_name,
      features = 25,
      dot_size = 3,
      color_by = "category",
      legend = TRUE
    ) +
      scale_fill_manual(values = group.colors) +
      theme(axis.text = element_text(size = 10))
    ggsave(
      plot = scatter_plot,
      filename = file.path(output_dir, "Scatter.png"),
      width = 30,
      height = 15
    )
  }
}

# Heatmap per view
plot_heatmap_per_view <- function(factor_name) {
  views <- names(MOFAobject.trained@expectations$W)
  for (view_name in views) {
    message(paste(
      "Generating heatmap for view:",
      view_name,
      "and factor:",
      factor_name
    ))
    output_dir <- ensure_output_dirs(factor_name, view_name)
    heatmap_plot <- plot_data_heatmap(
      MOFAobject.trained.cleanNames,
      factor = factor_name,
      view = view_name,
      denoise = TRUE,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_colnames = TRUE,
      show_rownames = TRUE,
      annotation_samples = "category",
      features = 25,
      annotation_colors = list("category" = group.colors),
      annotation_legend = FALSE,
      scale = "row"
    )
    ggsave(
      plot = heatmap_plot,
      filename = file.path(output_dir, "Heatmap.png"),
      width = 10,
      height = 10
    )
    tryCatch(while (dev.cur() > 1)
      dev.off())
  }
}

##### same as above but chatgpt ver #########
create_weight_outputs <- function(factor_name, MOFAobject.trained, view_groups) {
  # Ensure output directory exists
  output_dir <- paste0("Results/", factor_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Collect top features and weights across all views for this factor
  top_features_list <- list()
  all_weights_list <- list()
  
  for (view_name in names(get_weights(MOFAobject.trained))) {
    message(
      paste("Generating top features lists for view:", view_name, "and factor:", factor_name)
    )
    
    weight_df <- data.frame(
      Feature = gsub(paste('_', view_name, sep = ""), "", 
                     rownames(MOFAobject.trained.cleanNames@expectations$W[[view_name]])),
      View = view_name,
      Weight = MOFAobject.trained@expectations$W[[view_name]][, factor_name],
      AdjWeight = MOFAobject.trained@expectations$adjusted_W[[view_name]][, factor_name],
      MOFAobject.trained@expectations$Pvalues[[view_name]],
      MOFAobject.trained@expectations$FDRPvalues[[view_name]],
      MOFAobject.trained@data[[view_name]][["group1"]],
      MOFAobject.trained@expectations$logFC[[view_name]]
    )
    
    weight_df <- weight_df %>% arrange(desc(abs(Weight)))
    top_features_list[[view_name]] <- head(weight_df, 10)
    all_weights_list[[view_name]] <- weight_df
  }
  
  # Save all weights to CSV
  all_weights_df <- bind_rows(all_weights_list)
  write.csv(all_weights_df, paste0(output_dir, "/All_Weights.csv"), row.names = FALSE)
  
  # Group-level top features and plots
  for (group_name in names(view_groups)) {
    view_set <- view_groups[[group_name]]
    filtered_top_features <- top_features_list[view_set]
    top_features_df <- bind_rows(filtered_top_features)
    
    write.csv(
      top_features_df,
      file = paste0(output_dir, "/", group_name, "-Top_Features.csv"),
      row.names = FALSE
    )
    
    # Plot
    top_weights_df <- top_features_df %>% arrange(desc(abs(Weight)))
    top_weights_df$Feature <- substr(top_weights_df$Feature, 1, 40)
    
    top_weights_plot <- ggplot(top_weights_df, aes(
      x = reorder(Feature, -abs(Weight)),
      y = Weight,
      fill = View
    )) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(
        title = paste("Top Weights for", factor_name, "-", group_name),
        x = "Feature",
        y = "Weight"
      ) +
      theme_classic() +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
    
    ggsave(
      plot = top_weights_plot,
      filename = paste0(output_dir, "/", group_name, "-Top_Weights.png"),
      width = 10,
      height = 8
    )
  }
  
  return(invisible(NULL))
}


# Helper function to get GSEA input vectors for a given factor
get_gsea_input <- function(factor_name) {
  # Extract KEGG weights for the current factor across views
  kegg_weights <- MOFAobject.trained@expectations$KEGG_weights
  # Filter out any with less than 2 entries (if applicable)
  kegg_weights_filtered <- Filter(function(x)
    nrow(x) > 1, kegg_weights)
  
  # Extract named weight vectors for the specified factor from each available matrix/data frame
  extracted_list <- lapply(kegg_weights_filtered, function(df) {
    if ((is.matrix(df) ||
         is.data.frame(df)) && factor_name %in% colnames(df)) {
      vec <- df[, factor_name]
      names(vec) <- rownames(df)
      return(vec)
    } else {
      return(NULL)
    }
  })
  # Remove any NULL entries
  extracted_list <- extracted_list[!sapply(extracted_list, is.null)]
  
  # Deduplicate weights by taking mean for duplicated names, and sort in decreasing order
  dedup_vector <- function(x) {
    v <- tapply(x, names(x), mean)
    v <- as.vector(v)
    names(v) <- names(x <- tapply(x, names(x), mean))
    v[order(v, decreasing = TRUE)]
  }
  # Apply deduplication to all extracted vectors
  gsea_input <- lapply(extracted_list, dedup_vector)
  return(gsea_input)
}

# Function to perform KEGG enrichment analysis and output results for a given factor
perform_kegg_enrichment <- function(factor_name) {
  # Ensure output directory exists
  output_dir <- paste0("Results/", factor_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  message(paste("Generating FSEA for factor:", factor_name))
  
  # Prepare GSEA input data for this factor
  gsea_input <- get_gsea_input(factor_name)
  
  # Run GSEA for metabolomic data (if available)
  gsea_results <- list()
  
  if (length(metabolomic_sheets) >= 1) {
    subset_input <- gsea_input[metabolomic_sheets]
    names(subset_input) <- metabolomic_sheets
    met_res <- compareCluster(
      geneClusters = subset_input,
      fun = "GSEA",
      TERM2GENE = TERM2GENE_list[["metabolite"]],
      TERM2NAME = TERM2NAME,
      pvalueCutoff = 0.05,
      pAdjustMethod = "none",
      minGSSize = 10,
      maxGSSize = 500,
      nPermSimple = 10000,
      eps = 0
    )
    metabolite_gsea_results <- met_res
  }
  # Run GSEA for transcriptomic data (if available)
  if (length(transcriptomic_sheets) >= 1) {
    subset_input <- gsea_input[transcriptomic_sheets]
    names(subset_input) <- transcriptomic_sheets
    gsea_res <- compareCluster(
      geneClusters = subset_input,
      fun = "GSEA",
      TERM2GENE = TERM2GENE_list[["gene"]],
      TERM2NAME = TERM2NAME,
      pvalueCutoff = 0.05,
      pAdjustMethod = "none",
      minGSSize = 10,
      maxGSSize = 500,
      nPermSimple = 10000,
      eps = 0
    )
    transcript_gsea_results <- gsea_res
  }
  
  # Compile GSEA results if any are significant
  if (length(transcriptomic_sheets) >= 1) {
    if (exists("transcript_gsea_results") &&
        !is.null(transcript_gsea_results) &&
        nrow(transcript_gsea_results@compareClusterResult) > 0) {
      gsea_results[["gene"]] <- transcript_gsea_results
    }
  }
  if (length(metabolomic_sheets) >= 1) {
    if (exists("metabolite_gsea_results") &&
        !is.null(metabolite_gsea_results) &&
        nrow(metabolite_gsea_results@compareClusterResult) > 0) {
      gsea_results[["metabolite"]] <- metabolite_gsea_results
    }
  }
  
  
  if (length(gsea_results) > 0) {
    # Combine results into a single data frame and save
    enrichment_df <- do.call(rbind,
                             lapply(gsea_results, function(x)
                               x@compareClusterResult))
    write.csv(
      enrichment_df,
      paste0(output_dir, "/KEGG_Enrichment_Weighted.csv"),
      row.names = FALSE
    )
    # Plot enrichment results (dot plot of top pathways)
    options(
      enrichplot.colours = c("#ffe119", "#3a13eb"),
      enrichplot.color = "pvalue"
    )
    combined_result <- NULL
    if (exists("transcript_gsea_results") &&
        !is.null(transcript_gsea_results) &&
        exists("metabolite_gsea_results") &&
        !is.null(metabolite_gsea_results)) {
      combined_result <- transcript_gsea_results
      combined_result@compareClusterResult <- rbind(
        transcript_gsea_results@compareClusterResult,
        metabolite_gsea_results@compareClusterResult
      )
    } else if (exists("transcript_gsea_results") &&
               !is.null(transcript_gsea_results)) {
      combined_result <- transcript_gsea_results
    } else if (exists("metabolite_gsea_results") &&
               !is.null(metabolite_gsea_results)) {
      combined_result <- metabolite_gsea_results
    }
    
    KEGG_Enrichment_plot <- dotplot(
      combined_result,
      color = "pvalue",
      font.size = 8,
      showCategory = 5,
      title = paste(factor_name, ": Feature Set Enrichment Analysis"),
      split = ".sign"
    ) +
      facet_grid(. ~ .sign) +
      aes(color = enrichmentScore) +
      scale_color_viridis_c(name = "enrichmentScore") +
      guides(size = "none") +
      theme(legend.key.height = unit(0.8, "cm"))
    Sys.sleep(5)
    widthfactor <- (length(unique(enrichment_df$Cluster)) * 4) + 2
    ggsave(
      plot = KEGG_Enrichment_plot,
      filename = paste0(output_dir, "/KEGG_Enrichment.png"),
      width = widthfactor
    )
    
    # Extract the result slot
    enrichment_df <- combined_result@compareClusterResult
    
    # Define output file path
    output_file <- file.path(output_dir, "KEGG_Enrichment_Weighted.csv")
    
    # Write to CSV
    write.csv(enrichment_df, output_file, row.names = FALSE)
    
  } else {
    message("No significant GSEA results available.")
  }
  return(invisible(NULL))
}

# Function to create a network plot of enriched KEGG terms and features for a given factor
plot_enrichment_network <- function(factor_name) {
  # Ensure output directory exists
  output_dir <- paste0("Results/", factor_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  message(paste("Generating Igraph cnplot for factor:", factor_name))
  
  # Read enrichment results (requires perform_kegg_enrichment to have been run)
  enrichment_file <- paste0(output_dir, "/KEGG_Enrichment_Weighted.csv")
  if (!file.exists(enrichment_file)) {
    message("No enrichment results file found for factor ", factor_name)
    return(invisible(NULL))
  }
  enrichment_df <- read.csv(enrichment_file, stringsAsFactors = FALSE)
  
  # Create mapping from KEGG IDs to metabolite/gene names using available annotations
  kegg_to_name <- c()
  if (exists("metabolite_annotation") &&
      !is.null(metabolite_annotation)) {
    map <- metabolite_annotation[!is.na(metabolite_annotation$KEGG_ID), ]
    kegg_to_name <- c(kegg_to_name, setNames(map$Metabolite, gsub("^cpd:", "", map$KEGG_ID)))
  }
  if (exists("gene_annotation") && !is.null(gene_annotation)) {
    map <- gene_annotation[!is.na(gene_annotation$KEGG_ID), ]
    kegg_to_name <- c(kegg_to_name, setNames(map$gene, gsub("^\\w+:", "", map$KEGG_ID)))
  }
  # Replace KEGG IDs with names in the core_enrichment field
  replace_kegg_with_name <- function(KEGG_IDs, mapping) {
    ids <- unlist(strsplit(KEGG_IDs, "/", fixed = TRUE))
    names_vec <- mapping[ids]
    names_vec[is.na(names_vec)] <- ids[is.na(names_vec)]
    paste(names_vec, collapse = "/")
  }
  annotated_enrichment_df <- enrichment_df
  annotated_enrichment_df$core_enrichment <- sapply(enrichment_df$core_enrichment,
                                                    replace_kegg_with_name,
                                                    mapping = kegg_to_name)
  
  # Build an edge list linking enriched terms to each associated compound/gene
  edge_list <- list()
  for (i in seq_len(nrow(annotated_enrichment_df))) {
    cluster <- annotated_enrichment_df$Cluster[i]
    term <- annotated_enrichment_df$Description[i]
    compounds_str <- annotated_enrichment_df$core_enrichment[i]
    pvalue <- annotated_enrichment_df$pvalue[i]
    enrichmentScore <- annotated_enrichment_df$enrichmentScore[i]
    if (is.na(compounds_str) || compounds_str == "")
      next
    compounds <- strsplit(compounds_str, "/", fixed = TRUE)[[1]]
    edge_list[[i]] <- data.frame(
      Cluster = rep(cluster, length(compounds)),
      term = rep(term, length(compounds)),
      pvalue = rep(pvalue, length(compounds)),
      enrichmentScore = rep(enrichmentScore, length(compounds)),
      compound = compounds,
      stringsAsFactors = FALSE
    )
  }
  if (length(edge_list) == 0) {
    message("No enriched terms to visualize for factor ", factor_name)
    return(invisible(NULL))
  }
  edges_all <- do.call(rbind, edge_list)
  # Filter edges to only include significant pathways (p < 0.05)
  edges_all <- edges_all[edges_all$pvalue < 0.05, ]
  edges_all <- edges_all %>% mutate(
    compound_label = compound,
    compound_clustered = paste0(compound, " (", Cluster, ")")
  )
  # Select and limit to top 100 edges by absolute enrichment score
  # Select and arrange top 100 edges by enrichment score
  edges <- edges_all %>%
    dplyr::select(from = term,
                  to = compound_clustered,
                  enrichmentScore,
                  Cluster,
                  pvalue) %>%
    dplyr::arrange(desc(abs(enrichmentScore))) %>%
    head(100)
  
  if (nrow(edges) == 0) {
    message("No significant enriched connections to plot for factor ",
            factor_name)
    return(invisible(NULL))
  }
  
  # Create an igraph object and set attributes for nodes and edges
  g <- graph_from_data_frame(edges, directed = FALSE)
  terms <- unique(edges$from)
  V(g)$type <- ifelse(V(g)$name %in% terms, "term", "compound")
  # Aggregate attributes for compound nodes
  compound_attrs <- edges_all %>%
    dplyr::select(node = compound_clustered, label = compound_label, Cluster, pvalue) %>%
    dplyr::group_by(node) %>%
    dplyr::summarise(
      pvalue = mean(pvalue, na.rm = TRUE),
      Cluster = dplyr::first(as.character(Cluster)),
      label = dplyr::first(label),
      .groups = "drop"
    )
  # Initialize vertex attributes
  V(g)$pvalue <- NA_real_
  V(g)$Cluster <- NA_character_
  V(g)$label <- V(g)$name
  # Assign attributes to compound nodes
  compound_idx <- match(V(g)$name, compound_attrs$node)
  valid_idx <- !is.na(compound_idx)
  V(g)$pvalue[valid_idx] <- compound_attrs$pvalue[compound_idx[valid_idx]]
  V(g)$Cluster[valid_idx] <- compound_attrs$Cluster[compound_idx[valid_idx]]
  V(g)$label[valid_idx] <- compound_attrs$label[compound_idx[valid_idx]]
  # Set vertex size (term nodes get default size 10)
  V(g)$size <- -log10(V(g)$pvalue) * 3
  V(g)$size[is.na(V(g)$size)] <- 10
  V(g)$size[V(g)$size > 10] <- 10
  # Set vertex colors (terms: grey; compounds: view-specific color)
  V(g)$color <- ifelse(V(g)$type == "term", "grey95", view.colors[V(g)$Cluster])
  V(g)$color[is.na(V(g)$color)] <- "grey96"
  # Label size and color
  V(g)$label.cex <- ifelse(V(g)$type == "term", 0.6, 0.5)
  V(g)$label.color <- ifelse(V(g)$type == "term", "#e6194b", "black")
  # Edge width scaled by enrichment score
  E(g)$width <- rescale(abs(E(g)$enrichmentScore), to = c(1, 5))
  
  
  # Plot the network to a PNG file
  png(
    filename = file.path(output_dir, "Top_100_edges_KEGG_cnet.png"),
    width = 2000,
    height = 2000,
    res = 300
  )
  par(family = "Arial", mar = c(1, 1, 3, 1))
  set.seed(42)
  coords <- layout_with_fr(g)
  plot(
    g,
    layout = coords,
    vertex.label = V(g)$label,
    vertex.label.cex = V(g)$label.cex,
    vertex.label.color = V(g)$label.color,
    vertex.label.family = "Arial",
    vertex.label.dist = 0.5,
    edge.color = "grey80"
  )
  title(main = "CnetPlot of Top 100 Significantly Enriched (p < 0.05) Features", cex.main = 1)
  # Legend for data sources (views)
  legend(
    "bottomleft",
    legend = names(view.colors),
    col = view.colors,
    pch = 21,
    pt.bg = view.colors,
    pt.cex = 1,
    cex = 0.5,
    bty = "n",
    title = "Dataset View"
  )
  # Legend for p-value (node size)
  legend_pvals <- c(0.05, 0.01, 0.001)
  legend_sizes <- -log10(legend_pvals) * 3
  pt_cex_legend <- legend_sizes / 4
  legend(
    "topright",
    legend = c("    p = 0.05", "", "    p = 0.01", "", "    p = 0.001"),
    pt.cex = c(pt_cex_legend[1], NA, pt_cex_legend[2], NA, pt_cex_legend[3]),
    pch = 21,
    pt.bg = "white",
    col = "black",
    cex = 0.5,
    bty = "n",
    title = "Enrichment Pval"
  )
  # Legend for enrichment score (edge width)
  legend(
    "topleft",
    legend = c("Low enrichment", "High enrichment"),
    lwd = c(1, 5),
    col = "grey80",
    cex = 0.5,
    bty = "n",
    title = "Enrichment Score"
  )
  # Close the graphics device
  while (dev.cur() > 1)
    dev.off()
  return(invisible(NULL))
}

#generate pathway visualisation CHATGPT VER
generate_pathway_visualisations_all_views_for_factor <- function(factor_name) {
  message("Generating Pathview visualisations for factor: ", factor_name)
  
  # Load enrichment results
  enrichment_file <- file.path("Results", factor_name, "KEGG_Enrichment_Weighted.csv")
  if (!file.exists(enrichment_file)) {
    message("No enrichment file found for factor: ", factor_name)
    return(invisible(NULL))
  }
  
  enrichment_df <- read.csv(enrichment_file, stringsAsFactors = FALSE)
  enrichment_df <- enrichment_df[enrichment_df$pvalue < 0.05, ]
  if (nrow(enrichment_df) == 0) {
    message("No significant pathways for factor: ", factor_name)
    return(invisible(NULL))
  }
  
  enrichment_df$clean_id <- gsub("^map", "", enrichment_df$ID)
  enrichment_df$clean_id <- gsub(kegg_species, "", enrichment_df$clean_id)
  
  gsea_input <- get_gsea_input(factor_name)
  
  combine_gsea_named_list_to_matrix <- function(named_list) {
    df_long <- bind_rows(lapply(names(named_list), function(name) {
      vec <- named_list[[name]]
      if (length(vec) == 0) return(NULL)
      data.frame(sample = name, kegg_id = names(vec), value = as.numeric(vec))
    }))
    if (nrow(df_long) == 0) return(NULL)
    df_wide <- df_long %>% pivot_wider(names_from = sample, values_from = value)
    mat <- df_wide %>% column_to_rownames("kegg_id") %>% as.matrix()
    return(mat)
  }
  
  for (group_name in names(view_groups)) {
    message("🔎 Processing view group: ", group_name)
    views <- view_groups[[group_name]]
    
    gene_views <- intersect(views, transcriptomic_sheets)
    metabolite_views <- intersect(views, metabolomic_sheets)
    
    pathview_gene_data <- if (length(gene_views) >= 1) combine_gsea_named_list_to_matrix(gsea_input[gene_views]) else NULL
    pathview_cpd_data  <- if (length(metabolite_views) >= 1) combine_gsea_named_list_to_matrix(gsea_input[metabolite_views]) else NULL
    
    if (!is.null(pathview_gene_data)) pathview_gene_data[is.na(pathview_gene_data)] <- 0
    if (!is.null(pathview_cpd_data)) pathview_cpd_data[is.na(pathview_cpd_data)] <- 0
    
    output_dir <- file.path("Results", factor_name, paste0("FSEA-pathview_", group_name))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    enrichment_df_view <- enrichment_df[enrichment_df$Cluster %in% views, , drop = FALSE]
    
    for (i in seq_len(nrow(enrichment_df_view))) {
      raw_id <- enrichment_df_view$ID[i]
      clean_id <- enrichment_df_view$clean_id[i]
      message("Running Pathview for pathway: ", raw_id)
      
      pathview_call <- tryCatch({
        suppressMessages(suppressWarnings(
          pathview(
            gene.data   = pathview_gene_data,
            cpd.data    = pathview_cpd_data,
            pathway.id  = clean_id,
            species     = kegg_species,
            gene.idtype = "kegg",
            cpd.idtype  = "kegg",
            kegg.native = TRUE,
            multi.state = TRUE,
            match.data  = FALSE,
            same.layer  = FALSE,
            low         = list(gene = "#bcf60c", cpd = "#3a13eb"),
            mid         = list(gene = "grey", cpd = "grey"),
            high        = list(gene = "#e6194b", cpd = "#ffe119"),
            verbose     = FALSE
          )
        ))
      }, error = function(e) {
        message("pathview failed for ", raw_id, ": ", conditionMessage(e))
        return(NULL)
      })
      
      if (is.null(pathview_call)) next
      
      # Handle PNG output renaming if needed
      multi_png_path <- paste0(kegg_species, clean_id, ".pathview.multi.png")
      png_path <- paste0(kegg_species, clean_id, ".pathview.png")
      xml_path <- paste0(kegg_species, clean_id, ".xml")
      
      if (file.exists(multi_png_path)) file.rename(multi_png_path, png_path)
      
      if (!file.exists(png_path)) {
        message("PNG not generated for pathway: ", raw_id, " — skipping image step.")
        next
      }
      
      # Get pathway name from XML (fallback to raw_id)
      pathway_name <- tryCatch({
        ann <- xmlToList(xml_path)$.attrs[[4]]
        ann <- gsub("[^A-Za-z0-9_]", "_", ann)
        gsub("___+", "_", ann)
      }, error = function(e) raw_id)
      
      # Extract gene/compound data
      gene_ann <- tryCatch({
        if (!is.null(pathview_call$plot.data.gene)) {
          gene_data <- as.data.frame(pathview_call$plot.data.gene)[, c(1, 10)]
          names(gene_data)[1] <- "KEGG_ID"
          ann <- dplyr::left_join(gene_data, TERM2GENE_list$gene, by = "KEGG_ID")
          ann$pathway <- pathway_name
          ann$view <- "transcriptomic"
          ann
        } else NULL
      }, error = function(e) NULL)
      
      cpd_ann <- tryCatch({
        if (!is.null(pathview_call$plot.data.cpd)) {
          cpd_data <- as.data.frame(pathview_call$plot.data.cpd)[, c(1, 10)]
          names(cpd_data)[1] <- "KEGG_ID"
          ann <- dplyr::left_join(cpd_data, TERM2GENE_list$metabolite, by = "KEGG_ID")
          ann$pathway <- pathway_name
          ann$view <- "metabolomic"
          ann
        } else NULL
      }, error = function(e) NULL)
      
      all_ann <- dplyr::bind_rows(gene_ann, cpd_ann)
      write.table(all_ann, file = file.path(output_dir, paste0(pathway_name, ".tsv")),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      
      # Compose combined image
      summary_img_path <- "gsea_input_summary.png"
      png(summary_img_path, width = 800, height = 200)
      plot.new()
      title("GSEA Input Summary")
      if (!is.null(pathview_cpd_data)) {
        text(0.05, 0.9, "Compound sheet order:", adj = 0)
        text(0.05, 0.7, paste(colnames(pathview_cpd_data), collapse = " | "), adj = 0)
      }
      if (!is.null(pathview_gene_data)) {
        text(0.05, 0.3, "Gene sheet order:", adj = 0)
        text(0.05, 0.1, paste(colnames(pathview_gene_data), collapse = " | "), adj = 0)
      }
      dev.off()
      
      # Read and combine images
      tryCatch({
        pathview_img <- image_read(png_path)
        summary_img <- image_read(summary_img_path)
        summary_img <- image_resize(summary_img, geometry = paste0(image_info(pathview_img)$width, "x"))
        combined_img <- image_append(c(pathview_img, summary_img), stack = TRUE)
        image_write(combined_img, path = file.path(output_dir, paste0(pathway_name, ".png")))
      }, error = function(e) {
        message("Error combining image for ", pathway_name, ": ", conditionMessage(e))
      })
      
      # Clean up
      file.remove(list.files(".", pattern = "\\.(png|xml)$", full.names = TRUE))
    }
  }
}



run_fella_enrichment <- function(factor_name) {
  for (met_view in metabolomic_sheets) {
    message(paste("Generating FSEA for view: ", met_view, " and factor:" , factor_name))
    
    # Extract KEGG weights for this view
    kegg_weights <- MOFAobject.trained@expectations$KEGG_weights[[met_view]]
    
    if (is.null(kegg_weights) ||
        !factor_name %in% colnames(kegg_weights)) {
      message("Skipping ", met_view, ": no data or factor not found.")
      next
    }
    
    # Extract vector
    vec <- kegg_weights[, factor_name]
    names(vec) <- rownames(kegg_weights)
    
    # Deduplicate
    vec_dedup <- tapply(vec, names(vec), mean)
    vec_dedup <- vec_dedup[order(vec_dedup, decreasing = TRUE)]
    
    # Filter to valid KEGG compound IDs
    input_ids <- colnames(fella.data@diffusion@matrix)
    valid_ids <- names(vec_dedup)[names(vec_dedup) %in% input_ids]
    
    if (length(valid_ids) == 0) {
      message("No valid KEGG compound IDs in FELLA DB for ", met_view)
      next
    }
    
    # Run FELLA
    fella.user <- defineCompounds(compounds = valid_ids, data = fella.data)
    fella.user <- runDiffusion(object = fella.user, data = fella.data)
    
    # Generate table
    fella.table <- generateResultsTable(
      object = fella.user,
      data = fella.data,
      method = "diffusion",
      nlimit = 50
    )
    
    # Prepare output dir
    output_dir <- file.path("Results", factor_name, met_view, "FELLA")
    dir.create(output_dir,
               recursive = TRUE,
               showWarnings = FALSE)
    write.csv(fella.table,
              file.path(output_dir, "FELLA_Enrichment.csv"),
              row.names = FALSE)
    
    # Generate and subset enrichment graph
    graph <- generateResultsGraph(
      object = fella.user,
      method = "diffusion",
      data = fella.data,
      threshold = 0.05,
      nlimit = 1000
    )
    
    # Extract top nodes
    fella.table.full <- generateResultsTable(object = fella.user,
                                             data = fella.data,
                                             method = "diffusion")
    
    available_nodes <- V(graph)$name
    top_nodes <- head(fella.table.full$KEGG.id[fella.table.full$KEGG.id %in% available_nodes], 10)
    
    neighbours <- unlist(neighborhood(graph, order = 1, nodes = top_nodes))
    selected_nodes <- unique(c(top_nodes, V(graph)$name[neighbours]))
    graph_top10ctx <- induced_subgraph(graph, vids = selected_nodes)
    
    # Save graph plot
    png(
      file.path(
        output_dir,
        "FELLA_enrichment_graph_top_10_nodes_with_context.png"
      ),
      width = 2000,
      height = 2000,
      res = 300
    )
    plotGraph(graph_top10ctx, data = fella.data)
    dev.off()
    
    # Save input compound list
    sig_compounds <- getInput(fella.user)
    writeLines(sig_compounds,
               file.path(output_dir, "Input_Compounds.txt"))
  }
}


###########################################################
###########################################################
# CODE EXECUTION
###########################################################
###########################################################

# ==== Initialise the environment ====

# Sets up project environment, optionally with renv support
setup_project_environment(use_renv = FALSE)


# ==== Load and pre-process input data ====

# Load raw Excel sheets, clean up 'Name' columns, drop unused columns,
# harmonise column ordering across views
data_list <- load_and_clean_input_data(input_file, sheet_names, new_names)

# ==== Lipid annotation ====

# Load existing or generate new lipid annotations from lipidomic sheet names
#LIPID_ANNOTATIONS <- generate_lipid_annotations(lipidomic_sheets, data_list)


# ==== Generate KEGG pathway mappings ====

# Downloads or loads TERM2GENE and TERM2NAME for gene and metabolite KEGG pathways
kegg_results <- generate_kegg_annotation_lists(kegg_species)
TERM2GENE_list <- kegg_results$TERM2GENE
TERM2NAME <- kegg_results$TERM2NAME


# ==== Annotate transcriptomic data ====

# Maps transcriptomic features (e.g. Ensembl IDs) to KEGG gene IDs via Ensembl → Entrez conversion
gene_annotation <- annotate_transcriptomic_data_with_kegg(
  transcriptomic_sheets = transcriptomic_sheets,
  data_list = data_list,
  TERM2GENE_list = TERM2GENE_list,
  ensembl_to_entrez_mapping = ensembl_to_entrez_mapping,
  # data.frame with ensembl_gene_id and entrezgene_id
  kegg_species = "taes"
)

# ==== Annotate metabolomic data ====

# Matches compound names (post-pipe cleaned) to KEGG compound IDs using keggFind()
metabolite_annotation <- annotate_metabolomic_data_with_kegg(metabolomic_sheets = metabolomic_sheets, data_list = data_list)

# ==== Build or load FELLA database ====
# Builds or loads FELLA graph and diffusion matrix for pathway enrichment
fella.data <- load_or_build_fella_database(kegg_species)


# ==== Data normalisation and MOFA model training ====

# Normalise (impute, log-transform, zero out infs) and re-structure data for MOFA scale methods: "zscore", "minmax", "quantile"
mofa_list <- normalise_data_mofa(data_list, na_limit = 10-6)

# Create and train a new MOFA model, or load existing hdf5
MOFAobject.trained <- initialise_and_run_mofa(mofa_list)

# ==== Post-process MOFA model with KEGG annotations and per-feature analysis ====

MOFAobject.trained <- postprocess_mofa_model(
  MOFAobject.trained,
  kegg_species,
  metabolite_annotation = metabolite_annotation,
  gene_annotation = gene_annotation,
  transcriptomic_sheets = transcriptomic_sheets,
  metabolomic_sheets = metabolomic_sheets
)


# ==== Sample metadata and colour setup for visualisation ====

# Create metadata from sample names (e.g. categories), binary indicators, groupings
metadata <- generate_sample_metadata(colnames(mofa_list[[1]]))

# Attach sample metadata to MOFA object
MOFAobject.trained <- assign_sample_metadata_to_mofa(MOFAobject.trained, metadata)

# Define colour palettes per group and per view
group.colors <- define_group_colours(metadata)

view.colors <- define_view_colours(names(MOFAobject.trained@data))

# Create a MOFA object with standardised feature names (for plotting)
MOFAobject.trained.cleanNames <- standardise_feature_names(MOFAobject.trained)


# ==== Visualisation of MOFA global outputs ====

# Plot proportion of variance explained
generate_mofa_variance_plots(MOFAobject.trained, view.colors)

# Plot factor matrices (samples x factors)
generate_mofa_factor_matrix_plot(MOFAobject.trained, view.colors, group.colors)

# Plot comparisons of factor values across sample groups
generate_mofa_factor_comparison_plot(MOFAobject.trained, group.colors, view.colors)

# Plot factor–factor and factor–covariate correlation matrices
generate_mofa_factor_correlation_plot(MOFAobject.trained, view.colors)


# ==== Per-factor MOFA output generation ====

# Extract list of latent factor names from MOFA
#factors <- colnames(get_factors(MOFAobject.trained)[["group1"]])

factors <- colnames(get_factors(MOFAobject.trained, groups = "group1")[[1]]) #chatgpt ver

view_groups <- list(
  Transcriptomics = c("RNAseq"),
  Metabolomics = c("LCMS")
)

# Loop over factors to generate all downstream outputs
for (factor_name in factors) {
  message(paste("################################"))
  message(paste("Processing factor:", factor_name))
  message(paste("################################"))
  message("")
  
  # Generate and export weight tables for each view
  save_weight_table(factor_name)
  
  # Generate per-view volcano plots (logFC vs -log10 p-value)
  plot_volcano_per_pval(factor_name)
  
  # Plot top-ranked features per view
  plot_top_weights_view(factor_name)
  
  # Plot scatter plots of top feature weights across views
  plot_scatter_per_view(factor_name)
  
  # Plot heatmap of top-ranked features
  plot_heatmap_per_view(factor_name)
  
  # Save all weights in tidy format
  create_weight_outputs(factor_name, MOFAobject.trained, view_groups)
  
  # Perform KEGG pathway enrichment analysis
 perform_kegg_enrichment(factor_name)
  
  # Plot enrichment results as network
  plot_enrichment_network(factor_name)
  
  # Generate Pathview diagrams for enriched pathways
  generate_pathway_visualisations_all_views_for_factor(factor_name)
  
  # Perform FELLA-based enrichment analysis
  run_fella_enrichment(factor_name)
}
