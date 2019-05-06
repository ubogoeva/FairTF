#' Annotaion for terms taking into account True Path Rule
#'
#' @param rellist - list contains child terms as keys and nearest parents as values
#' @param annotlist - list contains terms as keys and direct annotations as values
#'
#' @return - list contains terms as keys and genes annotated according to True Path Rule as values
#' @export
#'
#' @examples
CreateAnnotationWithAncestry <- function(rellist, annotlist) {
  # extract terms that are leaves and roots in GO DAG to corresponding vectors
  leaves <- FindBoundaries(rellist)$leaves
  roots <- FindBoundaries(rellist)$roots
  
  # delete terms that are not presented in relation file (it must be a list entry)
  terms_not_in_rel <-  setdiff(names(annotlist), names(rellist))
  annotlist <- annotlist[-match(setdiff(terms_not_in_rel, roots), names(annotlist))]
  
  # traversing graph in down-up manner starting from leaves to the roots
  while (length(leaves) != 0) {
    #temporary vector for storing parental terms for next iteration of loop
    tempvec <- c()
    for (i in 1:length(leaves)) {
      if (!leaves[i] %in% roots) {
        # put all nearest parents for current term in parents vector
        # then iterate through parents vector and assign them child term annotation
        parents <- rellist[[leaves[i]]]
        for (j in 1:length(parents)) {
          annotlist[[parents[j]]] <-
            unique(c(annotlist[[parents[j]]], annotlist[[leaves[i]]]))
        }
        
        # append parents of current term to temporary vector
        tempvec <- unique(c(tempvec, parents))
      }
    }
    
    # reassign leaves vector with parents for next iteration of loop
    leaves <- tempvec
  }
  return(annotlist)
}

#' Intersection of annotation with custom background
#'
#' @param annotlist - list contains terms as keys and direct annotations as values
#' @param background_genes - vector contains background set of genes
#'
#' @return - list contains terms as keys and direct annotations as values with genes only from background
#' @export
#'
#' @examples
IntersectAnnotationWithBackground <- function(annotlist, background_genes) {
  term_ids <- names(annotlist)
  
  # intersect annotations for each term with background
  annotlist <-
    lapply(annotlist, function(x)
      intersect(x, background_genes))
  
  names(annotlist) <- term_ids
  
  # delete empty entries
  annotlist <- annotlist[lapply(annotlist, length) > 0]
  return(annotlist)
}

#' DAG roots and leaves extractor
#'
#' @param rellist - list contains child terms as keys and nearest parents as values
#'
#' @return - list contains root and leaf terms
#' @export
#'
#' @examples
FindBoundaries <- function(rellist) {
  outlist <- list()
  
  # extract child terms and parental terms
  children <- names(rellist)
  parents <- unlist(rellist)
  
  # define terms that are not  presented in the parental set as leaves
  outlist$leaves <- setdiff(children, parents)
  
  # define terms that are not presented in the children set as roots
  outlist$roots <- setdiff(parents, children)
  
  return(outlist)
}

#' Functional annotation test
#'
#' @param annotlist - list contains terms as keys and genes annotated according to True Path Rule as values
#' @param inputgeneset - vector of query genes
#' @param realnames - list contains terms as keys and corresponding real names as values
#' @param gene_amount_border - minimal number of genes annotated to a term (1 by default)
#' @param fisher.alternative - fisher exact test alternative,
#'                             possible values: "greater", "less", "two.sided" ("greater by default)
#' @param namespaces - list contains terms as keys and corresponding namespaces as values
#' @param p.adjust.method - method for multiple testing correction (to see all possible methods print: p.adjust.methods)
#'
#' @return - dataframe contains functional annotaion for gene set 
#' @export
#'
#' @examples
FuncAnnotTest <- function(annotlist,
                          inputgeneset,
                          realnames,
                          namespaces,
                          gene_amount_border = 1,
                          fisher.alternative = "greater",
                          p.adjust.method = "BH") {
  # extract background genes set from annotation list and calculate the amount
  bggenes <- unique(unlist(annotlist))
  bgtot <- length(bggenes)
  
  # delete genes from query genes set that are not presented in the bcakground and calculate the amount
  inputgeneset <- intersect(inputgeneset, bggenes)
  qtot <- length(inputgeneset)
  
  annotlist_names <- names(annotlist)
  
  qit_list <-
    lapply(annotlist, function(x)
      intersect(inputgeneset, x))
  
  names(qit_list) <- annotlist_names
  
  qit_ints <- lengths(qit_list)
  
  qit_list <- qit_list[qit_ints >= gene_amount_border]
  
  qit_ints <- lengths(qit_list)
  
  annotlist <- annotlist[match(names(qit_list), names(annotlist))]
  
  annotlist_names <- names(annotlist)
  
  bg_ints <- lengths(annotlist)
  
  values_matrix <-
    cbind(qit_ints, rep(qtot, length(qit_ints)), bg_ints, rep(bgtot, length(qit_ints)))
  
  fisher_matrix <-
    cbind(
      values_matrix[, 1],
      values_matrix[, 2] - values_matrix[, 1],
      values_matrix[, 3] - values_matrix[, 1],
      values_matrix[, 4] - values_matrix[, 2] - values_matrix[, 3] + values_matrix[, 1]
    )
  
  #pvalues <- c(double(0))
  #for(i in 1:nrow(fisher_matrix)){
  # pvalues <- c(pvalues, fisher.test(matrix(fisher_matrix[i,], nrow = 2), fisher.alternative)$p.value)
  #}
  
  pvalues <-
    apply(fisher_matrix, 1, function(x)
      fisher.test(matrix(x, nrow = 2), fisher.alternative)$p.value)
  
  real_names <-
    unlist(realnames[match(names(annotlist), names(realnames))])
  
  name_spaces <-
    unlist(namespaces[match(names(annotlist), names(namespaces))])
  
  genes <-
    unlist(lapply(qit_list, function(x)
      paste(x, collapse = "/")))
  
  padj <- p.adjust(pvalues, method = p.adjust.method)
  outdf <- data.frame(
    "GO_id" = annotlist_names,
    "namespace" = name_spaces,
    "name" = real_names,
    "qit" = values_matrix[, 1],
    "qtot" = values_matrix[, 2],
    "bgit" = values_matrix[, 3],
    "bgtot" = values_matrix[, 4],
    "pval" = pvalues,
    "padj" = padj,
    "qit_genes" = genes,
    stringsAsFactors = FALSE
  )
  
  # sort output df by p-values in ascending order
  outdf <- outdf[order(outdf$pval), ]
  # add row numbering to output data frame
  rownames(outdf) <- seq(nrow(outdf))
  return(outdf)
}


#' Functional annotation for multiple gene sets
#'
#' @param annotlist - list contains terms as keys and genes annotated according to True Path Rule as values
#' @param realnames - list contains terms as keys and corresponding real names as values
#' @param dir - directory contains only files with genes ids for annotation
#' @param namespaces - list contains terms as keys and corresponding namespaces as values
#' @param gene_amount_border - minimal number of genes annotated to a term (1 by default)
#' @param p.adjust.method -  method for multiple testing correction (to see all possible methods print: p.adjust.methods)
#'                          Benjamini-Hochberg by default      
#'
#' @return - list with filenames as keys and annotaton data frames as values
#' @export
#'
#' @examples
MultipleSetsAnnot <-
  function(annotlist,
           realnames,
           namespaces,
           dir,
           gene_amount_border = 1,
           p.adjust.method = "BH") {
    files <- list.files(dir)
    
    outlist <- lapply(files, function(file) {
      inputgenes <-
        read.table(paste(dir, file, sep = "/"), stringsAsFactors = FALSE)
      inputgenes <- unlist(inputgenes)
      
      # run functonal annotation
      return(
        FuncAnnotTest(
          annotlist,
          inputgenes,
          realnames,
          namespaces,
          gene_amount_border = gene_amount_border,
          p.adjust.method = p.adjust.method
        )
      )
      
    })
    
    files <- gsub(".txt", "", files)
    names(outlist) <- files
    return(outlist)
  }

#' Functional annotation for multiple gene sets. 
#'
#' @param annotlist - list contains terms as keys and genes annotated according to True Path Rule as values
#' @param realnames - list contains terms as keys and corresponding real names as values
#' @param namespaces - directory contains only files with genes ids for annotation
#' @param dir - directory contains only files with genes ids for annotation
#' @param gene_amount_border - minimal number of genes annotated to a term (1 by default)
#' @param p.adjust.method - method for multiple testing correction (to see all possible methods print: p.adjust.methods)
#'                          Benjamini-Hochberg by default
#' @param cl - cluster object
#'
#' @return
#' @export
#'
#' @examples
ParMultipleSetsAnnot <-
  function(annotlist,
           realnames,
           namespaces,
           dir,
           gene_amount_border = 1,
           p.adjust.method = "BH",
           cl) {
    files <- list.files(dir)
    
    clusterExport(
      cl = cl,
      varlist = c(
        "run.func.annot.test",
        "annotlist",
        "realnames",
        "namespaces",
        "dir",
        "gene_amount_border",
        "p.adjust.method"
      ),
      envir = environment()
    )
    
    outlist <- parLapply(cl, files, function(file) {
      inputgenes <-
        read.table(paste(dir, file, sep = "/"), stringsAsFactors = FALSE)
      inputgenes <- unlist(inputgenes)
      
      # run functonal annotation
      return(
        FuncAnnotTest(
          annotlist,
          inputgenes,
          realnames,
          namespaces,
          gene_amount_border = gene_amount_border,
          p.adjust.method = p.adjust.method
        )
      )
    })
    files <- gsub(".txt", "", files)
    names(outlist) <- files
    
    return(outlist)
}


#' Convert GAF format type annotation to list contains GO term id's as keys and Gene ID's as values
#'
#' @param data - GAF annotation table (output of read.gaf function)
#' @id_pattern - regex for Gene identifier
#'
#' @return - list with GO term id's as keys and Gene ID's as values
#' @export
#'
#' @examples
convert.annotation.wide <-
  function(data, id_pattern = patterns_enum$Arabidopsis_thaliana) {
    # extract GO id's and Gene id's annotated to them
    annotdf <-
      data.frame(
        "GOID" = data$GO_ID,
        "GeneID" = data$DB_Object_Name,
        stringsAsFactors = FALSE
      )
    
    # leave Gene ID's only of appropriate format if id_pattern is not NULL
    if (!is.null(id_pattern)) {
      annotdf <- annotdf[grep(id_pattern, annotdf[, 2]), ]
    }
    # convert data frame to list with GO term id's as keys and corresponding gene id's as values
    annotdf[["GOID"]] <- as.factor(annotdf[["GOID"]])
    outlist <- split(annotdf, annotdf[["GOID"]])
    outlist <- lapply(outlist, function(x)
      x[, -1])
    return(outlist)
  }

patterns_enum <- list(Arabidopsis_thaliana = "AT[1-5,M,C]G\\d{5}")



#' GAF file reader
#'
#' @param gafpath - full path to GAF file
#'
#' @return - list contains information about file and data accessible by $info and $data keys correspondingly
#' @export
#'
#' @examples
read.gaf2 <- function(gafpath){
  lines <- readLines(gafpath)
  
  # separate header and data parts of file
  inds <- grep("!.*",lines)
  header <- lines[inds]
  lines <- lines[-inds]
  
  # split lines and put them in data frame with specified col names
  df <- matrix(character(0),ncol = 17, nrow = length(lines))
  df <- t(sapply(lines, function(x) strsplit(x, "\t")[[1]]))
  df <- as.data.frame(df,stringsAsFactors = FALSE)
  rownames(df) <- NULL
  colnames(df) <- c("DB","DB_Object_ID","DB_Object_Symbol",
                    "Qualifier","GO_ID","DB:Reference",
                    "Evidence_Code","With_(or)_From","Aspect",
                    "DB_Object_Name","DB_Object_Synonym","DB_Object_Type",
                    "Taxon","Date","Assigned_By",
                    "Annotation_Extension","Gene_Product_Form_ID")
  outlist <- list(info = header, data = df)
  return(outlist)
}

#' OBO file reader
#'
#' @param obopath - path to the .obo format file (open biomedical ontologies)
#'
#' @return - lists that contain GO terms ancestry data, names for GO terms, and corresponding namespace for each GO term
#' @export
#'
#' @examples
read.obo2 <- function(obopath) {
  
  # read all lines from file
  lines <- readLines(obopath)
  
  # empty lists for data and id vector for storing GO ids
  outdata <- list()
  datalist <- list()
  nameslist <- list()
  namespacelist <- list()
  id_vector <- c(character(0))
  
  # counters: i - for iterating trough lines, ind - for creating indices for list entries
  i = 1
  list_ind = 1
  
  while (i != length(lines)) {
    
    id <- ""
    name <- ""
    namespace <- ""
    
    # flag that assigned with "true" value if current GO term marked as obsolete
    is_obsolete <- "false"
    
    # vector for storing parental GO terms
    parents <- c(character(0))
    
    # process piece of lines that contain informaton about distinct GO term
    if (lines[i] == "[Term]") {
      i = i + 1
      
      # extract GO term id
      if (startsWith(lines[i], "id:")) {
        id <- lines[i]
        id <- gsub("id:", "", id)
        id <- gsub("^\\s+|\\s+$", "", id)
        i = i + 1
        
      }
      
      # extract GO term real name
      if (startsWith(lines[i], "name:")) {
        name <- gsub("name:", "", lines[i])
        name <- gsub("^\\s+|\\s+$", "", name)
        i = i + 1
      }
      
      # extract GO term corresponding namespace
      if (startsWith(lines[i], "namespace:")) {
        namespace <- gsub("namespace:", "", lines[i])
        namespace <- gsub("^\\s+|\\s+$", "", namespace)
        i = i + 1
      }
      
      # process the rest part of lines that refer to current GO term and extract parental GO term ids 
      # by "is_a" relation keyword or "true" value for is_obsolete flag 
      while (lines[i] != "") {
        if (startsWith(lines[i], "is_a:")) {
          parent <- gsub("is_a:|!.*", "", lines[i])
          parent <- gsub("^\\s+|\\s+$", "", parent)
          parents <- c(parents, parent)
        } else if (startsWith(lines[i], "is_obsolete:")) {
          is_obsolete <- gsub("is_obsolete:", "", lines[i])
          is_obsolete <- gsub("^\\s+|\\s+$", "", is_obsolete)
        }
        
        i = i + 1
      }
      
      # if GO term is not obsolete fill in the lists with the appropriate values
      if (is_obsolete == "false") {
        id_vector[[list_ind]] <- id
        datalist[[list_ind]] <- parents
        nameslist[[list_ind]] <- name
        namespacelist[[list_ind]] <- namespace
        list_ind = list_ind + 1
      }
    }
    
    i = i + 1
    
    
  }
  
  # assign GO ids as names to lists
  names(datalist) <- id_vector
  names(nameslist) <- id_vector
  names(namespacelist) <- id_vector
  
  # remove entries without values from relation data list
  datalist <- datalist[lapply(datalist, length) > 0]
  
  outdata$relation_data <- datalist
  outdata$real_names <- nameslist
  outdata$namespaces <- namespacelist
  return(outdata)
}

prepare_annotation <- function(onto, annot, bg){
  annot <- convert.annotation.wide(annot)
  annot <- IntersectAnnotationWithBackground(annot, bg)
  annot <- CreateAnnotationWithAncestry(onto$relation_data, annot)
  return(annot)
}
