#!/usr/bin/env Rscript
#@author : Thomas NEFF
#@description : Wrapper for edgeR bulk RNAseq analysis

# lib ---------------------------------------------------------------------
library(edgeR)
library(xlsx)

# fun ---------------------------------------------------------------------

#' Create DGEList object for DE expression analysis with edgeR
#' 
#' @param counts reads counts file
#' @param metadata metadata file associated with reads counts file
#' @param norm normalization method apply on counts
#' @param prior prior apply with counts normalization
#' @param ... metadata vars file added to DGEList object
#' 
#' @return DGEList object
#' 
#' @example y <- list_dge(counts = counts, metadata = metadata, x1 = as.factor(x1) , x2 = as.factor(x2))
dge_list <- function(counts,metadata,norm="TMM",prior=NULL,...){
  
  y <- edgeR::DGEList(counts = counts)
  y <- y[edgeR::filterByExpr(y, keep.lib.sizes=FALSE),]
  
  if (all(sapply(list(...), function(x) is.factor(x)))) {
    y$samples <- list(...)
  } else {
    stop("the class attribut added to y$samples must be factor")
  }
  
  if (!is.null(prior)) {
    y <- edgeR::calcNormFactors(y,norm=norm,prior=prior)
  } else {
    y <- edgeR::calcNormFactors(y,norm=norm)
  }
  
  return(y)
}

#' Create design matrix
#' 
#' @param y DGEList
#' @param formula fornula to design model matrix
#' @param var_interest varaible of interest
#'
#' @return design model matrix
#' 
#' @example 
#' Glial Neuronal . . . pairX
#' 1      0       . . .      0
#' 2      1       . . .      0
#' .      .       . . .      .
#' 10     1       . . .      1
dge_design <- function(y,formula,var_interest=NULL){
  
  design <- model.matrix(eval(formula), data = y$samples)
  if (!is.null(var_interest)){
    colnames(design) <- gsub(var_interest,"",colnames(design))
  }
  return(design)
}

#' Differential expresion test 
#' 
#' @param y DGEList object
#' @param design design model matrix
#' @param contrast_by_design design model matrix for main var tested in DE see edgeR::makeContrasts
#' @param adj_method stat method to adjust p-value
#' @param model see edgeR::glmQLFit et edgeR::glmLRT
#' 
#' @return DGEList with new attributs (fit, lrt or qlf, and tt) 
dge_DE <- function(y,design,contrast_by_design,adj_method="BH",model='LRT'){
  
  #test attribut for model argument
  stopifnot((model!='LRT') | (model!='LTR'))
  
  y <- edgeR::estimateDisp(y = y, design = design,robust=TRUE)
  
  if (model=='LRT') {
    #regression
    y$fit <- edgeR::glmFit(y = y, design = design)
    #test
    y$lrt <- edgeR::glmLRT(glmfit = y$fit,contrast = contrast_by_design)
    #results
    y$tt <- edgeR::topTags(y$lrt,n=Inf,adjust.method = adj_method)$table
  } 
  
  if (model=='QLF') {
    #regression
    y$fit <- edgeR::glmQLFit(y = y, design = design)
    #test
    y$qlf <- edgeR::glmQLFTest(glmfit = y$fit,contrast = contrast_by_design)
    #result
    y$tt <- edgeR::topTags(y$qlf,n=Inf,adjust.method = adj_method)$table
  }
  
  return(y)
}

#' Save the results
#' 
#' @param y_lrt DGEList after dge_DE execution (model=LRT)
#' @param y_qlf DGEList after dge_DE execution (model=QLF)
#' @param file filename or file pathway
#' 
#' @return file xlsx with result
dge_xlsx <- function(y_lrt=NULL,y_qlf=NULL,file){

  #if file with same name already exist in folder, it will be deleted
  if (file.exists(file)) {
    system(paste0('rm ',file))
  }
  
  if ((!is.null(y_lrt)) & (!is.null(y_qlf))) {
    write.xlsx(x = y_lrt$tt, sheetName = "LRT",file = file)
    write.xlsx(x = y_qlf$tt, sheetName = "QLF", file = file, append = TRUE)
  } else {
    if (!is.null(y_lrt)) {
      write.xlsx(x = y_lrt$tt, sheetName = "LRT",file = file)
    }
    if (!is.null(y_qlf)) {
      write.xlsx(x = y_lrt$tt, sheetName = "QLF",file = file)
    }
  }
}