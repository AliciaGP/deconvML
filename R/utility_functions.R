#' List data sets included
#' @description getDataSetList shows data sets included as they are loom class variables.
#' @return vector of strings with the data sets included.
#' @example /home/aligo/deconvML/vignettes

getDataSetList <- function() {

  classVariables <-  eapply(.GlobalEnv, class)
  index <-  grep("loom", classVariables)
  datasets <- names(classVariables[index])

  if (length(datasets) != 0) {
    return(datasets)
  } else {
    return ("No dataset loaded")
  }

}



#' List main properties describing a data set
#' @description getDataSetInfo displays the following characteristics of a data
#'   set: version, number of cells, number of genes, number of type of tissues
#'   and number of type of cells.
#' @param dataset String with the name of the data set.
#' @return data frame of properties describing the dataset identified by the
#'   string.
#' @example /home/aligo/deconvML/vignettes

getDataSetInfo <- function(dataset) {

  properties <- as.data.frame(cbind(zheisel$version,
                                    zheisel$matrix$dims[1],
                                    zheisel$matrix$dims[2],
                                    nrow(table(zheisel$col.attrs$Tissue[])),
                                    nrow(table(zheisel$col.attrs$Class[]))))

  rownames(properties) <- paste0(deparse(substitute(zheisel)), "_", "Info")

  colnames(properties) <- c("Version", "NumberCells", "NumberGenes",
                            "NumberTypeTissues", "NumberTypeCells")

  return(properties)
}

#' List cell types included in a specific data set
#' @description getCellTypesInDataSet shows the different cell types for a
#'   specific data set.
#' @inheritParams getDataSetInfo
#' @return cell types considered in the data set identified by the string.
#' @example /home/aligo/deconvML/vignettes

getCellTypesInDataSet <- function(dataset) {

  if (!is.na(match("Class", names(dataset$col.attrs)))) {

    cellTypes <- names(table(dataset$col.attrs$Class[]))

    } else {

      cellTypes <- NULL
  }

  if (length(cellTypes)!=0) {
    return(cellTypes)

  } else {

    return("No cell types found")
  }
}



#' List tissue types included in a specific data set
#' @description getTissueTypesInDataSet collects all types of tissues included
#'   in the data set.
#' @inheritParams getDataSetInfo
#' @return the tissue types considered in the data set identified by the string.
#' @example /home/aligo/deconvML/vignettes

getTissueTypesInDataSet <- function(dataset) {

  if (!is.na(match("Tissue", names(dataset$col.attrs)))) {
    tissueTypes <-  names(table(dataset$col.attrs$Tissue[]))

  } else {

    tissueTypes <- NULL
  }

  if (length(tissueTypes)!=0) {

    return(tissueTypes)

  } else {


  }
  return("No tissue types found")
}



#' List available methods to estimate the active genes in a specific data set
#' @description getActiveGeneMethod returns the methods available to apply on a
#'   specific data set.
#' @inheritParams getDataSetInfo
#' @return available methods data.
#' @example /home/aligo/deconvML/vignettes

getActiveGeneMethod <- function(dataset) {

  if (!is.na(match("matrix", dataset$names))) {
    methods <-  c("VariationCoefficient", "RawExpression")

    } else {

      methods <-  NULL
      }
  if (length(methods)!=0) {
    return(methods)
  } else {
    return ("No methods found")
  }
}



#' List active genes in the tissue, cell type, of a specific data set for a
#' given method
#' @description getActiveGenesInCellTypes read the files containing the
#'   pre-calculated matrices to return ids of active genes in the tissue, cell
#'   type of a specific data set for a given method.
#' @inheritParams getDataSetInfo
#' @param region String with the name of the region in which we are interested.
#' @param celltype String with the name of the celltype in which we are
#'   interested.
#' @param method String with the name of the method in which we are interested.
#' @return gene ids of genes active in the region, cell type, for a given
#'   calculus method.
#' @example /home/aligo/deconvML/vignettes

getActiveGenesInCellType <- function (dataset, tissue, celltype, method){

  expr.data <- readRDS(paste0(dataset, "_", tissue, "_", method, ".rds"))

  expr.data <- expr.data[rownames(expr.data) == celltype, ]
  IDgenes  <- colnames(expr.data)[expr.data == 1]

  return(IDgenes)
}
