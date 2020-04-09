1. Introduction
---------------

The advent of single-cell genomics represents a turning point in cell
biology. For the first time, we can assay the expression level of every
gene in the genome across thousands of individual cells in a single
experiment (Trapnell [2015](#ref-trapnell2015defining)). This type of
study leads us to ask ourselves questions such as: what cell types does
your dataset have? how big are they? what are the most active genes in
cell type X?

To resolve these and many other questions, the deconvML package has been
developed. With this package we will be able to know the most relevant
information of our datasets, ask different questions that concern the
different cell types and tissues that compose it, know the methods that
we can apply to our dataset to estimate the active genes and,
ultimately, ask which ones are those genes for a certain cell type in a
tissue.

2. Materials
------------

The operation of this package is focused on datasets with a loom-type
structure, since these have the main advantage that loom objects are
merely connections to a file on disk, which enables scaling to massive
datasets with low memory consumption. The loom file is simply an HDF5
file with a strict structure imposed on it, consists of a container for
six sub-objects: one dataset (matrix) and five groups (layers,
row\_attrs, col\_attrs, row\_graphs, and col\_graphs). For more
information on this type of object, see:
<https://satijalab.org/loomR/loomR_tutorial.html>.

In this documentation, we will use as data set with this type of
structure, Zheisel, a mouse brain atlas of cell types from the
Linnarsson Lab (data available at:
<http://mousebrain.org/downloads.html>. For a first contact with this
repository and, in general with single-cell data, the documentation
“Processing raw data from Zheisel to get cell type specificity for
genes” by Juan A. Botía is recommended, where the data is examined of
single-cell focused on a single brain region.

3. Methodology
--------------

First of all, we have to install the deconvML package and load the
Zheisel dataset.

``` r
devtools::install_github('AliciaGP/deconvML')
#> Skipping install of 'deconvML' from a github remote, the SHA1 (145fdb84) has not changed since last install.
#>   Use `force = TRUE` to force installation
library(deconvML)
```

``` r
library("loomR")
#> Loading required package: R6
#> Loading required package: hdf5r
zheisel <- connect(filename = "/home/aligo/zheisel/l5_all.loom", mode = "r+") 
```

To show how to use the package we are going to show a typical script to
illustrate the process. In this case, the starting question for our
dataset will be: is the Gfap gene expressed in amygdala astrocytes?

To do this, the first thing we have to do is check the different loaded
datasets using the **getDataSetList** function. In this way, we ensure
that we are going to work with the one we want since with this package
we could study several data sets at the same time.

``` r
myDataset <- getDataSetList()
myDataset
#> [1] "zheisel"
```

In this case, we only have zheisel then there will be no problem but it
is always good to check it.

Also, if this is the first time you have dealt with a specific dataset,
it is fine to use the **getDataSetInfo** function for a first contact,
since it returns the main characteristics of a specific dataset.

``` r
getDataSetInfo(zheisel)
#>                 Version NumberCells NumberGenes NumberTypeTissues
#> zheisel_Info 0.2.1.9000      160796       27998                23
#>              NumberTypeCells
#> zheisel_Info               7
```

Once focused, we can proceed to solve the question posed. To do this, we
need to know the exact name by which the cell type and the tissue in
question has been denoted in this dataset since we will need it for
later. To do this we use the functions getCellTypesInDataSet (dataset)
and **getTissueTypesInDataSet**, which return a list of the cell types
and tissues available for said dataset, respectively.

``` r
cellTypes <- getCellTypesInDataSet(zheisel)
myCellType <- cellTypes[grep("[A-a]stro?", cellTypes)]
myCellType
#> [1] "Astrocytes"
```

``` r
tissueTypes <- getTissueTypesInDataSet(zheisel)
myTissueType <- tissueTypes[grep("[A-a]mygd?", tissueTypes)]
myTissueType
#> [1] "Amygd"
```

``` r
cat("In our case, we have to keep", myCellType, "and", myTissueType, "written exactly this way."
)
#> In our case, we have to keep Astrocytes and Amygd written exactly this way.
```

The next step is to know the type of method that we can apply to
evaluate if the Gfap gene is really expressed in the astrocyes within
the amygdala. So, we apply the function **getActiveGeneMethod**.

``` r
myMethods <- getActiveGeneMethod (zheisel)
myMethods
#> [1] "VariationCoefficient" "RawExpression"
```

The result is that we can apply both the raw expression and variation
coefficient methods to our dataset, so we are going to apply both as we
can (for more information on these methods consult the documentation:
“Processing raw data from Zheisel to get cell type specificity for
genes” by Juan A. Botía).

Once we have all the information we need, we can proceed to estimate
whether the Gfap gene is active in amygdala astrocytes with these two
available methods. In order to implement the function that performs it,
called *getActiveGenesInCellType*, we first pre-calculate the active
gene matrixes for each of the tissues. This process will be carried out
once per dataset in order to have all the information available to ask
all the questions we want, although we will return to our question
later. As we were saying, each matrix will be composed of as many rows
as cell types, each one composed of as many 1 as active genes and as
many 0 as inactive genes for this tissue. In this case, we will create
two matrixes for each of the tissues since we will execute each of the
available methods for this data set. Finally, each of these matrixes
will be saved as an rds file with the name “dataset\_tissue\_method.rds”
(these files are also available in the repository).

The steps to obtain these matrixes are detailed below:

**1.Collect the information necesary**. On the one hand, we will need
the expression data matrix, already available from the dataset. On the
other hand, we will need the tissue types, cell types and available
methods, which will be obtained by calling the previously described
functions. Obtained the information, we are going to verify that it is
correct, otherwise, the function would end here.

**2.Create each matrix**. From the main matrix of this dataset, we will
create one matrix composed of as many rows as cells belong to a tissue
and another one composed of as many rows as cells belong to a cell type.
Next, we will select the matching rows between them, that is, they
belong to a specific tissue and a specific cell type. This matrix will
be used to estimate the active genes using two methods: RawExpression
and VariationCoefficient. From each method, we will obtain the ids of
the active genes, that is, those that are expressed above 0.5 for a
visibility value of 0.05 (default value) or that have a coefficient of
variation above 0.05, depending on the method used.

**3.Build a boolean matrix.** At this point, we have the ids of the
active genes for a specific cell type within a specific tissue obtained
with each method for a specific dataset. But we want to build a boolean
matrix, then we have to represent the active genes by ones and zeros. To
do this, we will create a vector of the same length as the number of
columns in our smatrix, that is, the number of genes, and what we will
complete with zeros. Next, we will check the positions occupied by the
active genes calculated in the complete list of genes and replace those
positions with ones. Each created vector will be included as an
independent row in a dataframe, in this way, we will create a matrix
with as many rows as cell types and as many columns as genes. We carry
out this process for each method so we will simultaneously create two
dataframes for each tissue, one for each method.

**4.Repeat the process**. After evaluating all the cell types within a
tissue, we will save the files, initialize the variables again and
continue with the rest of the tissues following the same procedure.

``` r
calcActiveGenesInCellType <-  function(dataset){
 
# 1. Obtain the necessary data to execute the function.
  datasetName <- deparse(substitute(dataset))
  mymatrix <- dataset$matrix[, ]
 
  methods <- getActiveGeneMethod(dataset)
  tissues <- getTissueTypesInDataSet(dataset)
  cellTypes <- getCellTypesInDataSet(dataset)
  
  colnames(mymatrix) <- dataset$row.attrs$Gene[]
  rownames(mymatrix) <- dataset$col.attrs$CellID[]
  genes <- colnames(mymatrix)
  
  cv = function (x, na.rm = FALSE) {
        sd(x, na.rm = na.rm)/mean(x, na.rm = na.rm)}
  
  
# 2. Check if all data is available.
  
  if (is.null(methods)) stop("No methods available")
  if (is.null(tissues)) stop("No tissues available")
  if (is.null(cellTypes)) stop("No cell types available")
 
  
# 3. Access the expression data of each tissue.
  
  for (tissue in 1:length(tissues)) {
    
    tissueMatrix <- mymatrix[dataset$col.attrs$Tissue[] == tissues[tissue], ]
    example1 <- mymatrix[dataset$col.attrs$Tissue[] == tissues[1], ]
    
    # Two files will be created for each tissue, one for each method: RawExpression (RE) and VariationCoefficient (VC).
    activeGenesRETissue <- as.data.frame(matrix(numeric(), nrow = 0, ncol = length(genes)))
    activeGenesVCTissue <- as.data.frame(matrix(numeric(), nrow = 0, ncol = length(genes)))
       
    
# 3.1. And within each tissue, access the expression data of each cell type.
    
     for (type in 1:length(cellTypes)) {
       
      cellMatrix  <- mymatrix[dataset$col.attrs$Class[] == cellTypes[type], ]
      
      matchRows <- intersect(rownames(cellMatrix), rownames(tissueMatrix))
      tissueCellMatrix  <- as.matrix(cellMatrix[matchRows, ])
      
      
      # For each cell type, the active genes for each method will be saved (RE and VC).
      activeGenesRECellType <- c()
      activeGenesVCCellType <- c()
      
      
# 3.1.1 Apply RawExpression method.
      
      visibility  <- 0.05 # default value
      
      mask = colSums(tissueCellMatrix > 0.5)  > (visibility * nrow(tissueCellMatrix))
      IDactiveGenesRECellType <- genes[mask] 
      indexActiveGenesRECellTYpe <- match(IDactiveGenesRECellType, genes) 
      activeGenesRECellType <- replace(rep(0, length(genes)), indexActiveGenesRECellTYpe, 1)
      
     
# 3.1.2 Apply VariationCoefficient method.
      
      cvs <- apply(tissueCellMatrix, 2, cv, na.rm = T)
      IDactiveGenesVCCellType <- names(cvs[cvs > 0.05]) 
      indexActiveGenesVCCellTYpe  <- match(IDactiveGenesVCCellType, genes)
      activeGenesVCCellType <- replace(rep(0, length(genes)), indexActiveGenesVCCellTYpe, 1)
      
      
# 3.1.3. Add each vector as a row of its corresponding dataframe depending on the method.
      
      activeGenesRETissue <- rbind(activeGenesRETissue, activeGenesRECellType) 
      activeGenesVCTissue <- rbind(activeGenesVCTissue, activeGenesVCCellType)
     }
     
# 4. Save the two files created for each tissue with the corresponding name.
      
      rownames(activeGenesRETissue) <- cellTypes
      colnames(activeGenesRETissue) <- genes
      saveRDS(activeGenesRETissue, paste0(datasetName, "_", 
                                         tissues[tissue], "_", 
                                         "RawExpression.rds"))
      
      rownames(activeGenesVCTissue) <- cellTypes
      colnames(activeGenesVCTissue) <- genes
      saveRDS(activeGenesVCTissue, paste0(datasetName, "_", 
                                          tissues[tissue], "_", 
                                          "VariationCoefficient.rds"))
    }
  
  
# 5. Repeat the process for each tissue. 
  
  return("Completed process")  
}
```

``` r
calcActiveGenesInCellType(zheisel)
#> [1] "Completed process"
```

Therefore, once the matrixes are pre-calculated for each tissue, we will
be able to answer our question applying the function
**getActiveGenesInCellType**.

``` r
gene <- "Gfap"

for (method in myMethods) {
  
  results <- getActiveGenesInCellType(myDataset, myTissueType, myCellType, method)

  if (length(grep(gene, results
                  )) == 0) {
    cat("-", gene, "gene is not expressed in the", myCellType, "of the", myTissueType, 
        "of the data set", myDataset, "with the method", method, ".", "\n\n")
    
  } else {
    
     cat("-", gene, "gene is expressed in the", myCellType, "of the", myTissueType, 
         "of the data set", myDataset, "with the method", method, ".", "\n\n")
  }
}
#> - Gfap gene is expressed in the Astrocytes of the Amygd of the data set zheisel with the method VariationCoefficient . 
#> 
#> - Gfap gene is expressed in the Astrocytes of the Amygd of the data set zheisel with the method RawExpression .
```

In this way, questions of this type could be raised and resolved with
deconvML package for single-cell data sets, an increasingly used type of
data set.

References:
-----------

Trapnell, Cole. 2015. “Defining Cell Types and States with Single-Cell
Genomics.” *Genome Research* 25 (10). Cold Spring Harbor Lab: 1491–8.
