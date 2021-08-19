####################################################################################################
### loading required packages
####################################################################################################
rm(list = ls())
## Libraries and variables
{
  pkg<-c("ggplot2",    
         "dplyr",
         "tibble",    
         "factoextra", 
         "stringr",
         "caret",
         "plyr",
         "pryr",
         "gam",
         "randomForest")        
  lapply(pkg, function(x) {
    if(!require(x,character.only = TRUE)) install.packages(x,character.only = TRUE)
    library(x,character.only = TRUE)
  })
}



####################################################################################################
### Preparation of Auxiliary functions
####################################################################################################

multivariateFeatureElimination <- function(rf.ctrl) {
  preProc = c("center", "scale")
  rfProfile <- rfe(x=x.fs, y=y.fs,
                   sizes = c(1:(dim(dataTrain)[1]-1)),
                   preProc = preProc,
                   rfeControl = rf.ctrl)
  return(rfProfile)
}

functionSummary <- function (data, lev = NULL, model = NULL) {
  if (is.character(data$obs)) 
    data$obs <- factor(data$obs, levels = lev)
  postResample(data[, "pred"], data[, "obs"])
}

functionFit <- function(x, y, first, last, ...) {
  loadNamespace("randomForest")
  randomForest::randomForest(x, y, importance = first, ...)
}

functionPred <- function (object, x) {
  if (class(object) == "nullModel") {
    tmp <- predict(object, x)
    if (!is.null(object$levels)) {
      out <- cbind(data.frame(pred = tmp), as.data.frame(predict(object, 
                                                                 x, type = "prob")))
    }
    else out <- tmp
  }
  else {
    tmp <- predict(object, x)
    if (is.factor(object$y)) {
      out <- cbind(data.frame(pred = tmp), as.data.frame(predict(object, 
                                                                 x, type = "prob")))
    }
    else out <- tmp
  }
  out
}

functionScore <- function (x, y) {
  if (is.factor(y)) 
    anovaScores(x, y)
  else gamScores(x, y)
}

functionFilter <- function (score, x, y) { 
  score <= 0.05
}

functionRank <- function(object, x, y) {
  vimp <- varImp(object)
  if (is.factor(y)) {
    if (all(levels(y) %in% colnames(vimp))) {
      avImp <- apply(vimp[, levels(y), drop = TRUE], 1, 
                     mean)
      vimp$Overall <- avImp
    }
  }
  vimp <- vimp[order(vimp$Overall, decreasing = TRUE), , drop = FALSE]
  vimp$var <- rownames(vimp)
  vimp
}


####################################################################################################
### Loading data 
####################################################################################################

load("./data/datMeta_FeatureSelection.Rdata")

## Creating data matrix
mirna.data <- datMeta_FeatureSelection[,grep(x=colnames(datMeta_FeatureSelection), pattern="mmu", value=TRUE)]
phenotypic.data <- datMeta_FeatureSelection[,grep(x=colnames(datMeta_FeatureSelection), pattern="mmu", value=TRUE, invert=TRUE)]

### PCA
pca.phenotypic <- prcomp(phenotypic.data, center=TRUE, scale=TRUE)

## Making PC1 the dependent variable
dataset <- merge(x=mirna.data, y=subset(pca.phenotypic$x, select=PC1), by="row.names")
colnames(dataset) <- gsub(x=colnames(dataset), pattern="PC1", replacement="DepVar")
rownames(dataset) <- dataset$Row.names
dataset <- subset(dataset, select=-Row.names)

### Defining training data
dataTrain <- dataset 
x.fs <- subset(x=dataTrain, select=-DepVar)
y.fs <- unlist(subset(x=dataTrain, select=DepVar))
N <- nrow(x.fs)
D <- ncol(x.fs)

####################################################################################################
### Recursive Feature Elimination (RFE)
####################################################################################################

print(paste("[MAIN]: start recursive feature elimination", sep=""))

## bootstrap
set.seed(10)
message("performing RF multivariate bootstrapping...")

args.RFE <- list(functions = rfFuncs,
                 method = "boot",
                 number = 250,
                 repeats = 1,
                 p = 0.75, 
                 returnResamp = "all",
                 verbose = FALSE)
rf.ctrl <- do.call("rfeControl", args.RFE)
rf.ctrl$functions$fit <- functionFit
rf.ctrl$functions$rank <- functionRank

bootstrap.rfProfile <- multivariateFeatureElimination(rf.ctrl)
bootstrap.res <- varImp(bootstrap.rfProfile)

set.seed(10)
message("performing RF leave one out cross validation...")
  
args.RFE <- list(functions = rfFuncs,
                 method = "LOOCV",
                 number = 1,
                 repeats = 1,
                 returnResamp = "all",
                 verbose = FALSE)
args.RFE$functions$summary <- functionSummary
args.RFE$functions$fit <- functionFit
args.RFE$functions$pred <- functionPred
args.RFE$functions$score <- functionScore
args.RFE$functions$filter <- functionFilter

rf.ctrl <- do.call("rfeControl", args.RFE)

loocv.rfProfile <- multivariateFeatureElimination(rf.ctrl)
loocv.res <- varImp(loocv.rfProfile)


### SVM
set.seed(12345)

message("performing SVM...")


svmProfile <- rfe(x=x.fs, y=y.fs,
                   sizes = c(2, 5, 10, 20),
                   rfeControl = rfeControl(functions = caretFuncs, number = 100),
                   method = "svmRadial")

svm.res <- varImp(svmProfile)

#select top 10 microRNAs from each approach
bootstrap = rownames(bootstrap.res)[1:10]
loocv=rownames(loocv.res)[1:10]
svm=rownames(svm.res)[1:10]

list = list(bootstrap = bootstrap, loocv = loocv, svm = svm)

####################################################################
####################################################################
message("list of common microRNAs reported in Figure 2F")
#####################################################################
#####################################################################

print(Reduce(intersect, list))

sessionInfo()

