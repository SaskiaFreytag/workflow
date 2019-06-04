library(caret)
set.seed(3456)

ref_sce <- readRDS(here::here("data/meta/", "allen_mouse_atlas.rds"))

ref_train <- data.frame(colData(ref_sce)$subclass,
                        t(counts(ref_sce)[rowData(ref_sce)$hvg,]))
colnames(ref_train)[1] <- "subclass"

# filter stuff
shit <- c(names(table(ref_train$subclass))[table(ref_train$subclass)<10],
          "Batch Grouping", "Doublet", "High Intronic", "Low Quality",
          "No Class") 

ind_shit <- match(ref_train$subclass, shit)
ind_shit <- which(!is.na(ind_shit))

ref_train <- ref_train[-ind_shit,]
ref_train$subclass <- as.factor(as.character(ref_train$subclass ))

trainIndex <- createDataPartition(ref_train$subclass, p = .8, 
                                  list = FALSE, 
                                  times = 1)

refTrain <- ref_train[ trainIndex,]
refTest  <- ref_train[-trainIndex,]

descrCor <-  cor(refTest[,-1])
summary(descrCor[upper.tri(descrCor)])

highlyCorDescr <- findCorrelation(descrCor, cutoff = .75)
filteredrefTest <- refTest[,-(highlyCorDescr+1)]

fitControl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 1)

filteredrefTrain <- refTrain[, -c(1,highlyCorDescr+1)]

svmFit <- train(subclass ~ ., data = filteredrefTest, 
                method = "svmRadial", 
                trControl = fitControl, 
                preProc = c("center", "scale"),
                tuneLength = 10)

pred <- predict(svmFit, newdata = filteredrefTrain)

confusionMatrix(pred,refTrain$subclass)


