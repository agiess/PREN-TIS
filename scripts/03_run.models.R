#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library (h2o)
library (ggplot2)
library (data.table)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#start a h2o instance
h2o.init(nthreads = 1, max_mem_size="100g")
h2o.removeAll()


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

pos.raw<-read.csv(args[1], header=TRUE)
neg.raw<-read.csv(args[2], header=TRUE)

rownames(pos.raw) <- pos.raw[, 1]
pos.raw <- pos.raw[,-1]
pos.raw$seq_1[pos.raw$seq_1 == TRUE] <-"T"
pos.raw[,13:53] <- lapply(pos.raw[,13:53] , factor)
pos.raw[,54:ncol(pos.raw)] <- lapply(pos.raw[,54:ncol(pos.raw)] , as.numeric)

rownames(neg.raw) <- neg.raw[, 1]
neg.raw <- neg.raw[,-1]
neg.raw$seq_1[neg.raw$seq_1 == TRUE] <-"T"
neg.raw[,13:53] <- lapply(neg.raw[,13:53] , factor)
neg.raw[,54:ncol(pos.raw)] <- lapply(neg.raw[,54:ncol(pos.raw)] , as.numeric)

#the sets need to be randomised
set.seed(9999)
pos.ran <- pos.raw[sample(nrow(pos.raw)),]
neg.ran <- neg.raw[sample(nrow(neg.raw)),]

#set 1
pos.train.raw <- pos.ran[1:1500,] 
neg.train.raw <- neg.ran[1:6000,]
train.raw <- rbind(pos.train.raw, neg.train.raw)

pos.test.raw <- pos.ran[1501:1700,]  
neg.test.raw <- neg.ran[6001:6800,]
test.raw <- rbind(pos.test.raw, neg.test.raw)

train.raw$y=as.factor(train.raw$annotated_start_site)
levels(train.raw$y) =c('neg','pos')
test.raw$y=as.factor(test.raw$annotated_start_site)
levels(test.raw$y) =c('neg','pos')

train.hex <- as.h2o(train.raw, destination_frame="train.hex")
test.hex <- as.h2o(test.raw, destination_frame = "test.hex")


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#h2o without parameter scaling cross validation

alpha_opts = list(list(0), list(.2), list(.4), list(.6), list(.8), list(1))
hyper_parameters = list(alpha = alpha_opts)

grid.10fold <- h2o.grid(
  grid_id="glm_grid_10fold",
  algorithm="glm",
  family = "binomial",
  hyper_params = hyper_parameters, 
  y = ncol(train.hex),
  x = c(1, 3, 9:32, 36:(ncol(train.hex-1)) ),
  training_frame = train.hex, 
  standardize = TRUE,
  intercept = FALSE,
  lambda_search = TRUE, 
  nfolds=10,
  seed = 7777777,
  balance_classes = FALSE      
)

grid.models.10fold <- lapply(grid.10fold@model_ids, function(model_id) { model = h2o.getModel(model_id) })

grid.10fold.summary <- h2o.getGrid("glm_grid_10fold", sort_by = "auc", decreasing = TRUE)
model.10fold.ids <- grid.10fold.summary@model_ids
best.model.10fold <- h2o.getModel(model.10fold.ids[[1]])
best.10fold <- best.model.10fold@model$coefficients_table
sort.10fold <- as.data.frame(best.model.10fold@model$coefficients_table)
sorted.10fold <- sort.10fold[order(-sort.10fold$coefficients),]

#take non zero predictors
selected.10fold<-sorted.10fold[sorted.10fold$coefficients != 0, 1]
'%ni%'<-Negate('%in%')
selected.10fold<-selected.10fold[selected.10fold %ni% 'Intercept']

glm.xval <- h2o.performance(best.model.10fold, newdata = NULL, xval = TRUE, data = NULL)

#write model summary
glm.out = matrix( 
  c(glm.xval@metrics$MSE, 
    glm.xval@metrics$r2, 
    glm.xval@metrics$logloss, 
    glm.xval@metrics$mean_per_class_error, 
    glm.xval@metrics$AUC, glm.xval@metrics$Gini, 
    glm.xval@metrics$null_deviance, 
    glm.xval@metrics$residual_deviance, 
    glm.xval@metrics$AIC), 
  nrow=9, 
  ncol=1
)

rownames(glm.out) <- c("mse", "r2", "logloss", "mpce", "auc", "gini", "nulldev", "resdev", "aic")
colnames(glm.out) <- c("metric")
glm.out <- format(glm.out, scientific=F)

#output glm model summary
write.csv(glm.out, file=args[4])

#output glm coefficients
write.csv(sorted.10fold, file=args[5])


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#run models

hyper_params = list( ntrees = seq(50,1000,50) )

#tidy up the factor names from the glm so that they match the original column names 
tidy1 <- gsub("(codon).*", "\\1", selected.10fold)
tidy2 <- gsub("(seq_[0-9]+).*", "\\1", tidy1)
tidy3 <- gsub("(seq_.[0-9]+).*", "\\1", tidy2)
tidy4 <- unique(tidy3)

tree_grid_10fold <- h2o.grid(
  hyper_params = hyper_params,
  search_criteria = list(strategy = "Cartesian"),
  algorithm="randomForest",
  grid_id="tree_grid_10fold",
  x = c(tidy4), 
  y = ncol(train.hex), 
  training_frame = train.hex, 
  nfolds=10,
  seed = 7777777,                                                             
  balance_classes = FALSE                                             
)

grid.tree.10fold <- h2o.getGrid("tree_grid_10fold", sort_by = "auc", decreasing = TRUE) #check if this is validation AUC

model.ids.tree.10fold <- grid.tree.10fold@model_ids
best.model.tree.10fold <- h2o.getModel(model.ids.tree.10fold[[1]])
tree.10fold.var.imp <- as.data.frame(h2o.varimp(best.model.tree.10fold))

rf.xval <- h2o.performance(best.model.tree.10fold, newdata = NULL, xval = TRUE, data = NULL)

#write model summary
rf.out = matrix( 
  c(rf.xval@metrics$MSE, 
    rf.xval@metrics$r2, 
    rf.xval@metrics$logloss, 
    rf.xval@metrics$mean_per_class_error, 
    rf.xval@metrics$AUC, 
    rf.xval@metrics$Gini), 
  nrow=6, 
  ncol=1
)

rownames(rf.out) <- c("mse", "r2", "logloss", "mpce", "auc", "gini")
colnames(rf.out) <- c("metric")
rf.out <- format(rf.out, scientific=F)

#output random forest model summary
write.csv(rf.out,file=args[6])

#output random forest parameter important values
write.csv(tree.10fold.var.imp,file=args[7])


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#Evalulate the performance on the testing set
rf.test <- h2o.performance(best.model.tree.10fold, newdata = test.hex)

confusion.matrix <- h2o.confusionMatrix(rf.test)
score.matrix <- rf.test@metrics$thresholds_and_metric_scores

# area under the curve
auc <- h2o.auc(rf.test)
fpr <- h2o.fpr(rf.test)
tpr <- h2o.tpr(rf.test)
roc_plot<-ggplot( data.table(fpr = fpr, tpr = tpr), aes(fpr.fpr, tpr.tpr) ) + 
                  geom_line() + 
                  theme_bw() + 
                  ylab("tpr") +
                  xlab("fpr") +
                  ggtitle( sprintf('AUC: %f', auc) )

#write model summary
test.out = matrix( 
  c(rf.test@metrics$MSE, 
    rf.test@metrics$r2, 
    rf.test@metrics$logloss, 
    rf.test@metrics$mean_per_class_error, 
    rf.test@metrics$AUC, 
    rf.test@metrics$Gini), 
  nrow=6, 
  ncol=1
)

rownames(test.out) <- c("mse", "r2", "logloss", "mpce", "auc", "gini")
colnames(test.out) <- c("metric")
test.out <- format(test.out, scientific=F)

#output metrics on testing set
write.csv(as.data.frame(score.matrix), file=args[8])
write.csv(as.data.frame(confusion.matrix), file=args[9])
write.csv(test.out,file=args[10])


#output roc plot from testing set
ggsave(file=,args[11], roc_plot)


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
### Stop 2 stop

stop2stop<-read.csv(args[3], header=TRUE)

#rownames after filter
rownames(stop2stop) <- stop2stop[, 1]
stop2stop <- stop2stop[,-1]
stop2stop[,13:53] <- lapply(stop2stop[,13:53] , factor)
stop2stop[,54:ncol(pos.raw)] <- lapply(stop2stop[,54:ncol(pos.raw)] , as.numeric)

stop2stop$y=as.factor(stop2stop$annotated_start_site)
levels(stop2stop$y) =c('neg','pos')

stop2stop.hex <- as.h2o(stop2stop, destination_frame="stop2stop.hex")

stop2stop.RF.10fold <- h2o.predict(best.model.tree.10fold, stop2stop.hex)
stop2stop$pred       <- as.data.frame(stop2stop.RF.10fold$predict)[,1]
stop2stop$pred_pos   <- as.data.frame(stop2stop.RF.10fold$pos)[,1]

#turn off expontional notation
stop2stop$pred_pos<-format(stop2stop$pred_pos, scientific = FALSE)

#write predictions
write.csv(stop2stop,file=args[12])


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
