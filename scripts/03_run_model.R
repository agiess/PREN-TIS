#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library (h2o)
library (ggplot2)
library (data.table)

n.threads<-as.integer(args[15])
my_seed<-as.integer(args[16])

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#start a h2o instance
h2o.init(nthreads = n.threads)
h2o.removeAll()


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

training.pos.raw<-read.csv(args[1], header=TRUE)
training.neg.raw<-read.csv(args[2], header=TRUE)
testing.pos.raw<-read.csv(args[3], header=TRUE)
testing.neg.raw<-read.csv(args[4], header=TRUE)

rownames(training.pos.raw) <- training.pos.raw[, 1]
training.pos.raw <- training.pos.raw[,-1]
training.pos.raw$seq_1[training.pos.raw$seq_1 == TRUE] <-"T"
training.pos.raw[,13:53] <- lapply(training.pos.raw[,13:53] , factor)
training.pos.raw[,54:ncol(training.pos.raw)] <- lapply(training.pos.raw[,54:ncol(training.pos.raw)] , as.numeric)

rownames(training.neg.raw) <- training.neg.raw[, 1]
training.neg.raw <- training.neg.raw[,-1]
training.neg.raw$seq_1[training.neg.raw$seq_1 == TRUE] <-"T"
training.neg.raw[,13:53] <- lapply(training.neg.raw[,13:53] , factor)
training.neg.raw[,54:ncol(training.pos.raw)] <- lapply(training.neg.raw[,54:ncol(training.pos.raw)] , as.numeric)

rownames(testing.pos.raw) <- testing.pos.raw[, 1]
testing.pos.raw <- testing.pos.raw[,-1]
testing.pos.raw$seq_1[testing.pos.raw$seq_1 == TRUE] <-"T"
testing.pos.raw[,13:53] <- lapply(testing.pos.raw[,13:53] , factor)
testing.pos.raw[,54:ncol(testing.pos.raw)] <- lapply(testing.pos.raw[,54:ncol(testing.pos.raw)] , as.numeric)

rownames(testing.neg.raw) <- testing.neg.raw[, 1]
testing.neg.raw <- testing.neg.raw[,-1]
testing.neg.raw$seq_1[testing.neg.raw$seq_1 == TRUE] <-"T"
testing.neg.raw[,13:53] <- lapply(testing.neg.raw[,13:53] , factor)
testing.neg.raw[,54:ncol(testing.pos.raw)] <- lapply(testing.neg.raw[,54:ncol(testing.pos.raw)] , as.numeric)

train.raw <- rbind(training.pos.raw, training.neg.raw) 
test.raw <- rbind(testing.pos.raw, testing.neg.raw)

train.raw$y=as.factor(train.raw$annotated_start_site)
levels(train.raw$y) =c('neg','pos')
test.raw$y=as.factor(test.raw$annotated_start_site)
levels(test.raw$y) =c('neg','pos')

train.hex <- as.h2o(train.raw, destination_frame="train.hex")
test.hex <- as.h2o(test.raw, destination_frame = "test.hex")

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#run models

hyper_params = list( ntrees = seq(100, 1000, 100) )

tree_grid_10fold <- h2o.grid(
  hyper_params = hyper_params,
  search_criteria = list(strategy = "Cartesian"),
  algorithm="randomForest",
  grid_id="tree_grid_10fold",
  x = c(1, 3, 9:32, 36:(ncol(train.hex-1)) ),
  y = ncol(train.hex), 
  training_frame = train.hex, 
  nfolds=10,
  seed = my_seed,                                                             
  balance_classes = FALSE                                             
)

grid.tree.10fold <- h2o.getGrid("tree_grid_10fold", sort_by = "auc", decreasing = TRUE)

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
write.csv(rf.out,file=args[8])

#output random forest parameter important values
write.csv(tree.10fold.var.imp,file=args[9])


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
write.csv(as.data.frame(score.matrix), file=args[10])
write.csv(as.data.frame(confusion.matrix), file=args[11])
write.csv(test.out,file=args[12])

#output roc plot from testing set
ggsave(file=,args[13], roc_plot)


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#Stop 2 stop

stop2stop<-read.csv(args[5], header=TRUE)

#rownames after filter
rownames(stop2stop) <- stop2stop[, 1]
stop2stop <- stop2stop[,-1]
stop2stop[,13:53] <- lapply(stop2stop[,13:53] , factor)
stop2stop[,54:ncol(training.pos.raw)] <- lapply(stop2stop[,54:ncol(training.pos.raw)] , as.numeric)

#mark the positions that were used for model training
stop2stop$used_in_model_training <- 'no'
common <- intersect(rownames(train.raw), rownames(stop2stop))
stop2stop[common,]$used_in_model_training <- 'yes'

stop2stop$y=as.factor(stop2stop$annotated_start_site)
levels(stop2stop$y) =c('neg','pos')

stop2stop.hex <- as.h2o(stop2stop, destination_frame="stop2stop.hex")

stop2stop.RF.10fold <- h2o.predict(best.model.tree.10fold, stop2stop.hex)
stop2stop$pred       <- as.data.frame(stop2stop.RF.10fold$predict)[,1]
stop2stop$pred_pos   <- as.data.frame(stop2stop.RF.10fold$pos)[,1]

#turn off expontional notation
stop2stop$pred_pos<-format(stop2stop$pred_pos, scientific = FALSE)

#write predictions
write.csv(stop2stop,file=args[14])


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
