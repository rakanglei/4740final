Data <- read.csv(file = "/Users/minghanchang/Documents/Cornell/fall/DM and ML/assignment/final project/Kidney_2.csv")
names(Data)
summary(Data)
attach(Data)

#Kidney_2["Mapk1",]
#summary(Kidney_2)
install.packages('randomForest')
#set.seed(666)
#train = sample(40, 20) 
#Kidney_2.train = Kidney_2[,train]
#trans_Kidney_2.train = t(Kidney_2.train)
#Kidney_2.test = Kidney_2[,-train]
#trans_Kidney_2.test = t(Kidney_2.test)

library(randomForest)
#rf.fit = randomForest(Mapk1 ~ Cdc42+Pla2g6+Akt2+Plcg2+Rac2+Rik+Mapkapk2+Pik3cd+Pla2g5+Sphk2+Map2k1+Pik3r3+Ptk2+Nras+Nos3+Pik3r1+Pik3ca+Ppp3cb+Map2k2+Nfatc4+Mapk13+Rac1+Nfat5, data = trans_Kidney_2.train, importance=TRUE, ntree=2000)
#rf.fit = randomForest(Mapk1 ~ Cdc42+Pla2g6+Akt2+Plcg2+Rac2+Rik+Mapkapk2+Pik3cd+Pla2g5+Sphk2+Map2k1+Pik3r3+Ptk2+Nras+Nos3+Pik3r1+Pik3ca+Ppp3cb+Map2k2+Nfatc4+Mapk13+Rac1+Nfat5, data = Data, importance=TRUE, ntree=2000)
#

#varImpPlot(rf.fit) #look at what variables are important
#summary(rf.fit)
#prediction = predict(rf.fit, trans_Kidney_2.test)
#real = trans_Kidney_2.test[, "Mapk1"]
#mean((prediction-real)^2)

#exp 1
#rf1.fit = randomForest(Mapk1 ~ Cdc42+Plcg2+Rik+Mapkapk2+Pik3cd+Pla2g5+Map2k1+Pik3r3+Nras+Pik3r1+Pik3ca+Ppp3cb+Mapk13+Rac1+Nfat5, data = trans_Kidney_2.train, importance=TRUE, ntree=2000)
#varImpPlot(rf1.fit) #look at what variables are important
#prediction1 = predict(rf1.fit, trans_Kidney_2.test)
#real1 = trans_Kidney_2.test[, "Mapk1"]
#mean((prediction1-real1)^2)


# test variation of the restuls (set.seed from 1-9)
k = 10        # number of folds (try 5 too)： 10 is more stable
par(mfrow=c(3,3))
for (n in 1:9) {
  set.seed(n) 
  # Assign each observation to a single fold
  folds = sample(rep(1:k, dim(Data)[1]/k)) #equally divided
  nvar = 8
  # Create a matrix to store the results of our upcoming calculations
  cv_errors = matrix(NA, k, nvar, dimnames = list(NULL, paste(1:nvar)))
  
  for(j in 1:k){
    # The perform best subset selection on the full dataset, minus the jth fold
    rf.fit = randomForest(Data$Mapk1~., data = Data[folds!=j,], importance=TRUE, ntree=2000)
 
    # Inner loop iterates over each size i
    for(i in 1:nvar){
      # Predict the values of the current fold from the "best subset" model on i predictors
      pred = predict(rf.fit, Data[folds==j,], id=i)
      # Calculate the MSE, store it in the matrix we created above
      cv_errors[j,i] = mean((Data$Mapk1[folds==j]-pred)^2)
    }
  }
  
  # Take the mean of over all folds for each model size
  mean_cv_errors = apply(cv_errors, 2, mean)
  # Find the model size with the smallest cross-validation error
  min = which.min(mean_cv_errors)
  # Plot the cross-validation error for each model size, highlight the min
  plot(mean_cv_errors, type='b', main = paste('method: random forest, k=', k))
  points(min, mean_cv_errors[min][1], col = "red", cex = 2, pch = 20)

}



#"%IncMSE"
#if a predictor is important in your current model, then assigning other values for that predictor randomly but 'realistically' 
#(i.e.: permuting this predictor's values over your dataset), should have a negative influence on prediction, 
#i.e.: using the same model to predict from data that is the same except for the one variable, should give worse predictions.
#So, you take a predictive measure (MSE) with the original dataset and then with the 'permuted' dataset, and you compare them somehow. 
#One way, particularly since we expect the original MSE to always be smaller, the difference can be taken. 
#Finally, for making the values comparable over variables, these are scaled.
#IncNodePurity
#at each split, you can calculate how much this split reduces node impurity 
#(for regression trees, indeed, the difference between RSS before and after the split). 
#This is summed over all splits for that variable, over all trees.