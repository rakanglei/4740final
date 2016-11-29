Kidney_2["Mapk1",]
summary(Kidney_2)
install.packages('randomForest')
set.seed(666)
train = sample(40, 20) 
Kidney_2.train = Kidney_2[,train]
trans_Kidney_2.train = t(Kidney_2.train)
Kidney_2.test = Kidney_2[,-train]
trans_Kidney_2.test = t(Kidney_2.test)

library(randomForest)
rf.fit = randomForest(Mapk1 ~ Cdc42+Pla2g6+Akt2+Plcg2+Rac2+Rik+Mapkapk2+Pik3cd+Pla2g5+Sphk2+Map2k1+Pik3r3+Ptk2+Nras+Nos3+Pik3r1+Pik3ca+Ppp3cb+Map2k2+Nfatc4+Mapk13+Rac1+Nfat5, data = trans_Kidney_2.train, importance=TRUE, ntree=2000)
varImpPlot(rf.fit) #look at what variables are important
prediction = predict(rf.fit, trans_Kidney_2.test)
real = trans_Kidney_2.test[, "Mapk1"]
mean((prediction-real)^2)

#exp 1
rf1.fit = randomForest(Mapk1 ~ Cdc42+Plcg2+Rik+Mapkapk2+Pik3cd+Pla2g5+Map2k1+Pik3r3+Nras+Pik3r1+Pik3ca+Ppp3cb+Mapk13+Rac1+Nfat5, data = trans_Kidney_2.train, importance=TRUE, ntree=2000)
varImpPlot(rf1.fit) #look at what variables are important
prediction1 = predict(rf1.fit, trans_Kidney_2.test)
real1 = trans_Kidney_2.test[, "Mapk1"]
mean((prediction1-real1)^2)

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