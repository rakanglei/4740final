library(randomForest)
# name dataset
dataset = t(Kidney_2)
# set CV times as cv_times
cv_times = 10

# CV
cv_mse=c()
for (seed_num in 1:cv_times){
  set.seed(seed_num)
  print(seed_num)
  # One time cv
  # Randomly shuffle the data
  dataset<-dataset[sample(nrow(dataset)),]

  #Create 10 equally size folds
  folds <- cut(seq(1,nrow(dataset)),breaks=10,labels=FALSE)
  #Perform 10 fold cross validation
  mse = c()
  par(mfrow=c(5,2))
  for(i in 1:10){
    #target = dataset["Mapk1",]
    #features = dataset[rownames(Kidney_2)!="Mapk1",]
    #Segement your data by fold using the which() function
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData <- dataset[testIndexes, ]
    trainData <- dataset[-testIndexes, ]

    # train model
    rf.fit = randomForest(Mapk1~., data = trainData, importance=TRUE, ntree=100)
    varImpPlot(rf.fit)
    # test model
    prediction = predict(rf.fit, newdata = testData)
    real = testData[, "Mapk1"]
    mse=c(mse,mean((prediction-real)^2))
  }
  print(mean(mse))
  cv_mse = c(cv_mse,mean(mse))
}
print("Total error")
print(mean(cv_mse))
# 
# # One time cv
# # Randomly shuffle the data
# dataset<-dataset[sample(nrow(dataset)),]
# #Create 10 equally size folds
# folds <- cut(seq(1,nrow(dataset)),breaks=10,labels=FALSE)
# #Perform 10 fold cross validation
# mse = c()
# par(mfrow=c(5,2))
# for(i in 1:10){
#   #target = dataset["Mapk1",]
#   #features = dataset[rownames(Kidney_2)!="Mapk1",]
#   #Segement your data by fold using the which() function
#   testIndexes <- which(folds==i,arr.ind=TRUE)
#   testData <- dataset[testIndexes, ]
#   trainData <- dataset[-testIndexes, ]
# 
#   # train model
#   rf.fit = randomForest(Mapk1~., data = trainData, importance=TRUE, ntree=500)
#   varImpPlot(rf.fit)
#   # test model
#   prediction = predict(rf.fit, newdata = testData)
#   real = testData[, "Mapk1"]
#   mse=c(mse,mean((prediction-real)^2))
# }
# print(mean(mse))



# varImpPlot(rf.fit) #look at what variables are important
# prediction = predict(rf.fit, trans_Kidney_2.test)
# real = trans_Kidney_2.test[, "Mapk1"]
# mean((prediction-real)^2)
# 
# #exp 1
# rf1.fit = randomForest(Mapk1 ~ Cdc42+Plcg2+Rik+Mapkapk2+Pik3cd+Pla2g5+Map2k1+Pik3r3+Nras+Pik3r1+Pik3ca+Ppp3cb+Mapk13+Rac1+Nfat5, data = trans_Kidney_2.train, importance=TRUE, ntree=2000)
# varImpPlot(rf1.fit) #look at what variables are important
# prediction1 = predict(rf1.fit, trans_Kidney_2.test)
# real1 = trans_Kidney_2.test[, "Mapk1"]
# mean((prediction1-real1)^2)
