library(glmnet)
# name dataset
dataset = t(Kidney_2)
# set CV times as cv_times
cv_times = 10

# CV
# total mse in all cv
cv_mse = c()
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
    
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData <- dataset[testIndexes, ]
    trainData <- dataset[-testIndexes, ]
    
    # train model
    target = trainData[,"Mapk1"]
    features = trainData[,colnames(trainData)!="Mapk1"]
    cvfit_lasso = cv.glmnet(features, target, alpha=1,nfolds = 10)
    model = glmnet(features,target, alpha=1,lambda=cvfit_lasso$lambda.min)
    # test model
    target = testData[,"Mapk1"]
    features = testData[,colnames(testData)!="Mapk1"]
    prediction = predict(model, newx = features)
    mse=c(mse,mean((prediction-target)^2))
    }
  print(mean(mse))
  cv_mse = c(cv_mse,mean(mse))
  
}
print("total error")
print(mean(cv_mse))
