kidney <- read.csv(file = "/Users/yanyulong/Dropbox/1_Cornell/5th_semester/Machine_Learning/kidney_2.csv")
data = setNames(data.frame(t(kidney[,-1])), kidney[,1])

names(data)
summary(data)
attach(data)

# 1. distribution (boxplot) of all the variables
par(mfrow=c(3,8))
for (i in 1:24) {
  outlier_values <- boxplot.stats(eval(as.name(names(data)[i])))$out  # outlier values.
  boxplot(eval(as.name(names(data)[i])), main=names(data)[i])
  mtext(paste('outliers: ', length(outlier_values)), cex=0.6)
}


#--------------------------------------------------------------
library(leaps)
library(glmnet)

# split the data into n-fold, in each training set, the best model is selected; and test MSE is computed using n-fold CV
predict.regsubsets = function(object,newdata,id,...){
  form = as.formula(object$call[[2]]) # Extract the formula used when we called regsubsets()
  mat = model.matrix(form,newdata)    # Build the model matrix
  coefi = coef(object,id=id)          # Extract the coefficiants of the ith model
  xvars = names(coefi)                # Pull out the names of the predictors used in the ith model
  mat[,xvars]%*%coefi                 # Make predictions using matrix multiplication
}

par(mfrow=c(3,3))
nfold = 10; # 5, 8, 10, 20
seed_total = 9;
k = nfold  
model_total = 10 #number of models: 3 for best/fwd/back
mean_cv_errors = matrix(NA,seed_total,model_total)
min = rep(NA, seed_total)

for (nseed in 1:seed_total) {
set.seed(nseed) 
# Randomly shuffle the data
shuffle_ind = sample(nrow(data))
#Create 10 equally size folds
folds = cut(seq(1,nrow(data)),breaks=k,labels=FALSE)
# Create a matrix to store the results of our upcoming calculations
#cv_errors = matrix(NA, k, nvar, dimnames = list(NULL, paste(1:nvar)))
cv_errors = matrix(NA, k, model_total)
bestm_coef = rep(NA,6)

for(j in 1:k){ #every sliptting
  testIndexes = which(folds==j,arr.ind=TRUE)
  test_data   = data[shuffle_ind[testIndexes], ]
  train_data  = data[shuffle_ind[-testIndexes], ]
  
  titles = c('best','fwd','bwd')
  # par(mfrow=c(3,3))
  add = 3
  for (nmodel in 1:add) {
    if (nmodel == 1) {
    model = regsubsets(Mapk1~.,train_data)
    }
    else if (nmodel == 2) {
      model = regsubsets(Mapk1~.,train_data, method = 'forward')
    } else if (nmodel == 3) {
      model = regsubsets(Mapk1~.,train_data, method = 'backward')
      }
    model.summary = summary(model)
    
    # 1. use Cp and BIC to identify best model
    #plot(model.summary$cp,xlab="Number of Variables",ylab="Cp",sub = titles[nmodel], type='l')
    n = which.min(model.summary$cp)
    pred = predict(model, test_data, id = n)
    cv_errors[j,nmodel] = mean((test_data$Mapk1 - pred)^2)
    #points(n[nmodel],model.summary$cp[n[nmodel]],col="red",cex=2,pch=20)
    #bestm_coef[nmodel] = coef(model, n[nmodel]) #get the coefficients of the best model
    
    #plot(model.summary$bic,xlab="Number of Variables",ylab="BIC",sub = titles[nmodel], type='l')
    n = which.min(model.summary$bic)
    pred = predict(model, test_data, id = n)
    cv_errors[j,nmodel + add] = mean((test_data$Mapk1 - pred)^2)
    #points(n[nmodel+add],model.summary$bic[n[nmodel+add]],col="red",cex=2,pch=20)
    #bestm_coef[nmodel+add] = coef(model, n[nmodel+add]) #get the coefficients of the best model
    
    # 3. use CV to identify the best model
    k2 = nfold        
    set.seed(nseed) 
    shuffle_ind2 = sample(nrow(train_data))
    folds2 = cut(seq(1,nrow(train_data)),breaks=k2,labels=FALSE)
    nvar = 8 #number of models
    # Create a matrix to store the results of our upcoming calculations
    cv_errors2 = matrix(NA, k2, nvar)
    
    for(j2 in 1:k2){
      testIndexes2 = which(folds2==j2,arr.ind=TRUE)
      test_data2   = train_data[shuffle_ind2[testIndexes2], ]
      train_data2  = train_data[shuffle_ind2[-testIndexes2], ]
      
      if (nmodel == 1) {
        model2 = regsubsets(Mapk1~.,train_data2)
      }
      else if (nmodel == 2) {
        model2 = regsubsets(Mapk1~.,train_data2, method = 'forward')
      } else if (nmodel == 3) {
        model2 = regsubsets(Mapk1~.,train_data2, method = 'backward')
      }
      # Inner loop iterates over each size i
      for(i2 in 1:nvar){
        # Predict the values of the current fold from the "best subset" model on i predictors
        pred = predict(model2, test_data2, id=i2)
        # Calculate the MSE, store it in the matrix we created above
        cv_errors2[j2,i2] = mean((test_data2$Mapk1-pred)^2)
      }
    }
    
    # Take the mean of over all folds for each model size
    mean_cv_errors2 = apply(cv_errors2, 2, mean)
    # Find the model size with the smallest cross-validation error
    n = which.min(mean_cv_errors2)
    pred = predict(model2, test_data2, id = n)
    cv_errors[j,nmodel + add*2] = mean((test_data$Mapk1 - pred)^2)
    # Plot the cross-validation error for each model size, highlight the min
    #plot(mean_cv_errors2, xlab="Number of Variables",ylab="test MSE",sub = paste(titles[nmodel], 'CV'), type='l')
    #points(n, mean_cv_errors2[n][1], col = "red", cex = 2, pch = 20)
    
    #end of CV
    
  } # end of each model
  
  
  
  # 4. Lasso (the 10th method)
  # train model
  y = train_data[,"Mapk1"]
  x = model.matrix(Mapk1~., train_data)[,-1]
  cvfit_lasso = cv.glmnet(x, y, alpha=1, nfolds = k)
  model = glmnet(x,y, alpha=1,lambda=cvfit_lasso$lambda.min)
  #plot(model)
  # test model
  y = test_data[,"Mapk1"]
  x = model.matrix(Mapk1~., test_data)[,-1]
  pred = predict(model, newx = x)
  cv_errors[j,10]=mean((pred-y)^2)
  
  
  # vary nfold
  x=model.matrix(Mapk1~.,train_data)[,-1]  #remove the intercept
  y=train_data$Mapk1
  
  #par(mfrow=c(3,3))
  ns = 1:10
  fold = c(5, 10, 15, 40)
  #fold = 5:15
  nfold = length(fold)
  var.sd = rep(NA,nfold)
  var.med = rep(NA,nfold)
  ikk = 1
  for (kk in fold) {              
    var.length = rep(NA,length(ns))
    bestlam = rep(NA,length(ns))
    
    i = 1
    for (n in ns) {
      # use cross-validation to choose the bestlam
      set.seed(n) 
      cv.out=cv.glmnet(x, y, alpha=1, nfolds = kk) 
      bestlam[i]=cv.out$lambda.min
      
      out = glmnet(x, y, alpha = 1)  # the same as the lasso.mod
      lasso.coef = predict(out, type = "coefficients", s = bestlam[i])[1:24,]# Display coefficients using lambda chosen by CV
      #lasso.coef=coef(lasso.mod)[,which(grid == bestlam[i])]  #same results
      var.length[i] = length(lasso.coef[lasso.coef!=0]) #number of selected variable
      # print(lasso.coef[lasso.coef!=0])
      i = i+1
    }
    var.sd[ikk] = sd(var.length)
    var.med[ikk] = median(var.length)
    ikk = ikk+1
  }
  
  plot(fold, var.sd, xlab = 'nfold', ylab = 'sd of the number of var in model')
  plot(fold, var.med, xlab = 'nfold', ylab = 'median of the number of var in model')
  
  
  
  } # end of each split
  
  mean_cv_errors[nseed,] = apply(cv_errors, 2, mean) # test_MSE for the first four models
  
  # Find the best method with the smallest cross-validation error
  min[nseed] = which.min(mean_cv_errors[nseed,])  #when nfold = 10, nseed = 1, fwd_Cp give the best 
  # Plot the cross-validation error for each model size, highlight the min
  plot(mean_cv_errors[nseed,], type='b', main = paste('method comparison'), xlab = 'index of model', ylab = 'test MSE')
  points(min[nseed], mean_cv_errors[nseed,][min[nseed]][1], col = "red", cex = 2, pch = 20)  
}

min  

method_cv_errors = apply(mean_cv_errors, 1, mean) 
plot(method_cv_errors, type='b', main = paste('method comparison'), xlab = 'index of model', ylab = 'test MSE')
min2 = which.min(method_cv_errors) 
points(min2, method_cv_errors[min2][1], col = "red", cex = 2, pch = 20)  










#------------------------------------------------------
#check the model: linear assumptions
fit.model5 = lm(Mapk1 ~ Akt2 + Rik + Pik3r3 + Pik3r1 + Rac1, data = data)
plot(fit.model5)

#-----------------------------------------------------
# correlation between variables: colinearity
p_value = matrix(nrow = 24, ncol = 24)
r = matrix(nrow = 24, ncol = 24)
for (i in 1:24){
  for (k in 1:24) {
    if (cor.test(eval(as.name(names(data)[i])), eval(as.name(names(data)[k]))  )$p.value < 0.05) {
      p_value[i,k] = formatC(cor.test(eval(as.name(names(data)[i])), eval(as.name(names(data)[k]))  )$p.value, digit = 2) 
      r[i,k] = formatC(cor.test(eval(as.name(names(data)[i])), eval(as.name(names(data)[k]))  )$estimate, digit = 2)
    }
  }
}
p_value

#3, 7, 13, 17, 23
set = c(3, 7, 13, 17, 23)
cor = matrix(NA, length(set), length(set))
for (i in 1:length(set)) {
  for (j in 1:length(set)) {
    cor[i,j] = ( p_value [set[i],set[j]] )
  }
}
colnames(cor) = c(3, 7, 13, 17, 23)
rownames(cor) = c(3, 7, 13, 17, 23)
cor
# Rik and Pik3r3 is correlated




#-------------------------
# check correlation between each gene and Mapk1: give p and r value
par(mfrow=c(3,8))
for (i in 1:24) {
  plot(Mapk1, eval(as.name(names(data)[i])), ylab = names(data)[i], 
       sub = paste('p: ', formatC(cor.test(Mapk1, eval(as.name(names(data)[i])))$p.value , digit = 2), 
                   ',r: ', formatC(cor(Mapk1, eval(as.name(names(data)[i]))) , digit = 2)) )}


