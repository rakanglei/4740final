Data <- read.csv(file = "/Users/yanyulong/Dropbox/1_Cornell/5th_semester/Machine_Learning/kidney_1.csv")
names(Data)
summary(Data)
attach(Data)

# distribution (boxplot) of all the variables
par(mfrow=c(3,8))
for (i in 1:24) {
  boxplot(eval(as.name(names(Data)[i])), 
          xlab = names(Data)[i])}

# count outliers (may not matter, espeically in high dimensional space)
par(mfrow=c(3,8))
for (i in 1:24) {
outlier_values <- boxplot.stats(eval(as.name(names(Data)[i])))$out  # outlier values.
outlier_indice = numeric(length(outlier_values))
for (k in 1:length(outlier_values)) {
  outlier_indice[k] = which(eval(as.name(names(Data)[i])) == outlier_values[k], arr.ind = FALSE, useNames = TRUE)
}
boxplot(eval(as.name(names(Data)[i])), main=names(Data)[i])
mtext(paste('outliers: ', length(outlier_indice)), cex=0.6)
# if (sum(outlier_indice) > 0 )
# mtext(print(outlier_indice), cex=0.6)
}



# check correlation between each gene and Mapk1: give p and r value
par(mfrow=c(3,8))
for (i in 1:24) {
plot(Mapk1, eval(as.name(names(Data)[i])), ylab = names(Data)[i], 
     sub = paste('p: ', formatC(cor.test(Mapk1, eval(as.name(names(Data)[i])))$p.value , digit = 2), 
                 ',r: ', formatC(cor(Mapk1, eval(as.name(names(Data)[i]))) , digit = 2)) )}


# correlation between variables: colinearity
p_value = matrix(nrow = 24, ncol = 24)
r = matrix(nrow = 24, ncol = 24)
for (i in 1:24){
  for (k in 1:24) {
    if (cor.test(eval(as.name(names(Data)[i])), eval(as.name(names(Data)[k]))  )$p.value < 0.05) {
   p_value[i,k] = formatC(cor.test(eval(as.name(names(Data)[i])), eval(as.name(names(Data)[k]))  )$p.value, digit = 2) 
   r[i,k] = formatC(cor.test(eval(as.name(names(Data)[i])), eval(as.name(names(Data)[k]))  )$estimate, digit = 2)
    }
  }
}
p_value

#--------------------------------------------------------------
# preprocess
# need normalize the data (is done automatically in LASSO)
# centered and normalized
scale.data = data.frame(scale(Data))
# just normalized:
normalize.data = matrix(nrow = dim(Data)[1], ncol = dim(Data)[2] )
for (i in 1:24) {
  x =  eval(as.name(names(Data)[i])) 
  normalize.data[,i] = x / sd(x)
}
normalize.data = data.frame(normalize.data)
for (i in 1:24){
  names(normalize.data)[i] = names(Data)[i]
}
# check sd and mean
apply(normalize.data,2,sd)

# if no preprocess, use norm.data to save Data (not necessary, just I used norm.data later before I decided to use the un-preprocessed data)
norm.data = Data


#--------------------------------------------------------------
attach(norm.data)
# Best Subset Selection
library(leaps)
regfit.full=regsubsets(Mapk1~.,norm.data) #all subset selection 
# can constrain the number of included variable: regsubsets(...,nvmax = n)
coef(regfit.full,7)  # coefficient of best subset selection with 7 variables
summary(regfit.full) #mat
regfull.summary = summary(regfit.full)
names(regfull.summary) #reg.summary$rss

# Forward Stepwise Selection
regfit.fwd=regsubsets(Mapk1~.,norm.data,method="forward")
coef(regfit.fwd,7)
summary(regfit.fwd)
regfwd.summary = summary(regfit.fwd)

# Backward Stepwise Selection
regfit.bwd=regsubsets(Mapk1~.,norm.data,method="backward")
coef(regfit.bwd,7)
summary(regfit.bwd)
regbwd.summary =  summary(regfit.bwd)

#Compute AIC, BIC, adjusted R^2 in linear models and logistic regression (not used here)
lm1 = lm(Mapk1~.,norm.data) # the full model
AIC(lm1)
BIC(lm1)
summary(lm1)$adj.r.squared


# use Cp and BIC to identify best model
par(mfrow=c(3,2))
n = rep(NA,6)
#plot(reg.summary$rss,xlab="Number of Variables",ylab="RSS",sub = 'best selection', type="l")

#plot(reg.summary$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",sub = 'best selection', type="l")
#which.max(reg.summary$adjr2)
#points(8,reg.summary$adjr2[8], col="red",cex=2,pch=20)

plot(regfull.summary$cp,xlab="Number of Variables",ylab="Cp",sub = 'best selection', type='l')
n[1] = which.min(regfull.summary$cp)
points(n[1],regfull.summary$cp[n[1]],col="red",cex=2,pch=20)

reg_best = regsubsets(Mapk1~., data = norm.data)
best_cp.coef = coef(reg_best, n[1])

plot(regfull.summary$bic,xlab="Number of Variables",ylab="BIC",sub = 'best selection', type='l')
n[2]= which.min(regfull.summary$bic)
points(n[2],regfull.summary$bic[n[2]],col="red",cex=2,pch=20)

reg_best = regsubsets(Mapk1~., data = norm.data)
best_BIC.coef = coef(reg_best, n[2])

plot(regfwd.summary$cp,xlab="Number of Variables",ylab="Cp",sub = 'forward stepwise', type='l')
n[3] = which.min(regfwd.summary$cp)
points(n[3],regfwd.summary$cp[n[3]],col="red",cex=2,pch=20)

reg_best = regsubsets(Mapk1~., data = norm.data, method = 'forward')
forward_Cp.coef = coef(reg_best, n[3])

plot(regfwd.summary$bic,xlab="Number of Variables",ylab="BIC",sub = 'forward stepwise', type='l')
n[4] = which.min(regfwd.summary$bic)
points(n[4],regfwd.summary$bic[n[4]],col="red",cex=2,pch=20)

reg_best = regsubsets(Mapk1~., data = norm.data, method = 'forward')
forward_BIC.coef = coef(reg_best, n[4])

plot(regbwd.summary$cp,xlab="Number of Variables",ylab="Cp",sub = 'backward stepwise', type='l')
n[5] = which.min(regbwd.summary$cp)
points(n[5],regbwd.summary$cp[n[5]],col="red",cex=2,pch=20)

reg_best = regsubsets(Mapk1~., data = norm.data, method = 'backward')
backward_Cp.coef = coef(reg_best, n[5])

plot(regbwd.summary$bic,xlab="Number of Variables",ylab="BIC",sub = 'backward stepwise', type='l')
n[6]= which.min(regbwd.summary$bic)
points(n[6],regbwd.summary$bic[n[6]],col="red",cex=2,pch=20)

reg_best = regsubsets(Mapk1~., data = norm.data, method = 'backward')
backward_BIC.coef = coef(reg_best, n[6])


#--------------------------------------------------------------
# Choosing Among Models by CV
library(ISLR)
library(leaps)
library(dplyr)

# create a predict function for easier use
predict.regsubsets = function(object,newdata,id,...){
  form = as.formula(object$call[[2]]) # Extract the formula used when we called regsubsets()
  mat = model.matrix(form,newdata)    # Build the model matrix
  coefi = coef(object,id=id)          # Extract the coefficiants of the ith model
  xvars = names(coefi)                # Pull out the names of the predictors used in the ith model
  mat[,xvars]%*%coefi               # Make predictions using matrix multiplication
}



# evaluate best subset, forward and backward
k = 10        # number of folds (try 5 too)
par(mfrow=c(2,3))
  set.seed(2) 
  # Assign each observation to a single fold
  folds = sample(rep(1:k, dim(norm.data)[1]/k)) #equally divided
  #folds = sample(1:k, nrow(norm.data), replace = TRUE) #not equally divided
  nvar = 8
  # Create a matrix to store the results of our upcoming calculations
  cv_errors = matrix(NA, k, nvar, dimnames = list(NULL, paste(1:nvar)))
  
  # best subset
  for(j in 1:k){
    # The perform best subset selection on the full dataset, minus the jth fold
    best_fit = regsubsets(Mapk1~., data = norm.data[folds!=j,]) # method can be backward or nothing (best set)
    
    # Inner loop iterates over each size i
    for(i in 1:nvar){
      # Predict the values of the current fold from the "best subset" model on i predictors
      pred = predict(best_fit, norm.data[folds==j,], id=i)
      # Calculate the MSE, store it in the matrix we created above
      cv_errors[j,i] = mean((norm.data$Mapk1[folds==j]-pred)^2)
    }
  }
  
  # Take the mean of over all folds for each model size
  mean_cv_errors = apply(cv_errors, 2, mean)
  # Find the model size with the smallest cross-validation error
  min = which.min(mean_cv_errors)
  # Plot the cross-validation error for each model size, highlight the min
  plot(mean_cv_errors, type='b', main = paste('method: best subset, k=', k))
  points(min, mean_cv_errors[min][1], col = "red", cex = 2, pch = 20)
  
  # after obtaining the best model, refit using the full set to get the more correct coefficients
  reg_best = regsubsets(Mapk1~., data = norm.data)
  print(coef(reg_best, min))


  # forward
  for(j in 1:k){
    # The perform best subset selection on the full dataset, minus the jth fold
    best_fit = regsubsets(Mapk1~., data = norm.data[folds!=j,],method="forward")
    
    # Inner loop iterates over each size i
    for(i in 1:nvar){
      # Predict the values of the current fold from the "best subset" model on i predictors
      pred = predict(best_fit, norm.data[folds==j,], id=i)
      # Calculate the MSE, store it in the matrix we created above
      cv_errors[j,i] = mean((norm.data$Mapk1[folds==j]-pred)^2)
    }
  }
  
  # Take the mean of over all folds for each model size
  mean_cv_errors = apply(cv_errors, 2, mean)
  # Find the model size with the smallest cross-validation error
  min = which.min(mean_cv_errors)
  # Plot the cross-validation error for each model size, highlight the min
  plot(mean_cv_errors, type='b', main = paste('method: forward, k=', k))
  points(min, mean_cv_errors[min][1], col = "red", cex = 2, pch = 20)
  
  # after obtaining the best model, refit using the full set to get the more correct coefficients
  reg_best = regsubsets(Mapk1~., data = norm.data,method="forward")
  # reg_best = regsubsets(Mapk1~., data = norm.data)
  print(coef(reg_best, min))
  
  # backward
  for(j in 1:k){
    # The perform best subset selection on the full dataset, minus the jth fold
    best_fit = regsubsets(Mapk1~., data = norm.data[folds!=j,],method="backward")

    # Inner loop iterates over each size i
    for(i in 1:nvar){
      # Predict the values of the current fold from the "best subset" model on i predictors
      pred = predict(best_fit, norm.data[folds==j,], id=i)
      # Calculate the MSE, store it in the matrix we created above
      cv_errors[j,i] = mean((norm.data$Mapk1[folds==j]-pred)^2)
    }
  }
  
  # Take the mean of over all folds for each model size
  mean_cv_errors = apply(cv_errors, 2, mean)
  # Find the model size with the smallest cross-validation error
  min = which.min(mean_cv_errors)
  # Plot the cross-validation error for each model size, highlight the min
  plot(mean_cv_errors, type='b', main = paste('method: backward, k=', k))
  points(min, mean_cv_errors[min][1], col = "red", cex = 2, pch = 20)
  
  # after obtaining the best model, refit using the full set to get the more correct coefficients
  reg_best = regsubsets(Mapk1~., data = norm.data,method="backward")
  # reg_best = regsubsets(Mapk1~., data = norm.data)
  print(coef(reg_best, min))
  



# test variation of the restuls (set.seed from 1-9)
k = 10        # number of folds (try 5 too)： 10 is more stable
par(mfrow=c(3,3))
for (n in 1:9) {
set.seed(n) 
# Assign each observation to a single fold
folds = sample(rep(1:k, dim(norm.data)[1]/k)) #equally divided
nvar = 8
# Create a matrix to store the results of our upcoming calculations
cv_errors = matrix(NA, k, nvar, dimnames = list(NULL, paste(1:nvar)))

for(j in 1:k){
  # The perform best subset selection on the full dataset, minus the jth fold
  best_fit = regsubsets(Mapk1~., data = norm.data[folds!=j,])
  #best_fit = regsubsets(Mapk1~., data = norm.data[folds!=j,]) # method can be backward or nothing (best set)
  
  # Inner loop iterates over each size i
  for(i in 1:nvar){
    # Predict the values of the current fold from the "best subset" model on i predictors
    pred = predict(best_fit, norm.data[folds==j,], id=i)
    # Calculate the MSE, store it in the matrix we created above
    cv_errors[j,i] = mean((norm.data$Mapk1[folds==j]-pred)^2)
  }
}

# Take the mean of over all folds for each model size
mean_cv_errors = apply(cv_errors, 2, mean)
# Find the model size with the smallest cross-validation error
min = which.min(mean_cv_errors)
# Plot the cross-validation error for each model size, highlight the min
plot(mean_cv_errors, type='b', main = paste('method: bestset, k=', k))
points(min, mean_cv_errors[min][1], col = "red", cex = 2, pch = 20)

# after obtaining the best model, refit using the full set to get the more correct coefficients
reg_best = regsubsets(Mapk1~., data = norm.data)
# reg_best = regsubsets(Mapk1~., data = norm.data)
#backward_CV.coef = coef(reg_best, min)
#forward_CV.coef = coef(reg_best, min)
best_CV.coef = coef(reg_best, min)
}


#--------------------------------------------------------------
# validation set appraoch (not very good for small data set and results not stable - not reported): 
nrep  = 10 #different seed to see whether the results are stable
error.min = rep(NA,nrep)
for (k in 1:nrep) {
set.seed(k)
#train=sample(c(TRUE,FALSE), nrow(norm.data),rep=TRUE)
#you may also use train_index=sample(c(1:40),20) to create training data
train=sample(c(1:40),20)
#test=(!train)
regfit.best=regsubsets(Mapk1~.,norm.data[train,]) #method can also be forward/backward
test.mat=model.matrix(Mapk1~.,norm.data[-train,]) # create an X matrix of test data
val.errors=rep(NA,8)
for(i in 1:8){
  coefi=coef(regfit.best,id=i)
  pred=test.mat[,names(coefi)]%*%coefi
  val.errors[i]=mean((norm.data$Mapk1[-train]-pred)^2)
}
#val.errors
#which.min(val.errors) 
error.min[k] = which.min(val.errors) # use cross-validation: best is the 4th one
}
error.min

#We run the best subset selection with the entire data set and find the model with 4 predictors
regfit.best=regsubsets(Mapk1~.,norm.data)
coef(regfit.best,4)
#not the same as just using training data


 #--------------------------------------------------------------
# Ridge Regression and the Lasso (ridge will not be used)
x=model.matrix(Mapk1~.,norm.data)[,-1]  #remove the intercept
y=norm.data$Mapk1

# Ridge Regression
library(glmnet)
# some examples, no cross-validation
grid=10^seq(10,-2,length=100)   #create a grid for \lambda
ridge.mod=glmnet(x,y,alpha=0,lambda=grid) #alpha=0 is the ridge penalty, alpha=1 is the lasso penalty
dim(coef(ridge.mod))
ridge.mod$lambda[50] # grid[50]
grid[50]
coef(ridge.mod)[,50]
sqrt(sum(coef(ridge.mod)[-1,50]^2))  #calculuate L2 norm of beta
ridge.mod$lambda[60]
grid[60]
coef(ridge.mod)[,60]
sqrt(sum(coef(ridge.mod)[-1,60]^2)) #As we decrease lambda, L_2 norm of beta increase 

predict(ridge.mod,s=50,type="coefficients") #get the coefficient corresponding to lambda=50

# cross-validation
set.seed(1)
grid=10^seq(10,-2,length=100)   #create a grid for \lambda
cv.out=cv.glmnet(x,y,alpha=0, nfolds = 5, lambda = grid) #10 fold cross validation
#cv.out=cv.glmnet(x,y,alpha=0, nfolds = 5)
plot(cv.out) # Draw plot of training MSE as a function of lambda
ridge.mod=glmnet(x,y,alpha=0,lambda=grid) 
plot(ridge.mod)
#coef(ridge.mod)[,1]
bestlam=cv.out$lambda.min
bestlam #83
cv.out$lambda == bestlam
grid == bestlam

# for best lambdata
#ridge.mod=glmnet(x,y,alpha=0,lambda=grid)  # do not change whether lambda is
ridge.pred = predict(ridge.mod, s = bestlam, newx = x) # Use best lambda to predict test data
# use ridge.model to run the prediction; if use cv-out--> slight difference (not the full set)
mean((ridge.pred - y)^2) # Calculate test MSE

# plot test_error vs. log(lambda) not correct --> use plot(cv.out) directly
ridge.mod=glmnet(x,y,alpha=0,lambda=grid) 
te.error = rep(NA,100)
for (i in 1:100) {
  #for each lambda value, get the predict value, calculate test error
  ridge.pred = predict(ridge.mod, s=cv.out$lambda[i], newx = x)
  print(mean((ridge.pred- y)^2))
  te.error[i] = mean((y-ridge.pred)^2)
}
mean((mean(y) - y)^2) #largest MSE: 0.015

par(mfrow=c(3,2))
plot(log(cv.out$lambda), te.error, xlab = 'log(lambda)', ylab = 'test_error',  type = 'l')
plot((cv.out$lambda), te.error, xlab = 'lambda', ylab = 'test_error',  type = 'l')
# why lambda increase, test MSE always increase?

out = glmnet(x, y, alpha = 0, lambda=grid) # Fit ridge regression model on full dataset
predict(out, type = "coefficients", s = bestlam) # Display coefficients using lambda chosen by CV


#--------------------------------------------------------------
# The Lasso
# what can be change? grid (n1, n2, length); seed, nfolds
# bestlam
x=model.matrix(Mapk1~.,norm.data)[,-1]  #remove the intercept
y=norm.data$Mapk1

#how can n1, n2 and length affect test-error, lambda and number of selected variable
# n2 = -1:12; #when l = 1000, must [-1, 15], almost the same between -1 to 12, 13 is a little different (#var = 8)
# folder number is also a variable
par(mfrow=c(3,3))
n1 = 0 # must >= -1, stable between -1~4, better smaller
n2 = -4 #must <=-1, -2 has a little variation, between -8~-3 is stable
ns = 1:20
l = 200 #bigger will be less variable
nfold = 10
fold = 5 : (5+nfold-1)
#n2 = c(100, 200, 300, 400, 500, 600, 700, 800, 900)
var.sd = rep(NA,nfold)
var.med = rep(NA,nfold)
ikk = 1
for (kk in fold) { #kk: nfold
#for (ns in 11:16) {
  te.error = rep(NA,length(n2))
  var.length = rep(NA,length(n2))
  bestlam = rep(NA,length(n2))
  
i = 1
for (n in ns) {
grid=10^seq(n1,n2,length= l)  # n2 cannot be 0
# not doing cross-validation: use the full set
# use in prediction, value of grid does not matter
lasso.mod=glmnet(x,y,alpha=1) 
#plot(lasso.mod) # give the coefficient plot

# use cross-validation to choose the bestlam
set.seed(n) 
cv.out=cv.glmnet(x, y, alpha=1, nfolds = kk) #5 fold cross validation, #control grid
# plot(cv.out) #give the lambda vs. MSE plot
# cv.out=cv.glmnet(x, y, alpha=1, nfolds = kk) # use default lambda
bestlam[i]=cv.out$lambda.min


out = glmnet(x, y, alpha = 1)  # the same as the lasso.mod
lasso.coef = predict(out, type = "coefficients", s = bestlam[i])[1:24,]# Display coefficients using lambda chosen by CV, [1:24,] is just to transpose
#lasso.coef=coef(lasso.mod)[,which(grid == bestlam[i])]  #same results
var.length[i] = length(lasso.coef[lasso.coef!=0]) #number of selected variable
print(lasso.coef[lasso.coef!=0])


pred.lasso = predict(out, s = bestlam[i], newx = x) #not necessary
te.error[i] = mean((y-pred.lasso)^2) #not necessary
# te.error[i] = mean((y-predict(out, s = grid[which(grid == bestlam[i])], newx = x))^2) # same
# te.error always decrease, because as lambda decrease, te.error approximate LSE. It's different from the MSE in CV (is U-shaped)

i = i+1
}
# plot(ns, var.length, xlab = 'ns', ylab = '#variables',  main = paste('nfold_', kk), type = 'l')
var.sd[ikk] = sd(var.length)
var.med[ikk] = median(var.length)
ikk = ikk+1
}

plot(fold, var.sd, xlab = 'nfold', ylab = 'sd of the number of var in model')
plot(fold, var.med, xlab = 'nfold', ylab = 'median of the number of var in model')




grid=10^seq(-4,0,length=200)   # cannot be 0
set.seed(2)
cv.out=cv.glmnet(x,y,alpha=1,nfolds = 10) #5 fold cross validation
plot(cv.out)
lasso.mod = glmnet(x, y, alpha = 1) 
plot(lasso.mod)
bestlam=cv.out$lambda.min
print(bestlam)
#pred.lasso = predict(lasso.mod, s = bestlam, newx = x)
#mean((y-pred.lasso)^2)

out = glmnet(x, y, alpha = 1) # Fit lasso model on full dataset
lasso.coef = predict(out, type = "coefficients", s = bestlam)[1:24,]# Display coefficients using lambda chosen by CV
lasso.coef = lasso.coef[lasso.coef!=0] # Display only non-zero coefficients
# or use lambda = bestlam to build the out model, so it only contains 1 column
#out=glmnet(x,y,alpha=1,lambda=bestlam)
#lasso.coef=coef(out)[,1]
#print(lasso.coef[lasso.coef!=0])

#---------------------------------
#Compare models (not finished):
lasso.coef
best_cp.coef
best_BIC.coef
best_CV.coef
forward_Cp.coef
forward_BIC.coef
forward_CV.coef
backward_Cp.coef
backward_BIC.coef
backward_CV.coef