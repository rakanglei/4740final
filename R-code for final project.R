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
# 2. Best Subset Selection
library(leaps)
# Best subset
regfit.full=regsubsets(Mapk1~.,data) #all subset selection 
summary(regfit.full) #or summary(regfit.full)$outmat: check all the models
regfull.summary = summary(regfit.full)

# Forward Stepwise Selection
regfit.fwd=regsubsets(Mapk1~.,data,method="forward")
summary(regfit.fwd)
regfwd.summary = summary(regfit.fwd)

# Backward Stepwise Selection
regfit.bwd=regsubsets(Mapk1~.,data,method="backward")
summary(regfit.bwd)
regbwd.summary =  summary(regfit.bwd)


# use Cp and BIC to identify best model
par(mfrow=c(3,3))
n = rep(NA,6) #save the best models (number of variables)

# Best-Cp
plot(regfull.summary$cp,xlab="Number of Variables",ylab="Cp",sub = 'best selection', type='l')
n[1] = which.min(regfull.summary$cp)
points(n[1],regfull.summary$cp[n[1]],col="red",cex=2,pch=20)

reg_best = regsubsets(Mapk1~., data = data)
best_cp.coef = coef(reg_best, n[1]) #get the coefficients of the best model

# Forward-Cp
plot(regfwd.summary$cp,xlab="Number of Variables",ylab="Cp",sub = 'forward stepwise', type='l')
n[2] = which.min(regfwd.summary$cp)
points(n[2],regfwd.summary$cp[n[2]],col="red",cex=2,pch=20)

reg_best = regsubsets(Mapk1~., data = data, method = 'forward')
forward_Cp.coef = coef(reg_best, n[2])

# Backwar-Cp
plot(regbwd.summary$cp,xlab="Number of Variables",ylab="Cp",sub = 'backward stepwise', type='l')
n[3] = which.min(regbwd.summary$cp)
points(n[3],regbwd.summary$cp[n[3]],col="red",cex=2,pch=20)

reg_best = regsubsets(Mapk1~., data = data, method = 'backward')
backward_Cp.coef = coef(reg_best, n[3])

# Best-BIC
plot(regfull.summary$bic,xlab="Number of Variables",ylab="BIC",sub = 'best selection', type='l')
n[4]= which.min(regfull.summary$bic)
points(n[4],regfull.summary$bic[n[4]],col="red",cex=2,pch=20)

reg_best = regsubsets(Mapk1~., data = data)
best_BIC.coef = coef(reg_best, n[4])

# Forward-BIC
plot(regfwd.summary$bic,xlab="Number of Variables",ylab="BIC",sub = 'forward stepwise', type='l')
n[5] = which.min(regfwd.summary$bic)
points(n[5],regfwd.summary$bic[n[5]],col="red",cex=2,pch=20)

reg_best = regsubsets(Mapk1~., data = data, method = 'forward')
forward_BIC.coef = coef(reg_best, n[5])

# Backward-BIC
plot(regbwd.summary$bic,xlab="Number of Variables",ylab="BIC",sub = 'backward stepwise', type='l')
n[6]= which.min(regbwd.summary$bic)
points(n[6],regbwd.summary$bic[n[6]],col="red",cex=2,pch=20)

reg_best = regsubsets(Mapk1~., data = data, method = 'backward')
backward_BIC.coef = coef(reg_best, n[6])

#---------------------------------------
# use CV to identify best model
library(leaps)
library(dplyr)

# create a predict function for easier use
predict.regsubsets = function(object,newdata,id,...){
  form = as.formula(object$call[[2]]) # Extract the formula used when we called regsubsets()
  mat = model.matrix(form,newdata)    # Build the model matrix
  coefi = coef(object,id=id)          # Extract the coefficiants of the ith model
  xvars = names(coefi)                # Pull out the names of the predictors used in the ith model
  mat[,xvars]%*%coefi                 # Make predictions using matrix multiplication
}

# evaluate best subset, forward and backward
k = 10        
par(mfrow=c(2,3))
set.seed(2) 
# Assign each observation to a single fold
folds = sample(rep(1:k, dim(data)[1]/k)) #equally divided
nvar = 8 #number of models
# Create a matrix to store the results of our upcoming calculations
cv_errors = matrix(NA, k, nvar, dimnames = list(NULL, paste(1:nvar)))

# best subset
for(j in 1:k){
  # The perform best subset selection on the full dataset, minus the jth fold
  best_fit = regsubsets(Mapk1~., data = data[folds!=j,]) 
  # Inner loop iterates over each size i
  for(i in 1:nvar){
    # Predict the values of the current fold from the "best subset" model on i predictors
    pred = predict(best_fit, data[folds==j,], id=i)
    # Calculate the MSE, store it in the matrix we created above
    cv_errors[j,i] = mean((data$Mapk1[folds==j]-pred)^2)
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
reg_best = regsubsets(Mapk1~., data = data)
print(coef(reg_best, min))


# forward
for(j in 1:k){
  # The perform best subset selection on the full dataset, minus the jth fold
  best_fit = regsubsets(Mapk1~., data = data[folds!=j,],method="forward")
  # Inner loop iterates over each size i
  for(i in 1:nvar){
    # Predict the values of the current fold from the "best subset" model on i predictors
    pred = predict(best_fit, data[folds==j,], id=i)
    # Calculate the MSE, store it in the matrix we created above
    cv_errors[j,i] = mean((data$Mapk1[folds==j]-pred)^2)
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
reg_best = regsubsets(Mapk1~., data = data,method="forward")
# reg_best = regsubsets(Mapk1~., data = data)
print(coef(reg_best, min))

# backward
for(j in 1:k){
  # The perform best subset selection on the full dataset, minus the jth fold
  best_fit = regsubsets(Mapk1~., data = data[folds!=j,],method="backward")
  # Inner loop iterates over each size i
  for(i in 1:nvar){
    # Predict the values of the current fold from the "best subset" model on i predictors
    pred = predict(best_fit, data[folds==j,], id=i)
    # Calculate the MSE, store it in the matrix we created above
    cv_errors[j,i] = mean((data$Mapk1[folds==j]-pred)^2)
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
reg_best = regsubsets(Mapk1~., data = data,method="backward")
# reg_best = regsubsets(Mapk1~., data = data)
print(coef(reg_best, min))



# test variation of the restuls (set.seed from 1-9): need to hand change the method
k = 10        # number of folds (try 5 too)： 10 provides less variant results
par(mfrow=c(3,3))
for (n in 1:9) {
  set.seed(n) 
  # Assign each observation to a single fold
  folds = sample(rep(1:k, dim(data)[1]/k)) #equally divided
  nvar = 8
  # Create a matrix to store the results of our upcoming calculations
  cv_errors = matrix(NA, k, nvar, dimnames = list(NULL, paste(1:nvar)))
  
  for(j in 1:k){
    # The perform best subset selection on the full dataset, minus the jth fold
    best_fit = regsubsets(Mapk1~., data = data[folds!=j,])
    #best_fit = regsubsets(Mapk1~., data = data[folds!=j,], method = 'forward')
    #best_fit = regsubsets(Mapk1~., data = data[folds!=j,], method = 'backward') 
    # Inner loop iterates over each size i
    for(i in 1:nvar){
      # Predict the values of the current fold from the "best subset" model on i predictors
      pred = predict(best_fit, data[folds==j,], id=i)
      # Calculate the MSE, store it in the matrix we created above
      cv_errors[j,i] = mean((data$Mapk1[folds==j]-pred)^2)
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
  reg_best = regsubsets(Mapk1~., data = data)
  # reg_best = regsubsets(Mapk1~., data = data, method = 'forward')
  # reg_best = regsubsets(Mapk1~., data = data, method = 'backward')

  best_CV.coef = coef(reg_best, min)
  #forward_CV.coef = coef(reg_best, min)
  #backward_CV.coef = coef(reg_best, min)
}



#--------------------------------------------------------------
# Ridge Regression and the Lasso (ridge will not be used)
x=model.matrix(Mapk1~.,data)[,-1]  #remove the intercept
y=data$Mapk1

# Ridge Regression
library(glmnet)
# cross-validation
set.seed(1)
cv.out=cv.glmnet(x,y,alpha=0, nfolds = 10) #10 fold cross validation
plot(cv.out) # Draw plot of training MSE as a function of lambda
ridge.mod=glmnet(x,y,alpha=0) 
plot(ridge.mod)
bestlam=cv.out$lambda.min

out = glmnet(x, y, alpha = 0) # Fit ridge regression model on full dataset
coef = predict(out, type = "coefficients", s = bestlam)[1:24,] # Display coefficients using lambda chosen by CV
coef 


#--------------------------------------------------------------
# The Lasso
# what can be change? grid (n1, n2, length) --> use the default value; seed, nfolds
# bestlam
library(glmnet)
x=model.matrix(Mapk1~.,data)[,-1]  #remove the intercept
y=data$Mapk1

par(mfrow=c(3,3))
ns = 1:20
nfold = 10
fold = 5 : (5+nfold-1)
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
    print(lasso.coef[lasso.coef!=0])
    i = i+1
  }
  var.sd[ikk] = sd(var.length)
  var.med[ikk] = median(var.length)
  ikk = ikk+1
}

plot(fold, var.sd, xlab = 'nfold', ylab = 'sd of the number of var in model')
plot(fold, var.med, xlab = 'nfold', ylab = 'median of the number of var in model')



# use just one seed
set.seed(1)
cv.out=cv.glmnet(x,y,alpha=1,nfolds = 10) 
plot(cv.out)
lasso.mod = glmnet(x, y, alpha = 1) 
plot(lasso.mod)
bestlam=cv.out$lambda.min
print(bestlam)

out = glmnet(x, y, alpha = 1) # Fit lasso model on full dataset
lasso.coef = predict(out, type = "coefficients", s = bestlam)[1:24,]# Display coefficients using lambda chosen by CV
lasso.coef = lasso.coef[lasso.coef!=0] # Display only non-zero coefficients





#---------------------------------
#Compare models (not finished):
#lasso.coef
#best_cp.coef
#best_BIC.coef
#best_CV.coef
#forward_Cp.coef
#forward_BIC.coef
#forward_CV.coef
#backward_Cp.coef
#backward_BIC.coef
#backward_CV.coef


library(boot)
# vary nfold and set.seed
ns = 1:9
nfold=c(5, 10, 15, 20)
par(mfrow=c(3,3))
mins = matrix(NA,length(nfold),length(ns))
mins.med = rep(NA,length(nfold))
mins.sd = rep(NA,length(nfold))

kk = 1
for (k in nfold) { 
  cv.error = matrix(NA,length(ns),6)
  
  i = 1
  for (n in ns) {
    set.seed(n)
    glm.fit = glm(Mapk1 ~ Pik3r3 + Rac1, data = data) #best_CV
    cv.error[i,1] = cv.glm(data, glm.fit, K = k)$delta[1]
    
    set.seed(n)
    glm.fit = glm(Mapk1 ~ Ppp3cb + Rac1 , data = data) #back_CV
    cv.error[i,2] = cv.glm(data, glm.fit, K = k)$delta[1]
    
    set.seed(n)
    glm.fit = glm(Mapk1 ~ Akt2 + Rik + Pik3r3 + Rac1, data = data) #best_BIC
    cv.error[i,3] = cv.glm(data, glm.fit, K = k)$delta[1]
    
    
    set.seed(n)
    glm.fit = glm(Mapk1 ~ Akt2 + Rik + Pik3r3 + Pik3r1 + Rac1, data = data) #best_cp
    cv.error[i,1] = cv.glm(data, glm.fit, K = k)$delta[1]
    
    set.seed(n)
    glm.fit = glm(Mapk1 ~ Rik + Pik3cd + Pik3r3 + Rac1 + Nfat5, data = data) #lasso
    cv.error[i,5] = cv.glm(data, glm.fit, K = k)$delta[1]
    
    set.seed(n)
    glm.fit = glm(Mapk1 ~ Cdc42 + Akt2 + Plcg2 + Rac2 + Sphk2 + Ppp3cb + Rac1, data = data) #back_Cp/BIC
    cv.error[i,6] = cv.glm(data, glm.fit, K = k)$delta[1]
    
    set.seed(n)
    glm.fit = glm(Mapk1 ~ Rik + Pik3r3 + Pik3r1 + Rac1, data = data) #best_cp
    #cv.error[i,7] = cv.glm(data, glm.fit, K = k)$delta[1]
    
    set.seed(n)
    glm.fit = glm(Mapk1 ~ Akt2 + Rik + Pik3r3 + Rac1, data = data) #best_cp
    #cv.error[i,8] = cv.glm(data, glm.fit, K = k)$delta[1]
    
    set.seed(n)
    glm.fit = glm(Mapk1 ~ Rik + Pik3r3 + Rac1, data = data) #best_cp
    #cv.error[i,9] = cv.glm(data, glm.fit, K = k)$delta[1]
    
    set.seed(n)
    glm.fit = glm(Mapk1 ~ Akt2 + Rik + Pik3r3 + Rac1 + Nfat5, data = data) #best_cp
    #cv.error[i,10] = cv.glm(data, glm.fit, K = k)$delta[1]
    
    set.seed(n)
    glm.fit = glm(Mapk1 ~ Akt2 + Rik + Pik3r1 + Pik3r3 + Rac1 + Nfat5, data = data) #best_cp
    #cv.error[i,11] = cv.glm(data, glm.fit, K = k)$delta[1]
    
    plot(cv.error[i,], xlab = 'ind', ylab = 'MSE')
    min = which.min(cv.error[i,])
    points(min, cv.error[i,min][1], col = "red", cex = 2, pch = 20)
    
    i = i + 1
  }
  
  #cv.error
  mins[kk,] = apply(cv.error,1,which.min)
  mins.med[kk] = median(mins[kk,])
  mins.sd[kk] = sd(mins[kk,])
  kk = kk + 1
}


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


