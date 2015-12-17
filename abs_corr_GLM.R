##MAKE A PREDICTIVE MODEL OF LOCUS PRESENCE-ABSENCE IN PAIRWISE COMPARISONS
#Load in RF and GLMNET
library(glmnet)
library(randomForest)
library(bartMachine)

##Dataset created in 'train_builder.R' from pairwise shared matrix derived from ABS-PRES matrix in previous scripts (abs_pres_loci.py & abs_corr.R)
#Read in X and Y
X = read.table(file="~/Desktop/TRAIN",header=TRUE)
X = subset(X,select = -c(sample1,sample2))
PS = as.vector(t(read.table(file="~/Desktop/pairwise_shared/all_10k/PS_all")))
Y = PS

#Create training and test sets
n_total = length(Y)
train_sample = sample(1:n_total, size = 32000, replace = FALSE)
X_train = X[train_sample,]
X_test = X[-train_sample,]
Y_train = Y[train_sample]
Y_test = Y[-train_sample]

#Correlation matrix for understanding variable relationships
CORmatrix = cor(X)
plot(density(as.vector(CORmatrix),bw=.01),main="Variable Relationship: Density Distribution of Cor Matrix",xlab="Correlation")

##############################
##BUILD GLM WITH ELASTIC-NET##
##############################
#Turn X_train and X_test into matrices for GLMNET
X_train_mat = as.matrix(X_train)
X_test_mat = as.matrix(X_test)

#For variable interpretation
mat_names = colnames(X_train_mat)

#Build GLMs with alpha values (1, 0.95, 0.5, 0) LASSO to Ridge regression
GLM1 = glmnet(X_train_mat,Y_train,family="gaussian")
GLM2 = glmnet(X_train_mat,Y_train,alpha=0.95,family="gaussian")
GLM3 = glmnet(X_train_mat,Y_train,alpha=0.5,family="gaussian")
GLM4 = glmnet(X_train_mat,Y_train,alpha=0,family="gaussian")

#Let's plot them against lambda
plot(GLM1,xvar="lambda",label=TRUE)
plot(GLM2,xvar="lambda",label=TRUE)
plot(GLM3,xvar="lambda",label=TRUE)
plot(GLM4,xvar="lambda",label=TRUE)

#Build Cross-validated GLMs with alpha values (1, 0.95, 0.5, 0) LASSO to Ridge regression
GLMCV1 = cv.glmnet(X_train_mat,Y_train,family="gaussian")
GLMCV2 = cv.glmnet(X_train_mat,Y_train,alpha=0.95,family="gaussian")
GLMCV3 = cv.glmnet(X_train_mat,Y_train,alpha=0.5,family="gaussian")
GLMCV4 = cv.glmnet(X_train_mat,Y_train,alpha=0,family="gaussian")

#Let's plot them against lambda
plot(GLMCV1)
plot(GLMCV2)
plot(GLMCV3)
plot(GLMCV4)

#Let's use lambda.min and lambda.1se to predict
P1 = predict(GLMCV1, X_test_mat, s = "lambda.1se")
P2 = predict(GLMCV2, X_test_mat, s = "lambda.1se")
P3 = predict(GLMCV3, X_test_mat, s = "lambda.1se")
P4 = predict(GLMCV4, X_test_mat, s = "lambda.1se")
P5 = predict(GLMCV1, X_test_mat, s = "lambda.min")
P6 = predict(GLMCV2, X_test_mat, s = "lambda.min")
P7 = predict(GLMCV3, X_test_mat, s = "lambda.min")
P8 = predict(GLMCV4, X_test_mat, s = "lambda.min")

#Let's compare predicted vs. expected
plot(Y_test,P1)
abline(a = 0, b = 1, lty=2)
plot(Y_test,P2)
abline(a = 0, b = 1, lty=2)
plot(Y_test,P3)
abline(a = 0, b = 1, lty=2)
plot(Y_test,P4)
abline(a = 0, b = 1, lty=2)
plot(Y_test,P5,xlab="Observed Shared Loci",ylab="Expected Shared Loci (LASSO)")
abline(a = 0, b = 1, lty=2)
plot(Y_test,P6)
abline(a = 0, b = 1, lty=2)
plot(Y_test,P7)
abline(a = 0, b = 1, lty=2)
plot(Y_test,P8,xlab="Observed Shared Loci",ylab="Expected Shared Loci (RIDGE)")
abline(a = 0, b = 1, lty=2)

#Deviations
dev1 = Y_test - P1
dev2 = Y_test - P2
dev3 = Y_test - P3
dev4 = Y_test - P4
dev5 = Y_test - P5
dev6 = Y_test - P6
dev7 = Y_test - P7
dev8 = Y_test - P8
DEV = cbind(dev1,dev2,dev3,dev4,dev5,dev6,dev7,dev8)

#MSE and SD of devs
D1m = mean(dev1^2)
D2m = mean(dev2^2)
D3m = mean(dev3^2)
D4m = mean(dev4^2)
D5m = mean(dev5^2)
D6m = mean(dev6^2)
D7m = mean(dev7^2)
D8m = mean(dev8^2)
DEV_sm = c(D1m,D2m,D3m,D4m,D5m,D6m,D7m,D8m)

D1sd = sd(dev1)
D2sd = sd(dev2)
D3sd = sd(dev3)
D4sd = sd(dev4)
D5sd = sd(dev5)
D6sd = sd(dev6)
D7sd = sd(dev7)
D8sd = sd(dev8)
DEV_sd = c(D1sd,D2sd,D3sd,D4sd,D5sd,D6sd,D7sd,D8sd)

#Matching training (min & 1se) and test MSE
MSE1 = c(GLMCV1$cvm[match(GLMCV1$lambda.min,GLMCV1$lambda)],GLMCV1$cvm[match(GLMCV1$lambda.1se,GLMCV1$lambda)])
MSE2 = c(GLMCV2$cvm[match(GLMCV2$lambda.min,GLMCV2$lambda)],GLMCV2$cvm[match(GLMCV2$lambda.1se,GLMCV2$lambda)])
MSE3 = c(GLMCV3$cvm[match(GLMCV3$lambda.min,GLMCV3$lambda)],GLMCV3$cvm[match(GLMCV3$lambda.1se,GLMCV3$lambda)])
MSE4 = c(GLMCV4$cvm[match(GLMCV4$lambda.min,GLMCV4$lambda)],GLMCV4$cvm[match(GLMCV4$lambda.1se,GLMCV4$lambda)])

#Plotting training vs. test set MSE
#Checking for overfitting
plot(1:8,DEV_sm)
points(1:8,c(MSE1[2],MSE2[2],MSE3[2],MSE4[2],MSE1[1],MSE2[1],MSE3[1],MSE4[1]),col=2)

#######################
##BUILD RANDOM FOREST##
#######################
#Build random forest
RF1 = randomForest(Y_train ~ ., X_train,importance=TRUE,ntree=100)

#Check out Random Forest
plot(RF1,main="Random Forest: OOB Error")
varImpPlot(RF1,main="Random Forest: Variable Importance")
P9 = predict(RF1,X_test,type="response")

#Build PDM test set
PDM_n = 1000
PDMs = seq(0,0.05,5e-05)
PDM_L = length(PDMs)
X_PDM = as.data.frame(PDMs)
for (i in 2:length(colnames(X_test))){
    X_PDM = cbind(X_PDM,rep(mean(X_test[,i]),PDM_L))
}
colnames(X_PDM) = colnames(X_test)  

#Predict PDM test set
P_pdm = predict(RF1,X_PDM,type="response")
plot(PDMs,P_pdm,xlab="Phylogenetic Distance", ylab = "Expected Pairwise Shared Loci")

#Plot and Deviations
plot(Y_test,P9,xlab="Observed Shared Loci",ylab="Expected Shared Loci (Random Forest)", sub = "r = 0.99")
plot(log(Y_test),log(P9),xlab="Log Observed Shared Loci",ylab="Log Expected Shared Loci (Random Forest)", sub = "r = 0.99")
abline(a = 0, b = 1, lty=2)
dev9 = Y_test - P9
D9m = mean(dev9^2)

####################
##BUILD BART MODEL##
####################
#Build BART model
bart1 = bartMachine(X = X_train, y = Y_train, num_trees=5, run_in_sample = TRUE)

#Build BART CV model


#Predict with BART
P10 = predict(bart1,X_test)

#Deviations from P10
dev10 = Y_test - P10
D10m = mean(dev10^2)
D10sd = sd(dev10)

#Plot 
plot(Y_test,P10,xlab="Observed Shared Loci",ylab="Expected Shared Loci (BART)", sub = "r = 0.95; MSE = 1.1e7")
abline(a=0,b=1,lty=2)

#####################
##METHOD COMPARISON##
#####################
#Add RF to method comparision
DEV = cbind(DEV,dev9)
DEV_sma = c(DEV_sm,D9m)

#Plot methods (8 GLM, 1 RF)
plot(1:9,DEV_sma,pch=1:9,xlab="Model Number",ylab="Sum of Squared Errors",sub="Lasso (1,5); Elastic Net (2,6); Elastic Net (3,7); Ridge (4,8); Random Forest (9)", main="Model Comparison: Training Set Performance")

#Plot deviations and their density distributions
plot(density(dev9),main="Model Comparison: Deviation Density Distributions on Test Set")
lines(density(dev1),col=2)
lines(density(dev2),col=3)
lines(density(dev3),col=4)
lines(density(dev4),col=5)
lines(density(dev5),col=6)
lines(density(dev6),col=7)
lines(density(dev7),col=8)
lines(density(dev8),col=9)
