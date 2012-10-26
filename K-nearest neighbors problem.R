

options(width=65)
options(stringsAsFactors = FALSE)
library(gregmisc)
library(gmodels)
library(qvalue)
library(samr)
library(Biobase)
library(limma)

####  Problem 4  

# Download the data provided on Blackboard (dataHW3-breastcancer.zip), unzip
# the file and use \load" to read it into R. This will load the \breastdata" object,
# which is a dataset for 46 breast tumor samples where 23 are positive for an
# estrogen receptor (ER+) and 23 were negative (ER-) (West et al., PNAS 2001
# 98:11462-11467). There are expression levels for 7129 genes for each sample in
# this list (x) and class labels for each sample (y).

load("../Data/dataHW3-breastcancer.rdata")

dim(newpatients)   # newpatients is a matrix of 7129 expression values for  3 patients
length(trueclasses) #  trueclasses is a vector of classes for the 3 pts "ER+" "ER+" "ER-"  

str(breastcancer)   #  consists of an expression matrix x, a vector of genenames and a vector of class labels
class(breastcancer$x)
dim(breastcancer$x)   #   7129   46
head(breastcancer$x)
length(breastcancer$y)    #  46  consisting of ER+ or ER-  class labels

head(breastcancer$genenames)
head(breastcancer$geneid)   #  Affymentrix IDs?

#  Part a

#Code a function for k-nearest neighbor classification with an option for either
#Euclidean or correlation distance. Write your own function and please turn in all
#your code.
#############################################################################

#  this first section is "resubstitution"  we are testing samples on the same 
#  dataset used for classification.  This is generally bad practice as it is likely 
#  to be overfitted and perform poorly when used on new datasets

source("make_dist.r")
source("nw_knn.r")

KNN.class <- nw_knn(x = breastcancer$x, true.class = breastcancer$y, k = 7)


##  Part B
# Apply your function to the data and plot the error rate for different values of k
#  and both distance metrics.

#  So KNN.class is a vector of classifications based on KNN
#  compare to true classification
cbind(KNN.class, breastcancer$y)
CrossTable(KNN.class, breastcancer$y)
#  proportion misclassified 

#  Plot your results for different value of k and both distance metrics.
sum(KNN.class != breastcancer$y)/46   #   0.1521739 with k=2  euclidean
#   0.1304348 with k=3  euclidean
#   0.1304348 with k=4  euclidean
#   0.1304348 with k=5  euclidean
#   0.1521739 with k=6  euclidean
#   0.1304348 with k=7  euclidean
#   0.1521739 with k=8  euclidea


##   Part C
# Perform 5-fold cross validation to determine your error rate on 46 patients.
" Since 46/5 is not an integer, create 4 sets of 10 subjects and 1 set of 6 subjects. 
# Report the average error rate for the 4 sets of 10 subjects, after training on the 
# rest of the subjects".

#  setting up random groups of samples - 
Numbs <- c(1:46)
GroupA <- sample(Numbs, 10, replace=FALSE)
GroupNonA <- Numbs[-GroupA]
GroupB <- sample(Numbs[-GroupA],10, replace=FALSE)
GroupC <- sample(Numbs[-c(GroupA, GroupB)], 10, replace=FALSE)
GroupD <- sample(Numbs[-c(GroupA, GroupB, GroupC)],10,replace=FALSE)
GroupE <- sample(Numbs[-c(GroupA, GroupB, GroupC, GroupD)])

GROUPS <- rbind(GroupA, GroupB, GroupC, GroupD)

#  make a 
k <- 5
LL <- NULL
MM <- NULL
NN <- NULL
OO <- NULL
PP <- NULL
QQ <- NULL
RR <- NULL
SS <- NULL
PredictClass <- NULL
Prediction.vector <- NULL

for (j in 1:4){
  
  for (i in 1:10) {
    #  make a distance matrix with the first element of group A and all remaining samples
    #  This term indexes all the samples not in Group A: Numbs[-GROUPS[1, ]]
    LL <- make.dist(breastcancer$x[ ,c(GROUPS[j, i], Numbs[-GROUPS[j, ]])], "euclidean")
    #  select the first row of the new distance matrix - the "test" sample will always be the first
    #  Next select the k neighbors with the lowest distance from the test sample
    #  order the distnances of the other samples from the test sample and choose the k lowest
    MM <- LL[1, order(LL[1, ])[2:(k+1)]]
    #  some housekeeping
    NN <- labels(MM)  #  this gives a series of 3 labels "V3" or "V19"        
    OO <- substring(NN, 2,3)  #  this strips off the "V", leaving only a number 
    #  (but in character format)
    PP <- as.integer(OO)   #  this convert the character value to an integer 
    # so that we can use it to extract class info      
    RR <- breastcancer$y[PP]   # this gives the class labels for the k nearest neighbors, 
    # ie [1] "ER+" "ER+" "ER+"
    SS <- sum(RR == "ER+")     #  determines many of the nearest neighbors are ER+
    PredictClass[SS/k >.5] <- "ER+"     #  if >50% of nearest neighbors are ER+ makes T "ER+"
    PredictClass[SS/k <.5] <- "ER-"      #  if <50% of nearest neighbors are ER+ makes T "ER-"
    
    #  is this class label correct or incorrect?
    TrueClass <- breastcancer$y[GroupA[i]]
    Prediction <- TrueClass == PredictClass
    Prediction.vector <- c(Prediction.vector, Prediction)
  }
}

#  prediction accuracy
sum(Prediction.vector)/length(Prediction.vector)


##  Part D
# Select the best options for the value of k and the distance metric. Predict the
# tumor class for three new patients in \newpatients". The true classes are in
# \trueclasses". How well did you do?





#  from Coombs and Baggerly
library('class');
ds1 <- matrix(rnorm(40),20,2) + 1;
dim(ds1)
ds2 <- matrix(rnorm(40),20,2) + 2;
dim(ds1)

Train <- rbind(ds1[1:10,],ds2[1:10,])
Test <- rbind(ds1[11:20,],ds2[11:20,])
Class <- as.factor(c(rep(1,20),rep(2,20)))

my.knn <- knn(
  train = Train,
  test =  Test,
  k = 1,
  cl = Class)


#  from Class package

#  they are drawing these example data from iris3 which is an "array" object ?

train <- rbind(iris3[1:25,,1], iris3[1:25,,2], iris3[1:25,,3])
dim(train)   #  75  4
class(train)   #  a matrix
colnames(train)  #   "Sepal L." "Sepal W." "Petal L." "Petal W."


test <- rbind(iris3[26:50,,1], iris3[26:50,,2], iris3[26:50,,3])
dim(test)    #75  4
colnames(test)  #   "Sepal L." "Sepal W." "Petal L." "Petal W."

class <- factor(c(rep("s",25), rep("c",25), rep("v",25)))
class    #  these are the TRUE classifications of the training set

knn(train = train, 
    test = test, 
    cl = class, 
    k = 3, 
    prob=TRUE)   #  since prob = TRUE this gives us the proportion of votes 
#  for the winning class  

knn(train = train, 
    test = test, 
    cl = class, 
    k = 3, 
    prob=FALSE)


#  LOOK AT COOMBS AND BAGGERLY FOR LEAVE ONE OUT CROSS VALIDATION!!
