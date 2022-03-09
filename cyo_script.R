#Script:Exoplanet Hunting in Deep Space
#S. Kube 08/03/2022
# Description: this script performs all analysis and creates all figures as descripied in the report

#Libraries

if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(dplyr)) install.packages("dplyr", repos = "http://cran.us.r-project.org")
if(!require(readr)) install.packages("readr", repos = "http://cran.us.r-project.org")
if(!require(pracma)) install.packages("pracma", repos = "http://cran.us.r-project.org")
if(!require(matrixStats)) install.packages("matrixStats", repos = "http://cran.us.r-project.org")
if(!require(e1071)) install.packages("e1071", repos = "http://cran.us.r-project.org")
if(!require(googledrive)) install.packages("googledrive", repos = "http://cran.us.r-project.org")


library(tidyverse)
library(caret)
library(dplyr)
library(readr)
library(pracma)
library(matrixStats)
library(e1071)
library(googledrive)

#download data from google drive to current working directory
#caveat: google account is required
#original data is from kaggle 
#https://www.kaggle.com/keplersmachines/kepler-labelled-time-series-data/activity
#if download from google drive does not work, 
#download from kaggle manually and place them in your working directory
#then remove comments from next to lines
#train_data <- read_csv("exoTrain.csv")
#validation_data <- read_csv(exoTest.csv)
#then comment out the following section
# download validation and training data

tempf <- tempfile(fileext = ".zip")

#access google drive via API
dl <- drive_download(as_id("1FJ4aSozjirtxqAcw310w_tSGOeYeYNYB"), path = tempf, overwrite = TRUE)

#unzip archive
out <- unzip(tempf, exdir = tempdir())

#load data
train_data <- read_csv(out[3])

validation_data <- read_csv(out[1])
#section end
#data summary

# number of stars with exoplanets in the training data set
sum(train_data$LABEL == 2)

# probability of picking a star with exoplanets at random from training data
sum(train_data$LABEL == 2)/nrow(train_data)

# check if any "NA"s or "-Inf"s are in the data 
na_inf_check <- apply(train_data, 2, function(x) any(is.na(x) | is.infinite(x)))
sum(na_inf_check == TRUE)

# initial data wrangling

# add ID column and time interval column and convert wide to tidy

train_tidy <- train_data %>% 
  cbind(id=as.numeric(rownames(.)),.) %>% # add ID
  gather(time,flux, 'FLUX.1':'FLUX.3197') %>% # convert to tidy
  mutate(t=as.numeric(str_extract(time,"\\d{1,4}"))) %>% # extract numeric time interval t 
  select(id,LABEL,t,flux) # throw away column names

train_tidy_norm <- train_tidy %>%
  mutate_at(c("flux"), ~(scale(.) %>% as.vector))


# do the same to validation_data

validation_tidy <- validation_data %>% 
  cbind(id=rownames(.),.) %>% # add ID
  gather(time,flux, 'FLUX.1':'FLUX.3197') %>% # convert to tidy
  mutate(t=as.numeric(str_extract(time,"\\d{1,4}"))) %>% # extract numeric time interval t 
  select(id,LABEL,t,flux) # throw away column names

validation_tidy_norm <- validation_tidy %>%
  mutate_at(c("flux"), ~(scale(.) %>% as.vector))

# select three stars w/ exoplanets and three w/o exoplanets from train data as examples

# select ids of stars w/ exoplanets
ep <- train_data %>% 
  cbind(id=rownames(.),.) %>%
  filter(LABEL==2) %>%
  select(id)

# random select 3 of them

ep_3 <- sample(1:length(ep[,1]), 3)

# select ids of stars w/o exoplanets
no_ep <- train_data %>% 
  cbind(id=rownames(.),.) %>%
  filter(LABEL==1) %>%
  select(id)

# random select 3 of them

no_ep_3 <- sample(1:length(no_ep[,1]), 3)

# create example plots

# w/ exoplanets
fig_1 <- train_tidy_norm %>%
  filter(id %in% ep_3) %>%
  ggplot(.,aes(t,flux)) + 
  geom_line() +
  facet_grid(rows=vars(id)) + 
  theme_bw()

fig_2 <- train_tidy %>%
  filter(id %in% ep_3) %>%
  ggplot(.,aes(t,flux)) + 
  geom_line() +
  facet_grid(rows=vars(id)) + 
  theme_bw()

#w/o exoplanet
fig_3 <- train_tidy_norm %>%
  filter(id %in% no_ep_3) %>%
  ggplot(.,aes(t,flux)) + 
  geom_line() +
  facet_grid(rows=vars(id)) + 
  theme_bw()

fig_4 <- train_tidy %>%
  filter(id %in% no_ep_3) %>%
  ggplot(.,aes(t,flux)) + 
  geom_line() +
  facet_grid(rows=vars(id)) + 
  theme_bw()

#print figures
fig_1

fig_2

fig_3

fig_4

#pca of raw fluxes
# perform pca (normalized and center data in function, first 50 principal components)
pca_raw_flux <- prcomp(train_data[,2:3198],scale=TRUE,center=TRUE,rank.=50)
summary(pca_raw_flux)

# create scatter plots of first 10 principal components,
#limit data to first 1000 rows (planet w/ exoplanets are in the first 37 rows)
fig_5 <- data.frame(PC1 = pca_raw_flux$x[1:1000,1],
                    PC2 = pca_raw_flux$x[1:1000,2],
                    label=train_data[1:1000,1]) %>%
  ggplot(aes(PC1, PC2, fill=as.factor(LABEL)))+
  geom_point(cex=3, pch=21)

fig_6 <- data.frame(PC3 = pca_raw_flux$x[1:1000,3],
                    PC4 = pca_raw_flux$x[1:1000,4],
                    label=train_data[1:1000,1]) %>%
  ggplot(aes(PC3, PC4, fill=as.factor(LABEL)))+
  geom_point(cex=3, pch=21)

fig_7 <- data.frame(PC5 = pca_raw_flux$x[1:1000,5],
                    PC6 = pca_raw_flux$x[1:1000,6],
                    label=train_data[1:1000,1]) %>%
  ggplot(aes(PC5, PC6, fill=as.factor(LABEL)))+
  geom_point(cex=3, pch=21)

fig_8 <- data.frame(PC7 = pca_raw_flux$x[1:1000,7],
                    PC8 = pca_raw_flux$x[1:1000,8],
                    label=train_data[1:1000,1]) %>%
  ggplot(aes(PC7, PC8, fill=as.factor(LABEL)))+
  geom_point(cex=3, pch=21)

fig_9 <- data.frame(PC9 = pca_raw_flux$x[1:1000,9],
                    PC10 = pca_raw_flux$x[1:1000,10],
                    label=train_data[1:1000,1]) %>%
  ggplot(aes(PC9, PC10, fill=as.factor(LABEL)))+
  geom_point(cex=3, pch=21)

#print figures
fig_5

fig_6

fig_7

fig_8

fig_9

#calculation of periodograms
# calculation of periodograms
train_ps <- matrix(ncol = 1599, nrow = 0) #create empty ps matrix


for (i in 1:5087) { # extract time series row by row
  temp<-as.matrix(train_data[i,2:3197])
  ints2 <- abs(fft(temp))^2/3196 # calculate intensities from squared amplitudes
  scaled_ints <- (4/3196)*ints2[1:1599] # re-scale
  train_ps <- rbind(train_ps,scaled_ints) # create ps matrix
}

# rename rows and columns
rownames(train_ps)<- seq(1,5087,1)
colnames(train_ps)<- seq(1,1599,1)

# min-max scaling of each individual periodogram
train_ps_centered <- sweep(train_ps, 1, apply(train_ps,1,min))
train_ps_standardized <- sweep(train_ps_centered, 1, apply(train_ps,1,max), FUN = "/")

# add ID column, labels, and frequency column and convert wide to tidy

train_ps_tidy <- train_ps_standardized %>%
  as_tibble() %>%
  cbind(select(train_data,LABEL),.) %>% # add labels
  cbind(id=as.numeric(rownames(.)),.) %>% # add ID
  gather(key,ints, "1":"1599") %>% # convert to tidy
  mutate(f=(as.numeric(key)-1)/3196) %>% # add frequency
  select(id,LABEL,f,ints) # throw away column names

# extract periodograms for examples in fig. 1 and fig. 3
# w/ exoplanets
fig_10 <- train_ps_tidy %>%
  filter(id %in% ep_3) %>%
  ggplot(.,aes(f,ints)) + 
  geom_line() +
  facet_grid(rows=vars(id)) + 
  labs(x="frequency",y="intensity") +
  theme_bw()

# scale x to log10 to emphasize low-frequencies in visualization

fig_11 <- train_ps_tidy %>%
  filter(id %in% ep_3) %>%
  ggplot(.,aes(f,ints)) + 
  geom_line() +
  facet_grid(rows=vars(id)) + 
  labs(x="frequency",y="intensity") +
  scale_x_log10() +
  theme_bw()

#w/o exoplanet
fig_12 <- train_ps_tidy %>%
  filter(id %in% no_ep_3) %>%
  ggplot(.,aes(f,ints)) + 
  geom_line() +
  facet_grid(rows=vars(id)) +
  labs(x="frequency",y="intensity") +
  theme_bw()

# scale x to log10 to emphasize low-frequencies in visualization
fig_13 <- train_ps_tidy %>%
  filter(id %in% no_ep_3) %>%
  ggplot(.,aes(f,ints)) + 
  geom_line() +
  facet_grid(rows=vars(id)) +
  labs(x="frequency",y="intensity") +
  scale_x_log10() +
  theme_bw()

# print figures

fig_10

fig_11

fig_12

fig_13

# filter periodograms

# make empty matrix to hold filters
filters <-matrix(ncol = 1599, nrow = 0)

f <- 0:1598/3196

# calculate the average periodogram for a star w/ exoplanets

avg_ep <- train_ps_standardized[1:37,] %>%
  colMeans()

auc1 <- trapz(f,avg_ep)

avg_ep_norm <- avg_ep/auc1

filters <- rbind(filters,avg_ep_norm)

# calculate the average periodogram for a star w/ exoplanets

avg_no_ep <-  train_ps_standardized[38:5087,] %>%
  colMeans()

auc2 <- trapz(f,avg_no_ep)

avg_no_ep_norm <- avg_no_ep/auc2

filters <- rbind(filters,avg_no_ep_norm)

# calculate the difference between the two average and take the square amplitudes

diff <- (avg_no_ep-avg_ep)^2

auc3 <- trapz(f,diff)

diff_norm <- diff/auc3

filters <- rbind(filters,diff_norm)

# convert to tidy for visualization

filters_tidy <- filters %>%
  as_tibble() %>%
  cbind(name=c("avg_ep","avg_no_ep","diff"),.) %>%
  gather(key,ints,"1":"1599") %>% # convert to tidy
  mutate(f=(as.numeric(key)-1)/3196) %>% # add frequency
  select(-key)

# make figures for averages and diff periodograms

fig_14 <- filters_tidy %>%
  ggplot(aes(x=f,y=ints)) +
  geom_line()+
  facet_grid(rows=vars(name)) +
  labs(x="frequency",y="intensity") +
  theme_bw()

# create filter mask
sm <- with(filter(filters_tidy,name=="diff"),ksmooth(f,ints,kernel="normal",bandwidth = 63/3196))

filter_norm <- sm$y/trapz(f,sm$y)

# make overlay of diff and filter 

fig_15 <-filter(filters_tidy,name=="diff") %>% mutate(sm_ints=filter_norm) %>%
  ggplot(aes(x=f,y=ints)) +
  geom_line()+
  geom_line(aes(f,sm_ints),color="red") +
  labs(x="frequency",y="intensity") +
  annotate(geom ="text",x=0.04,y=50,label="region I") +
  annotate(geom ="text",x=0.1,y=50,label="region II") +
  annotate(geom ="text",x=0.18,y=50,label="region III") +
  annotate(geom ="text",x=0.27,y=50,label="region IV") +
  scale_y_log10() +
  theme_bw()

#print figures

fig_14

fig_15

# filter

train_ps_filt <- sweep(train_ps_standardized, 1, filter_norm, FUN = "*")

# convert to tidy

train_ps_filt_tidy <- train_ps_filt%>%
  as_tibble() %>%
  cbind(select(train_data,LABEL),.) %>% # add labels
  cbind(id=as.numeric(rownames(.)),.) %>% # add ID
  gather(key,ints, "1":"1599") %>% # convert to tidy
  mutate(f=(as.numeric(key)-1)/3196) %>% # add frequency
  select(id,LABEL,f,ints) # throw away column names

# extract periodograms for examples in fig. 1 and fig. 3
# w/ exoplanets
fig_16 <- train_ps_filt_tidy %>%
  filter(id %in% ep_3) %>%
  ggplot(.,aes(f,ints)) + 
  geom_line() +
  facet_grid(rows=vars(id)) + 
  labs(x="frequency",y="intensity") +
  theme_bw()

#w/o exoplanet
fig_17 <- train_ps_filt_tidy %>%
  filter(id %in% no_ep_3) %>%
  ggplot(.,aes(f,ints)) + 
  geom_line() +
  facet_grid(rows=vars(id)) +
  labs(x="frequency",y="intensity") +
  theme_bw()

#print figures

fig_16

fig_17

# perform pca on untreated periodograms and on the filtered ones and compare

# perform pca (normalized and center data in function,)
pca_ps <- prcomp(train_ps_standardized[,1:959],scale=FALSE,center=FALSE)
pca_ps_filt <- prcomp(train_ps_filt[,1:959],scale=FALSE,center=FALSE)

# extract cumulative proportion of variance from summarys

untreated <- summary(pca_ps)$importance[3,]

filtered <- summary(pca_ps_filt)$importance[3,]

fig_18 <- cbind(untreated,filtered) %>% 
  as_tibble() %>%
  ggplot(aes(x=1:959,y=untreated)) +
  geom_line()+
  geom_line(aes(x=1:959,y=filtered,color="red"))+
  labs(x="# principal components",y="cumulative proportion of explained variability") +
  theme_bw() +
  theme(legend.position = "none")

# create scatter plots of first 10 principal components,
#limit data to first 1000 rows (planet w/ exoplanets are in the first 37 rows)
fig_19 <- data.frame(PC1 = pca_ps_filt$x[1:1000,1],
                     PC2 = pca_ps_filt$x[1:1000,2],
                     label=train_data[1:1000,1]) %>%
  ggplot(aes(PC1, PC2, fill=as.factor(LABEL)))+
  geom_point(cex=3, pch=21)

fig_20 <- data.frame(PC1 = pca_ps$x[1:1000,1],
                     PC2 = pca_ps$x[1:100,2],
                     label=train_data[1:1000,1]) %>%
  ggplot(aes(PC1, PC2, fill=as.factor(LABEL)))+
  geom_point(cex=3, pch=21)

fig_21 <- data.frame(PC3 = pca_ps_filt$x[1:1000,3],
                     PC4 = pca_ps_filt$x[1:1000,4],
                     label=train_data[1:1000,1]) %>%
  ggplot(aes(PC3, PC4, fill=as.factor(LABEL)))+
  geom_point(cex=3, pch=21)

fig_22 <- data.frame(PC3 = pca_ps$x[1:1000,3],
                     PC4 = pca_ps$x[1:1000,4],
                     label=train_data[1:1000,1]) %>%
  ggplot(aes(PC3, PC4, fill=as.factor(LABEL)))+
  geom_point(cex=3, pch=21)

fig_23 <- data.frame(PC5 = pca_ps_filt$x[1:1000,5],
                     PC6 = pca_ps_filt$x[1:1000,6],
                     label=train_data[1:1000,1]) %>%
  ggplot(aes(PC5, PC6, fill=as.factor(LABEL)))+
  geom_point(cex=3, pch=21)

fig_24 <- data.frame(PC5 = pca_ps$x[1:1000,5],
                     PC6 = pca_ps$x[1:1000,6],
                     label=train_data[1:1000,1]) %>%
  ggplot(aes(PC5, PC6, fill=as.factor(LABEL)))+
  geom_point(cex=3, pch=21)

fig_25 <- data.frame(PC7 = pca_ps_filt$x[1:1000,7],
                     PC8 = pca_ps_filt$x[1:1000,8],
                     label=train_data[1:1000,1]) %>%
  ggplot(aes(PC7, PC8, fill=as.factor(LABEL)))+
  geom_point(cex=3, pch=21)

fig_26 <- data.frame(PC7 = pca_ps$x[1:1000,7],
                     PC8 = pca_ps$x[1:1000,8],
                     label=train_data[1:1000,1]) %>%
  ggplot(aes(PC7, PC8, fill=as.factor(LABEL)))+
  geom_point(cex=3, pch=21)

fig_27 <- data.frame(PC9 = pca_ps_filt$x[1:1000,9],
                     PC10 = pca_ps_filt$x[1:1000,10],
                     label=train_data[1:1000,1]) %>%
  ggplot(aes(PC9, PC10, fill=as.factor(LABEL)))+
  geom_point(cex=3, pch=21)

fig_28 <- data.frame(PC9 = pca_ps$x[1:1000,9],
                     PC10 = pca_ps$x[1:1000,10],
                     label=train_data[1:1000,1]) %>%
  ggplot(aes(PC9, PC10, fill=as.factor(LABEL)))+
  geom_point(cex=3, pch=21)

# print figures

fig_18

fig_19

fig_20

fig_21

fig_22

fig_23

fig_24

fig_25

fig_26

fig_27

fig_28

# create models
# downsampled methods
## define the settings for cross-validation.
ctrl_pars <- trainControl(method="repeatedcv",
                          classProbs = TRUE , #necessary for ROC curves
                          number=10, repeats = 5, # 10 fold cross-validation, 5 repeats
                          summaryFunction = twoClassSummary, 
                          sampling = "down") # downsampling
## define training inputs
TrainX <- pca_ps_filt$x
## define labels
# change labels to make original class "2" of stars w/ exoplanets the positive class "X1"
temp <- mutate(train_data,LABEL = ifelse(LABEL=="2","1","2")) 
# change labels to make original class "2" of stars w/ exoplanets the positive class "X1"
TrainY <- make.names(temp$LABEL)
# set up random forest model
set.seed(1, sample.kind="Rounding")
rf_model_down <- train(x=TrainX, y=TrainY,
                       method = "rf",
                       trControl = ctrl_pars,
                       metric = "ROC",
                       tuneGrid = data.frame(mtry =seq(3,99,3)))

## set up knn model
set.seed(1, sample.kind="Rounding")
knn_model_down <- train(x=TrainX, y=TrainY,
                        method = "knn",
                        trControl = ctrl_pars, 
                        metric = "ROC",
                        tuneGrid = data.frame(k =seq(1,100,2)))

## set up svm model
set.seed(1, sample.kind="Rounding")
svm_model_down <- train(x=TrainX, y=TrainY,
                        method = "svmRadial",
                        trControl = ctrl_pars,
                        metric = "ROC",
                        tuneLength = 10)

# create anomaly detector

# define the settings for cross-validation.
## define training inputs (stars w/ exoplanets only)
TrainX <- pca_ps_filt$x[1:37,]
## define labels
# change labels to make original class "2" of stars w/ exoplanets the positive class "X1"
temp <- mutate(train_data,LABEL = ifelse(LABEL=="2",TRUE,FALSE)) 
TrainY <- temp$LABEL[1:37]

# set up anomaly model in tune, no one-class predictor in caret available
set.seed(1, sample.kind="Rounding")
svm_model_anomaly <- tune(svm, train.x = TrainX, train.y = as.factor(TrainY),
                          kernel="radial",
                          type="one-classification",
                          ranges = list(gamma=c(0.1,0.5,1,2,4), 
                                        cost = c(0.1,1,10,100,1000)
                          ),
                          tunecontrol = tune.control(nrepeat = 10, #10 repeates
                                                     sampling = "cross",
                                                     cross = 10, #10-fold cross-validation
                                                     performances = TRUE)
)

#make predictions on the training data for all models
pred_rf <- predict(rf_model_down, pca_ps_filt$x)

pred_knn <- predict(knn_model_down, pca_ps_filt$x)

pred_svm <- predict(svm_model_down, pca_ps_filt$x)

pred_anomaly<- predict(svm_model_anomaly$best.model, pca_ps_filt$x)

#model selection
# prepare confusion matrices
# change labels to make original class "2" of stars w/ exoplanets the positive class "X1"
ref <- mutate(train_data,LABEL = ifelse(LABEL=="2","X1","X2"))
# change labels to make original class "2" of stars w/ exoplanets the positive class TRUE
ref_anom <- mutate(train_data,LABEL = ifelse(LABEL=="2",TRUE,FALSE))  

cm_rf <- table(Predicted_rf =pred_rf, Reference =ref$LABEL)

cm_knn <- table(Predicted_knn =pred_knn, Reference =ref$LABEL)

cm_svm <- table(Predicted_svm =pred_svm, Reference =ref$LABEL)

cm_anom <- table(Predicted_anom =pred_anomaly, Reference =ref_anom$LABEL)

# print matrices

cm_rf

cm_knn

cm_svm

cm_anom

# model validation

# process validation data

# calculation of power spectra
validation_ps <- matrix(ncol = 1599, nrow = 0) #create empty ps matrix

# create periograms
for (i in 1:570) { # extract time series row by row
  temp<-as.matrix(validation_data[i,2:3197])
  ints2 <- abs(fft(temp))^2/3196 # calculate intensities from squared amplitudes
  scaled_ints <- (4/3196)*ints2[1:1599] # re-scale
  validation_ps <- rbind(validation_ps,scaled_ints) # create ps matrix
}

# filter, min-max scaling
validation_ps_centered <- sweep(validation_ps, 1, apply(validation_ps,1,min))
validation_ps_standardized <- sweep(validation_ps_centered, 1,
                                    apply(validation_ps,1,max), FUN = "/")

validation_ps_filt <- sweep(validation_ps_standardized, 1,
                            filter_norm, FUN = "*")

# transform the validation dataset to obtain eigenvalues
pca_val <- validation_ps_filt[,1:959] %*% pca_ps_filt$rotation

#model 1: random forest predictor only

pred_val_1 <- predict(rf_model_down, pca_val)

#model 2: random forest predictor + anomaly detector

data_tier1 <- cbind(pca_val,
                    pred_val_1,
                    ind=1:length(validation_data$LABEL)) %>% 
  .[.[,"pred_val_1"] =="1",] # select "X1" from first prediction and add index ind

data_tier2 <- cbind(data_tier1,
                    pred_anom=predict(svm_model_anomaly$best.model,
                                      data_tier1[,1:959])) # combine with anomaly prediction

pred_val_2 <- ifelse(1:length(validation_data$LABEL) %in% data_tier2[data_tier2[,"ind"]],"TRUE","FALSE") # final prediction

#model 3: anomaly detector only

pred_val_3 <- predict(svm_model_anomaly$best.model, pca_val)

# evaluation of predictive performance
# change labels to make original class "2" of stars w/ exoplanets the positive class "X1"
ref_val <- mutate(validation_data,LABEL = ifelse(LABEL=="2","X1","X2")) 
# change labels to make original class "2" of stars w/ exoplanets the positive class TRUE
ref_val_anom <- mutate(validation_data,LABEL = ifelse(LABEL=="2",TRUE,FALSE)) 

cm_model_1 <- table(Predicted =pred_val_1, Reference =ref_val$LABEL)

cm_model_2 <- table(Predicted =pred_val_2, Reference =ref_val_anom$LABEL)

cm_model_3 <- table(Predicted =pred_val_3, Reference =ref_val_anom$LABEL)


#print rf matrix
print("Confusion matrix random forest model only")
cm_model_1

#print 2-tier matrix
print("Confusion matrix two tier model")
cm_model_2

#print anomaly detector matrix
print("Confusion matrix anomaly detector")
cm_model_3
