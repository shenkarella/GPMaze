#All Maze1 Data for the GP Manuscript:
#Began on July 8, 2025
#All the data analysis (written by Ella, and taken from Peiyao + Dan's code)
#To organize things a little bit better & make sure we're using accurate models + contrasts
#Heavily commented so we can refer back to ONE script. 


####SETUP####
library(dplyr)
library(ez)
library(lme4)
library(lmerTest)
library(optimx)
library(nloptr)
library(dfoptim)
library(parallel)
library(NbClust)
library(cluster)
library(class)
library(plyr)
#to not crash:
install.packages("RhpcBLASctl")
library(RhpcBLASctl)
library(lme4)

#Set wd
setwd("/Users/Ellashenkar_1/Desktop/Maze Organized")


#load data, pre-processed
mydata <- read.delim("/Users/Ellashenkar_1/Desktop/Maze Organized/GPMaze1dataForR_210202.tsv", header = TRUE, sep = "\t")

## bad subs remove based on the badsubs criteria
# Did not judge inanimates as 	significantly worse than plurals 
# 5cbcb382efda5c001a62a413 4.2 4.4 4.6 3.6 3.8 4 4.6 4 4.26667 0 1 2 1 0 23 Male    -0.06667 17
# 5df7d398591dbc000fad99b8 5.6 5.8 6.4 5.2 6.8 6.2 6.4 5.6 5.53333 0 7 1 2 2 37 Non-binary    0.06667 17
# 5d1f4bf21218f100017dad45 5.8 4.8 5.6 5.6 5.4 5.6 5.6 5.2 5.53333 0 8 0 3 1 40 Non-binary    0.26667 17
# 
# # Not a good subject number
# https://spellout.net/ibexexps/ellash729/BjorkmanCopy/experiment.html
# 
# # Don't have bjorkman data (should add this person back in when possible)
# #5e7b4ec9db5aba01700b311b	0	10	1	0	3	33	Female
# 
# # don't have bjorkman nor individual differences, also very low accuracy
# 5ee9f7b5979e5302a8d29301

subs_out <- c( "5cbcb382efda5c001a62a413", "5df7d398591dbc000fad99b8", "5d1f4bf21218f100017dad45", "5e7b4ec9db5aba01700b311b",  
               "5ee9f7b5979e5302a8d29301")
`%notin%` <- Negate(`%in%`)
mydata <- mydata[mydata$sID%notin%subs_out,]

cierr <- function(x) qnorm(0.975)*sd(x)/sqrt(length(x))
# data check
# summary(mydata)
# str(mydata)
# unique(mydata$sID)
# xtabs(~sID+match+plsg,mydata)
# xtabs(~corr, mydata)
z. <- function (x) scale (x)

## more play around
# filler <- mydata[mydata$cond == 'fil',]
# they <- mydata[mydata$word == 'they',]
# wrong <- mydata[mydata$corr == 'no',]

# 6) Describe exactly how outliers will be defined and handled, and your precise rule(s) for excluding observations.
# 1. If a participant makes an error prior to the critical region (pronoun and subsequent words), that trial will be omitted
# 2. If a participant?s accuracy is more than 2.5 SD lower than the median, their data will be excluded from the analyses.
# 3. If a participant?s average RT is more than 3 SD higher than the median, their data will be excluded from analyses.
# 4. If an individual word reaction time < 120 ms it will be excluded from the analyses.
# 5. If the reaction time on an individual word > 3 SD away from the condition mean for that word, it will be excluded..
# 6. Data from participants with <50% of their data remaining after the exclusions above will be discarded.

## create trial order
## exclude fillers
prodata <- mydata[mydata$cond != 'fil',]
critical <- unique(prodata[c(1,2,3)])
# total participants
unique(mydata$sID)
xtabs(~sID, critical) ## each participant has 120 critical trials and a total of 76 participants
critical$trialorder <- rep(1:120,76)


# data preprocessing
# do not exclude fillers
#nofiller <- mydata[mydata$cond != 'fil',]
xtabs(~corr, prodata) # check what's under corr
# check total trial in each participant, without filler 120
tt <- unique(prodata[c(1,2)])
ttcount <- data.frame(xtabs(~sID, tt))

## 2. If a participant's accuracy is more than 2.5 SD lower than the median, their data will be excluded from the analyses.
### based on erros in any position during a sentence, fillers included
erros <- prodata[prodata$corr == 'no', ]
error_trial <- unique(erros[c(1,2)])
error_count <- data.frame(xtabs(~sID,error_trial))
error_count$acc <- 1-(error_count$Freq/120)
error_count$lowcut <-  median(error_count$acc)-2.5*(sd(error_count$acc))
outlier_lowacc <- error_count[error_count$acc < error_count$lowcut,]


### endup with 2 outliers: 5e6b034bffc1a90c3fe8c234, 5ea84aa3583191095b292904,

subs_out <- c( "5e6b034bffc1a90c3fe8c234", "5ea84aa3583191095b292904")
`%notin%` <- Negate(`%in%`)
prepro2 <- prodata[prodata$sID%notin%subs_out,]

# 3. If a participant?s average RT is more than 3 SD higher than the median, their data will be excluded from analyses.
prepro2$resposetime <- as.numeric(as.character(prepro2$RT))
averRT <- ddply(prepro2, .(sID), summarise, rt=mean(resposetime, na.rm = TRUE)) 
averRT$highcut <-  median(averRT$rt)+3*(sd(averRT$rt))
outlier_highrt <- averRT[averRT$rt > averRT$highcut,]

## end up excluding one participant: 5eaa0e5b3b32cf1b13880d65
subs_out <- c("5eaa0e5b3b32cf1b13880d65")
`%notin%` <- Negate(`%in%`)

prepro3 <- prepro2[prepro2$sID%notin%subs_out,]

# 4. If an individual word reaction time < 120 ms it will be excluded from the analyses.
prepro4 <- subset(prepro3, resposetime >= 120 | is.na(resposetime))

#4a. no fillers
prepro4_nofiller <- subset(prepro4, cond != 'fil' | is.na(resposetime))


# 1. If a participant makes an error prior to the critical region (pronoun and subsequent words), that trial will be omitted
# critical region is associcated with 0 and 1 in word position (wdPos)
# figure out the sID and item number conbination
precri_error <- prepro4_nofiller[prepro4_nofiller$wdPos < 0 & prepro4_nofiller$corr == "no",]
totalerror <- unique(precri_error[c(1,2,3)])
## exclude error trials
prepro5 <- anti_join(prepro4_nofiller, totalerror, by=c("sID", "item", "cond"))


## 5. If the reaction time on an individual word > 3 SD away from the condition mean for that word, it will be excluded..
indiword <- ddply(prepro5, .(match, plsg, wdPos), summarise, meanrt=mean(resposetime, na.rm = TRUE), meansd=sd(resposetime, na.rm = TRUE)) 
indiword$highcut <- indiword$meanrt+3*indiword$meansd
prepro5.5 <- merge(prepro5,indiword, by = c('match','plsg','wdPos'), all.x=T, sort=FALSE)
prepro6 <- subset(prepro5.5, resposetime <= highcut |is.na(resposetime))



## 6. Data from participants with <50% of their data remaining after the exclusions above will be discarded.
# only based on critical trials
totalremain <- unique(prepro6[c(4,5)])
final_count <- data.frame(xtabs(~sID,totalremain))
final_count$accurate <- final_count$Freq/120
outlier_few <- final_count[final_count$accurate < 0.5,]
## removing no one in this step



# back to prepro5 which contains errors
final <- merge(prepro6, critical, by = c("sID", "item", "cond"), sort = FALSE)
unique(final$sID)
# final n =73


final$item <- as.factor(final$item)
final$sID <- as.factor(final$sID)
final$match <- as.factor(final$match)
final$plsg <- as.factor(final$plsg)
final$tPhobia <- as.numeric(as.character(final$tPhobia))
final$gEssential <- as.numeric(as.character(final$gEssential))
final$useOfThey <- as.numeric(as.character(final$useOfThey))
final$gIdentity <- as.numeric(as.character(final$gIdentity))
final$nBAcceptance <- as.numeric(as.character(final$nBAcceptance))
final$age <- as.numeric(as.character(final$age))


final$accuracy <- ifelse(final$corr == "yes", 1,0)

## descriptives
pronoun <- final[final$wdPos == 0,]
afterpro <- final[final$wdPos == 1,]

pronoun_acc <- ddply(pronoun, .(match, plsg), summarise, acc=mean(accuracy, na.rm = TRUE))

pronoun_rt <- ddply(pronoun[pronoun$accuracy == 1,], .(match, plsg), summarise, rt = mean(resposetime, na.rm = TRUE))
afterpro_acc <- ddply(afterpro, .(match, plsg), summarise, acc=mean(accuracy, na.rm = TRUE))
afterpro_rt <- ddply(afterpro[afterpro$accuracy == 1,], .(match, plsg), summarise, rt = mean(resposetime, na.rm = TRUE))


#Setting contrasts for pronoun
pronoun$matchness <-as.factor(pronoun$match)
pronoun$plurality <- as.factor(pronoun$plsg)
contrasts(pronoun$matchness) <- c(-0.5, 0.5) ## match is the first level, and the mismatch is the second level
contrasts(pronoun$plurality)  <- c(0.5, -0.5) ## pl is the first, sg is the second level

#Setting contrasts for after the pronoun
afterpro$matchness <-as.factor(afterpro$match)
afterpro$plurality <- as.factor(afterpro$plsg)
contrasts(afterpro$matchness) <- c(-0.5, 0.5) ## match is the first level, and the mismatch is the second level
contrasts(afterpro$plurality)  <- c(0.5, -0.5) ## pl is the first, sg is the second level

## descriptives
pronoun <- final[final$wdPos == 0,]
afterpro <- final[final$wdPos == 1,]


#--------------------------------------------------Pronoun Accuracy--------------------------------------------------
pronoun_acc <- ddply(pronoun, .(match, plsg), summarise, acc=mean(accuracy, na.rm = TRUE))
#  match plsg       acc
#1   mat   pl 0.9855228
#2   mat   sg 0.9933606
#3   mis   pl 0.9860176
#4   mis   sg 0.9899842
#----------GLMER models for accuracy----------

##Accuracy model 
#Maximal model, using nlminbwrap as optimizer
pro_acc6<- glmer(accuracy ~ matchness*plurality + (1+matchness*plurality|sID) + (1+matchness*plurality|item), data = pronoun, family = binomial, control = glmerControl(optimizer ='nlminbwrap'))
summary(pro_acc6)
#Fixed effects:
#Estimate Std. Error z value Pr(>|z|)    
#(Intercept)             5.4282     0.2921  18.585   <2e-16 ***
#  matchness1             -0.1034     0.6325  -0.164    0.870    
#plurality1             -0.7787     0.5648  -1.379    0.168    
#matchness1:plurality1   0.8389     1.1194   0.749    0.454  
df_pro_acc6  <- data.frame(
  effect = rownames(summary(pro_acc6)$coefficients),
  estimate = summary(pro_acc6)$coefficients[, "Estimate"],
  se = summary(pro_acc6)$coefficients[, "Std. Error"],
  df = summary(pro_acc6)$coefficients[, "df"],
  t_value = summary(pro_acc6)$coefficients[, "t value"],
  pvalue = summary(pro_acc6)$coefficients[, "Pr(>|t|)"]
)

##Using optimizer based on: 
diff_optims <- allFit(pro_acc4, maxfun = 1e5, parallel = 'multicore', ncpus = detectCores())
is.OK <- sapply(diff_optims, is, "merMod")
diff_optims.OK <- diff_optims[is.OK]
lapply(diff_optims.OK,function(x) x@optinfo$conv$lme4$messages)
#$bobyqa
#[1] "boundary (singular) fit: see ?isSingular"

#$Nelder_Mead
#[1] "unable to evaluate scaled gradient"                                       
#[2] "Model failed to converge: degenerate  Hessian with 1 negative eigenvalues"

#$nlminbwrap
#[1] "boundary (singular) fit: see ?isSingular"

#$nmkbw
#[1] "boundary (singular) fit: see ?isSingular"

#$`optimx.L-BFGS-B`
#[1] "boundary (singular) fit: see ?isSingular"

#$nloptwrap.NLOPT_LN_NELDERMEAD
#[1] "boundary (singular) fit: see ?isSingular"

#$nloptwrap.NLOPT_LN_BOBYQA
#[1] "boundary (singular) fit: see ?isSingular"

#instead, try gusing nlminbwrap
pro_acc5<- glmer(accuracy ~ matchness*plurality + (1+matchness+plurality|sID) + (1+matchness+plurality|item), data = pronoun, family = binomial, control = glmerControl(optimizer ='nlminbwrap'))
summary(pro_acc5)
#Fixed effects:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)             5.0972     0.3090  16.495   <2e-16 ***
#  matchness1              0.3677     0.4075   0.902   0.3669    
#plurality1             -1.0244     0.5056  -2.026   0.0427 *  
#  matchness1:plurality1  -0.4038     1.1951  -0.338   0.7355    


#--------------------------------------------------Pronoun Reaction Time--------------------------------------------------
pronoun_rt <- ddply(pronoun[pronoun$accuracy == 1,], .(match, plsg), summarise, rt = mean(resposetime, na.rm = TRUE))
#  match plsg       rt
#1   mat   pl 644.8134
#2   mat   sg 648.7342
#3   mis   pl 750.6239
#4   mis   sg 810.8962

#Maximal model of reaction time with optimizer (result is the same as non-maximal models)
pro_rt<- lmer(resposetime ~ matchness*plurality + (1+matchness*plurality|sID) + (1+matchness*plurality|item), data = pronoun[pronoun$accuracy ==1,], control = lmerControl(optimizer = "nloptwrap",optCtrl = list(algorithm = "NLOPT_LN_BOBYQA", maxit = 1e6)))
summary(pro_rt)
#Fixed effects:
#Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)            711.927     11.900   79.046  59.825  < 2e-16 ***
#  matchness1             135.114      7.455   72.946  18.125  < 2e-16 ***
#  plurality1             -31.022      6.153   66.150  -5.042 3.82e-06 ***
#  matchness1:plurality1  -57.061      9.121 7234.640  -6.256 4.17e-10 ***


#Checking best fit in terms of reaction time 
diff_optims <- allFit(pro_rt_sumcoding, maxfun = 1e5, parallel = 'multicore', ncpus = detectCores())
is.OK <- sapply(diff_optims, is, "merMod")
diff_optims.OK <- diff_optims[is.OK]
lapply(diff_optims.OK,function(x) x@optinfo$conv$lme4$messages)
#suggests that for pro_rt_sumcoding, $nloptwrap.NLOPT_LN_BOBYQA is the best fit. 

#follow-up analysis:
#singular
pro_rt_sg<- lmer(resposetime ~ matchness + (1+matchness|sID) + (1+matchness|item), data = pronoun[pronoun$accuracy ==1 & pronoun$plsg== 'sg',], control = lmerControl(optimizer = "nloptwrap",optCtrl = list(algorithm = "NLOPT_LN_BOBYQA", maxit = 1e6)))
summary(pro_rt_sg)
#Estimate Std. Error     df t value Pr(>|t|)    
#(Intercept)   727.88      12.97  81.59   56.13   <2e-16 ***
#  matchness1    162.47      11.10  88.86   14.64   <2e-16 ***

df_pro_rt_sg  <- data.frame(
  effect = rownames(summary(pro_rt_sg)$coefficients),
  estimate = summary(pro_rt_sg)$coefficients[, "Estimate"],
  se = summary(pro_rt_sg)$coefficients[, "Std. Error"],
  df = summary(pro_rt_sg)$coefficients[, "df"],
  t_value = summary(pro_rt_sg)$coefficients[, "t value"],
  pvalue = summary(pro_rt_sg)$coefficients[, "Pr(>|t|)"]
)
View(df_pro_rt_sg)


#plural 
pro_rt_pl<- lmer(resposetime ~ matchness + (1+matchness|sID) + (1+matchness|item), data = pronoun[pronoun$accuracy ==1 & pronoun$plsg== 'pl',], control = lmerControl(optimizer = "nloptwrap",optCtrl = list(algorithm = "NLOPT_LN_BOBYQA", maxit = 1e6)))
summary(pro_rt_pl)
#Fixed effects:
#  Estimate Std. Error      df t value Pr(>|t|)    
#(Intercept)  696.946     11.695  78.184   59.60   <2e-16 ***
#  matchness1   107.122      8.328  78.202   12.86   <2e-16 ***
df_pro_rt_pl  <- data.frame(
  effect = rownames(summary(pro_rt_pl)$coefficients),
  estimate = summary(pro_rt_pl)$coefficients[, "Estimate"],
  se = summary(pro_rt_pl)$coefficients[, "Std. Error"],
  df = summary(pro_rt_pl)$coefficients[, "df"],
  t_value = summary(pro_rt_pl)$coefficients[, "t value"],
  pvalue = summary(pro_rt_pl)$coefficients[, "Pr(>|t|)"]
)
View(df_pro_rt_pl)


#--------------------------------------------------After Pronoun Accuracy--------------------------------------------------

afterpro_acc <- ddply(afterpro, .(match, plsg), summarise, acc=mean(accuracy, na.rm = TRUE))
#  match plsg       acc
#1   mat   pl 0.9750927
#2   mat   sg 0.9759345
#3   mis   pl 0.9779487
#4   mis   sg 0.9779180

afterpro_acc<- glmer(accuracy ~ matchness*plurality + (1+matchness*plurality|sID) + (1+matchness*plurality|item), data = afterpro, family = binomial, control = glmerControl(optimizer ='nlminbwrap'))
summary(afterpro_acc)
#Fixed effects:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)             4.1000     0.1444  28.395   <2e-16 ***
#  matchness1              0.4364     0.2584   1.689   0.0913 .  
#plurality1             -0.1835     0.2548  -0.720   0.4714    
#matchness1:plurality1   0.1333     0.5086   0.262   0.7933  


#-------------------------------------------------- After Pronoun Reaction Time--------------------------------------------------

afterpro_rt <- ddply(afterpro[afterpro$accuracy == 1,], .(match, plsg), summarise, rt = mean(resposetime, na.rm = TRUE))
#match plsg       rt
#1   mat   pl 778.2967
#2   mat   sg 773.2539
#3   mis   pl 797.8668
#4   mis   sg 807.0172


#reaction time of afterpronoun
#models that failed to converge
afterpro_rt<- lmer(resposetime ~ match*plsg + (1+match*plsg|sID) + (1+match*plsg|item), data = afterpro[afterpro$accuracy ==1,], control = lmerControl(optimizer = "nloptwrap"))
afterpro_rt2<- lmer(resposetime ~ match*plsg + (1+match*plsg|sID) + (1+match*plsg|item), data = afterpro[afterpro$accuracy ==1,],control = lmerControl(optimizer ='Nelder_Mead'))

#checking better optimizers: 
diff_optims <- allFit(afterpro_rt, maxfun = 1e5, parallel = 'multicore', ncpus = detectCores())
is.OK <- sapply(diff_optims, is, "merMod")
diff_optims.OK <- diff_optims[is.OK]
lapply(diff_optims.OK,function(x) x@optinfo$conv$lme4$messages)

#maximal model:
afterpro_rt2<- lmer(resposetime ~ match*plurality + (1+match|sID) + (1+match|item), data = afterpro[afterpro$accuracy ==1,],control = lmerControl(optimizer ='Nelder_Mead'))
summary(afterpro_rt2)
#                Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)      777.709     15.406  116.340  50.482  < 2e-16 ***
#  matchmis          21.105      7.172 3183.485   2.943  0.00328 ** 
#  plsgsg            -7.198      7.128 7248.805  -1.010  0.31259    
#matchmis:plsgsg   14.515     10.028 7357.060   1.447  0.14783    

#follow-up analysis:

#SINGULAR
afterpro_rt3<- lmer(resposetime ~ matchness + (1+matchness|sID) + (1+matchness|item), data = afterpro[afterpro$accuracy ==1 & afterpro$plurality == "sg",],control = lmerControl(optimizer ='Nelder_Mead'))
summary(afterpro_rt3)
#Fixed effects:
#  Estimate Std. Error     df t value Pr(>|t|)    
#(Intercept)   788.03      14.09  95.54   55.92  < 2e-16 ***
#  matchness1     36.80       7.23 124.06    5.09 1.29e-06 ***
#PLURAL
afterpro_rt_pl<- lmer(resposetime ~ matchness + (1+matchness|sID) + (1+matchness|item), data = afterpro[afterpro$accuracy ==1 & afterpro$plurality == "pl",],control = lmerControl(optimizer ='nloptwrap'))
summary(afterpro_rt_pl)
# Estimate Std. Error      df t value Pr(>|t|)    
#(Intercept)  788.348     15.876 101.202   49.66   <2e-16 ***
#  matchness1    20.886      8.094  66.190    2.58   0.0121 *  





##_-------------------------------clustering results ----------------------------
d<-read.table("rating_Sadiethesis.txt", header=TRUE)
# use Hclust to label the old case, and split into training and testing data, and improve the accuracy before running on the new data
row.names(d) <- d$ID
d$ID <- NULL
d.scale <- scale(d)
nb <- NbClust(d.scale, distance = "euclidean", min.nc = 2, max.nc = 10, method = "complete", index ="all")
## we know that k=3, 3 clusters
distance <- get_dist(d.scale)
dm <-  dist(d.scale, method = "euclidean") # dissimilarity matrix
hc1 <- hclust(dm, method = "complete")
plot(hc1, cex = 0.6, hang = -1) 
sub_grp <- cutree(hc1, k=3)
table(sub_grp)
rect.hclust(hc1, k = 3, border = 2:4)
d$cluster <- sub_grp

## label the cluster

d$cluster <- as.factor(d$cluster)
d$user <- factor(d$cluster, levels = c("1", "2", "3"), labels = c("inno", "non", "super"))
table(d$user)

## normalizing numeric data

# function normalize 
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x))) }
d_normalized <- as.data.frame(lapply(d[1:9], normalize))
summary(d_normalized$plural)

## splitting into training vs. testing data by randomly selecting
set.seed(123)
dat.d <- sample(1:nrow(d_normalized),size=nrow(d_normalized)*0.7,replace = FALSE) #random selection of 70% data.

train.d <- d_normalized[dat.d,] # 70% training data
test.d <- d_normalized[-dat.d,] # remaining 30% test data

### create lables for each row in train and test sets
train.d_labels <- d[dat.d,11]
test.d_labels <- d[-dat.d,11]


## now train the model on the data, k here equals square root of total number
d_test_pred.1 <- knn(train = train.d, test = test.d,cl = train.d_labels, k=12)

## model evaluation
ACC.1 <- 100 * sum(test.d_labels == d_test_pred.1)/NROW(test.d_labels)
CrossTable(x=test.d_labels, y=d_test_pred.1, prop.chisq = FALSE)

## try k=1-13 the best is 11 & 9, 8 and 7, all 91.11, decrease at 21
## will go with 11, since it's closer to the 12
d_test_pred.2 <- knn(train = train.d, test = test.d,cl = train.d_labels, k=11)

## model evaluation
ACC.2 <- 100 * sum(test.d_labels == d_test_pred.2)/NROW(test.d_labels)
CrossTable(x=test.d_labels, y=d_test_pred.2, prop.chisq = FALSE)


## load new maze1 data 
cdata <- subset(mydata, select=c('sID','pl','qu','ngn','spk','gn','gf','nbname','gname','inanim'))
d_new <- unique(cdata)
row.names(d_new) <- d_new$sID
d_new$sID <- NULL


d_new_n <- as.data.frame(lapply(d_new[1:9], normalize))
# train on new data, with k=11

d_test_new <- knn(train = train.d, test = d_new_n,cl = train.d_labels, k=11)
table(d_test_new)
d_new$user <- d_test_new 
d_new$cluster[d_new$user =='non'] <- 1
d_new$cluster[d_new$user =='inno'] <- 2
d_new$cluster[d_new$user =='super'] <- 3
d_new <- d_new[, c(1:9,11)]

#add back ID
dclust2 <- group_by(d_new, cluster)
d2 <- unique(cdata)
dclust2$sID <- d2$sID

##plot the new 
dmeans<-aggregate(d_new,by=list(dclust2$cluster),mean) #get means for each var. for each cluster
head(dmeans) #see table for condition means by cluster
dmeans<-dmeans[,2:length(dmeans)] #remove Group.1 col
dmeans.long<-melt(dmeans,id.vars="cluster")

ggplot(dmeans.long,aes(x=variable,y=value,fill=factor(cluster)))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_discrete(name="cluster",
                      breaks=c(1, 2, 3),
                      labels=c("non-innovator", "innovater", "super-innovator"))+xlab("Condition")+ylab("Mean Naturalness Rating")

### merge the cluster information with the data
cluster_info <- dclust2[, c(10, 11)]
pronoun_cluster <- merge(pronoun, cluster_info, by = 'sID', sort = FALSE)
afterpro_cluster <- merge(afterpro, cluster_info, by = 'sID', sort = FALSE)
pronoun_cluster$cluster <- as.factor(pronoun_cluster$cluster)
afterpro_cluster$cluster <- as.factor(afterpro_cluster$cluster)

pronoun_cluster$user[pronoun_cluster$cluster == '1'] <- 'non-inno'
pronoun_cluster$user[pronoun_cluster$cluster == '2'] <- 'inno'
pronoun_cluster$user[pronoun_cluster$cluster == '3'] <- 'super'

afterpro_cluster$user[afterpro_cluster$cluster == '1'] <- 'non-inno'
afterpro_cluster$user[afterpro_cluster$cluster == '2'] <- 'inno'
afterpro_cluster$user[afterpro_cluster$cluster == '3'] <- 'super'


pronoun_cluster$clust.hC <- as.factor(pronoun_cluster$cluster)
contrasts(pronoun_cluster$clust.hC)<-  cbind(c(-1/2,1/2,0), c(1/3,1/3,-2/3))
colnames(contrasts(pronoun_cluster$clust.hC)) <- c('c1Vc2', 'c3Vothers')
contrasts(pronoun_cluster$clust.hC)

#setting other contrasts: 
pronoun_cluster$matchness <-as.factor(pronoun$match)
pronoun_cluster$plurality <- as.factor(pronoun$plsg)
contrasts(pronoun_cluster$matchness) <- c(-0.5, 0.5) ## match is the first level, and the mismatch is the second level
contrasts(pronoun_cluster$plurality)  <- c(0.5, -0.5) ## sg is the first, pl is the second

##Bad models b/c they don't use an optimizer
pro_rt_plu_cluster <- lmer(resposetime ~ match*clust.hC + (1+match|sID) + (1|item), data = pronoun_cluster[pronoun_cluster$accuracy == 1 & pronoun_cluster$plsg== 'pl',])
summary(pro_rt_plu_cluster)
#Fixed effects:
#Estimate Std. Error      df t value Pr(>|t|)    
#(Intercept)                 646.481     10.266  78.302  62.975   <2e-16 ***
#  matchmis                    107.470      7.859  71.688  13.676   <2e-16 ***
#  clust.hCc1Vc2               -40.225     24.405  70.533  -1.648   0.1038    
#clust.hCc3Vothers            39.953     21.086  69.729   1.895   0.0623 .  
#matchmis:clust.hCc1Vc2      -16.179     19.181  71.007  -0.843   0.4018    
#matchmis:clust.hCc3Vothers   39.563     16.523  69.388   2.394   0.0194 *  

#Same method as above 
pro_rt_sg_cluster <- lmer(resposetime ~ match*clust.hC + (1+match|sID) + (1+match|item), data = pronoun_cluster[pronoun_cluster$accuracy == 1 & pronoun_cluster$plsg== 'sg',])
summary(pro_rt_sg_cluster)
#Fixed effects:
#                           Estimate Std. Error     df t value Pr(>|t|)    
#(Intercept)                  651.16      10.31  70.61  63.178   <2e-16 ***
#matchmis                     160.49      11.81  89.86  13.586   <2e-16 ***
#clust.hCc1Vc2                -53.48      25.27  70.94  -2.116   0.0379 *  
#clust.hCc3Vothers             29.63      21.91  71.06   1.352   0.1806    
#matchmis:clust.hCc1Vc2        15.35      28.11  88.40   0.546   0.5864    
#matchmis:clust.hCc3Vothers    21.42      22.87  75.75   0.936   0.3521  

pro_rt_cluster <- lmer(resposetime ~ match*clust.hC*plurality + (1+match|sID) + (1+match|item), data = pronoun_cluster[pronoun_cluster$accuracy == 1 & pronoun_cluster$plsg== 'sg',])


#ALTERNATIVE METHOD:
pro_rt_sin_cluster <- lmer(resposetime ~ match*cluster + (1+match|sID) + (1|item), data = pronoun_cluster[pronoun_cluster$accuracy == 1 & pronoun_cluster$plsg== 'sg',])
summary(pro_rt_sin_cluster)
#Fixed effects:
#Estimate   Std. Error df    t value  Pr(>|t|)    
#(Intercept)            689.18      20.45  73.16  33.706    < 2e-16 ***
#  matchmis            158.47      20.11  70.57   7.879     2.91e-11 ***
#  cluster2            -55.55      25.31  70.48  -2.195     0.0315 *  
#  cluster3            -59.01      27.04  70.60  -2.182     0.0324 *  
#  matchmis:cluster2    17.07      25.12  70.39   0.680     0.4990    
#matchmis:cluster3     -10.98      26.80  70.03  -0.410     0.6832 

pro_rt_sg_cluster <- lmer(resposetime ~ match*clust.hC + (1+match|sID) + (1+match|item), data = pronoun_cluster[pronoun_cluster$accuracy == 1 & pronoun_cluster$plsg== 'sg',])
summary(pro_rt_sg_cluster)

#Less BIG
pro_rt_cluster.hC <- lmer(resposetime ~ match*clust.hC*plurality + (1+match*clust.hC*plurality|sID) + (1+match*clust.hC*plurality|item), data = pronoun_cluster[pronoun_cluster$accuracy == 1,])
summary(pro_rt_cluster.hC)


##What about a model to see if individual participants improved throughout the task? 
#(lmer using trial order as a predictor?) (treat trial number as a continuous variable)
pro_rt_trial_plura2l  <- lmer(resposetime ~ match*cluster*trialorder + (1+match|sID) + (1|item), data = pronoun_cluster[pronoun_cluster$accuracy == 1 & pronoun_cluster$plsg== 'sg',])
summary(pro_rt_trial_plura2l)
#Fixed effects:
#  Estimate Std. Error        df t value Pr(>|t|)    
#(Intercept)                   745.0437    25.7155  182.7107  28.973  < 2e-16 ***
#  matchmis                      220.8139    31.0185  320.6161   7.119 7.18e-12 ***
#  cluster2                      -80.7287    32.1855  183.8171  -2.508 0.012999 *  
#  cluster3                      -71.5523    34.6565  189.7416  -2.065 0.040319 *  
#  trialorder                     -0.9502     0.2664 3705.5087  -3.567 0.000366 ***
#  matchmis:cluster2              51.8031    39.0381  329.7397   1.327 0.185430    
#matchmis:cluster3              54.4508    41.9606  337.5482   1.298 0.195289    
#matchmis:trialorder            -1.0425     0.3862 3710.3839  -2.699 0.006986 ** 
#  cluster2:trialorder             0.4430     0.3364 3711.8353   1.317 0.187937    
#cluster3:trialorder             0.2363     0.3649 3711.6911   0.648 0.517292    
#matchmis:cluster2:trialorder   -0.5893     0.4861 3698.4917  -1.212 0.225536    
#matchmis:cluster3:trialorder   -1.0519     0.5219 3707.4115  -2.016 0.043923 *  

#Large lmer with trial order as a response variable (also using maximal effects)
#Model is good optimizer
pro_rt_trial1 <- lmer(resposetime ~ match*plsg*trialorder + (1+match*plsg*trialorder|sID) + (1+match*plsg*trialorder|item), data = pronoun_cluster[pronoun_cluster$accuracy == 1,], control = lmerControl(optimizer = "nloptwrap",optCtrl = list(algorithm = "NLOPT_LN_BOBYQA", maxit = 1e6)))
summary(pro_rt_trial1)
#Fixed effects:
#                           Estimate Std. Error         df   t value  #Pr(>|t|)    
#(Intercept)                 685.97866   13.70590   75.80417  50.050 < 2e-16 ***
#matchmis                    195.57718   17.14083  105.84180  11.410 < 2e-16 ***
#plsgsg                        1.51648   12.05789  505.88275   0.126    0.9000  
#trialorder                   -0.72095    0.12766  295.40378  -5.647 3.83e-08 *** (shows that people are getting faster over course of experiment)
#matchmis:plsgsg              61.27572   27.29366  143.89840   2.245  0.0263 *
#matchmis:trialorder          -1.43278    0.21419  124.36378  -6.689 6.81e-10 ***
#plsgsg:trialorder             0.03444    0.17092 2101.32020   0.201 0.8403 
#matchmis:plsgsg:trialorder   -0.13705    0.34854  126.60301  -0.3930.6948

pro_rt_trial_plural <- pro_rt_trial1 <- lmer(resposetime ~ match*trialorder + (1+match*trialorder|sID) + (1+match*trialorder|item), data = pronoun[pronoun$accuracy == 1 & pronoun$plsg == "pl",], control = lmerControl(optimizer = "nloptwrap"))
diff_optims <- allFit(pro_rt_trial_plural, maxfun = 1e5, parallel = 'multicore', ncpus = detectCores())
is.OK <- sapply(diff_optims, is, "merMod")
diff_optims.OK <- diff_optims[is.OK]
lapply(diff_optims.OK,function(x) x@optinfo$conv$lme4$messages)



#as trial order is growing, the negative effect means that match vs. mismatch effect is getting smaller
#should probably graph
#this means that, matchmis:trialorder, tells us that as match gets bigger,  
#matchmis: is positive because, the mismatch is positive, match is negative. This makes sense: intuitively, msimatch should have higher RT

#two trend lines: 1) match condition and misatch condition reaction as a function of trial order 
#item order on x ais, reaction times on y axis
#each trend line is a difference condition
library(ggplot2)
ggplot(pronoun, aes(trialorder, resposetime)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  geom_smooth(aes(y = mismat == "match" & ), method = "lm", se = FALSE, color = "red") +
  geom_smooth(aes(y = y3), method = "lm", se = FALSE, color = "green") +
  geom_smooth(aes(y = y4), method = "lm", se = FALSE, color = "purple") +
  labs(title = "Scatterplot with Four Trend Lines", x = "X-axis", y = "Y-axis") +
  theme_bw()

#follow-up analysis:
pro_rt_trial2_singular <- lmer(resposetime ~ match*trialorder + (1+match*trialorder|sID) + (1+match*trialorder|item), data = pronoun_cluster[pronoun_cluster$accuracy == 1 & pronoun_cluster$plsg == 'sg',], control = lmerControl(optimizer = "nloptwrap",optCtrl = list(algorithm = "NLOPT_LN_BOBYQA", maxit = 1e6)))
summary(pro_rt_trial2_singular)
#Fixed effects:
#Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)         688.2510    14.7270  72.0944  46.734  < 2e-16 ***
#  matchmis            255.1373    27.0920  80.0891   9.417 1.32e-14 ***
#  trialorder           -0.6932     0.1384 231.5910  -5.008 1.09e-06 ***
#  matchmis:trialorder  -1.5598     0.3151  79.0656  -4.950 4.10e-06 ***

#seeing which is the best optimizer:
diff_optims <- allFit(pro_rt_trial1, maxfun = 1e5, parallel = 'multicore', ncpus = detectCores())
is.OK <- sapply(diff_optims, is, "merMod")
diff_optims.OK <- diff_optims[is.OK]
lapply(diff_optims.OK,function(x) x@optinfo$conv$lme4$messages)


pro_rt_trial <- lmer(resposetime ~ match*plsg*trialorder + (1+match|sID) + (1|item), data = pronoun_cluster[pronoun_cluster$accuracy == 1,])
summary(pro_rt_trial)
#Fixed effects:
#  Estimate Std. Error         df t value Pr(>|t|)    
#(Intercept)                 684.95812   13.11149  224.83732  52.241  < 2e-16 ***
#  matchmis                    194.22443   13.97435  778.40895  13.899  < 2e-16 ***
#  plsgsg                        2.37362   12.32553 7440.08423   0.193 0.847295    
#trialorder                   -0.71174    0.12767 7423.83522  -5.575 2.57e-08 ***
#  matchmis:plsgsg              65.28527   17.54904 7423.69846   3.720 0.000201 ***
#  matchmis:trialorder          -1.42235    0.18049 7416.09513  -7.881 3.73e-15 ***
#  plsgsg:trialorder             0.02242    0.17808 7424.47960   0.126 0.899798    
#matchmis:plsgsg:trialorder   -0.17979    0.25308 7430.33057  -0.710 0.477468    
df_pro_rt_trial  <- data.frame(
  effect = rownames(summary(pro_rt_trial)$coefficients),
  estimate = summary(pro_rt_trial)$coefficients[, "Estimate"],
  se = summary(pro_rt_trial)$coefficients[, "Std. Error"],
  df = summary(pro_rt_trial)$coefficients[, "df"],
  t_value = summary(pro_rt_trial)$coefficients[, "t value"],
  pvalue = summary(pro_rt_trial)$coefficients[, "Pr(>|t|)"]
)
View(df_pro_rt_trial)

#checking best fit so that we can make a maximal model: 
diff_optims <- allFit(pro_rt_trial, maxfun = 1e5, parallel = 'multicore', ncpus = detectCores())
is.OK <- sapply(diff_optims, is, "merMod")
diff_optims.OK <- diff_optims[is.OK]
lapply(diff_optims.OK,function(x) x@optinfo$conv$lme4$messages)

#Maximal model:
pro_rt_trial_plural <- lmer(resposetime ~ match*trialorder + (1+match*trialorder|sID) + (1+match*trialorder|item), data = pronoun[pronoun$accuracy == 1 & pronoun$plsg== 'pl',])

pro_rt_trial_singular <- lmer(resposetime ~ match*trialorder + (1+match|sID) + (1+match|item), data = pronoun[pronoun$accuracy == 1 & pronoun$plsg== 'sg',])
#Fixed effects:
#  Estimate Std. Error        df t value Pr(>|t|)    
#(Intercept)          688.3565    13.0439  188.6001  52.772  < 2e-16 ***
#matchmis             260.6883    16.3558  358.6215  15.939  < 2e-16 ***
#  trialorder            -0.6937     0.1341 3626.5245  -5.172 2.44e-07 ***
#  matchmis:trialorder   -1.6409     0.1938 3689.9033  -8.469  < 2e-16 ***



#_______creating a scatterplot for the relationship between trial order and response time----
# filter data into four groups based on the conditions you want
group1 <- pronoun %>% filter(match == "mis" & plsg == "sg")
group2 <- pronoun %>% filter(match == "mat" & plsg == "pl")
group3 <- pronoun %>% filter(match == "mis" & plsg == "pl")
group4 <- pronoun %>% filter(match == "mat" & plsg == "sg")

# create a scatterplot with a trend line for each group
ggplot() +
  geom_smooth(data = group1, aes(x = trialorder, y = resposetime), method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
  geom_smooth(data = group2, aes(x = trialorder, y = resposetime), method = "lm", se = FALSE, color = "green", linetype = "solid") +
  geom_smooth(data = group3, aes(x = trialorder, y = resposetime), method = "lm", se = FALSE, color = "green", linetype = "dashed") +
  geom_smooth(data = group4, aes(x = trialorder, y = resposetime), method = "lm", se = FALSE, color = "blue", linetype = "solid") +
  labs(title = "Scatterplot with Four Trend Lines", x = "Trial Order", y = "Response Time") +
  theme_bw()

#anovabysub
anova_bysub_rt <- ddply(pronoun[pronoun$accuracy == 1,], .(sID, match, plsg), summarise,  rt = mean(resposetime, na.rm = TRUE))


# Calculate means and SDs by condition
means_by_cond <- ddply(anova_bysub_rt, c("plsg", "mat", "nb"), summarise,
                       N    = length(rt),
                       mean = meanrt,
                       sd   = sd(rt))

# Calculate CIs using t.test()
cis_by_cond <- ddply(anova_bysub, c("sgpl", "mismat", "nb"), summarise,
                     ci_lower = t.test(rt, conf.level = 0.95)$conf.int[1],
                     ci_upper = t.test(rt, conf.level = 0.95)$conf.int[2])

# Combine means, SDs, and CIs into one data frame
summary_anova_bysub <- merge(means_by_cond, cis_by_cond, by = c("sgpl", "mismat", "nb"))
summary_anova_bysub$se <- summary_anova_bysub$sd / sqrt(summary_anova_bysub$N)
