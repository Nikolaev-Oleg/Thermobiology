#| Description of variables:
#| Date -> date of collection, not used
#| Ter_N -> terrarium id
#| Sp_N -> animal id
#| Sp -> species
#| Time -> time of measurement
#| Time_initial -> time point when the light was switched off
#| Delta time -> Time - Time_initial
#| T_lizard -> body temperature
#| T_sub -> substratum temperature
#| T_air -> air temperature 1 cm above substratum
#| N_group -> number of individuals in a group
#| N_total -> number of individuals in a terrarium
#| N_ser -> series id
#| Sub_type -> type of substratum (rock[brick]/gravel)
#| Group -> is the individual in a group or not (True/False)
#| Season -> data collection before the hibernation (autumn) or after (spring)

#Import packages####

library(readxl)
library(ggplot2)
library(lme4)
library(lmerTest)
library(DrugClust)
library(dplyr)
library(tidyverse)
library(shiny)
library(vegan)
library(data.table)
library(wPerm)
library(cowplot)
library(splancs)
#Load and preprocess data####

data <- read.csv("data.csv")

Basking_a<-subset(data, !(N_group>3)) # Remove the only present group of four individuals
Basking_a<-subset(Basking_a, 
                  ((Sub_type=='gravel')|(Sub_type=='rock'))) # Double-check if strange substratum types are absent
Basking_a<-subset(Basking_a, Sp=='D.armeniaca') # Double-check if all individuals are D. armeniaca

#Remove individuals with few observations####

for (i in Basking_a$Sp_N){
  temp<-subset(Basking_a$Sp_N, Basking_a$Sp_N==i)
  if(length(temp)<5)Basking_a<-subset(Basking_a, !(Sp_N==i))}
for (i in Basking_a$Sp_N){
  temp<-subset(Basking_a, Sp_N==i)
  temp1<-subset(temp, N_group==1)
  temp2<-subset(temp, N_group==2)
  temp3<-subset(temp, N_group==3)
  if(nrow(temp1)<5)Basking_a<-subset(Basking_a, !(Sp_N==i&N_group==1))
  if(nrow(temp2)<5)Basking_a<-subset(Basking_a, !(Sp_N==i&N_group==2))
  if(nrow(temp3)<5)Basking_a<-subset(Basking_a, !(Sp_N==i&N_group==3))
}

#GLMM####
model2<-with(Basking_a, lmer(T_lizard~Sub_type+
                               Season+
                               as.character(N_group)+
                               N_total+
                               as.character(N_group)*T_air+
                               as.character(N_group)*T_sub+
                               as.character(N_group)*log(Delta_Time+1, exp(1))+
                               (1|Sp_N)+
                               (0+T_sub|Sp_N)+
                               (0+T_air|Sp_N)+
                               (0+log(Delta_Time+1, exp(1))|Sp_N), REML=FALSE)
)

# Summarize
summary(model2)
confint(model2)

# Calcuate R-squared
res<-residuals(model2)
res_sq<-res^2
ESS<-sum(res_sq)
mn<-mean(Basking_a$T_lizard)
res_sq<-(Basking_a$T_lizard-mn)^2
TSS<-sum(res_sq)
1-ESS/TSS

#Cross-validation#####

# p-values for each coefficient in every iteration
matrix_of_p_values<-matrix(nrow=1, ncol = 16) 
colnames(matrix_of_p_values)<-c(rownames(summary(model2)[['coefficients']]),'iter') # initiate empty matrix to write to

set.seed(179) # seed to generate 1000 pseudo-random seeds
new.seed<-runif(1000, 0, 10000) # pseudo-random seed for each of 1000 runs of 10-fold cross-validation

for(j in new.seed){
  set.seed(j)
  folds<-CreateFolds(Basking_a, 10) 
  Basking_a_folded<-cbind(Basking_a, folds) # randomly split the data into 10 batches
  
  fold_nums<-c(0:9)
  for(i in fold_nums){
    
    # Estimate model parameters on 9 out of 10 batches
    model_f<-with(subset(Basking_a_folded, !(folds==i)), lmer(T_lizard~
                                                                Sub_type+ 
                                                                Season+
                                                                as.character(N_group)+
                                                                N_total+
                                                                as.character(N_group)*T_air+
                                                                as.character(N_group)*T_sub+
                                                                as.character(N_group)*log(Delta_Time+1, exp(1))+
                                                                (1|Sp_N)+
                                                                (0+T_sub|Sp_N)+
                                                                (0+T_air|Sp_N)+
                                                                (0+log(Delta_Time+1, exp(1))|Sp_N), REML=FALSE)
    )
    
    predict_test<-predict(model_f, newdata=subset(Basking_a_folded, folds==i), allow.new.levels=TRUE) # Use the estimated parameters to predict Tb on the test batch
    
    signif.codes<-vector()
    for (i in summary(model_f)[['coefficients']][,5]){
      if(i>=0.1)signif.codes<-append(signif.codes, '_')
      if(i<0.1&i>=0.05)signif.codes<-append(signif.codes, '.')
      if(i<0.05&i>=0.01)signif.codes<-append(signif.codes, '*')
      if(i<0.01&i>=0.001)signif.codes<-append(signif.codes, '**')
      if(i<0.001)signif.codes<-append(signif.codes, '***')
    } # replace p-values with significance codes to meke the output human-readable
    signif.codes<-append(signif.codes, j) # use a random seed j as the iteration id
    
    matrix_of_p_values<- rbind(matrix_of_p_values, signif.codes) # add iteration info to the summary matrix
  }}

# Calculate R-squared 
matrix_of_R2<-matrix(ncol=3)
colnames(matrix_of_R2)<-c('train', 'test','iter') # initiate empty matrix to write to

for(j in new.seed){
  set.seed(j)
  folds<-CreateFolds(Basking_a, 10) 
  Basking_a_folded<-cbind(Basking_a, folds) # randomly split the data into 10 batches
  fold_nums<-c(0:9)
  
  for(i in fold_nums){
    
    # Estimate model parameters on 9 out of 10 batches
    model_f<-with(subset(Basking_a_folded, !(folds==i)), lmer(T_lizard~
                                                                Sub_type+
                                                                Season+
                                                                as.character(N_group)+
                                                                N_total+
                                                                as.character(N_group)*T_air+
                                                                as.character(N_group)*T_sub+
                                                                as.character(N_group)*log(Delta_Time+1, exp(1))+
                                                                (1|Sp_N)+
                                                                (0+T_sub|Sp_N)+
                                                                (0+T_air|Sp_N)+
                                                                (0+log(Delta_Time+1, exp(1))|Sp_N), REML=FALSE)
    )
    
    predict_test<-predict(model_f, newdata=subset(Basking_a_folded, folds==i), allow.new.levels=TRUE) # Use the estimated parameters to predict Tb on the test batch
    
    R2<-vector() 
    res<-residuals(model_f) 
    res_sq<-res^2
    ESS<-sum(res_sq)
    mn<-mean(subset(Basking_a_folded, !(folds==i))$T_lizard)
    res_sq<-(subset(Basking_a_folded, !(folds==i))$T_lizard-mn)^2
    TSS<-sum(res_sq)
    R2<-append(R2, 1-ESS/TSS) # calculate train-sample R-squared
    
    res<-predict_test-subset(Basking_a_folded, folds==i)$T_lizard
    ESS<-sum(res^2)
    mn<-mean(subset(Basking_a_folded, folds==i)$T_lizard)
    TSS<-sum((subset(Basking_a_folded, folds==i)$T_lizard-mn)^2)
    R2<-append(R2, 1-ESS/TSS) # calculate test-sample R-squared
    
    R2<-append(R2, j) # use seed as an iteration id
    matrix_of_R2<-rbind(matrix_of_R2,R2) # write iteration result to the summary matrix
  }}

matrix_of_R2<-as.data.frame(matrix_of_R2)
matrix_of_R2$delta<-matrix_of_R2$train-matrix_of_R2$test # calculate the difference between train and test R-squared

matrix_of_p_values<-as.data.frame(matrix_of_p_values)
matrix_of_p_values<-cbind(matrix_of_p_values, matrix_of_R2$test) # merge test-sample R-squared and significance codes for this model
matrix_of_p_values<-na.omit(matrix_of_p_values)

#best50<-matrix_of_p_values[matrix_of_p_values$`matrix_of_R2$test`>=median(na.omit(matrix_of_p_values$`matrix_of_R2$test`)),] #select 50% of the best models in order to check if there are different significant parameters in "good" and "bad" models

# Model coefficients for each iteration
model.coefs<-matrix(nrow=1, ncol = 16)
colnames(model.coefs)<-c(rownames(summary(model2)[['coefficients']]),'iter')
set.seed(179)
new.seed<-runif(1000, 0, 10000)
for(j in new.seed){
  set.seed(j)
  folds<-CreateFolds(Basking_a, 10)
  Basking_a_folded<-cbind(Basking_a, folds)
  fold_nums<-c(0:9)
  for(i in fold_nums){
    model_f<-with(subset(Basking_a_folded, !(folds==i)), lmer(T_lizard~
                                                                Sub_type+ 
                                                                Season+
                                                                as.character(N_group)+
                                                                N_total+
                                                                as.character(N_group)*T_air+
                                                                as.character(N_group)*T_sub+
                                                                as.character(N_group)*log(Delta_Time+1, exp(1))+
                                                                (1|Sp_N)+
                                                                (0+T_sub|Sp_N)+
                                                                (0+T_air|Sp_N)+
                                                                (0+log(Delta_Time+1, exp(1))|Sp_N), REML=FALSE)
    )
    
    predict_test<-predict(model_f, newdata=subset(Basking_a_folded, folds==i), allow.new.levels=TRUE)
    signif.codes<-summary(model_f)[['coefficients']][,1]
    signif.codes<-append(signif.codes, j)
    model.coefs<- rbind(model.coefs, signif.codes )
  }}
model.coefs<-na.omit(model.coefs)
model.coefs<-cbind(model.coefs, matrix_of_R2$test)
model.coefs<-as.data.frame(model.coefs)
#MW-ANOVA####
make_binary<-function(x)ifelse(x=='*'|x=='**'|x=='***', 1, 0) # significance codes to binary "significant or not" format

matrix_of_p_values_binary<- matrix_of_p_values[,1:15] %>% 
  mutate_all(make_binary) # transform p-values matrix to binary format

matrix_of_p_values_binary<-matrix_of_p_values_binary[,!(colMeans(matrix_of_p_values_binary)==1|colMeans(matrix_of_p_values_binary)==0)] # remove effects, which are always significant of always non-significant
matrix_of_p_values_binary<- matrix_of_p_values_binary %>% mutate(across(where(is.numeric), as.character)) # transform variables to discrete form
matrix_of_p_values_binary<- cbind(matrix_of_p_values_binary, matrix_of_p_values[,16:17]) # annotate the matrix

colnames(matrix_of_p_values_binary)<-c('(Intercept)',
                                 'Sub_typerock',
                                 'as.character(N_group)2', 
                                 'N_total',
                                 'as.character(N_group)2:T_air',
                                 'as.character(N_group)3:T_air',
                                 'as.character(N_group)2:T_sub',
                                 'as.character(N_group)3:T_sub',
                                 'as.character(N_group)3:log(Delta_Time + 1, exp(1))',
                                 'iter',
                                 'matrix_of_R2$test') # annotate the matrix

model.q.check<-with(matrix_of_p_values_binary,lm(`matrix_of_R2$test`~
                                             0+`(Intercept)`+
                                             `Sub_typeкирпич`+
                                             `as.character(N_group)2`+ 
                                             N_total +
                                             `as.character(N_group)2:T_air`+
                                             `as.character(N_group)3:T_air`+
                                             `as.character(N_group)2:T_sub`+
                                             `as.character(N_group)3:T_sub`+
                                             `as.character(N_group)3:log(Delta_Time + 1, exp(1))`)) 
summary(model.q.check) # estimate the effect of the model parameters on a test sample R-squared

#Calculate mean difference between train and test R-squared and CI####

matrix_of_R2<-matrix_of_R2[!is.na(matrix_of_R2$train),] # remove NA

set.seed(179)
new.seed<-runif(10000, 0, 10000) #here we use bootstapping to estimate mean and CI, new seed is generated for each of 10000 iterations

means_delta<-c() # initiate a vector containing means
for(j in new.seed){
  set.seed(j)
  
  folds<-CreateFolds(matrix_of_R2, 2) 
  delta_folded<-data.frame(delta=matrix_of_R2$delta, folds=folds) # select random 50% of values
  
  means_delta<-append(means_delta, mean(delta_folded[folds==1,]$delta)) # calculate mean delta on the batch and add it to a vector
}

means_delta<-data.frame(mn=means_delta) 

norm<-rnorm(
  n = nrow(means_delta),
  mean = mean(means_delta$mn),
  sd = sd(means_delta$mn)
) # generate normally distributed data for plotting using estimated parameters

means_delta<-cbind(means_delta, norm) 

mean(means_delta$mn)
sd(means_delta$mn) # CI is +/- 1.96*sd(means_delta$mn)
#Fig.1####
Tcloacal <- read.csv("CloacalSurface.csv") # load data

ggplot(Tcloacal, aes(T_clo, T_sur))+geom_smooth(method=lm, colour='grey40')+geom_point(colour='#11ee22', alpha=0.7,size=2.5)+geom_point(shape=1, size=3)+theme(panel.border = element_rect(size=1, fill=NA), panel.background = element_rect(fill=NA), axis.title = element_text(size=16))+xlab('Cloacal body temperature, ºC')+ylab('Surface body temperature, ºC')

#with(Tcloacal, cor.test(T_clo, T_sur, method = 'spearman'))
with(Tcloacal, cor.test(T_clo, T_sur, method = 'pearson'))

with(Tcloacal, lm(T_sur~T_clo))
#Fig.2####
Fig2a<-ggplot(matrix_of_R2)+
  geom_histogram(aes(train), fill='#ff3333', alpha=0.6, bins=100)+
  geom_histogram(aes(test), fill='#1177ff', alpha=0.6, bins = 100)+
  theme_classic()+
  theme(axis.title = element_text(size = 16),
        plot.title = element_text(size=16, face='bold', hjust=.97, vjust = -7))+
  xlab(expression(R^2))+
  ylab('N of iterations')+
  ggtitle('a')

Fig2b<-ggplot(matrix_of_R2)+
  geom_histogram(aes(delta), fill='#33ee77', alpha=0.6, bins=100)+
  theme_classic()+
  theme(axis.title = element_text(size = 16),
        plot.title = element_text(size=16, face='bold', hjust=.97, vjust = -7))+
  xlab(expression(Train ~R^2 - Test ~R^2))+
  ylab('N of iterations')+
  ggtitle('b')

Fig2d<-ggplot(matrix_of_R2, aes(train, test))+
  geom_point(shape = 21
  )+
  #geom_smooth(se=T, colour='#EE5533')+
  theme_classic()+
  theme(axis.title = element_text(size=16),
        plot.title = element_text(size=16, face='bold', hjust=.97, vjust = -7))+
  xlab(expression(Train ~R^2))+
  ylab(expression(Test ~R^2))+
  ggtitle('d')

means_delta$alpha<-'0'
for(i in 1:nrow(means_delta)){
  if(means_delta$mn[i]>mean(means_delta$mn)+1.96*sd(means_delta$mn) |
     means_delta$mn[i]<mean(means_delta$mn)-1.96*sd(means_delta$mn)
  ){
    means_delta$alpha[i]<-'1'
  }
}
Fig2c<- ggplot(means_delta)+
  geom_histogram(aes(mn, fill = alpha),
                 bins = 30,
                 colour = 'gray20')+
  geom_density(aes(norm), adjust =2)+
  theme_classic()+
  theme(axis.title = element_text(size=16),
        legend.position = 'none',
        plot.title = element_text(size=16, face='bold', hjust=.97, vjust = -7))+
  xlab('Mean difference')+
  ylab('Count')+
  scale_fill_manual(values = c('gray70', '#ff3333'))+
  ggtitle('c')

leg.data<-data.frame(val = rep(c(1,0), 5), grp = rep(c('Train', 'Test'), 5))
leg.plot<-ggplot(leg.data, aes(val, fill=grp))+
  geom_histogram()+
  theme_classic()+
  theme(legend.title = element_blank(), legend.text = element_text(size=14))+
  scale_fill_manual(values = c('#1177ff', '#ff3333'))
leg<-get_legend(leg.plot)

plot_grid(Fig2a, Fig2b, Fig2c, Fig2d, align = 'hv')+
  draw_plot(leg, x=0.02, y=0.8, width=1/6, height=0.2)
#Fig.6####
matrix_of_p_values_binary[matrix_of_p_values_binary == 1]<-'Significant'
matrix_of_p_values_binary[matrix_of_p_values_binary == '0']<-'Non-significant'
matrix_of_p_values_binary<-cbind(matrix_of_p_values_binary, matrix_of_R2)

Sup1.a<-ggplot(matrix_of_p_values_binary[,-1], 
               aes(train, test,
                   colour=as.character(X.Intercept.)))+
  geom_point(shape = 21)+
  theme_classic()+
  theme(axis.title = element_text(size=16), legend.title = element_blank(),
        plot.title = element_text(size=16, face='bold', hjust=.97, vjust = -7))+
  xlab('')+
  ylab(expression(Test ~R^2))+
  ggtitle('a')+
  guides(colour = 'none')

Sup1.b<-ggplot(matrix_of_p_values_binary[,-1], 
               aes(train, test,
                   colour=as.character(Sub_typeкирпич)))+
  geom_point(shape = 21)+
  theme_classic()+
  theme(axis.title = element_text(size=16), legend.title = element_blank(),
        plot.title = element_text(size=16, face='bold', hjust=.97, vjust = -7))+
  xlab('')+
  ylab('')+
  ggtitle('b')+
  guides(colour = 'none')

Sup1.c<-ggplot(matrix_of_p_values_binary[,-1], 
               aes(train, test,
                   colour=as.character(as.character.N_group.2)))+
  geom_point(shape = 21)+
  theme_classic()+
  theme(axis.title = element_text(size=16), legend.title = element_blank(),
        plot.title = element_text(size=16, face='bold', hjust=.97, vjust = -7))+
  xlab('')+
  ylab('')+
  ggtitle('c')+
  guides(colour = 'none')

Sup1.d<-ggplot(matrix_of_p_values_binary[,-1], 
               aes(train, test,
                   colour=as.character(N_total)))+
  geom_point(shape = 21)+
  theme_classic()+
  theme(axis.title = element_text(size=16), legend.title = element_blank(),
        plot.title = element_text(size=16, face='bold', hjust=.97, vjust = -7))+
  xlab('')+
  ylab(expression(Test ~R^2))+
  ggtitle('d')+
  guides(colour = 'none')

Sup1.e<-ggplot(matrix_of_p_values_binary[,-1], 
               aes(train, test,
                   colour=as.character(as.character.N_group.2.T_air)))+
  geom_point(shape = 21)+
  theme_classic()+
  theme(axis.title = element_text(size=16), legend.title = element_blank(),
        plot.title = element_text(size=16, face='bold', hjust=.97, vjust = -7))+
  xlab('')+
  ylab('')+
  ggtitle('e')+
  guides(colour = 'none')

Sup1.f<-ggplot(matrix_of_p_values_binary[,-1], 
               aes(train, test,
                   colour=as.character(as.character.N_group.3.T_air)))+
  geom_point(shape = 21)+
  theme_classic()+
  theme(axis.title = element_text(size=16), legend.title = element_blank(),
        plot.title = element_text(size=16, face='bold', hjust=.97, vjust = -7),
        legend.text = element_text(size=14))+
  xlab('')+
  ylab('')+
  ggtitle('f')+
  guides(color = guide_legend(override.aes = list(size = 5)))

leg<-get_legend(Sup1.f)
Sup1.f<-Sup1.f+guides(colour = 'none')

Sup1.g<-ggplot(matrix_of_p_values_binary[,-1], 
               aes(train, test,
                   colour=as.character(as.character.N_group.2.T_sub)))+
  geom_point(shape = 21)+
  theme_classic()+
  theme(axis.title = element_text(size=16), legend.title = element_blank(),
        plot.title = element_text(size=16, face='bold', hjust=.97, vjust = -7))+
  xlab(expression(Train ~R^2))+
  ylab(expression(Test ~R^2))+
  ggtitle('g')+
  guides(colour = 'none')

Sup1.h<-ggplot(matrix_of_p_values_binary[,-1], #!
               aes(train, test,
                   colour=as.character(as.character.N_group.3.T_sub)))+
  geom_point(shape = 21)+
  theme_classic()+
  theme(axis.title = element_text(size=16), legend.title = element_blank(),
        plot.title = element_text(size=16, face='bold', hjust=.97, vjust = -7))+
  xlab(expression(Train ~R^2))+
  ylab('')+
  ggtitle('h')+
  guides(colour = 'none')

Sup1.i<-ggplot(matrix_of_p_values_binary[,-1], #!
               aes(train, test,
                   colour=as.character(as.character.N_group.3.log.Delta_Time...1..exp.1..)))+
  geom_point(shape = 21)+
  theme_classic()+
  theme(axis.title = element_text(size=16), legend.title = element_blank(),
        plot.title = element_text(size=16, face='bold', hjust=.97, vjust = -7))+
  xlab(expression(Train ~R^2))+
  ylab('')+
  ggtitle('i')+
  guides(colour = 'none')

plot_grid(Sup1.a, Sup1.b, Sup1.c, NULL, 
          Sup1.d, Sup1.e, Sup1.f, NULL,
          Sup1.g, Sup1.h, Sup1.i, NULL,
          align = 'hv',
          rel_widths = rep(c(2,2,2, 1), 3),
          nrow = 3,
          ncol=4)+
  draw_plot(leg, x=6/7, y=1/3, width=1/7, height = 1/3)
#Table SI1####
cnames<-colnames(model.coefs)
cnames<-paste0(cnames, '_COEF')
colnames(model.coefs)<-cnames

make_binary<-function(x)ifelse(x=='*'|x=='**'|x=='***', 1, 0)
make_binary2<-function(x)ifelse(x>0, 1, 0)

matrix_of_p_values<-matrix_of_p_values[,-1]
matrix_of_p_values<- matrix_of_p_values[,1:15] %>% mutate_all(make_binary)

model.coefs<- model.coefs[,2:16] %>% mutate_all(make_binary2)

signif.pos <- rep(0, 15)
signif.neg <- rep(0, 15)
nonsignif.pos <- rep(0, 15)
nonsignif.neg <- rep(0, 15)

for(i in 1:nrow(model.coefs)){
  for(j in 1:length(model.coefs)){
    if(model.coefs[i,j] == 1 & matrix_of_p_values[i,j] == 1)signif.pos[j]<-signif.pos[j]+1
    ##//////////////##
    else if(model.coefs[i,j] == 0 & matrix_of_p_values[i,j] == 1)signif.neg[j]<-signif.neg[j]+1
    ##/////////////##
    if(model.coefs[i,j] == 1 & matrix_of_p_values[i,j] == 0)nonsignif.pos[j]<-nonsignif.pos[j]+1
    ##//////////////##
    else if(model.coefs[i,j] == 0 & matrix_of_p_values[i,j] == 0)nonsignif.neg[j]<-nonsignif.neg[j]+1
  }
}

tabS1<-rbind(signif.pos, signif.neg, nonsignif.pos, nonsignif.neg)