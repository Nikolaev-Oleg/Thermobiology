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

#| Some notes on series ids:
#| Series of measurements were not initially assigned with consequent ids. 
#| In addition, some series were removed as described in the code below. 
#| So, for consistency we renamed series in the main text:
#| [id in data.csv] -> [id in the main text]
#| i                -> i
#| ii               -> ii
#| iv               -> iii
#| v                -> iv
#| vi               -> v
#| vii              -> vi
#| ix               -> vii
#| xv               -> viii

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
library(ggpubr)
library(magick)
library(cowplot)
library(strucchange)
library(ggpp)
library(rstatix)
library(emmeans)

#Load data####
data <- read.csv("data.csv")

Basking_a_ser<-subset(data, !(is.na(N_ser))) # Extract series of consequent measurments
Basking_a_ser<-subset(Basking_a_ser, Sp=='D.armeniaca') # Check if data on other species is present

Basking_a_ser$N_ser[Basking_a_ser$N_ser=='x']<-'ix' # Trese are actually two distinct groups within the same series of measurements

Basking_a_ser$Dif<-Basking_a_ser$T_lizard-Basking_a_ser$T_sub
Basking_a_ser<-subset(Basking_a_ser, Dif>=-0.2) # Exclude measurements where Tb<Ts. Bias of 0.2ºC was included due to +/- 0.1ºC precision of the thermal imager


#Select series of measurements####
list_of_series<-list() # initiate empty list to write to

#fill the list with data subsets, containing measurements
for(i in Basking_a_ser$N_ser){
  list_of_series<-append(list_of_series, 
                         list(assign(paste0('Basking_a_ser_',i), subset(Basking_a_ser, N_ser==i)))) # name each element iteratively
}

# exclude repeated elements from the list
for (k in 1:length(list_of_series)){
  i<- 1
  while(i<=length(list_of_series)){
    j<- 1
    while (j <= length(list_of_series)&i<=length(list_of_series)){
      if(!(j==i)&list_of_series[[j]]$N_ser[1]==list_of_series[[i]]$N_ser[1]){list_of_series<-list_of_series[-j]
      }
      j<-j+1
    }
    i<- i+1}
}

#name elements of a list
names<-c() # initiate an empty vector to write to
for(i in list_of_series){
  names<-append(names,i$N_ser[1])
} # fill the vector
names(list_of_series)<-names

#exclude individuals with few observations
for(i in list_of_series){
  for(j in i$Sp_N){
    k<-subset(i, Sp_N==j)
    #print(nrow(k))
    if(nrow(k)<5){i<-subset(i, !(Sp_N==j))}
  }
}

#Plot Dif ~ log(deltaTime) for all series and store the plots in a list####

plotlist<-list() #initiate the list

for (i in list_of_series){
  assign(paste0('plot_',i$N_ser[1]), # save each plot as a separate variable
         (ggplot(i, aes(log(Delta_Time+1, exp(1)), Dif,    # Dif is Tb - Tsub
                        colour=Sp_N,                       # colour by individual
                        shape=as.character(N_group),       # shape by group size
                        linetype=as.character(N_group))))+ # line shape by group size
           geom_point()+
           geom_smooth(method='lm', se=FALSE)) # add linear regression
  
  plotlist<-append(plotlist, # add a plot to a list
                   list(( ggplot(i, aes(log(Delta_Time+1, exp(1)), Dif,    # Dif is Tb - Tsub
                                        colour=Sp_N,                       # colour by individual
                                        shape=as.character(N_group),       # shape by group size
                                        linetype=as.character(N_group))))+ # line shape by group size
                          geom_point()+
                          geom_smooth(method='lm', se=FALSE))) # add linear regression
}
names(plotlist)<-names # name each plot in the list
#Remove inappropriate series####
#| Exclude series iii: the only solitary lizard is heating 
#| Exclude series viii and xi: observations of lizards in a group of 2 only 

plotlist<-plotlist[names(plotlist) %in% c('viii','xi','iii')==FALSE] # remove from the list of plots
list_of_series<-list_of_series[names(list_of_series) %in% c('viii','xi','iii')==FALSE] # remove from the list of datasets

#Exclude animal with less than 5 observations within each group size class
for (i in 1:length(list_of_series)){
  j<-1
  while (j <=length(list_of_series[[i]]$Sp_N)){ # iterate over animals' ids
    k<-1
    while(k <= length(list_of_series[[i]]$N_group) & j <= length(list_of_series[[i]]$Sp_N)){ # iterate over group sizes within each id
      temp<-subset(list_of_series[[i]],
                   (Sp_N == list_of_series[[i]]$Sp_N[j] & N_group == list_of_series[[i]]$N_group[k]))
      if(nrow(temp)<5){
        list_of_series[[i]] <- subset(list_of_series[[i]], 
                                      !(Sp_N == list_of_series[[i]]$Sp_N[j] & N_group == list_of_series[[i]]$N_group[k]))
        }
      k<-k+1
    }
    j<-j+1}
}

names<-c()
for(i in list_of_series){
  names<-append(names,i$N_ser[1])
} # update names

#Update the list of plots in accordance with the updated list of series
plotlist<-list()
for (i in list_of_series){
  assign(paste0('plot_',i$N_ser[1]), # save each plot as a separate variable
         (ggplot(i, aes(log(Delta_Time+1, exp(1)), Dif,    # Dif is Tb - Tsub
                        colour=Sp_N,                       # colour by individual
                        shape=as.character(N_group),       # shape by group size
                        linetype=as.character(N_group))))+ # line shape by group size
           geom_point()+
           geom_smooth(method='lm', se=FALSE)) # add linear regression line
  plotlist<-append(plotlist, # add a plot to a list
                   list(( ggplot(i, aes(log(Delta_Time+1, exp(1)), Dif,              # Dif is Tb - Tsub
                                                  colour=Sp_N,                       # colour by individual
                                                  shape=as.character(N_group),       # shape by group size
                                                  linetype=as.character(N_group))))+ # line shape by group size
                                    geom_point()+
                                    geom_smooth(method='lm', se=FALSE))) # add linear regression line
}
names(plotlist)<-names
#Plot Tb ~ log(deltaTime) for all series and store the plots in a list####
plotlist2<-list()
for (i in list_of_series){
  plotlist2 <- append(plotlist2, 
                      list(( ggplot(i, aes(log(Delta_Time+1, exp(1)), T_lizard, # this time use non-adjusted Tb
                                           colour=Sp_N,                         # colour by individual
                                           shape=as.character(N_group),         # shape by group size
                                           linetype=as.character(N_group))))+   # line shape by group size
                             geom_point()+
                             geom_smooth(method='lm', se=FALSE))) # add linear regression line
}
names(plotlist2)<-names
#Regressions Tb ~ log(deltaTime)####
list_of_ci_intercept<-matrix(ncol=2) # initiate a matrix to write CIs for intercept to
list_of_ci_Time<-matrix(ncol=2) # initiate a matrix to write CIs for slope coefficient to

a<-c(0) # a vector containing intercepts
b<-c(0) # a vector containing slopes
signs<-c(0) # a vector containing a series id, an animal id and a number of individuals in a group
for (i in list_of_series){ # iterate over the series
  for(j in 1:nrow(i)){
    model<-with(subset(i, N_group == i$N_group[j] & Sp_N==i$Sp_N[j]), # iteratively subset data on each individual within each group (group composition may change)
                lm(T_lizard~log(1+Delta_Time, exp(1)))) # run linear regression
    assign(paste('lm', i$N_ser[1], i$Sp_N[j],i$N_group[j], sep='_'), model) # write the model to a variable
    ci<-confint(model) # calculate CI for the model coefficients
    
    signs<-append(signs, paste('lm', i$N_ser[1], i$Sp_N[j],i$N_group[j], sep='_')) #add the model name to a vector
    a<-append(a, coefficients(model)[1]) # add the intercept to a vector
    b<-append(b, coefficients(model)[2]) # add the slope to a vector
    
    list_of_ci_intercept<-rbind(list_of_ci_intercept, ci[1,]) # write the intercept CI to a vector
    list_of_ci_Time<-rbind(list_of_ci_Time, ci[2,]) # write the slope CI to a vector
  }
}

list_of_ci_intercept<-cbind(list_of_ci_intercept, signs) # add model ids to a matrix of intercepts
list_of_ci_Time<-cbind(list_of_ci_Time, signs) # add model ids to a matrix of slopes

list_of_ci_intercept<-cbind(list_of_ci_intercept, a) # merge intercepts and CIs
list_of_ci_Time<-cbind(list_of_ci_Time, b) # merge slopes and CIs

list_of_ci_intercept<-list_of_ci_intercept[!duplicated(list_of_ci_intercept),] # remove duplicates
list_of_ci_Time<-list_of_ci_Time[!duplicated(list_of_ci_Time),] #remove duplicates

colnames(list_of_ci_intercept)<-c('lower', 'upper', 'identifier', 'intercept') # set readable variable names
colnames(list_of_ci_Time)<-c('lower', 'upper', 'identifier', 'Time') # set readable variable names

list_of_ci_intercept<-as.data.frame(list_of_ci_intercept) # convert to a data.frame class
list_of_ci_Time<-as.data.frame(list_of_ci_Time) # convert to a data.frame class

# Split model ids into a series id + an animal id + number of individuals in a group
for (i in 1:nrow(list_of_ci_intercept)){
  splitted<- strsplit(list_of_ci_intercept$identifier[i], split='_') # split the model id
  list_of_ci_intercept$Series[i]<-splitted[[1]][2] # series id to Series variable
  list_of_ci_intercept$ID[i]<-splitted[[1]][3] # animal id to ID variable
  list_of_ci_intercept$N_group[i]<-splitted[[1]][4] # number of individuals in a group to N_group variable
}

for (i in 1:nrow(list_of_ci_Time)){
  splitted<- strsplit(list_of_ci_Time$identifier[i], split='_') # split the model id
  list_of_ci_Time$Series[i]<-splitted[[1]][2] # series id to Series variable
  list_of_ci_Time$ID[i]<-splitted[[1]][3] # animal id to ID variable
  list_of_ci_Time$N_group[i]<-splitted[[1]][4] # number of individuals in a group to N_group variable
}

list_of_ci_intercept<-list_of_ci_intercept[-1,] # remove the first row containing NAs 
list_of_ci_Time<-list_of_ci_Time[-1,] # remove the first row containing NAs

list_of_ci_intercept$lower<-as.numeric(list_of_ci_intercept$lower) # convert CI bonds to numeric format
list_of_ci_intercept$upper<-as.numeric(list_of_ci_intercept$upper) # convert CI bonds to numeric format

list_of_ci_Time$lower<-as.numeric(list_of_ci_Time$lower) # convert CI bonds to numeric format
list_of_ci_Time$upper<-as.numeric(list_of_ci_Time$upper) # convert CI bonds to numeric format

list_of_ci_intercept$intercept<-as.numeric(list_of_ci_intercept$intercept) # convert intercept to numeric format
list_of_ci_Time$Time<-as.numeric(list_of_ci_Time$Time) # convert slope to numeric format

list_of_ci_intercept_Tb<-list_of_ci_intercept
list_of_ci_Time_Tb<-list_of_ci_Time

# plot regressions coefficients
for(i in names(list_of_series)){
  assign(paste0('ci_noadj_plot_', i), # write a plot to a variable
         ggplot(list_of_ci_Time_Tb[list_of_ci_Time_Tb$Series == i,], 
                aes(ID, Time,        # plot slope coefficients vs id and group size
                    colour=ID,       # colour by animal id
                    fill=ID,   
                    group=N_group,   # group by groups
                    shape=N_group))+ # point shapes by groups
           geom_errorbar(aes(ymin=lower, ymax=upper), 
                         width=0.2, 
                         position=position_dodge(width=.9))+ #plot CI as an errorbar
           geom_point(size=2.5, alpha=0.7, position=position_dodge(width=.9))+ #plot estimated slope coefficient as a point
           theme_classic()+
           ggtitle(i) # add title
  )
}

#Regressions Dif ~ log(deltaTime)####
#| Dif is Tb - Tsub, as described earlier

list_of_ci_intercept<-matrix(ncol=2) # initiate a matrix to write CIs for intercept to
list_of_ci_Time<-matrix(ncol=2) # initiate a matrix to write CIs for slope coefficient to

a<-c(0) # a vector containing intercepts
b<-c(0) # a vector containing slopes
signs<-c(0) # a vector containing a series id, an animal id and a number of individuals in a group

for (i in list_of_series){ # iterate over the series 
  for(j in 1:nrow(i)){ # iterate over the individuals
    model<-with(subset(i, N_group ==i$N_group[j] & Sp_N == i$Sp_N[j]), # iteratively subset data on each individual within each group (group composition may change)
                lm(Dif~log(1+Delta_Time, exp(1)))) # run linear regression
    assign(paste('lm', i$N_ser[1], i$Sp_N[j],i$N_group[j], sep='_'), model) # write the model to a variable
    
    ci<-confint(model) # calculate CI for the model coefficients
    
    signs<-append(signs, paste('lm', i$N_ser[1], i$Sp_N[j],i$N_group[j], sep='_')) #add the model name to a vector
    a<-append(a, coefficients(model)[1]) # add the intercept to a vector
    b<-append(b, coefficients(model)[2]) # add the slope to a vector
    
    list_of_ci_intercept<-rbind(list_of_ci_intercept, ci[1,]) # write the intercept CI to a vector
    list_of_ci_Time<-rbind(list_of_ci_Time, ci[2,]) # write the slope CI to a vector
  }
}

list_of_ci_intercept<-cbind(list_of_ci_intercept, signs) # add model ids to a matrix of intercepts
list_of_ci_Time<-cbind(list_of_ci_Time, signs) # add model ids to a matrix of slopes

list_of_ci_intercept<-cbind(list_of_ci_intercept, a) # merge intercepts and CIs
list_of_ci_Time<-cbind(list_of_ci_Time, b) # merge slopes and CIs

list_of_ci_intercept<-list_of_ci_intercept[!duplicated(list_of_ci_intercept),] # remove duplicates
list_of_ci_Time<-list_of_ci_Time[!duplicated(list_of_ci_Time),] #remove duplicates

colnames(list_of_ci_intercept)<-c('lower', 'upper', 'identifier', 'intercept') # set readable variable names
colnames(list_of_ci_Time)<-c('lower', 'upper', 'identifier', 'Time') # set readable variable names

list_of_ci_intercept<-as.data.frame(list_of_ci_intercept) # convert to a data.frame object
list_of_ci_Time<-as.data.frame(list_of_ci_Time) # convert to a data.frame object

for (i in 1:nrow(list_of_ci_intercept)){
  splitted<- strsplit(list_of_ci_intercept$identifier[i], split='_') # split the model id
  
  list_of_ci_intercept$Series[i]<-splitted[[1]][2] # series id to Series variable
  list_of_ci_intercept$ID[i]<-splitted[[1]][3] # animal id to ID variable
  list_of_ci_intercept$N_group[i]<-splitted[[1]][4] # number of individuals in a group to N_group variable
}

for (i in 1:nrow(list_of_ci_Time)){
  splitted<- strsplit(list_of_ci_Time$identifier[i], split='_') # split the model id
  
  list_of_ci_Time$Series[i]<-splitted[[1]][2] # series id to Series variable
  list_of_ci_Time$ID[i]<-splitted[[1]][3] # animal id to ID variable
  list_of_ci_Time$N_group[i]<-splitted[[1]][4] # number of individuals in a group to N_group variable
}

list_of_ci_intercept<-list_of_ci_intercept[-1,] # remove the first row containing NAs
list_of_ci_Time<-list_of_ci_Time[-1,] # remove the first row containing NAs

list_of_ci_intercept$lower<-as.numeric(list_of_ci_intercept$lower) # convert intercept CI bonds to numeric format
list_of_ci_intercept$upper<-as.numeric(list_of_ci_intercept$upper) # convert intercept CI bonds to numeric format

list_of_ci_Time$lower<-as.numeric(list_of_ci_Time$lower) # convert slope CI bonds to numeric format
list_of_ci_Time$upper<-as.numeric(list_of_ci_Time$upper) # convert slope CI bonds to numeric format

list_of_ci_intercept$intercept<-as.numeric(list_of_ci_intercept$intercept) # convert intercept to numeric format
list_of_ci_Time$Time<-as.numeric(list_of_ci_Time$Time) # convert slope to numeric format

list_of_ci_intercept_Dif<-list_of_ci_intercept
list_of_ci_Time_Dif<-list_of_ci_Time

# plot regression coefficients
for(i in names(list_of_series)){
  assign(paste0('ci_adj_plot_', i), # write a plot to a variable
         ggplot(list_of_ci_Time_Dif[list_of_ci_Time_Dif$Series==i,], 
                aes(ID, Time, # plot slope coefficients vs id and group size
                    colour=ID, # colour by animal id
                    fill=ID, 
                    group=N_group, # grouping by groups
                    shape=N_group))+ # point shapes by groups
           geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2, position=position_dodge(width=.9))+ #plot CI as errorbars
           geom_point(size=2.5, alpha=0.7, position=position_dodge(width=.9))+ # plot estimated slope as a point
           theme_classic()+
           ggtitle(i) #add title
  )
}


#Nested ANOVA####
TukeyHSD(aov(with(list_of_ci_Time_Dif, lm(Time~N_group/Series)))) # Tukey test for the effect of a group size on the adjusted cooling rate
TukeyHSD(aov(with(list_of_ci_Time_Tb, lm(Time~N_group/Series)))) # Tukey test for the effect of a group size on the non-adjusted cooling rate

df<-cbind(list_of_ci_intercept_Tb, list_of_ci_Time_Tb) # merge data frames with intercepts and slopes (non-adjusted)
model<-lmer(data = df,
            Time~as.character(N_group)+intercept+(1|ID/Series)) # test if intercept and slope are correlated (non-adjusted)
summary(model)

emm.model<-emmeans(model, ~as.character(N_group)) # get marginal means for differen N_group levels (non-adjusted)
contrast(emm.model, interaction = "pairwise", adjust = 'holm') # Kenward-Roger test for comparison of slopes (i.e. non-adjusted cooling rates) between groups of different size

df<-cbind(list_of_ci_intercept_Dif, list_of_ci_Time_Dif) # merge data frames with intercepts and slopes (adjusted)
model<-lmer(data = df,
            Time~as.character(N_group)+intercept+(1|ID/Series)) # test if intercept and slope are correlated (adjusted)
summary(model)

emm.model<-emmeans(model, ~as.character(N_group)) # get marginal means for differen N_group levels (non-adjusted)
contrast(emm.model, interaction = "pairwise", adjust = 'holm') # Kenward-Roger test for comparison of slopes (i.e. adjusted cooling rates) between groups of different size
#Compare Tb of animals within groups####

df<-subset(Basking_a_ser_i, N_group > 1) # select data on huddling animals from the series i
df$logdt<-log(df$Delta_Time+1, exp(1)) # log-transform time difference
model<-lm(data = df,
          T_lizard~Sp_N+logdt+Sp_N:logdt) # run linear model
emm.model<-emmeans(model, ~Sp_N) # calculate marginal means
contrast(emm.model, interaction = "pairwise", adjust = 'holm') # Tukey test on marginal means (basically a KR-test, but without random effects)

df<-subset(Basking_a_ser_ii, N_group > 1) # select data on huddling animals from the series ii
df$logdt<-log(df$Delta_Time+1, exp(1)) # log-transform time difference
model<-lm(data = df,
          T_lizard~Sp_N+logdt+Sp_N:logdt) #run linear model
emm.model<-emmeans(model, ~Sp_N) # calculate marginal means
contrast(emm.model, interaction = "pairwise", adjust = 'holm') # Tukey test on marginal means (basically a KR-test, but without random effects)

df<-subset(Basking_a_ser_iv, N_group > 1) # select data on huddling animals from the series iv
df$logdt<-log(df$Delta_Time+1, exp(1)) # log-transform time difference
model<-lm(data = df,
          T_lizard~Sp_N+logdt+Sp_N:logdt) #run linear model
emm.model<-emmeans(model, ~Sp_N) # calculate marginal means
contrast(emm.model, interaction = "pairwise", adjust = 'holm') # Tukey test on marginal means (basically a KR-test, but without random effects)

df<-subset(Basking_a_ser_v, N_group > 1) # select data on huddling animals from the series v
df$logdt<-log(df$Delta_Time+1, exp(1)) # log-transform time difference
model<-lm(data = df,
          T_lizard~Sp_N+logdt+Sp_N:logdt) #run linear model
emm.model<-emmeans(model, ~Sp_N) # calculate marginal means
contrast(emm.model, interaction = "pairwise", adjust = 'holm') # Tukey test on marginal means (basically a KR-test, but without random effects)

df<-subset(Basking_a_ser_vi, N_group > 1) # select data on huddling animals from the series vi
df$logdt<-log(df$Delta_Time+1, exp(1)) # log-transform time difference
model<-lm(data = df,
          T_lizard~Sp_N+logdt+Sp_N:logdt) #run linear model
emm.model<-emmeans(model, ~Sp_N) # calculate marginal means
contrast(emm.model, interaction = "pairwise", adjust = 'holm') # Tukey test on marginal means (basically a KR-test, but without random effects)

df<-subset(Basking_a_ser_vii, N_group > 1) # select data on huddling animals from the series vi
df$logdt<-log(df$Delta_Time+1, exp(1)) # log-transform time difference
model<-lm(data = df,
          T_lizard~Sp_N+logdt+Sp_N:logdt) #run linear model
emm.model<-emmeans(model, ~Sp_N) # calculate marginal means
contrast(emm.model, interaction = "pairwise", adjust = 'holm') # Tukey test on marginal means (basically a KR-test, but without random effects)

df<-subset(Basking_a_ser_ix, N_group == 2) # select data on animals in a group of two from the series ix
df$logdt<-log(df$Delta_Time+1, exp(1)) # log-transform time difference
model<-lm(data = df,
          T_lizard~Sp_N+logdt+Sp_N:logdt) # run linear model
emm.model<-emmeans(model, ~Sp_N) # calculate marginal means
contrast(emm.model, interaction = "pairwise", adjust = 'holm') # Tukey test on marginal means (basically a KR-test, but without random effects)

df<-subset(Basking_a_ser_ix, N_group == 3) # select data on animals in a group of three from the series ix
df$logdt<-log(df$Delta_Time+1, exp(1)) # log-transform time difference
model<-lm(data = df,
          T_lizard~Sp_N+logdt+Sp_N:logdt) # run linear model
emm.model<-emmeans(model, ~Sp_N) # calculate marginal means
contrast(emm.model, interaction = "pairwise", adjust = 'holm') # Tukey test on marginal means (basically a KR-test, but without random effects)

df<-subset(Basking_a_ser_xv, N_group ==3) # select data on animals in a group of three from the series xv
df$logdt<-log(df$Delta_Time+1, exp(1)) # log-transform time difference
model<-lm(data = df,
          T_lizard~Sp_N+logdt+Sp_N:logdt) # run linear model
emm.model<-emmeans(model, ~Sp_N) # calculate marginal means
contrast(emm.model, interaction = "pairwise", adjust = 'holm') # Tukey test on marginal means (basically a KR-test, but without random effects)

df<-subset(Basking_a_ser_xv, N_group == 4) # select data on animals in a group of four from the series xv
df$logdt<-log(df$Delta_Time+1, exp(1)) # log-transform time difference

model<-lm(data = df,
          T_lizard~Sp_N+logdt+Sp_N:logdt) # run linear model
emm.model<-emmeans(model, ~Sp_N) # calculate marginal means
contrast(emm.model, interaction = "pairwise", adjust = 'holm') # Tukey test on marginal means (basically a KR-test, but without random effects)
#Fig.3####
a_b_cor_Dif<-cbind(list_of_ci_intercept_Dif$intercept, list_of_ci_Time_Dif$Time,list_of_ci_intercept_Dif$N_group)
a_b_cor_Tb<-cbind(list_of_ci_intercept_Tb$intercept, list_of_ci_Time_Tb$Time,list_of_ci_intercept_Tb$N_group)

a_b_cor_Dif<-as.data.frame(a_b_cor_Dif)
a_b_cor_Tb<-as.data.frame(a_b_cor_Tb)

a_b_cor_Dif$V1<-as.numeric(a_b_cor_Dif$V1)
a_b_cor_Dif$V2<-as.numeric(a_b_cor_Dif$V2)
a_b_cor_Tb$V1<-as.numeric(a_b_cor_Tb$V1)
a_b_cor_Tb$V2<-as.numeric(a_b_cor_Tb$V2)

Fig3A<-ggplot(a_b_cor_Tb, aes(V1, V2))+geom_smooth(method=lm, colour='grey40', fill='grey80')+geom_point(alpha=0.7,size=3,colour='black', shape=21, aes(fill=V3))+
  theme(panel.border = element_rect(size=1, fill=NA), 
        panel.background = element_rect(fill=NA), 
        legend.position = 'top', 
        axis.title =  element_text(size=16), 
        legend.title = element_text(size=18, hjust = .5), 
        legend.text = element_text(size=14),
        plot.title = element_text(size=14, hjust=.95, vjust=-7)
  )+
  xlab('Intercept\n\nNon-adjusted')+
  ylab('Cooling rate')+
  scale_fill_manual(values=c('grey40', "#00BE67","#F8766D", "#00A9FF"))+
  guides(fill=guide_legend(title = 'N of individuals in a group', title.position = 'top'))+
  ggtitle('a')

Fig3B<-ggplot(a_b_cor_Dif, aes(V1, V2))+geom_smooth(method=lm, colour='grey40', fill='grey80')+geom_point(alpha=0.7,size=3,colour='black', shape=21, aes(fill=V3))+
  theme(panel.border = element_rect(size=1, fill=NA), 
        panel.background = element_rect(fill=NA), 
        legend.position = 'top', 
        axis.title =  element_text(size=16), 
        legend.title = element_text(size=18, hjust = .5), 
        legend.text = element_text(size=14),
        plot.title = element_text(size=14, hjust=.95, vjust=-7)
  )+
  xlab('Intercept\n\nAdjusted')+
  ylab('Cooling rate')+
  scale_fill_manual(values=c('grey40',"#F8766D", "#00BE67", "#00A9FF"))+
  guides(fill=guide_legend(title = 'N of individuals in a group', title.position = 'top'))+
  ggtitle('b')

ggarrange(Fig3A, Fig3B, common.legend = T)
#Fig.4####
Fig4a<-ggplot(list_of_series[['i']], aes(log(Delta_Time +1), T_lizard))+
  geom_smooth(aes(fill = Sp_N,
                  colour= Sp_N),
              method='lm',
              alpha = .2)+
  geom_point(aes(fill = Sp_N, shape = as.character(N_group)),
             alpha = 0.7,
             size=2.5)+
  theme_classic()+
  theme(#axis.title.x = element_blank(),
    plot.title = element_text(size=16, face='bold', hjust=.97, vjust = -7),
    axis.title = element_text(size=16),
    legend.position = 'none')+
  xlab(expression(''))+
  ylab('Tb, ºC')+
  scale_fill_manual(values=c("#F8766D",'grey40',"#00A9FF", "#00BE67"))+
  scale_colour_manual(values=c("#F8766D",'grey40',"#00A9FF", "#00BE67"))+
  scale_shape_manual(values=c(21, 23))+
  ggtitle('i')

Fig4b<-ggplot(list_of_series[['ii']], aes(log(Delta_Time +1), T_lizard))+
  geom_smooth(aes(fill = Sp_N,
                  colour= Sp_N),
              method='lm',
              alpha = .2)+
  geom_point(aes(fill = Sp_N, shape = as.character(N_group)),
             alpha = 0.7,
             size=2.5)+
  theme_classic()+
  theme(#axis.title.x = element_blank(),
    plot.title = element_text(size=16, face='bold', hjust=.97, vjust = -7),
    axis.title = element_text(size=16),
    legend.position = 'none'
  )+
  xlab(expression(''))+
  ylab('')+
  scale_fill_manual(values=c("#F8766D",'grey40',"#00A9FF"))+
  scale_colour_manual(values=c("#F8766D",'grey40',"#00A9FF"))+
  scale_shape_manual(values=c(21, 22))+
  ggtitle('ii')

Fig4c<-ggplot(list_of_series[['iv']], aes(log(Delta_Time +1), T_lizard))+
  geom_smooth(aes(fill = Sp_N,
                  colour= Sp_N),
              method='lm',
              alpha = .2)+
  geom_point(aes(fill = Sp_N, shape = as.character(N_group)),
             alpha = 0.7,
             size=2.5)+
  theme_classic()+
  theme(#axis.title.x = element_blank(),
    plot.title = element_text(size=16, face='bold', hjust=.97, vjust = -7),
    axis.title = element_text(size=16),
    legend.position = 'none'
  )+
  xlab(expression(''))+
  ylab('Tb, ºC')+
  scale_fill_manual(values=c('grey40', "#F8766D","#00A9FF"))+
  scale_colour_manual(values=c('grey40', "#F8766D", "#00A9FF"))+
  scale_shape_manual(values=c(21, 22))+
  ggtitle('iii')

Fig4d<-ggplot(list_of_series[['v']], aes(log(Delta_Time +1), T_lizard))+
  geom_smooth(aes(fill = Sp_N,
                  colour= Sp_N),
              method='lm',
              alpha = .2)+
  geom_point(aes(fill = Sp_N, shape = as.character(N_group)),
             alpha = 0.7,
             size=2.5)+
  theme_classic()+
  theme(#axis.title.x = element_blank(),
    plot.title = element_text(size=16, face='bold', hjust=.97, vjust = -7),
    axis.title = element_text(size=16),
    legend.position = 'none'
  )+
  xlab(expression(''))+
  ylab('')+
  scale_fill_manual(values=c('grey40',"#F8766D", "#00A9FF"))+
  scale_colour_manual(values=c('grey40',"#F8766D","#00A9FF"))+
  scale_shape_manual(values=c(21, 22))+
  ggtitle('iv')

Fig4e<-ggplot(list_of_series[['vi']], aes(log(Delta_Time +1), T_lizard))+
  geom_smooth(aes(fill = Sp_N,
                  colour= Sp_N),
              method='lm',
              alpha = .2)+
  geom_point(aes(fill = Sp_N, shape = as.character(N_group)),
             alpha = 0.7,
             size=2.5)+
  theme_classic()+
  theme(#axis.title.x = element_blank(),
    plot.title = element_text(size=16, face='bold', hjust=.97, vjust = -7),
    axis.title = element_text(size=16),
    legend.position = 'none'
  )+
  xlab(expression(''))+
  ylab('Tb, ºC')+
  scale_fill_manual(values=c('grey40',"#F8766D","#00A9FF"))+
  scale_colour_manual(values=c('grey40',"#F8766D","#00A9FF"))+
  scale_shape_manual(values=c(21, 22))+
  ggtitle('v')

Fig4f<-ggplot(list_of_series[['vii']], aes(log(Delta_Time +1), T_lizard))+
  geom_smooth(aes(fill = Sp_N,
                  colour= Sp_N),
              method='lm',
              alpha = .2)+
  geom_point(aes(fill = Sp_N, shape = as.character(N_group)),
             alpha = 0.7,
             size=2.5)+
  theme_classic()+
  theme(#axis.title.x = element_blank(),
    plot.title = element_text(size=16, face='bold', hjust=.97, vjust = -7),
    axis.title = element_text(size=16),
    legend.position = 'none'
  )+
  xlab(expression(''))+
  ylab('')+
  scale_fill_manual(values=c("#F8766D","#00A9FF", 'grey40',"#00BE67"))+
  scale_colour_manual(values=c("#F8766D","#00A9FF",'grey40',"#00BE67"))+
  scale_shape_manual(values=c(21, 23))+
  ggtitle('vi')

Fig4g<-ggplot(list_of_series[['ix']], aes(log(Delta_Time +1), T_lizard))+
  geom_smooth(aes(fill = Sp_N,
                  colour= Sp_N),
              method='lm',
              alpha = .2)+
  geom_point(aes(fill = Sp_N, shape = as.character(N_group)),
             alpha = 0.7,
             size=2.5)+
  theme_classic()+
  theme(#axis.title.x = element_blank(),
    plot.title = element_text(size=16, face='bold', hjust=.97, vjust = -7),
    axis.title = element_text(size=16),
    legend.position = 'none'
  )+
  xlab(expression('ln('~Delta~'t+1)'))+
  ylab('Tb, ºC')+
  scale_fill_manual(values=c("#F8766D","#00A9FF", '#BB22BB', '#EEEE22',"#00BE67"))+
  scale_colour_manual(values=c("#F8766D","#00A9FF", '#BB22BB', '#EEEE22',"#00BE67"))+
  scale_shape_manual(values=c(22, 23))+
  ggtitle('vii')

Fig4h<-ggplot(list_of_series[['xv']], aes(log(Delta_Time +1), T_lizard))+
  geom_smooth(aes(fill = Sp_N,
                  colour= Sp_N),
              method='lm',
              alpha = .2)+
  geom_point(aes(fill = Sp_N, shape = as.character(N_group)),
             alpha = 0.7,
             size=2.5)+
  theme_classic()+
  theme(#axis.title.x = element_blank(),
    plot.title = element_text(size=16, face='bold', hjust=.97, vjust = -7),
    axis.title = element_text(size=16),
    legend.position = 'none'
  )+
  xlab(expression('ln('~Delta~'t+1)'))+
  ylab('')+
  scale_fill_manual(values=c("#F8766D","#00A9FF", 'gray40', '#EEEE22',"#00BE67"))+
  scale_colour_manual(values=c("#F8766D","#00A9FF", 'gray40', '#EEEE22',"#00BE67"))+
  scale_shape_manual(values=c(21, 23, 24))+
  ggtitle('viii')

leg.plot<-ggplot(full_join(list_of_series[['ii']], list_of_series[['xv']]),
                 aes(log(Delta_Time+1, exp(1)), Dif, shape=as.character(N_group)))+
  geom_point(size=5, fill='gray70')+
  scale_shape_manual(values = c(21, 22, 23, 24), name='Ng')+
  theme_classic()+
  theme(legend.title.position = 'top', legend.title = element_text(size=14, face = 'bold', hjust = .5), legend.text = element_text(size=12))+
  guides(shape=guide_legend(direction='horizontal'))

leg<-get_legend(leg.plot)

Fig4<-plot_grid(Fig4a, Fig4b, Fig4c, Fig4d, Fig4e, Fig4f, Fig4g, Fig4h, NULL, nrow=5, ncol=2, align = 'hv', rel_widths = c(rep(1, 8), 2))

Fig4+draw_plot(leg, x=1/4, y=0, width=1/2, height=0.4)
#Table SI####
SI1 <- cbind(list_of_ci_Time_Tb, list_of_ci_Time_Dif)
SI1 <- SI1[, -c(3, 5, 6, 7, 10)]
colnames(SI1) <- c('nonadjusted_lower_bond',
                   'nonadjusted_upper_bond',
                   'nonadjusted_rate',
                   'adjusted_lower_bond',
                   'adjusted_upper_bond',
                   'adjusted_rate',
                   'Series',
                   'ID',
                   'N_group')
rownames(SI1)<-NULL
