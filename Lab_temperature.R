#Packages####
library(gsheet)
library(ggplot2)
library(tidyverse)
library(ggpubr)
#Data####
df<-read.csv('logger_data.csv')
df_slice<-df %>%
  group_by(date) %>%
  top_n(1, wt=temp)
df_slice$date<-as.Date(df_slice$date, format = '%d.%m.%Y')
df_slice<-subset(df_slice, date > as.Date('01.05.2021', format = '%d.%m.%Y') &
                   date < as.Date('01.05.2022', format = '%d.%m.%Y'))
#Plot####
Sys.setlocale("LC_TIME", "English")
plot_temp<-ggplot(df_slice, aes(date, temp))+
  geom_point()+
  geom_path()+
  geom_smooth(span=0.20, se=F)+
  theme_classic()+
    xlab('')+
  ylab('Air temperature, °C')+
  scale_x_date(date_breaks = '2 month', date_labels = '%b %Y')+
  geom_rect(aes(xmin=as.Date('2021-11-16'), xmax=as.Date('2022-03-07'), ymin=0, ymax=40), fill='gray60', alpha=.008)+
  #geom_rect(aes(xmin=as.Date('2022-11-05'), xmax=as.Date('2023-02-23'), ymin=0, ymax=40), fill='gray60', alpha=.008)+
  theme(axis.title = element_text(size=20), plot.title = element_text(size=20))+
  ggtitle('a')

plot_hum<-ggplot(df_slice, aes(date, hum))+
  geom_point()+
  geom_path()+
  geom_smooth(span=0.2, se=F)+
  theme_classic()+
  xlab('Date')+
  ylab('Humidity, %')+
  scale_x_date(date_breaks = '2 month', date_labels = '%b %Y')+
  geom_rect(aes(xmin=as.Date('2021-11-16'), xmax=as.Date('2022-03-07'), ymin=0, ymax=95), fill='gray60', alpha=.008)+
  #geom_rect(aes(xmin=as.Date('2022-11-05'), xmax=as.Date('2023-02-23'), ymin=0, ymax=100), fill='gray60', alpha=.008)+
  theme(axis.title = element_text(size=20), plot.title = element_text(size=20))+
  ggtitle('b')


ggarrange(plot_temp, plot_hum, nrow=2)