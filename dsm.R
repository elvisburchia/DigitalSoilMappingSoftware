# DSM - Digital Soil Mapping
# UPDATE 04.02.2021
# Elvis Burchia & Giulio Genova
# Master thesis of Elvis Burchia 2021
# Program for mapping chemical soil properties
# Master of ecology and biodiversity (MSc)
# University of Innsbruck


# Libraries and functions ---------------------------------------------------------------
library(caret)
library(rgdal)
library(raster)
library(sp)
library(gstat)
library(randomForest)
library(ranger)
library(GSIF)
library(hydroGOF)
library(tidyr)
library(broom)
library(xlsx)
library(datasets)
library(ggplot2)
library(gridExtra)
library(e1071)
library(png)
library(corrplot)
library(data.table)
library(tidyverse)
library(ggpubr)
library(grid) 
library(gtable)
library(readr)
library(pdftools)
library(EBImage)
quad <- function(x){
  y <- x^2;
  return(y)}

inverse <- function(x){
  y <- 1/x;
  return(y)}

all_combn = function(vec){
  unlist(sapply(vec, function(y) combn(vec, y, paste, collapse = " and ")))
}
finderrors <- function(x){
  iscutoff<-2
  fehlererrors<-x
  fehlererrorsname<-names(fehlererrors)
  sequenza<-grepl("No convergence", fehlererrorsname, fixed = TRUE)
  for (sd in 1:length(sequenza)) {
    if (sequenza[sd]==TRUE) {
      iscutoff<-0
    } else {
      iscutoff<-1
    }
  }
  return(iscutoff)}

res <- 150
resolution <- paste0("",res,"")
predictors_combination <- "personalized"  # or "automatic"

# Define directory name and folders ---------------------------------------
script_dir<-dirname(rstudioapi::getActiveDocumentContext()$path)
script_path<-rstudioapi::getActiveDocumentContext()$path
predictors_dir_0<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/_predictors_independent_variables")
predictors_dir_1<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/_predictors_independent_variables/predictors_predictarea_1")
predictors_dir_2<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/_predictors_independent_variables/predictors_predictarea_2")
predictors_name_1<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/_predictors_independent_variables/names/names_predictarea_1")
predictors_name_2<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/_predictors_independent_variables/names/names_predictarea_2")
datapoints_dir<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/_masterfile_data_points")
predict_area_dir_0<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/_predict_area")
predict_area_dir_1<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/_predict_area/1_predict_area")
predict_area_dir_2<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/_predict_area/2_predict_area")
# Remove output folder with old outputs
setwd(script_dir)
unlink("./output", recursive=TRUE)
unlink("./output_overview", recursive=TRUE)
# Create output folder
dir.create("output")
out_dir<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/output")
dir.create("output_overview")

# Predict areas (input and output folders) --------------------------------
# Create output folders for different predict areas
setwd(predict_area_dir_0)
nrpredictareas<-length(dir()[file.info(dir())$isdir])+0
if (nrpredictareas==1){
  setwd(predict_area_dir_1)
  predictareas_name<-dir()[file.info(dir())$isdir]
} else if (nrpredictareas==2){
  predictareas_name<-c("predictarea1","predictarea2")
  setwd(predict_area_dir_1)
  predictareas_name[1]<-dir()[file.info(dir())$isdir]
  setwd(predict_area_dir_2)
  predictareas_name[2]<-dir()[file.info(dir())$isdir]
} else {}
setwd(out_dir)
for (i in 1:nrpredictareas)
{
  dir.create(predictareas_name[i])
  ff<-paste0("/output/",predictareas_name[i])
  assign(paste0("output_predictarea_dir",i),paste0(dirname(rstudioapi::getActiveDocumentContext()$path),ff))
}

for (i in 1:nrpredictareas) {
  if (i==1) {
    setwd(paste0(script_dir,"/output_overview"))
    dir.create("Overview_Predict_area_1")
  } else if (i==2) {
    setwd(paste0(script_dir,"/output_overview"))
    dir.create("Overview_Predict_area_2")
  } else {}
}
setwd(script_dir)

# Input of rasters of the predicted areas
if (nrpredictareas==1){
  assign(paste0("predictarea_1","_",predictareas_name[1]),raster(file.path(predict_area_dir_1,predictareas_name[1])))
  assign("predictarea_1",raster(file.path(predict_area_dir_1,predictareas_name[1])))
} else if (nrpredictareas==2){
  assign(paste0("predictarea_1","_",predictareas_name[1]),raster(file.path(predict_area_dir_1,predictareas_name[1])))
  assign("predictarea_1",raster(file.path(predict_area_dir_1,predictareas_name[1])))
  assign(paste0("predictarea_2","_",predictareas_name[2]),raster(file.path(predict_area_dir_2,predictareas_name[2])))
  assign("predictarea_2",raster(file.path(predict_area_dir_2,predictareas_name[2])))
} else {}

# Create output folder for different variables ----------------------------
# Variables from Masterfile (table with all data points)
obs<- read.csv(file = file.path(datapoints_dir,list.files(path = datapoints_dir)),
               header = TRUE,sep = ",",dec = ".",na.strings = c("-9999"))  
variables1<-names(obs)
lenvariables1<-length(variables1)+0
variables2<-variables1[-c((lenvariables1-1):lenvariables1)]
variables<-variables2[-c(1:8)]
lenvariables<-length(variables)+0
rm(variables1)
rm(variables2)
rm(lenvariables1)
# Create output folders for every variable
transformations<-c("raw","log","sqrt","quad","inverse")
lentransformations<-length(transformations)+0
for (i in 1:nrpredictareas){
  for (j in 1:lenvariables){
    ff<-paste0("/output/",predictareas_name[i])
    fz<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),ff)
    setwd(fz)
    dir.create(variables[j])
    fa<-paste0("/output/",predictareas_name[i],"/",variables[j])
    assign(paste0("output_",predictareas_name[i],"_",variables[j],"_dir"),paste0(dirname(rstudioapi::getActiveDocumentContext()$path),fa))
    if (i==1){
      setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_1"))
      dir.create(variables[j])
      setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_1/",variables[j]))
      dir.create("Overview_Tables")
      dir.create("Runs_Output")
      dir.create("Overview_Prediction_Ordinary_Kriging")
      dir.create("Overview_Exploratory_Statistic")
      dir.create("Final_Models")
      dir.create("Find_Best_Model")
      dir.create("PredictorsInfluence")
    } else if (i==2) {
      setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_2"))
      dir.create(variables[j])
      setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_2/",variables[j]))
      dir.create("Overview_Tables")
      dir.create("Runs_Output")
      dir.create("Overview_Prediction_Ordinary_Kriging")
      dir.create("Overview_Exploratory_Statistic")
      dir.create("Final_Models")
      dir.create("Find_Best_Model")
      dir.create("PredictorsInfluence")
    } else {}
    for (k in 1:lentransformations){
      fg<-paste0("/output/",predictareas_name[i],"/",variables[j])
      fh<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),fg)
      setwd(fh)
      dir.create(transformations[k])
      fk<-paste0("/output/",predictareas_name[i],"/",variables[j],"/",transformations[k])
      assign(paste0("output_",predictareas_name[i],"_",variables[j],transformations[k],"_dir"),paste0(dirname(rstudioapi::getActiveDocumentContext()$path),fk))
      
      # Clean data and output table to choose the right transformation
      var<-variables[j]
      obs<- read.csv(file = file.path(datapoints_dir,list.files(path = datapoints_dir)),
                     header = TRUE,sep = ",",dec = ".",na.strings = c("-9999"))  
      obs<- obs[!duplicated(obs[,c("x","y")]),]
      obs<-obs[which(!is.na(obs[,var])),]
      shapraw<-shapiro.test(obs[,var])
      shaplog<-shapiro.test(log(obs[,var]))
      shapsqrt<-shapiro.test(sqrt(obs[,var]))
      shapquad<-shapiro.test(quad(obs[,var]))
      shapinverse<-shapiro.test(inverse(obs[,var]))
      a <- c(transformations)
      b <- c(skewness(obs[,var]),skewness(log(obs[,var])), skewness(sqrt(obs[,var])), skewness(quad(obs[,var])), skewness(inverse(obs[,var])))
      c <- c(kurtosis(obs[,var]),kurtosis(log(obs[,var])), kurtosis(sqrt(obs[,var])), kurtosis(quad(obs[,var])), kurtosis(inverse(obs[,var])))
      d <- c(shapraw$p.value,shaplog$p.value,shapsqrt$p.value,shapquad$p.value,shapinverse$p.value)
      e <- c(shapraw$statistic,shaplog$statistic,shapsqrt$statistic,shapquad$statistic,shapinverse$statistic)
      f <- c("YES","YES","YES","YES","YES")
      df <- data.frame(a,b,c,d,e,f)
      names(df)<-c("transformation","skewness","kurtosis","shapiro_p.value","shapiro_statistic","Normal_distributed_data")
      if (df[1,5]<0.90){
        unlink("./raw", recursive=TRUE)
        df[1,6]<-"NO"
      }
      if (df[2,5]<0.90){
        unlink("./log", recursive=TRUE)
        df[2,6]<-"NO"
      }
      if (df[3,5]<0.90){
        unlink("./sqrt", recursive=TRUE)
        df[3,6]<-"NO"
      }
      if (df[4,5]<0.90){
        unlink("./quad", recursive=TRUE)
        df[4,6]<-"NO"
      }
      if (df[5,5]<0.90){
        unlink("./inverse", recursive=TRUE)
        df[5,6]<-"NO"
      }
      if (length(dir(all.files=TRUE))==0){
        dir.create("NoTransformationIsOk_NotNormallyDistributedData")
      }
      write.csv(df,file.path(fh,paste0("Choose_Transformation_",var,"_",predictareas_name[i],".csv")))
      write.xlsx(df,file.path(fh,paste0("Choose_Transformation_",var,"_",predictareas_name[i],".xlsx")))
      d <- head(df) 
      table <- tableGrob(d)
      title <- textGrob(paste0("Choose_Transformation_",var,"_",predictareas_name[i]),gp=gpar(fontsize=15)) 
      footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
      padding <- unit(2,"line") 
      table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
      table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
      table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
      png(paste0("Choose_Transformation_",var,"_",predictareas_name[i],".png"), height = 200, width = 630)
      grid.newpage() 
      grid.draw(table) 
      dev.off()
      # Save summary of every variable
      sumraw <- summary(obs[,var])
      yraw<-data.frame(sum=matrix(sumraw),row.names=names(sumraw))
      sumlog <- summary(log(obs[,var]))
      ylog<-data.frame(sum=matrix(sumlog),row.names=names(sumlog))
      sumsqrt <- summary(sqrt(obs[,var]))
      ysqrt<-data.frame(sum=matrix(sumsqrt),row.names=names(sumsqrt))
      sumquad <- summary(quad(obs[,var]))
      yquad<-data.frame(sum=matrix(sumquad),row.names=names(sumquad))
      suminverse <- summary(inverse(obs[,var]))
      yinverse<-data.frame(sum=matrix(suminverse),row.names=names(suminverse))
      a<-c(transformations)
      b<-c(yraw[1,1],ylog[1,1],ysqrt[1,1],yquad[1,1],yinverse[1,1])
      c<-c(yraw[2,1],ylog[2,1],ysqrt[2,1],yquad[2,1],yinverse[2,1])
      d<-c(yraw[3,1],ylog[3,1],ysqrt[3,1],yquad[3,1],yinverse[3,1])
      e<-c(yraw[4,1],ylog[4,1],ysqrt[4,1],yquad[4,1],yinverse[4,1])
      f<-c(yraw[5,1],ylog[5,1],ysqrt[5,1],yquad[5,1],yinverse[5,1])
      g<-c(yraw[6,1],ylog[6,1],ysqrt[6,1],yquad[6,1],yinverse[6,1])
      lendf<-length(df)+0
      h<-c(df[,6])
      dff<-data.frame(a,b,c,d,e,f,g,h)
      names(dff)<-c("transformation",row.names(yraw),"Normal_distributed_data")
      write.csv(dff,file.path(fh,paste0("Descriptive_Statistic_",var,"_",predictareas_name[i],".csv"))) 
      write.xlsx(dff,file.path(fh,paste0("Descriptive_Statistic_",var,"_",predictareas_name[i],".xlsx")))
      d <- head(dff) 
      table <- tableGrob(d)
      title <- textGrob(paste0("Descriptive_Statistic_",var,"_",predictareas_name[i]),gp=gpar(fontsize=15)) 
      footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
      padding <- unit(2,"line") 
      table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
      table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
      table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
      png(paste0("Descriptive_Statistic_",var,"_",predictareas_name[i],".png"), height = 200, width = 790)
      grid.newpage() 
      grid.draw(table) 
      dev.off()
    }
  }
}
setwd(script_dir)

# Output Summary Graphics -------------------------------------------------

for (i in 1:nrpredictareas) {
  print(predictareas_name[i])
  for (j in 1:lenvariables) {
    print(variables[j])
    ff<-paste0("/output/",predictareas_name[i],"/",variables[j])
    fz<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),ff)
    setwd(fz)
    ai<-dir()[file.info(dir())$isdir]
    lange<-length(ai)
    # Read Masterfile - Soil Data Points
    obs<- read.csv(file = file.path(datapoints_dir,list.files(path = datapoints_dir)),header = TRUE,sep = ",",dec = ".",na.strings = c("-9999"))  
    obs<- obs[!duplicated(obs[,c("x","y")]),]
    obs<-obs[which(!is.na(obs[,variables[j]])),]
    # Variable description for graphics
    if (variables[j]=="SOM"){
      varr = variables[j]
      einheit = "%"
      varbez = paste0(varr," ","[",einheit,"]")
    } else if (variables[j]=="pH"){
      varr = variables[j]
      einheit = " "
      varbez = varr
    } else if (variables[j]=="P_mg100g"){
      varr = "P"
      einheit = "mg/100g"
      varbez = paste0(varr," ","[",einheit,"]")
    } else if (variables[j]=="K_mg100g"){
      varr = "K"
      einheit = "mg/100g"
      varbez = paste0(varr," ","[",einheit,"]")
    } else if (variables[j]=="C"){
      varr = variables[j]
      einheit = "%"
      varbez = paste0(varr," ","[",einheit,"]")
    } else if (variables[j]=="N"){
      varr = variables[j]
      einheit = "%"
      varbez = paste0(varr," ","[",einheit,"]")
    } else if (variables[j]=="CN") {
      varr = variables[j]
      einheit = " "
      varbez = variables[j]
    } else {
      varbez = variables[j]
    }
    
    for (k in 1:lange) {
      print(ai[k])
      fg<-paste0("/output/",predictareas_name[i],"/",variables[j],"/",ai[k])
      fh<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),fg)
      setwd(fh)
      dir.create("Summary_Graphics")
      fgg<-paste0("/output/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","Summary_Graphics")
      fhh<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),fgg)
      setwd(fhh)
      
      if (ai[k]=="raw") {
        # Histogram___________________________
        #png("Histogram.png",res=600,width =4000, height=4000)
        #ggplot(data=obs, aes(obs[,variables[j]])) + geom_histogram() + 
        #scale_x_continuous(breaks=seq(0,max(obs[,variables[j]]),5), lim=c(0,max(obs[,variables[j]]))) + 
        #labs(x = paste0(varbez))+ggtitle(paste0("Histogram - ",variables[j],"  - (","raw - ",predictareas_name[i],")")) +
        #theme_minimal() +
        #theme(plot.title = element_text(hjust = 0.5))
        #dev.off()
        #png("Histogram.png",res=600,width =4000, height=4000)
        #histogram(obs[,variables[j]],xlab=paste0(varbez,"  - (","raw - ",predictareas_name[i],")"), col = "gray")
        #dev.off()
        png("Histogram_Normal_Distribution_Curve.png",res=600,width =4000, height=4000)
        m <- mean(obs[,variables[j]])
        std <- sqrt(var(obs[,variables[j]]))
        tmp <- density(obs[,variables[j]])
        hist(obs[,variables[j]],xlab=paste0(varbez,"  - (","raw - ",predictareas_name[i],")"), prob=TRUE, main=NULL, col="gray",xlim=c(0,max(tmp$x))) #,ylim=c(0,max(tmp$y)))
        curve(dnorm(x, mean=m, sd=std), col="red", lwd=2, add=TRUE, yaxt="n")
        dev.off()
        # QQ-Plot______________________________
        png("QQPlot.png",res=600,width =4000, height=4000)
        qqnorm(obs[,variables[j]],main=paste0("Normal Q-Q Plot - ",varbez," - (","raw - ",predictareas_name[i],")"))
        qqline(obs[,variables[j]])
        dev.off()
        # Boxplot and density__________________
        axis_text_size<-10
        axis_title_size<-10 
        legend_text_size<-10
        table_text_size<-10
        
        th_box<-theme(panel.background = element_rect(fill = NA),
                      panel.grid.major = element_line(linetype = "dashed",colour = "grey"),
                      legend.position = "none" ,
                      axis.title.y=element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks.y=element_blank(),
                      axis.title.x=element_blank(),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())
        
        th_dens<-theme(panel.background = element_rect(fill = NA),
                       panel.grid.major = element_line(linetype = "dashed",colour = "grey"),
                       axis.text.x=element_text(face = "plain",size = axis_text_size,colour = "black"),
                       axis.text.y=element_text(face = "plain",size = axis_text_size,colour = "black"),
                       #axis.title.y = element_text(size = axis_title_size,margin = margin(r=20,l=20)),
                       axis.title.x = element_text(size = axis_title_size),#,margin = margin(r=20,l=20)
                       legend.position = c(0.7, 0.7),
                       legend.title = element_blank(),
                       legend.text = element_text(size = legend_text_size,margin = 10),
                       #legend.key = element_rect(size = 2),
                       legend.key.size = unit(3, 'lines'),
                       axis.line.x = element_line(colour = "black",size=1.2))
        
        h <- ggplot(data=obs, aes_string(x=variables[j])) + geom_boxplot()+ th_box
        b <- ggplot(data=obs, aes_string(variables[j])) + geom_density(fill="blue",alpha=0.5) + th_dens + labs(x = paste0(varbez,"  - (","raw - ",predictareas_name[i],")"))
        plot_total<- ggarrange(h, b, heights = c(0.5,1), ncol = 1, nrow =2, align = "v")
        plot_total
        ggsave(filename = "Boxplot_Density.png",plot = plot_total,device ='png',width = 5, height = 4, dpi = 600, units = "in")
        
      } 
      
      else if (ai[k]=="log") {
        # Histogram___________________________
        #png("Histogram.png",res=600,width =4000, height=4000)
        #ggplot(data=obs, aes(log(obs[,variables[j]]))) + geom_histogram() + 
        #scale_x_continuous(breaks=seq(0,max(log(obs[,variables[j]])),5), lim=c(0,max(log(obs[,variables[j]])))) + 
        #labs(x = paste0(varbez))+ggtitle(paste0("Histogram - ",variables[j],"  - (","log - ",predictareas_name[i],")")) +
        #theme_minimal() +
        #theme(plot.title = element_text(hjust = 0.5))
        #dev.off()
        #png("Histogram.png",res=600,width =4000, height=4000)
        #histogram(log(obs[,variables[j]]),xlab=paste0(varbez,"  - (","log - ",predictareas_name[i],")"), col = "gray")
        #dev.off()
        png("Histogram_Normal_Distribution_Curve.png",res=600,width =4000, height=4000)
        m <- mean(log(obs[,variables[j]]))
        std <- sqrt(var(log(obs[,variables[j]])))
        tmp <- density(log(obs[,variables[j]]))
        hist(log(obs[,variables[j]]),xlab=paste0(varbez,"  - (","log - ",predictareas_name[i],")"), prob=TRUE, main=NULL, col="gray",xlim=c(0,max(tmp$x))) #,ylim=c(0,max(tmp$y)))
        curve(dnorm(x, mean=m, sd=std), col="red", lwd=2, add=TRUE, yaxt="n")
        dev.off()
        # QQ-Plot______________________________
        png("QQPlot.png",res=600,width =4000, height=4000)
        qqnorm(log(obs[,variables[j]]),main=paste0("Normal Q-Q Plot - ",varbez," - (","log - ",predictareas_name[i],")"))
        qqline(log(obs[,variables[j]]))
        dev.off()
        # Boxplot and density__________________
        axis_text_size<-10
        axis_title_size<-10 
        legend_text_size<-10
        table_text_size<-10
        
        th_box<-theme(panel.background = element_rect(fill = NA),
                      panel.grid.major = element_line(linetype = "dashed",colour = "grey"),
                      legend.position = "none" ,
                      axis.title.y=element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks.y=element_blank(),
                      axis.title.x=element_blank(),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())
        
        th_dens<-theme(panel.background = element_rect(fill = NA),
                       panel.grid.major = element_line(linetype = "dashed",colour = "grey"),
                       axis.text.x=element_text(face = "plain",size = axis_text_size,colour = "black"),
                       axis.text.y=element_text(face = "plain",size = axis_text_size,colour = "black"),
                       #axis.title.y = element_text(size = axis_title_size,margin = margin(r=20,l=20)),
                       axis.title.x = element_text(size = axis_title_size),#,margin = margin(r=20,l=20)
                       legend.position = c(0.7, 0.7),
                       legend.title = element_blank(),
                       legend.text = element_text(size = legend_text_size,margin = 10),
                       #legend.key = element_rect(size = 2),
                       legend.key.size = unit(3, 'lines'),
                       axis.line.x = element_line(colour = "black",size=1.2))
        
        o <- ggplot(data=obs, aes_string(paste0("x=log(",variables[j],")"))) + geom_boxplot()+ th_box
        a <- ggplot(data=obs, aes_string(paste0("log(",variables[j],")"))) + geom_density(fill="blue",alpha=0.5) + th_dens + labs(x = paste0(varbez,"  - (","log - ",predictareas_name[i],")"))
        plot_total_log<- ggarrange(o, a, heights = c(0.5,1), ncol = 1, nrow =2, align = "v")
        plot_total_log
        ggsave(filename="Boxplot_Density.png",plot=plot_total_log,device ='png',width = 5, height = 4, dpi = 600, units = "in")
        
      } 
      
      else if (ai[k]=="sqrt") {
        # Histogram___________________________
        #png("Histogram.png",res=600,width =4000, height=4000)
        #ggplot(data=obs, aes(sqrt(obs[,variables[j]]))) + geom_histogram() + 
        #scale_x_continuous(breaks=seq(0,max(sqrt(obs[,variables[j]])),5), lim=c(0,max(sqrt(obs[,variables[j]])))) + 
        #labs(x = paste0(varbez))+ggtitle(paste0("Histogram - ",variables[j],"  - (","sqrt - ",predictareas_name[i],")")) +
        #theme_minimal() +
        #theme(plot.title = element_text(hjust = 0.5))
        #dev.off()
        #png("Histogram.png",res=600,width =4000, height=4000)
        #histogram(sqrt(obs[,variables[j]]),xlab=paste0(varbez,"  - (","sqrt - ",predictareas_name[i],")"), col = "gray")
        #dev.off()
        png("Histogram_Normal_Distribution_Curve.png",res=600,width =4000, height=4000)
        m <- mean(sqrt(obs[,variables[j]]))
        std <- sqrt(var(sqrt(obs[,variables[j]])))
        tmp <- density(sqrt(obs[,variables[j]]))
        hist(sqrt(obs[,variables[j]]),xlab=paste0(varbez,"  - (","sqrt - ",predictareas_name[i],")"), prob=TRUE, main=NULL, col="gray",xlim=c(0,max(tmp$x))) #,ylim=c(0,max(tmp$y)))
        curve(dnorm(x, mean=m, sd=std), col="red", lwd=2, add=TRUE, yaxt="n")
        dev.off()
        # QQ-Plot______________________________
        png("QQPlot.png",res=600,width =4000, height=4000)
        qqnorm(sqrt(obs[,variables[j]]),main=paste0("Normal Q-Q Plot - ",varbez," - (","sqrt - ",predictareas_name[i],")"))
        qqline(sqrt(obs[,variables[j]]))
        dev.off()
        # Boxplot and density__________________
        axis_text_size<-10
        axis_title_size<-10 
        legend_text_size<-10
        table_text_size<-10
        
        th_box<-theme(panel.background = element_rect(fill = NA),
                      panel.grid.major = element_line(linetype = "dashed",colour = "grey"),
                      legend.position = "none" ,
                      axis.title.y=element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks.y=element_blank(),
                      axis.title.x=element_blank(),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())
        
        th_dens<-theme(panel.background = element_rect(fill = NA),
                       panel.grid.major = element_line(linetype = "dashed",colour = "grey"),
                       axis.text.x=element_text(face = "plain",size = axis_text_size,colour = "black"),
                       axis.text.y=element_text(face = "plain",size = axis_text_size,colour = "black"),
                       #axis.title.y = element_text(size = axis_title_size,margin = margin(r=20,l=20)),
                       axis.title.x = element_text(size = axis_title_size),#,margin = margin(r=20,l=20)
                       legend.position = c(0.7, 0.7),
                       legend.title = element_blank(),
                       legend.text = element_text(size = legend_text_size,margin = 10),
                       #legend.key = element_rect(size = 2),
                       legend.key.size = unit(3, 'lines'),
                       axis.line.x = element_line(colour = "black",size=1.2))
        
        o <- ggplot(data=obs, aes_string(paste0("x=sqrt(",variables[j],")"))) + geom_boxplot()+ th_box
        a <- ggplot(data=obs, aes_string(paste0("sqrt(",variables[j],")"))) + geom_density(fill="blue",alpha=0.5) + th_dens + labs(x = paste0(varbez,"  - (","sqrt - ",predictareas_name[i],")"))
        plot_total_sqrt<- ggarrange(o, a, heights = c(0.5,1), ncol = 1, nrow =2, align = "v")
        plot_total_sqrt
        ggsave(filename="Boxplot_Density.png",plot=plot_total_sqrt,device ='png',width = 5, height = 4, dpi = 600, units = "in")
        
      } 
      
      else if (ai[k]=="quad") {
        # Histogram___________________________
        #png("Histogram.png",res=600,width =4000, height=4000)
        #ggplot(data=obs, aes(quad(obs[,variables[j]]))) + geom_histogram() + 
        #scale_x_continuous(breaks=seq(0,max(quad(obs[,variables[j]])),5), lim=c(0,max(quad(obs[,variables[j]])))) + 
        #labs(x = paste0(varbez))+ggtitle(paste0("Histogram - ",variables[j],"  - (","quad - ",predictareas_name[i],")")) +
        #theme_minimal() +
        #theme(plot.title = element_text(hjust = 0.5))
        #dev.off()
        #png("Histogram.png",res=600,width =4000, height=4000)
        #histogram(quad(obs[,variables[j]]),xlab=paste0(varbez,"  - (","^2 - ",predictareas_name[i],")"), col = "gray")
        #dev.off()
        png("Histogram_Normal_Distribution_Curve.png",res=600,width =4000, height=4000)
        m <- mean(quad(obs[,variables[j]]))
        std <- sqrt(var(quad(obs[,variables[j]])))
        tmp <- density(quad(obs[,variables[j]]))
        hist(quad(obs[,variables[j]]),xlab=paste0(varbez,"  - (","^2 - ",predictareas_name[i],")"), prob=TRUE, main=NULL, col="gray",xlim=c(0,max(tmp$x))) #,ylim=c(0,max(tmp$y)))
        curve(dnorm(x, mean=m, sd=std), col="red", lwd=2, add=TRUE, yaxt="n")
        dev.off()
        # QQ-Plot______________________________
        png("QQPlot.png",res=600,width =4000, height=4000)
        qqnorm(quad(obs[,variables[j]]),main=paste0("Normal Q-Q Plot - ",varbez," - (","^2 - ",predictareas_name[i],")"))
        qqline(quad(obs[,variables[j]]))
        dev.off()
        # Boxplot and density__________________
        axis_text_size<-10
        axis_title_size<-10 
        legend_text_size<-10
        table_text_size<-10
        
        th_box<-theme(panel.background = element_rect(fill = NA),
                      panel.grid.major = element_line(linetype = "dashed",colour = "grey"),
                      legend.position = "none" ,
                      axis.title.y=element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks.y=element_blank(),
                      axis.title.x=element_blank(),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())
        
        th_dens<-theme(panel.background = element_rect(fill = NA),
                       panel.grid.major = element_line(linetype = "dashed",colour = "grey"),
                       axis.text.x=element_text(face = "plain",size = axis_text_size,colour = "black"),
                       axis.text.y=element_text(face = "plain",size = axis_text_size,colour = "black"),
                       #axis.title.y = element_text(size = axis_title_size,margin = margin(r=20,l=20)),
                       axis.title.x = element_text(size = axis_title_size),#,margin = margin(r=20,l=20)
                       legend.position = c(0.7, 0.7),
                       legend.title = element_blank(),
                       legend.text = element_text(size = legend_text_size,margin = 10),
                       #legend.key = element_rect(size = 2),
                       legend.key.size = unit(3, 'lines'),
                       axis.line.x = element_line(colour = "black",size=1.2))
        
        o <- ggplot(data=obs, aes_string(paste0("x=quad(",variables[j],")"))) + geom_boxplot()+ th_box
        a <- ggplot(data=obs, aes_string(paste0("quad(",variables[j],")"))) + geom_density(fill="blue",alpha=0.5) + th_dens + labs(x = paste0(varbez,"  - (","^2 - ",predictareas_name[i],")"))
        plot_total_quad<- ggarrange(o, a, heights = c(0.5,1), ncol = 1, nrow =2, align = "v")
        plot_total_quad
        ggsave(filename="Boxplot_Density.png",plot=plot_total_quad,device ='png',width = 5, height = 4, dpi = 600, units = "in")
        
      } 
      
      else if (ai[k]=="inverse") {
        # Histogram___________________________
        #png("Histogram.png",res=600,width =4000, height=4000)
        #ggplot(data=obs, aes(inverse(obs[,variables[j]]))) + geom_histogram() + 
        #scale_x_continuous(breaks=seq(0,max(inverse(obs[,variables[j]])),5), lim=c(0,max(inverse(obs[,variables[j]])))) + 
        #labs(x = paste0(varbez))+ggtitle(paste0("Histogram - ",variables[j],"  - (","inverse - ",predictareas_name[i],")")) +
        #theme_minimal() +
        #theme(plot.title = element_text(hjust = 0.5))
        #dev.off()
        #png("Histogram.png",res=600,width =4000, height=4000)
        #histogram(inverse(obs[,variables[j]]),xlab=paste0(varbez,"  - (","inverse - ",predictareas_name[i],")"), col = "gray")
        #dev.off()
        png("Histogram_Normal_Distribution_Curve.png",res=600,width =4000, height=4000)
        m <- mean(inverse(obs[,variables[j]]))
        std <- sqrt(var(inverse(obs[,variables[j]])))
        tmp <- density(inverse(obs[,variables[j]]))
        hist(inverse(obs[,variables[j]]),xlab=paste0(varbez,"  - (","inverse - ",predictareas_name[i],")"), prob=TRUE, main=NULL, col="gray",xlim=c(0,max(tmp$x))) #,ylim=c(0,max(tmp$y)))
        curve(dnorm(x, mean=m, sd=std), col="red", lwd=2, add=TRUE, yaxt="n")
        dev.off()
        # QQ-Plot______________________________
        png("QQPlot.png",res=600,width =4000, height=4000)
        qqnorm(inverse(obs[,variables[j]]),main=paste0("Normal Q-Q Plot - ",varbez," - (","inverse - ",predictareas_name[i],")"))
        qqline(inverse(obs[,variables[j]]))
        dev.off()
        # Boxplot and density__________________
        axis_text_size<-10
        axis_title_size<-10 
        legend_text_size<-10
        table_text_size<-10
        
        th_box<-theme(panel.background = element_rect(fill = NA),
                      panel.grid.major = element_line(linetype = "dashed",colour = "grey"),
                      legend.position = "none" ,
                      axis.title.y=element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks.y=element_blank(),
                      axis.title.x=element_blank(),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())
        
        th_dens<-theme(panel.background = element_rect(fill = NA),
                       panel.grid.major = element_line(linetype = "dashed",colour = "grey"),
                       axis.text.x=element_text(face = "plain",size = axis_text_size,colour = "black"),
                       axis.text.y=element_text(face = "plain",size = axis_text_size,colour = "black"),
                       #axis.title.y = element_text(size = axis_title_size,margin = margin(r=20,l=20)),
                       axis.title.x = element_text(size = axis_title_size),#,margin = margin(r=20,l=20)
                       legend.position = c(0.7, 0.7),
                       legend.title = element_blank(),
                       legend.text = element_text(size = legend_text_size,margin = 10),
                       #legend.key = element_rect(size = 2),
                       legend.key.size = unit(3, 'lines'),
                       axis.line.x = element_line(colour = "black",size=1.2))
        
        o <- ggplot(data=obs, aes_string(paste0("x=inverse(",variables[j],")"))) + geom_boxplot()+ th_box
        a <- ggplot(data=obs, aes_string(paste0("inverse(",variables[j],")"))) + geom_density(fill="blue",alpha=0.5) + th_dens + labs(x = paste0(varbez,"  - (","inverse - ",predictareas_name[i],")"))
        plot_total_inverse<- ggarrange(o, a, heights = c(0.5,1), ncol = 1, nrow =2, align = "v")
        plot_total_inverse
        ggsave(filename="Boxplot_Density.png",plot=plot_total_inverse,device ='png',width = 5, height = 4, dpi = 600, units = "in")
      }
      
      else {}
    }
  }
}


# Overview Summary_Graphics --------------------------------------------------------

# QQplot
for (i in 1:nrpredictareas){
  if (i==1){
    #setwd(output_predictarea_dir1)
    #hi<-dir()[file.info(dir())$isdir]
    for (j in 1:lenvariables){
      ff<-paste0("/output/",predictareas_name[i],"/",variables[j])
      fz<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),ff)
      setwd(fz)
      ho<-dir()[file.info(dir())$isdir]
      lungo<-length(ho)
      rest<-lungo%%2
      ganzdiv<-lungo%/%2
      nrrighe<-rest+ganzdiv
      if (lungo==1){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"QQPlot.png"))
      } else if (lungo==2){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"QQPlot.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"QQPlot.png"))
      } else if (lungo==3){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"QQPlot.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"QQPlot.png"))
        zrttt<-paste0(fz,"/",ho[3],"/Summary_Graphics")
        plot3 <- readPNG(file.path(zrttt,"QQPlot.png"))
      } else if (lungo==4){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"QQPlot.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"QQPlot.png"))
        zrttt<-paste0(fz,"/",ho[3],"/Summary_Graphics")
        plot3 <- readPNG(file.path(zrttt,"QQPlot.png"))
        zrtttt<-paste0(fz,"/",ho[4],"/Summary_Graphics")
        plot4 <- readPNG(file.path(zrtttt,"QQPlot.png"))
      } else if (lungo==5){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"QQPlot.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"QQPlot.png"))
        zrttt<-paste0(fz,"/",ho[3],"/Summary_Graphics")
        plot3 <- readPNG(file.path(zrttt,"QQPlot.png"))
        zrtttt<-paste0(fz,"/",ho[4],"/Summary_Graphics")
        plot4 <- readPNG(file.path(zrtttt,"QQPlot.png"))
        zrttttt<-paste0(fz,"/",ho[5],"/Summary_Graphics")
        plot5 <- readPNG(file.path(zrttttt,"QQPlot.png"))
      } else {}
      png(file.path(fz,"Overview_QQPlot.png"),res=600,width =8000, height=(4000*nrrighe))
      if (lungo==1){
        grid.arrange(rasterGrob(plot1),nrow=nrrighe, ncol=2)
      } else if (lungo==2){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2),nrow=nrrighe, ncol=2)
      } else if (lungo==3){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2), rasterGrob(plot3),nrow=nrrighe, ncol=2)
      } else if (lungo==4){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2), rasterGrob(plot3),rasterGrob(plot4),nrow=nrrighe, ncol=2)
      } else if (lungo==5){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2), rasterGrob(plot3),rasterGrob(plot4),rasterGrob(plot5),nrow=nrrighe, ncol=2)
      } else {}
      dev.off()
    }
  }
  if (i==2){
    #setwd(output_predictarea_dir1)
    #hi<-dir()[file.info(dir())$isdir]
    for (j in 1:lenvariables){
      ff<-paste0("/output/",predictareas_name[i],"/",variables[j])
      fz<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),ff)
      setwd(fz)
      ho<-dir()[file.info(dir())$isdir]
      lungo<-length(ho)
      rest<-lungo%%2
      ganzdiv<-lungo%/%2
      nrrighe<-rest+ganzdiv
      if (lungo==1){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"QQPlot.png"))
      } else if (lungo==2){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"QQPlot.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"QQPlot.png"))
      } else if (lungo==3){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"QQPlot.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"QQPlot.png"))
        zrttt<-paste0(fz,"/",ho[3],"/Summary_Graphics")
        plot3 <- readPNG(file.path(zrttt,"QQPlot.png"))
      } else if (lungo==4){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"QQPlot.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"QQPlot.png"))
        zrttt<-paste0(fz,"/",ho[3],"/Summary_Graphics")
        plot3 <- readPNG(file.path(zrttt,"QQPlot.png"))
        zrtttt<-paste0(fz,"/",ho[4],"/Summary_Graphics")
        plot4 <- readPNG(file.path(zrtttt,"QQPlot.png"))
      } else if (lungo==5){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"QQPlot.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"QQPlot.png"))
        zrttt<-paste0(fz,"/",ho[3],"/Summary_Graphics")
        plot3 <- readPNG(file.path(zrttt,"QQPlot.png"))
        zrtttt<-paste0(fz,"/",ho[4],"/Summary_Graphics")
        plot4 <- readPNG(file.path(zrtttt,"QQPlot.png"))
        zrttttt<-paste0(fz,"/",ho[5],"/Summary_Graphics")
        plot5 <- readPNG(file.path(zrttttt,"QQPlot.png"))
      } else {}
      png(file.path(fz,"Overview_QQPlot.png"),res=600,width =8000, height=(4000*nrrighe))
      if (lungo==1){
        grid.arrange(rasterGrob(plot1),nrow=nrrighe, ncol=2)
      } else if (lungo==2){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2),nrow=nrrighe, ncol=2)
      } else if (lungo==3){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2), rasterGrob(plot3),nrow=nrrighe, ncol=2)
      } else if (lungo==4){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2), rasterGrob(plot3),rasterGrob(plot4),nrow=nrrighe, ncol=2)
      } else if (lungo==5){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2), rasterGrob(plot3),rasterGrob(plot4),rasterGrob(plot5),nrow=nrrighe, ncol=2)
      } else {}
      dev.off()
    }
  }
}

# Histogram + Normal distribution curve
for (i in 1:nrpredictareas){
  if (i==1){
    #setwd(output_predictarea_dir1)
    #hi<-dir()[file.info(dir())$isdir]
    for (j in 1:lenvariables){
      ff<-paste0("/output/",predictareas_name[i],"/",variables[j])
      fz<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),ff)
      setwd(fz)
      ho<-dir()[file.info(dir())$isdir]
      lungo<-length(ho)
      rest<-lungo%%2
      ganzdiv<-lungo%/%2
      nrrighe<-rest+ganzdiv
      if (lungo==1){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"Histogram_Normal_Distribution_Curve.png"))
      } else if (lungo==2){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"Histogram_Normal_Distribution_Curve.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"Histogram_Normal_Distribution_Curve.png"))
      } else if (lungo==3){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"Histogram_Normal_Distribution_Curve.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"Histogram_Normal_Distribution_Curve.png"))
        zrttt<-paste0(fz,"/",ho[3],"/Summary_Graphics")
        plot3 <- readPNG(file.path(zrttt,"Histogram_Normal_Distribution_Curve.png"))
      } else if (lungo==4){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"Histogram_Normal_Distribution_Curve.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"Histogram_Normal_Distribution_Curve.png"))
        zrttt<-paste0(fz,"/",ho[3],"/Summary_Graphics")
        plot3 <- readPNG(file.path(zrttt,"Histogram_Normal_Distribution_Curve.png"))
        zrtttt<-paste0(fz,"/",ho[4],"/Summary_Graphics")
        plot4 <- readPNG(file.path(zrtttt,"Histogram_Normal_Distribution_Curve.png"))
      } else if (lungo==5){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"Histogram_Normal_Distribution_Curve.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"Histogram_Normal_Distribution_Curve.png"))
        zrttt<-paste0(fz,"/",ho[3],"/Summary_Graphics")
        plot3 <- readPNG(file.path(zrttt,"Histogram_Normal_Distribution_Curve.png"))
        zrtttt<-paste0(fz,"/",ho[4],"/Summary_Graphics")
        plot4 <- readPNG(file.path(zrtttt,"Histogram_Normal_Distribution_Curve.png"))
        zrttttt<-paste0(fz,"/",ho[5],"/Summary_Graphics")
        plot5 <- readPNG(file.path(zrttttt,"Histogram_Normal_Distribution_Curve.png"))
      } else {}
      png(file.path(fz,"Overview_Histogram.png"),res=600,width =8000, height=(4000*nrrighe))
      if (lungo==1){
        grid.arrange(rasterGrob(plot1),nrow=nrrighe, ncol=2)
      } else if (lungo==2){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2),nrow=nrrighe, ncol=2)
      } else if (lungo==3){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2), rasterGrob(plot3),nrow=nrrighe, ncol=2)
      } else if (lungo==4){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2), rasterGrob(plot3),rasterGrob(plot4),nrow=nrrighe, ncol=2)
      } else if (lungo==5){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2), rasterGrob(plot3),rasterGrob(plot4),rasterGrob(plot5),nrow=nrrighe, ncol=2)
      } else {}
      dev.off()
    }
  }
  if (i==2){
    #setwd(output_predictarea_dir1)
    #hi<-dir()[file.info(dir())$isdir]
    for (j in 1:lenvariables){
      ff<-paste0("/output/",predictareas_name[i],"/",variables[j])
      fz<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),ff)
      setwd(fz)
      ho<-dir()[file.info(dir())$isdir]
      lungo<-length(ho)
      rest<-lungo%%2
      ganzdiv<-lungo%/%2
      nrrighe<-rest+ganzdiv
      if (lungo==1){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"Histogram_Normal_Distribution_Curve.png"))
      } else if (lungo==2){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"Histogram_Normal_Distribution_Curve.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"Histogram_Normal_Distribution_Curve.png"))
      } else if (lungo==3){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"Histogram_Normal_Distribution_Curve.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"Histogram_Normal_Distribution_Curve.png"))
        zrttt<-paste0(fz,"/",ho[3],"/Summary_Graphics")
        plot3 <- readPNG(file.path(zrttt,"Histogram_Normal_Distribution_Curve.png"))
      } else if (lungo==4){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"Histogram_Normal_Distribution_Curve.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"Histogram_Normal_Distribution_Curve.png"))
        zrttt<-paste0(fz,"/",ho[3],"/Summary_Graphics")
        plot3 <- readPNG(file.path(zrttt,"Histogram_Normal_Distribution_Curve.png"))
        zrtttt<-paste0(fz,"/",ho[4],"/Summary_Graphics")
        plot4 <- readPNG(file.path(zrtttt,"Histogram_Normal_Distribution_Curve.png"))
      } else if (lungo==5){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"Histogram_Normal_Distribution_Curve.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"Histogram_Normal_Distribution_Curve.png"))
        zrttt<-paste0(fz,"/",ho[3],"/Summary_Graphics")
        plot3 <- readPNG(file.path(zrttt,"Histogram_Normal_Distribution_Curve.png"))
        zrtttt<-paste0(fz,"/",ho[4],"/Summary_Graphics")
        plot4 <- readPNG(file.path(zrtttt,"Histogram_Normal_Distribution_Curve.png"))
        zrttttt<-paste0(fz,"/",ho[5],"/Summary_Graphics")
        plot5 <- readPNG(file.path(zrttttt,"Histogram_Normal_Distribution_Curve.png"))
      } else {}
      png(file.path(fz,"Overview_Histogram.png"),res=600,width =8000, height=(4000*nrrighe))
      if (lungo==1){
        grid.arrange(rasterGrob(plot1),nrow=nrrighe, ncol=2)
      } else if (lungo==2){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2),nrow=nrrighe, ncol=2)
      } else if (lungo==3){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2), rasterGrob(plot3),nrow=nrrighe, ncol=2)
      } else if (lungo==4){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2), rasterGrob(plot3),rasterGrob(plot4),nrow=nrrighe, ncol=2)
      } else if (lungo==5){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2), rasterGrob(plot3),rasterGrob(plot4),rasterGrob(plot5),nrow=nrrighe, ncol=2)
      } else {}
      dev.off()
    }
  }
}

# Boxplot and density
for (i in 1:nrpredictareas){
  if (i==1){
    #setwd(output_predictarea_dir1)
    #hi<-dir()[file.info(dir())$isdir]
    for (j in 1:lenvariables){
      ff<-paste0("/output/",predictareas_name[i],"/",variables[j])
      fz<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),ff)
      setwd(fz)
      ho<-dir()[file.info(dir())$isdir]
      lungo<-length(ho)
      rest<-lungo%%2
      ganzdiv<-lungo%/%2
      nrrighe<-rest+ganzdiv
      if (lungo==1){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"Boxplot_Density.png"))
      } else if (lungo==2){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"Boxplot_Density.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"Boxplot_Density.png"))
      } else if (lungo==3){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"Boxplot_Density.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"Boxplot_Density.png"))
        zrttt<-paste0(fz,"/",ho[3],"/Summary_Graphics")
        plot3 <- readPNG(file.path(zrttt,"Boxplot_Density.png"))
      } else if (lungo==4){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"Boxplot_Density.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"Boxplot_Density.png"))
        zrttt<-paste0(fz,"/",ho[3],"/Summary_Graphics")
        plot3 <- readPNG(file.path(zrttt,"Boxplot_Density.png"))
        zrtttt<-paste0(fz,"/",ho[4],"/Summary_Graphics")
        plot4 <- readPNG(file.path(zrtttt,"Boxplot_Density.png"))
      } else if (lungo==5){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"Boxplot_Density.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"Boxplot_Density.png"))
        zrttt<-paste0(fz,"/",ho[3],"/Summary_Graphics")
        plot3 <- readPNG(file.path(zrttt,"Boxplot_Density.png"))
        zrtttt<-paste0(fz,"/",ho[4],"/Summary_Graphics")
        plot4 <- readPNG(file.path(zrtttt,"Boxplot_Density.png"))
        zrttttt<-paste0(fz,"/",ho[5],"/Summary_Graphics")
        plot5 <- readPNG(file.path(zrttttt,"Boxplot_Density.png"))
      } else {}
      png(file.path(fz,"Overview_BoxplotDensity.png"),res=600,width =6000, height=(2400*nrrighe))
      if (lungo==1){
        grid.arrange(rasterGrob(plot1),nrow=nrrighe, ncol=2)
      } else if (lungo==2){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2),nrow=nrrighe, ncol=2)
      } else if (lungo==3){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2), rasterGrob(plot3),nrow=nrrighe, ncol=2)
      } else if (lungo==4){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2), rasterGrob(plot3),rasterGrob(plot4),nrow=nrrighe, ncol=2)
      } else if (lungo==5){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2), rasterGrob(plot3),rasterGrob(plot4),rasterGrob(plot5),nrow=nrrighe, ncol=2)
      } else {}
      dev.off()
    }
  }
  if (i==2){
    #setwd(output_predictarea_dir1)
    #hi<-dir()[file.info(dir())$isdir]
    for (j in 1:lenvariables){
      ff<-paste0("/output/",predictareas_name[i],"/",variables[j])
      fz<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),ff)
      setwd(fz)
      ho<-dir()[file.info(dir())$isdir]
      lungo<-length(ho)
      rest<-lungo%%2
      ganzdiv<-lungo%/%2
      nrrighe<-rest+ganzdiv
      if (lungo==1){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"Boxplot_Density.png"))
      } else if (lungo==2){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"Boxplot_Density.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"Boxplot_Density.png"))
      } else if (lungo==3){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"Boxplot_Density.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"Boxplot_Density.png"))
        zrttt<-paste0(fz,"/",ho[3],"/Summary_Graphics")
        plot3 <- readPNG(file.path(zrttt,"Boxplot_Density.png"))
      } else if (lungo==4){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"Boxplot_Density.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"Boxplot_Density.png"))
        zrttt<-paste0(fz,"/",ho[3],"/Summary_Graphics")
        plot3 <- readPNG(file.path(zrttt,"Boxplot_Density.png"))
        zrtttt<-paste0(fz,"/",ho[4],"/Summary_Graphics")
        plot4 <- readPNG(file.path(zrtttt,"Boxplot_Density.png"))
      } else if (lungo==5){
        zrt<-paste0(fz,"/",ho[1],"/Summary_Graphics")
        plot1 <- readPNG(file.path(zrt,"Boxplot_Density.png"))
        zrtt<-paste0(fz,"/",ho[2],"/Summary_Graphics")
        plot2 <- readPNG(file.path(zrtt,"Boxplot_Density.png"))
        zrttt<-paste0(fz,"/",ho[3],"/Summary_Graphics")
        plot3 <- readPNG(file.path(zrttt,"Boxplot_Density.png"))
        zrtttt<-paste0(fz,"/",ho[4],"/Summary_Graphics")
        plot4 <- readPNG(file.path(zrtttt,"Boxplot_Density.png"))
        zrttttt<-paste0(fz,"/",ho[5],"/Summary_Graphics")
        plot5 <- readPNG(file.path(zrttttt,"Boxplot_Density.png"))
      } else {}
      png(file.path(fz,"Overview_BoxplotDensity.png"),res=600,width =6000, height=(2400*nrrighe))
      if (lungo==1){
        grid.arrange(rasterGrob(plot1),nrow=nrrighe, ncol=2)
      } else if (lungo==2){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2),nrow=nrrighe, ncol=2)
      } else if (lungo==3){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2), rasterGrob(plot3),nrow=nrrighe, ncol=2)
      } else if (lungo==4){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2), rasterGrob(plot3),rasterGrob(plot4),nrow=nrrighe, ncol=2)
      } else if (lungo==5){
        grid.arrange(rasterGrob(plot1), rasterGrob(plot2), rasterGrob(plot3),rasterGrob(plot4),rasterGrob(plot5),nrow=nrrighe, ncol=2)
      } else {}
      dev.off()
    }
  }
}    

# Correlation coefficients ------------------------------------------------

# Read Masterfile - Soil Data Points
obs<- read.csv(file = file.path(datapoints_dir,list.files(path = datapoints_dir)),header = TRUE,sep = ",",dec = ".",na.strings = c("-9999"))  
obs<- obs[!duplicated(obs[,c("x","y")]),]
obsAllCorr1<-obs
obsAllCorr2<-obs
coordinates(obs)= ~x+y
coordinates(obsAllCorr1)=~x+y
coordinates(obsAllCorr2)=~x+y


for (i in 1:nrpredictareas){
  if (i==1){
    predictors = list.dirs(path = predictors_dir_1)
    predictors = predictors[grepl(pattern = resolution,x = predictors)]
    zu<-paste0(predictareas_name[1])
    landuse<-raster(file.path(predict_area_dir_1,zu))  
    cropped_rasters = lapply(predictors, function(x) resample(raster(x),landuse, method="ngb"))
    r = raster::stack(cropped_rasters)
    r <- raster::mask(r,landuse)
    df=as.data.frame(r)
    df=df[complete.cases(df), ]
    gridCoords<-xyFromCell(object = r,which(complete.cases(as.data.frame(r))))
    grid=SpatialPixelsDataFrame(points = gridCoords,proj4string = CRS(proj4string(r)),data = df)
    proj4string(obs)= proj4string(grid)
    #proj4string(obsValid)= proj4string(grid)
    obs@data<-cbind(obs@data,raster::extract(r, obs))
    obsAllCorr1@data<-cbind(obsAllCorr1@data,raster::extract(r,obsAllCorr1))
    
    # ObsAllCorr@data => find numeric variables
    numericvariables1<-c()
    for (z in 9:length(obsAllCorr1@data)){
      prova<-unlist(obsAllCorr1@data[z])
      prova2<-as.numeric(prova)
      prova2<-prova2[which(!is.na(prova2))]
      prova2 <- prova2[ prova2 != 0]
      prova2 <- prova2[ prova2 != 1]
      if (is_empty(prova2)=="FALSE"){
        hilf<-length(numericvariables1)
        numericvariables1 <- append(numericvariables1,z,after=hilf) # Numeric columns
      } else {}
    }
    
    for (u in length(numericvariables1):1){
      if (u==length(numericvariables1)){
        corr_parameters1 <- obsAllCorr1@data[numericvariables1[length(numericvariables1)]]
      } else {
        corr_parameters1 <- cbind(obsAllCorr1@data[numericvariables1[u]],corr_parameters1)
      }
    }
    
    # Denomination of variables
    names(corr_parameters1)[names(corr_parameters1)== "P_1_aspect_ea_150"] <- "P1_aspect east"
    names(corr_parameters1)[names(corr_parameters1)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
    names(corr_parameters1)[names(corr_parameters1)== "P_3_aspect_no_150"] <- "P3_aspect north"
    names(corr_parameters1)[names(corr_parameters1)== "P_4_aspect_so_150"] <- "P4_aspect south"
    names(corr_parameters1)[names(corr_parameters1)== "P_5_aspect_we_150"] <- "P5_aspect west"
    names(corr_parameters1)[names(corr_parameters1)== "P_6_cl_grpre1_150"] <- "P6_land use - alpine meadow"
    names(corr_parameters1)[names(corr_parameters1)== "P_7_cl_grpre2_150"] <- "P7_land use - grassland"
    names(corr_parameters1)[names(corr_parameters1)== "P_8_cl_grpre3_150"] <- "P8_land use - grassland second shearing"
    names(corr_parameters1)[names(corr_parameters1)== "P_9_cl_grpre4_150"] <- "P9_land use - grassland special area"
    names(corr_parameters1)[names(corr_parameters1)== "P_10_dtm_150m"] <- "P10_dtm"
    names(corr_parameters1)[names(corr_parameters1)== "P_11_geo_bin_1_150"] <- "P11_geology - shist"
    names(corr_parameters1)[names(corr_parameters1)== "P_12_geo_bin_2_150"] <- "P12_geology - vulcanite"
    names(corr_parameters1)[names(corr_parameters1)== "P_13_geo_bin_3_150"] <- "P13_geology - vinschgau shear zone"
    names(corr_parameters1)[names(corr_parameters1)== "P_14_geo_bin_4_150"] <- "P14_geology - quaternary depositions"
    names(corr_parameters1)[names(corr_parameters1)== "P_15_geo_bin_5_150"] <- "P15_geology - sediment sequences"
    names(corr_parameters1)[names(corr_parameters1)== "P_16_geo_bin_6_150"] <- "P16_geology - plutons"
    names(corr_parameters1)[names(corr_parameters1)== "P_17_geo_bin_7_150"] <- "P17_geology - central gneiss"
    names(corr_parameters1)[names(corr_parameters1)== "P_18_geo_bin_8_150"] <- "P18_geology - quartz phyllite"
    names(corr_parameters1)[names(corr_parameters1)== "P_19_geo_bin_9_150"] <- "P19_geology - insignificant alpine metamorphosis"
    names(corr_parameters1)[names(corr_parameters1)== "P_20_geo_bi_10_150"] <- "P20_geology - medium-low alpine metamorphosis"
    names(corr_parameters1)[names(corr_parameters1)== "P_21_geo_bi_11_150"] <- "P21_geology - medium-high alpine metamorphosis"
    names(corr_parameters1)[names(corr_parameters1)== "P_22_slope_150m"] <- "P22_slope"
    names(corr_parameters1)[names(corr_parameters1)== "P_23_wetness_150m"] <- "P23_wetness_index"
    names(corr_parameters1)[names(corr_parameters1)== "P_24_gve_wgs84_150"] <- "P24_GVE"
    names(corr_parameters1)[names(corr_parameters1)== "P_25_geo_aeubin150"] <- "P25_Acidic soil"
    names(corr_parameters1)[names(corr_parameters1)== "P_26_geo_beubin150"] <- "P26_Basic soil"
    names(corr_parameters1)[names(corr_parameters1)== "P_27_geo_aelbin150"] <- "P27_Acidic soil"
    names(corr_parameters1)[names(corr_parameters1)== "P_28_geo_belbin150"] <- "P28_Basic soil"
    names(corr_parameters1)[names(corr_parameters1)== "P_29_geo_abelbi150"] <- "P29_Acidic-Basic soil"
    
    matrix_spearman <- cor(corr_parameters1,use="complete.obs",method="spearman")
    matrix_pearson <- cor(corr_parameters1,use="complete.obs",method="pearson")
    spearman = round(matrix_spearman, 2)
    write.csv(spearman,file.path(output_predictarea_dir1,"Spearman_Correlation_Predictarea1.csv"))      
    write.xlsx(spearman,file.path(output_predictarea_dir1,"Spearman_Correlation_Predictarea1.xlsx"))    
    png(file.path(output_predictarea_dir1,"Spearman_Correlation_Predictarea1.png"),res=600,width =7000, height=7000)
    corrplot(matrix_spearman, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
    dev.off()
    pearson = round(matrix_pearson, 2)
    write.csv(pearson,file.path(output_predictarea_dir1,"Pearson_Correlation_Predictarea1.csv"))      # ignore error
    write.xlsx(pearson,file.path(output_predictarea_dir1,"Pearson_Correlation_Predictarea1.xlsx"))    # ignore error
    png(file.path(output_predictarea_dir1,"Pearson_Correlation_Predictarea1.png"),res=600,width =7000, height=7000)
    corrplot(matrix_pearson, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
    dev.off()
    
    
  } else if (i==2){
    predictors = list.dirs(path = predictors_dir_2)
    predictors = predictors[grepl(pattern = resolution,x = predictors)]
    zu<-paste0(predictareas_name[2])
    landuse<-raster(file.path(predict_area_dir_2,zu))  
    cropped_rasters = lapply(predictors, function(x) resample(raster(x),landuse, method="ngb"))
    r = raster::stack(cropped_rasters)
    r <- raster::mask(r,landuse)
    df=as.data.frame(r)
    df=df[complete.cases(df), ]
    gridCoords<-xyFromCell(object = r,which(complete.cases(as.data.frame(r))))
    grid=SpatialPixelsDataFrame(points = gridCoords,proj4string = CRS(proj4string(r)),data = df)
    proj4string(obs)= proj4string(grid)
    #proj4string(obsValid)= proj4string(grid)
    obs@data<-cbind(obs@data,raster::extract(r, obs))
    obsAllCorr2@data<-cbind(obsAllCorr2@data,raster::extract(r,obsAllCorr2))
    
    # ObsAllCorr@data => find numeric variables
    numericvariables2<-c()
    for (z in 9:length(obsAllCorr2@data)){
      prova<-unlist(obsAllCorr2@data[z])
      prova2<-as.numeric(prova)
      prova2<-prova2[which(!is.na(prova2))]
      prova2 <- prova2[ prova2 != 0]
      prova2 <- prova2[ prova2 != 1]
      if (is_empty(prova2)=="FALSE"){
        hilf<-length(numericvariables2)
        numericvariables2 <- append(numericvariables2,z,after=hilf) # Numeric columns
      } else {}
    }
    
    for (u in length(numericvariables2):1){
      if (u==length(numericvariables2)){
        corr_parameters2 <- obsAllCorr2@data[numericvariables2[length(numericvariables2)]]
      } else {
        corr_parameters2 <- cbind(obsAllCorr2@data[numericvariables2[u]],corr_parameters2)
      }
    }
    
    # Denomination of variables
    names(corr_parameters2)[names(corr_parameters2)== "P_1_aspect_ea_150"] <- "P1_aspect east"
    names(corr_parameters2)[names(corr_parameters2)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
    names(corr_parameters2)[names(corr_parameters2)== "P_3_aspect_no_150"] <- "P3_aspect north"
    names(corr_parameters2)[names(corr_parameters2)== "P_4_aspect_so_150"] <- "P4_aspect south"
    names(corr_parameters2)[names(corr_parameters2)== "P_5_aspect_we_150"] <- "P5_aspect west"
    names(corr_parameters2)[names(corr_parameters2)== "P_6_cl_grarp1_150"] <- "P6_land use - alpine meadow"
    names(corr_parameters2)[names(corr_parameters2)== "P_7_cl_grarp2_150"] <- "P7_land use - grassland"
    names(corr_parameters2)[names(corr_parameters2)== "P_8_cl_grarp3_150"] <- "P8_land use - grassland second shearing"
    names(corr_parameters2)[names(corr_parameters2)== "P_9_cl_grarp4_150"] <- "P9land use - grassland special area"
    names(corr_parameters2)[names(corr_parameters2)== "P_10_cl_grarp5_150"] <- "P10_land use - crops"
    names(corr_parameters2)[names(corr_parameters2)== "P_11_dtm_150m"] <- "P11_dtm"
    names(corr_parameters2)[names(corr_parameters2)== "P_12_geo_bin_1_150"] <- "P12_geology - shist"
    names(corr_parameters2)[names(corr_parameters2)== "P_13_geo_bin_2_150"] <- "P13_geology - vulcanite"
    names(corr_parameters2)[names(corr_parameters2)== "P_14_geo_bin_3_150"] <- "P14_geology - vinschgau shear zone"
    names(corr_parameters2)[names(corr_parameters2)== "P_15_geo_bin_4_150"] <- "P15_geology - quaternary depositions"
    names(corr_parameters2)[names(corr_parameters2)== "P_16_geo_bin_5_150"] <- "P16_geology - sediment sequences"
    names(corr_parameters2)[names(corr_parameters2)== "P_17_geo_bin_6_150"] <- "P17_geology - plutons"
    names(corr_parameters2)[names(corr_parameters2)== "P_18_geo_bin_7_150"] <- "P18_geology - central gneiss"
    names(corr_parameters2)[names(corr_parameters2)== "P_19_geo_bin_8_150"] <- "P19_geology - quartz phyllite"
    names(corr_parameters2)[names(corr_parameters2)== "P_20_geo_bin_9_150"] <- "P20_geology - insignificant alpine metamorphosis"
    names(corr_parameters2)[names(corr_parameters2)== "P_21_geo_bi_10_150"] <- "P21_geology - medium-low alpine metamorphosis"
    names(corr_parameters2)[names(corr_parameters2)== "P_22_geo_bi_11_150"] <- "P22_geology - medium-high alpine metamorphosis"
    names(corr_parameters2)[names(corr_parameters2)== "P_23_slope_150m"] <- "P23_slope"
    names(corr_parameters2)[names(corr_parameters2)== "P_24_wetness_150m"] <- "P24_wetness_index"
    names(corr_parameters2)[names(corr_parameters2)== "P_25_gve_wgs84_150"] <- "P25_GVE"
    names(corr_parameters2)[names(corr_parameters2)== "P_26_geo_aeubin150"] <- "P26_Acidic soil"
    names(corr_parameters2)[names(corr_parameters2)== "P_27_geo_beubin150"] <- "P27_Basic soil"
    names(corr_parameters2)[names(corr_parameters2)== "P_28_geo_aelbin150"] <- "P28_Acidic soil"
    names(corr_parameters2)[names(corr_parameters2)== "P_29_geo_belbin150"] <- "P29_Basic soil"
    names(corr_parameters2)[names(corr_parameters2)== "P_30_geo_abelbi150"] <- "P30_Acidic-Basic soil"
    
    matrix_spearman <- cor(corr_parameters2,use="complete.obs",method="spearman")
    matrix_pearson <- cor(corr_parameters2,use="complete.obs",method="pearson")
    spearman = round(matrix_spearman, 2)
    write.csv(spearman,file.path(output_predictarea_dir2,"Spearman_Correlation_Predictarea2.csv"))      
    write.xlsx(spearman,file.path(output_predictarea_dir2,"Spearman_Correlation_Predictarea2.xlsx"))    
    png(file.path(output_predictarea_dir2,"Spearman_Correlation_Predictarea2.png"),res=600,width =7000, height=7000)
    corrplot(matrix_spearman, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
    dev.off()
    pearson = round(matrix_pearson, 2)
    write.csv(pearson,file.path(output_predictarea_dir2,"Pearson_Correlation_Predictarea2.csv"))      # ignore error
    write.xlsx(pearson,file.path(output_predictarea_dir2,"Pearson_Correlation_Predictarea2.xlsx"))    # ignore error
    png(file.path(output_predictarea_dir2,"Pearson_Correlation_Predictarea2.png"),res=600,width =7000, height=7000)
    corrplot(matrix_pearson, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
    dev.off()
    
  } else{}
}



# MAIN CODE (Random Forest/Regression - Ordinary Kriging/Prediction --------

nrun<-1
spalte7tot<-vector()
spalte10tot<-vector()

# Main Code
for (j in 1:lenvariables){   #(j in 1:lenvariables)
  print(j)
  print(variables[j])
  var<-variables[j]
  setwd(paste0(script_dir,"/output_overview"))
  dir.create(paste0("Overview_",variables[j]))
  for (i in 1:nrpredictareas){  #(i in 1:nrpredictareas)
    spalte1<-vector()
    spalte2<-vector()
    spalte3<-vector()
    spalte4<-vector()
    spalte5<-vector()
    spalte6<-vector()
    spalte7<-vector()
    spalte8<-vector()
    spalte9<-vector()
    spalte10<-vector()
    print(i)
    print(predictareas_name[i])
    ff<-paste0("/output/",predictareas_name[i],"/",variables[j])
    #fz<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),ff)
    fz<-paste0(script_dir,ff)
    setwd(fz)
    ai<-dir()[file.info(dir())$isdir]
    lange<-length(ai)
    # 2_Input_PredictArea_Predictors
    if (i==1){
      # PREDICTAREA 1 INPUT
      zu<-paste0(predictareas_name[1])
      landuse<-raster(file.path(predict_area_dir_1,zu))
      # PREDICTORS FOR PREDICTAREA 1
      predictors = list.dirs(path = predictors_dir_1)
      predictors = predictors[grepl(pattern = resolution,x = predictors)]
      if (predictors_combination=="personalized") {
        agj<-c(1:length(predictors))
        nrofcombinations<-26
        predictors_combi_1<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)
        predictors_combi_2<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)
        predictors_combi_3<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
        predictors_combi_4<-c(1,2,3,4,5,6,7,8,9,10,22,23,25,26)
        predictors_combi_5<-c(1,2,3,4,5,6,7,8,9,10,22,23,27,28,29)
        predictors_combi_6<-c(1,2,3,4,5,6,7,8,9,10,22,23,24,25,26)
        predictors_combi_7<-c(1,2,3,4,5,6,7,8,9,10,22,23,24,27,28,29)
        predictors_combi_8<-c(1,2,3,4,5,6,7,8,9,10,22,23,24)
        predictors_combi_9<-c(1,2,3,4,5,6,7,8,9,10,22,23)
        predictors_combi_10<-c(1,2,3,4,5,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
        predictors_combi_11<-c(1,2,3,4,5,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)
        predictors_combi_12<-c(1,2,3,4,5,10,22,23)
        predictors_combi_13<-c(1,2,3,4,5,10,22,23,24)
        predictors_combi_14<-c(1,3,4,5,6,7,9,10,11,12,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29)
        predictors_combi_15<-c(1,3,4,5,6,7,9,10,11,12,14,15,17,18,19,20,21,22,23,24)
        predictors_combi_16<-c(1,3,4,5,6,7,9,10,11,12,14,15,17,18,19,20,21,22,23) 
        predictors_combi_17<-c(1,3,4,5,6,7,9,10,22,23,25,26)
        predictors_combi_18<-c(1,3,4,5,6,7,9,10,22,23,27,28,29)
        predictors_combi_19<-c(1,3,4,5,6,7,9,10,22,23,24,25,26)
        predictors_combi_20<-c(1,3,4,5,6,7,9,10,22,23,24,27,28,29)
        predictors_combi_21<-c(1,3,4,5,6,7,9,10,22,23,24)
        predictors_combi_22<-c(1,3,4,5,6,7,9,10,22,23)
        predictors_combi_23<-c(1,3,4,5,10,11,12,14,15,17,18,19,20,21,22,23)
        predictors_combi_24<-c(1,3,4,5,10,11,12,14,15,17,18,19,20,21,22,23,24)
        predictors_combi_25<-c(1,3,4,5,10,22,23)
        predictors_combi_26<-c(1,3,4,5,10,22,23,24)
        my_different_combinations_of_predictors<-matrix(list(), nrow=1, ncol=nrofcombinations)
        for (ws in 1:nrofcombinations) {
          my_different_combinations_of_predictors[[1,ws]]<-get(paste0("predictors_combi_",ws))
        }
      } else if (predictors_combination=="automatic"){
        memory.limit(size = 9999999999999)
        stefanie<-c(1:length(predictors))
        ä<-all_combn(stefanie)
        resultelvis<-str_split(ä," and ")
        agj<-c(1:length(predictors))
        nrofcombinations<-length(ä)
      } else {}
    } else if (i==2){
      # PREDICTAREA 2 INPUT
      zu<-paste0(predictareas_name[2])
      landuse<-raster(file.path(predict_area_dir_2,zu))
      # PREDICTORS FOR PREDICTAREA 2
      predictors = list.dirs(path = predictors_dir_2)
      predictors = predictors[grepl(pattern = resolution,x = predictors)]
      if (predictors_combination=="personalized") {
        agj<-c(1:length(predictors))
        nrofcombinations<-26
        predictors_combi_1<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30)
        predictors_combi_2<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)
        predictors_combi_3<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)
        predictors_combi_4<-c(1,2,3,4,5,6,7,8,9,10,11,23,24,26,27)
        predictors_combi_5<-c(1,2,3,4,5,6,7,8,9,10,11,23,24,28,29,30)
        predictors_combi_6<-c(1,2,3,4,5,6,7,8,9,10,11,23,24,25,26,27)
        predictors_combi_7<-c(1,2,3,4,5,6,7,8,9,10,11,23,24,25,28,29,30)
        predictors_combi_8<-c(1,2,3,4,5,6,7,8,9,10,11,23,24,25)
        predictors_combi_9<-c(1,2,3,4,5,6,7,8,9,10,11,23,24)
        predictors_combi_10<-c(1,2,3,4,5,11,12,13,14,15,16,17,18,19,20,21,22,23,24)
        predictors_combi_11<-c(1,2,3,4,5,11,12,13,14,15,16,17,18,19,20,21,22,23,24)
        predictors_combi_12<-c(1,2,3,4,5,11,23,24)
        predictors_combi_13<-c(1,2,3,4,5,11,23,24,25)
        predictors_combi_14<-c(1,3,4,5,6,7,9,10,11,12,13,15,16,18,19,20,21,22,23,24,25,26,27,28,29,30)
        predictors_combi_15<-c(1,3,4,5,6,7,9,10,11,12,13,15,16,18,19,20,21,22,23,24,25)
        predictors_combi_16<-c(1,3,4,5,6,7,9,10,11,12,13,15,16,18,19,20,21,22,23,24)
        predictors_combi_17<-c(1,3,4,5,6,7,9,10,11,23,24,26,27)
        predictors_combi_18<-c(1,3,4,5,6,7,9,10,11,23,24,28,29,30)
        predictors_combi_19<-c(1,3,4,5,6,7,9,10,11,23,24,25,26,27)
        predictors_combi_20<-c(1,3,4,5,6,7,9,10,11,23,24,25,28,29,30)
        predictors_combi_21<-c(1,3,4,5,6,7,9,10,11,23,24,25)
        predictors_combi_22<-c(1,3,4,5,6,7,9,10,11,23,24)
        predictors_combi_23<-c(1,3,4,5,11,12,13,15,16,18,19,20,21,22,23,24)
        predictors_combi_24<-c(1,3,4,5,11,12,13,15,16,18,19,20,21,22,23,24,25)
        predictors_combi_25<-c(1,3,4,5,11,23,24)
        predictors_combi_26<-c(1,3,4,5,11,23,24,25)
        my_different_combinations_of_predictors<-matrix(list(), nrow=1, ncol=nrofcombinations)
        for (ws in 1:nrofcombinations) {
          my_different_combinations_of_predictors[[1,ws]]<-get(paste0("predictors_combi_",ws))
        }
      } else if (predictors_combination=="automatic"){
        memory.limit(size = 9999999999999)
        stefanie<-c(1:length(predictors))
        ä<-all_combn(stefanie)
        resultelvis<-str_split(ä," and ")
        agj<-c(1:length(predictors))
        nrofcombinations<-length(ä)
      } else {}
    } else{}
    for (k in 1:lange){ #(k in 1:lange)
      print(k)
      print(ai[k])
      for (s in 1:3) { # (s in 1:3)
        print(paste0("s=",s))
        # Read Masterfile - Soil Data Points
        obs<- read.csv(file = file.path(datapoints_dir,list.files(path = datapoints_dir)),header = TRUE,sep = ",",dec = ".",na.strings = c("-9999"))  
        obs<- obs[!duplicated(obs[,c("x","y")]),]
        # Clean Masterfile from missing values
        obs<-obs[which(!is.na(obs[,var])),]
        # Divide data - Training 80% vs. Validation 20%
        part<-createDataPartition(y = obs[,var],times = 1,p = 0.8,list = F)
        obsTrain<- obs[part,]
        obsValid<- obs[-part,]
        obsAll<-obs         
        obs<-obsTrain
        # Coordinates - Raster Transformation
        coordinates(obs)= ~x+y
        coordinates(obsValid)= ~x+y
        for (l in 1:nrofcombinations) { #(l in 1:nrofcombinations)
          #start_time <- Sys.time()
          print(l)
          print("Is the No. of the predictors combination")
          
          #_____________________________
          
          if (i==1){
            predictors = list.dirs(path = predictors_dir_1)
            predictors = predictors[grepl(pattern = resolution,x = predictors)]
            if (predictors_combination=="personalized"){
              forremoving<-my_different_combinations_of_predictors[[1,l]]
              assign(paste0("toremove"),agj [! agj %in% forremoving])
              for (opü in toremove) {
                predictors = predictors[!grepl(pattern=paste0("_",opü,"_"),x=predictors)]
              }
            } else if (predictors_combination=="automatic"){
              assign(paste0("combi_",l),as.numeric(resultelvis[[l]]))
              forremoving<-get(paste0("combi_",l))
              assign(paste0("toremove"),agj [! agj %in% forremoving])
              for (opü in toremove) {
                predictors = predictors[!grepl(pattern=paste0("_",opü,"_"),x=predictors)]
              }
            } else {}
            
          } else if (i==2) {
            predictors = list.dirs(path = predictors_dir_2)
            predictors = predictors[grepl(pattern = resolution,x = predictors)]
            if (predictors_combination=="personalized"){
              forremoving<-my_different_combinations_of_predictors[[1,l]]
              assign(paste0("toremove"),agj [! agj %in% forremoving])
              for (opü in toremove) {
                predictors = predictors[!grepl(pattern=paste0("_",opü,"_"),x=predictors)]
              }
            } else if (predictors_combination=="automatic") {
              assign(paste0("combi_",l),as.numeric(resultelvis[[l]]))
              forremoving<-get(paste0("combi_",l))
              assign(paste0("toremove"),agj [! agj %in% forremoving])
              for (opü in toremove) {
                predictors = predictors[!grepl(pattern=paste0("_",opü,"_"),x=predictors)]
              }
            } else {}
          } else {}
          
          # Cropped Rasters
          cropped_rasters = lapply(predictors, function(x) resample(raster(x),landuse, method="ngb"))
          r = raster::stack(cropped_rasters)
          r <- raster::mask(r,landuse)
          df=as.data.frame(r)
          df=df[complete.cases(df), ]
          gridCoords<-xyFromCell(object = r,which(complete.cases(as.data.frame(r))))
          grid=SpatialPixelsDataFrame(points = gridCoords,proj4string = CRS(proj4string(r)),data = df)
          proj4string(obs)= proj4string(grid)
          proj4string(obsValid)= proj4string(grid)
          obs@data<-cbind(obs@data,raster::extract(r, obs))
          
          print("Enter in RandomForest-Formula")
          if (ai[k]=="raw") {
            # RANDOM FOREST
            # Formula
            formula <- as.formula(paste(var," ~ ",  paste(names(r), collapse = "+")))
            formula_variabile = as.formula(paste(var," ~ 1"))
            formel<-formula
          } else if (ai[k]=="log") {
            # RANDOM FOREST
            # Formula
            obs@data$log = log(obs@data[[variables[j]]])
            formula <- as.formula(paste0("log"," ~ ",  paste(names(r), collapse = "+")))
            formula_variabile = as.formula(paste0("log"," ~ 1"))
            formel<-formula
          } else if (ai[k]=="sqrt") {
            # RANDOM FOREST
            # Formula
            obs@data$sqrt = sqrt(obs@data[[variables[j]]])
            formula <- as.formula(paste0("sqrt"," ~ ",  paste(names(r), collapse = "+")))
            formula_variabile = as.formula(paste0("sqrt"," ~ 1"))
            formel<-formula
          } else if (ai[k]=="quad") {
            # RANDOM FOREST
            # Formula
            obs@data$quad = quad(obs@data[[variables[j]]])
            formula <- as.formula(paste0("quad"," ~ ",  paste(names(r), collapse = "+")))
            formula_variabile = as.formula(paste0("quad"," ~ 1"))
            formel<-formula
          } else if (ai[k]=="inverse") {
            # RANDOM FOREST
            # Formula
            obs@data$inverse = inverse(obs@data[[variables[j]]])
            formula <- as.formula(paste0("inverse"," ~ ",  paste(names(r), collapse = "+")))
            formula_variabile = as.formula(paste0("inverse"," ~ 1"))
            formel<-formula
          } else {}
          
          #print("Enter in SearchCutoff")
          # SEARCH CUTOFF FILE GIUSTO
          cutoffs<-vector()
          searchcutoffsteps<-seq(from = 1000, to = 50000, by = 1000)
          for (lk in searchcutoffsteps) {
            try(modelRf<-fit.gstatModel(observations = obs,formulaString = formel,covariates = grid,method="randomForest",cutoff=lk))
            #modelRf<-fit.gstatModel(observations = obs,formulaString = formel,covariates = grid,method="randomForest",cutoff=lk)
            rangee<-modelRf@vgmModel$range[2]
            if (rangee<lk){
              sdf<-length(cutoffs)
              cutoffs<-append(cutoffs,lk,after=sdf)
              print(paste0(lk," is A cutoff"))
            } else if (rangee>=lk) {
              print(paste0(lk," is NO cutoff"))
            } else {}
          }
          #print(paste0("Fine_Search_Cutoff_",var,"_",ai[k],"_s=",s,"nrofcombi",l))
          
          #_____________________________
          
          if (length(cutoffs)==0) {
            cutoffs<-vector()
            searchcutoffsteps<-seq(from = 500, to = 50000, by = 500)
            grenz<-0
            
            while (length(cutoffs)==0 & grenz<10) {
              
              # Read Masterfile - Soil Data Points
              obs<- read.csv(file = file.path(datapoints_dir,list.files(path = datapoints_dir)),header = TRUE,sep = ",",dec = ".",na.strings = c("-9999"))  
              obs<- obs[!duplicated(obs[,c("x","y")]),]
              # Clean Masterfile from missing values
              obs<-obs[which(!is.na(obs[,var])),]
              # Divide data - Training 80% vs. Validation 20%
              part<-createDataPartition(y = obs[,var],times = 1,p = 0.8,list = F)
              obsTrain<- obs[part,]
              obsValid<- obs[-part,]
              obsAll<-obs         
              obs<-obsTrain
              # Coordinates - Raster Transformation
              coordinates(obs)= ~x+y
              coordinates(obsValid)= ~x+y
              
              if (i==1){
                predictors = list.dirs(path = predictors_dir_1)
                predictors = predictors[grepl(pattern = resolution,x = predictors)]
                if (predictors_combination=="personalized"){
                  forremoving<-my_different_combinations_of_predictors[[1,l]]
                  assign(paste0("toremove"),agj [! agj %in% forremoving])
                  for (opü in toremove) {
                    predictors = predictors[!grepl(pattern=paste0("_",opü,"_"),x=predictors)]
                  }
                } else if (predictors_combination=="automatic"){
                  assign(paste0("combi_",l),as.numeric(resultelvis[[l]]))
                  forremoving<-get(paste0("combi_",l))
                  assign(paste0("toremove"),agj [! agj %in% forremoving])
                  for (opü in toremove) {
                    predictors = predictors[!grepl(pattern=paste0("_",opü,"_"),x=predictors)]
                  }
                } else {}
                
              } else if (i==2) {
                predictors = list.dirs(path = predictors_dir_2)
                predictors = predictors[grepl(pattern = resolution,x = predictors)]
                if (predictors_combination=="personalized"){
                  forremoving<-my_different_combinations_of_predictors[[1,l]]
                  assign(paste0("toremove"),agj [! agj %in% forremoving])
                  for (opü in toremove) {
                    predictors = predictors[!grepl(pattern=paste0("_",opü,"_"),x=predictors)]
                  }
                } else if (predictors_combination=="automatic") {
                  assign(paste0("combi_",l),as.numeric(resultelvis[[l]]))
                  forremoving<-get(paste0("combi_",l))
                  assign(paste0("toremove"),agj [! agj %in% forremoving])
                  for (opü in toremove) {
                    predictors = predictors[!grepl(pattern=paste0("_",opü,"_"),x=predictors)]
                  }
                } else {}
              } else {}
              
              # Cropped Rasters
              cropped_rasters = lapply(predictors, function(x) resample(raster(x),landuse, method="ngb"))
              r = raster::stack(cropped_rasters)
              r <- raster::mask(r,landuse)
              df=as.data.frame(r)
              df=df[complete.cases(df), ]
              gridCoords<-xyFromCell(object = r,which(complete.cases(as.data.frame(r))))
              grid=SpatialPixelsDataFrame(points = gridCoords,proj4string = CRS(proj4string(r)),data = df)
              proj4string(obs)= proj4string(grid)
              proj4string(obsValid)= proj4string(grid)
              obs@data<-cbind(obs@data,raster::extract(r, obs))
              
              print("Enter in RandomForest-Formula")
              if (ai[k]=="raw") {
                # RANDOM FOREST
                # Formula
                formula <- as.formula(paste(var," ~ ",  paste(names(r), collapse = "+")))
                formula_variabile = as.formula(paste(var," ~ 1"))
                formel<-formula
              } else if (ai[k]=="log") {
                # RANDOM FOREST
                # Formula
                obs@data$log = log(obs@data[[variables[j]]])
                formula <- as.formula(paste0("log"," ~ ",  paste(names(r), collapse = "+")))
                formula_variabile = as.formula(paste0("log"," ~ 1"))
                formel<-formula
              } else if (ai[k]=="sqrt") {
                # RANDOM FOREST
                # Formula
                obs@data$sqrt = sqrt(obs@data[[variables[j]]])
                formula <- as.formula(paste0("sqrt"," ~ ",  paste(names(r), collapse = "+")))
                formula_variabile = as.formula(paste0("sqrt"," ~ 1"))
                formel<-formula
              } else if (ai[k]=="quad") {
                # RANDOM FOREST
                # Formula
                obs@data$quad = quad(obs@data[[variables[j]]])
                formula <- as.formula(paste0("quad"," ~ ",  paste(names(r), collapse = "+")))
                formula_variabile = as.formula(paste0("quad"," ~ 1"))
                formel<-formula
              } else if (ai[k]=="inverse") {
                # RANDOM FOREST
                # Formula
                obs@data$inverse = inverse(obs@data[[variables[j]]])
                formula <- as.formula(paste0("inverse"," ~ ",  paste(names(r), collapse = "+")))
                formula_variabile = as.formula(paste0("inverse"," ~ 1"))
                formel<-formula
              } else {}
              
              #print("Enter in SearchCutoff")
              # SEARCH CUTOFF FILE GIUSTO
              #searchcutoffsteps<-seq(from = 500, to = 50000, by = 500)
              for (lk in searchcutoffsteps) {
                try(modelRf<-fit.gstatModel(observations = obs,formulaString = formel,covariates = grid,method="randomForest",cutoff=lk))
                #modelRf<-fit.gstatModel(observations = obs,formulaString = formel,covariates = grid,method="randomForest",cutoff=lk)
                rangee<-modelRf@vgmModel$range[2]
                if (rangee<lk){
                  sdf<-length(cutoffs)
                  cutoffs<-append(cutoffs,lk,after=sdf)
                  print(paste0(lk," is A cutoff"))
                } else if (rangee>=lk) {
                  print(paste0(lk," is NO cutoff"))
                } else {}
              }
              #print(paste0("Fine_Search_Cutoff_",var,"_",ai[k],"_s=",s,"nrofcombi",l))
              
              grenz<-grenz+1
            }
            
          } else {}
          
          #___________________________________
          
          if (length(cutoffs)>=3) {
            print("length(cutoffs)>=3")
            # Output folders
            fg<-paste0("/output/",predictareas_name[i],"/",variables[j],"/",ai[k])
            #fh<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),fg)
            fh<-paste0(script_dir,fg)
            setwd(fh)
            dir.create(paste0("No.Run_",nrun))
            fgg<-paste0("/output/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun)
            #fhh<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),fgg)
            fhh<-paste0(script_dir,fgg)
            setwd(fhh)
            dir.create("Random_Forest_Regression")
            jkp<-paste0(out_dir,"/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun,"/Random_Forest_Regression")
            setwd(jkp)
            dir.create("1_Variogram")
            jkp1<-paste0(out_dir,"/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun,"/Random_Forest_Regression/1_Variogram")
            dir.create("2_RfModel")
            jkp2<-paste0(out_dir,"/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun,"/Random_Forest_Regression/2_RfModel")
            dir.create("3_VgmModel")
            jkp3<-paste0(out_dir,"/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun,"/Random_Forest_Regression/3_VgmModel")
            dir.create("4_PredictorsImportance")
            jkp4<-paste0(out_dir,"/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun,"/Random_Forest_Regression/4_PredictorsImportance")
            setwd(fhh)
            dir.create("Ordinary_Kriging_Prediction")
            jkpp<-paste0(out_dir,"/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun,"/Ordinary_Kriging_Prediction")
            dir.create("Final_Models")
            jkppp<-paste0(out_dir,"/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun,"/Final_Models")
            
            # Variogram
            #plot(variogram(formula_variabile, obs), plot.nu=T, pch="+")
            x=variogram(formula_variabile, obs)
            plott = ggplot(x,aes(dist,gamma))+
              geom_point()+
              theme_minimal()+
              ggtitle(paste0("Variogram - ","[",variables[j],"] ","- ",ai[k],"-","predictarea",i))+
              theme(plot.title = element_text(hjust = 0.5))
            jkpt<-paste0("Variogram - ","[",variables[j],"] ","- ",ai[k],"-","predictarea",i)
            #png(file.path(jkp1,paste0("1_Variogram1_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =4000, height=4000)
            #plot(variogram(formula_variabile, obs), plot.nu=T, pch="+",main=jkpt)
            #dev.off()
            ggsave(filename = file.path(jkp1,paste0("1_Variogram_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),plot = plott,device = "png",width=5, height=5)
            
            
            if (s==1) {
              print("length(cutoffs)>=3 / s==1")
              mincutoff<-min(cutoffs)
              try(modelRf<-fit.gstatModel(observations = obs,formulaString = formula,covariates = grid,method="randomForest",cutoff=mincutoff))
              taia<-mincutoff
              # MAIN CODE
              # Save RfModel & VgmModel
              #png(file.path(jkp2,paste0("2_RfModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".png")),res=300,width =8000, height=4000)
              #plot(modelRf)
              #dev.off()
              with(faithful,plot(modelRf)) 
              dev.copy(png,file.path(jkp2,paste0("2_RfModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".png")),res=300,width =3500, height=1750)
              dev.off()
              #dev.off(which = dev.cur())
              vgmModel<-modelRf@vgmModel
              #print(vgmModel)
              write.csv(vgmModel,file.path(jkp3,paste0("3_VgmModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".csv")))    
              write.xlsx(vgmModel,file.path(jkp3,paste0("3_VgmModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".xlsx")))    
              sink(file.path(jkp3,paste0("3_VgmModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".txt"))) 
              print(vgmModel)
              sink()
              
              # Save Importance of the predictors (%IncMse & IncNodePurity)
              # PredictorsImportance
              importance = modelRf@regModel$importance
              # %IncMse
              incmse<-sort(importance[,1],decreasing = TRUE)
              # Define names for predictors
              if (i==1) {
                names(incmse)[names(incmse)== "P_1_aspect_ea_150"] <- "P1__aspect east"
                names(incmse)[names(incmse)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(incmse)[names(incmse)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(incmse)[names(incmse)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(incmse)[names(incmse)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(incmse)[names(incmse)== "P_6_cl_grpre1_150"] <- "P6_land use - alpine meadow"
                names(incmse)[names(incmse)== "P_7_cl_grpre2_150"] <- "P7_land use - grassland"
                names(incmse)[names(incmse)== "P_8_cl_grpre3_150"] <- "P8_land use - grassland second shearing"
                names(incmse)[names(incmse)== "P_9_cl_grpre4_150"] <- "P9_land use - grassland special area"
                names(incmse)[names(incmse)== "P_10_dtm_150m"] <- "P10_dtm"
                names(incmse)[names(incmse)== "P_11_geo_bin_1_150"] <- "P11_geology - shist"
                names(incmse)[names(incmse)== "P_12_geo_bin_2_150"] <- "P12_geology - vulcanite"
                names(incmse)[names(incmse)== "P_13_geo_bin_3_150"] <- "P13_geology - vinschgau shear zone"
                names(incmse)[names(incmse)== "P_14_geo_bin_4_150"] <- "P14_geology - quaternary depositions"
                names(incmse)[names(incmse)== "P_15_geo_bin_5_150"] <- "P15_geology - sediment sequences"
                names(incmse)[names(incmse)== "P_16_geo_bin_6_150"] <- "P16_geology - plutons"
                names(incmse)[names(incmse)== "P_17_geo_bin_7_150"] <- "P17_geology - central gneiss"
                names(incmse)[names(incmse)== "P_18_geo_bin_8_150"] <- "P18_geology - quartz phyllite"
                names(incmse)[names(incmse)== "P_19_geo_bin_9_150"] <- "P19_geology - insignificant alpine metamorphosis"
                names(incmse)[names(incmse)== "P_20_geo_bi_10_150"] <- "P20_geology - medium-low alpine metamorphosis"
                names(incmse)[names(incmse)== "P_21_geo_bi_11_150"] <- "P21_geology - medium-high alpine metamorphosis"
                names(incmse)[names(incmse)== "P_22_slope_150m"] <- "P22_slope"
                names(incmse)[names(incmse)== "P_23_wetness_150m"] <- "P23_wetness_index"
                names(incmse)[names(incmse)== "P_24_gve_wgs84_150"] <- "P24_GVE"
                names(incmse)[names(incmse)== "P_25_geo_aeubin150"] <- "P25_Acidic soil"
                names(incmse)[names(incmse)== "P_26_geo_beubin150"] <- "P26_Basic soil"
                names(incmse)[names(incmse)== "P_27_geo_aelbin150"] <- "P27_Acidic soil"
                names(incmse)[names(incmse)== "P_28_geo_belbin150"] <- "P28_Basic soil"
                names(incmse)[names(incmse)== "P_29_geo_abelbi150"] <- "P29_Acidic-Basic soil"
              } else if (i==2) {
                names(incmse)[names(incmse)== "P_1_aspect_ea_150"] <- "P1_aspect east"
                names(incmse)[names(incmse)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(incmse)[names(incmse)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(incmse)[names(incmse)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(incmse)[names(incmse)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(incmse)[names(incmse)== "P_6_cl_grarp1_150"] <- "P6_land use - alpine meadow"
                names(incmse)[names(incmse)== "P_7_cl_grarp2_150"] <- "P7_land use - grassland"
                names(incmse)[names(incmse)== "P_8_cl_grarp3_150"] <- "P8_land use - grassland second shearing"
                names(incmse)[names(incmse)== "P_9_cl_grarp4_150"] <- "P9_land use - grassland special area"
                names(incmse)[names(incmse)== "P_10_cl_grarp5_150"] <- "P10_land use - crops"
                names(incmse)[names(incmse)== "P_11_dtm_150m"] <- "P11_dtm"
                names(incmse)[names(incmse)== "P_12_geo_bin_1_150"] <- "P12_geology - shist"
                names(incmse)[names(incmse)== "P_13_geo_bin_2_150"] <- "P13_geology - vulcanite"
                names(incmse)[names(incmse)== "P_14_geo_bin_3_150"] <- "P14_geology - vinschgau shear zone"
                names(incmse)[names(incmse)== "P_15_geo_bin_4_150"] <- "P15_geology - quaternary depositions"
                names(incmse)[names(incmse)== "P_16_geo_bin_5_150"] <- "P16_geology - sediment sequences"
                names(incmse)[names(incmse)== "P_17_geo_bin_6_150"] <- "P17_geology - plutons"
                names(incmse)[names(incmse)== "P_18_geo_bin_7_150"] <- "P18_geology - central gneiss"
                names(incmse)[names(incmse)== "P_19_geo_bin_8_150"] <- "P19_geology - quartz phyllite"
                names(incmse)[names(incmse)== "P_20_geo_bin_9_150"] <- "P20_geology - insignificant alpine metamorphosis"
                names(incmse)[names(incmse)== "P_21_geo_bi_10_150"] <- "P21_geology - medium-low alpine metamorphosis"
                names(incmse)[names(incmse)== "P_22_geo_bi_11_150"] <- "P22_geology - medium-high alpine metamorphosis"
                names(incmse)[names(incmse)== "P_23_slope_150m"] <- "P23_slope"
                names(incmse)[names(incmse)== "P_24_wetness_150m"] <- "P24_wetness_index"
                names(incmse)[names(incmse)== "P_25_gve_wgs84_150"] <- "P25_GVE"
                names(incmse)[names(incmse)== "P_26_geo_aeubin150"] <- "P26_Acidic soil"
                names(incmse)[names(incmse)== "P_27_geo_beubin150"] <- "P27_Basic soil"
                names(incmse)[names(incmse)== "P_28_geo_aelbin150"] <- "P28_Acidic soil"
                names(incmse)[names(incmse)== "P_29_geo_belbin150"] <- "P29_Basic soil"
                names(incmse)[names(incmse)== "P_30_geo_abelbi150"] <- "P30_Acidic-Basic soil"
              } else {}
              
              # Save PNG Table %IncMSE
              incmsedff <- data.frame(incmse)
              names(incmsedff) <- ("%IncMSE")
              table <- tableGrob(incmsedff)
              title <- textGrob(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia),gp=gpar(fontsize=8))
              footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
              padding <- unit(2,"line")
              table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
              table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
              table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
              png(file.path(jkp4,paste0("4_IncMse_Table_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")), units="cm", width=(ncol(table))*15, height=(nrow(table)), res=300) #, height = 200, width = 100)
              grid.newpage() 
              grid.draw(table) 
              dev.off()
              # Save VALUES Table %IncMSE
              write.csv(incmse,file.path(jkp4,paste0("4_IncMse_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".csv")))      
              write.xlsx(incmse,file.path(jkp4,paste0("4_IncMse_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".xlsx")))    
              # Save Plot %IncMSE
              x1<-incmsedff[,1]
              x<-rev(x1)
              y1<-names(incmse)
              y<-rev(y1)
              png(file.path(jkp4,paste0("4_IncMse_Plot_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              dotchart(x,y,xlab="%IncMse") 
              title(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun))
              dev.off()
              
              # IncNodePurity
              IncNodePurity<-round(sort(importance[,2],decreasing = TRUE),digits = 2) 
              # Define names for predictors
              if (i==1) {
                names(IncNodePurity)[names(IncNodePurity)== "P_1_aspect_ea_150"] <- "P1__aspect east"
                names(IncNodePurity)[names(IncNodePurity)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(IncNodePurity)[names(IncNodePurity)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(IncNodePurity)[names(IncNodePurity)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(IncNodePurity)[names(IncNodePurity)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(IncNodePurity)[names(IncNodePurity)== "P_6_cl_grpre1_150"] <- "P6_land use - alpine meadow"
                names(IncNodePurity)[names(IncNodePurity)== "P_7_cl_grpre2_150"] <- "P7_land use - grassland"
                names(IncNodePurity)[names(IncNodePurity)== "P_8_cl_grpre3_150"] <- "P8_land use - grassland second shearing"
                names(IncNodePurity)[names(IncNodePurity)== "P_9_cl_grpre4_150"] <- "P9_land use - grassland special area"
                names(IncNodePurity)[names(IncNodePurity)== "P_10_dtm_150m"] <- "P10_dtm"
                names(IncNodePurity)[names(IncNodePurity)== "P_11_geo_bin_1_150"] <- "P11_geology - shist"
                names(IncNodePurity)[names(IncNodePurity)== "P_12_geo_bin_2_150"] <- "P12_geology - vulcanite"
                names(IncNodePurity)[names(IncNodePurity)== "P_13_geo_bin_3_150"] <- "P13_geology - vinschgau shear zone"
                names(IncNodePurity)[names(IncNodePurity)== "P_14_geo_bin_4_150"] <- "P14_geology - quaternary depositions"
                names(IncNodePurity)[names(IncNodePurity)== "P_15_geo_bin_5_150"] <- "P15_geology - sediment sequences"
                names(IncNodePurity)[names(IncNodePurity)== "P_16_geo_bin_6_150"] <- "P16_geology - plutons"
                names(IncNodePurity)[names(IncNodePurity)== "P_17_geo_bin_7_150"] <- "P17_geology - central gneiss"
                names(IncNodePurity)[names(IncNodePurity)== "P_18_geo_bin_8_150"] <- "P18_geology - quartz phyllite"
                names(IncNodePurity)[names(IncNodePurity)== "P_19_geo_bin_9_150"] <- "P19_geology - insignificant alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_20_geo_bi_10_150"] <- "P20_geology - medium-low alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_21_geo_bi_11_150"] <- "P21_geology - medium-high alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_22_slope_150m"] <- "P22_slope"
                names(IncNodePurity)[names(IncNodePurity)== "P_23_wetness_150m"] <- "P23_wetness_index"
                names(IncNodePurity)[names(IncNodePurity)== "P_24_gve_wgs84_150"] <- "P24_GVE"
                names(IncNodePurity)[names(IncNodePurity)== "P_25_geo_aeubin150"] <- "P25_Acidic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_26_geo_beubin150"] <- "P26_Basic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_27_geo_aelbin150"] <- "P27_Acidic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_28_geo_belbin150"] <- "P28_Basic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_29_geo_abelbi150"] <- "P29_Acidic-Basic soil"
              } else if (i==2) {
                names(IncNodePurity)[names(IncNodePurity)== "P_1_aspect_ea_150"] <- "P1_aspect east"
                names(IncNodePurity)[names(IncNodePurity)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(IncNodePurity)[names(IncNodePurity)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(IncNodePurity)[names(IncNodePurity)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(IncNodePurity)[names(IncNodePurity)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(IncNodePurity)[names(IncNodePurity)== "P_6_cl_grarp1_150"] <- "P6_land use - alpine meadow"
                names(IncNodePurity)[names(IncNodePurity)== "P_7_cl_grarp2_150"] <- "P7_land use - grassland"
                names(IncNodePurity)[names(IncNodePurity)== "P_8_cl_grarp3_150"] <- "P8_land use - grassland second shearing"
                names(IncNodePurity)[names(IncNodePurity)== "P_9_cl_grarp4_150"] <- "P9land use - grassland special area"
                names(IncNodePurity)[names(IncNodePurity)== "P_10_cl_grarp5_150"] <- "P10_land use - crops"
                names(IncNodePurity)[names(IncNodePurity)== "P_11_dtm_150m"] <- "P11_dtm"
                names(IncNodePurity)[names(IncNodePurity)== "P_12_geo_bin_1_150"] <- "P12_geology - shist"
                names(IncNodePurity)[names(IncNodePurity)== "P_13_geo_bin_2_150"] <- "P13_geology - vulcanite"
                names(IncNodePurity)[names(IncNodePurity)== "P_14_geo_bin_3_150"] <- "P14_geology - vinschgau shear zone"
                names(IncNodePurity)[names(IncNodePurity)== "P_15_geo_bin_4_150"] <- "P15_geology - quaternary depositions"
                names(IncNodePurity)[names(IncNodePurity)== "P_16_geo_bin_5_150"] <- "P16_geology - sediment sequences"
                names(IncNodePurity)[names(IncNodePurity)== "P_17_geo_bin_6_150"] <- "P17_geology - plutons"
                names(IncNodePurity)[names(IncNodePurity)== "P_18_geo_bin_7_150"] <- "P18_geology - central gneiss"
                names(IncNodePurity)[names(IncNodePurity)== "P_19_geo_bin_8_150"] <- "P19_geology - quartz phyllite"
                names(IncNodePurity)[names(IncNodePurity)== "P_20_geo_bin_9_150"] <- "P20_geology - insignificant alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_21_geo_bi_10_150"] <- "P21_geology - medium-low alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_22_geo_bi_11_150"] <- "P22_geology - medium-high alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_23_slope_150m"] <- "P23_slope"
                names(IncNodePurity)[names(IncNodePurity)== "P_24_wetness_150m"] <- "P24_wetness_index"
                names(IncNodePurity)[names(IncNodePurity)== "P_25_gve_wgs84_150"] <- "P25_GVE"
                names(IncNodePurity)[names(IncNodePurity)== "P_26_geo_aeubin150"] <- "P26_Acidic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_27_geo_beubin150"] <- "P27_Basic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_28_geo_aelbin150"] <- "P28_Acidic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_29_geo_belbin150"] <- "P29_Basic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_30_geo_abelbi150"] <- "P30_Acidic-Basic soil"
              } else {}
              
              # Save PNG Table IncNodePurity
              incnodedff <- data.frame(IncNodePurity)
              names(incnodedff) <- ("IncNodePurity")
              table <- tableGrob(incnodedff)
              title <- textGrob(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia),gp=gpar(fontsize=8))
              footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
              padding <- unit(2,"line")
              table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
              table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
              table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
              png(file.path(jkp4,paste0("4_IncNodePurity_Table_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")), units="cm", width=(ncol(table))*15, height=(nrow(table)), res=300) #, height = 200, width = 100)
              grid.newpage() 
              grid.draw(table) 
              dev.off()
              # Save VALUES Table IncNodePurity
              write.csv(IncNodePurity,file.path(jkp4,paste0("4_IncNodePurity_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".csv")))      
              write.xlsx(IncNodePurity,file.path(jkp4,paste0("4_IncNodePurity_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".xlsx"))) 
              # Save Plot IncNodePurity
              x1<-incnodedff[,1]
              x<-rev(x1)
              y1<-names(IncNodePurity)
              y<-rev(y1)
              png(file.path(jkp4,paste0("4_IncNodePurity_Plot_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              dotchart(x,y,xlab="IncNodePurity") 
              title(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun))
              dev.off()
              
              # ImportanceSD (standard deviation)
              importancesdd<-modelRf@regModel$importanceSD
              importancesd<-sort(importancesdd,decreasing=TRUE)
              if (i==1) {
                names(importancesd)[names(importancesd)== "P_1_aspect_ea_150"] <- "P1__aspect east"
                names(importancesd)[names(importancesd)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(importancesd)[names(importancesd)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(importancesd)[names(importancesd)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(importancesd)[names(importancesd)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(importancesd)[names(importancesd)== "P_6_cl_grpre1_150"] <- "P6_land use - alpine meadow"
                names(importancesd)[names(importancesd)== "P_7_cl_grpre2_150"] <- "P7_land use - grassland"
                names(importancesd)[names(importancesd)== "P_8_cl_grpre3_150"] <- "P8_land use - grassland second shearing"
                names(importancesd)[names(importancesd)== "P_9_cl_grpre4_150"] <- "P9_land use - grassland special area"
                names(importancesd)[names(importancesd)== "P_10_dtm_150m"] <- "P10_dtm"
                names(importancesd)[names(importancesd)== "P_11_geo_bin_1_150"] <- "P11_geology - shist"
                names(importancesd)[names(importancesd)== "P_12_geo_bin_2_150"] <- "P12_geology - vulcanite"
                names(importancesd)[names(importancesd)== "P_13_geo_bin_3_150"] <- "P13_geology - vinschgau shear zone"
                names(importancesd)[names(importancesd)== "P_14_geo_bin_4_150"] <- "P14_geology - quaternary depositions"
                names(importancesd)[names(importancesd)== "P_15_geo_bin_5_150"] <- "P15_geology - sediment sequences"
                names(importancesd)[names(importancesd)== "P_16_geo_bin_6_150"] <- "P16_geology - plutons"
                names(importancesd)[names(importancesd)== "P_17_geo_bin_7_150"] <- "P17_geology - central gneiss"
                names(importancesd)[names(importancesd)== "P_18_geo_bin_8_150"] <- "P18_geology - quartz phyllite"
                names(importancesd)[names(importancesd)== "P_19_geo_bin_9_150"] <- "P19_geology - insignificant alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_20_geo_bi_10_150"] <- "P20_geology - medium-low alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_21_geo_bi_11_150"] <- "P21_geology - medium-high alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_22_slope_150m"] <- "P22_slope"
                names(importancesd)[names(importancesd)== "P_23_wetness_150m"] <- "P23_wetness_index"
                names(importancesd)[names(importancesd)== "P_24_gve_wgs84_150"] <- "P24_GVE"
                names(importancesd)[names(importancesd)== "P_25_geo_aeubin150"] <- "P25_Acidic soil"
                names(importancesd)[names(importancesd)== "P_26_geo_beubin150"] <- "P26_Basic soil"
                names(importancesd)[names(importancesd)== "P_27_geo_aelbin150"] <- "P27_Acidic soil"
                names(importancesd)[names(importancesd)== "P_28_geo_belbin150"] <- "P28_Basic soil"
                names(importancesd)[names(importancesd)== "P_29_geo_abelbi150"] <- "P29_Acidic-Basic soil"
              } else if (i==2) {
                names(importancesd)[names(importancesd)== "P_1_aspect_ea_150"] <- "P1_aspect east"
                names(importancesd)[names(importancesd)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(importancesd)[names(importancesd)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(importancesd)[names(importancesd)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(importancesd)[names(importancesd)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(importancesd)[names(importancesd)== "P_6_cl_grarp1_150"] <- "P6_land use - alpine meadow"
                names(importancesd)[names(importancesd)== "P_7_cl_grarp2_150"] <- "P7_land use - grassland"
                names(importancesd)[names(importancesd)== "P_8_cl_grarp3_150"] <- "P8_land use - grassland second shearing"
                names(importancesd)[names(importancesd)== "P_9_cl_grarp4_150"] <- "P9land use - grassland special area"
                names(importancesd)[names(importancesd)== "P_10_cl_grarp5_150"] <- "P10_land use - crops"
                names(importancesd)[names(importancesd)== "P_11_dtm_150m"] <- "P11_dtm"
                names(importancesd)[names(importancesd)== "P_12_geo_bin_1_150"] <- "P12_geology - shist"
                names(importancesd)[names(importancesd)== "P_13_geo_bin_2_150"] <- "P13_geology - vulcanite"
                names(importancesd)[names(importancesd)== "P_14_geo_bin_3_150"] <- "P14_geology - vinschgau shear zone"
                names(importancesd)[names(importancesd)== "P_15_geo_bin_4_150"] <- "P15_geology - quaternary depositions"
                names(importancesd)[names(importancesd)== "P_16_geo_bin_5_150"] <- "P16_geology - sediment sequences"
                names(importancesd)[names(importancesd)== "P_17_geo_bin_6_150"] <- "P17_geology - plutons"
                names(importancesd)[names(importancesd)== "P_18_geo_bin_7_150"] <- "P18_geology - central gneiss"
                names(importancesd)[names(importancesd)== "P_19_geo_bin_8_150"] <- "P19_geology - quartz phyllite"
                names(importancesd)[names(importancesd)== "P_20_geo_bin_9_150"] <- "P20_geology - insignificant alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_21_geo_bi_10_150"] <- "P21_geology - medium-low alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_22_geo_bi_11_150"] <- "P22_geology - medium-high alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_23_slope_150m"] <- "P23_slope"
                names(importancesd)[names(importancesd)== "P_24_wetness_150m"] <- "P24_wetness_index"
                names(importancesd)[names(importancesd)== "P_25_gve_wgs84_150"] <- "P25_GVE"
                names(importancesd)[names(importancesd)== "P_26_geo_aeubin150"] <- "P26_Acidic soil"
                names(importancesd)[names(importancesd)== "P_27_geo_beubin150"] <- "P27_Basic soil"
                names(importancesd)[names(importancesd)== "P_28_geo_aelbin150"] <- "P28_Acidic soil"
                names(importancesd)[names(importancesd)== "P_29_geo_belbin150"] <- "P29_Basic soil"
                names(importancesd)[names(importancesd)== "P_30_geo_abelbi150"] <- "P30_Acidic-Basic soil"
              } else {}
              
              # Save PNG Table importancesd
              importancesddff <- data.frame(importancesd)
              names(importancesddff) <- ("ImportanceSD")
              table <- tableGrob(importancesddff)
              title <- textGrob(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia),gp=gpar(fontsize=8))
              footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
              padding <- unit(2,"line")
              table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
              table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
              table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
              png(file.path(jkp4,paste0("4_ImportanceSD_Table_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")), units="cm", width=(ncol(table))*15, height=(nrow(table)), res=300) #, height = 200, width = 100)
              grid.newpage() 
              grid.draw(table) 
              dev.off()
              # Save VALUES Table importancesd
              write.csv(importancesd,file.path(jkp4,paste0("4_ImportanceSD_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".csv")))      
              write.xlsx(importancesd,file.path(jkp4,paste0("4_ImportanceSD_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".xlsx"))) 
              # Save Plot importancesd
              x1<-importancesddff[,1]
              x<-rev(x1)
              y1<-names(importancesd)
              y<-rev(y1)
              png(file.path(jkp4,paste0("4_ImportanceSD_Plot_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              dotchart(x,y,xlab="ImportanceSD") 
              title(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun))
              dev.off()
              
              # Mean of square residuals / % Var. explained
              printregmodel<-modelRf@regModel
              sink(file.path(jkp2,paste0("2_RfModel_Values",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".txt")))
              print(printregmodel)
              sink()
              setwd(paste0(out_dir,"/",predictareas_name[i],"/",var,"/",ai[k],"/","No.Run_",nrun,"/Random_Forest_Regression/2_RfModel"))
              armin<-read.delim(paste0("2_RfModel_Values",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".txt"),stringsAsFactors = FALSE)
              meanofsquareresiduals<-readr::parse_number(armin[5,])
              varexplained<-readr::parse_number(armin[6,])
              
              
              
              # Predictions ORDINARY KRIGING (with 80% of the data)
              predictionsRf<-GSIF::predict.gstatModel(object = modelRf, predictionLocations = grid, predict.method = "RK",verbose=T,nfold = 5)#vgmmodel = model@vgmModel,
              predictRf<-raster(predictionsRf@predicted)
              if (ai[k]=="raw") {
                
              } else if (ai[k]=="log") {
                predictRf_transf_log<-exp(predictRf)
              } else if (ai[k]=="sqrt") {
                predictRf_transf_sqrt<-quad(predictRf)
              } else if (ai[k]=="quad") {
                predictRf_transf_quad<-sqrt(predictRf)
              } else if (ai[k]=="inverse") {
                predictRf_transf_inverse<-inverse(predictRf)
              } else {}
              # Save Plots
              png(file.path(jkpp,paste0("5_Prediction_OrdinaryKriging_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              plot(predictionsRf) #paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              dev.off()
              
              
              
              # Plot Raster of predicted values
              png(file.path(jkppp,paste0("6_Predict_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              if (ai[k]=="raw") {
                raster::plot(predictRf,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else if (ai[k]=="log") {
                raster::plot(predictRf_transf_log,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else if (ai[k]=="sqrt") {
                raster::plot(predictRf_transf_sqrt,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else if (ai[k]=="quad") {
                raster::plot(predictRf_transf_quad,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else if (ai[k]=="inverse") {
                raster::plot(predictRf_transf_inverse,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else {}
              dev.off()
              
              
              
              
              # Overlay validation and prediction
              obsValidRf<-sp::over(obsValid,predictionsRf@predicted)
              if (ai[k]=="raw") {
                
              } else if (ai[k]=="log") {
                obsValidRf$original<-exp(obsValidRf$log)
              } else if (ai[k]=="sqrt") {
                obsValidRf$original<-quad(obsValidRf$sqrt)
              } else if (ai[k]=="quad") {
                obsValidRf$original<-sqrt(obsValidRf$quad)
              } else if (ai[k]=="inverse") {
                obsValidRf$original<-inverse(obsValidRf$inverse)
              } else {}
              
              
              
              # RMSE & NRMSE & R-SQUARED
              if (ai[k]=="raw") {
                rmseRf <- RMSE(pred = obsValidRf[,var], obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf[,var], obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf[,var], obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf[,var], obs = obsValid@data[[var]], na.rm = T)
              } else if (ai[k]=="log") {
                rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              } else if (ai[k]=="sqrt") {
                rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              } else if (ai[k]=="quad") {
                rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              } else if (ai[k]=="inverse") {
                rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              } else {}
              
              
              
              # Maps - Final rasters (GTiff)
              if (ai[k]=="raw") {
                writeRaster(x = predictRf, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else if (ai[k]=="log") {
                writeRaster(x = predictRf_transf_log, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else if (ai[k]=="sqrt") {
                writeRaster(x = predictRf_transf_sqrt, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else if (ai[k]=="quad") {
                writeRaster(x = predictRf_transf_quad, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else if (ai[k]=="inverse") {
                writeRaster(x = predictRf_transf_inverse,filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else {}
              
              # spalte 1 - No.Run
              sdf1<-length(spalte1)
              spalte1<-append(spalte1,nrun,after=sdf1)
              nrun<-nrun+1
              # Spalte 2 - Transformation
              sdf2<-length(spalte2)
              spalte2<-append(spalte2,ai[k],after=sdf2)
              # Spalte 3 - Cutoff
              sdf3<-length(spalte3)
              spalte3<-append(spalte3,taia,after=sdf3)
              # Spalte 4 - Predictors
              sdf4<-length(spalte4)
              spalte4<-append(spalte4,paste(get(paste0("predictors_combi_",l)), collapse = "/"))
              # Spalte 5 - Mean of squared residuals
              sdf5<-length(spalte5)
              spalte5<-append(spalte5,meanofsquareresiduals,after=sdf5)
              # Spalte 6 - % Var explained
              sdf6<-length(spalte6)
              spalte6<-append(spalte6,varexplained,after=sdf6)
              # Spalte 7 - RMSE
              sdf7<-length(spalte7)
              spalte7<-append(spalte7,rmseRf,after=sdf7)
              sdf7tot<-length(spalte7tot)
              spalte7tot<-append(spalte7tot,rmseRf,after=sdf7tot)
              # Spalte 8 - NRMSE
              sdf8<-length(spalte8)
              spalte8<-append(spalte8,nrmseRf,after=sdf8)
              # Spalte 9 - NRMSEMAXMIN
              sdf9<-length(spalte9)
              spalte9<-append(spalte9,nrmseMaxMinRf,after=sdf9)
              # Spalte 10 - R-SQUARED
              sdf10<-length(spalte10)
              spalte10<-append(spalte10,r2Rf,after=sdf10)
              sdf10tot<-length(spalte10tot)
              spalte10tot<-append(spalte10tot,r2Rf,after=sdf10tot)
              
            } else if (s==2) {
              print("length(cutoffs)>=3 / s==2")
              meancutoff<-cutoffs[length(cutoffs)/2]
              try(modelRf<-fit.gstatModel(observations = obs,formulaString = formula,covariates = grid,method="randomForest",cutoff=meancutoff))
              taia<-meancutoff
              # MAIN CODE
              # Save RfModel & VgmModel
              with(faithful,plot(modelRf)) 
              dev.copy(png,file.path(jkp2,paste0("2_RfModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".png")),res=300,width =3500, height=1750)
              dev.off()
              vgmModel = modelRf@vgmModel
              print(vgmModel)
              write.csv(vgmModel,file.path(jkp3,paste0("3_VgmModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".csv")))    
              write.xlsx(vgmModel,file.path(jkp3,paste0("3_VgmModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".xlsx")))    
              sink(file.path(jkp3,paste0("3_VgmModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".txt"))) 
              print(vgmModel)
              sink()
              
              # Save Importance of the predictors (%IncMse & IncNodePurity)
              # PredictorsImportance
              importance = modelRf@regModel$importance
              # %IncMse
              incmse<-sort(importance[,1],decreasing = TRUE)
              # Define names for predictors
              if (i==1) {
                names(incmse)[names(incmse)== "P_1_aspect_ea_150"] <- "P1__aspect east"
                names(incmse)[names(incmse)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(incmse)[names(incmse)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(incmse)[names(incmse)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(incmse)[names(incmse)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(incmse)[names(incmse)== "P_6_cl_grpre1_150"] <- "P6_land use - alpine meadow"
                names(incmse)[names(incmse)== "P_7_cl_grpre2_150"] <- "P7_land use - grassland"
                names(incmse)[names(incmse)== "P_8_cl_grpre3_150"] <- "P8_land use - grassland second shearing"
                names(incmse)[names(incmse)== "P_9_cl_grpre4_150"] <- "P9_land use - grassland special area"
                names(incmse)[names(incmse)== "P_10_dtm_150m"] <- "P10_dtm"
                names(incmse)[names(incmse)== "P_11_geo_bin_1_150"] <- "P11_geology - shist"
                names(incmse)[names(incmse)== "P_12_geo_bin_2_150"] <- "P12_geology - vulcanite"
                names(incmse)[names(incmse)== "P_13_geo_bin_3_150"] <- "P13_geology - vinschgau shear zone"
                names(incmse)[names(incmse)== "P_14_geo_bin_4_150"] <- "P14_geology - quaternary depositions"
                names(incmse)[names(incmse)== "P_15_geo_bin_5_150"] <- "P15_geology - sediment sequences"
                names(incmse)[names(incmse)== "P_16_geo_bin_6_150"] <- "P16_geology - plutons"
                names(incmse)[names(incmse)== "P_17_geo_bin_7_150"] <- "P17_geology - central gneiss"
                names(incmse)[names(incmse)== "P_18_geo_bin_8_150"] <- "P18_geology - quartz phyllite"
                names(incmse)[names(incmse)== "P_19_geo_bin_9_150"] <- "P19_geology - insignificant alpine metamorphosis"
                names(incmse)[names(incmse)== "P_20_geo_bi_10_150"] <- "P20_geology - medium-low alpine metamorphosis"
                names(incmse)[names(incmse)== "P_21_geo_bi_11_150"] <- "P21_geology - medium-high alpine metamorphosis"
                names(incmse)[names(incmse)== "P_22_slope_150m"] <- "P22_slope"
                names(incmse)[names(incmse)== "P_23_wetness_150m"] <- "P23_wetness_index"
                names(incmse)[names(incmse)== "P_24_gve_wgs84_150"] <- "P24_GVE"
                names(incmse)[names(incmse)== "P_25_geo_aeubin150"] <- "P25_Acidic soil"
                names(incmse)[names(incmse)== "P_26_geo_beubin150"] <- "P26_Basic soil"
                names(incmse)[names(incmse)== "P_27_geo_aelbin150"] <- "P27_Acidic soil"
                names(incmse)[names(incmse)== "P_28_geo_belbin150"] <- "P28_Basic soil"
                names(incmse)[names(incmse)== "P_29_geo_abelbi150"] <- "P29_Acidic-Basic soil"
              } else if (i==2) {
                names(incmse)[names(incmse)== "P_1_aspect_ea_150"] <- "P1_aspect east"
                names(incmse)[names(incmse)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(incmse)[names(incmse)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(incmse)[names(incmse)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(incmse)[names(incmse)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(incmse)[names(incmse)== "P_6_cl_grarp1_150"] <- "P6_land use - alpine meadow"
                names(incmse)[names(incmse)== "P_7_cl_grarp2_150"] <- "P7_land use - grassland"
                names(incmse)[names(incmse)== "P_8_cl_grarp3_150"] <- "P8_land use - grassland second shearing"
                names(incmse)[names(incmse)== "P_9_cl_grarp4_150"] <- "P9_land use - grassland special area"
                names(incmse)[names(incmse)== "P_10_cl_grarp5_150"] <- "P10_land use - crops"
                names(incmse)[names(incmse)== "P_11_dtm_150m"] <- "P11_dtm"
                names(incmse)[names(incmse)== "P_12_geo_bin_1_150"] <- "P12_geology - shist"
                names(incmse)[names(incmse)== "P_13_geo_bin_2_150"] <- "P13_geology - vulcanite"
                names(incmse)[names(incmse)== "P_14_geo_bin_3_150"] <- "P14_geology - vinschgau shear zone"
                names(incmse)[names(incmse)== "P_15_geo_bin_4_150"] <- "P15_geology - quaternary depositions"
                names(incmse)[names(incmse)== "P_16_geo_bin_5_150"] <- "P16_geology - sediment sequences"
                names(incmse)[names(incmse)== "P_17_geo_bin_6_150"] <- "P17_geology - plutons"
                names(incmse)[names(incmse)== "P_18_geo_bin_7_150"] <- "P18_geology - central gneiss"
                names(incmse)[names(incmse)== "P_19_geo_bin_8_150"] <- "P19_geology - quartz phyllite"
                names(incmse)[names(incmse)== "P_20_geo_bin_9_150"] <- "P20_geology - insignificant alpine metamorphosis"
                names(incmse)[names(incmse)== "P_21_geo_bi_10_150"] <- "P21_geology - medium-low alpine metamorphosis"
                names(incmse)[names(incmse)== "P_22_geo_bi_11_150"] <- "P22_geology - medium-high alpine metamorphosis"
                names(incmse)[names(incmse)== "P_23_slope_150m"] <- "P23_slope"
                names(incmse)[names(incmse)== "P_24_wetness_150m"] <- "P24_wetness_index"
                names(incmse)[names(incmse)== "P_25_gve_wgs84_150"] <- "P25_GVE"
                names(incmse)[names(incmse)== "P_26_geo_aeubin150"] <- "P26_Acidic soil"
                names(incmse)[names(incmse)== "P_27_geo_beubin150"] <- "P27_Basic soil"
                names(incmse)[names(incmse)== "P_28_geo_aelbin150"] <- "P28_Acidic soil"
                names(incmse)[names(incmse)== "P_29_geo_belbin150"] <- "P29_Basic soil"
                names(incmse)[names(incmse)== "P_30_geo_abelbi150"] <- "P30_Acidic-Basic soil"
              } else {}
              
              # Save PNG Table %IncMSE
              incmsedff <- data.frame(incmse)
              names(incmsedff) <- ("%IncMSE")
              table <- tableGrob(incmsedff)
              title <- textGrob(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia),gp=gpar(fontsize=8))
              footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
              padding <- unit(2,"line")
              table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
              table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
              table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
              png(file.path(jkp4,paste0("4_IncMse_Table_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")), units="cm", width=(ncol(table))*15, height=(nrow(table)), res=300) #, height = 200, width = 100)
              grid.newpage() 
              grid.draw(table) 
              dev.off()
              # Save VALUES Table %IncMSE
              write.csv(incmse,file.path(jkp4,paste0("4_IncMse_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".csv")))      
              write.xlsx(incmse,file.path(jkp4,paste0("4_IncMse_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".xlsx")))    
              # Save Plot %IncMSE
              x1<-incmsedff[,1]
              x<-rev(x1)
              y1<-names(incmse)
              y<-rev(y1)
              png(file.path(jkp4,paste0("4_IncMse_Plot_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              dotchart(x,y,xlab="%IncMse") 
              title(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun))
              dev.off()
              
              # IncNodePurity
              IncNodePurity<-round(sort(importance[,2],decreasing = TRUE),digits = 2) 
              # Define names for predictors
              if (i==1) {
                names(IncNodePurity)[names(IncNodePurity)== "P_1_aspect_ea_150"] <- "P1__aspect east"
                names(IncNodePurity)[names(IncNodePurity)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(IncNodePurity)[names(IncNodePurity)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(IncNodePurity)[names(IncNodePurity)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(IncNodePurity)[names(IncNodePurity)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(IncNodePurity)[names(IncNodePurity)== "P_6_cl_grpre1_150"] <- "P6_land use - alpine meadow"
                names(IncNodePurity)[names(IncNodePurity)== "P_7_cl_grpre2_150"] <- "P7_land use - grassland"
                names(IncNodePurity)[names(IncNodePurity)== "P_8_cl_grpre3_150"] <- "P8_land use - grassland second shearing"
                names(IncNodePurity)[names(IncNodePurity)== "P_9_cl_grpre4_150"] <- "P9_land use - grassland special area"
                names(IncNodePurity)[names(IncNodePurity)== "P_10_dtm_150m"] <- "P10_dtm"
                names(IncNodePurity)[names(IncNodePurity)== "P_11_geo_bin_1_150"] <- "P11_geology - shist"
                names(IncNodePurity)[names(IncNodePurity)== "P_12_geo_bin_2_150"] <- "P12_geology - vulcanite"
                names(IncNodePurity)[names(IncNodePurity)== "P_13_geo_bin_3_150"] <- "P13_geology - vinschgau shear zone"
                names(IncNodePurity)[names(IncNodePurity)== "P_14_geo_bin_4_150"] <- "P14_geology - quaternary depositions"
                names(IncNodePurity)[names(IncNodePurity)== "P_15_geo_bin_5_150"] <- "P15_geology - sediment sequences"
                names(IncNodePurity)[names(IncNodePurity)== "P_16_geo_bin_6_150"] <- "P16_geology - plutons"
                names(IncNodePurity)[names(IncNodePurity)== "P_17_geo_bin_7_150"] <- "P17_geology - central gneiss"
                names(IncNodePurity)[names(IncNodePurity)== "P_18_geo_bin_8_150"] <- "P18_geology - quartz phyllite"
                names(IncNodePurity)[names(IncNodePurity)== "P_19_geo_bin_9_150"] <- "P19_geology - insignificant alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_20_geo_bi_10_150"] <- "P20_geology - medium-low alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_21_geo_bi_11_150"] <- "P21_geology - medium-high alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_22_slope_150m"] <- "P22_slope"
                names(IncNodePurity)[names(IncNodePurity)== "P_23_wetness_150m"] <- "P23_wetness_index"
                names(IncNodePurity)[names(IncNodePurity)== "P_24_gve_wgs84_150"] <- "P24_GVE"
                names(IncNodePurity)[names(IncNodePurity)== "P_25_geo_aeubin150"] <- "P25_Acidic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_26_geo_beubin150"] <- "P26_Basic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_27_geo_aelbin150"] <- "P27_Acidic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_28_geo_belbin150"] <- "P28_Basic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_29_geo_abelbi150"] <- "P29_Acidic-Basic soil"
              } else if (i==2) {
                names(IncNodePurity)[names(IncNodePurity)== "P_1_aspect_ea_150"] <- "P1_aspect east"
                names(IncNodePurity)[names(IncNodePurity)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(IncNodePurity)[names(IncNodePurity)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(IncNodePurity)[names(IncNodePurity)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(IncNodePurity)[names(IncNodePurity)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(IncNodePurity)[names(IncNodePurity)== "P_6_cl_grarp1_150"] <- "P6_land use - alpine meadow"
                names(IncNodePurity)[names(IncNodePurity)== "P_7_cl_grarp2_150"] <- "P7_land use - grassland"
                names(IncNodePurity)[names(IncNodePurity)== "P_8_cl_grarp3_150"] <- "P8_land use - grassland second shearing"
                names(IncNodePurity)[names(IncNodePurity)== "P_9_cl_grarp4_150"] <- "P9land use - grassland special area"
                names(IncNodePurity)[names(IncNodePurity)== "P_10_cl_grarp5_150"] <- "P10_land use - crops"
                names(IncNodePurity)[names(IncNodePurity)== "P_11_dtm_150m"] <- "P11_dtm"
                names(IncNodePurity)[names(IncNodePurity)== "P_12_geo_bin_1_150"] <- "P12_geology - shist"
                names(IncNodePurity)[names(IncNodePurity)== "P_13_geo_bin_2_150"] <- "P13_geology - vulcanite"
                names(IncNodePurity)[names(IncNodePurity)== "P_14_geo_bin_3_150"] <- "P14_geology - vinschgau shear zone"
                names(IncNodePurity)[names(IncNodePurity)== "P_15_geo_bin_4_150"] <- "P15_geology - quaternary depositions"
                names(IncNodePurity)[names(IncNodePurity)== "P_16_geo_bin_5_150"] <- "P16_geology - sediment sequences"
                names(IncNodePurity)[names(IncNodePurity)== "P_17_geo_bin_6_150"] <- "P17_geology - plutons"
                names(IncNodePurity)[names(IncNodePurity)== "P_18_geo_bin_7_150"] <- "P18_geology - central gneiss"
                names(IncNodePurity)[names(IncNodePurity)== "P_19_geo_bin_8_150"] <- "P19_geology - quartz phyllite"
                names(IncNodePurity)[names(IncNodePurity)== "P_20_geo_bin_9_150"] <- "P20_geology - insignificant alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_21_geo_bi_10_150"] <- "P21_geology - medium-low alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_22_geo_bi_11_150"] <- "P22_geology - medium-high alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_23_slope_150m"] <- "P23_slope"
                names(IncNodePurity)[names(IncNodePurity)== "P_24_wetness_150m"] <- "P24_wetness_index"
                names(IncNodePurity)[names(IncNodePurity)== "P_25_gve_wgs84_150"] <- "P25_GVE"
                names(IncNodePurity)[names(IncNodePurity)== "P_26_geo_aeubin150"] <- "P26_Acidic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_27_geo_beubin150"] <- "P27_Basic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_28_geo_aelbin150"] <- "P28_Acidic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_29_geo_belbin150"] <- "P29_Basic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_30_geo_abelbi150"] <- "P30_Acidic-Basic soil"
              } else {}
              
              # Save PNG Table IncNodePurity
              incnodedff <- data.frame(IncNodePurity)
              names(incnodedff) <- ("IncNodePurity")
              table <- tableGrob(incnodedff)
              title <- textGrob(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia),gp=gpar(fontsize=8))
              footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
              padding <- unit(2,"line")
              table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
              table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
              table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
              png(file.path(jkp4,paste0("4_IncNodePurity_Table_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")), units="cm", width=(ncol(table))*15, height=(nrow(table)), res=300) #, height = 200, width = 100)
              grid.newpage() 
              grid.draw(table) 
              dev.off()
              # Save VALUES Table IncNodePurity
              write.csv(IncNodePurity,file.path(jkp4,paste0("4_IncNodePurity_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".csv")))      
              write.xlsx(IncNodePurity,file.path(jkp4,paste0("4_IncNodePurity_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".xlsx"))) 
              # Save Plot IncNodePurity
              x1<-incnodedff[,1]
              x<-rev(x1)
              y1<-names(IncNodePurity)
              y<-rev(y1)
              png(file.path(jkp4,paste0("4_IncNodePurity_Plot_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              dotchart(x,y,xlab="IncNodePurity") 
              title(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun))
              dev.off()
              
              # ImportanceSD (standard deviation)
              importancesdd<-modelRf@regModel$importanceSD
              importancesd<-sort(importancesdd,decreasing=TRUE)
              if (i==1) {
                names(importancesd)[names(importancesd)== "P_1_aspect_ea_150"] <- "P1__aspect east"
                names(importancesd)[names(importancesd)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(importancesd)[names(importancesd)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(importancesd)[names(importancesd)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(importancesd)[names(importancesd)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(importancesd)[names(importancesd)== "P_6_cl_grpre1_150"] <- "P6_land use - alpine meadow"
                names(importancesd)[names(importancesd)== "P_7_cl_grpre2_150"] <- "P7_land use - grassland"
                names(importancesd)[names(importancesd)== "P_8_cl_grpre3_150"] <- "P8_land use - grassland second shearing"
                names(importancesd)[names(importancesd)== "P_9_cl_grpre4_150"] <- "P9_land use - grassland special area"
                names(importancesd)[names(importancesd)== "P_10_dtm_150m"] <- "P10_dtm"
                names(importancesd)[names(importancesd)== "P_11_geo_bin_1_150"] <- "P11_geology - shist"
                names(importancesd)[names(importancesd)== "P_12_geo_bin_2_150"] <- "P12_geology - vulcanite"
                names(importancesd)[names(importancesd)== "P_13_geo_bin_3_150"] <- "P13_geology - vinschgau shear zone"
                names(importancesd)[names(importancesd)== "P_14_geo_bin_4_150"] <- "P14_geology - quaternary depositions"
                names(importancesd)[names(importancesd)== "P_15_geo_bin_5_150"] <- "P15_geology - sediment sequences"
                names(importancesd)[names(importancesd)== "P_16_geo_bin_6_150"] <- "P16_geology - plutons"
                names(importancesd)[names(importancesd)== "P_17_geo_bin_7_150"] <- "P17_geology - central gneiss"
                names(importancesd)[names(importancesd)== "P_18_geo_bin_8_150"] <- "P18_geology - quartz phyllite"
                names(importancesd)[names(importancesd)== "P_19_geo_bin_9_150"] <- "P19_geology - insignificant alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_20_geo_bi_10_150"] <- "P20_geology - medium-low alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_21_geo_bi_11_150"] <- "P21_geology - medium-high alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_22_slope_150m"] <- "P22_slope"
                names(importancesd)[names(importancesd)== "P_23_wetness_150m"] <- "P23_wetness_index"
                names(importancesd)[names(importancesd)== "P_24_gve_wgs84_150"] <- "P24_GVE"
                names(importancesd)[names(importancesd)== "P_25_geo_aeubin150"] <- "P25_Acidic soil"
                names(importancesd)[names(importancesd)== "P_26_geo_beubin150"] <- "P26_Basic soil"
                names(importancesd)[names(importancesd)== "P_27_geo_aelbin150"] <- "P27_Acidic soil"
                names(importancesd)[names(importancesd)== "P_28_geo_belbin150"] <- "P28_Basic soil"
                names(importancesd)[names(importancesd)== "P_29_geo_abelbi150"] <- "P29_Acidic-Basic soil"
              } else if (i==2) {
                names(importancesd)[names(importancesd)== "P_1_aspect_ea_150"] <- "P1_aspect east"
                names(importancesd)[names(importancesd)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(importancesd)[names(importancesd)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(importancesd)[names(importancesd)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(importancesd)[names(importancesd)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(importancesd)[names(importancesd)== "P_6_cl_grarp1_150"] <- "P6_land use - alpine meadow"
                names(importancesd)[names(importancesd)== "P_7_cl_grarp2_150"] <- "P7_land use - grassland"
                names(importancesd)[names(importancesd)== "P_8_cl_grarp3_150"] <- "P8_land use - grassland second shearing"
                names(importancesd)[names(importancesd)== "P_9_cl_grarp4_150"] <- "P9land use - grassland special area"
                names(importancesd)[names(importancesd)== "P_10_cl_grarp5_150"] <- "P10_land use - crops"
                names(importancesd)[names(importancesd)== "P_11_dtm_150m"] <- "P11_dtm"
                names(importancesd)[names(importancesd)== "P_12_geo_bin_1_150"] <- "P12_geology - shist"
                names(importancesd)[names(importancesd)== "P_13_geo_bin_2_150"] <- "P13_geology - vulcanite"
                names(importancesd)[names(importancesd)== "P_14_geo_bin_3_150"] <- "P14_geology - vinschgau shear zone"
                names(importancesd)[names(importancesd)== "P_15_geo_bin_4_150"] <- "P15_geology - quaternary depositions"
                names(importancesd)[names(importancesd)== "P_16_geo_bin_5_150"] <- "P16_geology - sediment sequences"
                names(importancesd)[names(importancesd)== "P_17_geo_bin_6_150"] <- "P17_geology - plutons"
                names(importancesd)[names(importancesd)== "P_18_geo_bin_7_150"] <- "P18_geology - central gneiss"
                names(importancesd)[names(importancesd)== "P_19_geo_bin_8_150"] <- "P19_geology - quartz phyllite"
                names(importancesd)[names(importancesd)== "P_20_geo_bin_9_150"] <- "P20_geology - insignificant alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_21_geo_bi_10_150"] <- "P21_geology - medium-low alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_22_geo_bi_11_150"] <- "P22_geology - medium-high alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_23_slope_150m"] <- "P23_slope"
                names(importancesd)[names(importancesd)== "P_24_wetness_150m"] <- "P24_wetness_index"
                names(importancesd)[names(importancesd)== "P_25_gve_wgs84_150"] <- "P25_GVE"
                names(importancesd)[names(importancesd)== "P_26_geo_aeubin150"] <- "P26_Acidic soil"
                names(importancesd)[names(importancesd)== "P_27_geo_beubin150"] <- "P27_Basic soil"
                names(importancesd)[names(importancesd)== "P_28_geo_aelbin150"] <- "P28_Acidic soil"
                names(importancesd)[names(importancesd)== "P_29_geo_belbin150"] <- "P29_Basic soil"
                names(importancesd)[names(importancesd)== "P_30_geo_abelbi150"] <- "P30_Acidic-Basic soil"
              } else {}
              
              # Save PNG Table importancesd
              importancesddff <- data.frame(importancesd)
              names(importancesddff) <- ("ImportanceSD")
              table <- tableGrob(importancesddff)
              title <- textGrob(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia),gp=gpar(fontsize=8))
              footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
              padding <- unit(2,"line")
              table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
              table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
              table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
              png(file.path(jkp4,paste0("4_ImportanceSD_Table_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")), units="cm", width=(ncol(table))*15, height=(nrow(table)), res=300) #, height = 200, width = 100)
              grid.newpage() 
              grid.draw(table) 
              dev.off()
              # Save VALUES Table importancesd
              write.csv(importancesd,file.path(jkp4,paste0("4_ImportanceSD_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".csv")))      
              write.xlsx(importancesd,file.path(jkp4,paste0("4_ImportanceSD_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".xlsx"))) 
              # Save Plot importancesd
              x1<-importancesddff[,1]
              x<-rev(x1)
              y1<-names(importancesd)
              y<-rev(y1)
              png(file.path(jkp4,paste0("4_ImportanceSD_Plot_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              dotchart(x,y,xlab="ImportanceSD") 
              title(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun))
              dev.off()
              
              # Mean of square residuals / % Var. explained
              printregmodel<-modelRf@regModel
              sink(file.path(jkp2,paste0("2_RfModel_Values",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".txt")))
              print(printregmodel)
              sink()
              setwd(paste0(out_dir,"/",predictareas_name[i],"/",var,"/",ai[k],"/","No.Run_",nrun,"/Random_Forest_Regression/2_RfModel"))
              armin<-read.delim(paste0("2_RfModel_Values",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".txt"),stringsAsFactors = FALSE)
              meanofsquareresiduals<-readr::parse_number(armin[5,])
              varexplained<-readr::parse_number(armin[6,])
              
              
              
              # Predictions ORDINARY KRIGING (with 80% of the data)
              predictionsRf<-GSIF::predict.gstatModel(object = modelRf, predictionLocations = grid, predict.method = "RK",verbose=T,nfold = 5)#vgmmodel = model@vgmModel,
              predictRf<-raster(predictionsRf@predicted)
              if (ai[k]=="raw") {
                
              } else if (ai[k]=="log") {
                predictRf_transf_log<-exp(predictRf)
              } else if (ai[k]=="sqrt") {
                predictRf_transf_sqrt<-quad(predictRf)
              } else if (ai[k]=="quad") {
                predictRf_transf_quad<-sqrt(predictRf)
              } else if (ai[k]=="inverse") {
                predictRf_transf_inverse<-inverse(predictRf)
              } else {}
              # Save Plots
              png(file.path(jkpp,paste0("5_Prediction_OrdinaryKriging_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              plot(predictionsRf) #paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              dev.off()
              
              
              
              # Plot Raster of predicted values
              png(file.path(jkppp,paste0("6_Predict_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              if (ai[k]=="raw") {
                raster::plot(predictRf,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else if (ai[k]=="log") {
                raster::plot(predictRf_transf_log,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else if (ai[k]=="sqrt") {
                raster::plot(predictRf_transf_sqrt,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else if (ai[k]=="quad") {
                raster::plot(predictRf_transf_quad,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else if (ai[k]=="inverse") {
                raster::plot(predictRf_transf_inverse,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else {}
              dev.off()
              
              
              
              
              # Overlay validation and prediction
              obsValidRf<-sp::over(obsValid,predictionsRf@predicted)
              if (ai[k]=="raw") {
                
              } else if (ai[k]=="log") {
                obsValidRf$original<-exp(obsValidRf$log)
              } else if (ai[k]=="sqrt") {
                obsValidRf$original<-quad(obsValidRf$sqrt)
              } else if (ai[k]=="quad") {
                obsValidRf$original<-sqrt(obsValidRf$quad)
              } else if (ai[k]=="inverse") {
                obsValidRf$original<-inverse(obsValidRf$inverse)
              } else {}
              
              
              
              # RMSE & NRMSE & R-SQUARED
              if (ai[k]=="raw") {
                rmseRf <- RMSE(pred = obsValidRf[,4], obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf[,4], obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf[,4], obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf[,4], obs = obsValid@data[[var]], na.rm = T)
              } else if (ai[k]=="log") {
                rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              } else if (ai[k]=="sqrt") {
                rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              } else if (ai[k]=="quad") {
                rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              } else if (ai[k]=="inverse") {
                rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              } else {}
              
              
              
              # Maps - Final rasters (GTiff)
              if (ai[k]=="raw") {
                writeRaster(x = predictRf, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else if (ai[k]=="log") {
                writeRaster(x = predictRf_transf_log, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else if (ai[k]=="sqrt") {
                writeRaster(x = predictRf_transf_sqrt, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else if (ai[k]=="quad") {
                writeRaster(x = predictRf_transf_quad, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else if (ai[k]=="inverse") {
                writeRaster(x = predictRf_transf_inverse,filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else {}
              
              # spalte 1 - No.Run
              sdf1<-length(spalte1)
              spalte1<-append(spalte1,nrun,after=sdf1)
              nrun<-nrun+1
              # Spalte 2 - Transformation
              sdf2<-length(spalte2)
              spalte2<-append(spalte2,ai[k],after=sdf2)
              # Spalte 3 - Cutoff
              sdf3<-length(spalte3)
              spalte3<-append(spalte3,taia,after=sdf3)
              # Spalte 4 - Predictors
              sdf4<-length(spalte4)
              spalte4<-append(spalte4,paste(get(paste0("predictors_combi_",l)), collapse = "/"))
              # Spalte 5 - Mean of squared residuals
              sdf5<-length(spalte5)
              spalte5<-append(spalte5,meanofsquareresiduals,after=sdf5)
              # Spalte 6 - % Var explained
              sdf6<-length(spalte6)
              spalte6<-append(spalte6,varexplained,after=sdf6)
              # Spalte 7 - RMSE
              sdf7<-length(spalte7)
              spalte7<-append(spalte7,rmseRf,after=sdf7)
              sdf7tot<-length(spalte7tot)
              spalte7tot<-append(spalte7tot,rmseRf,after=sdf7tot)
              # Spalte 8 - NRMSE
              sdf8<-length(spalte8)
              spalte8<-append(spalte8,nrmseRf,after=sdf8)
              # Spalte 9 - NRMSEMAXMIN
              sdf9<-length(spalte9)
              spalte9<-append(spalte9,nrmseMaxMinRf,after=sdf9)
              # Spalte 10 - R-SQUARED
              sdf10<-length(spalte10)
              spalte10<-append(spalte10,r2Rf,after=sdf10)
              sdf10tot<-length(spalte10tot)
              spalte10tot<-append(spalte10tot,r2Rf,after=sdf10tot)
              
            } else if (s==3) {
              print("length(cutoffs)>=3 / s==3")
              maxcutoff<-max(cutoffs)
              try(modelRf<-fit.gstatModel(observations = obs,formulaString = formula,covariates = grid,method="randomForest",cutoff=maxcutoff))
              taia<-maxcutoff
              # MAIN CODE
              # Save RfModel & VgmModel
              with(faithful,plot(modelRf)) 
              dev.copy(png,file.path(jkp2,paste0("2_RfModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".png")),res=300,width =3500, height=1750)
              dev.off()
              vgmModel = modelRf@vgmModel
              print(vgmModel)
              write.csv(vgmModel,file.path(jkp3,paste0("3_VgmModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".csv")))    
              write.xlsx(vgmModel,file.path(jkp3,paste0("3_VgmModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".xlsx")))    
              sink(file.path(jkp3,paste0("3_VgmModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".txt"))) 
              print(vgmModel)
              sink()
              
              # Save Importance of the predictors (%IncMse & IncNodePurity)
              # PredictorsImportance
              importance = modelRf@regModel$importance
              # %IncMse
              incmse<-sort(importance[,1],decreasing = TRUE)
              # Define names for predictors
              if (i==1) {
                names(incmse)[names(incmse)== "P_1_aspect_ea_150"] <- "P1__aspect east"
                names(incmse)[names(incmse)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(incmse)[names(incmse)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(incmse)[names(incmse)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(incmse)[names(incmse)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(incmse)[names(incmse)== "P_6_cl_grpre1_150"] <- "P6_land use - alpine meadow"
                names(incmse)[names(incmse)== "P_7_cl_grpre2_150"] <- "P7_land use - grassland"
                names(incmse)[names(incmse)== "P_8_cl_grpre3_150"] <- "P8_land use - grassland second shearing"
                names(incmse)[names(incmse)== "P_9_cl_grpre4_150"] <- "P9_land use - grassland special area"
                names(incmse)[names(incmse)== "P_10_dtm_150m"] <- "P10_dtm"
                names(incmse)[names(incmse)== "P_11_geo_bin_1_150"] <- "P11_geology - shist"
                names(incmse)[names(incmse)== "P_12_geo_bin_2_150"] <- "P12_geology - vulcanite"
                names(incmse)[names(incmse)== "P_13_geo_bin_3_150"] <- "P13_geology - vinschgau shear zone"
                names(incmse)[names(incmse)== "P_14_geo_bin_4_150"] <- "P14_geology - quaternary depositions"
                names(incmse)[names(incmse)== "P_15_geo_bin_5_150"] <- "P15_geology - sediment sequences"
                names(incmse)[names(incmse)== "P_16_geo_bin_6_150"] <- "P16_geology - plutons"
                names(incmse)[names(incmse)== "P_17_geo_bin_7_150"] <- "P17_geology - central gneiss"
                names(incmse)[names(incmse)== "P_18_geo_bin_8_150"] <- "P18_geology - quartz phyllite"
                names(incmse)[names(incmse)== "P_19_geo_bin_9_150"] <- "P19_geology - insignificant alpine metamorphosis"
                names(incmse)[names(incmse)== "P_20_geo_bi_10_150"] <- "P20_geology - medium-low alpine metamorphosis"
                names(incmse)[names(incmse)== "P_21_geo_bi_11_150"] <- "P21_geology - medium-high alpine metamorphosis"
                names(incmse)[names(incmse)== "P_22_slope_150m"] <- "P22_slope"
                names(incmse)[names(incmse)== "P_23_wetness_150m"] <- "P23_wetness_index"
                names(incmse)[names(incmse)== "P_24_gve_wgs84_150"] <- "P24_GVE"
                names(incmse)[names(incmse)== "P_25_geo_aeubin150"] <- "P25_Acidic soil"
                names(incmse)[names(incmse)== "P_26_geo_beubin150"] <- "P26_Basic soil"
                names(incmse)[names(incmse)== "P_27_geo_aelbin150"] <- "P27_Acidic soil"
                names(incmse)[names(incmse)== "P_28_geo_belbin150"] <- "P28_Basic soil"
                names(incmse)[names(incmse)== "P_29_geo_abelbi150"] <- "P29_Acidic-Basic soil"
              } else if (i==2) {
                names(incmse)[names(incmse)== "P_1_aspect_ea_150"] <- "P1_aspect east"
                names(incmse)[names(incmse)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(incmse)[names(incmse)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(incmse)[names(incmse)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(incmse)[names(incmse)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(incmse)[names(incmse)== "P_6_cl_grarp1_150"] <- "P6_land use - alpine meadow"
                names(incmse)[names(incmse)== "P_7_cl_grarp2_150"] <- "P7_land use - grassland"
                names(incmse)[names(incmse)== "P_8_cl_grarp3_150"] <- "P8_land use - grassland second shearing"
                names(incmse)[names(incmse)== "P_9_cl_grarp4_150"] <- "P9_land use - grassland special area"
                names(incmse)[names(incmse)== "P_10_cl_grarp5_150"] <- "P10_land use - crops"
                names(incmse)[names(incmse)== "P_11_dtm_150m"] <- "P11_dtm"
                names(incmse)[names(incmse)== "P_12_geo_bin_1_150"] <- "P12_geology - shist"
                names(incmse)[names(incmse)== "P_13_geo_bin_2_150"] <- "P13_geology - vulcanite"
                names(incmse)[names(incmse)== "P_14_geo_bin_3_150"] <- "P14_geology - vinschgau shear zone"
                names(incmse)[names(incmse)== "P_15_geo_bin_4_150"] <- "P15_geology - quaternary depositions"
                names(incmse)[names(incmse)== "P_16_geo_bin_5_150"] <- "P16_geology - sediment sequences"
                names(incmse)[names(incmse)== "P_17_geo_bin_6_150"] <- "P17_geology - plutons"
                names(incmse)[names(incmse)== "P_18_geo_bin_7_150"] <- "P18_geology - central gneiss"
                names(incmse)[names(incmse)== "P_19_geo_bin_8_150"] <- "P19_geology - quartz phyllite"
                names(incmse)[names(incmse)== "P_20_geo_bin_9_150"] <- "P20_geology - insignificant alpine metamorphosis"
                names(incmse)[names(incmse)== "P_21_geo_bi_10_150"] <- "P21_geology - medium-low alpine metamorphosis"
                names(incmse)[names(incmse)== "P_22_geo_bi_11_150"] <- "P22_geology - medium-high alpine metamorphosis"
                names(incmse)[names(incmse)== "P_23_slope_150m"] <- "P23_slope"
                names(incmse)[names(incmse)== "P_24_wetness_150m"] <- "P24_wetness_index"
                names(incmse)[names(incmse)== "P_25_gve_wgs84_150"] <- "P25_GVE"
                names(incmse)[names(incmse)== "P_26_geo_aeubin150"] <- "P26_Acidic soil"
                names(incmse)[names(incmse)== "P_27_geo_beubin150"] <- "P27_Basic soil"
                names(incmse)[names(incmse)== "P_28_geo_aelbin150"] <- "P28_Acidic soil"
                names(incmse)[names(incmse)== "P_29_geo_belbin150"] <- "P29_Basic soil"
                names(incmse)[names(incmse)== "P_30_geo_abelbi150"] <- "P30_Acidic-Basic soil"
              } else {}
              
              # Save PNG Table %IncMSE
              incmsedff <- data.frame(incmse)
              names(incmsedff) <- ("%IncMSE")
              table <- tableGrob(incmsedff)
              title <- textGrob(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia),gp=gpar(fontsize=8))
              footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
              padding <- unit(2,"line")
              table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
              table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
              table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
              png(file.path(jkp4,paste0("4_IncMse_Table_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")), units="cm", width=(ncol(table))*15, height=(nrow(table)), res=300) #, height = 200, width = 100)
              grid.newpage() 
              grid.draw(table) 
              dev.off()
              # Save VALUES Table %IncMSE
              write.csv(incmse,file.path(jkp4,paste0("4_IncMse_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".csv")))      
              write.xlsx(incmse,file.path(jkp4,paste0("4_IncMse_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".xlsx")))    
              # Save Plot %IncMSE
              x1<-incmsedff[,1]
              x<-rev(x1)
              y1<-names(incmse)
              y<-rev(y1)
              png(file.path(jkp4,paste0("4_IncMse_Plot_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              dotchart(x,y,xlab="%IncMse") 
              title(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun))
              dev.off()
              
              # IncNodePurity
              IncNodePurity<-round(sort(importance[,2],decreasing = TRUE),digits = 2) 
              # Define names for predictors
              if (i==1) {
                names(IncNodePurity)[names(IncNodePurity)== "P_1_aspect_ea_150"] <- "P1__aspect east"
                names(IncNodePurity)[names(IncNodePurity)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(IncNodePurity)[names(IncNodePurity)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(IncNodePurity)[names(IncNodePurity)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(IncNodePurity)[names(IncNodePurity)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(IncNodePurity)[names(IncNodePurity)== "P_6_cl_grpre1_150"] <- "P6_land use - alpine meadow"
                names(IncNodePurity)[names(IncNodePurity)== "P_7_cl_grpre2_150"] <- "P7_land use - grassland"
                names(IncNodePurity)[names(IncNodePurity)== "P_8_cl_grpre3_150"] <- "P8_land use - grassland second shearing"
                names(IncNodePurity)[names(IncNodePurity)== "P_9_cl_grpre4_150"] <- "P9_land use - grassland special area"
                names(IncNodePurity)[names(IncNodePurity)== "P_10_dtm_150m"] <- "P10_dtm"
                names(IncNodePurity)[names(IncNodePurity)== "P_11_geo_bin_1_150"] <- "P11_geology - shist"
                names(IncNodePurity)[names(IncNodePurity)== "P_12_geo_bin_2_150"] <- "P12_geology - vulcanite"
                names(IncNodePurity)[names(IncNodePurity)== "P_13_geo_bin_3_150"] <- "P13_geology - vinschgau shear zone"
                names(IncNodePurity)[names(IncNodePurity)== "P_14_geo_bin_4_150"] <- "P14_geology - quaternary depositions"
                names(IncNodePurity)[names(IncNodePurity)== "P_15_geo_bin_5_150"] <- "P15_geology - sediment sequences"
                names(IncNodePurity)[names(IncNodePurity)== "P_16_geo_bin_6_150"] <- "P16_geology - plutons"
                names(IncNodePurity)[names(IncNodePurity)== "P_17_geo_bin_7_150"] <- "P17_geology - central gneiss"
                names(IncNodePurity)[names(IncNodePurity)== "P_18_geo_bin_8_150"] <- "P18_geology - quartz phyllite"
                names(IncNodePurity)[names(IncNodePurity)== "P_19_geo_bin_9_150"] <- "P19_geology - insignificant alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_20_geo_bi_10_150"] <- "P20_geology - medium-low alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_21_geo_bi_11_150"] <- "P21_geology - medium-high alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_22_slope_150m"] <- "P22_slope"
                names(IncNodePurity)[names(IncNodePurity)== "P_23_wetness_150m"] <- "P23_wetness_index"
                names(IncNodePurity)[names(IncNodePurity)== "P_24_gve_wgs84_150"] <- "P24_GVE"
                names(IncNodePurity)[names(IncNodePurity)== "P_25_geo_aeubin150"] <- "P25_Acidic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_26_geo_beubin150"] <- "P26_Basic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_27_geo_aelbin150"] <- "P27_Acidic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_28_geo_belbin150"] <- "P28_Basic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_29_geo_abelbi150"] <- "P29_Acidic-Basic soil"
              } else if (i==2) {
                names(IncNodePurity)[names(IncNodePurity)== "P_1_aspect_ea_150"] <- "P1_aspect east"
                names(IncNodePurity)[names(IncNodePurity)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(IncNodePurity)[names(IncNodePurity)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(IncNodePurity)[names(IncNodePurity)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(IncNodePurity)[names(IncNodePurity)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(IncNodePurity)[names(IncNodePurity)== "P_6_cl_grarp1_150"] <- "P6_land use - alpine meadow"
                names(IncNodePurity)[names(IncNodePurity)== "P_7_cl_grarp2_150"] <- "P7_land use - grassland"
                names(IncNodePurity)[names(IncNodePurity)== "P_8_cl_grarp3_150"] <- "P8_land use - grassland second shearing"
                names(IncNodePurity)[names(IncNodePurity)== "P_9_cl_grarp4_150"] <- "P9land use - grassland special area"
                names(IncNodePurity)[names(IncNodePurity)== "P_10_cl_grarp5_150"] <- "P10_land use - crops"
                names(IncNodePurity)[names(IncNodePurity)== "P_11_dtm_150m"] <- "P11_dtm"
                names(IncNodePurity)[names(IncNodePurity)== "P_12_geo_bin_1_150"] <- "P12_geology - shist"
                names(IncNodePurity)[names(IncNodePurity)== "P_13_geo_bin_2_150"] <- "P13_geology - vulcanite"
                names(IncNodePurity)[names(IncNodePurity)== "P_14_geo_bin_3_150"] <- "P14_geology - vinschgau shear zone"
                names(IncNodePurity)[names(IncNodePurity)== "P_15_geo_bin_4_150"] <- "P15_geology - quaternary depositions"
                names(IncNodePurity)[names(IncNodePurity)== "P_16_geo_bin_5_150"] <- "P16_geology - sediment sequences"
                names(IncNodePurity)[names(IncNodePurity)== "P_17_geo_bin_6_150"] <- "P17_geology - plutons"
                names(IncNodePurity)[names(IncNodePurity)== "P_18_geo_bin_7_150"] <- "P18_geology - central gneiss"
                names(IncNodePurity)[names(IncNodePurity)== "P_19_geo_bin_8_150"] <- "P19_geology - quartz phyllite"
                names(IncNodePurity)[names(IncNodePurity)== "P_20_geo_bin_9_150"] <- "P20_geology - insignificant alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_21_geo_bi_10_150"] <- "P21_geology - medium-low alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_22_geo_bi_11_150"] <- "P22_geology - medium-high alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_23_slope_150m"] <- "P23_slope"
                names(IncNodePurity)[names(IncNodePurity)== "P_24_wetness_150m"] <- "P24_wetness_index"
                names(IncNodePurity)[names(IncNodePurity)== "P_25_gve_wgs84_150"] <- "P25_GVE"
                names(IncNodePurity)[names(IncNodePurity)== "P_26_geo_aeubin150"] <- "P26_Acidic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_27_geo_beubin150"] <- "P27_Basic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_28_geo_aelbin150"] <- "P28_Acidic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_29_geo_belbin150"] <- "P29_Basic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_30_geo_abelbi150"] <- "P30_Acidic-Basic soil"
              } else {}
              
              # Save PNG Table IncNodePurity
              incnodedff <- data.frame(IncNodePurity)
              names(incnodedff) <- ("IncNodePurity")
              table <- tableGrob(incnodedff)
              title <- textGrob(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia),gp=gpar(fontsize=8))
              footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
              padding <- unit(2,"line")
              table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
              table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
              table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
              png(file.path(jkp4,paste0("4_IncNodePurity_Table_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")), units="cm", width=(ncol(table))*15, height=(nrow(table)), res=300) #, height = 200, width = 100)
              grid.newpage() 
              grid.draw(table) 
              dev.off()
              # Save VALUES Table IncNodePurity
              write.csv(IncNodePurity,file.path(jkp4,paste0("4_IncNodePurity_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".csv")))      
              write.xlsx(IncNodePurity,file.path(jkp4,paste0("4_IncNodePurity_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".xlsx"))) 
              # Save Plot IncNodePurity
              x1<-incnodedff[,1]
              x<-rev(x1)
              y1<-names(IncNodePurity)
              y<-rev(y1)
              png(file.path(jkp4,paste0("4_IncNodePurity_Plot_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              dotchart(x,y,xlab="IncNodePurity") 
              title(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun))
              dev.off()
              
              # ImportanceSD (standard deviation)
              importancesdd<-modelRf@regModel$importanceSD
              importancesd<-sort(importancesdd,decreasing=TRUE)
              if (i==1) {
                names(importancesd)[names(importancesd)== "P_1_aspect_ea_150"] <- "P1__aspect east"
                names(importancesd)[names(importancesd)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(importancesd)[names(importancesd)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(importancesd)[names(importancesd)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(importancesd)[names(importancesd)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(importancesd)[names(importancesd)== "P_6_cl_grpre1_150"] <- "P6_land use - alpine meadow"
                names(importancesd)[names(importancesd)== "P_7_cl_grpre2_150"] <- "P7_land use - grassland"
                names(importancesd)[names(importancesd)== "P_8_cl_grpre3_150"] <- "P8_land use - grassland second shearing"
                names(importancesd)[names(importancesd)== "P_9_cl_grpre4_150"] <- "P9_land use - grassland special area"
                names(importancesd)[names(importancesd)== "P_10_dtm_150m"] <- "P10_dtm"
                names(importancesd)[names(importancesd)== "P_11_geo_bin_1_150"] <- "P11_geology - shist"
                names(importancesd)[names(importancesd)== "P_12_geo_bin_2_150"] <- "P12_geology - vulcanite"
                names(importancesd)[names(importancesd)== "P_13_geo_bin_3_150"] <- "P13_geology - vinschgau shear zone"
                names(importancesd)[names(importancesd)== "P_14_geo_bin_4_150"] <- "P14_geology - quaternary depositions"
                names(importancesd)[names(importancesd)== "P_15_geo_bin_5_150"] <- "P15_geology - sediment sequences"
                names(importancesd)[names(importancesd)== "P_16_geo_bin_6_150"] <- "P16_geology - plutons"
                names(importancesd)[names(importancesd)== "P_17_geo_bin_7_150"] <- "P17_geology - central gneiss"
                names(importancesd)[names(importancesd)== "P_18_geo_bin_8_150"] <- "P18_geology - quartz phyllite"
                names(importancesd)[names(importancesd)== "P_19_geo_bin_9_150"] <- "P19_geology - insignificant alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_20_geo_bi_10_150"] <- "P20_geology - medium-low alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_21_geo_bi_11_150"] <- "P21_geology - medium-high alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_22_slope_150m"] <- "P22_slope"
                names(importancesd)[names(importancesd)== "P_23_wetness_150m"] <- "P23_wetness_index"
                names(importancesd)[names(importancesd)== "P_24_gve_wgs84_150"] <- "P24_GVE"
                names(importancesd)[names(importancesd)== "P_25_geo_aeubin150"] <- "P25_Acidic soil"
                names(importancesd)[names(importancesd)== "P_26_geo_beubin150"] <- "P26_Basic soil"
                names(importancesd)[names(importancesd)== "P_27_geo_aelbin150"] <- "P27_Acidic soil"
                names(importancesd)[names(importancesd)== "P_28_geo_belbin150"] <- "P28_Basic soil"
                names(importancesd)[names(importancesd)== "P_29_geo_abelbi150"] <- "P29_Acidic-Basic soil"
              } else if (i==2) {
                names(importancesd)[names(importancesd)== "P_1_aspect_ea_150"] <- "P1_aspect east"
                names(importancesd)[names(importancesd)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(importancesd)[names(importancesd)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(importancesd)[names(importancesd)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(importancesd)[names(importancesd)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(importancesd)[names(importancesd)== "P_6_cl_grarp1_150"] <- "P6_land use - alpine meadow"
                names(importancesd)[names(importancesd)== "P_7_cl_grarp2_150"] <- "P7_land use - grassland"
                names(importancesd)[names(importancesd)== "P_8_cl_grarp3_150"] <- "P8_land use - grassland second shearing"
                names(importancesd)[names(importancesd)== "P_9_cl_grarp4_150"] <- "P9land use - grassland special area"
                names(importancesd)[names(importancesd)== "P_10_cl_grarp5_150"] <- "P10_land use - crops"
                names(importancesd)[names(importancesd)== "P_11_dtm_150m"] <- "P11_dtm"
                names(importancesd)[names(importancesd)== "P_12_geo_bin_1_150"] <- "P12_geology - shist"
                names(importancesd)[names(importancesd)== "P_13_geo_bin_2_150"] <- "P13_geology - vulcanite"
                names(importancesd)[names(importancesd)== "P_14_geo_bin_3_150"] <- "P14_geology - vinschgau shear zone"
                names(importancesd)[names(importancesd)== "P_15_geo_bin_4_150"] <- "P15_geology - quaternary depositions"
                names(importancesd)[names(importancesd)== "P_16_geo_bin_5_150"] <- "P16_geology - sediment sequences"
                names(importancesd)[names(importancesd)== "P_17_geo_bin_6_150"] <- "P17_geology - plutons"
                names(importancesd)[names(importancesd)== "P_18_geo_bin_7_150"] <- "P18_geology - central gneiss"
                names(importancesd)[names(importancesd)== "P_19_geo_bin_8_150"] <- "P19_geology - quartz phyllite"
                names(importancesd)[names(importancesd)== "P_20_geo_bin_9_150"] <- "P20_geology - insignificant alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_21_geo_bi_10_150"] <- "P21_geology - medium-low alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_22_geo_bi_11_150"] <- "P22_geology - medium-high alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_23_slope_150m"] <- "P23_slope"
                names(importancesd)[names(importancesd)== "P_24_wetness_150m"] <- "P24_wetness_index"
                names(importancesd)[names(importancesd)== "P_25_gve_wgs84_150"] <- "P25_GVE"
                names(importancesd)[names(importancesd)== "P_26_geo_aeubin150"] <- "P26_Acidic soil"
                names(importancesd)[names(importancesd)== "P_27_geo_beubin150"] <- "P27_Basic soil"
                names(importancesd)[names(importancesd)== "P_28_geo_aelbin150"] <- "P28_Acidic soil"
                names(importancesd)[names(importancesd)== "P_29_geo_belbin150"] <- "P29_Basic soil"
                names(importancesd)[names(importancesd)== "P_30_geo_abelbi150"] <- "P30_Acidic-Basic soil"
              } else {}
              
              # Save PNG Table importancesd
              importancesddff <- data.frame(importancesd)
              names(importancesddff) <- ("ImportanceSD")
              table <- tableGrob(importancesddff)
              title <- textGrob(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia),gp=gpar(fontsize=8))
              footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
              padding <- unit(2,"line")
              table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
              table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
              table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
              png(file.path(jkp4,paste0("4_ImportanceSD_Table_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")), units="cm", width=(ncol(table))*15, height=(nrow(table)), res=300) #, height = 200, width = 100)
              grid.newpage() 
              grid.draw(table) 
              dev.off()
              # Save VALUES Table importancesd
              write.csv(importancesd,file.path(jkp4,paste0("4_ImportanceSD_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".csv")))      
              write.xlsx(importancesd,file.path(jkp4,paste0("4_ImportanceSD_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".xlsx"))) 
              # Save Plot importancesd
              x1<-importancesddff[,1]
              x<-rev(x1)
              y1<-names(importancesd)
              y<-rev(y1)
              png(file.path(jkp4,paste0("4_ImportanceSD_Plot_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              dotchart(x,y,xlab="ImportanceSD") 
              title(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun))
              dev.off()
              
              # Mean of square residuals / % Var. explained
              printregmodel<-modelRf@regModel
              sink(file.path(jkp2,paste0("2_RfModel_Values",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".txt")))
              print(printregmodel)
              sink()
              setwd(paste0(out_dir,"/",predictareas_name[i],"/",var,"/",ai[k],"/","No.Run_",nrun,"/Random_Forest_Regression/2_RfModel"))
              armin<-read.delim(paste0("2_RfModel_Values",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".txt"),stringsAsFactors = FALSE)
              meanofsquareresiduals<-readr::parse_number(armin[5,])
              varexplained<-readr::parse_number(armin[6,])
              
              
              
              # Predictions ORDINARY KRIGING 
              predictionsRf<-GSIF::predict.gstatModel(object = modelRf, predictionLocations = grid, predict.method = "RK",verbose=T,nfold = 5)#vgmmodel = model@vgmModel,
              predictRf<-raster(predictionsRf@predicted)
              if (ai[k]=="raw") {
                
              } else if (ai[k]=="log") {
                predictRf_transf_log<-exp(predictRf)
              } else if (ai[k]=="sqrt") {
                predictRf_transf_sqrt<-quad(predictRf)
              } else if (ai[k]=="quad") {
                predictRf_transf_quad<-sqrt(predictRf)
              } else if (ai[k]=="inverse") {
                predictRf_transf_inverse<-inverse(predictRf)
              } else {}
              # Save Plots
              png(file.path(jkpp,paste0("5_Prediction_OrdinaryKriging_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              plot(predictionsRf) #paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              dev.off()
              
              
              
              # Plot Raster of predicted values
              png(file.path(jkppp,paste0("6_Predict_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              if (ai[k]=="raw") {
                raster::plot(predictRf,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else if (ai[k]=="log") {
                raster::plot(predictRf_transf_log,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else if (ai[k]=="sqrt") {
                raster::plot(predictRf_transf_sqrt,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else if (ai[k]=="quad") {
                raster::plot(predictRf_transf_quad,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else if (ai[k]=="inverse") {
                raster::plot(predictRf_transf_inverse,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else {}
              dev.off()
              
              
              
              
              # Overlay validation and prediction
              obsValidRf<-sp::over(obsValid,predictionsRf@predicted)
              if (ai[k]=="raw") {
                
              } else if (ai[k]=="log") {
                obsValidRf$original<-exp(obsValidRf$log)
              } else if (ai[k]=="sqrt") {
                obsValidRf$original<-quad(obsValidRf$sqrt)
              } else if (ai[k]=="quad") {
                obsValidRf$original<-sqrt(obsValidRf$quad)
              } else if (ai[k]=="inverse") {
                obsValidRf$original<-inverse(obsValidRf$inverse)
              } else {}
              
              
              
              # RMSE & NRMSE & R-SQUARED
              if (ai[k]=="raw") {
                rmseRf <- RMSE(pred = obsValidRf[,4], obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf[,4], obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf[,4], obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf[,4], obs = obsValid@data[[var]], na.rm = T)
              } else if (ai[k]=="log") {
                rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              } else if (ai[k]=="sqrt") {
                rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              } else if (ai[k]=="quad") {
                rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              } else if (ai[k]=="inverse") {
                rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              } else {}
              
              
              
              # Maps - Final rasters (GTiff)
              if (ai[k]=="raw") {
                writeRaster(x = predictRf, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else if (ai[k]=="log") {
                writeRaster(x = predictRf_transf_log, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else if (ai[k]=="sqrt") {
                writeRaster(x = predictRf_transf_sqrt, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else if (ai[k]=="quad") {
                writeRaster(x = predictRf_transf_quad, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else if (ai[k]=="inverse") {
                writeRaster(x = predictRf_transf_inverse,filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else {}
              
              # spalte 1 - No.Run
              sdf1<-length(spalte1)
              spalte1<-append(spalte1,nrun,after=sdf1)
              nrun<-nrun+1
              # Spalte 2 - Transformation
              sdf2<-length(spalte2)
              spalte2<-append(spalte2,ai[k],after=sdf2)
              # Spalte 3 - Cutoff
              sdf3<-length(spalte3)
              spalte3<-append(spalte3,taia,after=sdf3)
              # Spalte 4 - Predictors
              sdf4<-length(spalte4)
              spalte4<-append(spalte4,paste(get(paste0("predictors_combi_",l)), collapse = "/"))
              # Spalte 5 - Mean of squared residuals
              sdf5<-length(spalte5)
              spalte5<-append(spalte5,meanofsquareresiduals,after=sdf5)
              # Spalte 6 - % Var explained
              sdf6<-length(spalte6)
              spalte6<-append(spalte6,varexplained,after=sdf6)
              # Spalte 7 - RMSE
              sdf7<-length(spalte7)
              spalte7<-append(spalte7,rmseRf,after=sdf7)
              sdf7tot<-length(spalte7tot)
              spalte7tot<-append(spalte7tot,rmseRf,after=sdf7tot)
              # Spalte 8 - NRMSE
              sdf8<-length(spalte8)
              spalte8<-append(spalte8,nrmseRf,after=sdf8)
              # Spalte 9 - NRMSEMAXMIN
              sdf9<-length(spalte9)
              spalte9<-append(spalte9,nrmseMaxMinRf,after=sdf9)
              # Spalte 10 - R-SQUARED
              sdf10<-length(spalte10)
              spalte10<-append(spalte10,r2Rf,after=sdf10)
              sdf10tot<-length(spalte10tot)
              spalte10tot<-append(spalte10tot,r2Rf,after=sdf10tot)
              
            } else {}
            
            
          } else if (length(cutoffs)==2) {
            print("length(cutoffs)=2")
            # Output folders
            fg<-paste0("/output/",predictareas_name[i],"/",variables[j],"/",ai[k])
            #fh<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),fg)
            fh<-paste0(script_dir,fg)
            setwd(fh)
            dir.create(paste0("No.Run_",nrun))
            fgg<-paste0("/output/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun)
            #fhh<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),fgg)
            fhh<-paste0(script_dir,fgg)
            setwd(fhh)
            dir.create("Random_Forest_Regression")
            jkp<-paste0(out_dir,"/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun,"/Random_Forest_Regression")
            setwd(jkp)
            dir.create("1_Variogram")
            jkp1<-paste0(out_dir,"/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun,"/Random_Forest_Regression/1_Variogram")
            dir.create("2_RfModel")
            jkp2<-paste0(out_dir,"/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun,"/Random_Forest_Regression/2_RfModel")
            dir.create("3_VgmModel")
            jkp3<-paste0(out_dir,"/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun,"/Random_Forest_Regression/3_VgmModel")
            dir.create("4_PredictorsImportance")
            jkp4<-paste0(out_dir,"/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun,"/Random_Forest_Regression/4_PredictorsImportance")
            setwd(fhh)
            dir.create("Ordinary_Kriging_Prediction")
            jkpp<-paste0(out_dir,"/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun,"/Ordinary_Kriging_Prediction")
            dir.create("Final_Models")
            jkppp<-paste0(out_dir,"/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun,"/Final_Models")
            
            # Variogram
            #plot(variogram(formula_variabile, obs), plot.nu=T, pch="+")
            x=variogram(formula_variabile, obs)
            plott = ggplot(x,aes(dist,gamma))+
              geom_point()+
              theme_minimal()+
              ggtitle(paste0("Variogram - ","[",variables[j],"] ","- ",ai[k],"-","predictarea",i))+
              theme(plot.title = element_text(hjust = 0.5))
            jkpt<-paste0("Variogram - ","[",variables[j],"] ","- ",ai[k],"-","predictarea",i)
            #png(file.path(jkp1,paste0("1_Variogram1_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =4000, height=4000)
            #plot(variogram(formula_variabile, obs), plot.nu=T, pch="+",main=jkpt)
            #dev.off()
            ggsave(filename = file.path(jkp1,paste0("1_Variogram_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),plot = plott,device = "png",width=5, height=5)
            
            
            if (s==1) {
              print("length(cutoffs)=2 / s==1")
              mincutoff<-min(cutoffs)
              try(modelRf<-fit.gstatModel(observations = obs,formulaString = formula,covariates = grid,method="randomForest",cutoff=mincutoff))
              taia<-mincutoff
              # MAIN CODE
              # Save RfModel & VgmModel
              with(faithful,plot(modelRf)) 
              dev.copy(png,file.path(jkp2,paste0("2_RfModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".png")),res=300,width =3500, height=1750)
              dev.off()
              vgmModel = modelRf@vgmModel
              print(vgmModel)
              write.csv(vgmModel,file.path(jkp3,paste0("3_VgmModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".csv")))    
              write.xlsx(vgmModel,file.path(jkp3,paste0("3_VgmModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".xlsx")))    
              sink(file.path(jkp3,paste0("3_VgmModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".txt"))) 
              print(vgmModel)
              sink()
              
              # Save Importance of the predictors (%IncMse & IncNodePurity)
              # PredictorsImportance
              importance = modelRf@regModel$importance
              # %IncMse
              incmse<-sort(importance[,1],decreasing = TRUE)
              # Define names for predictors
              if (i==1) {
                names(incmse)[names(incmse)== "P_1_aspect_ea_150"] <- "P1__aspect east"
                names(incmse)[names(incmse)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(incmse)[names(incmse)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(incmse)[names(incmse)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(incmse)[names(incmse)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(incmse)[names(incmse)== "P_6_cl_grpre1_150"] <- "P6_land use - alpine meadow"
                names(incmse)[names(incmse)== "P_7_cl_grpre2_150"] <- "P7_land use - grassland"
                names(incmse)[names(incmse)== "P_8_cl_grpre3_150"] <- "P8_land use - grassland second shearing"
                names(incmse)[names(incmse)== "P_9_cl_grpre4_150"] <- "P9_land use - grassland special area"
                names(incmse)[names(incmse)== "P_10_dtm_150m"] <- "P10_dtm"
                names(incmse)[names(incmse)== "P_11_geo_bin_1_150"] <- "P11_geology - shist"
                names(incmse)[names(incmse)== "P_12_geo_bin_2_150"] <- "P12_geology - vulcanite"
                names(incmse)[names(incmse)== "P_13_geo_bin_3_150"] <- "P13_geology - vinschgau shear zone"
                names(incmse)[names(incmse)== "P_14_geo_bin_4_150"] <- "P14_geology - quaternary depositions"
                names(incmse)[names(incmse)== "P_15_geo_bin_5_150"] <- "P15_geology - sediment sequences"
                names(incmse)[names(incmse)== "P_16_geo_bin_6_150"] <- "P16_geology - plutons"
                names(incmse)[names(incmse)== "P_17_geo_bin_7_150"] <- "P17_geology - central gneiss"
                names(incmse)[names(incmse)== "P_18_geo_bin_8_150"] <- "P18_geology - quartz phyllite"
                names(incmse)[names(incmse)== "P_19_geo_bin_9_150"] <- "P19_geology - insignificant alpine metamorphosis"
                names(incmse)[names(incmse)== "P_20_geo_bi_10_150"] <- "P20_geology - medium-low alpine metamorphosis"
                names(incmse)[names(incmse)== "P_21_geo_bi_11_150"] <- "P21_geology - medium-high alpine metamorphosis"
                names(incmse)[names(incmse)== "P_22_slope_150m"] <- "P22_slope"
                names(incmse)[names(incmse)== "P_23_wetness_150m"] <- "P23_wetness_index"
                names(incmse)[names(incmse)== "P_24_gve_wgs84_150"] <- "P24_GVE"
                names(incmse)[names(incmse)== "P_25_geo_aeubin150"] <- "P25_Acidic soil"
                names(incmse)[names(incmse)== "P_26_geo_beubin150"] <- "P26_Basic soil"
                names(incmse)[names(incmse)== "P_27_geo_aelbin150"] <- "P27_Acidic soil"
                names(incmse)[names(incmse)== "P_28_geo_belbin150"] <- "P28_Basic soil"
                names(incmse)[names(incmse)== "P_29_geo_abelbi150"] <- "P29_Acidic-Basic soil"
              } else if (i==2) {
                names(incmse)[names(incmse)== "P_1_aspect_ea_150"] <- "P1_aspect east"
                names(incmse)[names(incmse)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(incmse)[names(incmse)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(incmse)[names(incmse)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(incmse)[names(incmse)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(incmse)[names(incmse)== "P_6_cl_grarp1_150"] <- "P6_land use - alpine meadow"
                names(incmse)[names(incmse)== "P_7_cl_grarp2_150"] <- "P7_land use - grassland"
                names(incmse)[names(incmse)== "P_8_cl_grarp3_150"] <- "P8_land use - grassland second shearing"
                names(incmse)[names(incmse)== "P_9_cl_grarp4_150"] <- "P9_land use - grassland special area"
                names(incmse)[names(incmse)== "P_10_cl_grarp5_150"] <- "P10_land use - crops"
                names(incmse)[names(incmse)== "P_11_dtm_150m"] <- "P11_dtm"
                names(incmse)[names(incmse)== "P_12_geo_bin_1_150"] <- "P12_geology - shist"
                names(incmse)[names(incmse)== "P_13_geo_bin_2_150"] <- "P13_geology - vulcanite"
                names(incmse)[names(incmse)== "P_14_geo_bin_3_150"] <- "P14_geology - vinschgau shear zone"
                names(incmse)[names(incmse)== "P_15_geo_bin_4_150"] <- "P15_geology - quaternary depositions"
                names(incmse)[names(incmse)== "P_16_geo_bin_5_150"] <- "P16_geology - sediment sequences"
                names(incmse)[names(incmse)== "P_17_geo_bin_6_150"] <- "P17_geology - plutons"
                names(incmse)[names(incmse)== "P_18_geo_bin_7_150"] <- "P18_geology - central gneiss"
                names(incmse)[names(incmse)== "P_19_geo_bin_8_150"] <- "P19_geology - quartz phyllite"
                names(incmse)[names(incmse)== "P_20_geo_bin_9_150"] <- "P20_geology - insignificant alpine metamorphosis"
                names(incmse)[names(incmse)== "P_21_geo_bi_10_150"] <- "P21_geology - medium-low alpine metamorphosis"
                names(incmse)[names(incmse)== "P_22_geo_bi_11_150"] <- "P22_geology - medium-high alpine metamorphosis"
                names(incmse)[names(incmse)== "P_23_slope_150m"] <- "P23_slope"
                names(incmse)[names(incmse)== "P_24_wetness_150m"] <- "P24_wetness_index"
                names(incmse)[names(incmse)== "P_25_gve_wgs84_150"] <- "P25_GVE"
                names(incmse)[names(incmse)== "P_26_geo_aeubin150"] <- "P26_Acidic soil"
                names(incmse)[names(incmse)== "P_27_geo_beubin150"] <- "P27_Basic soil"
                names(incmse)[names(incmse)== "P_28_geo_aelbin150"] <- "P28_Acidic soil"
                names(incmse)[names(incmse)== "P_29_geo_belbin150"] <- "P29_Basic soil"
                names(incmse)[names(incmse)== "P_30_geo_abelbi150"] <- "P30_Acidic-Basic soil"
              } else {}
              
              # Save PNG Table %IncMSE
              incmsedff <- data.frame(incmse)
              names(incmsedff) <- ("%IncMSE")
              table <- tableGrob(incmsedff)
              title <- textGrob(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia),gp=gpar(fontsize=8))
              footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
              padding <- unit(2,"line")
              table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
              table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
              table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
              png(file.path(jkp4,paste0("4_IncMse_Table_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")), units="cm", width=(ncol(table))*15, height=(nrow(table)), res=300) #, height = 200, width = 100)
              grid.newpage() 
              grid.draw(table) 
              dev.off()
              # Save VALUES Table %IncMSE
              write.csv(incmse,file.path(jkp4,paste0("4_IncMse_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".csv")))      
              write.xlsx(incmse,file.path(jkp4,paste0("4_IncMse_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".xlsx")))    
              # Save Plot %IncMSE
              x1<-incmsedff[,1]
              x<-rev(x1)
              y1<-names(incmse)
              y<-rev(y1)
              png(file.path(jkp4,paste0("4_IncMse_Plot_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              dotchart(x,y,xlab="%IncMse") 
              title(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun))
              dev.off()
              
              # IncNodePurity
              IncNodePurity<-round(sort(importance[,2],decreasing = TRUE),digits = 2) 
              # Define names for predictors
              if (i==1) {
                names(IncNodePurity)[names(IncNodePurity)== "P_1_aspect_ea_150"] <- "P1__aspect east"
                names(IncNodePurity)[names(IncNodePurity)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(IncNodePurity)[names(IncNodePurity)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(IncNodePurity)[names(IncNodePurity)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(IncNodePurity)[names(IncNodePurity)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(IncNodePurity)[names(IncNodePurity)== "P_6_cl_grpre1_150"] <- "P6_land use - alpine meadow"
                names(IncNodePurity)[names(IncNodePurity)== "P_7_cl_grpre2_150"] <- "P7_land use - grassland"
                names(IncNodePurity)[names(IncNodePurity)== "P_8_cl_grpre3_150"] <- "P8_land use - grassland second shearing"
                names(IncNodePurity)[names(IncNodePurity)== "P_9_cl_grpre4_150"] <- "P9_land use - grassland special area"
                names(IncNodePurity)[names(IncNodePurity)== "P_10_dtm_150m"] <- "P10_dtm"
                names(IncNodePurity)[names(IncNodePurity)== "P_11_geo_bin_1_150"] <- "P11_geology - shist"
                names(IncNodePurity)[names(IncNodePurity)== "P_12_geo_bin_2_150"] <- "P12_geology - vulcanite"
                names(IncNodePurity)[names(IncNodePurity)== "P_13_geo_bin_3_150"] <- "P13_geology - vinschgau shear zone"
                names(IncNodePurity)[names(IncNodePurity)== "P_14_geo_bin_4_150"] <- "P14_geology - quaternary depositions"
                names(IncNodePurity)[names(IncNodePurity)== "P_15_geo_bin_5_150"] <- "P15_geology - sediment sequences"
                names(IncNodePurity)[names(IncNodePurity)== "P_16_geo_bin_6_150"] <- "P16_geology - plutons"
                names(IncNodePurity)[names(IncNodePurity)== "P_17_geo_bin_7_150"] <- "P17_geology - central gneiss"
                names(IncNodePurity)[names(IncNodePurity)== "P_18_geo_bin_8_150"] <- "P18_geology - quartz phyllite"
                names(IncNodePurity)[names(IncNodePurity)== "P_19_geo_bin_9_150"] <- "P19_geology - insignificant alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_20_geo_bi_10_150"] <- "P20_geology - medium-low alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_21_geo_bi_11_150"] <- "P21_geology - medium-high alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_22_slope_150m"] <- "P22_slope"
                names(IncNodePurity)[names(IncNodePurity)== "P_23_wetness_150m"] <- "P23_wetness_index"
                names(IncNodePurity)[names(IncNodePurity)== "P_24_gve_wgs84_150"] <- "P24_GVE"
                names(IncNodePurity)[names(IncNodePurity)== "P_25_geo_aeubin150"] <- "P25_Acidic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_26_geo_beubin150"] <- "P26_Basic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_27_geo_aelbin150"] <- "P27_Acidic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_28_geo_belbin150"] <- "P28_Basic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_29_geo_abelbi150"] <- "P29_Acidic-Basic soil"
              } else if (i==2) {
                names(IncNodePurity)[names(IncNodePurity)== "P_1_aspect_ea_150"] <- "P1_aspect east"
                names(IncNodePurity)[names(IncNodePurity)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(IncNodePurity)[names(IncNodePurity)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(IncNodePurity)[names(IncNodePurity)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(IncNodePurity)[names(IncNodePurity)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(IncNodePurity)[names(IncNodePurity)== "P_6_cl_grarp1_150"] <- "P6_land use - alpine meadow"
                names(IncNodePurity)[names(IncNodePurity)== "P_7_cl_grarp2_150"] <- "P7_land use - grassland"
                names(IncNodePurity)[names(IncNodePurity)== "P_8_cl_grarp3_150"] <- "P8_land use - grassland second shearing"
                names(IncNodePurity)[names(IncNodePurity)== "P_9_cl_grarp4_150"] <- "P9land use - grassland special area"
                names(IncNodePurity)[names(IncNodePurity)== "P_10_cl_grarp5_150"] <- "P10_land use - crops"
                names(IncNodePurity)[names(IncNodePurity)== "P_11_dtm_150m"] <- "P11_dtm"
                names(IncNodePurity)[names(IncNodePurity)== "P_12_geo_bin_1_150"] <- "P12_geology - shist"
                names(IncNodePurity)[names(IncNodePurity)== "P_13_geo_bin_2_150"] <- "P13_geology - vulcanite"
                names(IncNodePurity)[names(IncNodePurity)== "P_14_geo_bin_3_150"] <- "P14_geology - vinschgau shear zone"
                names(IncNodePurity)[names(IncNodePurity)== "P_15_geo_bin_4_150"] <- "P15_geology - quaternary depositions"
                names(IncNodePurity)[names(IncNodePurity)== "P_16_geo_bin_5_150"] <- "P16_geology - sediment sequences"
                names(IncNodePurity)[names(IncNodePurity)== "P_17_geo_bin_6_150"] <- "P17_geology - plutons"
                names(IncNodePurity)[names(IncNodePurity)== "P_18_geo_bin_7_150"] <- "P18_geology - central gneiss"
                names(IncNodePurity)[names(IncNodePurity)== "P_19_geo_bin_8_150"] <- "P19_geology - quartz phyllite"
                names(IncNodePurity)[names(IncNodePurity)== "P_20_geo_bin_9_150"] <- "P20_geology - insignificant alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_21_geo_bi_10_150"] <- "P21_geology - medium-low alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_22_geo_bi_11_150"] <- "P22_geology - medium-high alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_23_slope_150m"] <- "P23_slope"
                names(IncNodePurity)[names(IncNodePurity)== "P_24_wetness_150m"] <- "P24_wetness_index"
                names(IncNodePurity)[names(IncNodePurity)== "P_25_gve_wgs84_150"] <- "P25_GVE"
                names(IncNodePurity)[names(IncNodePurity)== "P_26_geo_aeubin150"] <- "P26_Acidic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_27_geo_beubin150"] <- "P27_Basic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_28_geo_aelbin150"] <- "P28_Acidic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_29_geo_belbin150"] <- "P29_Basic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_30_geo_abelbi150"] <- "P30_Acidic-Basic soil"
              } else {}
              
              # Save PNG Table IncNodePurity
              incnodedff <- data.frame(IncNodePurity)
              names(incnodedff) <- ("IncNodePurity")
              table <- tableGrob(incnodedff)
              title <- textGrob(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia),gp=gpar(fontsize=8))
              footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
              padding <- unit(2,"line")
              table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
              table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
              table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
              png(file.path(jkp4,paste0("4_IncNodePurity_Table_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")), units="cm", width=(ncol(table))*15, height=(nrow(table)), res=300) #, height = 200, width = 100)
              grid.newpage() 
              grid.draw(table) 
              dev.off()
              # Save VALUES Table IncNodePurity
              write.csv(IncNodePurity,file.path(jkp4,paste0("4_IncNodePurity_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".csv")))      
              write.xlsx(IncNodePurity,file.path(jkp4,paste0("4_IncNodePurity_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".xlsx"))) 
              # Save Plot IncNodePurity
              x1<-incnodedff[,1]
              x<-rev(x1)
              y1<-names(IncNodePurity)
              y<-rev(y1)
              png(file.path(jkp4,paste0("4_IncNodePurity_Plot_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              dotchart(x,y,xlab="IncNodePurity") 
              title(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun))
              dev.off()
              
              # ImportanceSD (standard deviation)
              importancesdd<-modelRf@regModel$importanceSD
              importancesd<-sort(importancesdd,decreasing=TRUE)
              if (i==1) {
                names(importancesd)[names(importancesd)== "P_1_aspect_ea_150"] <- "P1__aspect east"
                names(importancesd)[names(importancesd)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(importancesd)[names(importancesd)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(importancesd)[names(importancesd)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(importancesd)[names(importancesd)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(importancesd)[names(importancesd)== "P_6_cl_grpre1_150"] <- "P6_land use - alpine meadow"
                names(importancesd)[names(importancesd)== "P_7_cl_grpre2_150"] <- "P7_land use - grassland"
                names(importancesd)[names(importancesd)== "P_8_cl_grpre3_150"] <- "P8_land use - grassland second shearing"
                names(importancesd)[names(importancesd)== "P_9_cl_grpre4_150"] <- "P9_land use - grassland special area"
                names(importancesd)[names(importancesd)== "P_10_dtm_150m"] <- "P10_dtm"
                names(importancesd)[names(importancesd)== "P_11_geo_bin_1_150"] <- "P11_geology - shist"
                names(importancesd)[names(importancesd)== "P_12_geo_bin_2_150"] <- "P12_geology - vulcanite"
                names(importancesd)[names(importancesd)== "P_13_geo_bin_3_150"] <- "P13_geology - vinschgau shear zone"
                names(importancesd)[names(importancesd)== "P_14_geo_bin_4_150"] <- "P14_geology - quaternary depositions"
                names(importancesd)[names(importancesd)== "P_15_geo_bin_5_150"] <- "P15_geology - sediment sequences"
                names(importancesd)[names(importancesd)== "P_16_geo_bin_6_150"] <- "P16_geology - plutons"
                names(importancesd)[names(importancesd)== "P_17_geo_bin_7_150"] <- "P17_geology - central gneiss"
                names(importancesd)[names(importancesd)== "P_18_geo_bin_8_150"] <- "P18_geology - quartz phyllite"
                names(importancesd)[names(importancesd)== "P_19_geo_bin_9_150"] <- "P19_geology - insignificant alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_20_geo_bi_10_150"] <- "P20_geology - medium-low alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_21_geo_bi_11_150"] <- "P21_geology - medium-high alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_22_slope_150m"] <- "P22_slope"
                names(importancesd)[names(importancesd)== "P_23_wetness_150m"] <- "P23_wetness_index"
                names(importancesd)[names(importancesd)== "P_24_gve_wgs84_150"] <- "P24_GVE"
                names(importancesd)[names(importancesd)== "P_25_geo_aeubin150"] <- "P25_Acidic soil"
                names(importancesd)[names(importancesd)== "P_26_geo_beubin150"] <- "P26_Basic soil"
                names(importancesd)[names(importancesd)== "P_27_geo_aelbin150"] <- "P27_Acidic soil"
                names(importancesd)[names(importancesd)== "P_28_geo_belbin150"] <- "P28_Basic soil"
                names(importancesd)[names(importancesd)== "P_29_geo_abelbi150"] <- "P29_Acidic-Basic soil"
              } else if (i==2) {
                names(importancesd)[names(importancesd)== "P_1_aspect_ea_150"] <- "P1_aspect east"
                names(importancesd)[names(importancesd)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(importancesd)[names(importancesd)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(importancesd)[names(importancesd)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(importancesd)[names(importancesd)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(importancesd)[names(importancesd)== "P_6_cl_grarp1_150"] <- "P6_land use - alpine meadow"
                names(importancesd)[names(importancesd)== "P_7_cl_grarp2_150"] <- "P7_land use - grassland"
                names(importancesd)[names(importancesd)== "P_8_cl_grarp3_150"] <- "P8_land use - grassland second shearing"
                names(importancesd)[names(importancesd)== "P_9_cl_grarp4_150"] <- "P9land use - grassland special area"
                names(importancesd)[names(importancesd)== "P_10_cl_grarp5_150"] <- "P10_land use - crops"
                names(importancesd)[names(importancesd)== "P_11_dtm_150m"] <- "P11_dtm"
                names(importancesd)[names(importancesd)== "P_12_geo_bin_1_150"] <- "P12_geology - shist"
                names(importancesd)[names(importancesd)== "P_13_geo_bin_2_150"] <- "P13_geology - vulcanite"
                names(importancesd)[names(importancesd)== "P_14_geo_bin_3_150"] <- "P14_geology - vinschgau shear zone"
                names(importancesd)[names(importancesd)== "P_15_geo_bin_4_150"] <- "P15_geology - quaternary depositions"
                names(importancesd)[names(importancesd)== "P_16_geo_bin_5_150"] <- "P16_geology - sediment sequences"
                names(importancesd)[names(importancesd)== "P_17_geo_bin_6_150"] <- "P17_geology - plutons"
                names(importancesd)[names(importancesd)== "P_18_geo_bin_7_150"] <- "P18_geology - central gneiss"
                names(importancesd)[names(importancesd)== "P_19_geo_bin_8_150"] <- "P19_geology - quartz phyllite"
                names(importancesd)[names(importancesd)== "P_20_geo_bin_9_150"] <- "P20_geology - insignificant alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_21_geo_bi_10_150"] <- "P21_geology - medium-low alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_22_geo_bi_11_150"] <- "P22_geology - medium-high alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_23_slope_150m"] <- "P23_slope"
                names(importancesd)[names(importancesd)== "P_24_wetness_150m"] <- "P24_wetness_index"
                names(importancesd)[names(importancesd)== "P_25_gve_wgs84_150"] <- "P25_GVE"
                names(importancesd)[names(importancesd)== "P_26_geo_aeubin150"] <- "P26_Acidic soil"
                names(importancesd)[names(importancesd)== "P_27_geo_beubin150"] <- "P27_Basic soil"
                names(importancesd)[names(importancesd)== "P_28_geo_aelbin150"] <- "P28_Acidic soil"
                names(importancesd)[names(importancesd)== "P_29_geo_belbin150"] <- "P29_Basic soil"
                names(importancesd)[names(importancesd)== "P_30_geo_abelbi150"] <- "P30_Acidic-Basic soil"
              } else {}
              
              # Save PNG Table importancesd
              importancesddff <- data.frame(importancesd)
              names(importancesddff) <- ("ImportanceSD")
              table <- tableGrob(importancesddff)
              title <- textGrob(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia),gp=gpar(fontsize=8))
              footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
              padding <- unit(2,"line")
              table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
              table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
              table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
              png(file.path(jkp4,paste0("4_ImportanceSD_Table_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")), units="cm", width=(ncol(table))*15, height=(nrow(table)), res=300) #, height = 200, width = 100)
              grid.newpage() 
              grid.draw(table) 
              dev.off()
              # Save VALUES Table importancesd
              write.csv(importancesd,file.path(jkp4,paste0("4_ImportanceSD_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".csv")))      
              write.xlsx(importancesd,file.path(jkp4,paste0("4_ImportanceSD_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".xlsx"))) 
              # Save Plot importancesd
              x1<-importancesddff[,1]
              x<-rev(x1)
              y1<-names(importancesd)
              y<-rev(y1)
              png(file.path(jkp4,paste0("4_ImportanceSD_Plot_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              dotchart(x,y,xlab="ImportanceSD") 
              title(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun))
              dev.off()
              
              # Mean of square residuals / % Var. explained
              printregmodel<-modelRf@regModel
              sink(file.path(jkp2,paste0("2_RfModel_Values",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".txt")))
              print(printregmodel)
              sink()
              setwd(paste0(out_dir,"/",predictareas_name[i],"/",var,"/",ai[k],"/","No.Run_",nrun,"/Random_Forest_Regression/2_RfModel"))
              armin<-read.delim(paste0("2_RfModel_Values",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".txt"),stringsAsFactors = FALSE)
              meanofsquareresiduals<-readr::parse_number(armin[5,])
              varexplained<-readr::parse_number(armin[6,])
              
              
              
              # Predictions ORDINARY KRIGING (with 80% of the data)
              predictionsRf<-GSIF::predict.gstatModel(object = modelRf, predictionLocations = grid, predict.method = "RK",verbose=T,nfold = 5)#vgmmodel = model@vgmModel,
              predictRf<-raster(predictionsRf@predicted)
              if (ai[k]=="raw") {
                
              } else if (ai[k]=="log") {
                predictRf_transf_log<-exp(predictRf)
              } else if (ai[k]=="sqrt") {
                predictRf_transf_sqrt<-quad(predictRf)
              } else if (ai[k]=="quad") {
                predictRf_transf_quad<-sqrt(predictRf)
              } else if (ai[k]=="inverse") {
                predictRf_transf_inverse<-inverse(predictRf)
              } else {}
              # Save Plots
              png(file.path(jkpp,paste0("5_Prediction_OrdinaryKriging_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              plot(predictionsRf) #paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              dev.off()
              
              
              
              # Plot Raster of predicted values
              png(file.path(jkppp,paste0("6_Predict_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              if (ai[k]=="raw") {
                raster::plot(predictRf,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else if (ai[k]=="log") {
                raster::plot(predictRf_transf_log,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else if (ai[k]=="sqrt") {
                raster::plot(predictRf_transf_sqrt,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else if (ai[k]=="quad") {
                raster::plot(predictRf_transf_quad,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else if (ai[k]=="inverse") {
                raster::plot(predictRf_transf_inverse,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else {}
              dev.off()
              
              
              
              
              # Overlay validation and prediction
              obsValidRf<-sp::over(obsValid,predictionsRf@predicted)
              if (ai[k]=="raw") {
                
              } else if (ai[k]=="log") {
                obsValidRf$original<-exp(obsValidRf$log)
              } else if (ai[k]=="sqrt") {
                obsValidRf$original<-quad(obsValidRf$sqrt)
              } else if (ai[k]=="quad") {
                obsValidRf$original<-sqrt(obsValidRf$quad)
              } else if (ai[k]=="inverse") {
                obsValidRf$original<-inverse(obsValidRf$inverse)
              } else {}
              
              
              
              # RMSE & NRMSE & R-SQUARED
              if (ai[k]=="raw") {
                rmseRf <- RMSE(pred = obsValidRf[,4], obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf[,4], obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf[,4], obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf[,4], obs = obsValid@data[[var]], na.rm = T)
              } else if (ai[k]=="log") {
                rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              } else if (ai[k]=="sqrt") {
                rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              } else if (ai[k]=="quad") {
                rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              } else if (ai[k]=="inverse") {
                rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              } else {}
              
              
              
              # Maps - Final rasters (GTiff)
              if (ai[k]=="raw") {
                writeRaster(x = predictRf, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else if (ai[k]=="log") {
                writeRaster(x = predictRf_transf_log, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else if (ai[k]=="sqrt") {
                writeRaster(x = predictRf_transf_sqrt, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else if (ai[k]=="quad") {
                writeRaster(x = predictRf_transf_quad, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else if (ai[k]=="inverse") {
                writeRaster(x = predictRf_transf_inverse,filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else {}
              
              # spalte 1 - No.Run
              sdf1<-length(spalte1)
              spalte1<-append(spalte1,nrun,after=sdf1)
              nrun<-nrun+1
              # Spalte 2 - Transformation
              sdf2<-length(spalte2)
              spalte2<-append(spalte2,ai[k],after=sdf2)
              # Spalte 3 - Cutoff
              sdf3<-length(spalte3)
              spalte3<-append(spalte3,taia,after=sdf3)
              # Spalte 4 - Predictors
              sdf4<-length(spalte4)
              spalte4<-append(spalte4,paste(get(paste0("predictors_combi_",l)), collapse = "/"))
              # Spalte 5 - Mean of squared residuals
              sdf5<-length(spalte5)
              spalte5<-append(spalte5,meanofsquareresiduals,after=sdf5)
              # Spalte 6 - % Var explained
              sdf6<-length(spalte6)
              spalte6<-append(spalte6,varexplained,after=sdf6)
              # Spalte 7 - RMSE
              sdf7<-length(spalte7)
              spalte7<-append(spalte7,rmseRf,after=sdf7)
              sdf7tot<-length(spalte7tot)
              spalte7tot<-append(spalte7tot,rmseRf,after=sdf7tot)
              # Spalte 8 - NRMSE
              sdf8<-length(spalte8)
              spalte8<-append(spalte8,nrmseRf,after=sdf8)
              # Spalte 9 - NRMSEMAXMIN
              sdf9<-length(spalte9)
              spalte9<-append(spalte9,nrmseMaxMinRf,after=sdf9)
              # Spalte 10 - R-SQUARED
              sdf10<-length(spalte10)
              spalte10<-append(spalte10,r2Rf,after=sdf10)
              sdf10tot<-length(spalte10tot)
              spalte10tot<-append(spalte10tot,r2Rf,after=sdf10tot)
              
              
            } else if (s==2 | s==3) {
              print("length(cutoffs)=2 / s==2 | s==3")
              maxcutoff<-max(cutoffs)
              try(modelRf<-fit.gstatModel(observations = obs,formulaString = formula,covariates = grid,method="randomForest",cutoff=maxcutoff))
              taia<-maxcutoff
              # MAIN CODE
              # Save RfModel & VgmModel
              with(faithful,plot(modelRf)) 
              dev.copy(png,file.path(jkp2,paste0("2_RfModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".png")),res=300,width =3500, height=1750)
              dev.off()
              vgmModel = modelRf@vgmModel
              print(vgmModel)
              write.csv(vgmModel,file.path(jkp3,paste0("3_VgmModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".csv")))    
              write.xlsx(vgmModel,file.path(jkp3,paste0("3_VgmModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".xlsx")))    
              sink(file.path(jkp3,paste0("3_VgmModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".txt"))) 
              print(vgmModel)
              sink()
              
              # Save Importance of the predictors (%IncMse & IncNodePurity)
              # PredictorsImportance
              importance = modelRf@regModel$importance
              # %IncMse
              incmse<-sort(importance[,1],decreasing = TRUE)
              # Define names for predictors
              if (i==1) {
                names(incmse)[names(incmse)== "P_1_aspect_ea_150"] <- "P1__aspect east"
                names(incmse)[names(incmse)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(incmse)[names(incmse)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(incmse)[names(incmse)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(incmse)[names(incmse)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(incmse)[names(incmse)== "P_6_cl_grpre1_150"] <- "P6_land use - alpine meadow"
                names(incmse)[names(incmse)== "P_7_cl_grpre2_150"] <- "P7_land use - grassland"
                names(incmse)[names(incmse)== "P_8_cl_grpre3_150"] <- "P8_land use - grassland second shearing"
                names(incmse)[names(incmse)== "P_9_cl_grpre4_150"] <- "P9_land use - grassland special area"
                names(incmse)[names(incmse)== "P_10_dtm_150m"] <- "P10_dtm"
                names(incmse)[names(incmse)== "P_11_geo_bin_1_150"] <- "P11_geology - shist"
                names(incmse)[names(incmse)== "P_12_geo_bin_2_150"] <- "P12_geology - vulcanite"
                names(incmse)[names(incmse)== "P_13_geo_bin_3_150"] <- "P13_geology - vinschgau shear zone"
                names(incmse)[names(incmse)== "P_14_geo_bin_4_150"] <- "P14_geology - quaternary depositions"
                names(incmse)[names(incmse)== "P_15_geo_bin_5_150"] <- "P15_geology - sediment sequences"
                names(incmse)[names(incmse)== "P_16_geo_bin_6_150"] <- "P16_geology - plutons"
                names(incmse)[names(incmse)== "P_17_geo_bin_7_150"] <- "P17_geology - central gneiss"
                names(incmse)[names(incmse)== "P_18_geo_bin_8_150"] <- "P18_geology - quartz phyllite"
                names(incmse)[names(incmse)== "P_19_geo_bin_9_150"] <- "P19_geology - insignificant alpine metamorphosis"
                names(incmse)[names(incmse)== "P_20_geo_bi_10_150"] <- "P20_geology - medium-low alpine metamorphosis"
                names(incmse)[names(incmse)== "P_21_geo_bi_11_150"] <- "P21_geology - medium-high alpine metamorphosis"
                names(incmse)[names(incmse)== "P_22_slope_150m"] <- "P22_slope"
                names(incmse)[names(incmse)== "P_23_wetness_150m"] <- "P23_wetness_index"
                names(incmse)[names(incmse)== "P_24_gve_wgs84_150"] <- "P24_GVE"
                names(incmse)[names(incmse)== "P_25_geo_aeubin150"] <- "P25_Acidic soil"
                names(incmse)[names(incmse)== "P_26_geo_beubin150"] <- "P26_Basic soil"
                names(incmse)[names(incmse)== "P_27_geo_aelbin150"] <- "P27_Acidic soil"
                names(incmse)[names(incmse)== "P_28_geo_belbin150"] <- "P28_Basic soil"
                names(incmse)[names(incmse)== "P_29_geo_abelbi150"] <- "P29_Acidic-Basic soil"
              } else if (i==2) {
                names(incmse)[names(incmse)== "P_1_aspect_ea_150"] <- "P1_aspect east"
                names(incmse)[names(incmse)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(incmse)[names(incmse)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(incmse)[names(incmse)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(incmse)[names(incmse)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(incmse)[names(incmse)== "P_6_cl_grarp1_150"] <- "P6_land use - alpine meadow"
                names(incmse)[names(incmse)== "P_7_cl_grarp2_150"] <- "P7_land use - grassland"
                names(incmse)[names(incmse)== "P_8_cl_grarp3_150"] <- "P8_land use - grassland second shearing"
                names(incmse)[names(incmse)== "P_9_cl_grarp4_150"] <- "P9_land use - grassland special area"
                names(incmse)[names(incmse)== "P_10_cl_grarp5_150"] <- "P10_land use - crops"
                names(incmse)[names(incmse)== "P_11_dtm_150m"] <- "P11_dtm"
                names(incmse)[names(incmse)== "P_12_geo_bin_1_150"] <- "P12_geology - shist"
                names(incmse)[names(incmse)== "P_13_geo_bin_2_150"] <- "P13_geology - vulcanite"
                names(incmse)[names(incmse)== "P_14_geo_bin_3_150"] <- "P14_geology - vinschgau shear zone"
                names(incmse)[names(incmse)== "P_15_geo_bin_4_150"] <- "P15_geology - quaternary depositions"
                names(incmse)[names(incmse)== "P_16_geo_bin_5_150"] <- "P16_geology - sediment sequences"
                names(incmse)[names(incmse)== "P_17_geo_bin_6_150"] <- "P17_geology - plutons"
                names(incmse)[names(incmse)== "P_18_geo_bin_7_150"] <- "P18_geology - central gneiss"
                names(incmse)[names(incmse)== "P_19_geo_bin_8_150"] <- "P19_geology - quartz phyllite"
                names(incmse)[names(incmse)== "P_20_geo_bin_9_150"] <- "P20_geology - insignificant alpine metamorphosis"
                names(incmse)[names(incmse)== "P_21_geo_bi_10_150"] <- "P21_geology - medium-low alpine metamorphosis"
                names(incmse)[names(incmse)== "P_22_geo_bi_11_150"] <- "P22_geology - medium-high alpine metamorphosis"
                names(incmse)[names(incmse)== "P_23_slope_150m"] <- "P23_slope"
                names(incmse)[names(incmse)== "P_24_wetness_150m"] <- "P24_wetness_index"
                names(incmse)[names(incmse)== "P_25_gve_wgs84_150"] <- "P25_GVE"
                names(incmse)[names(incmse)== "P_26_geo_aeubin150"] <- "P26_Acidic soil"
                names(incmse)[names(incmse)== "P_27_geo_beubin150"] <- "P27_Basic soil"
                names(incmse)[names(incmse)== "P_28_geo_aelbin150"] <- "P28_Acidic soil"
                names(incmse)[names(incmse)== "P_29_geo_belbin150"] <- "P29_Basic soil"
                names(incmse)[names(incmse)== "P_30_geo_abelbi150"] <- "P30_Acidic-Basic soil"
              } else {}
              
              # Save PNG Table %IncMSE
              incmsedff <- data.frame(incmse)
              names(incmsedff) <- ("%IncMSE")
              table <- tableGrob(incmsedff)
              title <- textGrob(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia),gp=gpar(fontsize=8))
              footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
              padding <- unit(2,"line")
              table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
              table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
              table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
              png(file.path(jkp4,paste0("4_IncMse_Table_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")), units="cm", width=(ncol(table))*15, height=(nrow(table)), res=300) #, height = 200, width = 100)
              grid.newpage() 
              grid.draw(table) 
              dev.off()
              # Save VALUES Table %IncMSE
              write.csv(incmse,file.path(jkp4,paste0("4_IncMse_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".csv")))      
              write.xlsx(incmse,file.path(jkp4,paste0("4_IncMse_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".xlsx")))    
              # Save Plot %IncMSE
              x1<-incmsedff[,1]
              x<-rev(x1)
              y1<-names(incmse)
              y<-rev(y1)
              png(file.path(jkp4,paste0("4_IncMse_Plot_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              dotchart(x,y,xlab="%IncMse") 
              title(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun))
              dev.off()
              
              # IncNodePurity
              IncNodePurity<-round(sort(importance[,2],decreasing = TRUE),digits = 2) 
              # Define names for predictors
              if (i==1) {
                names(IncNodePurity)[names(IncNodePurity)== "P_1_aspect_ea_150"] <- "P1__aspect east"
                names(IncNodePurity)[names(IncNodePurity)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(IncNodePurity)[names(IncNodePurity)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(IncNodePurity)[names(IncNodePurity)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(IncNodePurity)[names(IncNodePurity)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(IncNodePurity)[names(IncNodePurity)== "P_6_cl_grpre1_150"] <- "P6_land use - alpine meadow"
                names(IncNodePurity)[names(IncNodePurity)== "P_7_cl_grpre2_150"] <- "P7_land use - grassland"
                names(IncNodePurity)[names(IncNodePurity)== "P_8_cl_grpre3_150"] <- "P8_land use - grassland second shearing"
                names(IncNodePurity)[names(IncNodePurity)== "P_9_cl_grpre4_150"] <- "P9_land use - grassland special area"
                names(IncNodePurity)[names(IncNodePurity)== "P_10_dtm_150m"] <- "P10_dtm"
                names(IncNodePurity)[names(IncNodePurity)== "P_11_geo_bin_1_150"] <- "P11_geology - shist"
                names(IncNodePurity)[names(IncNodePurity)== "P_12_geo_bin_2_150"] <- "P12_geology - vulcanite"
                names(IncNodePurity)[names(IncNodePurity)== "P_13_geo_bin_3_150"] <- "P13_geology - vinschgau shear zone"
                names(IncNodePurity)[names(IncNodePurity)== "P_14_geo_bin_4_150"] <- "P14_geology - quaternary depositions"
                names(IncNodePurity)[names(IncNodePurity)== "P_15_geo_bin_5_150"] <- "P15_geology - sediment sequences"
                names(IncNodePurity)[names(IncNodePurity)== "P_16_geo_bin_6_150"] <- "P16_geology - plutons"
                names(IncNodePurity)[names(IncNodePurity)== "P_17_geo_bin_7_150"] <- "P17_geology - central gneiss"
                names(IncNodePurity)[names(IncNodePurity)== "P_18_geo_bin_8_150"] <- "P18_geology - quartz phyllite"
                names(IncNodePurity)[names(IncNodePurity)== "P_19_geo_bin_9_150"] <- "P19_geology - insignificant alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_20_geo_bi_10_150"] <- "P20_geology - medium-low alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_21_geo_bi_11_150"] <- "P21_geology - medium-high alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_22_slope_150m"] <- "P22_slope"
                names(IncNodePurity)[names(IncNodePurity)== "P_23_wetness_150m"] <- "P23_wetness_index"
                names(IncNodePurity)[names(IncNodePurity)== "P_24_gve_wgs84_150"] <- "P24_GVE"
                names(IncNodePurity)[names(IncNodePurity)== "P_25_geo_aeubin150"] <- "P25_Acidic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_26_geo_beubin150"] <- "P26_Basic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_27_geo_aelbin150"] <- "P27_Acidic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_28_geo_belbin150"] <- "P28_Basic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_29_geo_abelbi150"] <- "P29_Acidic-Basic soil"
              } else if (i==2) {
                names(IncNodePurity)[names(IncNodePurity)== "P_1_aspect_ea_150"] <- "P1_aspect east"
                names(IncNodePurity)[names(IncNodePurity)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(IncNodePurity)[names(IncNodePurity)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(IncNodePurity)[names(IncNodePurity)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(IncNodePurity)[names(IncNodePurity)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(IncNodePurity)[names(IncNodePurity)== "P_6_cl_grarp1_150"] <- "P6_land use - alpine meadow"
                names(IncNodePurity)[names(IncNodePurity)== "P_7_cl_grarp2_150"] <- "P7_land use - grassland"
                names(IncNodePurity)[names(IncNodePurity)== "P_8_cl_grarp3_150"] <- "P8_land use - grassland second shearing"
                names(IncNodePurity)[names(IncNodePurity)== "P_9_cl_grarp4_150"] <- "P9land use - grassland special area"
                names(IncNodePurity)[names(IncNodePurity)== "P_10_cl_grarp5_150"] <- "P10_land use - crops"
                names(IncNodePurity)[names(IncNodePurity)== "P_11_dtm_150m"] <- "P11_dtm"
                names(IncNodePurity)[names(IncNodePurity)== "P_12_geo_bin_1_150"] <- "P12_geology - shist"
                names(IncNodePurity)[names(IncNodePurity)== "P_13_geo_bin_2_150"] <- "P13_geology - vulcanite"
                names(IncNodePurity)[names(IncNodePurity)== "P_14_geo_bin_3_150"] <- "P14_geology - vinschgau shear zone"
                names(IncNodePurity)[names(IncNodePurity)== "P_15_geo_bin_4_150"] <- "P15_geology - quaternary depositions"
                names(IncNodePurity)[names(IncNodePurity)== "P_16_geo_bin_5_150"] <- "P16_geology - sediment sequences"
                names(IncNodePurity)[names(IncNodePurity)== "P_17_geo_bin_6_150"] <- "P17_geology - plutons"
                names(IncNodePurity)[names(IncNodePurity)== "P_18_geo_bin_7_150"] <- "P18_geology - central gneiss"
                names(IncNodePurity)[names(IncNodePurity)== "P_19_geo_bin_8_150"] <- "P19_geology - quartz phyllite"
                names(IncNodePurity)[names(IncNodePurity)== "P_20_geo_bin_9_150"] <- "P20_geology - insignificant alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_21_geo_bi_10_150"] <- "P21_geology - medium-low alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_22_geo_bi_11_150"] <- "P22_geology - medium-high alpine metamorphosis"
                names(IncNodePurity)[names(IncNodePurity)== "P_23_slope_150m"] <- "P23_slope"
                names(IncNodePurity)[names(IncNodePurity)== "P_24_wetness_150m"] <- "P24_wetness_index"
                names(IncNodePurity)[names(IncNodePurity)== "P_25_gve_wgs84_150"] <- "P25_GVE"
                names(IncNodePurity)[names(IncNodePurity)== "P_26_geo_aeubin150"] <- "P26_Acidic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_27_geo_beubin150"] <- "P27_Basic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_28_geo_aelbin150"] <- "P28_Acidic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_29_geo_belbin150"] <- "P29_Basic soil"
                names(IncNodePurity)[names(IncNodePurity)== "P_30_geo_abelbi150"] <- "P30_Acidic-Basic soil"
              } else {}
              
              # Save PNG Table IncNodePurity
              incnodedff <- data.frame(IncNodePurity)
              names(incnodedff) <- ("IncNodePurity")
              table <- tableGrob(incnodedff)
              title <- textGrob(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia),gp=gpar(fontsize=8))
              footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
              padding <- unit(2,"line")
              table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
              table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
              table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
              png(file.path(jkp4,paste0("4_IncNodePurity_Table_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")), units="cm", width=(ncol(table))*15, height=(nrow(table)), res=300) #, height = 200, width = 100)
              grid.newpage() 
              grid.draw(table) 
              dev.off()
              # Save VALUES Table IncNodePurity
              write.csv(IncNodePurity,file.path(jkp4,paste0("4_IncNodePurity_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".csv")))      
              write.xlsx(IncNodePurity,file.path(jkp4,paste0("4_IncNodePurity_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".xlsx"))) 
              # Save Plot IncNodePurity
              x1<-incnodedff[,1]
              x<-rev(x1)
              y1<-names(IncNodePurity)
              y<-rev(y1)
              png(file.path(jkp4,paste0("4_IncNodePurity_Plot_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              dotchart(x,y,xlab="IncNodePurity") 
              title(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun))
              dev.off()
              
              # ImportanceSD (standard deviation)
              importancesdd<-modelRf@regModel$importanceSD
              importancesd<-sort(importancesdd,decreasing=TRUE)
              if (i==1) {
                names(importancesd)[names(importancesd)== "P_1_aspect_ea_150"] <- "P1__aspect east"
                names(importancesd)[names(importancesd)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(importancesd)[names(importancesd)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(importancesd)[names(importancesd)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(importancesd)[names(importancesd)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(importancesd)[names(importancesd)== "P_6_cl_grpre1_150"] <- "P6_land use - alpine meadow"
                names(importancesd)[names(importancesd)== "P_7_cl_grpre2_150"] <- "P7_land use - grassland"
                names(importancesd)[names(importancesd)== "P_8_cl_grpre3_150"] <- "P8_land use - grassland second shearing"
                names(importancesd)[names(importancesd)== "P_9_cl_grpre4_150"] <- "P9_land use - grassland special area"
                names(importancesd)[names(importancesd)== "P_10_dtm_150m"] <- "P10_dtm"
                names(importancesd)[names(importancesd)== "P_11_geo_bin_1_150"] <- "P11_geology - shist"
                names(importancesd)[names(importancesd)== "P_12_geo_bin_2_150"] <- "P12_geology - vulcanite"
                names(importancesd)[names(importancesd)== "P_13_geo_bin_3_150"] <- "P13_geology - vinschgau shear zone"
                names(importancesd)[names(importancesd)== "P_14_geo_bin_4_150"] <- "P14_geology - quaternary depositions"
                names(importancesd)[names(importancesd)== "P_15_geo_bin_5_150"] <- "P15_geology - sediment sequences"
                names(importancesd)[names(importancesd)== "P_16_geo_bin_6_150"] <- "P16_geology - plutons"
                names(importancesd)[names(importancesd)== "P_17_geo_bin_7_150"] <- "P17_geology - central gneiss"
                names(importancesd)[names(importancesd)== "P_18_geo_bin_8_150"] <- "P18_geology - quartz phyllite"
                names(importancesd)[names(importancesd)== "P_19_geo_bin_9_150"] <- "P19_geology - insignificant alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_20_geo_bi_10_150"] <- "P20_geology - medium-low alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_21_geo_bi_11_150"] <- "P21_geology - medium-high alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_22_slope_150m"] <- "P22_slope"
                names(importancesd)[names(importancesd)== "P_23_wetness_150m"] <- "P23_wetness_index"
                names(importancesd)[names(importancesd)== "P_24_gve_wgs84_150"] <- "P24_GVE"
                names(importancesd)[names(importancesd)== "P_25_geo_aeubin150"] <- "P25_Acidic soil"
                names(importancesd)[names(importancesd)== "P_26_geo_beubin150"] <- "P26_Basic soil"
                names(importancesd)[names(importancesd)== "P_27_geo_aelbin150"] <- "P27_Acidic soil"
                names(importancesd)[names(importancesd)== "P_28_geo_belbin150"] <- "P28_Basic soil"
                names(importancesd)[names(importancesd)== "P_29_geo_abelbi150"] <- "P29_Acidic-Basic soil"
              } else if (i==2) {
                names(importancesd)[names(importancesd)== "P_1_aspect_ea_150"] <- "P1_aspect east"
                names(importancesd)[names(importancesd)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
                names(importancesd)[names(importancesd)== "P_3_aspect_no_150"] <- "P3_aspect north"
                names(importancesd)[names(importancesd)== "P_4_aspect_so_150"] <- "P4_aspect south"
                names(importancesd)[names(importancesd)== "P_5_aspect_we_150"] <- "P5_aspect west"
                names(importancesd)[names(importancesd)== "P_6_cl_grarp1_150"] <- "P6_land use - alpine meadow"
                names(importancesd)[names(importancesd)== "P_7_cl_grarp2_150"] <- "P7_land use - grassland"
                names(importancesd)[names(importancesd)== "P_8_cl_grarp3_150"] <- "P8_land use - grassland second shearing"
                names(importancesd)[names(importancesd)== "P_9_cl_grarp4_150"] <- "P9land use - grassland special area"
                names(importancesd)[names(importancesd)== "P_10_cl_grarp5_150"] <- "P10_land use - crops"
                names(importancesd)[names(importancesd)== "P_11_dtm_150m"] <- "P11_dtm"
                names(importancesd)[names(importancesd)== "P_12_geo_bin_1_150"] <- "P12_geology - shist"
                names(importancesd)[names(importancesd)== "P_13_geo_bin_2_150"] <- "P13_geology - vulcanite"
                names(importancesd)[names(importancesd)== "P_14_geo_bin_3_150"] <- "P14_geology - vinschgau shear zone"
                names(importancesd)[names(importancesd)== "P_15_geo_bin_4_150"] <- "P15_geology - quaternary depositions"
                names(importancesd)[names(importancesd)== "P_16_geo_bin_5_150"] <- "P16_geology - sediment sequences"
                names(importancesd)[names(importancesd)== "P_17_geo_bin_6_150"] <- "P17_geology - plutons"
                names(importancesd)[names(importancesd)== "P_18_geo_bin_7_150"] <- "P18_geology - central gneiss"
                names(importancesd)[names(importancesd)== "P_19_geo_bin_8_150"] <- "P19_geology - quartz phyllite"
                names(importancesd)[names(importancesd)== "P_20_geo_bin_9_150"] <- "P20_geology - insignificant alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_21_geo_bi_10_150"] <- "P21_geology - medium-low alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_22_geo_bi_11_150"] <- "P22_geology - medium-high alpine metamorphosis"
                names(importancesd)[names(importancesd)== "P_23_slope_150m"] <- "P23_slope"
                names(importancesd)[names(importancesd)== "P_24_wetness_150m"] <- "P24_wetness_index"
                names(importancesd)[names(importancesd)== "P_25_gve_wgs84_150"] <- "P25_GVE"
                names(importancesd)[names(importancesd)== "P_26_geo_aeubin150"] <- "P26_Acidic soil"
                names(importancesd)[names(importancesd)== "P_27_geo_beubin150"] <- "P27_Basic soil"
                names(importancesd)[names(importancesd)== "P_28_geo_aelbin150"] <- "P28_Acidic soil"
                names(importancesd)[names(importancesd)== "P_29_geo_belbin150"] <- "P29_Basic soil"
                names(importancesd)[names(importancesd)== "P_30_geo_abelbi150"] <- "P30_Acidic-Basic soil"
              } else {}
              
              # Save PNG Table importancesd
              importancesddff <- data.frame(importancesd)
              names(importancesddff) <- ("ImportanceSD")
              table <- tableGrob(importancesddff)
              title <- textGrob(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia),gp=gpar(fontsize=8))
              footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
              padding <- unit(2,"line")
              table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
              table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
              table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
              png(file.path(jkp4,paste0("4_ImportanceSD_Table_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")), units="cm", width=(ncol(table))*15, height=(nrow(table)), res=300) #, height = 200, width = 100)
              grid.newpage() 
              grid.draw(table) 
              dev.off()
              # Save VALUES Table importancesd
              write.csv(importancesd,file.path(jkp4,paste0("4_ImportanceSD_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".csv")))      
              write.xlsx(importancesd,file.path(jkp4,paste0("4_ImportanceSD_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".xlsx"))) 
              # Save Plot importancesd
              x1<-importancesddff[,1]
              x<-rev(x1)
              y1<-names(importancesd)
              y<-rev(y1)
              png(file.path(jkp4,paste0("4_ImportanceSD_Plot_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              dotchart(x,y,xlab="ImportanceSD") 
              title(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun))
              dev.off()
              
              # Mean of square residuals / % Var. explained
              printregmodel<-modelRf@regModel
              sink(file.path(jkp2,paste0("2_RfModel_Values",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".txt")))
              print(printregmodel)
              sink()
              setwd(paste0(out_dir,"/",predictareas_name[i],"/",var,"/",ai[k],"/","No.Run_",nrun,"/Random_Forest_Regression/2_RfModel"))
              armin<-read.delim(paste0("2_RfModel_Values",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".txt"),stringsAsFactors = FALSE)
              meanofsquareresiduals<-readr::parse_number(armin[5,])
              varexplained<-readr::parse_number(armin[6,])
              
              
              
              # Predictions ORDINARY KRIGING (with 80% of the data)
              predictionsRf<-GSIF::predict.gstatModel(object = modelRf, predictionLocations = grid, predict.method = "RK",verbose=T,nfold = 5)#vgmmodel = model@vgmModel,
              predictRf<-raster(predictionsRf@predicted)
              if (ai[k]=="raw") {
                
              } else if (ai[k]=="log") {
                predictRf_transf_log<-exp(predictRf)
              } else if (ai[k]=="sqrt") {
                predictRf_transf_sqrt<-quad(predictRf)
              } else if (ai[k]=="quad") {
                predictRf_transf_quad<-sqrt(predictRf)
              } else if (ai[k]=="inverse") {
                predictRf_transf_inverse<-inverse(predictRf)
              } else {}
              # Save Plots
              png(file.path(jkpp,paste0("5_Prediction_OrdinaryKriging_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              plot(predictionsRf) #paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              dev.off()
              
              
              
              # Plot Raster of predicted values
              png(file.path(jkppp,paste0("6_Predict_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
              if (ai[k]=="raw") {
                raster::plot(predictRf,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else if (ai[k]=="log") {
                raster::plot(predictRf_transf_log,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else if (ai[k]=="sqrt") {
                raster::plot(predictRf_transf_sqrt,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else if (ai[k]=="quad") {
                raster::plot(predictRf_transf_quad,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else if (ai[k]=="inverse") {
                raster::plot(predictRf_transf_inverse,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
              } else {}
              dev.off()
              
              
              
              
              # Overlay validation and prediction
              obsValidRf<-sp::over(obsValid,predictionsRf@predicted)
              if (ai[k]=="raw") {
                
              } else if (ai[k]=="log") {
                obsValidRf$original<-exp(obsValidRf$log)
              } else if (ai[k]=="sqrt") {
                obsValidRf$original<-quad(obsValidRf$sqrt)
              } else if (ai[k]=="quad") {
                obsValidRf$original<-sqrt(obsValidRf$quad)
              } else if (ai[k]=="inverse") {
                obsValidRf$original<-inverse(obsValidRf$inverse)
              } else {}
              
              
              
              # RMSE & NRMSE & R-SQUARED
              if (ai[k]=="raw") {
                rmseRf <- RMSE(pred = obsValidRf[,4], obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf[,4], obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf[,4], obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf[,4], obs = obsValid@data[[var]], na.rm = T)
              } else if (ai[k]=="log") {
                rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              } else if (ai[k]=="sqrt") {
                rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              } else if (ai[k]=="quad") {
                rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              } else if (ai[k]=="inverse") {
                rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
                nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
                r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              } else {}
              
              
              
              # Maps - Final rasters (GTiff)
              if (ai[k]=="raw") {
                writeRaster(x = predictRf, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else if (ai[k]=="log") {
                writeRaster(x = predictRf_transf_log, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else if (ai[k]=="sqrt") {
                writeRaster(x = predictRf_transf_sqrt, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else if (ai[k]=="quad") {
                writeRaster(x = predictRf_transf_quad, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else if (ai[k]=="inverse") {
                writeRaster(x = predictRf_transf_inverse,filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                            overwrite = TRUE, format = "GTiff")
                save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                     list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                     file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
                
              } else {}
              
              # spalte 1 - No.Run
              sdf1<-length(spalte1)
              spalte1<-append(spalte1,nrun,after=sdf1)
              nrun<-nrun+1
              # Spalte 2 - Transformation
              sdf2<-length(spalte2)
              spalte2<-append(spalte2,ai[k],after=sdf2)
              # Spalte 3 - Cutoff
              sdf3<-length(spalte3)
              spalte3<-append(spalte3,taia,after=sdf3)
              # Spalte 4 - Predictors
              sdf4<-length(spalte4)
              spalte4<-append(spalte4,paste(get(paste0("predictors_combi_",l)), collapse = "/"))
              # Spalte 5 - Mean of squared residuals
              sdf5<-length(spalte5)
              spalte5<-append(spalte5,meanofsquareresiduals,after=sdf5)
              # Spalte 6 - % Var explained
              sdf6<-length(spalte6)
              spalte6<-append(spalte6,varexplained,after=sdf6)
              # Spalte 7 - RMSE
              sdf7<-length(spalte7)
              spalte7<-append(spalte7,rmseRf,after=sdf7)
              sdf7tot<-length(spalte7tot)
              spalte7tot<-append(spalte7tot,rmseRf,after=sdf7tot)
              # Spalte 8 - NRMSE
              sdf8<-length(spalte8)
              spalte8<-append(spalte8,nrmseRf,after=sdf8)
              # Spalte 9 - NRMSEMAXMIN
              sdf9<-length(spalte9)
              spalte9<-append(spalte9,nrmseMaxMinRf,after=sdf9)
              # Spalte 10 - R-SQUARED
              sdf10<-length(spalte10)
              spalte10<-append(spalte10,r2Rf,after=sdf10)
              sdf10tot<-length(spalte10tot)
              spalte10tot<-append(spalte10tot,r2Rf,after=sdf10tot)
              
            } else {}
            
            
          } else if (length(cutoffs)==1) {
            print("length(cutoffs)=1")
            # Output folders
            fg<-paste0("/output/",predictareas_name[i],"/",variables[j],"/",ai[k])
            #fh<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),fg)
            fh<-paste0(script_dir,fg)
            setwd(fh)
            dir.create(paste0("No.Run_",nrun))
            fgg<-paste0("/output/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun)
            #fhh<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),fgg)
            fhh<-paste0(script_dir,fgg)
            setwd(fhh)
            dir.create("Random_Forest_Regression")
            jkp<-paste0(out_dir,"/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun,"/Random_Forest_Regression")
            setwd(jkp)
            dir.create("1_Variogram")
            jkp1<-paste0(out_dir,"/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun,"/Random_Forest_Regression/1_Variogram")
            dir.create("2_RfModel")
            jkp2<-paste0(out_dir,"/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun,"/Random_Forest_Regression/2_RfModel")
            dir.create("3_VgmModel")
            jkp3<-paste0(out_dir,"/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun,"/Random_Forest_Regression/3_VgmModel")
            dir.create("4_PredictorsImportance")
            jkp4<-paste0(out_dir,"/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun,"/Random_Forest_Regression/4_PredictorsImportance")
            setwd(fhh)
            dir.create("Ordinary_Kriging_Prediction")
            jkpp<-paste0(out_dir,"/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun,"/Ordinary_Kriging_Prediction")
            dir.create("Final_Models")
            jkppp<-paste0(out_dir,"/",predictareas_name[i],"/",variables[j],"/",ai[k],"/","No.Run_",nrun,"/Final_Models")
            
            # Variogram
            #plot(variogram(formula_variabile, obs), plot.nu=T, pch="+")
            x=variogram(formula_variabile, obs)
            plott = ggplot(x,aes(dist,gamma))+
              geom_point()+
              theme_minimal()+
              ggtitle(paste0("Variogram - ","[",variables[j],"] ","- ",ai[k],"-","predictarea",i))+
              theme(plot.title = element_text(hjust = 0.5))
            jkpt<-paste0("Variogram - ","[",variables[j],"] ","- ",ai[k],"-","predictarea",i)
            #png(file.path(jkp1,paste0("1_Variogram1_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =4000, height=4000)
            #plot(variogram(formula_variabile, obs), plot.nu=T, pch="+",main=jkpt)
            #dev.off()
            ggsave(filename = file.path(jkp1,paste0("1_Variogram_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),plot = plott,device = "png",width=5, height=5)
            
            mincutoff<-min(cutoffs)
            try(modelRf<-fit.gstatModel(observations = obs,formulaString = formula,covariates = grid,method="randomForest",cutoff=mincutoff))
            taia<-mincutoff
            # MAIN CODE
            # Save RfModel & VgmModel
            with(faithful,plot(modelRf)) 
            dev.copy(png,file.path(jkp2,paste0("2_RfModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".png")),res=300,width =3500, height=1750)
            dev.off()
            vgmModel = modelRf@vgmModel
            print(vgmModel)
            write.csv(vgmModel,file.path(jkp3,paste0("3_VgmModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".csv")))    
            write.xlsx(vgmModel,file.path(jkp3,paste0("3_VgmModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".xlsx")))    
            sink(file.path(jkp3,paste0("3_VgmModel_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".txt"))) 
            print(vgmModel)
            sink()
            
            # Save Importance of the predictors (%IncMse & IncNodePurity)
            # PredictorsImportance
            importance = modelRf@regModel$importance
            # %IncMse
            incmse<-sort(importance[,1],decreasing = TRUE)
            # Define names for predictors
            if (i==1) {
              names(incmse)[names(incmse)== "P_1_aspect_ea_150"] <- "P1__aspect east"
              names(incmse)[names(incmse)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
              names(incmse)[names(incmse)== "P_3_aspect_no_150"] <- "P3_aspect north"
              names(incmse)[names(incmse)== "P_4_aspect_so_150"] <- "P4_aspect south"
              names(incmse)[names(incmse)== "P_5_aspect_we_150"] <- "P5_aspect west"
              names(incmse)[names(incmse)== "P_6_cl_grpre1_150"] <- "P6_land use - alpine meadow"
              names(incmse)[names(incmse)== "P_7_cl_grpre2_150"] <- "P7_land use - grassland"
              names(incmse)[names(incmse)== "P_8_cl_grpre3_150"] <- "P8_land use - grassland second shearing"
              names(incmse)[names(incmse)== "P_9_cl_grpre4_150"] <- "P9_land use - grassland special area"
              names(incmse)[names(incmse)== "P_10_dtm_150m"] <- "P10_dtm"
              names(incmse)[names(incmse)== "P_11_geo_bin_1_150"] <- "P11_geology - shist"
              names(incmse)[names(incmse)== "P_12_geo_bin_2_150"] <- "P12_geology - vulcanite"
              names(incmse)[names(incmse)== "P_13_geo_bin_3_150"] <- "P13_geology - vinschgau shear zone"
              names(incmse)[names(incmse)== "P_14_geo_bin_4_150"] <- "P14_geology - quaternary depositions"
              names(incmse)[names(incmse)== "P_15_geo_bin_5_150"] <- "P15_geology - sediment sequences"
              names(incmse)[names(incmse)== "P_16_geo_bin_6_150"] <- "P16_geology - plutons"
              names(incmse)[names(incmse)== "P_17_geo_bin_7_150"] <- "P17_geology - central gneiss"
              names(incmse)[names(incmse)== "P_18_geo_bin_8_150"] <- "P18_geology - quartz phyllite"
              names(incmse)[names(incmse)== "P_19_geo_bin_9_150"] <- "P19_geology - insignificant alpine metamorphosis"
              names(incmse)[names(incmse)== "P_20_geo_bi_10_150"] <- "P20_geology - medium-low alpine metamorphosis"
              names(incmse)[names(incmse)== "P_21_geo_bi_11_150"] <- "P21_geology - medium-high alpine metamorphosis"
              names(incmse)[names(incmse)== "P_22_slope_150m"] <- "P22_slope"
              names(incmse)[names(incmse)== "P_23_wetness_150m"] <- "P23_wetness_index"
              names(incmse)[names(incmse)== "P_24_gve_wgs84_150"] <- "P24_GVE"
              names(incmse)[names(incmse)== "P_25_geo_aeubin150"] <- "P25_Acidic soil"
              names(incmse)[names(incmse)== "P_26_geo_beubin150"] <- "P26_Basic soil"
              names(incmse)[names(incmse)== "P_27_geo_aelbin150"] <- "P27_Acidic soil"
              names(incmse)[names(incmse)== "P_28_geo_belbin150"] <- "P28_Basic soil"
              names(incmse)[names(incmse)== "P_29_geo_abelbi150"] <- "P29_Acidic-Basic soil"
            } else if (i==2) {
              names(incmse)[names(incmse)== "P_1_aspect_ea_150"] <- "P1_aspect east"
              names(incmse)[names(incmse)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
              names(incmse)[names(incmse)== "P_3_aspect_no_150"] <- "P3_aspect north"
              names(incmse)[names(incmse)== "P_4_aspect_so_150"] <- "P4_aspect south"
              names(incmse)[names(incmse)== "P_5_aspect_we_150"] <- "P5_aspect west"
              names(incmse)[names(incmse)== "P_6_cl_grarp1_150"] <- "P6_land use - alpine meadow"
              names(incmse)[names(incmse)== "P_7_cl_grarp2_150"] <- "P7_land use - grassland"
              names(incmse)[names(incmse)== "P_8_cl_grarp3_150"] <- "P8_land use - grassland second shearing"
              names(incmse)[names(incmse)== "P_9_cl_grarp4_150"] <- "P9_land use - grassland special area"
              names(incmse)[names(incmse)== "P_10_cl_grarp5_150"] <- "P10_land use - crops"
              names(incmse)[names(incmse)== "P_11_dtm_150m"] <- "P11_dtm"
              names(incmse)[names(incmse)== "P_12_geo_bin_1_150"] <- "P12_geology - shist"
              names(incmse)[names(incmse)== "P_13_geo_bin_2_150"] <- "P13_geology - vulcanite"
              names(incmse)[names(incmse)== "P_14_geo_bin_3_150"] <- "P14_geology - vinschgau shear zone"
              names(incmse)[names(incmse)== "P_15_geo_bin_4_150"] <- "P15_geology - quaternary depositions"
              names(incmse)[names(incmse)== "P_16_geo_bin_5_150"] <- "P16_geology - sediment sequences"
              names(incmse)[names(incmse)== "P_17_geo_bin_6_150"] <- "P17_geology - plutons"
              names(incmse)[names(incmse)== "P_18_geo_bin_7_150"] <- "P18_geology - central gneiss"
              names(incmse)[names(incmse)== "P_19_geo_bin_8_150"] <- "P19_geology - quartz phyllite"
              names(incmse)[names(incmse)== "P_20_geo_bin_9_150"] <- "P20_geology - insignificant alpine metamorphosis"
              names(incmse)[names(incmse)== "P_21_geo_bi_10_150"] <- "P21_geology - medium-low alpine metamorphosis"
              names(incmse)[names(incmse)== "P_22_geo_bi_11_150"] <- "P22_geology - medium-high alpine metamorphosis"
              names(incmse)[names(incmse)== "P_23_slope_150m"] <- "P23_slope"
              names(incmse)[names(incmse)== "P_24_wetness_150m"] <- "P24_wetness_index"
              names(incmse)[names(incmse)== "P_25_gve_wgs84_150"] <- "P25_GVE"
              names(incmse)[names(incmse)== "P_26_geo_aeubin150"] <- "P26_Acidic soil"
              names(incmse)[names(incmse)== "P_27_geo_beubin150"] <- "P27_Basic soil"
              names(incmse)[names(incmse)== "P_28_geo_aelbin150"] <- "P28_Acidic soil"
              names(incmse)[names(incmse)== "P_29_geo_belbin150"] <- "P29_Basic soil"
              names(incmse)[names(incmse)== "P_30_geo_abelbi150"] <- "P30_Acidic-Basic soil"
            } else {}
            
            # Save PNG Table %IncMSE
            incmsedff <- data.frame(incmse)
            names(incmsedff) <- ("%IncMSE")
            table <- tableGrob(incmsedff)
            title <- textGrob(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia),gp=gpar(fontsize=8))
            footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
            padding <- unit(2,"line")
            table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
            table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
            table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
            png(file.path(jkp4,paste0("4_IncMse_Table_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")), units="cm", width=(ncol(table))*15, height=(nrow(table)), res=300) #, height = 200, width = 100)
            grid.newpage() 
            grid.draw(table) 
            dev.off()
            # Save VALUES Table %IncMSE
            write.csv(incmse,file.path(jkp4,paste0("4_IncMse_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".csv")))      
            write.xlsx(incmse,file.path(jkp4,paste0("4_IncMse_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".xlsx")))    
            # Save Plot %IncMSE
            x1<-incmsedff[,1]
            x<-rev(x1)
            y1<-names(incmse)
            y<-rev(y1)
            png(file.path(jkp4,paste0("4_IncMse_Plot_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
            dotchart(x,y,xlab="%IncMse") 
            title(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun))
            dev.off()
            
            # IncNodePurity
            IncNodePurity<-round(sort(importance[,2],decreasing = TRUE),digits = 2) 
            # Define names for predictors
            if (i==1) {
              names(IncNodePurity)[names(IncNodePurity)== "P_1_aspect_ea_150"] <- "P1__aspect east"
              names(IncNodePurity)[names(IncNodePurity)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
              names(IncNodePurity)[names(IncNodePurity)== "P_3_aspect_no_150"] <- "P3_aspect north"
              names(IncNodePurity)[names(IncNodePurity)== "P_4_aspect_so_150"] <- "P4_aspect south"
              names(IncNodePurity)[names(IncNodePurity)== "P_5_aspect_we_150"] <- "P5_aspect west"
              names(IncNodePurity)[names(IncNodePurity)== "P_6_cl_grpre1_150"] <- "P6_land use - alpine meadow"
              names(IncNodePurity)[names(IncNodePurity)== "P_7_cl_grpre2_150"] <- "P7_land use - grassland"
              names(IncNodePurity)[names(IncNodePurity)== "P_8_cl_grpre3_150"] <- "P8_land use - grassland second shearing"
              names(IncNodePurity)[names(IncNodePurity)== "P_9_cl_grpre4_150"] <- "P9_land use - grassland special area"
              names(IncNodePurity)[names(IncNodePurity)== "P_10_dtm_150m"] <- "P10_dtm"
              names(IncNodePurity)[names(IncNodePurity)== "P_11_geo_bin_1_150"] <- "P11_geology - shist"
              names(IncNodePurity)[names(IncNodePurity)== "P_12_geo_bin_2_150"] <- "P12_geology - vulcanite"
              names(IncNodePurity)[names(IncNodePurity)== "P_13_geo_bin_3_150"] <- "P13_geology - vinschgau shear zone"
              names(IncNodePurity)[names(IncNodePurity)== "P_14_geo_bin_4_150"] <- "P14_geology - quaternary depositions"
              names(IncNodePurity)[names(IncNodePurity)== "P_15_geo_bin_5_150"] <- "P15_geology - sediment sequences"
              names(IncNodePurity)[names(IncNodePurity)== "P_16_geo_bin_6_150"] <- "P16_geology - plutons"
              names(IncNodePurity)[names(IncNodePurity)== "P_17_geo_bin_7_150"] <- "P17_geology - central gneiss"
              names(IncNodePurity)[names(IncNodePurity)== "P_18_geo_bin_8_150"] <- "P18_geology - quartz phyllite"
              names(IncNodePurity)[names(IncNodePurity)== "P_19_geo_bin_9_150"] <- "P19_geology - insignificant alpine metamorphosis"
              names(IncNodePurity)[names(IncNodePurity)== "P_20_geo_bi_10_150"] <- "P20_geology - medium-low alpine metamorphosis"
              names(IncNodePurity)[names(IncNodePurity)== "P_21_geo_bi_11_150"] <- "P21_geology - medium-high alpine metamorphosis"
              names(IncNodePurity)[names(IncNodePurity)== "P_22_slope_150m"] <- "P22_slope"
              names(IncNodePurity)[names(IncNodePurity)== "P_23_wetness_150m"] <- "P23_wetness_index"
              names(IncNodePurity)[names(IncNodePurity)== "P_24_gve_wgs84_150"] <- "P24_GVE"
              names(IncNodePurity)[names(IncNodePurity)== "P_25_geo_aeubin150"] <- "P25_Acidic soil"
              names(IncNodePurity)[names(IncNodePurity)== "P_26_geo_beubin150"] <- "P26_Basic soil"
              names(IncNodePurity)[names(IncNodePurity)== "P_27_geo_aelbin150"] <- "P27_Acidic soil"
              names(IncNodePurity)[names(IncNodePurity)== "P_28_geo_belbin150"] <- "P28_Basic soil"
              names(IncNodePurity)[names(IncNodePurity)== "P_29_geo_abelbi150"] <- "P29_Acidic-Basic soil"
            } else if (i==2) {
              names(IncNodePurity)[names(IncNodePurity)== "P_1_aspect_ea_150"] <- "P1_aspect east"
              names(IncNodePurity)[names(IncNodePurity)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
              names(IncNodePurity)[names(IncNodePurity)== "P_3_aspect_no_150"] <- "P3_aspect north"
              names(IncNodePurity)[names(IncNodePurity)== "P_4_aspect_so_150"] <- "P4_aspect south"
              names(IncNodePurity)[names(IncNodePurity)== "P_5_aspect_we_150"] <- "P5_aspect west"
              names(IncNodePurity)[names(IncNodePurity)== "P_6_cl_grarp1_150"] <- "P6_land use - alpine meadow"
              names(IncNodePurity)[names(IncNodePurity)== "P_7_cl_grarp2_150"] <- "P7_land use - grassland"
              names(IncNodePurity)[names(IncNodePurity)== "P_8_cl_grarp3_150"] <- "P8_land use - grassland second shearing"
              names(IncNodePurity)[names(IncNodePurity)== "P_9_cl_grarp4_150"] <- "P9land use - grassland special area"
              names(IncNodePurity)[names(IncNodePurity)== "P_10_cl_grarp5_150"] <- "P10_land use - crops"
              names(IncNodePurity)[names(IncNodePurity)== "P_11_dtm_150m"] <- "P11_dtm"
              names(IncNodePurity)[names(IncNodePurity)== "P_12_geo_bin_1_150"] <- "P12_geology - shist"
              names(IncNodePurity)[names(IncNodePurity)== "P_13_geo_bin_2_150"] <- "P13_geology - vulcanite"
              names(IncNodePurity)[names(IncNodePurity)== "P_14_geo_bin_3_150"] <- "P14_geology - vinschgau shear zone"
              names(IncNodePurity)[names(IncNodePurity)== "P_15_geo_bin_4_150"] <- "P15_geology - quaternary depositions"
              names(IncNodePurity)[names(IncNodePurity)== "P_16_geo_bin_5_150"] <- "P16_geology - sediment sequences"
              names(IncNodePurity)[names(IncNodePurity)== "P_17_geo_bin_6_150"] <- "P17_geology - plutons"
              names(IncNodePurity)[names(IncNodePurity)== "P_18_geo_bin_7_150"] <- "P18_geology - central gneiss"
              names(IncNodePurity)[names(IncNodePurity)== "P_19_geo_bin_8_150"] <- "P19_geology - quartz phyllite"
              names(IncNodePurity)[names(IncNodePurity)== "P_20_geo_bin_9_150"] <- "P20_geology - insignificant alpine metamorphosis"
              names(IncNodePurity)[names(IncNodePurity)== "P_21_geo_bi_10_150"] <- "P21_geology - medium-low alpine metamorphosis"
              names(IncNodePurity)[names(IncNodePurity)== "P_22_geo_bi_11_150"] <- "P22_geology - medium-high alpine metamorphosis"
              names(IncNodePurity)[names(IncNodePurity)== "P_23_slope_150m"] <- "P23_slope"
              names(IncNodePurity)[names(IncNodePurity)== "P_24_wetness_150m"] <- "P24_wetness_index"
              names(IncNodePurity)[names(IncNodePurity)== "P_25_gve_wgs84_150"] <- "P25_GVE"
              names(IncNodePurity)[names(IncNodePurity)== "P_26_geo_aeubin150"] <- "P26_Acidic soil"
              names(IncNodePurity)[names(IncNodePurity)== "P_27_geo_beubin150"] <- "P27_Basic soil"
              names(IncNodePurity)[names(IncNodePurity)== "P_28_geo_aelbin150"] <- "P28_Acidic soil"
              names(IncNodePurity)[names(IncNodePurity)== "P_29_geo_belbin150"] <- "P29_Basic soil"
              names(IncNodePurity)[names(IncNodePurity)== "P_30_geo_abelbi150"] <- "P30_Acidic-Basic soil"
            } else {}
            
            # Save PNG Table IncNodePurity
            incnodedff <- data.frame(IncNodePurity)
            names(incnodedff) <- ("IncNodePurity")
            table <- tableGrob(incnodedff)
            title <- textGrob(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia),gp=gpar(fontsize=8))
            footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
            padding <- unit(2,"line")
            table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
            table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
            table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
            png(file.path(jkp4,paste0("4_IncNodePurity_Table_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")), units="cm", width=(ncol(table))*15, height=(nrow(table)), res=300) #, height = 200, width = 100)
            grid.newpage() 
            grid.draw(table) 
            dev.off()
            # Save VALUES Table IncNodePurity
            write.csv(IncNodePurity,file.path(jkp4,paste0("4_IncNodePurity_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".csv")))      
            write.xlsx(IncNodePurity,file.path(jkp4,paste0("4_IncNodePurity_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".xlsx"))) 
            # Save Plot IncNodePurity
            x1<-incnodedff[,1]
            x<-rev(x1)
            y1<-names(IncNodePurity)
            y<-rev(y1)
            png(file.path(jkp4,paste0("4_IncNodePurity_Plot_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
            dotchart(x,y,xlab="IncNodePurity") 
            title(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun))
            dev.off()
            
            # ImportanceSD (standard deviation)
            importancesdd<-modelRf@regModel$importanceSD
            importancesd<-sort(importancesdd,decreasing=TRUE)
            if (i==1) {
              names(importancesd)[names(importancesd)== "P_1_aspect_ea_150"] <- "P1__aspect east"
              names(importancesd)[names(importancesd)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
              names(importancesd)[names(importancesd)== "P_3_aspect_no_150"] <- "P3_aspect north"
              names(importancesd)[names(importancesd)== "P_4_aspect_so_150"] <- "P4_aspect south"
              names(importancesd)[names(importancesd)== "P_5_aspect_we_150"] <- "P5_aspect west"
              names(importancesd)[names(importancesd)== "P_6_cl_grpre1_150"] <- "P6_land use - alpine meadow"
              names(importancesd)[names(importancesd)== "P_7_cl_grpre2_150"] <- "P7_land use - grassland"
              names(importancesd)[names(importancesd)== "P_8_cl_grpre3_150"] <- "P8_land use - grassland second shearing"
              names(importancesd)[names(importancesd)== "P_9_cl_grpre4_150"] <- "P9_land use - grassland special area"
              names(importancesd)[names(importancesd)== "P_10_dtm_150m"] <- "P10_dtm"
              names(importancesd)[names(importancesd)== "P_11_geo_bin_1_150"] <- "P11_geology - shist"
              names(importancesd)[names(importancesd)== "P_12_geo_bin_2_150"] <- "P12_geology - vulcanite"
              names(importancesd)[names(importancesd)== "P_13_geo_bin_3_150"] <- "P13_geology - vinschgau shear zone"
              names(importancesd)[names(importancesd)== "P_14_geo_bin_4_150"] <- "P14_geology - quaternary depositions"
              names(importancesd)[names(importancesd)== "P_15_geo_bin_5_150"] <- "P15_geology - sediment sequences"
              names(importancesd)[names(importancesd)== "P_16_geo_bin_6_150"] <- "P16_geology - plutons"
              names(importancesd)[names(importancesd)== "P_17_geo_bin_7_150"] <- "P17_geology - central gneiss"
              names(importancesd)[names(importancesd)== "P_18_geo_bin_8_150"] <- "P18_geology - quartz phyllite"
              names(importancesd)[names(importancesd)== "P_19_geo_bin_9_150"] <- "P19_geology - insignificant alpine metamorphosis"
              names(importancesd)[names(importancesd)== "P_20_geo_bi_10_150"] <- "P20_geology - medium-low alpine metamorphosis"
              names(importancesd)[names(importancesd)== "P_21_geo_bi_11_150"] <- "P21_geology - medium-high alpine metamorphosis"
              names(importancesd)[names(importancesd)== "P_22_slope_150m"] <- "P22_slope"
              names(importancesd)[names(importancesd)== "P_23_wetness_150m"] <- "P23_wetness_index"
              names(importancesd)[names(importancesd)== "P_24_gve_wgs84_150"] <- "P24_GVE"
              names(importancesd)[names(importancesd)== "P_25_geo_aeubin150"] <- "P25_Acidic soil"
              names(importancesd)[names(importancesd)== "P_26_geo_beubin150"] <- "P26_Basic soil"
              names(importancesd)[names(importancesd)== "P_27_geo_aelbin150"] <- "P27_Acidic soil"
              names(importancesd)[names(importancesd)== "P_28_geo_belbin150"] <- "P28_Basic soil"
              names(importancesd)[names(importancesd)== "P_29_geo_abelbi150"] <- "P29_Acidic-Basic soil"
            } else if (i==2) {
              names(importancesd)[names(importancesd)== "P_1_aspect_ea_150"] <- "P1_aspect east"
              names(importancesd)[names(importancesd)== "P_2_aspect_fl_150"] <- "P2_aspect flat"
              names(importancesd)[names(importancesd)== "P_3_aspect_no_150"] <- "P3_aspect north"
              names(importancesd)[names(importancesd)== "P_4_aspect_so_150"] <- "P4_aspect south"
              names(importancesd)[names(importancesd)== "P_5_aspect_we_150"] <- "P5_aspect west"
              names(importancesd)[names(importancesd)== "P_6_cl_grarp1_150"] <- "P6_land use - alpine meadow"
              names(importancesd)[names(importancesd)== "P_7_cl_grarp2_150"] <- "P7_land use - grassland"
              names(importancesd)[names(importancesd)== "P_8_cl_grarp3_150"] <- "P8_land use - grassland second shearing"
              names(importancesd)[names(importancesd)== "P_9_cl_grarp4_150"] <- "P9land use - grassland special area"
              names(importancesd)[names(importancesd)== "P_10_cl_grarp5_150"] <- "P10_land use - crops"
              names(importancesd)[names(importancesd)== "P_11_dtm_150m"] <- "P11_dtm"
              names(importancesd)[names(importancesd)== "P_12_geo_bin_1_150"] <- "P12_geology - shist"
              names(importancesd)[names(importancesd)== "P_13_geo_bin_2_150"] <- "P13_geology - vulcanite"
              names(importancesd)[names(importancesd)== "P_14_geo_bin_3_150"] <- "P14_geology - vinschgau shear zone"
              names(importancesd)[names(importancesd)== "P_15_geo_bin_4_150"] <- "P15_geology - quaternary depositions"
              names(importancesd)[names(importancesd)== "P_16_geo_bin_5_150"] <- "P16_geology - sediment sequences"
              names(importancesd)[names(importancesd)== "P_17_geo_bin_6_150"] <- "P17_geology - plutons"
              names(importancesd)[names(importancesd)== "P_18_geo_bin_7_150"] <- "P18_geology - central gneiss"
              names(importancesd)[names(importancesd)== "P_19_geo_bin_8_150"] <- "P19_geology - quartz phyllite"
              names(importancesd)[names(importancesd)== "P_20_geo_bin_9_150"] <- "P20_geology - insignificant alpine metamorphosis"
              names(importancesd)[names(importancesd)== "P_21_geo_bi_10_150"] <- "P21_geology - medium-low alpine metamorphosis"
              names(importancesd)[names(importancesd)== "P_22_geo_bi_11_150"] <- "P22_geology - medium-high alpine metamorphosis"
              names(importancesd)[names(importancesd)== "P_23_slope_150m"] <- "P23_slope"
              names(importancesd)[names(importancesd)== "P_24_wetness_150m"] <- "P24_wetness_index"
              names(importancesd)[names(importancesd)== "P_25_gve_wgs84_150"] <- "P25_GVE"
              names(importancesd)[names(importancesd)== "P_26_geo_aeubin150"] <- "P26_Acidic soil"
              names(importancesd)[names(importancesd)== "P_27_geo_beubin150"] <- "P27_Basic soil"
              names(importancesd)[names(importancesd)== "P_28_geo_aelbin150"] <- "P28_Acidic soil"
              names(importancesd)[names(importancesd)== "P_29_geo_belbin150"] <- "P29_Basic soil"
              names(importancesd)[names(importancesd)== "P_30_geo_abelbi150"] <- "P30_Acidic-Basic soil"
            } else {}
            
            # Save PNG Table importancesd
            importancesddff <- data.frame(importancesd)
            names(importancesddff) <- ("ImportanceSD")
            table <- tableGrob(importancesddff)
            title <- textGrob(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia),gp=gpar(fontsize=8))
            footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
            padding <- unit(2,"line")
            table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
            table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
            table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
            png(file.path(jkp4,paste0("4_ImportanceSD_Table_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")), units="cm", width=(ncol(table))*15, height=(nrow(table)), res=300) #, height = 200, width = 100)
            grid.newpage() 
            grid.draw(table) 
            dev.off()
            # Save VALUES Table importancesd
            write.csv(importancesd,file.path(jkp4,paste0("4_ImportanceSD_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".csv")))      
            write.xlsx(importancesd,file.path(jkp4,paste0("4_ImportanceSD_ValuesTable_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".xlsx"))) 
            # Save Plot importancesd
            x1<-importancesddff[,1]
            x<-rev(x1)
            y1<-names(importancesd)
            y<-rev(y1)
            png(file.path(jkp4,paste0("4_ImportanceSD_Plot_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
            dotchart(x,y,xlab="ImportanceSD") 
            title(paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun))
            dev.off()
            
            # Mean of square residuals / % Var. explained
            printregmodel<-modelRf@regModel
            sink(file.path(jkp2,paste0("2_RfModel_Values",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".txt")))
            print(printregmodel)
            sink()
            setwd(paste0(out_dir,"/",predictareas_name[i],"/",var,"/",ai[k],"/","No.Run_",nrun,"/Random_Forest_Regression/2_RfModel"))
            armin<-read.delim(paste0("2_RfModel_Values",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia,".txt"),stringsAsFactors = FALSE)
            meanofsquareresiduals<-readr::parse_number(armin[5,])
            varexplained<-readr::parse_number(armin[6,])
            
            
            
            # Predictions ORDINARY KRIGING (with 80% of the data)
            predictionsRf<-GSIF::predict.gstatModel(object = modelRf, predictionLocations = grid, predict.method = "RK",verbose=T,nfold = 5)#vgmmodel = model@vgmModel,
            predictRf<-raster(predictionsRf@predicted)
            if (ai[k]=="raw") {
              
            } else if (ai[k]=="log") {
              predictRf_transf_log<-exp(predictRf)
            } else if (ai[k]=="sqrt") {
              predictRf_transf_sqrt<-quad(predictRf)
            } else if (ai[k]=="quad") {
              predictRf_transf_quad<-sqrt(predictRf)
            } else if (ai[k]=="inverse") {
              predictRf_transf_inverse<-inverse(predictRf)
            } else {}
            # Save Plots
            png(file.path(jkpp,paste0("5_Prediction_OrdinaryKriging_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
            plot(predictionsRf) #paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
            dev.off()
            
            
            
            # Plot Raster of predicted values
            png(file.path(jkppp,paste0("6_Predict_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".png")),res=300,width =2000, height=2000)
            if (ai[k]=="raw") {
              raster::plot(predictRf,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
            } else if (ai[k]=="log") {
              raster::plot(predictRf_transf_log,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
            } else if (ai[k]=="sqrt") {
              raster::plot(predictRf_transf_sqrt,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
            } else if (ai[k]=="quad") {
              raster::plot(predictRf_transf_quad,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
            } else if (ai[k]=="inverse") {
              raster::plot(predictRf_transf_inverse,main=paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
            } else {}
            dev.off()
            
            
            
            
            # Overlay validation and prediction
            obsValidRf<-sp::over(obsValid,predictionsRf@predicted)
            if (ai[k]=="raw") {
              
            } else if (ai[k]=="log") {
              obsValidRf$original<-exp(obsValidRf$log)
            } else if (ai[k]=="sqrt") {
              obsValidRf$original<-quad(obsValidRf$sqrt)
            } else if (ai[k]=="quad") {
              obsValidRf$original<-sqrt(obsValidRf$quad)
            } else if (ai[k]=="inverse") {
              obsValidRf$original<-inverse(obsValidRf$inverse)
            } else {}
            
            
            
            # RMSE & NRMSE & R-SQUARED
            if (ai[k]=="raw") {
              rmseRf <- RMSE(pred = obsValidRf[,4], obs = obsValid@data[[var]], na.rm = T)
              nrmseRf <- nrmse(sim = obsValidRf[,4], obs = obsValid@data[[var]], na.rm = T)
              nrmseMaxMinRf <- nrmse(sim = obsValidRf[,4], obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
              r2Rf <- R2(pred = obsValidRf[,4], obs = obsValid@data[[var]], na.rm = T)
            } else if (ai[k]=="log") {
              rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
              r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
            } else if (ai[k]=="sqrt") {
              rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
              r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
            } else if (ai[k]=="quad") {
              rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
              r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
            } else if (ai[k]=="inverse") {
              rmseRf <- RMSE(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              nrmseRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
              nrmseMaxMinRf <- nrmse(sim = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T, norm = "maxmin")
              r2Rf <- R2(pred = obsValidRf$original, obs = obsValid@data[[var]], na.rm = T)
            } else {}
            
            
            
            # Maps - Final rasters (GTiff)
            if (ai[k]=="raw") {
              writeRaster(x = predictRf, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                          overwrite = TRUE, format = "GTiff")
              save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                   list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                   file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
              
            } else if (ai[k]=="log") {
              writeRaster(x = predictRf_transf_log, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                          overwrite = TRUE, format = "GTiff")
              save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                   list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                   file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
              
            } else if (ai[k]=="sqrt") {
              writeRaster(x = predictRf_transf_sqrt, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                          overwrite = TRUE, format = "GTiff")
              save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                   list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                   file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
              
            } else if (ai[k]=="quad") {
              writeRaster(x = predictRf_transf_quad, filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                          overwrite = TRUE, format = "GTiff")
              save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                   list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                   file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
              
            } else if (ai[k]=="inverse") {
              writeRaster(x = predictRf_transf_inverse,filename = file.path(jkppp,paste0("7_Predicted_FinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".tif")),
                          overwrite = TRUE, format = "GTiff")
              save(modelRf, obs, obsValid, predictionsRf, rmseRf, r2Rf, nrmseRf, nrmseMaxMinRf,
                   list = c("modelRf", "obs", "obsValid", "predictionsRf", "rmseRf", "r2Rf", "nrmseRf", "nrmseMaxMinRf"),
                   file = file.path(jkppp, paste0("7_Values_ModelFinalRaster_",var,"_","predictarea",i,"_",ai[k],"_","run",nrun,".RData")))
              
            } else {}
            
            # spalte 1 - No.Run
            sdf1<-length(spalte1)
            spalte1<-append(spalte1,nrun,after=sdf1)
            nrun<-nrun+1
            # Spalte 2 - Transformation
            sdf2<-length(spalte2)
            spalte2<-append(spalte2,ai[k],after=sdf2)
            # Spalte 3 - Cutoff
            sdf3<-length(spalte3)
            spalte3<-append(spalte3,taia,after=sdf3)
            # Spalte 4 - Predictors
            sdf4<-length(spalte4)
            spalte4<-append(spalte4,paste(get(paste0("predictors_combi_",l)), collapse = "/"))
            # Spalte 5 - Mean of squared residuals
            sdf5<-length(spalte5)
            spalte5<-append(spalte5,meanofsquareresiduals,after=sdf5)
            # Spalte 6 - % Var explained
            sdf6<-length(spalte6)
            spalte6<-append(spalte6,varexplained,after=sdf6)
            # Spalte 7 - RMSE
            sdf7<-length(spalte7)
            spalte7<-append(spalte7,rmseRf,after=sdf7)
            sdf7tot<-length(spalte7tot)
            spalte7tot<-append(spalte7tot,rmseRf,after=sdf7tot)
            # Spalte 8 - NRMSE
            sdf8<-length(spalte8)
            spalte8<-append(spalte8,nrmseRf,after=sdf8)
            # Spalte 9 - NRMSEMAXMIN
            sdf9<-length(spalte9)
            spalte9<-append(spalte9,nrmseMaxMinRf,after=sdf9)
            # Spalte 10 - R-SQUARED
            sdf10<-length(spalte10)
            spalte10<-append(spalte10,r2Rf,after=sdf10)
            sdf10tot<-length(spalte10tot)
            spalte10tot<-append(spalte10tot,r2Rf,after=sdf10tot)
            
            
          } else {}
          
          try(dev.off(which = dev.cur()))
          #end_time <- Sys.time()
        }
        
        #time.taken <- end_time - start_time
        #time.taken
      }
    }
    
    # Copy plots describing the correlation coefficients 
    memory.limit(size=9999999999999)
    for (b in 1:nrpredictareas) {
      pre <- list.files(paste0(out_dir,"/",predictareas_name[b]),paste0("Correlation"),recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
      prenew<-paste0(script_dir,"/output_overview")
      file.copy(pre,prenew,recursive=TRUE)
    }
    
    
    if (length(list.files(paste0(out_dir,"/",predictareas_name[i],"/",variables[j]),"No.Run",recursive=TRUE, full.names=FALSE, include.dirs=TRUE))==0) {
      unlink(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/*"), recursive=TRUE)
      setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var))
      fileConn<-file("Read_Me.txt")
      writeLines(c(paste0("The source data of the dependent variable ","'",variables[j],"'"," is not good enough."),"Therefore no prediction and output data could be generated.","One cause could be that there are too few available soil data points of this variable."), fileConn)
      close(fileConn)
      
    } else {
      
      if (i==1) {
        overviewdir<-paste0(script_dir,"/output_overview/Overview_Predict_area_1/",variables[j],"/Overview_Tables")
      } else if (i==2) {
        overviewdir<-paste0(script_dir,"/output_overview/Overview_Predict_area_2/",variables[j],"/Overview_Tables/")
      } else {}
      
      # Total overview table
      memory.limit(size=999999999)
      finaltabledf <- data.frame(spalte1,spalte2,spalte3,spalte4,spalte5,spalte6,spalte7,spalte8,spalte9,spalte10)
      names(finaltabledf)<-c("No.Run","Transformation","Cutoff","Predictors","Mean_Of_Squared_Residuals","%VarExplained","RMSE","NRMSE","NRMSEMAXMIN","R_SQUARED")
      table <- tableGrob(finaltabledf)
      title <- textGrob(paste0(variables[j]," - predictarea",i),gp=gpar(fontsize=15))
      footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
      padding <- unit(2,"line")
      table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
      table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
      table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
      #png(file.path(overviewdir,paste0("8_Output_Overview_Total_predictarea",i,".png")), units="cm", width=(ncol(table))*4, height=(nrow(table)), res=300) #, height = 200, width = 100)
      #grid.newpage() 
      #grid.draw(table) 
      #dev.off()
      write.csv(finaltabledf,file.path(overviewdir,paste0("8_Output_Overview_Total_predictarea",i,".csv")))      
      write.xlsx(finaltabledf,file.path(overviewdir,paste0("8_Output_Overview_Total_predictarea",i,".xlsx")))
      # Total overview table (order by highest R2)
      memory.limit(size=999999999)
      finaltabledfhighestr2 <- data.frame(spalte1,spalte2,spalte3,spalte4,spalte5,spalte6,spalte7,spalte8,spalte9,spalte10)
      finaltabledfhighestr2 <- finaltabledfhighestr2[order(-finaltabledfhighestr2[,10]),]
      names(finaltabledfhighestr2)<-c("No.Run","Transformation","Cutoff","Predictors","Mean_Of_Squared_Residuals","%VarExplained","RMSE","NRMSE","NRMSEMAXMIN","R_SQUARED")
      table <- tableGrob(finaltabledfhighestr2)
      title <- textGrob(paste0(variables[j]," - predictarea",i," - R-Squared (descending order)"),gp=gpar(fontsize=15))
      footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
      padding <- unit(2,"line")
      table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
      table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
      table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
      #png(file.path(overviewdir,paste0("8_Output_Overview_Total_predictarea",i,"_HighestR2.png")), units="cm", width=(ncol(table))*4, height=(nrow(table)), res=300) #, height = 200, width = 100)
      #grid.newpage() 
      #grid.draw(table) 
      #dev.off()
      write.csv(finaltabledfhighestr2,file.path(overviewdir,paste0("8_Output_Overview_Total_predictarea",i,"_HighestR2.csv")))      
      write.xlsx(finaltabledfhighestr2,file.path(overviewdir,paste0("8_Output_Overview_Total_predictarea",i,"_HighestR2.xlsx")))
      # Total overview table (order by lowest RMSE)
      memory.limit(size=999999999)
      finaltabledflowestRMSE <- data.frame(spalte1,spalte2,spalte3,spalte4,spalte5,spalte6,spalte7,spalte8,spalte9,spalte10)
      finaltabledflowestRMSE <- finaltabledflowestRMSE[order(finaltabledflowestRMSE[,7]),]
      names(finaltabledflowestRMSE)<-c("No.Run","Transformation","Cutoff","Predictors","Mean_Of_Squared_Residuals","%VarExplained","RMSE","NRMSE","NRMSEMAXMIN","R_SQUARED")
      table <- tableGrob(finaltabledflowestRMSE)
      title <- textGrob(paste0(variables[j]," - predictarea",i," - RMSE (ascending order)"),gp=gpar(fontsize=15))
      footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
      padding <- unit(2,"line")
      table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
      table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
      table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
      #png(file.path(overviewdir,paste0("8_Output_Overview_Total_predictarea",i,"_LowestRMSE.png")), units="cm", width=(ncol(table))*4, height=(nrow(table)), res=300) #, height = 200, width = 100)
      #grid.newpage() 
      #grid.draw(table) 
      #dev.off()
      write.csv(finaltabledflowestRMSE,file.path(overviewdir,paste0("8_Output_Overview_Total_predictarea",i,"_LowestRMSE.csv")))      
      write.xlsx(finaltabledflowestRMSE,file.path(overviewdir,paste0("8_Output_Overview_Total_predictarea",i,"_LowestRMSE.xlsx")))
      
      if (length(spalte10)>0 & length(spalte10)<=5) {
        
        if (length(spalte10)==1) {
          copyspalte10<-spalte10
          copyspalte10tot<-spalte10tot
          highestr2data<-sort(spalte10,decreasing = TRUE) 
          index_highest_r2_11<-which(copyspalte10==highestr2data[1]) 
          index_highest_r2_1<-which(copyspalte10tot==highestr2data[1])
          copyspalte7<-spalte7
          copyspalte7tot<-spalte7tot
          lowestrmsedata<-sort(spalte7,decreasing = FALSE)
          index_lowest_rmse_11<-which(copyspalte7==lowestrmsedata[1])
          index_lowest_rmse_1<-which(copyspalte7tot==lowestrmsedata[1])
          # Copy Runs Output
          bestindex<-vector()
          bestindexwithoutduplicates<-vector()
          bestindex<-c(index_highest_r2_1,
                       index_lowest_rmse_1)
          bestindexwithoutduplicates<-unique(bestindex)
          setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output"))
          for (b in 1:length(bestindexwithoutduplicates)) {
            dir.create(paste0("No_Run_",bestindexwithoutduplicates[b]))
            pre <- list.files(paste0(out_dir,"/",predictareas_name[i],"/",var), paste0("No.Run_",bestindexwithoutduplicates[b]) , recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
            pree <- paste0(pre,"/Final_Models")
            pree <- pree[1]
            prenew<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output","/No_Run_",bestindexwithoutduplicates[b])
            file.copy(pree,prenew,recursive=TRUE)
            pree <- paste0(pre,"/Ordinary_Kriging_Prediction")
            pree <- pree[1]
            file.copy(pree,prenew,recursive=TRUE)
            pree <- paste0(pre,"/Random_Forest_Regression")
            pree <- pree[1]
            file.copy(pree,prenew,recursive=TRUE)
          }
          
          # Overview Ordinary Kriging Output
          setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output"))
          aii<-dir()[file.info(dir())$isdir]
          for (b in 1:length(aii)) {
            rn<-extract_numeric(aii[b])
            pre2 <- list.files(path=paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output/",aii[b],"/Ordinary_Kriging_Prediction"), full.names = TRUE, recursive = TRUE)
            #pre2 <- paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output/",aii[b],"/Ordinary_Kriging_Prediction/","5_Prediction_OrdinaryKriging_",var,"_","predictarea",i,"_",ai[k],"_","run",rn,".png")
            prenew2<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
            file.copy(pre2,prenew2,recursive=TRUE)
          }
          
          zr<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
          filesvector <- list.files(zr)
          for (b in 1:length(aii)) {
            memory.limit(size=9999999999999)
            rn<-extract_numeric(aii[b])
            zr<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
            img1 <- readImage(file.path(list.files(zr,paste0("run",rn,".png"),full.names = TRUE, recursive = TRUE)))
            index_transformation <- which(spalte1==rn)
            transformationhere <- spalte2[index_transformation]
            #img1 <- readImage(file.path(zr,paste0("5_Prediction_OrdinaryKriging_",var,"_predictarea",i,"_",ai[k],"_run",rn,".png")))
            img2 <- img1[1018:2000,1:1000,]
            png(file.path(zr,paste0("No.Run_",rn,".png")),res=300,width = 2000,height = 2000)
            plot(img2) #paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
            title(paste0("No.Run",rn," - ",transformationhere))
            dev.off()
            #assign(paste0("plot",rn),readImage(file.path(list.files(zr,paste0("No.Run_",rn,".png"),full.names = TRUE, recursive = TRUE))))
            #assign(paste0("plot",b),readPNG(file.path(zr,paste0("No.Run_",rn,".png"))))
          }
          
          memory.limit(size=9999999999999)
          index_lowest_rmse_total <- c(index_lowest_rmse_1)
          index_highest_r2_total <- c(index_highest_r2_1)
          fz<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
          fx<-paste0("5_Overview_Ordinary_Kriging_",var,"_predictarea",i,"_LowestRMSE.png")
          fxx<-paste0("5_Overview_Ordinary_Kriging_",var,"_predictarea",i,"_HighestR2.png")
          plotted1 <- readImage(file.path(paste0(fz,"/No.Run_",index_lowest_rmse_total[1],".png")))
          png(file.path(fz,fx),res=300,width =10000, height=2000)
          grid.arrange(rasterGrob(plotted1), nrow=1, ncol=5)
          dev.off()
          plottted1 <- readImage(file.path(paste0(fz,"/No.Run_",index_highest_r2_total[1],".png")))
          png(file.path(fz,fxx),res=300,width =10000, height=2000)
          grid.arrange(rasterGrob(plottted1), nrow=1, ncol=5)
          dev.off()
          
        } else if (length(spalte10)==2) {
          copyspalte10<-spalte10
          copyspalte10tot<-spalte10tot
          highestr2data<-sort(spalte10,decreasing = TRUE)  
          index_highest_r2_11<-which(copyspalte10==highestr2data[1]) 
          index_highest_r2_22<-which(copyspalte10==highestr2data[2]) 
          index_highest_r2_1<-which(copyspalte10tot==highestr2data[1])  
          index_highest_r2_2<-which(copyspalte10tot==highestr2data[2])
          copyspalte7<-spalte7
          copyspalte7tot<-spalte7tot
          lowestrmsedata<-sort(spalte7,decreasing = FALSE)
          index_lowest_rmse_11<-which(copyspalte7==lowestrmsedata[1])
          index_lowest_rmse_22<-which(copyspalte7==lowestrmsedata[2])
          index_lowest_rmse_1<-which(copyspalte7tot==lowestrmsedata[1])
          index_lowest_rmse_2<-which(copyspalte7tot==lowestrmsedata[2])
          # Copy Runs Output
          bestindex<-vector()
          bestindexwithoutduplicates<-vector()
          bestindex<-c(index_highest_r2_1,index_highest_r2_2,
                       index_lowest_rmse_1,index_lowest_rmse_2)
          bestindexwithoutduplicates<-unique(bestindex)
          setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output"))
          for (b in 1:length(bestindexwithoutduplicates)) {
            dir.create(paste0("No_Run_",bestindexwithoutduplicates[b]))
            pre <- list.files(paste0(out_dir,"/",predictareas_name[i],"/",var), paste0("No.Run_",bestindexwithoutduplicates[b]) , recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
            pree <- paste0(pre,"/Final_Models")
            pree <- pree[1]
            prenew<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output","/No_Run_",bestindexwithoutduplicates[b])
            file.copy(pree,prenew,recursive=TRUE)
            pree <- paste0(pre,"/Ordinary_Kriging_Prediction")
            pree <- pree[1]
            file.copy(pree,prenew,recursive=TRUE)
            pree <- paste0(pre,"/Random_Forest_Regression")
            pree <- pree[1]
            file.copy(pree,prenew,recursive=TRUE)
          }
          
          # Overview Ordinary Kriging Output
          setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output"))
          aii<-dir()[file.info(dir())$isdir]
          for (b in 1:length(aii)) {
            rn<-extract_numeric(aii[b])
            pre2 <- list.files(path=paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output/",aii[b],"/Ordinary_Kriging_Prediction"), full.names = TRUE, recursive = TRUE)
            #pre2 <- paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output/",aii[b],"/Ordinary_Kriging_Prediction/","5_Prediction_OrdinaryKriging_",var,"_","predictarea",i,"_",ai[k],"_","run",rn,".png")
            prenew2<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
            file.copy(pre2,prenew2,recursive=TRUE)
          }
          
          zr<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
          filesvector <- list.files(zr)
          for (b in 1:length(aii)) {
            memory.limit(size=9999999999999)
            rn<-extract_numeric(aii[b])
            zr<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
            img1 <- readImage(file.path(list.files(zr,paste0("run",rn,".png"),full.names = TRUE, recursive = TRUE)))
            index_transformation <- which(spalte1==rn)
            transformationhere <- spalte2[index_transformation]
            #img1 <- readImage(file.path(zr,paste0("5_Prediction_OrdinaryKriging_",var,"_predictarea",i,"_",ai[k],"_run",rn,".png")))
            img2 <- img1[1018:2000,1:1000,]
            png(file.path(zr,paste0("No.Run_",rn,".png")),res=300,width = 2000,height = 2000)
            plot(img2) #paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
            title(paste0("No.Run",rn," - ",transformationhere))
            dev.off()
            #assign(paste0("plot",rn),readImage(file.path(list.files(zr,paste0("No.Run_",rn,".png"),full.names = TRUE, recursive = TRUE))))
            #assign(paste0("plot",b),readPNG(file.path(zr,paste0("No.Run_",rn,".png"))))
          }
          
          memory.limit(size=9999999999999)
          index_lowest_rmse_total <- c(index_lowest_rmse_1,index_lowest_rmse_2)
          index_highest_r2_total <- c(index_highest_r2_1,index_highest_r2_2)
          fz<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
          fx<-paste0("5_Overview_Ordinary_Kriging_",var,"_predictarea",i,"_LowestRMSE.png")
          fxx<-paste0("5_Overview_Ordinary_Kriging_",var,"_predictarea",i,"_HighestR2.png")
          plotted1 <- readImage(file.path(paste0(fz,"/No.Run_",index_lowest_rmse_total[1],".png")))
          plotted2 <- readImage(file.path(paste0(fz,"/No.Run_",index_lowest_rmse_total[2],".png")))
          png(file.path(fz,fx),res=300,width =10000, height=2000)
          grid.arrange(rasterGrob(plotted1),rasterGrob(plotted2), nrow=1, ncol=5)
          dev.off()
          plottted1 <- readImage(file.path(paste0(fz,"/No.Run_",index_highest_r2_total[1],".png")))
          plottted2 <- readImage(file.path(paste0(fz,"/No.Run_",index_highest_r2_total[2],".png")))
          png(file.path(fz,fxx),res=300,width =10000, height=2000)
          grid.arrange(rasterGrob(plottted1),rasterGrob(plottted2), nrow=1, ncol=5)
          dev.off()
          
        } else if (length(spalte10)==3) {
          copyspalte10<-spalte10
          copyspalte10tot<-spalte10tot
          highestr2data<-sort(spalte10,decreasing = TRUE)  
          index_highest_r2_11<-which(copyspalte10==highestr2data[1]) 
          index_highest_r2_22<-which(copyspalte10==highestr2data[2]) 
          index_highest_r2_33<-which(copyspalte10==highestr2data[3]) 
          index_highest_r2_1<-which(copyspalte10tot==highestr2data[1])  
          index_highest_r2_2<-which(copyspalte10tot==highestr2data[2])
          index_highest_r2_3<-which(copyspalte10tot==highestr2data[3])
          copyspalte7<-spalte7
          copyspalte7tot<-spalte7tot
          lowestrmsedata<-sort(spalte7,decreasing = FALSE)
          index_lowest_rmse_11<-which(copyspalte7==lowestrmsedata[1])
          index_lowest_rmse_22<-which(copyspalte7==lowestrmsedata[2])
          index_lowest_rmse_33<-which(copyspalte7==lowestrmsedata[3])
          index_lowest_rmse_1<-which(copyspalte7tot==lowestrmsedata[1])
          index_lowest_rmse_2<-which(copyspalte7tot==lowestrmsedata[2])
          index_lowest_rmse_3<-which(copyspalte7tot==lowestrmsedata[3])
          # Copy Runs Output
          bestindex<-vector()
          bestindexwithoutduplicates<-vector()
          bestindex<-c(index_highest_r2_1,index_highest_r2_2,index_highest_r2_3,
                       index_lowest_rmse_1,index_lowest_rmse_2,index_lowest_rmse_3)
          bestindexwithoutduplicates<-unique(bestindex)
          setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output"))
          for (b in 1:length(bestindexwithoutduplicates)) {
            dir.create(paste0("No_Run_",bestindexwithoutduplicates[b]))
            pre <- list.files(paste0(out_dir,"/",predictareas_name[i],"/",var), paste0("No.Run_",bestindexwithoutduplicates[b]) , recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
            pree <- paste0(pre,"/Final_Models")
            pree <- pree[1]
            prenew<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output","/No_Run_",bestindexwithoutduplicates[b])
            file.copy(pree,prenew,recursive=TRUE)
            pree <- paste0(pre,"/Ordinary_Kriging_Prediction")
            pree <- pree[1]
            file.copy(pree,prenew,recursive=TRUE)
            pree <- paste0(pre,"/Random_Forest_Regression")
            pree <- pree[1]
            file.copy(pree,prenew,recursive=TRUE)
          }
          
          # Overview Ordinary Kriging Output
          setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output"))
          aii<-dir()[file.info(dir())$isdir]
          for (b in 1:length(aii)) {
            rn<-extract_numeric(aii[b])
            pre2 <- list.files(path=paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output/",aii[b],"/Ordinary_Kriging_Prediction"), full.names = TRUE, recursive = TRUE)
            #pre2 <- paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output/",aii[b],"/Ordinary_Kriging_Prediction/","5_Prediction_OrdinaryKriging_",var,"_","predictarea",i,"_",ai[k],"_","run",rn,".png")
            prenew2<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
            file.copy(pre2,prenew2,recursive=TRUE)
          }
          
          zr<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
          filesvector <- list.files(zr)
          for (b in 1:length(aii)) {
            memory.limit(size=9999999999999)
            rn<-extract_numeric(aii[b])
            zr<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
            img1 <- readImage(file.path(list.files(zr,paste0("run",rn,".png"),full.names = TRUE, recursive = TRUE)))
            index_transformation <- which(spalte1==rn)
            transformationhere <- spalte2[index_transformation]
            #img1 <- readImage(file.path(zr,paste0("5_Prediction_OrdinaryKriging_",var,"_predictarea",i,"_",ai[k],"_run",rn,".png")))
            img2 <- img1[1018:2000,1:1000,]
            png(file.path(zr,paste0("No.Run_",rn,".png")),res=300,width = 2000,height = 2000)
            plot(img2) #paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
            title(paste0("No.Run",rn," - ",transformationhere))
            dev.off()
            #assign(paste0("plot",rn),readImage(file.path(list.files(zr,paste0("No.Run_",rn,".png"),full.names = TRUE, recursive = TRUE))))
            #assign(paste0("plot",b),readPNG(file.path(zr,paste0("No.Run_",rn,".png"))))
          }
          
          memory.limit(size=9999999999999)
          index_lowest_rmse_total <- c(index_lowest_rmse_1,index_lowest_rmse_2,index_lowest_rmse_3)
          index_highest_r2_total <- c(index_highest_r2_1,index_highest_r2_2,index_highest_r2_3)
          fz<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
          fx<-paste0("5_Overview_Ordinary_Kriging_",var,"_predictarea",i,"_LowestRMSE.png")
          fxx<-paste0("5_Overview_Ordinary_Kriging_",var,"_predictarea",i,"_HighestR2.png")
          plotted1 <- readImage(file.path(paste0(fz,"/No.Run_",index_lowest_rmse_total[1],".png")))
          plotted2 <- readImage(file.path(paste0(fz,"/No.Run_",index_lowest_rmse_total[2],".png")))
          plotted3 <- readImage(file.path(paste0(fz,"/No.Run_",index_lowest_rmse_total[3],".png")))
          png(file.path(fz,fx),res=300,width =10000, height=2000)
          grid.arrange(rasterGrob(plotted1),rasterGrob(plotted2),rasterGrob(plotted3), nrow=1, ncol=5)
          dev.off()
          plottted1 <- readImage(file.path(paste0(fz,"/No.Run_",index_highest_r2_total[1],".png")))
          plottted2 <- readImage(file.path(paste0(fz,"/No.Run_",index_highest_r2_total[2],".png")))
          plottted3 <- readImage(file.path(paste0(fz,"/No.Run_",index_highest_r2_total[3],".png")))
          png(file.path(fz,fxx),res=300,width =10000, height=2000)
          grid.arrange(rasterGrob(plottted1),rasterGrob(plottted2),rasterGrob(plottted3), nrow=1, ncol=5)
          dev.off()
          
        } else if (length(spalte10)==4) {
          copyspalte10<-spalte10
          copyspalte10tot<-spalte10tot
          highestr2data<-sort(spalte10,decreasing = TRUE)  
          index_highest_r2_11<-which(copyspalte10==highestr2data[1]) 
          index_highest_r2_22<-which(copyspalte10==highestr2data[2]) 
          index_highest_r2_33<-which(copyspalte10==highestr2data[3]) 
          index_highest_r2_44<-which(copyspalte10==highestr2data[4]) 
          index_highest_r2_1<-which(copyspalte10tot==highestr2data[1])  
          index_highest_r2_2<-which(copyspalte10tot==highestr2data[2])
          index_highest_r2_3<-which(copyspalte10tot==highestr2data[3])
          index_highest_r2_4<-which(copyspalte10tot==highestr2data[4])
          copyspalte7<-spalte7
          copyspalte7tot<-spalte7tot
          lowestrmsedata<-sort(spalte7,decreasing = FALSE)
          index_lowest_rmse_11<-which(copyspalte7==lowestrmsedata[1])
          index_lowest_rmse_22<-which(copyspalte7==lowestrmsedata[2])
          index_lowest_rmse_33<-which(copyspalte7==lowestrmsedata[3])
          index_lowest_rmse_44<-which(copyspalte7==lowestrmsedata[4])
          index_lowest_rmse_1<-which(copyspalte7tot==lowestrmsedata[1])
          index_lowest_rmse_2<-which(copyspalte7tot==lowestrmsedata[2])
          index_lowest_rmse_3<-which(copyspalte7tot==lowestrmsedata[3])
          index_lowest_rmse_4<-which(copyspalte7tot==lowestrmsedata[4])
          # Copy Runs Output
          bestindex<-vector()
          bestindexwithoutduplicates<-vector()
          bestindex<-c(index_highest_r2_1,index_highest_r2_2,index_highest_r2_3,index_highest_r2_4,
                       index_lowest_rmse_1,index_lowest_rmse_2,index_lowest_rmse_3,index_lowest_rmse_4)
          bestindexwithoutduplicates<-unique(bestindex)
          setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output"))
          for (b in 1:length(bestindexwithoutduplicates)) {
            dir.create(paste0("No_Run_",bestindexwithoutduplicates[b]))
            pre <- list.files(paste0(out_dir,"/",predictareas_name[i],"/",var), paste0("No.Run_",bestindexwithoutduplicates[b]) , recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
            pree <- paste0(pre,"/Final_Models")
            pree <- pree[1]
            prenew<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output","/No_Run_",bestindexwithoutduplicates[b])
            file.copy(pree,prenew,recursive=TRUE)
            pree <- paste0(pre,"/Ordinary_Kriging_Prediction")
            pree <- pree[1]
            file.copy(pree,prenew,recursive=TRUE)
            pree <- paste0(pre,"/Random_Forest_Regression")
            pree <- pree[1]
            file.copy(pree,prenew,recursive=TRUE)
          }
          
          # Overview Ordinary Kriging Output
          setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output"))
          aii<-dir()[file.info(dir())$isdir]
          for (b in 1:length(aii)) {
            rn<-extract_numeric(aii[b])
            pre2 <- list.files(path=paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output/",aii[b],"/Ordinary_Kriging_Prediction"), full.names = TRUE, recursive = TRUE)
            #pre2 <- paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output/",aii[b],"/Ordinary_Kriging_Prediction/","5_Prediction_OrdinaryKriging_",var,"_","predictarea",i,"_",ai[k],"_","run",rn,".png")
            prenew2<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
            file.copy(pre2,prenew2,recursive=TRUE)
          }
          
          zr<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
          filesvector <- list.files(zr)
          for (b in 1:length(aii)) {
            memory.limit(size=9999999999999)
            rn<-extract_numeric(aii[b])
            zr<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
            img1 <- readImage(file.path(list.files(zr,paste0("run",rn,".png"),full.names = TRUE, recursive = TRUE)))
            index_transformation <- which(spalte1==rn)
            transformationhere <- spalte2[index_transformation]
            #img1 <- readImage(file.path(zr,paste0("5_Prediction_OrdinaryKriging_",var,"_predictarea",i,"_",ai[k],"_run",rn,".png")))
            img2 <- img1[1018:2000,1:1000,]
            png(file.path(zr,paste0("No.Run_",rn,".png")),res=300,width = 2000,height = 2000)
            plot(img2) #paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
            title(paste0("No.Run",rn," - ",transformationhere))
            dev.off()
            #assign(paste0("plot",rn),readImage(file.path(list.files(zr,paste0("No.Run_",rn,".png"),full.names = TRUE, recursive = TRUE))))
            #assign(paste0("plot",b),readPNG(file.path(zr,paste0("No.Run_",rn,".png"))))
          }
          
          memory.limit(size=9999999999999)
          index_lowest_rmse_total <- c(index_lowest_rmse_1,index_lowest_rmse_2,index_lowest_rmse_3,index_lowest_rmse_4)
          index_highest_r2_total <- c(index_highest_r2_1,index_highest_r2_2,index_highest_r2_3,index_highest_r2_4)
          fz<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
          fx<-paste0("5_Overview_Ordinary_Kriging_",var,"_predictarea",i,"_LowestRMSE.png")
          fxx<-paste0("5_Overview_Ordinary_Kriging_",var,"_predictarea",i,"_HighestR2.png")
          plotted1 <- readImage(file.path(paste0(fz,"/No.Run_",index_lowest_rmse_total[1],".png")))
          plotted2 <- readImage(file.path(paste0(fz,"/No.Run_",index_lowest_rmse_total[2],".png")))
          plotted3 <- readImage(file.path(paste0(fz,"/No.Run_",index_lowest_rmse_total[3],".png")))
          plotted4 <- readImage(file.path(paste0(fz,"/No.Run_",index_lowest_rmse_total[4],".png")))
          png(file.path(fz,fx),res=300,width =10000, height=2000)
          grid.arrange(rasterGrob(plotted1),rasterGrob(plotted2),rasterGrob(plotted3),rasterGrob(plotted4), nrow=1, ncol=5)
          dev.off()
          plottted1 <- readImage(file.path(paste0(fz,"/No.Run_",index_highest_r2_total[1],".png")))
          plottted2 <- readImage(file.path(paste0(fz,"/No.Run_",index_highest_r2_total[2],".png")))
          plottted3 <- readImage(file.path(paste0(fz,"/No.Run_",index_highest_r2_total[3],".png")))
          plottted4 <- readImage(file.path(paste0(fz,"/No.Run_",index_highest_r2_total[4],".png")))
          png(file.path(fz,fxx),res=300,width =10000, height=2000)
          grid.arrange(rasterGrob(plottted1),rasterGrob(plottted2),rasterGrob(plottted3),rasterGrob(plottted4), nrow=1, ncol=5)
          dev.off()
          
        } else if (length(spalte10)==5) {
          copyspalte10<-spalte10
          copyspalte10tot<-spalte10tot
          highestr2data<-sort(spalte10,decreasing = TRUE)  
          index_highest_r2_11<-which(copyspalte10==highestr2data[1]) 
          index_highest_r2_22<-which(copyspalte10==highestr2data[2]) 
          index_highest_r2_33<-which(copyspalte10==highestr2data[3]) 
          index_highest_r2_44<-which(copyspalte10==highestr2data[4])
          index_highest_r2_55<-which(copyspalte10==highestr2data[5])
          index_highest_r2_1<-which(copyspalte10tot==highestr2data[1])  
          index_highest_r2_2<-which(copyspalte10tot==highestr2data[2])
          index_highest_r2_3<-which(copyspalte10tot==highestr2data[3])
          index_highest_r2_4<-which(copyspalte10tot==highestr2data[4])
          index_highest_r2_5<-which(copyspalte10tot==highestr2data[5])
          copyspalte7<-spalte7
          copyspalte7tot<-spalte7tot
          lowestrmsedata<-sort(spalte7,decreasing = FALSE)
          index_lowest_rmse_11<-which(copyspalte7==lowestrmsedata[1])
          index_lowest_rmse_22<-which(copyspalte7==lowestrmsedata[2])
          index_lowest_rmse_33<-which(copyspalte7==lowestrmsedata[3])
          index_lowest_rmse_44<-which(copyspalte7==lowestrmsedata[4])
          index_lowest_rmse_55<-which(copyspalte7==lowestrmsedata[5])
          index_lowest_rmse_1<-which(copyspalte7tot==lowestrmsedata[1])
          index_lowest_rmse_2<-which(copyspalte7tot==lowestrmsedata[2])
          index_lowest_rmse_3<-which(copyspalte7tot==lowestrmsedata[3])
          index_lowest_rmse_4<-which(copyspalte7tot==lowestrmsedata[4])
          index_lowest_rmse_5<-which(copyspalte7tot==lowestrmsedata[5])
          # Copy Runs Output
          bestindex<-vector()
          bestindexwithoutduplicates<-vector()
          bestindex<-c(index_highest_r2_1,index_highest_r2_2,index_highest_r2_3,index_highest_r2_4,index_highest_r2_5,
                       index_lowest_rmse_1,index_lowest_rmse_2,index_lowest_rmse_3,index_lowest_rmse_4,index_lowest_rmse_5)
          bestindexwithoutduplicates<-unique(bestindex)
          setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output"))
          for (b in 1:length(bestindexwithoutduplicates)) {
            dir.create(paste0("No_Run_",bestindexwithoutduplicates[b]))
            pre <- list.files(paste0(out_dir,"/",predictareas_name[i],"/",var), paste0("No.Run_",bestindexwithoutduplicates[b]) , recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
            pree <- paste0(pre,"/Final_Models")
            pree <- pree[1]
            prenew<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output","/No_Run_",bestindexwithoutduplicates[b])
            file.copy(pree,prenew,recursive=TRUE)
            pree <- paste0(pre,"/Ordinary_Kriging_Prediction")
            pree <- pree[1]
            file.copy(pree,prenew,recursive=TRUE)
            pree <- paste0(pre,"/Random_Forest_Regression")
            pree <- pree[1]
            file.copy(pree,prenew,recursive=TRUE)
          }
          
          # Overview Ordinary Kriging Output
          setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output"))
          aii<-dir()[file.info(dir())$isdir]
          for (b in 1:length(aii)) {
            rn<-extract_numeric(aii[b])
            pre2 <- list.files(path=paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output/",aii[b],"/Ordinary_Kriging_Prediction"), full.names = TRUE, recursive = TRUE)
            #pre2 <- paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output/",aii[b],"/Ordinary_Kriging_Prediction/","5_Prediction_OrdinaryKriging_",var,"_","predictarea",i,"_",ai[k],"_","run",rn,".png")
            prenew2<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
            file.copy(pre2,prenew2,recursive=TRUE)
          }
          
          zr<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
          filesvector <- list.files(zr)
          for (b in 1:length(aii)) {
            memory.limit(size=9999999999999)
            rn<-extract_numeric(aii[b])
            zr<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
            img1 <- readImage(file.path(list.files(zr,paste0("run",rn,".png"),full.names = TRUE, recursive = TRUE)))
            index_transformation <- which(spalte1==rn)
            transformationhere <- spalte2[index_transformation]
            #img1 <- readImage(file.path(zr,paste0("5_Prediction_OrdinaryKriging_",var,"_predictarea",i,"_",ai[k],"_run",rn,".png")))
            img2 <- img1[1018:2000,1:1000,]
            png(file.path(zr,paste0("No.Run_",rn,".png")),res=300,width = 2000,height = 2000)
            plot(img2) #paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
            title(paste0("No.Run",rn," - ",transformationhere))
            dev.off()
            #assign(paste0("plot",rn),readImage(file.path(list.files(zr,paste0("No.Run_",rn,".png"),full.names = TRUE, recursive = TRUE))))
            #assign(paste0("plot",b),readPNG(file.path(zr,paste0("No.Run_",rn,".png"))))
          }
          
          memory.limit(size=9999999999999)
          index_lowest_rmse_total <- c(index_lowest_rmse_1,index_lowest_rmse_2,index_lowest_rmse_3,index_lowest_rmse_4,index_lowest_rmse_5)
          index_highest_r2_total <- c(index_highest_r2_1,index_highest_r2_2,index_highest_r2_3,index_highest_r2_4,index_highest_r2_5)
          fz<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
          fx<-paste0("5_Overview_Ordinary_Kriging_",var,"_predictarea",i,"_LowestRMSE.png")
          fxx<-paste0("5_Overview_Ordinary_Kriging_",var,"_predictarea",i,"_HighestR2.png")
          plotted1 <- readImage(file.path(paste0(fz,"/No.Run_",index_lowest_rmse_total[1],".png")))
          plotted2 <- readImage(file.path(paste0(fz,"/No.Run_",index_lowest_rmse_total[2],".png")))
          plotted3 <- readImage(file.path(paste0(fz,"/No.Run_",index_lowest_rmse_total[3],".png")))
          plotted4 <- readImage(file.path(paste0(fz,"/No.Run_",index_lowest_rmse_total[4],".png")))
          plotted5 <- readImage(file.path(paste0(fz,"/No.Run_",index_lowest_rmse_total[5],".png")))
          png(file.path(fz,fx),res=300,width =10000, height=2000)
          grid.arrange(rasterGrob(plotted1),rasterGrob(plotted2),rasterGrob(plotted3),rasterGrob(plotted4),rasterGrob(plotted5), nrow=1, ncol=5)
          dev.off()
          plottted1 <- readImage(file.path(paste0(fz,"/No.Run_",index_highest_r2_total[1],".png")))
          plottted2 <- readImage(file.path(paste0(fz,"/No.Run_",index_highest_r2_total[2],".png")))
          plottted3 <- readImage(file.path(paste0(fz,"/No.Run_",index_highest_r2_total[3],".png")))
          plottted4 <- readImage(file.path(paste0(fz,"/No.Run_",index_highest_r2_total[4],".png")))
          plottted5 <- readImage(file.path(paste0(fz,"/No.Run_",index_highest_r2_total[5],".png")))
          png(file.path(fz,fxx),res=300,width =10000, height=2000)
          grid.arrange(rasterGrob(plottted1),rasterGrob(plottted2),rasterGrob(plottted3),rasterGrob(plottted4),rasterGrob(plotted5), nrow=1, ncol=5)
          dev.off()
          
        } else {}
        
        # 5 runs with lowest RMSE
        memory.limit(size=999999999)
        finaltabledflowestRMSE <- data.frame(spalte1,spalte2,spalte3,spalte4,spalte5,spalte6,spalte7,spalte8,spalte9,spalte10)
        finaltabledflowestRMSE <- finaltabledflowestRMSE[order(finaltabledflowestRMSE[,7]),]
        names(finaltabledflowestRMSE)<-c("No.Run","Transformation","Cutoff","Predictors","Mean_Of_Squared_Residuals","%VarExplained","RMSE","NRMSE","NRMSEMAXMIN","R_SQUARED")
        table <- tableGrob(finaltabledflowestRMSE)
        title <- textGrob(paste0(variables[j]," - predictarea",i," - RMSE (ascending order)"),gp=gpar(fontsize=15))
        footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
        padding <- unit(2,"line")
        table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
        table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
        table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
        png(file.path(overviewdir,paste0("8_Output_Overview_5LowestRMSE_predictarea",i,".png")),width=5200, height=710, res=300)
        grid.newpage() 
        grid.draw(table) 
        dev.off()
        write.csv(finaltabledflowestRMSE,file.path(overviewdir,paste0("8_Output_Overview_5LowestRMSE_predictarea",i,".csv")))      
        write.xlsx(finaltabledflowestRMSE,file.path(overviewdir,paste0("8_Output_Overview_5LowestRMSE_predictarea",i,".xlsx")))
        # 5 runs with highest R2
        memory.limit(size=999999999)
        finaltabledfhighestr2 <- data.frame(spalte1,spalte2,spalte3,spalte4,spalte5,spalte6,spalte7,spalte8,spalte9,spalte10)
        finaltabledfhighestr2 <- finaltabledfhighestr2[order(-finaltabledfhighestr2[,10]),]
        names(finaltabledfhighestr2)<-c("No.Run","Transformation","Cutoff","Predictors","Mean_Of_Squared_Residuals","%VarExplained","RMSE","NRMSE","NRMSEMAXMIN","R_SQUARED")
        table <- tableGrob(finaltabledfhighestr2)
        title <- textGrob(paste0(variables[j]," - predictarea",i," - R-Squared (descending order)"),gp=gpar(fontsize=15))
        footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
        padding <- unit(2,"line")
        table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
        table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
        table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
        png(file.path(overviewdir,paste0("8_Output_Overview_5HighestR2_predictarea",i,".png")),width=5200, height=710, res=300)
        grid.newpage() 
        grid.draw(table) 
        dev.off()
        write.csv(finaltabledfhighestr2,file.path(overviewdir,paste0("8_Output_Overview_5HighestR2_predictarea",i,".csv")))      
        write.xlsx(finaltabledfhighestr2,file.path(overviewdir,paste0("8_Output_Overview_5HighestR2_predictarea",i,".xlsx")))
        
        
        
      } else if (length(spalte10)>5) {
        
        # Find 5 runs with highest R^2
        copyspalte10<-spalte10
        copyspalte10tot<-spalte10tot
        highestr2data<-sort(spalte10,decreasing = TRUE)  
        index_highest_r2_11<-which(copyspalte10==highestr2data[1]) 
        index_highest_r2_22<-which(copyspalte10==highestr2data[2]) 
        index_highest_r2_33<-which(copyspalte10==highestr2data[3]) 
        index_highest_r2_44<-which(copyspalte10==highestr2data[4]) 
        index_highest_r2_55<-which(copyspalte10==highestr2data[5]) 
        index_highest_r2_1<-which(copyspalte10tot==highestr2data[1])  
        index_highest_r2_2<-which(copyspalte10tot==highestr2data[2])
        index_highest_r2_3<-which(copyspalte10tot==highestr2data[3])
        index_highest_r2_4<-which(copyspalte10tot==highestr2data[4])
        index_highest_r2_5<-which(copyspalte10tot==highestr2data[5])
        finaltabledf2 <- finaltabledf
        hilfsvektor<-c(1:length(spalte10))
        hilfsvektor <- hilfsvektor[!hilfsvektor %in% c(index_highest_r2_11,index_highest_r2_22,index_highest_r2_33,index_highest_r2_44,index_highest_r2_55)]
        finaltabledf2 <- finaltabledf2[-hilfsvektor,]
        finaltabledf2 <- finaltabledf2[order(-finaltabledf2[,10]),]
        names(finaltabledf2)<-c("No.Run","Transformation","Cutoff","Predictors","Mean_Of_Squared_Residuals","%VarExplained","RMSE","NRMSE","NRMSEMAXMIN","R_SQUARED")
        table <- tableGrob(finaltabledf2)
        title <- textGrob(paste0(variables[j]," - predictarea",i," - 5 runs with highest R-Squared"),gp=gpar(fontsize=15))
        footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
        padding <- unit(2,"line")
        table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
        table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
        table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
        #png(file.path(overviewdir,paste0("8_Output_Overview_HighestR2_predictarea",i,".png")), units="cm", width=(ncol(table))*4, height=(nrow(table)), res=300) #, height = 200, width = 100)
        png(file.path(overviewdir,paste0("8_Output_Overview_5HighestR2_predictarea",i,".png")), width=5200, height=710, res=300)
        grid.newpage() 
        grid.draw(table) 
        dev.off()
        write.csv(finaltabledf2,file.path(overviewdir,paste0("8_Output_Overview_5HighestR2_predictarea",i,".csv")))      
        write.xlsx(finaltabledf2,file.path(overviewdir,paste0("8_Output_Overview_5HighestR2_predictarea",i,".xlsx")))
        # Find 5 runs with lowest RMSE
        copyspalte7<-spalte7
        copyspalte7tot<-spalte7tot
        lowestrmsedata<-sort(spalte7,decreasing = FALSE)
        index_lowest_rmse_11<-which(copyspalte7==lowestrmsedata[1])
        index_lowest_rmse_22<-which(copyspalte7==lowestrmsedata[2])
        index_lowest_rmse_33<-which(copyspalte7==lowestrmsedata[3])
        index_lowest_rmse_44<-which(copyspalte7==lowestrmsedata[4])
        index_lowest_rmse_55<-which(copyspalte7==lowestrmsedata[5])
        index_lowest_rmse_1<-which(copyspalte7tot==lowestrmsedata[1])
        index_lowest_rmse_2<-which(copyspalte7tot==lowestrmsedata[2])
        index_lowest_rmse_3<-which(copyspalte7tot==lowestrmsedata[3])
        index_lowest_rmse_4<-which(copyspalte7tot==lowestrmsedata[4])
        index_lowest_rmse_5<-which(copyspalte7tot==lowestrmsedata[5])
        finaltabledf3 <- finaltabledf
        hilfsvektor<-c(1:length(spalte7))
        hilfsvektor <- hilfsvektor[!hilfsvektor %in% c(index_lowest_rmse_11,index_lowest_rmse_22,index_lowest_rmse_33,index_lowest_rmse_44,index_lowest_rmse_55)]
        finaltabledf3 <- finaltabledf3[-hilfsvektor,]
        finaltabledf3 <- finaltabledf3[order(finaltabledf3[,7]),]
        names(finaltabledf3)<-c("No.Run","Transformation","Cutoff","Predictors","Mean_Of_Squared_Residuals","%VarExplained","RMSE","NRMSE","NRMSEMAXMIN","R_SQUARED")
        table <- tableGrob(finaltabledf3)
        title <- textGrob(paste0(variables[j]," - predictarea",i," - 5 runs with lowest RMSE"),gp=gpar(fontsize=15))
        footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
        padding <- unit(2,"line")
        table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
        table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
        table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
        #png(file.path(overviewdir,paste0("8_Output_Overview_LowestRMSE_predictarea",i,".png")), units="cm", width=(ncol(table))*4, height=(nrow(table)), res=300) #, height = 200, width = 100)
        png(file.path(overviewdir,paste0("8_Output_Overview_5LowestRMSE_predictarea",i,".png")),width=5200, height=710, res=300)
        grid.newpage() 
        grid.draw(table) 
        dev.off()
        write.csv(finaltabledf3,file.path(overviewdir,paste0("8_Output_Overview_5LowestRMSE_predictarea",i,".csv")))      
        write.xlsx(finaltabledf3,file.path(overviewdir,paste0("8_Output_Overview_5LowestRMSE_predictarea",i,".xlsx")))
        # Copy Runs Output
        bestindex<-vector()
        bestindexwithoutduplicates<-vector()
        bestindex<-c(index_highest_r2_1,index_highest_r2_2,index_highest_r2_3,index_highest_r2_4,index_highest_r2_5,
                     index_lowest_rmse_1,index_lowest_rmse_2,index_lowest_rmse_3,index_lowest_rmse_4,index_lowest_rmse_5)
        bestindexwithoutduplicates<-unique(bestindex)
        #for (b in 1:length(bestindexwithoutduplicates)) {
        #pre <- list.files(paste0(out_dir,"/",predictareas_name[i],"/",var), paste0("No.Run_",bestindexwithoutduplicates[b]) , recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
        #prenew<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output")
        #file.copy(pre,prenew,recursive=TRUE)
        #}
        setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output"))
        for (b in 1:length(bestindexwithoutduplicates)) {
          dir.create(paste0("No_Run_",bestindexwithoutduplicates[b]))
          pre <- list.files(paste0(out_dir,"/",predictareas_name[i],"/",var), paste0("No.Run_",bestindexwithoutduplicates[b]) , recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
          pree <- paste0(pre,"/Final_Models")
          pree <- pree[1]
          prenew<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output","/No_Run_",bestindexwithoutduplicates[b])
          file.copy(pree,prenew,recursive=TRUE)
          pree <- paste0(pre,"/Ordinary_Kriging_Prediction")
          pree <- pree[1]
          file.copy(pree,prenew,recursive=TRUE)
          pree <- paste0(pre,"/Random_Forest_Regression")
          pree <- pree[1]
          file.copy(pree,prenew,recursive=TRUE)
        }
        
        # Overview Ordinary Kriging Output
        setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output"))
        aii<-dir()[file.info(dir())$isdir]
        for (b in 1:length(aii)) {
          rn<-extract_numeric(aii[b])
          pre2 <- list.files(path=paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output/",aii[b],"/Ordinary_Kriging_Prediction"), full.names = TRUE, recursive = TRUE)
          #pre2 <- paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output/",aii[b],"/Ordinary_Kriging_Prediction/","5_Prediction_OrdinaryKriging_",var,"_","predictarea",i,"_",ai[k],"_","run",rn,".png")
          prenew2<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
          file.copy(pre2,prenew2,recursive=TRUE)
        }
        
        zr<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
        filesvector <- list.files(zr)
        for (b in 1:length(aii)) {
          memory.limit(size=9999999999999)
          rn<-extract_numeric(aii[b])
          zr<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
          img1 <- readImage(file.path(list.files(zr,paste0("run",rn,".png"),full.names = TRUE, recursive = TRUE)))
          index_transformation <- which(spalte1==rn)
          transformationhere <- spalte2[index_transformation]
          #img1 <- readImage(file.path(zr,paste0("5_Prediction_OrdinaryKriging_",var,"_predictarea",i,"_",ai[k],"_run",rn,".png")))
          img2 <- img1[1018:2000,1:1000,]
          png(file.path(zr,paste0("No.Run_",rn,".png")),res=300,width = 2000,height = 2000)
          plot(img2) #paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
          title(paste0("No.Run",rn," - ",transformationhere))
          dev.off()
          #assign(paste0("plot",rn),readImage(file.path(list.files(zr,paste0("No.Run_",rn,".png"),full.names = TRUE, recursive = TRUE))))
          #assign(paste0("plot",b),readPNG(file.path(zr,paste0("No.Run_",rn,".png"))))
        }
        
        memory.limit(size=9999999999999)
        index_lowest_rmse_total <- c(index_lowest_rmse_1,index_lowest_rmse_2,index_lowest_rmse_3,index_lowest_rmse_4,index_lowest_rmse_5)
        index_highest_r2_total <- c(index_highest_r2_1,index_highest_r2_2,index_highest_r2_3,index_highest_r2_4,index_highest_r2_5)
        fz<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
        fx<-paste0("5_Overview_Ordinary_Kriging_",var,"_predictarea",i,"_LowestRMSE.png")
        fxx<-paste0("5_Overview_Ordinary_Kriging_",var,"_predictarea",i,"_HighestR2.png")
        plotted1 <- readImage(file.path(paste0(fz,"/No.Run_",index_lowest_rmse_total[1],".png")))
        plotted2 <- readImage(file.path(paste0(fz,"/No.Run_",index_lowest_rmse_total[2],".png")))
        plotted3 <- readImage(file.path(paste0(fz,"/No.Run_",index_lowest_rmse_total[3],".png")))
        plotted4 <- readImage(file.path(paste0(fz,"/No.Run_",index_lowest_rmse_total[4],".png")))
        plotted5 <- readImage(file.path(paste0(fz,"/No.Run_",index_lowest_rmse_total[5],".png")))
        png(file.path(fz,fx),res=300,width =10000, height=2000)
        grid.arrange(rasterGrob(plotted1),rasterGrob(plotted2),rasterGrob(plotted3),rasterGrob(plotted4),rasterGrob(plotted5), nrow=1, ncol=5)
        dev.off()
        plottted1 <- readImage(file.path(paste0(fz,"/No.Run_",index_highest_r2_total[1],".png")))
        plottted2 <- readImage(file.path(paste0(fz,"/No.Run_",index_highest_r2_total[2],".png")))
        plottted3 <- readImage(file.path(paste0(fz,"/No.Run_",index_highest_r2_total[3],".png")))
        plottted4 <- readImage(file.path(paste0(fz,"/No.Run_",index_highest_r2_total[4],".png")))
        plottted5 <- readImage(file.path(paste0(fz,"/No.Run_",index_highest_r2_total[5],".png")))
        png(file.path(fz,fxx),res=300,width =10000, height=2000)
        grid.arrange(rasterGrob(plottted1),rasterGrob(plottted2),rasterGrob(plottted3),rasterGrob(plottted4),rasterGrob(plottted5), nrow=1, ncol=5)
        dev.off()
        
      } else {}
      
      
      # Copy Final models
      memory.limit(size=9999999999999)
      for (b in 1:length(aii)) {
        rn<-extract_numeric(aii[b])
        pre<-list.files(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output/No_Run_",rn,"/Final_Models"), recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
        setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Final_Models"))
        dir.create(paste0("No_Run_",rn))
        prenew<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Final_Models/No_Run_",rn)
        file.copy(pre,prenew,recursive=TRUE)
      }
      
      # Copy Summary Graphics
      memory.limit(size=9999999999999)
      for (b in 1:length(aii)) {
        pre <- list.files(paste0(out_dir,"/",predictareas_name[i],"/",var),paste0("Overview"),recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
        prenew<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Exploratory_Statistic")
        file.copy(pre,prenew,recursive=TRUE)
        pre2 <- list.files(paste0(out_dir,"/",predictareas_name[i],"/",var),paste0("Choose_Transformation_",var,"_",predictareas_name[i],".png"),recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
        prenew2<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Exploratory_Statistic")
        file.copy(pre2,prenew2,recursive=TRUE)
        pre3 <- list.files(paste0(out_dir,"/",predictareas_name[i],"/",var),paste0("Descriptive_Statistic_",var,"_",predictareas_name[i],".png"),recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
        prenew3<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Exploratory_Statistic")
        file.copy(pre3,prenew3,recursive=TRUE)
      }
      
      # Find Best Models (OrdinaryKriging output + Table with 5 runs with lowest RMSE)
      memory.limit(size=9999999999999)
      findout <- paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Find_Best_Model")
      plot1 <- readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Tables/8_Output_Overview_5LowestRMSE_predictarea",i,".png")))
      plot2 <- readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging/",fx)))
      memory.limit(size=9999999999999)
      png(file.path(findout,paste0("8_Find_Best_Modell_",var,"_","predictarea",i,"_LowestRMSE.png")),res=300,width = 4000, height = 2355) 
      grid.arrange(rasterGrob(plot1), rasterGrob(plot2),nrow=2, ncol=1)
      dev.off()
      
      # Find Best Models (OrdinaryKriging output + Table with 5 runs with highest R-SQUARED)
      memory.limit(size=9999999999999)
      findout <- paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Find_Best_Model")
      plot1 <- readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Tables/8_Output_Overview_5HighestR2_predictarea",i,".png")))
      plot2 <- readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging/",fxx)))
      memory.limit(size=9999999999999)
      png(file.path(findout,paste0("8_Find_Best_Modell_",var,"_","predictarea",i,"_HighestR2.png")),res=300,width = 4000, height = 2355) 
      grid.arrange(rasterGrob(plot1), rasterGrob(plot2),nrow=2, ncol=1)
      dev.off()
      
      # Predictors' importance (from the regression part: Random Forest (RF))
      memory.limit(size=9999999999999)
      for (b in 1:length(aii)) {
        rn<-extract_numeric(aii[b])
        pre<-list.files(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output/No_Run_",rn,"/Random_Forest_Regression/4_PredictorsImportance"), recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
        setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/PredictorsInfluence"))
        dir.create(paste0("No_Run_",rn))
        prenew<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/PredictorsInfluence/No_Run_",rn)
        file.copy(pre,prenew,recursive=TRUE)
      }
    }
    
  }
  
  if (length(list.files(paste0(out_dir,"/",predictareas_name[i],"/",variables[j]),"No.Run",recursive=TRUE, full.names=FALSE, include.dirs=TRUE))==0) {
    setwd(paste0(script_dir,"/output_overview/Overview_",variables[j]))
    fileConn<-file("Read_Me.txt")
    writeLines(c(paste0("The source data of the dependent variable ","'",variables[j],"'"," is not good enough."),"Therefore no prediction and output data could be generated.","One cause could be that there are too few available soil data points of this variable."), fileConn)
    close(fileConn)
  } else {
    # Save total overview tables for each dependent variable (lowest RMSE)
    if (nrpredictareas==1){
      fg<-paste0("/output/",predictareas_name[1],"/",variables[j])
      fh<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),fg)
      findout <- paste0(script_dir,"/output_overview/Overview_",variables[j])
      plot00 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_1/",variables[j],"/Overview_Exploratory_Statistic/Descriptive_Statistic_",variables[j],"_",predictareas_name[1],".png"))) #w:1420 h:200
      plot11 <-  readImage(file.path(paste0(fh,"/Choose_Transformation_",variables[j],"_",predictareas_name[1],".png"))) #w:630 h:200
      png(file.path(findout,paste0("Overview_Exp_Des_Statistic_Tot_",variables[j],"_",predictareas_name[1],".png")),width = 1420, height = 200)
      grid.arrange(rasterGrob(plot00),rasterGrob(plot11),nrow=1, ncol=2)
      dev.off()
      plot0 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_",variables[j],"/Overview_Exp_Des_Statistic_Tot_",variables[j],"_",predictareas_name[1],".png")))
      plot1 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_1/",variables[j],"/Find_Best_Model/8_Find_Best_Modell_",variables[j],"_","predictarea1_LowestRMSE.png")))
      img1 <- plot1[,300:810,]
      plot2 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_1/",variables[j],"/Find_Best_Model/8_Find_Best_Modell_",variables[j],"_","predictarea1_LowestRMSE.png")))
      img2 <- plot2[,1380:2185,]
      memory.limit(size=9999999999999)
      findout <- paste0(script_dir,"/output_overview/Overview_",variables[j]) 
      png(file.path(findout,paste0("Total_Overview_",variables[j],"_LowestRMSE.png")),res=300,width = 4000, height = 1515)
      grid.arrange(rasterGrob(plot0),rasterGrob(img1), rasterGrob(img2) ,top = paste0(variables[j]," - OVERVIEW - predictarea 1"),nrow=3, ncol=1)
      dev.off()
      setwd(findout)
      pdf(paste0("Total_Overview_",variables[j],"_LowestRMSE.pdf"),paper="a4r",width = 4000, height = 1515)
      grid.arrange(rasterGrob(plot0),rasterGrob(img1), rasterGrob(img2) ,top = paste0(variables[j]," - OVERVIEW - predictarea 1"),nrow=3, ncol=1)
      dev.off()
      findout2 <- paste0(script_dir,"/output_overview")
      setwd(findout2)
      pdf(paste0("Total_Overview_",variables[j],"_LowestRMSE.pdf"),paper="a4r",width = 4000, height = 1515)
      grid.arrange(rasterGrob(plot0),rasterGrob(img1), rasterGrob(img2) ,top = paste0(variables[j]," - OVERVIEW - predictarea 1"),nrow=3, ncol=1)
      dev.off()
    } else if (nrpredictareas==2) {
      fg<-paste0("/output/",predictareas_name[1],"/",variables[j])
      fh<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),fg)
      findout <- paste0(script_dir,"/output_overview/Overview_",variables[j])
      plot00 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_1/",variables[j],"/Overview_Exploratory_Statistic/Descriptive_Statistic_",variables[j],"_",predictareas_name[1],".png"))) #w:1420 h:200
      plot11 <-  readImage(file.path(paste0(fh,"/Choose_Transformation_",variables[j],"_",predictareas_name[1],".png"))) #w:630 h:200
      png(file.path(findout,paste0("Overview_Exp_Des_Statistic_Tot_",variables[j],"_",predictareas_name[1],".png")),width = 1420, height = 200)
      grid.arrange(rasterGrob(plot00),rasterGrob(plot11),nrow=1, ncol=2)
      dev.off()
      plot0 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_",variables[j],"/Overview_Exp_Des_Statistic_Tot_",variables[j],"_",predictareas_name[1],".png")))
      plot1 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_1/",variables[j],"/Find_Best_Model/8_Find_Best_Modell_",variables[j],"_","predictarea1_LowestRMSE.png")))
      img1 <- plot1[,300:810,]
      plot2 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_1/",variables[j],"/Find_Best_Model/8_Find_Best_Modell_",variables[j],"_","predictarea1_LowestRMSE.png")))
      img2 <- plot2[,1380:2185,]
      memory.limit(size=9999999999999)
      findout <- paste0(script_dir,"/output_overview/Overview_",variables[j]) 
      png(file.path(findout,paste0("Overview_predictarea1_",variables[j],"_LowestRMSE.png")),res=300,width = 4000, height = 1515)
      grid.arrange(rasterGrob(plot0),rasterGrob(img1), rasterGrob(img2) ,top = paste0(variables[j]," - OVERVIEW - predictarea 1"),nrow=3, ncol=1)
      dev.off()
      plot3 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_2/",variables[j],"/Find_Best_Model/8_Find_Best_Modell_",variables[j],"_","predictarea2_LowestRMSE.png")))
      img3 <- plot3[,300:810,]
      plot4 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_2/",variables[j],"/Find_Best_Model/8_Find_Best_Modell_",variables[j],"_","predictarea2_LowestRMSE.png")))
      img4 <- plot4[,1380:2185,]
      memory.limit(size=9999999999999)
      findout <- paste0(script_dir,"/output_overview/Overview_",variables[j]) 
      png(file.path(findout,paste0("Overview_predictarea2_",variables[j],"_LowestRMSE.png")),res=300,width = 4000, height = 1515)
      grid.arrange(rasterGrob(plot0),rasterGrob(img3), rasterGrob(img4) ,top = paste0(variables[j]," - OVERVIEW - predictarea 2"),nrow=3, ncol=1)
      dev.off()
      plot11 <- readImage(file.path(findout,paste0("Overview_predictarea1_",variables[j],"_LowestRMSE.png")))
      plot22 <- readImage(file.path(findout,paste0("Overview_predictarea2_",variables[j],"_LowestRMSE.png")))
      memory.limit(size=9999999999999)
      png(file.path(findout,paste0("Total_Overview_",variables[j],"_LowestRMSE.png")),res=300,width = 4000, height = 3030)
      grid.arrange(rasterGrob(plot11), rasterGrob(plot22),top = paste0(variables[j]," - TOTAL OVERVIEW"),nrow=2, ncol=1)
      dev.off()
      setwd(findout)
      pdf(paste0("Total_Overview_",variables[j],"_LowestRMSE.pdf"),paper="a4r",width = 4000, height = 3030)
      grid.arrange(rasterGrob(plot11), rasterGrob(plot22),top = paste0(variables[j]," - TOTAL OVERVIEW"),nrow=2, ncol=1)
      dev.off()
      findout2 <- paste0(script_dir,"/output_overview")
      setwd(findout2)
      pdf(paste0("Total_Overview_",variables[j],"_LowestRMSE.pdf"),paper="a4r",width = 4000, height = 3030)
      grid.arrange(rasterGrob(plot11), rasterGrob(plot22),top = paste0(variables[j]," - TOTAL OVERVIEW"),nrow=2, ncol=1)
      dev.off()
    } else {}
    
    
    # Save total overview tables for each dependent variable (highest R-Squared)
    if (nrpredictareas==1){
      fg<-paste0("/output/",predictareas_name[1],"/",variables[j])
      fh<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),fg)
      findout <- paste0(script_dir,"/output_overview/Overview_",variables[j])
      plot00 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_1/",variables[j],"/Overview_Exploratory_Statistic/Descriptive_Statistic_",variables[j],"_",predictareas_name[1],".png"))) #w:1420 h:200
      plot11 <-  readImage(file.path(paste0(fh,"/Choose_Transformation_",variables[j],"_",predictareas_name[1],".png"))) #w:630 h:200
      png(file.path(findout,paste0("Overview_Exp_Des_Statistic_Tot_",variables[j],"_",predictareas_name[1],".png")),width = 1420, height = 200)
      grid.arrange(rasterGrob(plot00),rasterGrob(plot11),nrow=1, ncol=2)
      dev.off()
      plot0 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_",variables[j],"/Overview_Exp_Des_Statistic_Tot_",variables[j],"_",predictareas_name[1],".png")))
      plot1 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_1/",variables[j],"/Find_Best_Model/8_Find_Best_Modell_",variables[j],"_","predictarea1_HighestR2.png")))
      img1 <- plot1[,300:810,]
      plot2 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_1/",variables[j],"/Find_Best_Model/8_Find_Best_Modell_",variables[j],"_","predictarea1_HighestR2.png")))
      img2 <- plot2[,1380:2185,]
      memory.limit(size=9999999999999)
      findout <- paste0(script_dir,"/output_overview/Overview_",variables[j]) 
      png(file.path(findout,paste0("Total_Overview_",variables[j],"_HighestR2.png")),res=300,width = 4000, height = 1515)
      grid.arrange(rasterGrob(plot0),rasterGrob(img1), rasterGrob(img2) ,top = paste0(variables[j]," - OVERVIEW - predictarea 1"),nrow=3, ncol=1)
      dev.off()
      setwd(findout)
      pdf(paste0("Total_Overview_",variables[j],"_HighestR2.pdf"),paper="a4r",width = 4000, height = 1515)
      grid.arrange(rasterGrob(plot0),rasterGrob(img1), rasterGrob(img2) ,top = paste0(variables[j]," - OVERVIEW - predictarea 1"),nrow=3, ncol=1)
      dev.off()
      findout2 <- paste0(script_dir,"/output_overview")
      setwd(findout2)
      pdf(paste0("Total_Overview_",variables[j],"_HighestR2.pdf"),paper="a4r",width = 4000, height = 1515)
      grid.arrange(rasterGrob(plot0),rasterGrob(img1), rasterGrob(img2) ,top = paste0(variables[j]," - OVERVIEW - predictarea 1"),nrow=3, ncol=1)
      dev.off()
    } else if (nrpredictareas==2) {
      fg<-paste0("/output/",predictareas_name[1],"/",variables[j])
      fh<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),fg)
      findout <- paste0(script_dir,"/output_overview/Overview_",variables[j])
      plot00 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_1/",variables[j],"/Overview_Exploratory_Statistic/Descriptive_Statistic_",variables[j],"_",predictareas_name[1],".png"))) #w:1420 h:200
      plot11 <-  readImage(file.path(paste0(fh,"/Choose_Transformation_",variables[j],"_",predictareas_name[1],".png"))) #w:630 h:200
      png(file.path(findout,paste0("Overview_Exp_Des_Statistic_Tot_",variables[j],"_",predictareas_name[1],".png")),width = 1420, height = 200)
      grid.arrange(rasterGrob(plot00),rasterGrob(plot11),nrow=1, ncol=2)
      dev.off()
      plot0 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_",variables[j],"/Overview_Exp_Des_Statistic_Tot_",variables[j],"_",predictareas_name[1],".png")))
      plot1 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_1/",variables[j],"/Find_Best_Model/8_Find_Best_Modell_",variables[j],"_","predictarea1_HighestR2.png")))
      img1 <- plot1[,300:810,]
      plot2 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_1/",variables[j],"/Find_Best_Model/8_Find_Best_Modell_",variables[j],"_","predictarea1_HighestR2.png")))
      img2 <- plot2[,1380:2185,]
      memory.limit(size=9999999999999)
      findout <- paste0(script_dir,"/output_overview/Overview_",variables[j]) 
      png(file.path(findout,paste0("Overview_predictarea1_",variables[j],"_HighestR2.png")),res=300,width = 4000, height = 1515)
      grid.arrange(rasterGrob(plot0),rasterGrob(img1), rasterGrob(img2) ,top = paste0(variables[j]," - OVERVIEW - predictarea 1"),nrow=3, ncol=1)
      dev.off()
      plot3 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_2/",variables[j],"/Find_Best_Model/8_Find_Best_Modell_",variables[j],"_","predictarea2_HighestR2.png")))
      img3 <- plot3[,300:810,]
      plot4 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_2/",variables[j],"/Find_Best_Model/8_Find_Best_Modell_",variables[j],"_","predictarea2_HighestR2.png")))
      img4 <- plot4[,1380:2185,]
      memory.limit(size=9999999999999)
      findout <- paste0(script_dir,"/output_overview/Overview_",variables[j]) 
      png(file.path(findout,paste0("Overview_predictarea2_",variables[j],"_HighestR2.png")),res=300,width = 4000, height = 1515)
      grid.arrange(rasterGrob(plot0),rasterGrob(img3), rasterGrob(img4) ,top = paste0(variables[j]," - OVERVIEW - predictarea 2"),nrow=3, ncol=1)
      dev.off()
      plot11 <- readImage(file.path(findout,paste0("Overview_predictarea1_",variables[j],"_HighestR2.png")))
      plot22 <- readImage(file.path(findout,paste0("Overview_predictarea2_",variables[j],"_HighestR2.png")))
      memory.limit(size=9999999999999)
      png(file.path(findout,paste0("Total_Overview_",variables[j],"_HighestR2.png")),res=300,width = 4000, height = 3030)
      grid.arrange(rasterGrob(plot11), rasterGrob(plot22),top = paste0(variables[j]," - TOTAL OVERVIEW"),nrow=2, ncol=1)
      dev.off()
      setwd(findout)
      pdf(paste0("Total_Overview_",variables[j],"_HighestR2.pdf"),paper="a4r",width = 4000, height = 3030)
      grid.arrange(rasterGrob(plot11), rasterGrob(plot22),top = paste0(variables[j]," - TOTAL OVERVIEW"),nrow=2, ncol=1)
      dev.off()
      findout2 <- paste0(script_dir,"/output_overview")
      setwd(findout2)
      pdf(paste0("Total_Overview_",variables[j],"_HighestR2.pdf"),paper="a4r",width = 4000, height = 3030)
      grid.arrange(rasterGrob(plot11), rasterGrob(plot22),top = paste0(variables[j]," - TOTAL OVERVIEW"),nrow=2, ncol=1)
      dev.off()
    } else {}
  }
  
}

for (j in 1:lenvariables){
  var<-variables[j]
  for (i in 1:nrpredictareas){
    if (i==1) {
      overviewdir<-paste0(script_dir,"/output_overview/Overview_Predict_area_1/",var,"/Overview_Tables")
    } else if (i==2) {
      overviewdir<-paste0(script_dir,"/output_overview/Overview_Predict_area_2/",var,"/Overview_Tables/")
    } else {}
    filepath <- paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Tables/")
    filepath2 <-paste0(script_dir,"/output/",predictareas_name[i],"/",var,"/")
    som_1 <- read.csv(file=paste0(filepath,"8_Output_Overview_Total_predictarea",i,".csv"),header=TRUE)
    tot_1 <- read.csv(file=paste0(filepath2,"Descriptive_Statistic_",var,"_",predictareas_name[i],".csv"),header=TRUE)
    noruns <- vector()
    results <- vector()
    jc<-1
    while (jc<=5){
      result2 <- 999999
      norun <-999
      indexlinenorun<-1
      for (ii in 1:length(som_1[,8])) { 
        print(ii)
        if (som_1[ii,3]=="raw") {
          diviso <- tot_1[1,8] - tot_1[1,3]
        } else if (som_1[ii,3]=="log") {
          diviso <- tot_1[2,8] - tot_1[2,3]
        } else if (som_1[ii,3]=="sqrt") {
          diviso <- tot_1[3,8] - tot_1[3,3]
        } else if (som_1[ii,3]=="quad") {
          diviso <- tot_1[4,8] - tot_1[4,3]
        } else if (som_1[ii,3]=="inverse") {
          diviso <- tot_1[5,8] - tot_1[5,3]
        } else {}
        min_RMSE <- min(som_1[,8])
        max_RMSE <- max(som_1[,8])
        max_R2 <- max(som_1[,11])
        min_R2 <- min(som_1[,11])
        yy <- abs(som_1[ii,8]-min_RMSE) / diviso   # Distanz zum min RMSE
        y <- yy
        zz <- abs(som_1[ii,11]-max_R2)   # Distanz zum max R2
        z <- zz
        result <- y+z
        if (result < result2) {
          norun <- som_1[ii,2]
          indexlinenorun <- ii
          result2<-y+z
        } else {}
      }
      sdf<-length(noruns)
      noruns <- append(noruns,norun,after=sdf)
      som_1 <- som_1[-indexlinenorun,]
      jc<-jc+1
      sdf2<-length(results)
      results <- append(results,result2,after=sdf2)
    }
    som_2 <- read.csv(file=paste0(filepath,"8_Output_Overview_Total_predictarea",i,".csv"),header=TRUE)
    hilfsvektor<-som_2[,2]
    hilfsvektor <- hilfsvektor[!hilfsvektor %in% noruns]
    som_2 <- som_2[-hilfsvektor,]
    som_3 <- data.frame()
    spalt <- som_2[,2]
    rank1 <- match(noruns[1],spalt)
    rank2 <- match(noruns[2],spalt)
    rank3 <- match(noruns[3],spalt)
    rank4 <- match(noruns[4],spalt)
    rank5 <- match(noruns[5],spalt)
    som_3 <- rbind(som_2[rank5,],som_3)
    som_3 <- rbind(som_2[rank4,],som_3)
    som_3 <- rbind(som_2[rank3,],som_3)
    som_3 <- rbind(som_2[rank2,],som_3)
    som_3 <- rbind(som_2[rank1,],som_3)
    som_3 <- som_3[,-1]
    # Save Table with 5 runs with best combination of R2 and R-Squared
    table <- tableGrob(som_3)
    title <- textGrob(paste0(var," - predictarea",i," - 5 runs with best RMSE and R-Squared combi"),gp=gpar(fontsize=15))
    footnote <- textGrob("", x=0, hjust=0, gp=gpar(fontface="italic")) 
    padding <- unit(2,"line")
    table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0) 
    table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding) 
    table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))
    png(file.path(overviewdir,paste0("8_Output_Overview_5RMSER2Combi_predictarea",i,".png")), width=5200, height=710, res=300)
    grid.newpage() 
    grid.draw(table) 
    dev.off()
    write.csv(som_3,file.path(overviewdir,paste0("8_Output_Overview_5RMSER2Combi_predictarea_predictarea",i,".csv")))      
    write.xlsx(som_3,file.path(overviewdir,paste0("8_Output_Overview_5RMSER2Combi_predictarea_predictarea",i,".xlsx")))
    
    # Copy Runs Output
    setwd(script_dir)
    dir.create("BestCombi")
    out_dir <- paste0(script_dir,"/output")
    setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output"))
    for (b in 1:length(noruns)) {
      dir.create(paste0("No_Run_",noruns[b]))
      pre <- list.files(paste0(out_dir,"/",predictareas_name[i],"/",var), paste0("No.Run_",noruns[b]) , recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
      pree <- paste0(pre,"/Final_Models")
      pree <- pree[1]
      prenew<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output","/No_Run_",noruns[b])
      file.copy(pree,prenew,recursive=TRUE)
      pree <- paste0(pre,"/Ordinary_Kriging_Prediction")
      pree <- pree[1]
      file.copy(pree,prenew,recursive=TRUE)
      pree <- paste0(pre,"/Random_Forest_Regression")
      pree <- pree[1]
      file.copy(pree,prenew,recursive=TRUE)
    }
    
    # Overview Ordinary Kriging Output
    aii<-c(paste0("No_Run_",noruns[1]),paste0("No_Run_",noruns[2]),paste0("No_Run_",noruns[3]),paste0("No_Run_",noruns[4]),paste0("No_Run_",noruns[5]))
    for (b in 1:length(aii)) {
      rn<-extract_numeric(aii[b])
      pre2 <- list.files(path=paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output/",aii[b],"/Ordinary_Kriging_Prediction"), full.names = TRUE, recursive = TRUE)
      prenew2<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
      file.copy(pre2,prenew2,recursive=TRUE)
    }
    zr<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
    filesvector <- list.files(zr)
    for (b in 1:length(aii)) {
      memory.limit(size=9999999999999)
      rn<-extract_numeric(aii[b])
      zr<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
      img1 <- readImage(file.path(list.files(zr,paste0("run",rn,".png"),full.names = TRUE, recursive = TRUE)))
      index_transformation <- which(som_3[,1]==rn)
      transformationhere <- som_3[index_transformation,2]
      img2 <- img1[1018:2000,1:1000,]
      png(file.path(zr,paste0("No.Run_",rn,".png")),res=300,width = 2000,height = 2000)
      plot(img2) #paste0(var,"_","predictarea",i,"_",ai[k],"_","run",nrun,"_cutoff_",taia))
      title(paste0("No.Run",rn," - ",transformationhere))
      dev.off()
    }
    memory.limit(size=9999999999999)
    fz<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging")
    fx<-paste0("5_Overview_Ordinary_Kriging_",var,"_predictarea",i,"_BestCombiRMSER2.png")
    plotted1 <- readImage(file.path(paste0(fz,"/No.Run_",noruns[1],".png")))
    plotted2 <- readImage(file.path(paste0(fz,"/No.Run_",noruns[2],".png")))
    plotted3 <- readImage(file.path(paste0(fz,"/No.Run_",noruns[3],".png")))
    plotted4 <- readImage(file.path(paste0(fz,"/No.Run_",noruns[4],".png")))
    plotted5 <- readImage(file.path(paste0(fz,"/No.Run_",noruns[5],".png")))
    png(file.path(fz,fx),res=300,width =10000, height=2000)
    grid.arrange(rasterGrob(plotted1),rasterGrob(plotted2),rasterGrob(plotted3),rasterGrob(plotted4),rasterGrob(plotted5), nrow=1, ncol=5)
    dev.off()
    
    # Copy Final models
    memory.limit(size=9999999999999)
    for (b in 1:length(aii)) {
      rn<-extract_numeric(aii[b])
      pre<-list.files(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output/No_Run_",rn,"/Final_Models"), recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
      setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Final_Models"))
      dir.create(paste0("No_Run_",rn))
      prenew<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Final_Models/No_Run_",rn)
      file.copy(pre,prenew,recursive=TRUE)
    }
    
    # Find Best Models (OrdinaryKriging output + Table with 5 runs with lowest RMSE)
    memory.limit(size=9999999999999)
    findout <- paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Find_Best_Model")
    plot1 <- readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Tables/8_Output_Overview_5RMSER2Combi_predictarea",i,".png")))
    plot2 <- readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Overview_Prediction_Ordinary_Kriging/",fx)))
    memory.limit(size=9999999999999)
    png(file.path(findout,paste0("8_Find_Best_Modell_",var,"_","predictarea",i,"_bestRMSER2Combi.png")),res=300,width = 4000, height = 2355) 
    grid.arrange(rasterGrob(plot1), rasterGrob(plot2),nrow=2, ncol=1)
    dev.off()
    
    # Predictors' importance (from the regression part: Random Forest (RF))
    memory.limit(size=9999999999999)
    for (b in 1:length(aii)) {
      rn<-extract_numeric(aii[b])
      pre<-list.files(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/Runs_Output/No_Run_",rn,"/Random_Forest_Regression/4_PredictorsImportance"), recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
      setwd(paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/PredictorsInfluence"))
      dir.create(paste0("No_Run_",rn))
      prenew<-paste0(script_dir,"/output_overview/Overview_Predict_area_",i,"/",var,"/PredictorsInfluence/No_Run_",rn)
      file.copy(pre,prenew,recursive=TRUE)
    }
    
  }
  
  # Save total overview tables for each dependent variable (best RMSE & R-Squared combination)
  if (nrpredictareas==1){
    fg<-paste0("/output/",predictareas_name[1],"/",var)
    fh<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),fg)
    findout <- paste0(script_dir,"/output_overview/Overview_",var)
    plot0 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_",var,"/Overview_Exp_Des_Statistic_Tot_",var,"_",predictareas_name[1],".png")))
    plot1 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_1/",var,"/Find_Best_Model/8_Find_Best_Modell_",var,"_","predictarea1_bestRMSER2Combi.png")))
    img1 <- plot1[,300:810,]
    plot2 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_1/",var,"/Find_Best_Model/8_Find_Best_Modell_",var,"_","predictarea1_bestRMSER2Combi.png")))
    img2 <- plot2[,1380:2185,]
    memory.limit(size=9999999999999)
    findout <- paste0(script_dir,"/output_overview/Overview_",var) 
    png(file.path(findout,paste0("Total_Overview_",var,"_bestRMSER2Combi.png")),res=300,width = 4000, height = 1515)
    grid.arrange(rasterGrob(plot0),rasterGrob(img1), rasterGrob(img2) ,top = paste0(var," - OVERVIEW - predictarea 1"),nrow=3, ncol=1)
    dev.off()
    setwd(findout)
    pdf(paste0("Total_Overview_",var,"_bestRMSER2Combi.pdf"),paper="a4r",width = 4000, height = 1515)
    grid.arrange(rasterGrob(plot0),rasterGrob(img1), rasterGrob(img2) ,top = paste0(var," - OVERVIEW - predictarea 1"),nrow=3, ncol=1)
    dev.off()
    findout2 <- paste0(script_dir,"/output_overview")
    setwd(findout2)
    pdf(paste0("Total_Overview_",var,"_bestRMSER2Combi.pdf"),paper="a4r",width = 4000, height = 1515)
    grid.arrange(rasterGrob(plot0),rasterGrob(img1), rasterGrob(img2) ,top = paste0(var," - OVERVIEW - predictarea 1"),nrow=3, ncol=1)
    dev.off()
  } else if (nrpredictareas==2) {
    fg<-paste0("/output/",predictareas_name[1],"/",var)
    fh<-paste0(dirname(rstudioapi::getActiveDocumentContext()$path),fg)
    findout <- paste0(script_dir,"/output_overview/Overview_",var)
    plot0 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_",var,"/Overview_Exp_Des_Statistic_Tot_",var,"_",predictareas_name[1],".png")))
    plot1 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_1/",var,"/Find_Best_Model/8_Find_Best_Modell_",var,"_","predictarea1_bestRMSER2Combi.png")))
    img1 <- plot1[,300:810,]
    plot2 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_1/",var,"/Find_Best_Model/8_Find_Best_Modell_",var,"_","predictarea1_bestRMSER2Combi.png")))
    img2 <- plot2[,1380:2185,]
    memory.limit(size=9999999999999)
    findout <- paste0(script_dir,"/output_overview/Overview_",var) 
    png(file.path(findout,paste0("Overview_predictarea1_",var,"_bestRMSER2Combi.png")),res=300,width = 4000, height = 1515)
    grid.arrange(rasterGrob(plot0),rasterGrob(img1), rasterGrob(img2) ,top = paste0(var," - OVERVIEW - predictarea 1"),nrow=3, ncol=1)
    dev.off()
    plot3 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_2/",var,"/Find_Best_Model/8_Find_Best_Modell_",var,"_","predictarea2_bestRMSER2Combi.png")))
    img3 <- plot3[,300:810,]
    plot4 <-  readImage(file.path(paste0(script_dir,"/output_overview/Overview_Predict_area_2/",var,"/Find_Best_Model/8_Find_Best_Modell_",var,"_","predictarea2_bestRMSER2Combi.png")))
    img4 <- plot4[,1380:2185,]
    memory.limit(size=9999999999999)
    findout <- paste0(script_dir,"/output_overview/Overview_",var) 
    png(file.path(findout,paste0("Overview_predictarea2_",var,"_bestRMSER2Combi.png")),res=300,width = 4000, height = 1515)
    grid.arrange(rasterGrob(plot0),rasterGrob(img3), rasterGrob(img4) ,top = paste0(var," - OVERVIEW - predictarea 2"),nrow=3, ncol=1)
    dev.off()
    plot11 <- readImage(file.path(findout,paste0("Overview_predictarea1_",var,"_bestRMSER2Combi.png")))
    plot22 <- readImage(file.path(findout,paste0("Overview_predictarea2_",var,"_bestRMSER2Combi.png")))
    memory.limit(size=9999999999999)
    png(file.path(findout,paste0("Total_Overview_",var,"_bestRMSER2Combi.png")),res=300,width = 4000, height = 3030)
    grid.arrange(rasterGrob(plot11), rasterGrob(plot22),top = paste0(var," - TOTAL OVERVIEW"),nrow=2, ncol=1)
    dev.off()
    setwd(findout)
    pdf(paste0("Total_Overview_",var,"_bestRMSER2Combi.pdf"),paper="a4r",width = 4000, height = 3030)
    grid.arrange(rasterGrob(plot11), rasterGrob(plot22),top = paste0(var," - TOTAL OVERVIEW"),nrow=2, ncol=1)
    dev.off()
    findout2 <- paste0(script_dir,"/output_overview")
    setwd(findout2)
    pdf(paste0("Total_Overview_",var,"_bestRMSER2Combi.pdf"),paper="a4r",width = 4000, height = 3030)
    grid.arrange(rasterGrob(plot11), rasterGrob(plot22),top = paste0(var," - TOTAL OVERVIEW"),nrow=2, ncol=1)
    dev.off()
  } else {}
  
  
}

findout2<-paste0(script_dir,"/output_overview")
setwd(findout2)
mergevector <- list.files (findout2, "*.pdf", full.names = FALSE)
pdf_combine(mergevector,output="Total_Overview_Print.pdf")
file.remove(mergevector)