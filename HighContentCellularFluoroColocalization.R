library(data.table)
library(stringr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(rstatix)
library(plyr)
library(gridExtra)
library(viridis)

##Functions
sem <- function(x) sd(x) / sqrt(length(x))

##View the data
data_cells<- fread("V:/Laboratoires/Lécuyer/Juan-Carlos/HCS/20230802_MCF7_Ab_DMSOGWFTY-2/20231006_CellProfiler-Analysis/MyExpt_Cells.csv")
data_cyto<- fread("V:/Laboratoires/Lécuyer/Juan-Carlos/HCS/20230802_MCF7_Ab_DMSOGWFTY-2/20231006_CellProfiler-Analysis/MyExpt_Cytoplasm.csv")
data_image<- fread("V:/Laboratoires/Lécuyer/Juan-Carlos/HCS/20230802_MCF7_Ab_DMSOGWFTY-2/20231006_CellProfiler-Analysis/MyExpt_Image.csv")
data_nuclei<- fread("V:/Laboratoires/Lécuyer/Juan-Carlos/HCS/20230802_MCF7_Ab_DMSOGWFTY-2/20231006_CellProfiler-Analysis/MyExpt_Nuclei.csv")

#Read data
data_cells<- fread("V:/Laboratoires/Lécuyer/Juan-Carlos/HCS/20230802_MCF7_Ab_DMSOGWFTY-2/20231006_CellProfiler-Analysis/MyExpt_Nuclei.csv", 
             select = c("ImageNumber",
                        "Metadata_Plate",
                        "Metadata_Site",
                        "Metadata_Well",
                        "Intensity_IntegratedIntensity_Cy5",
                        "Intensity_IntegratedIntensity_GFP",
                        "Intensity_MeanIntensity_Cy5",
                        "Intensity_MeanIntensity_GFP",
                        "Intensity_MedianIntensity_Cy5",
                        "Intensity_MedianIntensity_GFP",
                        "Correlation_Correlation_Cy5_GFP",      
                        "Correlation_Costes_Cy5_GFP",           
                        "Correlation_Costes_GFP_Cy5",           
                        "Correlation_K_Cy5_GFP",                
                        "Correlation_K_GFP_Cy5",                
                        "Correlation_Manders_Cy5_GFP",          
                        "Correlation_Manders_GFP_Cy5",          
                        "Correlation_Overlap_Cy5_GFP",          
                        "Correlation_RWC_Cy5_GFP",              
                        "Correlation_RWC_GFP_Cy5",
                        "RadialDistribution_FracAtD_Cy5_1of3",  
                        "RadialDistribution_FracAtD_Cy5_2of3",  
                        "RadialDistribution_FracAtD_Cy5_3of3",  
                        "RadialDistribution_FracAtD_GFP_1of3",  
                        "RadialDistribution_FracAtD_GFP_2of3",  
                        "RadialDistribution_FracAtD_GFP_3of3",  
                        "RadialDistribution_MeanFrac_Cy5_1of3", 
                        "RadialDistribution_MeanFrac_Cy5_2of3", 
                        "RadialDistribution_MeanFrac_Cy5_3of3", 
                        "RadialDistribution_MeanFrac_GFP_1of3", 
                        "RadialDistribution_MeanFrac_GFP_2of3", 
                        "RadialDistribution_MeanFrac_GFP_3of3", 
                        "RadialDistribution_RadialCV_Cy5_1of3", 
                        "RadialDistribution_RadialCV_Cy5_2of3", 
                        "RadialDistribution_RadialCV_Cy5_3of3", 
                        "RadialDistribution_RadialCV_GFP_1of3", 
                        "RadialDistribution_RadialCV_GFP_2of3", 
                        "RadialDistribution_RadialCV_GFP_3of3"))
#remove data from H02
#data_cells<-data_cells[which(!data_cells$Metadata_Well == "H02")]
#data_cells<-data_cells[which(!data_cells$Metadata_Well == "H01")]

# Replace multiple strings at a time for condition
rep_str_condition = c("DMSO"="DMSO", "GW"="GW4869", "FTY"="FTY720")

data_cells$Metadata_Plate <- str_replace_all(data_cells$Metadata_Plate, rep_str_condition)

# Replace multiple strings at a time for GFP
rep_str_gfp = c("A01"="NSM", "A02"="NSM", "A03"="NSM", "A04"="NSM", "A05"="NSM", "A06"="NSM", "A07"="NSM", "A08"="NSM", "A09"="NSM", "A10"="NSM", "A11"="NSM", "A12"="NSM",
                "B01"="ASM", "B02"="ASM", "B03"="ASM", "B04"="ASM", "B05"="ASM", "B06"="ASM", "B07"="ASM", "B08"="ASM", "B09"="ASM", "B10"="ASM", "B11"="ASM", "B12"="ASM",
                "D01"="hnRNPA2B1", "D02"="hnRNPA2B1", "D03"="hnRNPA2B1", "D04"="hnRNPA2B1", "D05"="hnRNPA2B1", "D06"="hnRNPA2B1", "D07"="hnRNPA2B1", "D08"="hnRNPA2B1", "D09"="hnRNPA2B1", "D10"="hnRNPA2B1", "D11"="hnRNPA2B1", "D12"="hnRNPA2B1",
                "E01"="HNRNPD", "E02"="HNRNPD", "E03"="HNRNPD", "E04"="HNRNPD", "E05"="HNRNPD", "E06"="HNRNPD", "E07"="HNRNPD", "E08"="HNRNPD", "E09"="HNRNPD", "E10"="HNRNPD", "E11"="HNRNPD", "E12"="HNRNPD",
                "F01"="G3BP2", "F02"="G3BP2", "F03"="G3BP2", "F04"="G3BP2", "F05"="G3BP2", "F06"="G3BP2", "F07"="G3BP2", "F08"="G3BP2", "F09"="G3BP2", "F10"="G3BP2", "F11"="G3BP2", "F12"="G3BP2",
                "G01"="ANXA2", "G02"="ANXA2","G03"="ANXA2","G04"="ANXA2","G05"="ANXA2","G06"="ANXA2","G07"="ANXA2","G08"="ANXA2","G09"="ANXA2","G10"="ANXA2","G11"="ANXA2", "G12"="ANXA2",
                "H01"="PTBP1", "H02"="PTBP1","H03"="PTBP1","H04"="PTBP1","H05"="PTBP1","H06"="PTBP1","H07"="PTBP1","H08"="PTBP1","H09"="PTBP1","H10"="PTBP1","H11"="PTBP1", "H12"="PTBP1",
                "C01"="TSG101", "C02"="TSG101", "C03"="TSG101", "C04"="TSG101", "C05"="TSG101", "C06"="TSG101", "C07"="TSG101", "C08"="TSG101", "C09"="TSG101", "C10"="TSG101", "C11"="TSG101", "C12"="TSG101"
                )

data_cells$gfp <- str_replace_all(data_cells$Metadata_Well, rep_str_gfp)

# Replace multiple strings at a time for Cy5
rep_str_cy5 = c("A01"="Ceramide", "A02"="Ceramide", "A03"="CD63", "A04"="CD63", "A05"="EEA1", "A06"="EEA1", "A07"="LC3B", "A08"="LC3B", "A09"="CD9", "A10"="CD9", "A11"="LAMP1", "A12"="LAMP1",
                "B01"="Ceramide", "B02"="Ceramide", "B03"="CD63", "B04"="CD63", "B05"="EEA1", "B06"="EEA1", "B07"="LC3B", "B08"="LC3B", "B09"="CD9", "B10"="CD9", "B11"="LAMP1", "B12"="LAMP1",
                "C01"="Ceramide", "C02"="Ceramide", "C03"="CD63", "C04"="CD63", "C05"="EEA1", "C06"="EEA1", "C07"="LC3B", "C08"="LC3B", "C09"="CD9", "C10"="CD9", "C11"="LAMP1", "C12"="LAMP1",
                "D01"="Ceramide", "D02"="Ceramide", "D03"="CD63", "D04"="CD63", "D05"="EEA1", "D06"="EEA1", "D07"="LC3B", "D08"="LC3B", "D09"="CD9", "D10"="CD9", "D11"="LAMP1", "D12"="LAMP1",
                "E01"="Ceramide", "E02"="Ceramide", "E03"="CD63", "E04"="CD63", "E05"="EEA1", "E06"="EEA1", "E07"="LC3B", "E08"="LC3B", "E09"="CD9", "E10"="CD9", "E11"="LAMP1", "E12"="LAMP1",
                "F01"="Ceramide", "F02"="Ceramide", "F03"="CD63", "F04"="CD63", "F05"="EEA1", "F06"="EEA1", "F07"="LC3B", "F08"="LC3B", "F09"="CD9", "F10"="CD9", "F11"="LAMP1", "F12"="LAMP1",
                "G01"="Ceramide", "G02"="Ceramide", "G03"="CD63", "G04"="CD63", "G05"="EEA1", "G06"="EEA1", "G07"="LC3B", "G08"="LC3B", "G09"="CD9", "G10"="CD9", "G11"="LAMP1", "G12"="LAMP1",
                "H01"="Ceramide", "H02"="Ceramide", "H03"="CD63", "H04"="CD63", "H05"="EEA1", "H06"="EEA1", "H07"="LC3B", "H08"="LC3B", "H09"="CD9", "H10"="CD9", "H11"="LAMP1", "H12"="LAMP1")                
data_cells$cy5 <- str_replace_all(data_cells$Metadata_Well, rep_str_cy5)

##creating identifiers
data_cells$plate_Gfp_cy5 <- paste(data_cells$Metadata_Plate,"_",data_cells$gfp,"_",data_cells$cy5, sep = "")
data_cells$Gfp_cy5 <- paste(data_cells$gfp,"_",data_cells$cy5, sep = "")

#making the identifier as factors
#data_cells$plate_Gfp_cy5<-as.factor(data_cells$plate_Gfp_cy5)
#data_cells$Gfp_cy5<-as.factor(data_cells$Gfp_cy5)

##removing the outlier
remove_outliers <- function(x, low, upp) {
  lower_bound <- quantile(x, low, na.rm = T)
  upper_bound <- quantile(x, upp, na.rm = T)
  x[x < lower_bound | x > upper_bound] <- NA
  #x_filtered <- x[x >= lower_bound & x <= upper_bound]
  return(x)
}

#remove outliers from intergrated or mean intensity
data_cells$Intensity_IntegratedIntensity_Cy5 <- remove_outliers(data_cells$Intensity_IntegratedIntensity_Cy5, 0.1,0.9)
data_cells$Intensity_IntegratedIntensity_GFP <- remove_outliers(data_cells$Intensity_IntegratedIntensity_GFP, 0.1,0.9)
data_cells$Intensity_MeanIntensity_Cy5 <- remove_outliers(data_cells$Intensity_MeanIntensity_Cy5,0.1,0.9)
data_cells$Intensity_MeanIntensity_GFP<- remove_outliers(data_cells$Intensity_MeanIntensity_GFP,0.1,0.9)

#normalizing the dataset Normalized value = (x – x bar) / sd
##where:
### x = data value
### x bar = mean of dataset
### sd = standard deviation of dataset



#### Loops for plotting####


#### Loop for plotting normalized integrated intensity GFP in cell images####

pdf(file = "Normalized integrated intensity GFP_Boxplot.pdf",onefile = TRUE)

for(i in 1:length(unique(data_cells$gfp))){
  fitdata <- data_cells[data_cells$gfp == unique(data_cells$gfp)[i],]
  fitdata[["Metadata_Plate"]] = factor(  fitdata[["Metadata_Plate"]], levels = c("DMSO", "GW4869", "FTY720"))
  ##Calculate the Mean Intensity for Each Condition:
  mean_IntegratedIntensity <- aggregate(Intensity_IntegratedIntensity_GFP ~ Metadata_Plate, fitdata, mean)
  ##Calculate the Fold Change Relative to DMSO:
  mean_FC_DMSO <- mean_IntegratedIntensity$Intensity_IntegratedIntensity_GFP[mean_IntegratedIntensity$Metadata_Plate == "DMSO"]
  fitdata$normalized_DMSO_MeanIntegratedIntensity_GFP<-fitdata$Intensity_IntegratedIntensity_GFP / mean_FC_DMSO

  fit <- fitdata %>% 
    anova_test(normalized_DMSO_MeanIntegratedIntensity_GFP ~ Metadata_Plate) %>% 
    add_significance()
  ##fit2 <- fitdata %>% 
    ##kruskal_test(normalized_DMSO_MeanIntegratedIntensity_GFP ~ Metadata_Plate) %>% 
    ##add_significance()
  
  #### Run Dunn ###
  ##dunn <- fitdata %>% 
    ##dunn_test(normalized_DMSO_MeanIntegratedIntensity_GFP ~ Metadata_Plate) %>% 
    ##add_significance() %>% 
    ##add_xy_position()
  
  #### Run Tukey ###
  tukey <- fitdata %>% 
    tukey_hsd(normalized_DMSO_MeanIntegratedIntensity_GFP ~ Metadata_Plate) %>% 
    add_significance() %>% 
    add_xy_position()
  
  
  A<-ggboxplot(fitdata,
              "Metadata_Plate",
              "normalized_DMSO_MeanIntegratedIntensity_GFP", color = "#000004",
              fill="Metadata_Plate",palette = c("#000004", "#3B0F70", "#8C2981"),
              width = 0.2, ##control the width of the boxes
              size = 0.5, ##controls the size of the outlines and points
              #ylim = c(0, 3),
              title = c(unique(fitdata$gfp)))+
    scale_x_discrete(name = "Condition") +
    scale_y_continuous(name = "Normalized Integrated Intensity") +
    stat_pvalue_manual(tukey,
                       hide.ns = F, y.position = tukey$y.position+0.05, step.increase = 0.05, label = "p.adj.signif", tip.length = 0.01, bracket.nudge.y = 0.0)+
    theme(aspect.ratio=9/3, axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = .9), legend.position="right", legend.title = element_blank(),axis.line=element_line(size=0.5))+
    ##TO add mean label on top
    ##stat_summary(fun.data = function(x) data.frame(y=(max(tukey$y.position)+1), label = paste("Mean=",round(mean(x),2))), geom="text")+
    stat_summary(geom = "crossbar", width=1, fatten=0, color="white", 
                 fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))})
    ##labs(subtitle = get_test_label(fit,detailed = TRUE), caption = get_pwc_label(dunn))
  grid.arrange(A)
}
dev.off()

#### Loop for plotting normalized integrated intensity Cy5 in cell images####

pdf(file = "Normalized integrated intensity Cy5_Boxplot.pdf",onefile = TRUE)

for(i in 1:length(unique(data_cells$cy5))){
  fitdata <- data_cells[data_cells$cy5 == unique(data_cells$cy5)[i],]
  fitdata[["Metadata_Plate"]] = factor(  fitdata[["Metadata_Plate"]], levels = c("DMSO", "GW4869", "FTY720"))
  ##Calculate the Mean Intensity for Each Condition:
  mean_IntegratedIntensity <- aggregate(Intensity_IntegratedIntensity_Cy5 ~ Metadata_Plate, fitdata, mean)
  ##Calculate the Fold Change Relative to DMSO:
  mean_FC_DMSO <- mean_IntegratedIntensity$Intensity_IntegratedIntensity_Cy5[mean_IntegratedIntensity$Metadata_Plate == "DMSO"]
  fitdata$normalized_DMSO_MeanIntegratedIntensity_Cy5<-fitdata$Intensity_IntegratedIntensity_Cy5 / mean_FC_DMSO
  
  fit <- fitdata %>% 
    anova_test(normalized_DMSO_MeanIntegratedIntensity_Cy5 ~ Metadata_Plate) %>% 
    add_significance()
  ##fit2 <- fitdata %>% 
  ##kruskal_test(normalized_DMSO_MeanIntegratedIntensity_Cy5 ~ Metadata_Plate) %>% 
  ##add_significance()
  
  #### Run Dunn ###
  ##dunn <- fitdata %>% 
  ##dunn_test(normalized_DMSO_MeanIntegratedIntensity_Cy5 ~ Metadata_Plate) %>% 
  ##add_significance() %>% 
  ##add_xy_position()
  
  #### Run Tukey ###
  tukey <- fitdata %>% 
    tukey_hsd(normalized_DMSO_MeanIntegratedIntensity_Cy5 ~ Metadata_Plate) %>% 
    add_significance() %>% 
    add_xy_position()
  
  
  A<-ggboxplot(fitdata,
               "Metadata_Plate",
               "normalized_DMSO_MeanIntegratedIntensity_Cy5", color = "#000004",
               fill="Metadata_Plate",palette = c("#000004", "#3B0F70", "#8C2981"),
               #order = c("DMSO", "GW", "FTY"),
               width = 0.2, ##control the width of the boxes
               size = 0.5, ##controls the size of the outlines and points
               #ylim = c(0, 3),
               title = c(unique(fitdata$cy5)))+
    scale_x_discrete(name = "Condition") +
    scale_y_continuous(name = "Normalized Integrated Intensity") +
    stat_pvalue_manual(tukey,
                       hide.ns = F, y.position = tukey$y.position+0.05, step.increase = 0.05, label = "p.adj.signif", tip.length = 0.01, bracket.nudge.y = 0.0)+
    theme(aspect.ratio=9/3, axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = .9), legend.position="right", legend.title = element_blank(),axis.line=element_line(size=0.5))+
    ##TO add mean label on top
    ##stat_summary(fun.data = function(x) data.frame(y=(max(tukey$y.position)+1), label = paste("Mean=",round(mean(x),2))), geom="text")+
    stat_summary(geom = "crossbar", width=1, fatten=0, color="white", 
                 fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))})
  ##labs(subtitle = get_test_label(fit,detailed = TRUE), caption = get_pwc_label(dunn))
  grid.arrange(A)
}
dev.off()

#### Loop for plotting normalized mean intensity GFP in cell images####

pdf(file = "Normalized mean intensity GFP_Boxplot.pdf",onefile = TRUE)

for(i in 1:length(unique(data_cells$gfp))){
  fitdata <- data_cells[data_cells$gfp == unique(data_cells$gfp)[i],]
  fitdata[["Metadata_Plate"]] = factor(  fitdata[["Metadata_Plate"]], levels = c("DMSO", "GW4869", "FTY720"))
  ##Calculate the Mean Intensity for Each Condition:
  mean_MeanIntensity <- aggregate(Intensity_MeanIntensity_GFP ~ Metadata_Plate, fitdata, mean)
  ##Calculate the Fold Change Relative to DMSO:
  mean_FC_DMSO <- mean_MeanIntensity$Intensity_MeanIntensity_GFP[mean_MeanIntensity$Metadata_Plate == "DMSO"]
  fitdata$normalized_DMSO_MeanMeanIntensity_GFP<-fitdata$Intensity_MeanIntensity_GFP / mean_FC_DMSO
  
  fit <- fitdata %>% 
    anova_test(normalized_DMSO_MeanMeanIntensity_GFP~ Metadata_Plate) %>% 
    add_significance()
  ##fit2 <- fitdata %>% 
  ##kruskal_test(normalized_DMSO_MeanIntensity_GFP ~ Metadata_Plate) %>% 
  ##add_significance()
  
  #### Run Dunn ###
  ##dunn <- fitdata %>% 
  ##dunn_test(normalized_DMSO_MeanIntensity_GFP ~ Metadata_Plate) %>% 
  ##add_significance() %>% 
  ##add_xy_position()
  
  #### Run Tukey ###
  tukey <- fitdata %>% 
    tukey_hsd(normalized_DMSO_MeanMeanIntensity_GFP ~ Metadata_Plate) %>% 
    add_significance() %>% 
    add_xy_position()
  
  
  A<-ggboxplot(fitdata,
               "Metadata_Plate",
               "normalized_DMSO_MeanMeanIntensity_GFP", color = "#000004",
               fill="Metadata_Plate",palette = c("#000004", "#3B0F70", "#8C2981"),
               width = 0.2, ##control the width of the boxes
               size = 0.5, ##controls the size of the outlines and points
               #ylim = c(0, 3),
               title = c(unique(fitdata$gfp)))+
    scale_x_discrete(name = "Condition") +
    scale_y_continuous(name = "Normalized Mean Intensity") +
    stat_pvalue_manual(tukey,
                       hide.ns = F, y.position = tukey$y.position+0.05, step.increase = 0.05, label = "p.adj.signif", tip.length = 0.01, bracket.nudge.y = 0.0)+
    theme(aspect.ratio=9/3, axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = .9), legend.position="right", legend.title = element_blank(),axis.line=element_line(size=0.5))+
    ##TO add mean label on top
    ##stat_summary(fun.data = function(x) data.frame(y=(max(tukey$y.position)+1), label = paste("Mean=",round(mean(x),2))), geom="text")+
    stat_summary(geom = "crossbar", width=1, fatten=0, color="white", 
                 fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))})
  ##labs(subtitle = get_test_label(fit,detailed = TRUE), caption = get_pwc_label(dunn))
  grid.arrange(A)
}
dev.off()

#### Loop for plotting normalized mean intensity Cy5 in cell images####

pdf(file = "Normalized mean intensity Cy5_Boxplot.pdf",onefile = TRUE)

for(i in 1:length(unique(data_cells$cy5))){
  fitdata <- data_cells[data_cells$cy5 == unique(data_cells$cy5)[i],]
  fitdata[["Metadata_Plate"]] = factor(  fitdata[["Metadata_Plate"]], levels = c("DMSO", "GW4869", "FTY720"))
  ##Calculate the Mean Intensity for Each Condition:
  mean_MeanIntensity <- aggregate(Intensity_MeanIntensity_Cy5 ~ Metadata_Plate, fitdata, mean)
  ##Calculate the Fold Change Relative to DMSO:
  mean_FC_DMSO <- mean_MeanIntensity$Intensity_MeanIntensity_Cy5[mean_MeanIntensity$Metadata_Plate == "DMSO"]
  fitdata$normalized_DMSO_MeanMeanIntensity_Cy5<-fitdata$Intensity_MeanIntensity_Cy5 / mean_FC_DMSO
  
  fit <- fitdata %>% 
    anova_test(normalized_DMSO_MeanMeanIntensity_Cy5 ~ Metadata_Plate) %>% 
    add_significance()
  ##fit2 <- fitdata %>% 
  ##kruskal_test(normalized_DMSO_MeanMeanIntensity_Cy5 ~ Metadata_Plate) %>% 
  ##add_significance()
  
  #### Run Dunn ###
  ##dunn <- fitdata %>% 
  ##dunn_test(normalized_DMSO_MeanMeanIntensity_Cy5 ~ Metadata_Plate) %>% 
  ##add_significance() %>% 
  ##add_xy_position()
  
  #### Run Tukey ###
  tukey <- fitdata %>% 
    tukey_hsd(normalized_DMSO_MeanMeanIntensity_Cy5 ~ Metadata_Plate) %>% 
    add_significance() %>% 
    add_xy_position()
  
  
  A<-ggboxplot(fitdata,
               "Metadata_Plate",
               "normalized_DMSO_MeanMeanIntensity_Cy5", color = "#000004",
               fill="Metadata_Plate",palette = c("#000004", "#3B0F70", "#8C2981"),
               #order = c("DMSO", "GW", "FTY"),
               width = 0.2, ##control the width of the boxes
               size = 0.5, ##controls the size of the outlines and points
               #ylim = c(0, 3),
               title = c(unique(fitdata$cy5)))+
    scale_x_discrete(name = "Condition") +
    scale_y_continuous(name = "Normalized Mean Intensity") +
    stat_pvalue_manual(tukey,
                       hide.ns = F, y.position = tukey$y.position+0.05, step.increase = 0.05, label = "p.adj.signif", tip.length = 0.01, bracket.nudge.y = 0.0)+
    theme(aspect.ratio=9/3, 
          axis.line=element_line(size=0.2),
          axis.ticks=element_line(size=0.2),
          legend.box.background = element_rect(color = NA),
          axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = .9), legend.position="right", legend.title = element_blank(),axis.line=element_line(size=0.5))+
    ##TO add mean label on top
    ##stat_summary(fun.data = function(x) data.frame(y=(max(tukey$y.position)+1), label = paste("Mean=",round(mean(x),2))), geom="text")+
    stat_summary(geom = "crossbar", width=1, fatten=0, color="white", 
                 fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))})
  ##labs(subtitle = get_test_label(fit,detailed = TRUE), caption = get_pwc_label(dunn))
  grid.arrange(A)
}
dev.off()


#### Loop for plotting correlation pearsons in cell images####

corColNames <- c(       "Correlation_Correlation_Cy5_GFP",      
                        "Correlation_Costes_Cy5_GFP",           
                        "Correlation_Costes_GFP_Cy5",           
                        "Correlation_K_Cy5_GFP",                
                        "Correlation_K_GFP_Cy5",                
                        "Correlation_Manders_Cy5_GFP",          
                        "Correlation_Manders_GFP_Cy5",          
                        "Correlation_Overlap_Cy5_GFP",          
                        "Correlation_RWC_Cy5_GFP",              
                        "Correlation_RWC_GFP_Cy5")


for(m in 1:length(corColNames)){
  
  pdf(file = paste0(corColNames[m],".pdf"),onefile = TRUE)
  
  for(i in 1:length(unique(data_cells$Gfp_cy5))){
    fitdata <- data_cells[data_cells$Gfp_cy5 == unique(data_cells$Gfp_cy5)[i],]
    fitdata[["Metadata_Plate"]] = factor(  fitdata[["Metadata_Plate"]], levels = c("DMSO", "GW4869", "FTY720"))
    fitdata[[corColNames[m]]] <- remove_outliers(fitdata[[corColNames[m]]], 0.025, 0.975)
    
    fit <- fitdata %>% 
      kruskal_test(as.formula(paste(corColNames[m],"Metadata_Plate", sep=" ~ "))) %>% 
      add_significance()
    
    #### Run Dunn ###
    if(any(is.na(tryCatch(fitdata %>% 
               dunn_test(as.formula(paste(corColNames[m],"Metadata_Plate", sep=" ~ "))) %>% 
               add_significance() %>% 
               add_xy_position(), 
             error=function(e) NA)))== TRUE){
      print(paste(unique(data_cells$Gfp_cy5)[i],"failed for",corColNames[m]))
      fitdata <- data_cells[data_cells$Gfp_cy5 == unique(data_cells$Gfp_cy5)[i],]
      fitdata[["Metadata_Plate"]] = factor(  fitdata[["Metadata_Plate"]], levels = c("DMSO", "GW4869", "FTY720"))
      dunn <- fitdata %>% 
        dunn_test(as.formula(paste(corColNames[m],"Metadata_Plate", sep=" ~ "))) %>% 
        add_significance() %>% 
        add_xy_position()
      }else {
        dunn <- fitdata %>% 
      dunn_test(as.formula(paste(corColNames[m],"Metadata_Plate", sep=" ~ "))) %>% 
      add_significance() %>% 
      add_xy_position()
        }
    
    A<-ggviolin(fitdata,
                 "Metadata_Plate",
                 corColNames[m], color = "#000004",
                 fill="Metadata_Plate",palette = c("#000004", "#3B0F70", "#8C2981"),
                 #order = c("DMSO", "GW", "FTY"),
                 width = 0.7, ##control the width of the boxes
                 size = 0.5, ##controls the size of the outlines and points
                 #ylim = c(0, 1.2),
                 title = c(unique(fitdata$Gfp_cy5)))+
      scale_x_discrete(name = "Condition") +
      scale_y_continuous(name = corColNames[m]) +
      stat_pvalue_manual(dunn,
                         hide.ns = F, y.position = dunn$y.position+0.2, step.increase = 0.05, label = "p.adj.signif", tip.length = 0.01, bracket.nudge.y = 0.0)+
      theme(aspect.ratio=9/3, 
            axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = .9), 
            legend.position="right", 
            legend.title = element_blank(),
            axis.line=element_line(size=0.2),
            axis.ticks=element_line(size=0.2),
            legend.box.background = element_rect(color = NA))+
      ##TO add mean label on top
      ##stat_summary(fun.data = function(x) data.frame(y=(max(dunn$y.position)+1), label = paste("Mean=",round(mean(x),2))), geom="text")+
      stat_summary(geom = "crossbar", width=1, fatten=0, color="white", 
                   fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))})
    ##labs(subtitle = get_test_label(fit,detailed = TRUE), caption = get_pwc_label(dunn))
    grid.arrange(A)
  }
  dev.off()
}  
