# This script is for deciding on best matched internal standards in quality control (QC) samples run on TSQ for Global Metabolite F. cylindrus experiment 

# Load required packages into environment
library(ggplot2)
library(dplyr)
library(raster)
library(MetBrewer)
library(gridExtra)
library(tidyr)

# set working directory 
setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/Life/School/MSc/1_CHII/FC_GM_Data/TSQ_Frag_Metab/TSQ_Frag_Metab_25042022/TSQ_Frag_Metab_250422_RawData")


# load QC data and only include QC treatments
QC.Data <- read.csv("ER3_149_Frag_Catalina_01_Cal03022022_output.csv") %>% filter( grepl("QC",Replicate.Name))


# Plot QC peaks per compound
QC.Data$Replicate.Name <- factor(QC.Data$Replicate.Name, levels = c("ER3_149_2fdQC_03", "ER3_149_2fdQC_04" , "ER3_149_2fdQC_08" , "ER3_149_2fdQC_09", "ER3_149_2fdQC_44" , "ER3_149_2fdQC_79" ,
 "ER3_149_2fdQC_114", "ER3_149_2fdQC_142"))

# Fix DMB name 
QC.Data <- QC.Data %>% 
  mutate(Molecule.Name = replace(Molecule.Name, Molecule.Name == "Dimethyl-benzimidazole (DMB)", "DMB"))

QC.Peaks <- ggplot() + 
  geom_bar(data = QC.Data, aes(y = Total.Area, x = Replicate.Name), stat = "identity") + 
  ggtitle("QC Peak Areas (No Normalization)") + 
  ylab("Peak Area") +
  xlab("QC Injection") +
  facet_wrap( ~ Molecule.Name, scales = "free") +
  theme(axis.text.x = element_text(angle = 90,   hjust = 1)) + 
  scale_x_discrete(breaks = levels(QC.Data$Replicate.Name), 
                   labels=c("03", "04", "08", "09", "44", "79", "114", "142" ))
  

# Get mean, standard deviation, and cv of all QC's per compound
QC.Data.SumStats <- QC.Data %>%
  group_by(Molecule.Name) %>%
  dplyr::summarise(Mean.Peak.Area = mean(Total.Area), SD.Peak.Area = sd(Total.Area), CV.Peak.Area = raster::cv(Total.Area))



# Get Normalized Peaks (divide each light QC by each heavy QC) ------------


# Divide QC peaks by heavy B1, B2, CN-B12, B7
BMIS.Data <- QC.Data %>% filter(!grepl('heavy', Molecule.Name))


# Get Heavy B1 Normalized peaks (divide peaks by heavy B1)
Heavy.B1.Peaks <-  QC.Data %>% filter(Molecule.Name == "B1-heavy")
BMIS.Data$Heavy.B1.Peaks <- Heavy.B1.Peaks$Total.Area
BMIS.Data$Heavy.B1.Norm <- BMIS.Data$Total.Area /BMIS.Data$Heavy.B1.Peaks

# Get heavy B2 normalized peaks (divide peaks by heavy B2)
Heavy.B2.Peaks <-  QC.Data %>% filter(Molecule.Name == "B2-heavy")
BMIS.Data$Heavy.B2.Peaks <- Heavy.B2.Peaks$Total.Area
BMIS.Data$Heavy.B2.Norm <- BMIS.Data$Total.Area /BMIS.Data$Heavy.B2.Peaks

# Get heavy CN-B12 normalized peaks (divide peaks by heavy CN-B12)
Heavy.CNB12.Peaks <-  QC.Data %>% filter(Molecule.Name == "B12-CN-heavy")
BMIS.Data$Heavy.CNB12.Peaks <- Heavy.CNB12.Peaks$Total.Area
BMIS.Data$Heavy.CNB12.Norm <- BMIS.Data$Total.Area /BMIS.Data$Heavy.CNB12.Peaks

# Get heavy B7 normalized peaks (divide peaks by heavy CN-B12)
Heavy.B7.Peaks <-  QC.Data %>% filter(Molecule.Name == "B12-CN-heavy")
BMIS.Data$Heavy.B7.Peaks <- Heavy.B7.Peaks$Total.Area
BMIS.Data$Heavy.B7.Norm <- BMIS.Data$Total.Area /BMIS.Data$Heavy.B7.Peaks



# Compare Heavy-Normalized QC to non-normed QC ----------------------------
HeavyNormQC.Data.SumStats <- BMIS.Data %>%
  group_by(Molecule.Name) %>%
  dplyr::summarise(Mean.Peak.Area = mean(Total.Area), SD.Peak.Area = sd(Total.Area), CV.Peak.Area = raster::cv(Total.Area), 
            Mean.Peak.Area.B1.Norm = mean(Heavy.B1.Norm),  SD.Peak.Area.B1.Norm = sd(Heavy.B1.Norm), CV.Peak.Area.B1.Norm = raster::cv(Heavy.B1.Norm), 
            Mean.Peak.Area.B2.Norm = mean(Heavy.B2.Norm),  SD.Peak.Area.B2.Norm = sd(Heavy.B2.Norm), CV.Peak.Area.B2.Norm = raster::cv(Heavy.B2.Norm), 
            Mean.Peak.Area.CNB12.Norm = mean(Heavy.CNB12.Norm),  SD.Peak.Area.CNB12.Norm = sd(Heavy.CNB12.Norm), CV.Peak.Area.CNB12.Norm = raster::cv(Heavy.CNB12.Norm), 
            Mean.Peak.Area.B7.Norm = mean(Heavy.B7.Norm),  SD.Peak.Area.B7.Norm = sd(Heavy.B7.Norm), CV.Peak.Area.B7.Norm = raster::cv(Heavy.B7.Norm))


# Compare CV's per compound and normalization type -------------------------------
HeavyNormQC.Data.SumStats1 <- HeavyNormQC.Data.SumStats %>% dplyr::select(Molecule.Name, CV.Peak.Area, CV.Peak.Area.B1.Norm, CV.Peak.Area.B2.Norm, CV.Peak.Area.CNB12.Norm, CV.Peak.Area.B7.Norm)

colnames(HeavyNormQC.Data.SumStats1)[2:6] <- c("No Normalization", "Heavy B1", "Heavy B2", "Heavy CN-B12", "Heavy B7")


gather_cols <- c("No Normalization", "Heavy B1", "Heavy B2", "Heavy CN-B12", "Heavy B7")


# change to long format for plotting
HeavyNormQC.Data.Comp <- gather(HeavyNormQC.Data.SumStats1, Norm.Compound, Norm.Peak.CV,  gather_cols, factor_key = TRUE)

ggplot() + 
  geom_bar(data = HeavyNormQC.Data.Comp, aes(y = Norm.Peak.CV, x = Norm.Compound, fill = Norm.Compound), stat = "identity") + 
  ylab("CV") +
  xlab("Normalization Compound") +
  facet_wrap( ~ Molecule.Name, scales = "fixed") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 15)) +
  scale_fill_manual(values=met.brewer("Cross", 5), name = "Normalization Compound")



# Choose BMIS For Each Light Compound -------------------------------------

# Calculate improvement (or decrease) in CV by normalization
# subtract each normalization CV from non-normalized cv
HeavyNormQC.Data.SumStats$B1 <- HeavyNormQC.Data.SumStats$CV.Peak.Area - HeavyNormQC.Data.SumStats$CV.Peak.Area.B1.Norm 
HeavyNormQC.Data.SumStats$B2 <- HeavyNormQC.Data.SumStats$CV.Peak.Area - HeavyNormQC.Data.SumStats$CV.Peak.Area.B2.Norm 
HeavyNormQC.Data.SumStats$CNB12 <- HeavyNormQC.Data.SumStats$CV.Peak.Area - HeavyNormQC.Data.SumStats$CV.Peak.Area.CNB12.Norm 
HeavyNormQC.Data.SumStats$B7 <- HeavyNormQC.Data.SumStats$CV.Peak.Area - HeavyNormQC.Data.SumStats$CV.Peak.Area.B7.Norm 

# create a list of highest changes and which heavy standard norm led to them
HeavyNormQC.Data.SumStats.deltaCV <- HeavyNormQC.Data.SumStats %>% dplyr::select(B1, B2, CNB12, B7)
BMIS_list <- data.frame(colnames(HeavyNormQC.Data.SumStats.deltaCV)[apply(HeavyNormQC.Data.SumStats.deltaCV,1,which.max)])
colnames(BMIS_list) <- "BMIS"
BMIS_list$deltaCV <- apply(HeavyNormQC.Data.SumStats.deltaCV, 1, max)

# Add back to dataframe
HeavyNormQC.Data.BMIS <- cbind(HeavyNormQC.Data.SumStats, BMIS_list)

# Change column nam
HeavyNormQC.Data.BMIS$BMIS_col_name <- paste("Mean.Peak.Area.",HeavyNormQC.Data.BMIS$BMIS,".Norm", sep = "")




# Make loop that sees if max change is larger than 30%
# If so, make final norm column show the normalization that caused improvement 
# Start indexing vector at 7 to ignore B1, the B12's, and B2 for BMIS selection (they have corresponding internal standards). B7 not included because B7:heavy B7 norm increases cv rather than reduces

# Make df with compounds with corresponding internal standards
HeavyNormQC.Data.BMIS.corresp <- HeavyNormQC.Data.BMIS[1:6,]


for (i in 1:nrow(HeavyNormQC.Data.BMIS.corresp)){
  HeavyNormQC.Data.BMIS.corresp$Final.BMIS.Norm.Peak[i] <- HeavyNormQC.Data.BMIS.corresp[i,as.character(HeavyNormQC.Data.BMIS.corresp$BMIS_col_name[i])]
  HeavyNormQC.Data.BMIS.corresp$BMIS_used[[i]] <- HeavyNormQC.Data.BMIS.corresp$BMIS[i]
}


# Make df for other compounds 
HeavyNormQC.Data.BMIS.NOcorresp <- HeavyNormQC.Data.BMIS[7:26,]

# Do loop to choose BMIS
for (i in 1:nrow(HeavyNormQC.Data.BMIS.NOcorresp)){
  if (HeavyNormQC.Data.BMIS.NOcorresp$deltaCV[i] < 30){
    HeavyNormQC.Data.BMIS.NOcorresp$Final.BMIS.Norm.Peak[i] <- HeavyNormQC.Data.BMIS.NOcorresp$Mean.Peak.Area[i]
    HeavyNormQC.Data.BMIS.NOcorresp$BMIS_used[i] <- "none"
  } else {
    HeavyNormQC.Data.BMIS.NOcorresp$Final.BMIS.Norm.Peak[i] <- HeavyNormQC.Data.BMIS.NOcorresp[i,as.character(HeavyNormQC.Data.BMIS.NOcorresp$BMIS_col_name[i])]
    HeavyNormQC.Data.BMIS.NOcorresp$BMIS_used[[i]] <- HeavyNormQC.Data.BMIS.NOcorresp$BMIS[i]
  }
}

HeavyNormQC.Data.BMIS <- rbind(HeavyNormQC.Data.BMIS.corresp, HeavyNormQC.Data.BMIS.NOcorresp)


QC_Norm_Export <- HeavyNormQC.Data.BMIS %>% dplyr::select(Molecule.Name, Mean.Peak.Area, SD.Peak.Area, CV.Peak.Area, 
                                                          Mean.Peak.Area.B1.Norm, SD.Peak.Area.B1.Norm, CV.Peak.Area.B1.Norm, 
                                                          Mean.Peak.Area.B2.Norm, SD.Peak.Area.B2.Norm, CV.Peak.Area.B2.Norm,
                                                          Mean.Peak.Area.CNB12.Norm, SD.Peak.Area.CNB12.Norm, CV.Peak.Area.CNB12.Norm,
                                                          Mean.Peak.Area.B7.Norm, SD.Peak.Area.B7.Norm, CV.Peak.Area.B7.Norm,BMIS, deltaCV, BMIS_used, Final.BMIS.Norm.Peak)

QC_Norm_Export_Sum <- QC_Norm_Export %>% dplyr::select(Molecule.Name, BMIS, deltaCV, BMIS_used, Final.BMIS.Norm.Peak)

# Fix weird encoding error
QC_Norm_Export <- apply(QC_Norm_Export,2,as.character)
QC_Norm_Export_Sum <- apply(QC_Norm_Export_Sum,2,as.character)


# Change working directory to move to BMIS folder
setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/Life/School/MSc/1_CHII/FC_GM_Data/TSQ_Frag_Metab/TSQ_Frag_Metab_25042022/TSQ_Frag_Metab_250422_BMIS")

# Export results and summary from BMIS analysis
write.csv(QC_Norm_Export, file = "QC_BMIS_results.csv")
write.csv(QC_Norm_Export_Sum, file = "QC_BMIS_results_sum.csv")

# Create a table from export
grid.table(QC_Norm_Export_Sum)

# Plot QC's before and after normalization for those that get a normalization

# Get before and after CV's
QC_norm_comp_df <- HeavyNormQC.Data.BMIS %>% filter(BMIS_used != "none") %>% dplyr::select(Molecule.Name, CV.Peak.Area, deltaCV) 

# get columns for before and after normalization
QC_norm_comp_df$After_norm <- QC_norm_comp_df$CV.Peak.Area - QC_norm_comp_df$deltaCV

# change to long form data with pivot (note: use of amend gather above)
colnames(QC_norm_comp_df)[c(2,4)] <- c("Before Normalization", "After Normalization") 

QC_norm_comp_df <- QC_norm_comp_df %>% dplyr::select(Molecule.Name, "Before Normalization", "After Normalization")

HeavyNormQC.Data.Comp <- QC_norm_comp_df %>%
  pivot_longer(!Molecule.Name, names_to = "Norm", values_to = "CV")

# filter out methinonine (cv too high)
HeavyNormQC.Data.Comp<- HeavyNormQC.Data.Comp %>% dplyr::filter(Molecule.Name != "Methionine")


HeavyNormQC.Data.Comp$Norm <- factor(HeavyNormQC.Data.Comp$Norm, levels = c("Before Normalization", "After Normalization") )

# Open a pdf file
pdf("Frag_BMIS.pdf") 

# make plot comparing effects of BMIS normalization
ggplot(HeavyNormQC.Data.Comp, aes(fill=Norm, y=CV, x=Molecule.Name)) + 
  geom_bar(position="dodge", stat="identity") +
  xlab("Molecule Name") +
  ylab("Quality Control CV") +
  scale_fill_manual(values=met.brewer("Cross", 2), name = "Normalization")  + 
  theme_classic() +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90))

# Close the pdf file
dev.off() 





  
  
