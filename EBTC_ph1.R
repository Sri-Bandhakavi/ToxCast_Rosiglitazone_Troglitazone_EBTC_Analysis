
library(ggplot2)
#library(ggthemes)
library(tidyr)
library(ComplexHeatmap)
#source("https://bioconductor.org/biocLite.R")
#biocLite("ComplexHeatmap")
getwd()

setwd('C:/Users/Sri Bandhakavi/Desktop/EBTC/EBTC_Phase1_3')
getwd()

#############################################################################################################################
#######1.Load datasets#########################################################################################################
#############################################################################################################################
AssayEndpointResults <-read.csv('Tox21_assay_endpoint_results.csv') 
AssayInformation <-read.csv('Tox21_assay_information.csv')

colnames(AssayEndpointResults)

####################################################################################################################################
########2.Compare two drugs - Rosiglitazone Maleate and Toglitazone for common cell-based assays and positive results#################
####################################################################################################################################
#extract data for two drugs - Troglitazone and Rosiglitazone Maleate into one common data frame
RosMalTrog<-data.frame(AssayEndpointResults$aeid, AssayEndpointResults$Set.2...cpd2, AssayEndpointResults$Set.2...cpd3) 
colnames(RosMalTrog)[1]<- 'AssayID'
colnames(RosMalTrog)[2] <- 'RosiglitazoneMaleate'
colnames(RosMalTrog)[3] <- 'Troglitazone'

#Remove any entries with NA values once each for Rosiglitazone Maleate and then for Troglitazone
RosiglitazoneMal_RT <- na.omit(data.frame(RosMalTrog$AssayID, RosMalTrog$RosiglitazoneMaleate)) # RT suffix in dataframe title indicates use of dataframe for Rosiglitazone and Troglitazone
Troglitazone_RT <- na.omit(data.frame(RosMalTrog$AssayID, RosMalTrog$Troglitazone))

#Identify "common tests" for both Rosiglitazone and Rosiglitazone Maleate
Commontests_RMT<-merge(RosiglitazoneMal_RT, Troglitazone_RT, by = "RosMalTrog.AssayID", all.x = TRUE)

#parse "common tests" into individual results for each drug and rename columns for convenience
RosiglitazoneMal_RMT_Tests<-data.frame(Commontests_RMT$RosMalTrog.AssayID, Commontests_RMT$RosMalTrog.RosiglitazoneMaleate)
colnames(RosiglitazoneMal_RMT_Tests)[1] <- 'AssayID'
colnames(RosiglitazoneMal_RMT_Tests)[2] <- 'Results'

Troglitazone_RMT_Tests <-data.frame(Commontests_RMT$RosMalTrog.AssayID, Commontests_RMT$RosMalTrog.Troglitazone)
colnames(Troglitazone_RMT_Tests)[1] <- 'AssayID'
colnames(Troglitazone_RMT_Tests)[2] <- 'Results'

# Select for each drug (Rosiglitazone Maleate and Troglitazone), "positive" tests/assay IDs (based on result < 1000000) and visualize total count of postiive tests/drug
RosiglitazoneMal_RMT_PosTests<- RosiglitazoneMal_RMT_Tests[which(RosiglitazoneMal_RMT_Tests$Results < 1.000000e+06), ]
Troglitazone_RMT_PosTests<- Troglitazone_RMT_Tests[which(Troglitazone_RMT_Tests$Results < 1.000000e+06), ]

#Summary table of positive tests for each drug & generate bar plot
PosTests_RMT_Table <- data.frame(Drug = factor(c("RosiglitazoneMaleate", "Troglitazone")), PosTests = c(nrow(RosiglitazoneMal_RMT_PosTests), nrow(Troglitazone_RMT_PosTests)))
ggplot(data=PosTests_RMT_Table, aes(x=Drug, y=PosTests, fill = PosTests)) + geom_bar(stat="identity") + guides(fill = FALSE) + geom_text(label=PosTests_RMT_Table$PosTests, vjust=1.6, color = "white", size = 3.5)


#Merge positive tests from each drug into new dataframe & write to excel
PosRosMalTrog <- merge(RosiglitazoneMal_RMT_PosTests, Troglitazone_RMT_PosTests, by='AssayID', all=TRUE)
colnames(PosRosMalTrog)[1]<- 'AssayID'
colnames(PosRosMalTrog)[2]<- 'RosiglitazoneMaleate Positive'
colnames(PosRosMalTrog)[3]<- 'Troglitazone Positive' 

write.csv(PosRosMalTrog, "PosTestsRosMalTrog.csv", row.names = FALSE)

#Extract class I and Class II "tests" from PosRosMalTrog dataframe
# Class I tests = "positive in Rosiglitazone Maleate"
# Class II tests = "positive in Troglitazone only"

PosRosmalTrog <-PosRosMalTrog[!is.na(PosRosMalTrog$`RosiglitazoneMaleate Positive`) & !is.na(PosRosMalTrog$`Troglitazone Positive`), ]

ClassI <-PosRosMalTrog[!is.na(PosRosMalTrog$`RosiglitazoneMaleate Positive`) & !is.na(PosRosMalTrog$`Troglitazone Positive`) | !is.na(PosRosMalTrog$`RosiglitazoneMaleate Positive`), ]
ClassI$Class <- 'Class 1'

ClassII <-PosRosMalTrog[is.na(PosRosMalTrog$`RosiglitazoneMaleate Positive`) & !is.na(PosRosMalTrog$`Troglitazone Positive`), ]

ClassII$Class <- 'Class 2'

#merge AssayInformation table to ClassI assays, extract data of interest, rename columns and save as ClassIAssays2 (cleaned dataframe for further analysis)

ClassIAssays<-merge(ClassI, AssayInformation, by.x='AssayID', by.y ='aeid', all.x=TRUE)

ClassIAssays2 <-data.frame(ClassIAssays$AssayID, ClassIAssays$`RosiglitazoneMaleate Positive`, 
                           ClassIAssays$`Troglitazone Positive`, ClassIAssays$assay_component_name, 
                           ClassIAssays$biological_process_target, ClassIAssays$intended_target_family,
                           ClassIAssays$intended_target_family_sub, ClassIAssays$intended_target_official_full_name, 
                           ClassIAssays$intended_target_gene_name)


colnames(ClassIAssays2)[1] <- 'AssayID'
colnames(ClassIAssays2)[2] <- 'EC50(Rosiglitazone Maleate)' 
colnames(ClassIAssays2)[3] <- 'EC50(Troglidazone)'        
colnames(ClassIAssays2)[4] <- 'Assay Component Name'
colnames(ClassIAssays2)[5] <- 'Biological Process'
colnames(ClassIAssays2)[6] <- 'Target Family'
colnames(ClassIAssays2)[7] <- 'Target SubFamily'
colnames(ClassIAssays2)[8] <- 'Target Full Name'
colnames(ClassIAssays2)[9] <- 'Target Gene Name'
ClassIAssays2$Class <- 'Class 1'

write.csv(ClassIAssays2, "ClassIAssays.csv", row.names = FALSE)
  
#merge AssayInformation table to ClassII assays, extract data of interest, rename columns, and save as ClassIIAssays2 (cleaned dataframe for further analysis)

ClassIIAssays<-merge(ClassII, AssayInformation, by.x='AssayID', by.y ='aeid', all.x=TRUE)

ClassIIAssays2 <-data.frame(ClassIIAssays$AssayID, ClassIIAssays$`RosiglitazoneMaleate Positive`, 
                           ClassIIAssays$`Troglitazone Positive`, ClassIIAssays$assay_component_name, 
                           ClassIIAssays$biological_process_target, ClassIIAssays$intended_target_family,
                           ClassIIAssays$intended_target_family_sub, ClassIIAssays$intended_target_official_full_name, 
                           ClassIIAssays$intended_target_gene_name)


colnames(ClassIIAssays2)[1] <- 'AssayID'
colnames(ClassIIAssays2)[2] <- 'EC50(Rosiglitazone Maleate)'
colnames(ClassIIAssays2)[3] <- 'EC50(Troglidazone)'
colnames(ClassIIAssays2)[4] <- 'Assay Component Name'
colnames(ClassIIAssays2)[5] <- 'Biological Process'
colnames(ClassIIAssays2)[6] <- 'Target Family'
colnames(ClassIIAssays2)[7] <- 'Target SubFamily'
colnames(ClassIIAssays2)[8] <- 'Target Full Name'
colnames(ClassIIAssays2)[9] <- 'Target Gene Name'
ClassIIAssays2$Class <- 'Class 2' 

write.csv(ClassIIAssays2, "ClassIIAssays.csv", row.names=FALSE)
  
##rbind ClassIassays2 and ClassIIassays2 for data exploration 
AllClassesPosAssays<- as.data.frame(rbind(ClassIAssays2, ClassIIAssays2))
write.csv(AllClassesPosAssays, "PosAssaysBothClasses.csv", row.names = FALSE)


#################################################################################################################################################################################################################
#######################################################################################3. Biological processes per each drug######################################################################################
#################################################################################################################################################################################################################

  # Generate stacked bar plot and table of all biological process based on Rosiglitazone alone or Troglitazone alone or both drugs

BiolProcByDrug<-read.csv("BiolProcByDrug.csv")
colnames(BiolProcByDrug)[1] <- "Class"
colnames(BiolProcByDrug)[1] <- "cell cycle or cell morphology"
colnames(BiolProcByDrug)[1] <- "regulator of transcription factor activity"
colnames(BiolProcByDrug)[1] <- "regulator of gene expression"


LongBiolProcByDrug <- gather(BiolProcByDrug, BiologicalProcess, Count, 2:4)  # convert data to long format for ggplot - stacked bar plot
colnames(LongBiolProcByDrug)[1] <- 'Class'

ggplot(data=LongBiolProcByDrug, aes(x=Class, y= Count, fill= BiologicalProcess)) + 
  geom_bar(stat="identity")  + theme_classic(base_size = 18) +
  theme(legend.position = "bottom", legend.title = element_blank()) + labs (y = "Count", x = " ") # stacked bar plot of all affected biological processes for both drugs


                           

#################################################################################################################################################################################################################
#####4.#Visualize targets affected by biological processes,  # tests per target protein, & EC50 values for common test for each drug (import excel cleaned dataframe based on ClassIAssays 2 Plus ClassIIAssays2)
#################################################################################################################################################################################################################

###################################
#####A1. Transcriptional Targets####
###################################

###Transcriptional Targets (Panoramic map of all targets activated by each drug type and grouped by their target familes)
UniqueTRNTargets<-read.csv("Panoramic_Trn_Targets.csv")

ggplot(data=UniqueTRNTargets, aes(x=Target, y = Class)) + 
  geom_tile(aes(fill = Target.Family), color = "white", size = 0.1) + coord_flip() + 
  theme_classic(base_size = 16) + scale_fill_discrete(name="Target Family") 


### Transcriptional Targets (Number of "positive tests for each targeted protein per drug) 
TranscriptionalTargets<-read.csv("Trn_Targets_By_Drug.csv")

TranscriptionTargetsClassTable<-as.data.frame(table(TranscriptionalTargets$Target, TranscriptionalTargets$Class))
colnames(TranscriptionTargetsClassTable)[1] <- 'Target'
colnames(TranscriptionTargetsClassTable)[2] <- 'Class'
colnames(TranscriptionTargetsClassTable)[3] <- 'Count'

TranscriptionTargetsClassTable <- TranscriptionTargetsClassTable[-c(1:4, 6:12, 16:18, 25:37, 40:42,44:48, 50:58, 62:67, 69:82, 84:86, 92, 99,101:111, 124:126, 129), ]

TranscriptionTargetsClassTable <- TranscriptionTargetsClassTable [order(-TranscriptionTargetsClassTable$Count), ]

ggplot(data=TranscriptionTargetsClassTable, aes(x=reorder(Target, Count), y= Count, fill= Class)) + 
  geom_bar(stat="identity") + coord_flip() + theme_classic(base_size = 16) + labs (x= "Target", y = 
 'Count of distinct positive assays per target') + theme(legend.position = "bottom", legend.title = element_blank())


#################################################################################################################################
############Activation Score stratification of Transcriptional targets for each drug#############################################
############Cmax and EC50 are expressed in micromolar units & based on Avandia/Rezulin product leaflets & pharmapendium data#####
##############Rosiglitazone Maleate Cmax = 1.26 micromolar, Troglitazone Cmax = 6.38 micromolar #################################
#################################################################################################################################

###########################
##########A2.Troglitazone - NAS based stratification of transcriptional targets
###########################
TRN<-read.csv("C1C2_TRN.csv", na.strings=c("#N/A", "NA", " "))
str(TRN)

##Troglidazone - transcriptional targets activation score stratification

TrogTRN <-TRN[!(is.na(TRN$EC50.Troglidazone)), ] # select all tests positive with Troglidazone

TrogTRN<-data.frame(TrogTRN$EC50.Troglidazone, TrogTRN$Cmax.Troglidazone, 
                    TrogTRN$Target.Assay., TrogTRN$Target) #subset selected data 
colnames(TrogTRN)[1]<- "EC50"
colnames(TrogTRN)[2]<- "Cmax(micromolar)"
colnames(TrogTRN)[3]<- "AssayTarget"
colnames(TrogTRN)[4]<- "Target"

TrogTRN <-TrogTRN[order(TrogTRN$EC50), ] # order Trog+ve assay targets based on their EC50 values

TrogTRN$AssayTarget <- factor(TrogTRN$AssayTarget, levels = TrogTRN$AssayTarget) # retain order based on EC50 for ggplot

temp <- c()
for(i in TrogTRN$EC50){
if (i < 6.387316) {             # < Cmax (at highest dose = 6.387316 micromolar) set as threshold for "high" activation score
  temp <- append(temp, "high")
  } else if (i< 31.93657984) {  # > Cmax but < 5*Cmax set as threshold for "medium" activation score 
  temp <- append(temp, "medium")
  } else {
  temp <- append(temp, "low")   # > 5* Cmax set as threshold for "low" activation score
  }
} # create a temp vector with high/medium/low expression label for appending to TrogTRN (based on Cmax versus EC50 value)

TrogTRN$Activn <- temp # append temp vector with labels for ggplot purposes

TrogTRN$ActivnScore <- (TrogTRN$Cmax-TrogTRN$EC50)/(TrogTRN$Cmax) # create normalized activation score based on distance of EC50 from Cmax; more positive indicates higher activation potential

TrogTRN$ActivnScore <- round((TrogTRN$ActivnScore), 2) # round Activation score to 2 digit

## plot all targets/assay combinations against their activation scores (includes duplicated targets that may be positive in multiple assays)

library("RColorBrewer")
ggplot(TrogTRN, aes(x=TrogTRN$AssayTarget, y=TrogTRN$ActivnScore)) + 
  geom_point(stat='identity', aes(col=TrogTRN$Activn), size = 8, alpha = 0.8) + coord_flip() + 
  scale_color_manual(values = c(brewer.pal(n = 3, name = "Accent"))) +
  geom_text(label = TrogTRN$ActivnScore, color="black", size=3)+ 
  geom_segment(aes(y = 0, x = TrogTRN$AssayTarget, yend = TrogTRN$ActivnScore, 
  xend = TrogTRN$AssayTarget), color = "honeydew4") + theme_classic()+ theme(axis.text=element_text(size=12)) +
  labs(y= "*(Normalized) Activation Score", x = "Target(Assay#)")

write.csv(TrogTRN, "TrogTRN.csv", row.names=FALSE) # write to excel/csv; remove duplicate targets from Trog TRN- retaining entry for each target with highest activation score

TrogTRN2 <-read.csv("TrogTRN_nodup_targets.csv")

TrogTRN2 <-TrogTRN2[order(TrogTRN2$EC50), ] # order Trog+ve assay targets based on their EC50 values

TrogTRN2$TargetAssay <- factor(TrogTRN2$AssayTarget, levels = TrogTRN2$AssayTarget) # retain order based on EC50 for ggplot

## plot all unique targets against their activation scores (highest activation score retained for each target if duplicated across multiple assays)

ggplot(TrogTRN2, aes(x=TrogTRN2$TargetAssay, y=TrogTRN2$ActivnScore)) + 
  geom_point(stat='identity', aes(col=TrogTRN2$Activn), size = 8, alpha = 0.8) + coord_flip() + 
  scale_color_manual(values = c("tomato", "steelblue", "orange")) +
  geom_text(label = TrogTRN2$ActivnScore, color="black", size=3)+ 
  geom_segment(aes(y = 0, x = TrogTRN2$TargetAssay, yend = TrogTRN2$ActivnScore, 
                   xend = TrogTRN2$TargetAssay), color = "honeydew4") + theme_classic()+ 
  theme(axis.text=element_text(size=12)) + labs(y= "*(Normalized) Activation Score", x = "Target(Assay#)")+
  theme(legend.position = "none")


#######################################
##########A3. Rosiglitazone Maleate - NAS based transcriptional targets stratification
######################################
RosMalTRN <-TRN[!is.na(TRN$EC50.Rosiglitazone.Maleate), ] # select all tests positive with Rosiglitazone Maleate

RosMalTRN<-data.frame(RosMalTRN$EC50.Rosiglitazone.Maleate, RosMalTRN$Cmax.Rosiglitazone.Maleate, 
                      RosMalTRN$Target.Assay., RosMalTRN$Target) #subset selected data 
colnames(RosMalTRN)[1]<- "EC50"
colnames(RosMalTRN)[2]<- "Cmax(micromolar)"
colnames(RosMalTRN)[3]<- "Target(Assay)"
colnames(RosMalTRN)[4]<- "Target"

RosMalTRN <-RosMalTRN[order(RosMalTRN$EC50), ] # order Rosiglitazone +ve assay targets based on their EC50 values

RosMalTRN$`Target(Assay)` <- factor(RosMalTRN$`Target(Assay)`, levels = RosMalTRN$`Target(Assay)`) # retain order based on EC50 for ggplot

temp <- c()
for(i in RosMalTRN$EC50){
  if (i < 1.262882) {             # < Cmax at highest dose set as threshold for "high" activation score
    temp <- append(temp, "high")
  } else if (i< 6.31441) {  # > Cmax but < 5*Cmax set as threshold for "medium" activation score 
    temp <- append(temp, "medium")
  } else {
    temp <- append(temp, "low")   # > 5* Cmax set as threshold for "low" activation score
  }
} # create a temp vector with high/medium/low expression label for appending to RosMalTRN (based on Cmax/EC50 based "Activation Score")

RosMalTRN$Activn <- temp # append temp vector with labels for ggplot purposes

RosMalTRN$ActivnScore <- (RosMalTRN$Cmax-RosMalTRN$EC50)/(RosMalTRN$Cmax) # create normalized activation score based on distance of EC50 from Cmax; more positive indicates higher activation potential

RosMalTRN$ActivnScore <- round((RosMalTRN$ActivnScore), 2) # round Activation score to 2 digit

## plot all targets/assay combinations against their activation scores (includes duplicated targets that may be positive in multiple assays)

library("RColorBrewer")
ggplot(RosMalTRN, aes(x=RosMalTRN$`Target(Assay)`, y=RosMalTRN$ActivnScore)) + 
  geom_point(stat='identity', aes(col=RosMalTRN$Activn), size = 8, alpha = 0.8) + coord_flip() + 
  scale_color_manual(values = c(brewer.pal(n = 3, name = "Accent"))) +
  geom_text(label = RosMalTRN$ActivnScore, color="black", size=3)+ 
  geom_segment(aes(y = 0, x = RosMalTRN$`Target(Assay)`, yend = RosMalTRN$ActivnScore, 
                   xend = RosMalTRN$`Target(Assay)`), color = "honeydew4") + theme_classic()+ theme(axis.text=element_text(size=12)) +
  labs(y= "*(Normalized) Activation Score", x = "Target(Assay#)")

write.csv(RosMalTRN, "RosMalTRN.csv", row.names=FALSE) # write to excel/csv; remove duplicate targets from RosMal TRN- retaining entry for each target with highest activation score

RosMalTRN2 <-read.csv("RosMalTRN_nodup_targets.csv")

RosMalTRN2 <-RosMalTRN2[order(RosMalTRN2$EC50), ] # order Ros+ve assay targets based on their EC50 values

colnames(RosMalTRN2)[2]<-"Cmax"
colnames(RosMalTRN2)[3]<-"TargetAssay"

RosMalTRN2$TargetAssay <- factor(RosMalTRN2$TargetAssay, levels = RosMalTRN2$TargetAssay) # retain order based on EC50 for ggplot

## plot all unique targets against their activation scores (highest activation score retained for each target if duplicated across multiple assays)

ggplot(RosMalTRN2, aes(x=RosMalTRN2$TargetAssay, y=RosMalTRN2$ActivnScore)) + 
  geom_point(stat='identity', aes(col=RosMalTRN2$Activn), size = 8, alpha = 0.8) + coord_flip() + 
  scale_color_manual(values = c("tomato", "steelblue", "orange")) +
  geom_text(label = RosMalTRN2$ActivnScore, color="black", size=3)+ 
  geom_segment(aes(y = 0, x = RosMalTRN2$TargetAssay, yend = RosMalTRN2$ActivnScore, 
                   xend = RosMalTRN2$TargetAssay), color = "honeydew4") + theme_classic()+ 
  theme(axis.text=element_text(size=12)) + labs(y= "*(Normalized) Activation Score", x = "Target(Assay#)")+
  theme(legend.position = "none")


##################################################
#####B1. Gene expression regulatory Targets#######
##################################################


###Gene expression regulation targets (Panoramic Map of all targets activated by each drug type and grouped by their target familes)
UniqueGXPTargets<-read.csv("Panoramic_Gxp_Targets.csv")


ggplot(data=UniqueGXPTargets, aes(x=Target, y = Class)) + 
  geom_tile(aes(fill = Target.Family), color = "white", size = 0.1) + coord_flip() + 
  theme_classic(base_size = 16) + scale_fill_discrete(name="Target Family") 

### Gene Expresion Targets (Number of "positive tests for each targeted protein per drug) 
GeneExpRegulnTargets<-read.csv("Gxp_Targets_By_Drug.csv")

GeneExpTargetsClassTable<-as.data.frame(table(GeneExpRegulnTargets$Target, GeneExpRegulnTargets$Class))
colnames(GeneExpTargetsClassTable)[1] <- 'Target'
colnames(GeneExpTargetsClassTable)[2] <- 'Class'
colnames(GeneExpTargetsClassTable)[3] <- 'Count'

GeneExpTargetsClassTable <- GeneExpTargetsClassTable[which(GeneExpTargetsClassTable$Count !=0), ]

GeneExpTargetsClassTable <- GeneExpTargetsClassTable [order(-GeneExpTargetsClassTable$Count), ]

ggplot(data=GeneExpTargetsClassTable, aes(x=reorder(Target, Count), y= Count, fill= Class)) + 
  geom_bar(stat="identity") + coord_flip() + theme_classic(base_size = 16) + labs (x= "Target", y = 
  'Count of distinct positive assays per target') + theme(legend.position = "bottom", legend.title = element_blank()) 


#################################################################################################################################
############Activation Score stratification of gene expression regulatory targetss for each drug#################################
############Cmax and EC50 are expressed in micromolar units & based on Avandia/Rezulin product leaflets & pharmapendium data#####
##############Rosiglitazone Maleate Cmax = 1.26 micromolar, Troglitazone Cmax = 6.38 micromolar #################################
#################################################################################################################################


GXP<-read.csv("C1C2_GXP.csv", na.strings=c("#N/A", "NA", " "))

###########################
##########B2.Troglitazone - NAS based stratification of gene expression regulatory targets
###########################


TrogGXP <-GXP[!is.na(GXP$EC50.Troglidazone), ] # select all tests positive with Troglitazone

TrogGXP<-data.frame(TrogGXP$EC50.Troglidazone, TrogGXP$Cmax.Troglidazone, 
                    TrogGXP$Target.Assay., TrogGXP$Target) #subset selected data 
colnames(TrogGXP)[1]<- "EC50"
colnames(TrogGXP)[2]<- "Cmax(micromolar)"
colnames(TrogGXP)[3]<- "Target(Assay)"
colnames(TrogGXP)[4]<- "Target"

TrogGXP <-TrogGXP[order(TrogGXP$EC50), ] # order Trog+ve assay targets based on their EC50 values

TrogGXP$`Target(Assay)` <- factor(TrogGXP$`Target(Assay)`, levels = TrogGXP$`Target(Assay)`) # retain order based on EC50 for ggplot

temp <- c()
for(i in TrogGXP$EC50){
  if (i < 6.387316) {             # < Cmax (at highest dose = 6.387316 micromolar) set as threshold for "high" activation score
    temp <- append(temp, "high")
  } else if (i< 31.93657984) {  # > Cmax but < 5*Cmax set as threshold for "medium" activation score 
    temp <- append(temp, "medium")
  } else {
    temp <- append(temp, "low")   # > 5* Cmax set as threshold for "low" activation score
  }
} # create a temp vector with high/medium/low expression label for appending to TrogTRN (based on Cmax versus EC50 value)

TrogGXP$Activn <- temp # append temp vector with labels for ggplot purposes

TrogGXP$ActivnScore <- (TrogGXP$Cmax-TrogGXP$EC50)/(TrogGXP$Cmax) # create normalized activation score based on distance of EC50 from Cmax; more positive indicates higher activation potential

TrogGXP$ActivnScore <- round((TrogGXP$ActivnScore), 2) # round Activation score to 2 digit

## plot all targets/assay combinations against their activation scores (includes duplicated targets that may be positive in multiple assays)

library("RColorBrewer")
ggplot(TrogGXP, aes(x=TrogGXP$`Target(Assay)`, y=TrogGXP$ActivnScore)) + 
  geom_point(stat='identity', aes(col=TrogGXP$Activn), size = 8, alpha = 0.8) + coord_flip() + 
  scale_color_manual(values = c(brewer.pal(n = 3, name = "Accent"))) +
  geom_text(label = TrogGXP$ActivnScore, color="black", size=3)+ 
  geom_segment(aes(y = 0, x = TrogGXP$`Target(Assay)`, yend = TrogGXP$ActivnScore, 
                   xend = TrogGXP$`Target(Assay)`), color = "honeydew4") + theme_classic()+ 
  theme(axis.text=element_text(size=12)) +labs(y= "*(Normalized) Activation Score", x = "Target(Assay#)")

write.csv(TrogGXP, "TrogGXP.csv", row.names = FALSE) # write to excel/csv; remove duplicate targets from TrogGXP- retaining entry for each target with highest activation score

TrogGXP2 <-read.csv("TrogGXP_nodup_targets.csv") # read file from excel/csv with each target represented only once with its highest activation score (lower activation scores removed)

TrogGXP2 <-TrogGXP2[order(TrogGXP2$EC50), ] # order Trog+ve assay targets based on their EC50 values

TrogGXP2$Target.Assay. <- factor(TrogGXP2$Target.Assay., levels = TrogGXP2$Target.Assay.) # retain order based on EC50 for ggplot

## plot all unique targets against their activation scores (highest activation score retained for each target if duplicated across multiple assays)

ggplot(TrogGXP2, aes(x=TrogGXP2$Target.Assay., y=TrogGXP2$ActivnScore)) + 
  geom_point(stat='identity', aes(col=TrogGXP2$Activn), size = 8, alpha = 0.8) + coord_flip() + 
  scale_color_manual(values = c("tomato", "steelblue", "orange")) +
  geom_text(label = TrogGXP2$ActivnScore, color="black", size=3)+ 
  geom_segment(aes(y = 0, x = TrogGXP2$Target.Assay., yend = TrogGXP2$ActivnScore, 
                   xend = TrogGXP2$Target.Assay.), color = "honeydew4") + theme_classic()+ 
  theme(axis.text=element_text(size=12)) + labs(y= "*(Normalized) Activation Score", x = "Target(Assay#)")+
  theme(legend.position = "none")

###################################
##########B3.Rosiglitazone Maleate - NAS based stratification of gene expression regulatory targets
###################################

RosMalGXP <-GXP[!is.na(GXP$EC50.Rosiglitazone.Maleate), ] # select all tests positive with Rosiglitazone Maleate

RosMalGXP<-data.frame(RosMalGXP$EC50.Rosiglitazone.Maleate, RosMalGXP$Cmax.Rosiglitazone.Maleate, 
                      RosMalGXP$Target.Assay., RosMalGXP$Target) #subset selected data 
colnames(RosMalGXP)[1]<- "EC50"
colnames(RosMalGXP)[2]<- "Cmax(micromolar)"
colnames(RosMalGXP)[3]<- "Target(Assay)"
colnames(RosMalGXP)[4]<- "Target"

RosMalGXP <-RosMalGXP[order(RosMalGXP$EC50), ] # order Rosiglitazone +ve assay targets based on their EC50 values

RosMalGXP$`Target(Assay)` <- factor(RosMalGXP$`Target(Assay)`, levels = RosMalGXP$`Target(Assay)`) # retain order based on EC50 for ggplot

temp <- c()
for(i in RosMalGXP$EC50){
  if (i < 1.262882) {             # < Cmax at highest dose set as threshold for "high" activation score
    temp <- append(temp, "high")
  } else if (i< 6.31441) {  # > Cmax but < 5*Cmax set as threshold for "medium" activation score 
    temp <- append(temp, "medium")
  } else {
    temp <- append(temp, "low")   # > 5* Cmax set as threshold for "low" activation score
  }
} # create a temp vector with high/medium/low expression label for appending to RosMalTRN (based on Cmax/EC50 based "Activation Score")

RosMalGXP$Activn <- temp # append temp vector with labels for ggplot purposes

RosMalGXP$ActivnScore <- (RosMalGXP$Cmax-RosMalGXP$EC50)/(RosMalGXP$Cmax) # create normalized activation score based on distance of EC50 from Cmax; more positive indicates higher activation potential

RosMalGXP$ActivnScore <- round((RosMalGXP$ActivnScore), 2) # round Activation score to 2 digit

## plot all targets/assay combinations against their activation scores (includes duplicated targets that may be positive in multiple assays)

library("RColorBrewer")
ggplot(RosMalGXP, aes(x=RosMalGXP$`Target(Assay)`, y=RosMalGXP$ActivnScore)) + 
  geom_point(stat='identity', aes(col=RosMalGXP$Activn), size = 8, alpha = 0.8) + coord_flip() + 
  scale_color_manual(values = c(brewer.pal(n = 3, name = "Accent"))) +
  geom_text(label = RosMalGXP$ActivnScore, color="black", size=3)+ 
  geom_segment(aes(y = 0, x = RosMalGXP$`Target(Assay)`, yend = RosMalGXP$ActivnScore, 
                   xend = RosMalGXP$`Target(Assay)`), color = "honeydew4") + theme_classic()+ theme(axis.text=element_text(size=12)) +
  labs(y= "*(Normalized) Activation Score", x = "Target(Assay#)")

write.csv(RosMalGXP, "RosMalGXP.csv", row.names = FALSE) # write to excel/csv; remove duplicate targets from RosMal GXP- retaining entry for each target with highest activation score

RosMalGXP2 <-read.csv("RosMalGXP_nodup_targets.csv") # read in csv file with duplicate targets removed

colnames(RosMalGXP2)[2]<-"Cmax"
colnames(RosMalGXP2)[3]<-"TargetAssay"

RosMalGXP2 <-RosMalGXP2[order(RosMalGXP2$EC50), ] # order Rosaglitazone Maleate +ve assay targets based on their EC50 values

RosMalGXP2$TargetAssay <- factor(RosMalGXP2$TargetAssay, levels = RosMalGXP2$TargetAssay) # retain order based on EC50 for ggplot

## plot all unique targets against their activation scores (highest activation score retained for each target if duplicated across multiple assays)

ggplot(RosMalGXP2, aes(x=RosMalGXP2$TargetAssay, y=RosMalGXP2$ActivnScore)) + 
  geom_point(stat='identity', aes(col=RosMalGXP2$Activn), size = 8, alpha = 0.8) + coord_flip() + 
  scale_color_manual(values = c("tomato", "steelblue", "orange")) +
  geom_text(label = RosMalGXP2$ActivnScore, color="black", size=3)+ 
  geom_segment(aes(y = 0, x = RosMalGXP2$TargetAssay, yend = RosMalGXP2$ActivnScore, 
                   xend = RosMalGXP2$TargetAssay), color = "honeydew4") + theme_classic()+ 
  theme(axis.text=element_text(size=12)) + labs(y= "*(Normalized) Activation Score", x = "Target(Assay#)")+
  theme(legend.position = "none")

##################################################
#####C1. Cell Cycle/Cell Morphology Targets#######
##################################################

###Cell Cycle/Morphology targets (Panoramic Map of all targets activated by each drug type and grouped by their target familes)
UniqueCellCycleMorphTargets<-read.csv("Panoramic_CCM_Targets.csv")


ggplot(data=UniqueCellCycleMorphTargets, aes(x=Target, y = Class)) + 
  geom_tile(aes(fill = Target.Family), color = "white", size = 0.1) + coord_flip() + 
  theme_classic(base_size = 16) + scale_fill_discrete(name="Target Family") 

### Cell Cycle/Morphology targets (Number of "positive tests for each targeted protein per drug) 
CellCyclMorphTargets<-read.csv("CCM_Targets_By_Drug.csv")

CellCycleMorphClassTable<-as.data.frame(table(CellCyclMorphTargets$Target, CellCyclMorphTargets$Class))
colnames(CellCycleMorphClassTable)[1] <- 'Target'
colnames(CellCycleMorphClassTable)[2] <- 'Class'
colnames(CellCycleMorphClassTable)[3] <- 'Count'


CellCycleMorphClassTable <- CellCycleMorphClassTable[which(CellCycleMorphClassTable$Count !=0), ]

CellCycleMorphClassTable<- CellCycleMorphClassTable [order(-CellCycleMorphClassTable$Count), ]

ggplot(data=CellCycleMorphClassTable, aes(x=reorder(Target, Count), y= Count, fill= Class)) + 
  geom_bar(stat="identity") + coord_flip() + theme_classic(base_size = 16) + labs (x= "Target", y = 
                                                                                     'Count of distinct positive assays per target') +
  scale_y_continuous(breaks = seq(0,1)) + theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_fill_manual(values=c("steel blue", "pink"))



###########################
##########C2.Troglitazone - NAS based stratification of cell cycle/cell morphology targets
###########################

CCM<-read.csv("C1C2_CCM.csv", na.strings=c("#N/A", "NA", " "))


TrogCCM <-CCM[!is.na(CCM$EC50.Troglidazone), ] # select all tests positive with Troglitazone

TrogCCM<-data.frame(TrogCCM$EC50.Troglidazone, TrogCCM$Cmax.Troglidazone, 
                    TrogCCM$Target.Assay., TrogCCM$Target) #subset selected data 
colnames(TrogCCM)[1]<- "EC50"
colnames(TrogCCM)[2]<- "Cmax(micromolar)"
colnames(TrogCCM)[3]<- "Target(Assay)"
colnames(TrogCCM)[4]<- "Target"

TrogCCM <-TrogCCM[order(TrogCCM$EC50), ] # order Trog+ve assay targets based on their EC50 values

TrogCCM$`Target(Assay)` <- factor(TrogCCM$`Target(Assay)`, levels = TrogCCM$`Target(Assay)`) # retain order based on EC50 for ggplot

temp <- c()
for(i in TrogCCM$EC50){
  if (i < 6.387316) {             # < Cmax (at highest dose = 6.387316 micromolar) set as threshold for "high" activation score
    temp <- append(temp, "high")
  } else if (i< 31.93657984) {  # > Cmax but < 5*Cmax set as threshold for "medium" activation score 
    temp <- append(temp, "medium")
  } else {
    temp <- append(temp, "low")   # > 5* Cmax set as threshold for "low" activation score
  }
} # create a temp vector with high/medium/low expression label for appending to TrogCCM (based on Cmax versus EC50 value)

TrogCCM$Activn <- temp # append temp vector with labels for ggplot purposes

TrogCCM$ActivnScore <- (TrogCCM$Cmax-TrogCCM$EC50)/(TrogCCM$Cmax) # create normalized activation score based on distance of EC50 from Cmax; more positive indicates higher activation potential

TrogCCM$ActivnScore <- round((TrogCCM$ActivnScore), 2) # round Activation score to 2 digit

## plot all targets/assay combinations against their activation scores (includes duplicated targets that may be positive in multiple assays)

library("RColorBrewer")
ggplot(TrogCCM, aes(x=TrogCCM$`Target(Assay)`, y=TrogCCM$ActivnScore)) + 
  geom_point(stat='identity', aes(col=TrogCCM$Activn), size = 8, alpha = 0.8) + coord_flip() + 
  scale_color_manual(values = c(brewer.pal(n = 3, name = "Accent"))) +
  geom_text(label = TrogCCM$ActivnScore, color="black", size=3)+ 
  geom_segment(aes(y = 0, x = TrogCCM$`Target(Assay)`, yend = TrogCCM$ActivnScore, 
                   xend = TrogCCM$`Target(Assay)`), color = "honeydew4") + theme_classic()+ 
  theme(axis.text=element_text(size=12)) +labs(y= "*(Normalized) Activation Score", x = "Target(Assay#)")

write.csv(TrogCCM, "TrogCCM.csv", row.names = FALSE) # write to excel/csv; remove duplicate targets from Trog CCM- retaining entry for each target with highest activation score

TrogCCM2 <-read.csv("TrogCCM_nodup_targets.csv") # read file from excel/csv with each target represented only once with its highest activation score (lower activation scores removed)

TrogCCM2 <-TrogCCM2[order(TrogCCM2$EC50), ] # order Trog+ve assay targets based on their EC50 values
colnames(TrogCCM2)[2]<-"Cmax"
colnames(TrogCCM2)[3]<-"TargetAssay"

TrogCCM2$TargetAssay <- factor(TrogCCM2$TargetAssay, levels = TrogCCM2$TargetAssay) # retain order based on EC50 for ggplot

## plot all unique targets against their activation scores (highest activation score retained for each target if duplicated across multiple assays)

ggplot(TrogCCM2, aes(x=TrogCCM2$TargetAssay, y=TrogCCM2$ActivnScore)) + 
  geom_point(stat='identity', aes(col=TrogCCM2$Activn), size = 8, alpha = 0.8) + coord_flip() + 
  scale_color_manual(values = c("tomato", "steelblue", "orange")) +
  geom_text(label = TrogCCM2$ActivnScore, color="black", size=4)+ 
  geom_segment(aes(y = 0, x = TrogCCM2$TargetAssay, yend = TrogCCM2$ActivnScore, 
                   xend = TrogCCM2$TargetAssay), color = "honeydew4") + theme_classic()+ 
  theme(axis.text=element_text(size=12)) + labs(y= "*(Normalized) Activation Score", x = "Target(Assay#)")+
  theme(legend.position = "none")


###########################
##########C3.Rosiglitazone Maleate - NAS based stratification of cell cycle/cell morphology targets
###########################

RosMalCCM <-CCM[!is.na(CCM$EC50.Rosiglitazone.Maleate), ] # select all tests positive with Rosiglitazone Maleate

RosMalCCM<-data.frame(RosMalCCM$EC50.Rosiglitazone.Maleate, RosMalCCM$Cmax.Rosiglitazone.Maleate, 
                      RosMalCCM$Target.Assay., RosMalCCM$Target) #subset selected data 
colnames(RosMalCCM)[1]<- "EC50"
colnames(RosMalCCM)[2]<- "Cmax"
colnames(RosMalCCM)[3]<- "Target(Assay)"
colnames(RosMalCCM)[4]<- "Target"

RosMalCCM <-RosMalCCM[order(RosMalCCM$EC50), ] # order Rosiglitazone +ve assay targets based on their EC50 values

RosMalCCM$`Target(Assay)` <- factor(RosMalCCM$`Target(Assay)`, levels = RosMalCCM$`Target(Assay)`) # retain order based on EC50 for ggplot

temp <- c()
for(i in RosMalCCM$EC50){
  if (i < 1.262882) {             # < Cmax at highest dose set as threshold for "high" activation score
    temp <- append(temp, "high")
  } else if (i< 6.31441) {  # > Cmax but < 5*Cmax set as threshold for "medium" activation score 
    temp <- append(temp, "medium")
  } else {
    temp <- append(temp, "low")   # > 5* Cmax set as threshold for "low" activation score
  }
} # create a temp vector with high/medium/low expression label for appending to RosMalTRN (based on Cmax/EC50 based "Activation Score")

RosMalCCM$Activn <- temp # append temp vector with labels for ggplot purposes

RosMalCCM$ActivnScore <- (RosMalCCM$Cmax-RosMalCCM$EC50)/(RosMalCCM$Cmax) # create normalized activation score based on distance of EC50 from Cmax; more positive indicates higher activation potential

RosMalCCM$ActivnScore <- round((RosMalCCM$ActivnScore), 2) # round Activation score to 2 digit

## plot all targets/assay combinations against their activation scores (includes duplicated targets that may be positive in multiple assays)

library("RColorBrewer")
ggplot(RosMalCCM, aes(x=RosMalCCM$`Target(Assay)`, y=RosMalCCM$ActivnScore)) + 
  geom_point(stat='identity', aes(col=RosMalCCM$Activn), size = 8, alpha = 0.8) + coord_flip() + 
  scale_color_manual(values = c(brewer.pal(n = 3, name = "Accent"))) +
  geom_text(label = RosMalCCM$ActivnScore, color="black", size=3)+ 
  geom_segment(aes(y = 0, x = RosMalCCM$`Target(Assay)`, yend = RosMalCCM$ActivnScore, 
                   xend = RosMalCCM$`Target(Assay)`), color = "honeydew4") + theme_classic()+ 
  theme(axis.text=element_text(size=12)) + labs(y= "*(Normalized) Activation Score", x = "Target(Assay#)") +
  theme(legend.position = "none")

write.csv(RosMalCCM, "RosMalCCM.csv", row.names=FALSE) # write to excel; remove duplicate targets from RosMal CCM- retaining entry for each target with highest activation score

RosMalCCM2 <-read.csv("RosMalCCM_nodup_targets.csv")

RosMalCCM2 <-RosMalCCM2[order(RosMalCCM2$EC50), ] # order Rosaglitazone Maleate +ve assay targets based on their EC50 values
colnames(RosMalCCM2)[2]<-"Cmax"
colnames(RosMalCCM2)[3]<-"TargetAssay"


RosMalCCM2$TargetAssay <- factor(RosMalCCM2$TargetAssay, levels = RosMalCCM2$TargetAssay) # retain order based on EC50 for ggplot

## plot all unique targets against their activation scores (highest activation score retained for each target if duplicated across multiple assays)

ggplot(RosMalCCM2, aes(x=RosMalCCM2$TargetAssay, y=RosMalCCM2$ActivnScore)) + 
  geom_point(stat='identity', aes(col=RosMalCCM2$Activn), size = 8, alpha = 0.8) + coord_flip() + 
  scale_color_manual(values = c("steelblue")) +
  geom_text(label = RosMalCCM2$ActivnScore, color="black", size=4)+ 
  geom_segment(aes(y = 0, x = RosMalCCM2$TargetAssay, yend = RosMalCCM2$ActivnScore, 
                   xend = RosMalCCM2$TargetAssay), color = "honeydew4") + theme_classic()+ 
  theme(axis.text=element_text(size=12)) + labs(y= "*(Normalized) Activation Score", x = "Target(Assay#)")+
  theme(legend.position = "none")


#################################################################################################################################
#################################################################################################################################
#####5. Prepare data to generate a heat map of all tests/targets affected by each drug with their NAS scores#######################
######This will be done in following steps:
#####A. For both drugs, generate single table with NAS values for all targets with TRN category; remove extraneous columns 
#####B. For both drugs, generate single table with NAS values for all targets with GXP category; remove extraneous columns
#####C. For both drugs, generate single table with NAS values for all targets with CCM category; remove extraneous columns
#####D. rbind all tables and proceed to heat mapping
#########impute any missing values created by merges to -1000000 (indicates test/target is not activated)########################
#################################################################################################################################
#################################################################################################################################

#A.Generate common table with NAS values for all targets within Transcription (TRN) regulators category
RosMalTRN2$RosMalNAS <-RosMalTRN2$ActivnScore
TrogTRN2$TrogNAS<-TrogTRN2$ActivnScore
TRNTargets<-merge(TrogTRN2, RosMalTRN2, by="Target", all = TRUE)
colnames(TRNTargets)

TRNTargets$AssayTarget<-as.character(TRNTargets$AssayTarget)
TRNTargets[38:43, 4]<-c("forkhead box protein (1425)", "nuclear receptor subfamily 1, group I, member 3 (101)",
                        "peroxisome proliferator-activated receptor alpha (132)",
                        "peroxisome proliferator-activated receptor delta (1124)",
                        "TOX21_ERa_BLA_Antagonist_ch1 (1189)","TOX21_VDR_BLA_Agonist_ch2 (1130)")

TRNTargets$AssayTarget<-factor(TRNTargets$AssayTarget, levels = TRNTargets$AssayTarget) 
str(TRNTargets)

TRNTargets<-TRNTargets[, c(4,8,14)]

#B.Generate common table with NAS values for all targets within gene expression (GXP) regulators category
RosMalGXP2$RosMalNAS <-RosMalGXP2$ActivnScore
TrogGXP2$TrogNAS<-TrogGXP2$ActivnScore
GXPTargets<-merge(TrogGXP2, RosMalGXP2, by="Target", all = TRUE)
colnames(GXPTargets)
    #fill in assay/target column for any entries missed
GXPTargets$Target.Assay.<-as.character(GXPTargets$Target.Assay.)
GXPTargets[32:33, 4]<- c("low density lipoprotein receptor (214)", "prostaglandin E receptor 2 (290)")
GXPTargets$Target.Assay.<-factor(GXPTargets$Target.Assay., levels=GXPTargets$Target.Assay.)

colnames(GXPTargets)[4]<-"AssayTarget"

GXPTargets<-GXPTargets[, c(4,7,13)]

#C.Generate common table with NAS values for all targets within cell cycle & morphology (CCM) targets category
RosMalCCM2$RosMalNAS <-RosMalCCM2$ActivnScore
TrogCCM2$TrogNAS<-TrogCCM2$ActivnScore
CCMTargets<-merge(TrogCCM2, RosMalCCM2, by="Target", all = TRUE)

colnames(CCMTargets)[4]<-"AssayTarget"

CCMTargets<-CCMTargets[, c(4,7,13)]

#D. rbind all tests/targets for both drugs and impute missing values created by above merges
HMDAll<-rbind(TRNTargets, GXPTargets, CCMTargets)
HMDAll[is.na(HMDAll)]<--1000000

data.frame(table(HMDAll$AssayTarget))


#########################################################################################################
###################6. Heatmap all targets across all biological processes################################
#########################################################################################################

set.seed(10)

row.names(HMDAll) <- HMDAll$Target # Converts target names to "rows" of heatmap
AllTargData <- as.matrix(HMDAll[ , c(2:3)]) # extract heat map "values" i.e., activation scores into a matrix format

fontsize <- 0.2
require(circlize)
coul = colorRamp2(c(-1000000, -40, -20, -10, -5, -2, 0, 1), c("gray14", "purple", "blue", "steelblue", "salmon", "orange", "tomato", "firebrick"))
Heatmap(AllTargData, cluster_columns = TRUE,
        row_names_side = "left",
        row_dend_side = "right",
        col=coul,
        row_names_gp = gpar(cex=fontsize),
        row_dend_width = unit(1, "cm"),
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D",
        km= 3)

##Simplify heat map from above by comparison of gene/pathway symbols and their NAS values across drugs
##merge gene or pathway symbols (into a column called as "GPS") if gene targets not defined for some tests ias an additional layer (to simplify "Target" column) and faciliate removal of duplicate gene names
##Remove gene names if duplicated by retaining those with highest NAS score
##regenerate heat map

write.csv(HMDAll, "HMDAll.csv", row.names = FALSE) # export for curation of gene/pathway symbols and removal of duplciated gene/pathway targets

HMDAll_GT<-read.csv("HMDAll_GeneTgts.csv") #in 2 cases (PTGER2, SERPINE1 all entries were removed as it was not possible to choose highest NAS of each replicate across both drugs - they moved orthogonally in two tests across both drugs making this an ambiguous selection)
#require(circlize)

set.seed(104)

row.names(HMDAll_GT) <- HMDAll_GT$GPS # Converts gene names to "rows" of heatmap
GeneTargets <- as.matrix(HMDAll_GT[ , c(2:3)]) # extract heat map "values" i.e., activation scores into a matrix format

fontsize <- 0.77

coul = colorRamp2(c(-1000000, -40, -20, -10, -5, -2, 0, 1), c("gray14", "purple", "blue", "steelblue", "salmon", "orange", "tomato", "firebrick"))
Heatmap(GeneTargets, cluster_columns = TRUE,
        row_names_side = "left",
        row_dend_side = "right",
        col=coul,
        row_names_gp = gpar(cex=fontsize),
        row_dend_width = unit(1, "cm"),
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D",
        km= 3)


############################################################################################################################
#################################################THE END####################################################################
############################################################################################################################


