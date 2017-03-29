#########################################################################################################################
#########################################################################################################################
#
# TRANSMISION RATIO ANALYSIS
#
#########################################################################################################################
#########################################################################################################################
#                                                                                                                       #
# README:                                                                                                               #
# This script  have some comments with explanations about the code.                                                     #
# All the orders that are after a #, are there just to see how the output and input files are. For running these        #
# orders in R, you will need to erase the # sign present at the beginning of the line.                                  #
#                                                                                                                       #
# This script was generated in March 2017 by Mar Muniz, PhD student IGBMC : Herault Lab.                                #
# Last modification of the Script: 080317. Mantained by Herault Lab.                                                    #
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#############################################  PACKAGES NEEDED  #########################################################


source("http://bioconductor.org/biocLite.R")
biocLite("tidyr")
biocLite("dplyr")
biocLite("lubridate")
biocLite("data.table")

library("tidyr")
library("dplyr")
library("lubridate")
library("data.table")
##############################################################################################################

# wd: change the path accordinly to where you wanna work with these files.
wd= "/Users/marmmoreno/Desktop/YH_Lab/transmision_ratio/"  #mac
wd= "/home/munizmom/Documents/YH_Lab/transmision_ratio/" #ubuntu
setwd(wd)
getwd()

#All the files generated should be keep on output_folder
folder <- "R_Analysis/"
folder_input <- "input/" #folder name where we are gonna put the output files of this script 
folder_to_analyse <- "060317_mice_line_kinasDead/"
output_folder <- "results/"
output_data <- "data/"
output_plots <- "plots/"
output_mice_line <- "per_MiceLine/"
output_dyrk1b <- "dyrk1b/"
output_kd <- "kd/"
output_dyrk1ae2 <- "dyrk1ae2/"

output_all <- "total/"

dir.create(file.path(getwd (), folder), showWarnings = FALSE)
dir.create(file.path(getwd (), folder,), showWarnings = FALSE)
dir.create(file.path(getwd (), folder, output_folder), showWarnings = FALSE)
dir.create(file.path(getwd (), folder, output_folder,output_data), showWarnings = FALSE)
dir.create(file.path(getwd (), folder, output_folder,output_data,output_mice_line), showWarnings = FALSE)
dir.create(file.path(getwd (), folder, output_folder,output_data,output_mice_line,output_dyrk1b,folder_to_analyse), showWarnings = FALSE)
dir.create(file.path(getwd (), folder, output_folder,output_data,output_mice_line, output_kd,folder_to_analyse), showWarnings = FALSE)
dir.create(file.path(getwd (), folder, output_folder,output_data,output_mice_line, output_dyrk1ae2,folder_to_analyse), showWarnings = FALSE)

dir.create(file.path(getwd (), folder, output_folder,output_plots,output_all), showWarnings = FALSE)
dir.create(file.path(getwd (), folder, output_folder,output_data,output_all), showWarnings = FALSE)


dir.create(file.path(getwd (), folder, output_folder,output_plots,output_mice_line), showWarnings = FALSE)
dir.create(file.path(getwd (), folder, output_folder,output_plots,output_mice_line, output_dyrk1b,folder_to_analyse), showWarnings = FALSE)
dir.create(file.path(getwd (), folder, output_folder,output_plots,output_mice_line, output_kd,folder_to_analyse), showWarnings = FALSE)
dir.create(file.path(getwd (), folder, output_folder,output_plots,output_mice_line, output_dyrk1ae2,folder_to_analyse), showWarnings = FALSE)



getwd()
options(stringsAsFactors = FALSE)

##########################################################################################################################
##########################################################################################################################
######                                                                                 ################################### 
######    Reading in R the working files: mice lines breeding files                    ###################################
######                                                                                 ################################### 
##########################################################################################################################
##########################################################################################################################

####
#########################
## 1)   KinaseDead      #
#########################
#first formating step
#formating the .xls files from the beginning in bash
#awk 'BEGIN { FS=":" } {print $2}' Dyrk1B_pheno.csv | sed -e 's/\?/\-/g'> Kinase_Dead_Strain_F.csv
#awk 'BEGIN { FS=":" } {print $2}' Dyrk1B_pheno.csv | sed -e 's/\?/\-/g'|sed -e 's/\;/\ /g' > Kinase_Dead_Strain_F.txt
#awk 'BEGIN {FS=";"} $1~/^[0-9]+/ {print $0}' Dyrk1B_ics.csv |  sed -e 's/\?/\-/g'| sed -e 's/\;/ /g' > Kinase_Dead_PhenoF.txt
#awk 'BEGIN {FS=";"} $1~/^[0-9]+/ {print $0}' Dyrk1B_ics.csv |  sed -e 's/\?/\-/g'> Kinase_Dead_PhenoF.csv


kd_file1 <- read.table('./input/kinaseDead/Kinase_Dead_Strain_F.csv', header=F, skip=1, sep=";" , 
                       colClasses=c("factor","factor","factor","Date","factor", "factor", "character", "factor", "character", "numeric"), dec=",") #cannot use date as the colclasse because i lose the info of the times as it keeps only the date

colnames(kd_file1) <- c("Origin", "Id", "Born_week", "Born_date", "Sex", "Genotype", "Location", "Status", "Dead_date", "Age_Decimal")
head(kd_file1)

kd_file2 <- read.table('./input/kinaseDead/Kinase_Dead_PhenoF.csv', header=F, skip=1, sep=";", 
                       colClasses=c( "factor","factor","factor","Date","factor", "factor", "character", "factor", "numeric"), dec=",") #cannot use date as the colclasse because i lose the info of the times as it keeps only the date
colnames(kd_file2) <- c("Origin", "Id", "Born_week", "Born_date", "Sex", "Genotype", "Location", "Status",  "Age_Decimal")


kd_file2$Dead_date <- rep("-", nrow(kd_file2)) 
#kd_file2$Dead_date <- as.Date( rep("-", nrow(kd_file2)) ) # is ok even if there is an error: character string is not in a standard unambiguous format as the format instead of being date is character
kd_file2 <- kd_file2[, c(1,2,3,4,5,6,7,8,10,9)] 
head(kd_file2)
#for the next one, ignore the erorrs, lets combine both in one file
kd <- full_join(kd_file1, kd_file2, by=c("Origin","Id","Born_week", "Born_date", "Sex", "Genotype", "Location", "Status", "Dead_date", "Age_Decimal"), sep=";")
head(kd) #seems ok
#dim(kd_file1) #195
#dim(kd_file2) #34
dim(kd) #229 #ook!
#arrange(kd, Id) #to check the ones with 0

#to see duplicated ids numbers, take con mind to do this step after 700 animals they began counting again or even before to make sure
duplicated_Ids <- data.frame(table(kd$Id))
rows_num <- duplicated_Ids[duplicated_Ids$Freq > 1,]
kd[kd$id %in% duplicated_Ids$Var1[duplicated_Ids$Freq > 1],]

###
#########################
#####   statistics      #
#########################

After_wheaning <- ifelse(as.numeric(kd$Age_Decimal) > 4.00, "Alive_weaning", "death_before_weaning")
kd$After_wheaning <- After_wheaning

##filter by status: to see the mortality rates
Status_DC <- filter(kd, Status == "DC") 
dim(Status_DC) #24 11
#defining variables:
number_animals <- nrow(kd)
#total number of females/males born
sexes_total <- as.data.frame(summarise(group_by(kd, Sex),count =n())) 
num_fem <- nrow(filter(kd, Sex=="F" )) 
num_males <-  nrow(filter(kd, Sex=="M" )) 



# A) Mortality
  mortality <- as.data.frame(summarise(group_by(Status_DC, After_wheaning),count =n()))
  deaths_before_weaning <- as.numeric(mortality[1,2]) 
  mortality_sex <- as.data.frame(summarise(group_by(Status_DC, Sex, After_wheaning),count =n())) 
  mortality_sex_simple <- as.data.frame(summarise(group_by(Status_DC, Sex, After_wheaning),count =n()))  #same than above, but i will use this one later
  mortality_sex$Percentage <- mortality_sex$count/deaths_before_weaning*100
  mortality_sex$all_dead_Percentage <- max(cumsum(mortality_sex$count)/number_animals)*100
  write.table(mortality_sex, file=paste0(wd,folder, output_folder,output_data,output_mice_line, output_kd,folder_to_analyse,"stats_KD_total_mortality.txt"),sep = "\t", col.names=TRUE, row.names=T)


#mortality per genotype and sex

mortality_sex_genotype <- as.data.frame(summarise(group_by(Status_DC, Sex, Genotype, After_wheaning),count =n())) 
dead_animals_sex_geno <- inner_join(mortality_sex_genotype,sexes_total, by="Sex")
colnames(dead_animals_sex_geno) <- c("Sex", "Genotype", "After_wheaning", "dead_animals_per_sex_geno", "total_animals_per_sex")

dead_animals_sex_geno$Percentage_geno_sex <- dead_animals_sex_geno$dead_animals_per_sex_geno/dead_animals_sex_geno$total_animals_per_sex*100
dead_animals_sex_geno$Percentage_sex <- dead_animals_sex_geno$dead_animals_per_sex_geno/dead_animals_sex_geno$total_animals_per_sex*100
dead_animals_sex_M <- filter(dead_animals_sex_geno, Sex=="M")
dead_animals_sex_M$dead_animals_sex <- sum(dead_animals_sex_M$dead_animals_per_sex_geno)
dead_animals_sex_M$Percentage_per_sex <- dead_animals_sex_M$dead_animals_sex/dead_animals_sex_M$total_animals_per_sex*100

dead_animals_sex_F <- filter(dead_animals_sex_geno, Sex=="F")
dead_animals_sex_F$dead_animals_sex <- sum(dead_animals_sex_F$dead_animals_per_sex_geno)
dead_animals_sex_F$Percentage_per_sex <- dead_animals_sex_F$dead_animals_sex/dead_animals_sex_F$total_animals_per_sex*100

all_dead_animals_geno_sex <- full_join(dead_animals_sex_M,dead_animals_sex_F, by=c("Sex", "Genotype", "After_wheaning", "dead_animals_per_sex_geno", "total_animals_per_sex","Percentage_geno_sex","Percentage_sex","dead_animals_sex","Percentage_per_sex"))
write.table(all_dead_animals_geno_sex, file=paste0(wd,folder, output_folder,output_data,output_mice_line, output_kd,folder_to_analyse, "stats_KD_mortality_genotype_sex.txt"),sep = "\t", col.names=TRUE, row.names=T)

# we dont really know, as only 2 of 24 animals were genotyped... so


#sumup vector
mortality_by_sexes <- c(total_mortality_rate,fem/num_fem*100,males/num_males*100, all_dead/number_animals*100)


total_born<- as.data.frame(summarise(group_by(kd, Status, After_wheaning),count =n())) #
total_born_sex<- as.data.frame(summarise(group_by(kd,  After_wheaning, Sex),count =n())) #
write.table(total_born_sex, file=paste0(wd,folder, output_folder,output_data,output_mice_line, output_kd,folder_to_analyse, "stats_KD_born_genotype_sex.txt"),sep = "\t", col.names=TRUE, row.names=T)


# B) Transmission

dim(kd) #229 -dead(24)= 205
transmission <- filter(kd, !Genotype == "-") 
dim(transmission) #176 11


### B1 TOTAL TRANSMISSION PER GENOTYPE
  total_transmision<- as.data.frame(summarise(group_by(transmission,  Genotype),count =n())) #
  n_genotyped <- max(cumsum(total_transmision$count))  #176
  total_transmision$Percentage <- total_transmision$count/n_genotyped*100
  total_transmision$Mean <-mean(total_transmision$Percentage)
  write.table(total_transmision, file=paste0(wd,folder, output_folder,output_data,output_mice_line, output_kd,folder_to_analyse, "stats_KD_genotype.txt"),sep = "\t", col.names=TRUE, row.names=T)




### B2 TRANSMISSION  PER SEX PER GENOTYPE
total_transmision_sex<- as.data.frame(summarise(group_by(transmission,  Genotype, Sex),count =n())) #
total_transmision_sex_male <- filter(total_transmision_sex, Sex=="M")
males_genotyped <- max(cumsum(total_transmision_sex_male$count))
total_transmision_sex_male$Percentage <- total_transmision_sex_male$count/males_genotyped*100
total_transmision_sex_male$Mean <- mean(total_transmision_sex_male$Percentage)


total_transmision_sex_fem <- filter(total_transmision_sex, Sex=="F")
fem_genotyped <- max(cumsum(total_transmision_sex_fem$count))
total_transmision_sex_fem$Percentage <- total_transmision_sex_fem$count/fem_genotyped*100
total_transmision_sex_fem$Mean <- mean(total_transmision_sex_fem$Percentage)
transmision_sex_all <- full_join(total_transmision_sex_fem, total_transmision_sex_male, by=c("Genotype", "Sex", "count", "Percentage","Mean"))
write.table(transmision_sex_all, file=paste0(wd,folder, output_folder,output_data,output_mice_line, output_kd,folder_to_analyse, "stats_KD_genotype_sex.txt"),sep = "\t", col.names=TRUE, row.names=T)




### B3 TOTAL TRANSMISSION PER BREEDING CAGE
#list_breeding_cages <- as.list(setNames(total_transmision_bcage$Breeding_couple_Id,total_transmision_bcage$Total_Number_mice_born_from_that_couple ))
#a <- mapply(c, total_transmision_bcage$Breeding_couple_Id, total_transmision_bcage$Total_Number_mice_born_from_that_couple, SIMPLIFY = FALSE)
#vals["key42"]
#######################################################################################################################################

total_transmision_bcage<- as.data.frame(summarise(group_by(transmission,  Origin),count =n()))
colnames(total_transmision_bcage) <- c("Origin", "Number_mice_born")
total_transmision_bcage_geno<- as.data.frame(summarise(group_by(transmission,  Origin, Genotype),count =n())) #
total_transmision_bcage_geno_f <- inner_join(total_transmision_bcage_geno,total_transmision_bcage, by ="Origin" )
total_transmision_bcage_geno_f$Percentages <- total_transmision_bcage_geno_f$count * 100/total_transmision_bcage_geno_f$Number_mice_born


geno_wt <- filter(total_transmision_bcage_geno_f, !Genotype == "R5580wt") 
geno_wt$Mean <- mean(geno_wt$Percentages)
total_transmision_bcage_geno$Mean<- total_transmision_bcage_geno$count * 100/n_genotyped
geno_kd <- filter(total_transmision_bcage_geno_f, !Genotype == "R5580L-/WT") 
geno_kd$Mean <- mean(geno_kd$Percentages)
geno_all <- full_join(geno_wt, geno_kd, by=c("Origin", "Genotype", "count", "Number_mice_born", "Percentages","Mean"))
write.table(geno_all, file=paste0(wd,folder, output_folder,output_data,output_mice_line, output_kd,folder_to_analyse, "stats_KD_breedingCage_genotype.txt"),sep = "\t", col.names=TRUE, row.names=T)



#taking out 11 and 8 breeding cages
total_transmision_bcage_geno_f_11_out <- filter(total_transmision_bcage_geno_f, !Origin == c("11") )
total_transmision_bcage_geno_f_8_11_out <- filter(total_transmision_bcage_geno_f_11_out, !Origin == c("8") )


geno_wt <- filter(total_transmision_bcage_geno_f_8_11_out, !Genotype == "R5580wt") 
geno_wt$Mean <- mean(geno_wt$Percentages)
total_transmision_bcage_geno_f_8_11_out$Mean<- total_transmision_bcage_geno_f_8_11_out$count * 100/n_genotyped
geno_kd <- filter(total_transmision_bcage_geno_f_8_11_out, !Genotype == "R5580L-/WT") 
geno_kd$Mean <- mean(geno_kd$Percentages)
geno_all <- full_join(geno_wt, geno_kd, by=c("Origin", "Genotype", "count", "Number_mice_born", "Percentages","Mean"))
write.table(geno_all, file=paste0(wd,folder, output_folder,output_data,output_mice_line, output_kd,folder_to_analyse, "stats_KD_breedingCage_genotypewo11_8.txt"),sep = "\t", col.names=TRUE, row.names=T)



#grouped by sex
total_transmision_bcage_sex<- as.data.frame(summarise(group_by(transmission,  Origin, Sex),count =n())) #

total_transmision_bcage_sex_f <- inner_join(total_transmision_bcage_sex,total_transmision_bcage, by ="Origin" )
total_transmision_bcage_sex_f$percentage <- total_transmision_bcage_sex_f$count/total_transmision_bcage_sex_f$Number_mice_born*100
male_total_transmision_bcage_sex_f <- filter(total_transmision_bcage_sex_f, Sex=="M")
male_total_transmision_bcage_sex_f$Mean <- mean(male_total_transmision_bcage_sex_f$percentage)

fem_total_transmision_bcage_sex_f <- filter(total_transmision_bcage_sex_f, Sex=="F")
fem_total_transmision_bcage_sex_f$Mean <- mean(fem_total_transmision_bcage_sex_f$percentage)
all_total_transmision_bcage_sex_f <- full_join(male_total_transmision_bcage_sex_f, fem_total_transmision_bcage_sex_f, by= c("Origin", "Sex", "count", "Number_mice_born", "percentage","Mean"))
write.table(all_total_transmision_bcage_sex_f, file=paste0(wd,folder, output_folder,output_data,output_mice_line, output_kd,folder_to_analyse, "stats_KD_breedingCage_sex.txt"),sep = "\t", col.names=TRUE, row.names=T)



#grouped by genotype and sex

total_transmision_bcage_sex_geno<- as.data.frame(summarise(group_by(transmission,  Origin,Genotype, Sex),count =n())) #
total_transmision_bcage_sex_geno_f <- inner_join(total_transmision_bcage_sex_geno,total_transmision_bcage, by ="Origin" )
total_transmision_bcage_sex_geno_f$Percentage <- total_transmision_bcage_sex_geno_f$count/ total_transmision_bcage_sex_geno_f$Number_mice_born*100


wt<- filter(total_transmision_bcage_sex_geno_f, Genotype=="R5580wt")
wt_fem <- filter(wt, Sex=="F")
wt_fem$Mean <- mean(wt_fem$Percentage)
wt_male<- filter(wt, Sex=="M")
wt_male$Mean <- mean(wt_male$Percentage)
all_sex_wt <- full_join(wt_fem, wt_male, by=c("Origin", "Genotype", "Sex", "count", "Number_mice_born", "Percentage","Mean"))


kds <- filter(total_transmision_bcage_sex_geno_f, Genotype=="R5580L-/WT")
kd_fem <- filter(kds, Sex=="F")
kd_fem$Mean <- mean(kd_fem$Percentage)
kd_male <- filter(kds, Sex=="M")
kd_male$Mean <- mean(kd_male$Percentage)
all_sex_kd <- full_join(kd_fem, kd_male, by=c("Origin", "Genotype", "Sex", "count", "Number_mice_born", "Percentage","Mean"))


all_sex_gen <- full_join(all_sex_wt, all_sex_kd, by=c("Origin", "Genotype", "Sex", "count", "Number_mice_born", "Percentage","Mean"))
write.table(all_sex_gen, file=paste0(wd,folder, output_folder,output_data,output_mice_line, output_kd,folder_to_analyse, "stats_KD_breedingCage_genotype_sex.txt"),sep = "\t", col.names=TRUE, row.names=T)


### end for KD  ###

##




####
#########################
## 2)   Dyrk1B          #
#########################

#first formating step
#final orders
#formating the .xls files from the beginning in bash
#awk 'BEGIN { FS=":" } {print $2}' Dyrk1B_pheno.csv | sed -e 's/\?/\-/g'> Dyrk1B_pheno_F.csv
#awk 'BEGIN { FS=":" } {print $2}' Dyrk1B_pheno.csv | sed -e 's/\?/\-/g'|sed -e 's/\;/\ /g' > Dyrk1B_pheno_F.txt
#awk 'BEGIN {FS=";"} $1~/^[0-9]+/ {print $0}' Dyrk1B_ics.csv |  sed -e 's/\?/\-/g'| sed -e 's/\;/ /g' > Dyrk1B_ics_F.txt
#awk 'BEGIN {FS=";"} $1~/^[0-9]+/ {print $0}' Dyrk1B_ics.csv |  sed -e 's/\?/\-/g'> Dyrk1B_ics_F.csv




kd_file1 <- read.table('./input/dyrk1b/Dyrk1B_ics_F.csv', header=F, skip=0, sep=";" , 
                       colClasses=c("factor","factor","factor","Date","factor", "character", "factor", "character", "numeric"), dec=",") #cannot use date as the colclasse because i lose the info of the times as it keeps only the date
colnames(kd_file1) <- c("Origin", "Id", "Born_week", "Born_date", "Sex","Genotype", "Location", "Status", "Age_Decimal")
head(kd_file1)

kd_file2 <- read.table('./input/dyrk1b/Dyrk1B_pheno_F.csv', header=F, skip=0, sep=";", 
                       colClasses=c( "factor","factor","factor","Date","factor", "factor", "character", "factor", "numeric"), dec=",") #cannot use date as the colclasse because i lose the info of the times as it keeps only the date
colnames(kd_file2) <- c("Origin", "Id", "Born_week", "Born_date", "Sex", "Genotype", "Location", "Status",  "Age_Decimal")

head(kd_file2)
#for the next one, ignore the erorrs, lets combine both in one file
kd <- full_join(kd_file1, kd_file2, by=c("Origin","Id","Born_week", "Born_date", "Sex", "Genotype", "Location", "Status", "Age_Decimal"), sep=";")
head(kd) #seems ok
dim(kd_file1) #658
dim(kd_file2) #23
dim(kd) #681 #ook!
#arrange(kd, Id) #to check the ones with 0

#to see duplicated ids numbers, take con mind to do this step after 700 animals they began counting again or even before to make sure
duplicated_Ids <- data.frame(table(kd$Id))
rows_num <- duplicated_Ids[duplicated_Ids$Freq > 1,] # 2 animals have the same ID
#kd[kd$id %in% duplicated_Ids$Var1[duplicated_Ids$Freq > 1],]
kd[kd$Id=="321",] 
#12       8 321     27/16 0006-07-16   M     dead        0     DC           2
#416     10 321     31/16 0003-08-16   F    WT/WT     1004     SC          10
###
#########################
#####   statistics      #
#########################
After_wheaning <- ifelse(as.numeric(kd$Genotype) > L3/WT, "Alive_weaning", "death_before_weaning")

kd$Genotype<- sub("P4181L3/WT","L3/WT", kd$Genotype )
kd$Genotype<- sub("P4181WT/WT","WT/WT", kd$Genotype )

After_wheaning <- ifelse(as.numeric(kd$Age_Decimal) > 4.00, "Alive_weaning", "death_before_weaning")
kd$After_wheaning <- After_wheaning

##filter by status: to see the mortality rates
Status_DC <- filter(kd, Status == "DC") 
dim(Status_DC) #81 10
#defining variables:
number_animals <- nrow(kd) #681
#total number of females/males born
sexes_total <- as.data.frame(summarise(group_by(kd, Sex),count =n())) 
num_fem <- nrow(filter(kd, Sex=="F" )) 
num_males <-  nrow(filter(kd, Sex=="M" )) 
num_unk <-  nrow(filter(kd, Sex=="-" )) 


# A) Mortality
mortality <- as.data.frame(summarise(group_by(Status_DC, After_wheaning),count =n()))
deaths_before_weaning <- as.numeric(filter(mortality, After_wheaning=="death_before_weaning" )[1,2]) #72
mortality_sex <- as.data.frame(summarise(group_by(Status_DC, Sex, After_wheaning),count =n())) 
mortality_sex_simple <- as.data.frame(summarise(group_by(Status_DC, Sex, After_wheaning),count =n()))  #same than above, but i will use this one later
mortality_sex$Percentage <- mortality_sex$count/deaths_before_weaning*100
mortality_sex$all_dead_Percentage <- max(cumsum(mortality_sex$count)/number_animals)*100
write.table(mortality_sex, file=paste0(wd,folder,output_folder,output_data,output_mice_line, output_dyrk1b,folder_to_analyse, "stats_dyrk1b_total_mortality.txt"),sep = "\t", col.names=TRUE, row.names=T)


#mortality per genotype and sex

mortality_sex_genotype <- as.data.frame(summarise(group_by(Status_DC, Sex, Genotype, After_wheaning),count =n())) 
dead_animals_sex_geno <- inner_join(mortality_sex_genotype,sexes_total, by="Sex")
colnames(dead_animals_sex_geno) <- c("Sex", "Genotype", "After_wheaning", "dead_animals_per_sex_geno", "total_animals_per_sex")

dead_animals_sex_geno$Percentage_geno_sex <- dead_animals_sex_geno$dead_animals_per_sex_geno/dead_animals_sex_geno$total_animals_per_sex*100
dead_animals_sex_geno$Percentage_sex <- dead_animals_sex_geno$dead_animals_per_sex_geno/dead_animals_sex_geno$total_animals_per_sex*100
dead_animals_sex_M <- filter(dead_animals_sex_geno, Sex=="M")
dead_animals_sex_M$dead_animals_sex <- sum(dead_animals_sex_M$dead_animals_per_sex_geno)
dead_animals_sex_M$Percentage_per_sex <- dead_animals_sex_M$dead_animals_sex/dead_animals_sex_M$total_animals_per_sex*100

dead_animals_sex_F <- filter(dead_animals_sex_geno, Sex=="F")
dead_animals_sex_F$dead_animals_sex <- sum(dead_animals_sex_F$dead_animals_per_sex_geno)
dead_animals_sex_F$Percentage_per_sex <- dead_animals_sex_F$dead_animals_sex/dead_animals_sex_F$total_animals_per_sex*100

all_dead_animals_geno_sex <- full_join(dead_animals_sex_M,dead_animals_sex_F, by=c("Sex", "Genotype", "After_wheaning", "dead_animals_per_sex_geno", "total_animals_per_sex","Percentage_geno_sex","Percentage_sex","dead_animals_sex","Percentage_per_sex"))
write.table(all_dead_animals_geno_sex, file=paste0(wd,folder, folder_to_analyse,output_folder,output_data,output_mice_line, output_dyrk1b, "stats_dyrk1b_mortality_genotype_sex.txt"),sep = "\t", col.names=TRUE, row.names=T)


total_born<- as.data.frame(summarise(group_by(kd, Status, After_wheaning),count =n())) #
total_born_sex<- as.data.frame(summarise(group_by(kd,  After_wheaning, Sex),count =n())) #
write.table(total_born_sex, file=paste0(wd,folder,output_folder,output_data,output_mice_line, output_dyrk1b,folder_to_analyse,"stats_dyrk1b_born_genotype_sex.txt"),sep = "\t", col.names=TRUE, row.names=T)



# B) Transmission

dim(kd)  #681
transmission <- filter(kd, !Genotype == "-") 
dim(transmission) #471 11


### B1 TOTAL TRANSMISSION PER GENOTYPE
total_transmision<- as.data.frame(summarise(group_by(transmission,  Genotype),count =n())) #
off <- c(1:3,5) #careful with this step... depends on the data
total_transmision <- total_transmision[-off,] 
n_genotyped <- max(cumsum(total_transmision$count))  #176
total_transmision$Percentage <- total_transmision$count/n_genotyped*100
total_transmision$Mean <-mean(total_transmision$Percentage)
write.table(total_transmision, file=paste0(wd,folder,output_folder,output_data,output_mice_line, output_dyrk1b,folder_to_analyse, "stats_dyrk1b_genotype.txt"),sep = "\t", col.names=TRUE, row.names=T)




### B2 TRANSMISSION  PER SEX PER GENOTYPE
total_transmision_sex<- as.data.frame(summarise(group_by(transmission,  Genotype, Sex),count =n())) #
off <- c(1:5,8) 
total_transmision_sex <- total_transmision_sex[-off,] 
total_transmision_sex_male <- filter(total_transmision_sex, Sex=="M")
males_genotyped <- max(cumsum(total_transmision_sex_male$count))
total_transmision_sex_male$Percentage <- total_transmision_sex_male$count/males_genotyped*100
total_transmision_sex_male$Mean <- mean(total_transmision_sex_male$Percentage)


total_transmision_sex_fem <- filter(total_transmision_sex, Sex=="F")
fem_genotyped <- max(cumsum(total_transmision_sex_fem$count))
total_transmision_sex_fem$Percentage <- total_transmision_sex_fem$count/fem_genotyped*100
total_transmision_sex_fem$Mean <- mean(total_transmision_sex_fem$Percentage)
transmision_sex_all <- full_join(total_transmision_sex_fem, total_transmision_sex_male, by=c("Genotype", "Sex", "count", "Percentage","Mean"))
write.table(transmision_sex_all, file=paste0(wd,folder,output_folder,output_data,output_mice_line, output_dyrk1b,folder_to_analyse, "stats_dyrk1b_genotype_sex.txt"),sep = "\t", col.names=TRUE, row.names=T)




### B3 TOTAL TRANSMISSION PER BREEDING CAGE
#list_breeding_cages <- as.list(setNames(total_transmision_bcage$Breeding_couple_Id,total_transmision_bcage$Total_Number_mice_born_from_that_couple ))
#a <- mapply(c, total_transmision_bcage$Breeding_couple_Id, total_transmision_bcage$Total_Number_mice_born_from_that_couple, SIMPLIFY = FALSE)
#vals["key42"]
#######################################################################################################################################

total_transmision_bcage<- as.data.frame(summarise(group_by(transmission,  Origin),count =n()))
colnames(total_transmision_bcage) <- c("Origin", "Number_mice_born")
total_transmision_bcage_geno<- as.data.frame(summarise(group_by(transmission,  Origin, Genotype),count =n())) #
total_transmision_bcage_geno_f <- inner_join(total_transmision_bcage_geno,total_transmision_bcage, by ="Origin" )
total_transmision_bcage_geno_f$Percentages <- total_transmision_bcage_geno_f$count * 100/total_transmision_bcage_geno_f$Number_mice_born


geno_wt <- filter(total_transmision_bcage_geno_f, !Genotype == "R5580wt") 
geno_wt$Mean <- mean(geno_wt$Percentages)
total_transmision_bcage_geno$Mean<- total_transmision_bcage_geno$count * 100/n_genotyped
geno_kd <- filter(total_transmision_bcage_geno_f, !Genotype == "R5580L-/WT") 
geno_kd$Mean <- mean(geno_kd$Percentages)
geno_all <- full_join(geno_wt, geno_kd, by=c("Origin", "Genotype", "count", "Number_mice_born", "Percentages","Mean"))
write.table(geno_all, file=paste0(wd,folder,output_folder,output_data,output_mice_line, output_dyrk1b,folder_to_analyse, "stats_dyrk1b_breedingCage_genotype.txt"),sep = "\t", col.names=TRUE, row.names=T)



#taking out 11 and 8 breeding cages
total_transmision_bcage_geno_f_11_out <- filter(total_transmision_bcage_geno_f, !Origin == c("11") )
total_transmision_bcage_geno_f_8_11_out <- filter(total_transmision_bcage_geno_f_11_out, !Origin == c("8") )


geno_wt <- filter(total_transmision_bcage_geno_f_8_11_out, !Genotype == "R5580wt") 
geno_wt$Mean <- mean(geno_wt$Percentages)
total_transmision_bcage_geno_f_8_11_out$Mean<- total_transmision_bcage_geno_f_8_11_out$count * 100/n_genotyped
geno_kd <- filter(total_transmision_bcage_geno_f_8_11_out, !Genotype == "R5580L-/WT") 
geno_kd$Mean <- mean(geno_kd$Percentages)
geno_all <- full_join(geno_wt, geno_kd, by=c("Origin", "Genotype", "count", "Number_mice_born", "Percentages","Mean"))
write.table(geno_all, file=paste0(wd,folder,output_folder,output_data,output_mice_line, output_dyrk1b,folder_to_analyse, "stats_dyrk1b_breedingCage_genotypewo11_8.txt"),sep = "\t", col.names=TRUE, row.names=T)



#grouped by sex
total_transmision_bcage_sex<- as.data.frame(summarise(group_by(transmission,  Origin, Sex),count =n())) #

total_transmision_bcage_sex_f <- inner_join(total_transmision_bcage_sex,total_transmision_bcage, by ="Origin" )
total_transmision_bcage_sex_f$percentage <- total_transmision_bcage_sex_f$count/total_transmision_bcage_sex_f$Number_mice_born*100
male_total_transmision_bcage_sex_f <- filter(total_transmision_bcage_sex_f, Sex=="M")
male_total_transmision_bcage_sex_f$Mean <- mean(male_total_transmision_bcage_sex_f$percentage)

fem_total_transmision_bcage_sex_f <- filter(total_transmision_bcage_sex_f, Sex=="F")
fem_total_transmision_bcage_sex_f$Mean <- mean(fem_total_transmision_bcage_sex_f$percentage)
all_total_transmision_bcage_sex_f <- full_join(male_total_transmision_bcage_sex_f, fem_total_transmision_bcage_sex_f, by= c("Origin", "Sex", "count", "Number_mice_born", "percentage","Mean"))
write.table(all_total_transmision_bcage_sex_f, file=paste0(wd,folder,output_folder,output_data,output_mice_line, output_dyrk1b,folder_to_analyse,"stats_dyrk1b_breedingCage_sex.txt"),sep = "\t", col.names=TRUE, row.names=T)



#grouped by genotype and sex

total_transmision_bcage_sex_geno<- as.data.frame(summarise(group_by(transmission,  Origin,Genotype, Sex),count =n())) #
total_transmision_bcage_sex_geno_f <- inner_join(total_transmision_bcage_sex_geno,total_transmision_bcage, by ="Origin" )
total_transmision_bcage_sex_geno_f$Percentage <- total_transmision_bcage_sex_geno_f$count/ total_transmision_bcage_sex_geno_f$Number_mice_born*100


wt<- filter(total_transmision_bcage_sex_geno_f, Genotype=="R5580wt")
wt_fem <- filter(wt, Sex=="F")
wt_fem$Mean <- mean(wt_fem$Percentage)
wt_male<- filter(wt, Sex=="M")
wt_male$Mean <- mean(wt_male$Percentage)
all_sex_wt <- full_join(wt_fem, wt_male, by=c("Origin", "Genotype", "Sex", "count", "Number_mice_born", "Percentage","Mean"))


kds <- filter(total_transmision_bcage_sex_geno_f, Genotype=="R5580L-/WT")
kd_fem <- filter(kds, Sex=="F")
kd_fem$Mean <- mean(kd_fem$Percentage)
kd_male <- filter(kds, Sex=="M")
kd_male$Mean <- mean(kd_male$Percentage)
all_sex_kd <- full_join(kd_fem, kd_male, by=c("Origin", "Genotype", "Sex", "count", "Number_mice_born", "Percentage","Mean"))


all_sex_gen <- full_join(all_sex_wt, all_sex_kd, by=c("Origin", "Genotype", "Sex", "count", "Number_mice_born", "Percentage","Mean"))
write.table(all_sex_gen, file=paste0(wd,folder,output_folder,output_data,output_mice_line, output_dyrk1b,folder_to_analyse, "stats_dyrk1b_breedingCage_genotype_sex.txt"),sep = "\t", col.names=TRUE, row.names=T)


### end for KD  ###

##

#final orders

awk 'BEGIN { FS=":" } {print $2}' Dyrk1B_pheno.csv | sed -e 's/\?/\-/g'> Dyrk1B_pheno_F.csv
awk 'BEGIN { FS=":" } {print $2}' Dyrk1B_pheno.csv | sed -e 's/\?/\-/g'|sed -e 's/\;/\ /g' > Dyrk1B_pheno_F.txt
awk 'BEGIN {FS=";"} $1~/^[0-9]+/ {print $0}' Dyrk1B_ics.csv |  sed -e 's/\?/\-/g'| sed -e 's/\;/ /g' > Dyrk1B_ics_F.txt
awk 'BEGIN {FS=";"} $1~/^[0-9]+/ {print $0}' Dyrk1B_ics.csv |  sed -e 's/\?/\-/g'> Dyrk1B_ics_F.csv



####
#########################
## 1)   DYRK1a2         #
#########################
#1 0  06/10  10/02/10 -  0 DC 1
#1;34; 12/10 ;24/03/10;F;189N3Tg/+;0;SC;15,71
#Origin  #	Born the	Born the (date)	Sex	Genotype	Location	Status	Age

kd <- read.table('./input/dyrk1b//189N3_mouse_strain_TgDyrk1A2_F.csv', header=F, skip=0, sep=";" , 
                       colClasses=c("factor","factor","factor","Date","factor", "factor", "factor", "factor", "numeric"), dec=",") #cannot use date as the colclasse because i lose the info of the times as it keeps only the date
colnames(kd) <- c("Origin", "Id", "Born_week", "Born_date", "Sex", "Genotype", "Location", "Status", "Age_Decimal")
head(kd)
dim(kd) #2785    9 #ook!

#to see duplicated ids numbers, take con mind to do this step after 700 animals they began counting again or even before to make sure
duplicated_Ids <- data.frame(table(kd$Id))
rows_num <- duplicated_Ids[duplicated_Ids$Freq > 1,] #on this case, the line has going on for many years, there are many repeated ids, we need to create an unique id
kd[kd$id %in% duplicated_Ids$Var1[duplicated_Ids$Freq > 1],]

kd$uniq_id <- make.names(kd[,2], unique=T)
rownames(kd) <- kd$uniq_id

head(kd)
###
#########################
#####   statistics      #
#########################

After_wheaning <- ifelse(as.numeric(kd$Age_Decimal) > 4.00, "Alive_weaning", "death_before_weaning")
kd$After_wheaning <- After_wheaning

##filter by status: to see the mortality rates
Status_DC <- filter(kd, Status == "DC") 
dim(Status_DC) #669  11
#defining variables:
number_animals <- nrow(kd) #2785
#total number of females/males born
sexes_total <- as.data.frame(summarise(group_by(kd, Sex),count =n())) 
num_fem <- nrow(filter(kd, Sex=="F" )) 
num_males <-  nrow(filter(kd, Sex=="M" )) 


# A) Mortality

####cricricri
mortality <- as.data.frame(summarise(group_by(Status_DC, After_wheaning),count =n()))
deaths_before_weaning <- as.numeric(mortality[1,2])
mortality_sex <- as.data.frame(summarise(group_by(Status_DC, Sex, After_wheaning),count =n())) 
mortality_sex_simple <- as.data.frame(summarise(group_by(Status_DC, Sex, After_wheaning),count =n()))  #same than above, but i will use this one later
mortality_sex$Percentage <- mortality_sex$count/deaths_before_weaning*100
mortality_sex$all_dead_Percentage <- max(cumsum(mortality_sex$count)/number_animals)*100
write.table(mortality_sex, file=paste0(wd,folder, output_folder,output_data,output_mice_line,output_dyrk1ae2,folder_to_analyse, "stats_tgDyrk1a2_total_mortality.txt"),sep = "\t", col.names=TRUE, row.names=T)


#mortality per genotype and sex

mortality_sex_genotype <- as.data.frame(summarise(group_by(Status_DC, Sex, Genotype, After_wheaning),count =n())) 
dead_animals_sex_geno <- inner_join(mortality_sex_genotype,sexes_total, by="Sex")
colnames(dead_animals_sex_geno) <- c("Sex", "Genotype", "After_wheaning", "dead_animals_per_sex_geno", "total_animals_per_sex")

dead_animals_sex_geno$Percentage_geno_sex <- dead_animals_sex_geno$dead_animals_per_sex_geno/dead_animals_sex_geno$total_animals_per_sex*100
dead_animals_sex_geno$Percentage_sex <- dead_animals_sex_geno$dead_animals_per_sex_geno/dead_animals_sex_geno$total_animals_per_sex*100
dead_animals_sex_M <- filter(dead_animals_sex_geno, Sex=="M")
dead_animals_sex_M$dead_animals_sex <- sum(dead_animals_sex_M$dead_animals_per_sex_geno)
dead_animals_sex_M$Percentage_per_sex <- dead_animals_sex_M$dead_animals_sex/dead_animals_sex_M$total_animals_per_sex*100

dead_animals_sex_F <- filter(dead_animals_sex_geno, Sex=="F")
dead_animals_sex_F$dead_animals_sex <- sum(dead_animals_sex_F$dead_animals_per_sex_geno)
dead_animals_sex_F$Percentage_per_sex <- dead_animals_sex_F$dead_animals_sex/dead_animals_sex_F$total_animals_per_sex*100

all_dead_animals_geno_sex <- full_join(dead_animals_sex_M,dead_animals_sex_F, by=c("Sex", "Genotype", "After_wheaning", "dead_animals_per_sex_geno", "total_animals_per_sex","Percentage_geno_sex","Percentage_sex","dead_animals_sex","Percentage_per_sex"))
write.table(all_dead_animals_geno_sex, file=paste0(wd,folder, output_folder,output_data,output_mice_line,output_dyrk1ae2,folder_to_analyse,  "stats_tgDyrk1a2_mortality_genotype_sex.txt"),sep = "\t", col.names=TRUE, row.names=T)

# we dont really know, as only 2 of 24 animals were genotyped... so


#sumup vector
mortality_by_sexes <- c(total_mortality_rate,fem/num_fem*100,males/num_males*100, all_dead/number_animals*100)


total_born<- as.data.frame(summarise(group_by(kd, Status, After_wheaning),count =n())) #
total_born_sex<- as.data.frame(summarise(group_by(kd,  After_wheaning, Sex),count =n())) #
write.table(total_born_sex, file=paste0(wd,folder, output_folder,output_data,output_mice_line,output_dyrk1ae2,folder_to_analyse,  "stats_tgDyrk1a2_born_genotype_sex.txt"),sep = "\t", col.names=TRUE, row.names=T)



# B) Transmission

dim(kd) #229 -dead(24)= 205
transmission <- filter(kd, !Genotype == "-") 
dim(transmission) #176 11


### B1 TOTAL TRANSMISSION PER GENOTYPE
total_transmision<- as.data.frame(summarise(group_by(transmission,  Genotype),count =n())) #
n_genotyped <- max(cumsum(total_transmision$count))  #176
total_transmision$Percentage <- total_transmision$count/n_genotyped*100
total_transmision$Mean <-mean(total_transmision$Percentage)
write.table(total_transmision, file=paste0(wd,folder, output_folder,output_data,output_mice_line,output_dyrk1ae2,folder_to_analyse,  "stats_tgDyrk1a2_genotype.txt"),sep = "\t", col.names=TRUE, row.names=T)




### B2 TRANSMISSION  PER SEX PER GENOTYPE
total_transmision_sex<- as.data.frame(summarise(group_by(transmission,  Genotype, Sex),count =n())) #
total_transmision_sex_male <- filter(total_transmision_sex, Sex=="M")
males_genotyped <- max(cumsum(total_transmision_sex_male$count))
total_transmision_sex_male$Percentage <- total_transmision_sex_male$count/males_genotyped*100
total_transmision_sex_male$Mean <- mean(total_transmision_sex_male$Percentage)


total_transmision_sex_fem <- filter(total_transmision_sex, Sex=="F")
fem_genotyped <- max(cumsum(total_transmision_sex_fem$count))
total_transmision_sex_fem$Percentage <- total_transmision_sex_fem$count/fem_genotyped*100
total_transmision_sex_fem$Mean <- mean(total_transmision_sex_fem$Percentage)
transmision_sex_all <- full_join(total_transmision_sex_fem, total_transmision_sex_male, by=c("Genotype", "Sex", "count", "Percentage","Mean"))
write.table(transmision_sex_all, file=paste0(wd,folder, output_folder,output_data,output_mice_line,output_dyrk1ae2,folder_to_analyse,  "stats_tgDyrk1a2_genotype_sex.txt"),sep = "\t", col.names=TRUE, row.names=T)




### B3 TOTAL TRANSMISSION PER BREEDING CAGE
#######################################################################################################################################

total_transmision_bcage<- as.data.frame(summarise(group_by(transmission,  Origin),count =n()))
colnames(total_transmision_bcage) <- c("Origin", "Number_mice_born")
total_transmision_bcage_geno<- as.data.frame(summarise(group_by(transmission,  Origin, Genotype),count =n())) #
total_transmision_bcage_geno_f <- inner_join(total_transmision_bcage_geno,total_transmision_bcage, by ="Origin" )
total_transmision_bcage_geno_f$Percentages <- total_transmision_bcage_geno_f$count * 100/total_transmision_bcage_geno_f$Number_mice_born


geno_wt <- filter(total_transmision_bcage_geno_f, !Genotype == "R5580wt") 
geno_wt$Mean <- mean(geno_wt$Percentages)
total_transmision_bcage_geno$Mean<- total_transmision_bcage_geno$count * 100/n_genotyped
geno_kd <- filter(total_transmision_bcage_geno_f, !Genotype == "R5580L-/WT") 
geno_kd$Mean <- mean(geno_kd$Percentages)
geno_all <- full_join(geno_wt, geno_kd, by=c("Origin", "Genotype", "count", "Number_mice_born", "Percentages","Mean"))
write.table(geno_all, file=paste0(wd,folder, output_folder,output_data,output_mice_line,output_dyrk1ae2,folder_to_analyse,  "stats_tgDyrk1a2_breedingCage_genotype.txt"),sep = "\t", col.names=TRUE, row.names=T)



#taking out 11 and 8 breeding cages
total_transmision_bcage_geno_f_11_out <- filter(total_transmision_bcage_geno_f, !Origin == c("11") )
total_transmision_bcage_geno_f_8_11_out <- filter(total_transmision_bcage_geno_f_11_out, !Origin == c("8") )


geno_wt <- filter(total_transmision_bcage_geno_f_8_11_out, !Genotype == "R5580wt") 
geno_wt$Mean <- mean(geno_wt$Percentages)
total_transmision_bcage_geno_f_8_11_out$Mean<- total_transmision_bcage_geno_f_8_11_out$count * 100/n_genotyped
geno_kd <- filter(total_transmision_bcage_geno_f_8_11_out, !Genotype == "R5580L-/WT") 
geno_kd$Mean <- mean(geno_kd$Percentages)
geno_all <- full_join(geno_wt, geno_kd, by=c("Origin", "Genotype", "count", "Number_mice_born", "Percentages","Mean"))
write.table(geno_all, file=paste0(wd,folder, output_folder,output_data,output_mice_line,output_dyrk1ae2,folder_to_analyse,  "stats_tgDyrk1a2_breedingCage_genotypewo11_8.txt"),sep = "\t", col.names=TRUE, row.names=T)



#grouped by sex
total_transmision_bcage_sex<- as.data.frame(summarise(group_by(transmission,  Origin, Sex),count =n())) #

total_transmision_bcage_sex_f <- inner_join(total_transmision_bcage_sex,total_transmision_bcage, by ="Origin" )
total_transmision_bcage_sex_f$percentage <- total_transmision_bcage_sex_f$count/total_transmision_bcage_sex_f$Number_mice_born*100
male_total_transmision_bcage_sex_f <- filter(total_transmision_bcage_sex_f, Sex=="M")
male_total_transmision_bcage_sex_f$Mean <- mean(male_total_transmision_bcage_sex_f$percentage)

fem_total_transmision_bcage_sex_f <- filter(total_transmision_bcage_sex_f, Sex=="F")
fem_total_transmision_bcage_sex_f$Mean <- mean(fem_total_transmision_bcage_sex_f$percentage)
all_total_transmision_bcage_sex_f <- full_join(male_total_transmision_bcage_sex_f, fem_total_transmision_bcage_sex_f, by= c("Origin", "Sex", "count", "Number_mice_born", "percentage","Mean"))
write.table(all_total_transmision_bcage_sex_f, file=paste0(wd,folder, output_folder,output_data,output_mice_line,output_dyrk1ae2,folder_to_analyse,  "stats_tgDyrk1a2_breedingCage_sex.txt"),sep = "\t", col.names=TRUE, row.names=T)



#grouped by genotype and sex

total_transmision_bcage_sex_geno<- as.data.frame(summarise(group_by(transmission,  Origin,Genotype, Sex),count =n())) #
total_transmision_bcage_sex_geno_f <- inner_join(total_transmision_bcage_sex_geno,total_transmision_bcage, by ="Origin" )
total_transmision_bcage_sex_geno_f$Percentage <- total_transmision_bcage_sex_geno_f$count/ total_transmision_bcage_sex_geno_f$Number_mice_born*100


wt<- filter(total_transmision_bcage_sex_geno_f, Genotype=="R5580wt")
wt_fem <- filter(wt, Sex=="F")
wt_fem$Mean <- mean(wt_fem$Percentage)
wt_male<- filter(wt, Sex=="M")
wt_male$Mean <- mean(wt_male$Percentage)
all_sex_wt <- full_join(wt_fem, wt_male, by=c("Origin", "Genotype", "Sex", "count", "Number_mice_born", "Percentage","Mean"))


kds <- filter(total_transmision_bcage_sex_geno_f, Genotype=="R5580L-/WT")
kd_fem <- filter(kds, Sex=="F")
kd_fem$Mean <- mean(kd_fem$Percentage)
kd_male <- filter(kds, Sex=="M")
kd_male$Mean <- mean(kd_male$Percentage)
all_sex_kd <- full_join(kd_fem, kd_male, by=c("Origin", "Genotype", "Sex", "count", "Number_mice_born", "Percentage","Mean"))


all_sex_gen <- full_join(all_sex_wt, all_sex_kd, by=c("Origin", "Genotype", "Sex", "count", "Number_mice_born", "Percentage","Mean"))
write.table(all_sex_gen, file=paste0(wd,folder, output_folder,output_data,output_mice_line,output_dyrk1ae2,folder_to_analyse,  "stats_tgDyrk1a2_breedingCage_genotype_sex.txt"),sep = "\t", col.names=TRUE, row.names=T)


### end for dyrk1a2 ###

###Keep going GuiGui! :)

























