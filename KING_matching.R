# #Run the KING matching
library(gtools)
library(data.table)
library(tidyverse)


#Recode everything to be chr and position

setwd("/lustre/scratch123/hgi/mdt2/projects/ibdgwas_bioresource/mt27/CSI/Scratch/")

cohorts = c("WES2","GWAS1","GWAS2","GWAS3","WGS","IBDBR","b25","b27")
cohorts_combinations = combinations(n = length(cohorts), r = 2, v = cohorts )
#
#Run through all the pairs
for (i in 1:length(cohorts_combinations)) {
  print(i)
  #Regenerate the PLINK file
  system(sprintf("/lustre/scratch123/hgi/mdt2/projects/ibdgwas_bioresource/qz2/software/king -b KING_Plink_Cohorts/%s.bed,KING_Plink_Cohorts/%s.bed --related --degree 1  --prefix KING_Cohort_Results/%s_%s", cohorts_combinations[i,1],cohorts_combinations[i,2],cohorts_combinations[i,1],cohorts_combinations[i,2]))
}
#
#Load in all the matches:
files   = list.files(path = "/lustre/scratch123/hgi/projects/ibdgwas_bioresource/mt27/CSI/Scratch/KING_Cohort_Results/",pattern = "kin0",full.names = T)
files   = files[which(!str_detect(files, "X"))]
matches = rbindlist(lapply(files, fread),fill = T)

duplicate_matches = matches[which(matches$Kinship > 0.34),]

#Remove duplicate rows
duplicate_matches = duplicate_matches[!duplicated(duplicate_matches[,c(1,3)]),]
# Add to the table of the matches (to both the FID1 and FID2) which Cohort they originated from
# Read in all the FAM files for this:
WES_fam   = fread("/lustre/scratch123/hgi/projects/ibdgwas_bioresource/mt27/CSI/Scratch/KING_Plink_Cohorts/WES2.fam")
GWAS1_fam = fread("/lustre/scratch123/hgi/projects/ibdgwas_bioresource/mt27/CSI/Scratch/KING_Plink_Cohorts/GWAS1.fam")
GWAS2_fam = fread("/lustre/scratch123/hgi/projects/ibdgwas_bioresource/mt27/CSI/Scratch/KING_Plink_Cohorts/GWAS2.fam")
GWAS3_fam = fread("/lustre/scratch123/hgi/projects/ibdgwas_bioresource/mt27/CSI/Scratch/KING_Plink_Cohorts/GWAS3.fam")
IBDBR_fam = fread("/lustre/scratch123/hgi/projects/ibdgwas_bioresource/mt27/CSI/Scratch/KING_Plink_Cohorts/IBDBR.fam")
WGS_fam   = fread("/lustre/scratch123/hgi/projects/ibdgwas_bioresource/mt27/CSI/Scratch/KING_Plink_Cohorts/WGS.fam")
WGS_fam   = fread("/lustre/scratch123/hgi/projects/ibdgwas_bioresource/mt27/CSI/Scratch/KING_Plink_Cohorts/WGS.fam")
b25_fam = fread("/lustre/scratch123/hgi/projects/ibdgwas_bioresource/mt27/CSI/Scratch/KING_Plink_Cohorts/b25.fam")
b27_fam = fread("/lustre/scratch123/hgi/projects/ibdgwas_bioresource/mt27/CSI/Scratch/KING_Plink_Cohorts/b27.fam")




#Go through the seperate cohorts for both FID1 and FID2 and assign a cohort ID to those indexes
duplicate_matches$FID2_Cohort_B                   = NA
duplicate_matches$FID1_Cohort_A                   = NA

#WES
WES_indexes_FID1                                  =  which(duplicate_matches$FID1 %in% WES_fam$V1)
duplicate_matches$FID1_Cohort_A[WES_indexes_FID1] = "WES"
WES_indexes_FID2                                  =  which(duplicate_matches$FID2 %in% WES_fam$V1)
duplicate_matches$FID2_Cohort_B[WES_indexes_FID2] = "WES"

#GWAS1
GWAS1_indexes_FID1                                    =  which(duplicate_matches$FID1 %in% GWAS1_fam$V1)
duplicate_matches$FID1_Cohort_A[GWAS1_indexes_FID1]   = "GWAS1"
GWAS1_indexes_FID2                                    =  which(duplicate_matches$FID2 %in% GWAS1_fam$V1)
duplicate_matches$FID2_Cohort_B[GWAS1_indexes_FID2]   = "GWAS1"

#GWAS2
GWAS2_indexes_FID1                                    =  which(duplicate_matches$FID1 %in% GWAS2_fam$V1)
duplicate_matches$FID1_Cohort_A[GWAS2_indexes_FID1]   = "GWAS2"
GWAS2_indexes_FID2                                    =  which(duplicate_matches$FID2 %in% GWAS2_fam$V1)
duplicate_matches$FID2_Cohort_B[GWAS2_indexes_FID2]   = "GWAS2"

#GWAS3
GWAS3_indexes_FID1                                    =  which(duplicate_matches$FID1 %in% GWAS3_fam$V1)
duplicate_matches$FID1_Cohort_A[GWAS3_indexes_FID1]   = "GWAS3"
GWAS3_indexes_FID2                                    =  which(duplicate_matches$FID2 %in% GWAS3_fam$V1)
duplicate_matches$FID2_Cohort_B[GWAS3_indexes_FID2]   = "GWAS3"

#IBDBR
IBDBR_indexes_FID1                                    =  which(duplicate_matches$FID1 %in% IBDBR_fam$V1)
duplicate_matches$FID1_Cohort_A[IBDBR_indexes_FID1]   = "IBDBR"
IBDBR_indexes_FID2                                    =  which(duplicate_matches$FID2 %in% IBDBR_fam$V1)
duplicate_matches$FID2_Cohort_B[IBDBR_indexes_FID2]   = "IBDBR"

#WGS
WGS_indexes_FID1                                    =  which(duplicate_matches$FID1 %in% WGS_fam$V1)
duplicate_matches$FID1_Cohort_A[WGS_indexes_FID1]   = "WGS"
WGS_indexes_FID2                                    =  which(duplicate_matches$FID2 %in% WGS_fam$V1)
duplicate_matches$FID2_Cohort_B[WGS_indexes_FID2]   = "WGS"

#b25
b25_indexes_FID1                                    =  which(duplicate_matches$FID1 %in% b25_fam$V1)
duplicate_matches$FID1_Cohort_A[b25_indexes_FID1]   = "b25"
b25_indexes_FID2                                    =  which(duplicate_matches$FID2 %in% b25_fam$V1)
duplicate_matches$FID2_Cohort_B[b25_indexes_FID2]   = "b25"

#b27
b27_indexes_FID1                                    =  which(duplicate_matches$FID1 %in% b27_fam$V1)
duplicate_matches$FID1_Cohort_A[b27_indexes_FID1]   = "b27"
b27_indexes_FID2                                    =  which(duplicate_matches$FID2 %in% b27_fam$V1)
duplicate_matches$FID2_Cohort_B[b27_indexes_FID2]   = "b27"

#Print out a contigency table of the matches
table(duplicate_matches$FID1_Cohort_A, duplicate_matches$FID2_Cohort_B)



#
write.table(duplicate_matches, "matches", quote = F, col.names = T, row.names = F)
# #
# # #Extract out the IBDBR samples from the WES:
# wes_summary = fread("/nfs/users/nfs_m/mt27/all_wes_updated.tsv")
#
# for (i in 1:nrow(duplicate_matches)) {
#   print(i)
#   egan_track = NA
#   if (duplicate_matches[i,]$FID1 %in% wes_summary$EGA_accession_number) {
#     index = which(wes_summary$EGA_accession_number %in% duplicate_matches[i,]$FID1)
#     egan_track = wes_summary$EGA_accession_number[index]
#     if (wes_summary[which(wes_summary$EGA_accession_number == egan_track),]$source == "ibd_bio")
#       if(duplicate_matches[i,]$FID1_Cohort_A == "WES"){
#         duplicate_matches[i,]$FID1_Cohort_A = "WES_IBDBR"
#       }
#       if(duplicate_matches[i,]$FID1_Cohort_A == "WGS"){
#         duplicate_matches[i,]$FID1_Cohort_A = "WGS_IBDBR"
#     }
#
#   }
#   if (duplicate_matches[i,]$FID2 %in% wes_summary$EGA_accession_number) {
#     index = which(wes_summary$EGA_accession_number %in% duplicate_matches[i,]$FID2)
#     egan_track = wes_summary$EGA_accession_number[index]
#     if (wes_summary[which(wes_summary$EGA_accession_number == egan_track),]$source == "ibd_bio")
#       if(duplicate_matches[i,]$FID2_Cohort_B == "WES"){
#         duplicate_matches[i,]$FID2_Cohort_B = "WES_IBDBR"
#       }
#       if(duplicate_matches[i,]$FID2_Cohort_B == "WGS"){
#         duplicate_matches[i,]$FID2_Cohort_B = "WGS_IBDBR"
#     }
#   }
# }
#
# # ## Add things to Carolyn's table:
# # for (i in 1:nrow(spid_linkage_table)) {
# #   print(i)
# #   index = which(wes_summary$Sanger.sample == spid_linkage_table$sanger_sample_id[i])
# #   egan = wes_summary[index,]$EGA_accession_number
# #   duplicate_egan_index = which(duplicate_matches$FID1 == egan)
# #   if(length(duplicate_egan_index) == 0 ){
# #     spid_linkage_table$link_found[i] = "Unsequenced"
# #   }
# #   else{
# #     if(duplicate_matches$FID2_Cohort_B[duplicate_egan_index] == "IBDBR"){
# #       spid_linkage_table$link_found[i] = duplicate_matches$FID2[duplicate_egan_index]
# #     }
# #   }
# #   duplicate_egan_index = which(duplicate_matches$FID2 == egan)
# #   if(length(duplicate_egan_index) == 0 ){
# #     spid_linkage_table$link_found[i] = "Unsequenced"
# #   }
# #   else{
# #     if(duplicate_matches$FID1_Cohort_A[duplicate_egan_index] == "IBDBR"){
# #       spid_linkage_table$link_found[i] = duplicate_matches$ID1[duplicate_egan_index]
# #     }
# #   }
# # }
# # dna_hand_high_confidence = fread("/lustre/scratch123/hgi/mdt2/projects/ibdgwas_bioresource/mt27/CSI/dnahand_out/kinship/kinship_results_minkin_0.45_nsnps_15.csv")
# # dna_hand_high_confidence = dna_hand_high_confidence[which(dna_hand_high_confidence$V4 > 15 & dna_hand_high_confidence$V3 == 0.5),]
# #
# # containsSTDY1 <- grepl("stdy", dna_hand_high_confidence$V1)
# # df_filtered1 <- subset(dna_hand_high_confidence, containsSTDY1)
# #
# # containsSTDY2 <- grepl("stdy", dna_hand_high_confidence$V2)
# # df_filtered2 <- subset(dna_hand_high_confidence, containsSTDY2)
# #
# # sanger = rbind(df_filtered1,df_filtered2)
# # #Add to these the supplier names, sources, EGANs for samples_1 and samples_2
# # sanger = as.data.frame(cbind(sanger$V1,sanger$V2))
# #
# # sanger$sanger_sample_1_supplier = NA
# # sanger$sanger_sample_2_supplier = NA
# #
# # sanger$sanger_sample_1_source = NA
# # sanger$sanger_sample_2_source = NA
# #
# # sanger$sanger_sample_1_EGAN = NA
# # sanger$sanger_sample_2_EGAN = NA
# #
# #
# # for (i in 1:nrow(sanger)) {
# #   print(i)
# #   index = which(wes_summary$Sanger.sample_ed %in% sanger$V1[i])
# #   if(length(index) == 0){
# #     sanger$sanger_sample_1_supplier[i] = "Unsequenced"
# #     sanger$sanger_sample_1_source[i] = "Unsequenced"
# #     sanger$sanger_sample_1_EGAN[i] = "Unsequenced"
# #   }
# #   else{
# #     sanger$sanger_sample_1_supplier[i] = wes_summary$Supplier.name[index]
# #     sanger$sanger_sample_1_source[i] = wes_summary$source[index]
# #     sanger$sanger_sample_1_EGAN[i] = wes_summary$EGA_accession_number[index]
# #   }
# # }
# #
# #
# # for (i in 1:nrow(sanger)) {
# #   print(i)
# #   index = which(wes_summary$Sanger.sample_ed %in% sanger$V2[i])
# #   if(length(index) == 0){
# #     sanger$sanger_sample_2_supplier[i] = "Unsequenced"
# #     sanger$sanger_sample_2_source[i] = "Unsequenced"
# #     sanger$sanger_sample_2_EGAN[i] = "Unsequenced"
# #   }
# #   else{
# #     sanger$sanger_sample_2_supplier[i] = wes_summary$Supplier.name[index]
# #     sanger$sanger_sample_2_source[i] = wes_summary$source[index]
# #     sanger$sanger_sample_2_EGAN[i] = wes_summary$EGA_accession_number[index]
# #   }
# # }
# #
