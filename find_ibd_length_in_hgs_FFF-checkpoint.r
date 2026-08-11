################################
#
#Claudia Moreau
#
#aout 2023
#
#Find IBD length of segments in HGA and B
#
################################


rm(list = ls());
library("GENLIB")
library(stringr)
library(pals)
library(ggplot2)
library(tidyr)
library(reshape2)
library(ape)
library(data.tree)
library("ggdendro")
library("grid")
library("ggpubr")
library(gridExtra)
library(dplyr)

setwd('/lustre03/project/6033529/ARSACS/results/haplotypes')
getwd()

## Your ht file with hgs
ht_file<-'patients_ctrls_hts_clusters_chr13:23790642-24531781.txt'
hts_df <-  read.table(ht_file,header=TRUE,sep="")

carriers_list <- '/lustre03/project/6033529/ARSACS/data/Jasmin/carriers_list.asc'
carriers_df <-  read.table(carriers_list,header=F,sep="")

## fams and mixed up are already replaced
hts_df$sample_no_ht <- gsub('_1','',hts_df$sample)
hts_df$sample_no_ht <- gsub('_2','',hts_df$sample_no_ht)
hts_df$chr <- sapply(str_split(hts_df$sample, '_'), function(x) if(length(x) >= 3) x[[3]] else NA)

hts_df$fam <- sapply(str_split(hts_df$sample,'_'),'[[',1)
## for ctrls
# Ensure there are no NA values in sample_no_ht
hts_df$sample_no_ht[is.na(hts_df$sample_no_ht)] <- ''

# Perform the assignment
hts_df$sample_no_ht <- ifelse(
  substring(hts_df$sample, 1, 3) == '111',
  paste(
    sapply(str_split(hts_df$sample[substring(hts_df$sample, 1, 3) == '111'], '_'), '[[', 1),
    sapply(str_split(hts_df$sample[substring(hts_df$sample, 1, 3) == '111'], '_'), '[[', 1),
    sep = '_'
  ),
  hts_df$sample_no_ht
)

length(hts_df$sample)

## Do we have samples with A and B??
#samples_HGA <- hts_df$sample[hts_df$hg=='hgA']
#samples_HGB <- hts_df$sample[hts_df$hg=='hgB']
#samples_HGU <- hts_df$sample[hts_df$hg=='hgU']

## ibd.gz file
raw_ibd_file<-'/lustre03/project/6033529/ARSACS/data/genotypes2/ARSACS_CaG_Sag_commonsnps_geno0.05_noATGC_refinedIBD.calculated_chr13_window40_scale10_lod3_length0.5.ibd.gz'
raw_ibd_df <- read.table(gzfile(raw_ibd_file))
dimnames(raw_ibd_df)[[2]] <- c('ind1','chr1','ind2','chr2','chr','start','end','LOD','length')
## check
#raw_ibd_df[raw_ibd_df$ind1=='FDM001_JM1051' & raw_ibd_df$ind2=='FDM397_JM1609',]
inds_ibd_file2<-'/lustre03/project/6033529/ARSACS/data/genotypes2/ARSACS_CaG_Sag_commonsnps_geno0.05_noATGC_refinedIBD.calculated_chr13_window40_scale10_lod3_length0.5.hbd.gz'
inds_ibd_df2 <- read.table(gzfile(inds_ibd_file2))
dimnames(inds_ibd_df2)[[2]] <- c('ind1','chr1','ind2','chr2','chr','start','end','LOD','length')

raw_ibd_df2 <- rbind(inds_ibd_df2,raw_ibd_df)

## 2cM
length <- 0.5
raw_ibd_df2 <- raw_ibd_df2[raw_ibd_df2$length>=0.5,]

raw_ibd_df2$new_col1_2 <- paste(raw_ibd_df2$ind1, raw_ibd_df2$chr1, sep = "_")

raw_ibd_df2$new_col3_4 <- paste(raw_ibd_df2$ind2, raw_ibd_df2$chr2, sep = "_")

## here define your haplogroups
raw_ibd_df2$hg1 <- hts_df$ht[match(raw_ibd_df2$new_col1_2,hts_df$sample)]
raw_ibd_df2$hg2 <- hts_df$ht[match(raw_ibd_df2$new_col3_4,hts_df$sample)]
#raw_ibd_df$ht1 <- hts_df$H2N[match(raw_ibd_df$ind1,hts_df$sample_no_ht)]
#raw_ibd_df$ht2 <- hts_df$H2N[match(raw_ibd_df$ind2,hts_df$sample_no_ht)]
raw_ibd_df2[raw_ibd_df2$length<0.5,]
######## IBD files ######### 
ARSACS_str_pos <- 23909171



## put fam num (maybe not for ARSACS...)
#ARSACS_raw_ibd_df$fam1 <- sapply(str_split(ARSACS_raw_ibd_df$ind1,'_'),'[[',1)
#ARSACS_raw_ibd_df$fam2 <- sapply(str_split(ARSACS_raw_ibd_df$ind2,'_'),'[[',1)


## See if hg mean sharing is driven by fams
#fam_intrahg <-  ARSACS_raw_ibd_df[ARSACS_raw_ibd_df$hg1==ARSACS_raw_ibd_df$hg2,]
#fam_intrahg$intrafam <- 0
#fam_intrahg$intrafam[fam_intrahg$fam1==fam_intrahg$fam2 ] <- 1

## mean length in HG



unique_values <- c("0.8MBhaplotypes001", "0.8MBhaplotypes012", "0.8MBhaplotypes096","0.8MBhaplotypes769","0.8MBhaplotypes351")
                  #"0.8MBhaplotypes1308","0.8MBhaplotypes1917", "0.8MBhaplotypes1920", "0.8MBhaplotypes1930", 
                   #"0.8MBhaplotypes1932", "0.8MBhaplotypes1931", "0.8MBhaplotypes1982", "0.8MBhaplotypes1903",
                   #"0.8MBhaplotypes1698", "0.8MBhaplotypes1515", "0.8MBhaplotypes1918")



#unique_values <- c("2MBhaplotypes001", "2MBhaplotypes064", "2MBhaplotypes065", "2MBhaplotypes087",
                   #"2MBhaplotypes160", "2MBhaplotypes211", "2MBhaplotypes277", "2MBhaplotypes399",
                   #"2MBhaplotypes467", "2MBhaplotypes596", "2MBhaplotypes610")
#"2MBhaplotypes1002","2MBhaplotypes1004", "2MBhaplotypes1005", "2MBhaplotypes1010", "2MBhaplotypes1011",
                  #"2MBhaplotypes1012", "2MBhaplotypes1200", "2MBhaplotypes1318", "2MBhaplotypes1520",
                  #"2MBhaplotypes1802", "2MBhaplotypes1848", "2MBhaplotypes2019", "2MBhaplotypes2022",
                  #"2MBhaplotypes2147", "2MBhaplotypes2209", "2MBhaplotypes2544", "2MBhaplotypes2637",
                  #"2MBhaplotypes2836", "2MBhaplotypes2878", "2MBhaplotypes3187", "2MBhaplotypes3258",
                  #"2MBhaplotypes3317", "2MBhaplotypes3631", "2MBhaplotypes610", "2MBhaplotypes938",
                  #"2MBhaplotypes982", "2MBhaplotypes983", "2MBhaplotypes985", "2MBhaplotypes987",
                  #"2MBhaplotypes989", "2MBhaplotypes993", "2MBhaplotypes995")


#unique_values <- c("4MBhaplotypes001", "4MBhaplotypes002", "4MBhaplotypes022", "4MBhaplotypes047",
                   #"4MBhaplotypes057", "4MBhaplotypes135", "4MBhaplotypes169", "4MBhaplotypes187",
                   #"4MBhaplotypes198", "4MBhaplotypes199", "4MBhaplotypes269", "4MBhaplotypes308",
                   #"4MBhaplotypes398", "4MBhaplotypes405", "4MBhaplotypes133", "4MBhaplotypes206",
                   #"4MBhaplotypes237", "4MBhaplotypes282", "4MBhaplotypes364", "4MBhaplotypes474")
                   #"4MBhaplotypes1587", "4MBhaplotypes1773", "4MBhaplotypes1835", "4MBhaplotypes1891",
                   #"4MBhaplotypes2011", "4MBhaplotypes1086", "4MBhaplotypes2173", "4MBhaplotypes2352",
                   #"4MBhaplotypes1101", "4MBhaplotypes2615", "4MBhaplotypes2632", "4MBhaplotypes2730",
                   #"4MBhaplotypes2732", "4MBhaplotypes2735", "4MBhaplotypes2737", "4MBhaplotypes2751",
                   #"4MBhaplotypes2790", "4MBhaplotypes2791", "4MBhaplotypes1001", "4MBhaplotypes2937",
                   #"4MBhaplotypes3036", "4MBhaplotypes3060", "4MBhaplotypes3410","4MBhaplotypes3446",
                   #"4MBhaplotypes3492", "4MBhaplotypes3493", "4MBhaplotypes3494","4MBhaplotypes3503",
                   #"4MBhaplotypes3504", "4MBhaplotypes3507", "4MBhaplotypes3510","4MBhaplotypes1378",
                   #"4MBhaplotypes3761", "4MBhaplotypes3762", "4MBhaplotypes3763","4MBhaplotypes3765",
                   #"4MBhaplotypes3766", "4MBhaplotypes3767", "4MBhaplotypes3771","4MBhaplotypes3772",
                   #"4MBhaplotypes3774", "4MBhaplotypes3775", "4MBhaplotypes3777","4MBhaplotypes3780",
                   #"4MBhaplotypes3782", "4MBhaplotypes4284", "4MBhaplotypes453","4MBhaplotypes4662",
                   #"4MBhaplotypes1530", "4MBhaplotypes632")


# Filter the rows
filtered_table <- raw_ibd_df2 %>%
  filter(hg1 %in% unique_values & hg2 %in% unique_values)

# View the filtered table
print(filtered_table)


# Filter rows based on the list values
filtered_table2 <- filtered_table[
  filtered_table$ind1 %in% carriers_df$V1 | filtered_table$ind2 %in% carriers_df$V1, 
]

## make it square
ARSACS_raw_ibd_df <- filtered_table2[,c(3,4,1,2,5:9,11,10,13,12)]
dimnames(ARSACS_raw_ibd_df)[[2]] <- c('ind2','chr2','ind1','chr1','chr','start','end','LOD','length','new_col3_4','new-col1_2','hg2','hg1')


colnames(ARSACS_raw_ibd_df) <- colnames(filtered_table2)


ARSACS_raw_ibd_df <- rbind(filtered_table2,ARSACS_raw_ibd_df)
#rm(ARSACS_raw2_ibd_df)

## Calculate generation num
ARSACS_raw_ibd_df$gen <- 3/(2*((filtered_table2$length)/100))

## Make mean
mean_IBD_df <- as.data.frame(xtabs(length~hg1+hg2,aggregate(length~hg1+hg2,ARSACS_raw_ibd_df,mean)))
mean_IBD_df$gen <- 3/(2*((mean_IBD_df$Freq)/100))


myfile <- paste('mean_IBD_chr13_window40_scale10_lod3_length',length,'_sharing_hgs_mean0.8Mb(5haplotypes_with_themselves).jpg',sep='')
jpeg(filename=myfile,width=15,height=11,units="in",res=300)
ggplot(mean_IBD_df, aes(hg1,hg2, fill= Freq)) +
  geom_tile() +
  theme(axis.text.x = element_text( angle = 45, hjust = 1,size = 16 )) +
  theme(axis.text.y = element_text( size = 16 )) +
  # scale_fill_gradientn(colors = c("#006228", "#85BC47",  "#FEFDBE"  ,"#EB9C00", "#A51122")) +
  scale_fill_gradientn(colors = c('#1B1465', '#5DB9FA',"#F0A0FF", '#8759F0','#CCCCFF')) +
  geom_text(aes(label=round(gen)),cex=3,color='white') +
  labs(title = paste('IBD sharing length (',length,'cM) between haplogroups',sep=''),x='',y='',fill='Mean length (cM)')
dev.off()



