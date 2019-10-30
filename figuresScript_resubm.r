library(MASS)
library(dplyr)
library(ggplot2)
library(mltools)
library(reshape2)
library(Hmisc)
library(tidyr)
library(stringr)
library(ggpubr)
library(lme4)
library(lsmeans)
library(contrast)
library(multcomp)
library(BiocManager)
library(HTqPCR)
library(nondetects)


setwd("/media/nvakirlis/hd1/projects/PCaP1_paper/pcap_article_bpi/final_csv/")

confirm_contrasts <- function(factors, plus, minus)
{
  res = c()
  found=0
  for (i in factors)
  {
    if (i %in% plus)
    {
      res <- c(res, 1)
      found=found+1
    }
    else if (i %in% minus)
    {
      res <- c(res, -1)
      found=found+1
    }
    else
    {
      res <- c(res, 0)
    }
  }
  
  return(res)
}

lmer_df_method = "Kenward-Roger"


####################
##### Figure 1 #####
####################

VIGS_d_data <- read.csv("resubmission/data/Figure1_data.csv")

total_contrast_df <- data.frame()

for (gene in c("PCaP", "PVY"))
{
  df <- select(VIGS_d_data,
               Sample,
               Treatment,
               UBI,
               gene) %>%
    rename(Ct_ref="UBI", Ct_target=gene) %>%
    filter(!is.na(Ct_target)) %>%
    gather(gene,
           "Ct",
           "Ct_ref",
           "Ct_target")
  
  
  df$Treatment <- factor(df$Treatment, 
                         levels = unique(df$Treatment),
                         ordered=TRUE)
  
  df$gene <- factor(df$gene,
                    levels = c("Ct_ref", "Ct_target"),
                    ordered=TRUE)
                    
  for (i in unique(df$Treatment))
  {
	var <- filter(df, Treatment == i)
	temp <- var(-(var[which(var$gene=="Ct_ref"), "Ct"] - var[which(var$gene=="Ct_target"), "Ct"]))
	print(c(gene, i, temp))
  }                  
  
  df_lmm <- lmer(Ct ~ Treatment*gene + (1|Sample),
                 data=df)
  
  lsmod <- lsmeans(df_lmm,
                   pairwise ~ Treatment*gene,
                   lmer.df = lmer_df_method,
                   adjust="tukey")
  
  combs <- paste(summary(lsmod$lsmeans)[,1],
                 summary(lsmod$lsmeans)[,2],
                 sep=",")
  
  contrastVIGS <- confirm_contrasts(combs,
                                    c("TRV/PVY,Ct_target", "VIGS,Ct_ref"),
                                    c("TRV/PVY,Ct_ref", "VIGS,Ct_target"))
  contrastH <- confirm_contrasts(combs,
                                 c("TRV/PVY,Ct_target", "H,Ct_ref"),
                                 c("TRV/PVY,Ct_ref", "H,Ct_target"))
  contrastPVY <- confirm_contrasts(combs,
                                   c("TRV/PVY,Ct_target", "PVY,Ct_ref"),
                                   c("TRV/PVY,Ct_ref", "PVY,Ct_target"))
  
  if (gene=="PCaP")
  {
    contrast_list <- list(VIGS = contrastVIGS,
                          H = contrastH,
                          PVY = contrastPVY)
    pcap_fig2_contrasts <- as.data.frame(emmeans::contrast(lsmod[[1]],
                                                           contrast_list,
                                                           by=NULL,
                                                           adjust="tukey"))
    
    pcap_fig2_contrasts$contrast <- factor(pcap_fig2_contrasts$contrast,
                                           levels=c("H", "PVY", "VIGS"),
                                           ordered = TRUE)
  }
  else if (gene=="PVY")
  {
    contrast_list <- list(PVY = contrastPVY,
                          VIGS = contrastVIGS)
    pcap_fig2_contrasts <- as.data.frame(emmeans::contrast(lsmod[[1]],
                                                           contrast_list,
                                                           by=NULL,
                                                           adjust="tukey"))
    
    pcap_fig2_contrasts$contrast <- factor(pcap_fig2_contrasts$contrast,
                                           levels=c("PVY", "VIGS"),
                                           ordered = TRUE)
  }
  
  pcap_fig2_contrasts$Target <- gene
  total_contrast_df <- rbind(total_contrast_df,
                             pcap_fig2_contrasts)
  
}

total_contrast_df <- filter(total_contrast_df, !is.na(p.value))
Fig1_pvals <- filter(total_contrast_df, !is.na(p.value))

sampl_num <- as.data.frame(group_by(VIGS_d_data, Treatment) %>%
				dplyr::summarize(Freq = length(unique(Sample))))
				
Fig1_pvals$no_plants_left <- sampl_num[match(Fig1_pvals$contrast, sampl_num$Treatment), "Freq"]
Fig1_pvals$no_plants_right <- sampl_num[which(sampl_num$Treatment=="TRV/PVY"), "Freq"]

fig1_both <- ggplot(data = total_contrast_df,
                   aes(x=contrast,
                       y=estimate,
                       fill=Target)) +
  geom_bar(stat="identity",
           colour="black") +
  geom_errorbar(aes(ymin=estimate-2*SE,
                    ymax=estimate+2*SE),
                width=0.3) +
  geom_hline(yintercept = 0) +
  facet_grid(~Target,
             scales="free_x",
             space = "free") +
  labs(x = "Treatment",
       y = expression(Log[2]*"FC")) +
  scale_fill_manual(values=c("white", "gray80", "gray38")) +
  theme_bw() +
  theme(text = element_text(size=14),
        axis.title = element_text(size=16)) +
  ggtitle("VIGS both")

ggsave("resubmission/figures/Figure1.pdf")



####################
##### Figure 3 #####
####################

fig2_data <- read.csv("resubmission/data/Figure3_data.csv")


total_contrast_df <- data.frame()

for (gene in c("PcaP", "P3_PVY", "TRV"))
{
  df <- select(fig2_data,
               Sample,
               Treatment,
               TIP41,
               gene) %>%
    rename(Ct_ref="TIP41", Ct_target=gene) %>%
    filter(!is.na(Ct_target)) %>%
    gather(gene,
           "Ct",
           "Ct_ref",
           "Ct_target")
  
  
  df$Treatment <- factor(df$Treatment, 
                         levels = unique(df$Treatment),
                         ordered=TRUE)
  
  df$gene <- factor(df$gene,
                    levels = c("Ct_ref", "Ct_target"),
                    ordered=TRUE)
  
  df_lmm <- lmer(Ct ~ Treatment*gene + (1|Sample),
                 data=df)          
  
  lsmod <- lsmeans(df_lmm,
                   pairwise ~ Treatment*gene,
                   lmer.df = lmer_df_method,
                   adjust="tukey")
  
  combs <- paste(summary(lsmod$lsmeans)[,1],
                 summary(lsmod$lsmeans)[,2],
                 sep=",")
  
  contrastTRV <- confirm_contrasts(combs,
                                   c("TRV_PVY,Ct_target", "TRV,Ct_ref"),
                                   c("TRV_PVY,Ct_ref", "TRV,Ct_target"))
  contrastH <- confirm_contrasts(combs,
                                 c("TRV_PVY,Ct_target", "H,Ct_ref"),
                                 c("TRV_PVY,Ct_ref", "H,Ct_target"))
  contrastPVY <- confirm_contrasts(combs,
                                   c("TRV_PVY,Ct_target", "PVY,Ct_ref"),
                                   c("TRV_PVY,Ct_ref", "PVY,Ct_target"))
  contrastVIGS_PVY <- confirm_contrasts(combs,
                                        c("TRV_PVY,Ct_target", "VIGS_PVY,Ct_ref"),
                                        c("TRV_PVY,Ct_ref", "VIGS_PVY,Ct_target"))
  contrastVIGS <- confirm_contrasts(combs,
                                    c("TRV_PVY,Ct_target", "VIGS,Ct_ref"),
                                    c("TRV_PVY,Ct_ref", "VIGS,Ct_target"))
  
  if (gene=="PcaP")
  {
    contrast_list <- list(TRV = contrastTRV,
                          H = contrastH,
                          PVY = contrastPVY,
                          VIGS_PVY = contrastVIGS_PVY,
                          VIGS = contrastVIGS)
    pcap_fig2_contrasts <- as.data.frame(emmeans::contrast(lsmod[[1]],
                                                           contrast_list,
                                                           by=NULL,
                                                           adjust="tukey"))
    
    pcap_fig2_contrasts$contrast <- factor(pcap_fig2_contrasts$contrast,
                                           levels=c("H", "PVY", "TRV", "VIGS", "VIGS_PVY"),
                                           ordered = TRUE)
                                        
  }
  else if (gene=="P3_PVY")
  {
    contrast_list <- list(PVY = contrastPVY,
                          VIGS_PVY = contrastVIGS_PVY)
    pcap_fig2_contrasts <- as.data.frame(emmeans::contrast(lsmod[[1]],
                                                           contrast_list,
                                                           by=NULL,
                                                           adjust="tukey"))
                                                           
    
    
    pcap_fig2_contrasts$contrast <- factor(pcap_fig2_contrasts$contrast,
                                           levels=c("PVY", "VIGS_PVY"),
                                           ordered = TRUE)
                                           
  }
  else if (gene=="TRV")
  {
    contrast_list <- list(TRV = contrastTRV,
                          VIGS_PVY = contrastVIGS_PVY,
                          VIGS = contrastVIGS)
    pcap_fig2_contrasts <- as.data.frame(emmeans::contrast(lsmod[[1]],
                                                           contrast_list,
                                                           by=NULL,
                                                           adjust="tukey"))
                                                           
                                                          
    
    pcap_fig2_contrasts$contrast <- factor(pcap_fig2_contrasts$contrast,
                                           levels=c("TRV", "VIGS", "VIGS_PVY"),
                                           ordered = TRUE)
  }
  
  pcap_fig2_contrasts$Target <- gene
  total_contrast_df <- rbind(total_contrast_df,
                             pcap_fig2_contrasts)
  
}

total_contrast_df <- filter(total_contrast_df, !is.na(p.value))
Fig3_pvals <- filter(total_contrast_df, !is.na(p.value))

sampl_num <- as.data.frame(group_by(fig2_data, Treatment) %>%
				dplyr::summarize(Freq = length(unique(Sample))))
				
Fig3_pvals$no_plants_left <- sampl_num[match(Fig3_pvals$contrast, sampl_num$Treatment), "Freq"]
Fig3_pvals$no_plants_right <- sampl_num[which(sampl_num$Treatment=="TRV_PVY"), "Freq"]

fig2 <- ggplot(data = total_contrast_df,
               aes(x=contrast,
                   y=estimate,
                   fill=Target)) +
  geom_bar(stat="identity",
           colour="black") +
  geom_errorbar(aes(ymin=estimate-2*SE,
                    ymax=estimate+2*SE),
                width=0.3) +
  geom_hline(yintercept = 0) +
  facet_grid(~Target, scales = "free_x", space="free") +
  labs(x = "Treatment",
       y = expression(Log[2]*"FC")) +
  theme_bw() +
  scale_fill_manual(values=c("white", "gray80", "gray38")) +
  scale_y_continuous(limits=c(-5, 7),
                     breaks = seq(-5, 7, by = 2)) +
  theme(legend.position = "None",
        text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.text.x = element_text(angle=60, hjust=1))

ggsave("resubmission/figures/Figure3.pdf")



####################
##### Figure 4 #####
####################

fig3_data <- read.csv("resubmission/data/Figure4_data.csv")


total_contrast_df <- data.frame()


for (gene in c("PcaP", "P3_PVY", "TRV"))
{
  df <- select(fig3_data,
                Sample,
                Treatment,
                UBI,
                gene) %>%
    rename(Ct_ref="UBI", Ct_target=gene) %>%
    filter(!is.na(Ct_target)) %>%
    gather(gene,
           "Ct",
           "Ct_ref",
           "Ct_target")
  
  
  df$Treatment <- factor(df$Treatment, 
                         levels = unique(df$Treatment),
                         ordered=TRUE)
  
  df$gene <- factor(df$gene,
                    levels = c("Ct_ref", "Ct_target"),
                    ordered=TRUE)
  
  df_lmm <- lmer(Ct ~ Treatment*gene + (1|Sample),
                 data=df)
  
  lsmod <- lsmeans(df_lmm,
                   pairwise ~ Treatment*gene,
                   lmer.df = lmer_df_method,
                   adjust="tukey")
  
  combs <- paste(summary(lsmod$lsmeans)[,1],
                 summary(lsmod$lsmeans)[,2],
                 sep=",")
  
  contrastTRV <- confirm_contrasts(combs,
                                   c("TRV_PVY,Ct_target", "TRV,Ct_ref"),
                                   c("TRV_PVY,Ct_ref", "TRV,Ct_target"))
  contrastH <- confirm_contrasts(combs,
                                 c("TRV_PVY,Ct_target", "H,Ct_ref"),
                                 c("TRV_PVY,Ct_ref", "H,Ct_target"))
  contrastPVY <- confirm_contrasts(combs,
                                   c("TRV_PVY,Ct_target", "PVY,Ct_ref"),
                                   c("TRV_PVY,Ct_ref", "PVY,Ct_target"))
  contrastOEX_PVY <- confirm_contrasts(combs,
                                       c("TRV_PVY,Ct_target", "OEX_PVY,Ct_ref"),
                                       c("TRV_PVY,Ct_ref", "OEX_PVY,Ct_target"))
  contrastOEX <- confirm_contrasts(combs,
                                   c("TRV_PVY,Ct_target", "OEX,Ct_ref"),
                                   c("TRV_PVY,Ct_ref", "OEX,Ct_target"))
  
  if (gene=="PcaP")
  {
    contrast_list <- list(TRV = contrastTRV,
                          H = contrastH,
                          PVY = contrastPVY,
                          OEX_PVY = contrastOEX_PVY,
                          OEX = contrastOEX)
    pcap_fig2_contrasts <- as.data.frame(emmeans::contrast(lsmod[[1]],
                                                           contrast_list,
                                                           by=NULL,
                                                           adjust="tukey"))
    
    pcap_fig2_contrasts$contrast <- factor(pcap_fig2_contrasts$contrast,
                                           levels=c("H", "PVY", "TRV", "OEX", "OEX_PVY"),
                                           ordered = TRUE)
  }
  else if (gene=="P3_PVY")
  {
    contrast_list <- list(PVY = contrastPVY,
                          OEX_PVY = contrastOEX_PVY)
    pcap_fig2_contrasts <- as.data.frame(emmeans::contrast(lsmod[[1]],
                                                           contrast_list,
                                                           by=NULL,
                                                           adjust="tukey"))
    
    pcap_fig2_contrasts$contrast <- factor(pcap_fig2_contrasts$contrast,
                                           levels=c("PVY", "OEX_PVY"),
                                           ordered = TRUE)
  }
  else if (gene=="TRV")
  {
    contrast_list <- list(TRV = contrastTRV,
                          OEX_PVY = contrastOEX_PVY,
                          OEX = contrastOEX)
    pcap_fig2_contrasts <- as.data.frame(emmeans::contrast(lsmod[[1]],
                                                           contrast_list,
                                                           by=NULL,
                                                           adjust="tukey"))
    
    pcap_fig2_contrasts$contrast <- factor(pcap_fig2_contrasts$contrast,
                                           levels=c("TRV", "OEX", "OEX_PVY"),
                                           ordered = TRUE)
  }
  
  pcap_fig2_contrasts$Target <- gene
  total_contrast_df <- rbind(total_contrast_df,
                             pcap_fig2_contrasts)
  
}

total_contrast_df <- filter(total_contrast_df, !is.na(p.value))
Fig4_pvals <- filter(total_contrast_df, !is.na(p.value))

sampl_num <- as.data.frame(group_by(fig3_data, Treatment) %>%
				dplyr::summarize(Freq = length(unique(Sample))))
				
Fig4_pvals$no_plants_left <- sampl_num[match(Fig4_pvals$contrast, sampl_num$Treatment), "Freq"]
Fig4_pvals$no_plants_right <- sampl_num[which(sampl_num$Treatment=="TRV_PVY"), "Freq"]

fig3_ge1 <- ggplot(data = total_contrast_df,
                   aes(x=contrast,
                       y=estimate,
                       fill=Target)) +
  geom_bar(stat="identity",
           colour="black") +
  geom_errorbar(aes(ymin=estimate-2*SE,
                    ymax=estimate+2*SE),
                width=0.3) +
  geom_hline(yintercept = 0) +
  facet_grid(~Target,
             scales="free_x",
             space = "free") +
  labs(x = "Treatment",
       y = expression(Log[2]*"FC")) +
  scale_fill_manual(values=c("white", "gray80", "gray38")) +
  scale_y_continuous(limits=c(-4, 8),
                     breaks = seq(-4, 8, by = 2)) +
  theme_bw() +
  theme(text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "None")

ggsave("resubmission/figures/Figure4.pdf")



#######################
##### Supp. Fig 2 #####
#######################


rv_data <- read.csv("resubmission/data/SupFigure2_data.csv")

total_contrast_df <- data.frame()
total_contrast_df_2 <- data.frame()

for (gene in c("PcaP", "PVY"))
{
  
  df <- select(rv_data,
               Sample,
               Treatment,
               UBI,
               gene) %>%
    dplyr::rename(Ct_ref="UBI", Ct_target=gene) %>%
    filter(!is.na(Ct_target)) %>%
    gather(gene,
           "Ct",
           "Ct_ref",
           "Ct_target")
  
  
  df$Treatment <- factor(df$Treatment, 
                         levels = unique(df$Treatment),
                         ordered=TRUE)
  
  df$gene <- factor(df$gene,
                    levels = c("Ct_ref", "Ct_target"),
                    ordered=TRUE)
  
  df_lmm <- lmer(Ct ~ Treatment*gene + (1|Sample),
                 data=df)
  
  lsmod <- lsmeans(df_lmm,
                   pairwise ~ Treatment*gene,
                   lmer.df = lmer_df_method,
                   adjust="tukey")
  
  combs <- paste(summary(lsmod$lsmeans)[,1],
                 summary(lsmod$lsmeans)[,2],
                 sep=",")
  
  contrastVIGS <- confirm_contrasts(combs,
                                    c("PVY,Ct_target", "VIGS,Ct_ref"),
                                    c("PVY,Ct_ref", "VIGS,Ct_target"))
  contrastVIGS_p19 <- confirm_contrasts(combs,
                                 c("PVY,Ct_target", "VIGS_p19,Ct_ref"),
                                 c("PVY,Ct_ref", "VIGS_p19,Ct_target"))
  contrastPVY_p19 <- confirm_contrasts(combs,
                                        c("PVY,Ct_target", "PVY_p19,Ct_ref"),
                                        c("PVY,Ct_ref", "PVY_p19,Ct_target"))
                                        
  contrastVIGS_p19_VIGS <- confirm_contrasts(combs,
                                 c("VIGS,Ct_target", "VIGS_p19,Ct_ref"),
                                 c("VIGS,Ct_ref", "VIGS_p19,Ct_target"))										 
  
  
  contrast_list <- list(VIGS = contrastVIGS,
                        VIGS_p19 = contrastVIGS_p19,
                        PVY_p19 = contrastPVY_p19)
                        
  contrast_list_2 <- list(VIGSp19_VIGS = contrastVIGS_p19_VIGS)                        
  
  pcap_fig2_contrasts <- as.data.frame(emmeans::contrast(lsmod[[1]],
                                                         contrast_list,
                                                         by=NULL,
                                                         adjust="tukey"))
                                                         
  pcap_fig2_contrasts_2 <- as.data.frame(emmeans::contrast(lsmod[[1]],
                                                           contrast_list_2,
                                                           by=NULL,
                                                           adjust="tukey"))                                                         
  
  pcap_fig2_contrasts$contrast <- factor(pcap_fig2_contrasts$contrast,
                                         levels=c("PVY_p19", "VIGS", "VIGS_p19"),
                                         ordered = TRUE)
  
  pcap_fig2_contrasts$Target <- gene
  
  pcap_fig2_contrasts_2$Target <- gene
  
  total_contrast_df <- rbind(total_contrast_df,
                             pcap_fig2_contrasts)
                             
  total_contrast_df_2 <- rbind(total_contrast_df_2,
                               pcap_fig2_contrasts_2)                             
  
}

#> total_contrast_df_2 this is a contrast not shown in the figure but commented in the article
#VIGSp19_VIGS	3.152974	0.4330127	4	7.281481	0.001890368	PcaP	2	2	Supp. Fig. 2
#VIGSp19_VIGS	3.671151	0.503736	4	7.287847	0.001884159	PVY	2	2	Supp. Fig. 2

add_sfig2_pvals <- total_contrast_df_2
add_sfig2_pvals$no_plants_left <- 2
add_sfig2_pvals$no_plants_right <- 2

total_contrast_df <- filter(total_contrast_df, !is.na(p.value))
total_contrast_df$Target <- as.character(total_contrast_df$Target)
sfig2_pvals <- filter(total_contrast_df, !is.na(p.value))

sampl_num <- as.data.frame(group_by(rv_data, Treatment) %>%
				dplyr::summarize(Freq = length(unique(Sample))))
				
sfig2_pvals$no_plants_left <- sampl_num[match(sfig2_pvals$contrast, sampl_num$Treatment), "Freq"]
sfig2_pvals$no_plants_right <- sampl_num[which(sampl_num$Treatment=="PVY"), "Freq"]

reverse_VIGS_figure <- ggplot(data = total_contrast_df,
                            aes(x=contrast,
                                y=estimate,
                                fill=Target)) +
  geom_bar(stat="identity",
           colour="black",
           position = "dodge") +
  geom_errorbar(aes(ymin=estimate-2*SE,
                    ymax=estimate+2*SE),
                width=0.3,
                position = position_dodge(width=0.9)) +
  geom_hline(yintercept = 0) +
  labs(x = "Treatment",
       y = expression(Log[2]*"FC")) +
  scale_fill_manual(values=c("white", "gray80")) +
  scale_y_continuous(limits=c(-8, 2),
                     breaks = seq(-8, 2, by = 2)) +
  theme_bw() +
  theme(text = element_text(size=14),
        axis.title = element_text(size=16))

ggsave("resubmission/figures/SupFigure2.pdf")



###############################
######## Sup. Figure 3 ########
###############################

fig2_data <- read.csv("resubmission/data/SupFigure3_data.csv")

total_contrast_df <- data.frame()

for (gene in c("PcaP", "P3_PVY", "TRV"))
{
  g <- gene  	
  fig2_data$gene <- as.character(fig2_data$gene)
  df <- filter(fig2_data,
          !is.na(Ct),
          (gene==g | gene=="UBI"))
          
   		
  df[which(df$gene!="UBI"), "gene"] = "Ct_target"
  df[which(df$gene=="UBI"), "gene"] = "Ct_ref"
  
  df$Treatment <- factor(df$Treatment, 
                         levels = unique(df$Treatment),
                         ordered=TRUE)
  
  df$gene <- factor(df$gene,
                    levels = c("Ct_ref", "Ct_target"),
                    ordered=TRUE)
                    

  
  df_lmm <- lmer(Ct ~ Treatment*gene + (1|Sample),
                 data=df)          
  
  lsmod <- lsmeans(df_lmm,
                   pairwise ~ Treatment*gene,
                   lmer.df = lmer_df_method,
                   adjust="tukey")
  
  combs <- paste(summary(lsmod$lsmeans)[,1],
                 summary(lsmod$lsmeans)[,2],
                 sep=",")
  
  contrastTRV <- confirm_contrasts(combs,
                                   c("TRV_PVY,Ct_target", "TRV,Ct_ref"),
                                   c("TRV_PVY,Ct_ref", "TRV,Ct_target"))
  contrastH <- confirm_contrasts(combs,
                                 c("TRV_PVY,Ct_target", "H,Ct_ref"),
                                 c("TRV_PVY,Ct_ref", "H,Ct_target"))
  contrastPVY <- confirm_contrasts(combs,
                                   c("TRV_PVY,Ct_target", "PVY,Ct_ref"),
                                   c("TRV_PVY,Ct_ref", "PVY,Ct_target"))
  contrastVIGS_PVY <- confirm_contrasts(combs,
                                        c("TRV_PVY,Ct_target", "VIGS_PVY,Ct_ref"),
                                        c("TRV_PVY,Ct_ref", "VIGS_PVY,Ct_target"))
  contrastVIGS <- confirm_contrasts(combs,
                                    c("TRV_PVY,Ct_target", "VIGS,Ct_ref"),
                                    c("TRV_PVY,Ct_ref", "VIGS,Ct_target"))
  
  if (gene=="PcaP")
  {
    contrast_list <- list(TRV = contrastTRV,
                          H = contrastH,
                          PVY = contrastPVY,
                          VIGS_PVY = contrastVIGS_PVY,
                          VIGS = contrastVIGS)
    pcap_fig2_contrasts <- as.data.frame(emmeans::contrast(lsmod[[1]],
                                                           contrast_list,
                                                           by=NULL,
                                                           adjust="tukey"))
    
    pcap_fig2_contrasts$contrast <- factor(pcap_fig2_contrasts$contrast,
                                           levels=c("H", "PVY", "TRV", "VIGS", "VIGS_PVY"),
                                           ordered = TRUE)
  }
  else if (gene=="P3_PVY")
  {
    contrast_list <- list(PVY = contrastPVY,
                          VIGS_PVY = contrastVIGS_PVY)
    pcap_fig2_contrasts <- as.data.frame(emmeans::contrast(lsmod[[1]],
                                                           contrast_list,
                                                           by=NULL,
                                                           adjust="tukey"))
    
    pcap_fig2_contrasts$contrast <- factor(pcap_fig2_contrasts$contrast,
                                           levels=c("PVY", "VIGS_PVY"),
                                           ordered = TRUE)
  }
  else if (gene=="TRV")
  {
    contrast_list <- list(TRV = contrastTRV,
                          VIGS_PVY = contrastVIGS_PVY,
                          VIGS = contrastVIGS)
    pcap_fig2_contrasts <- as.data.frame(emmeans::contrast(lsmod[[1]],
                                                           contrast_list,
                                                           by=NULL,
                                                           adjust="tukey"))
    
    pcap_fig2_contrasts$contrast <- factor(pcap_fig2_contrasts$contrast,
                                           levels=c("TRV", "VIGS", "VIGS_PVY"),
                                           ordered = TRUE)
  }
  
  pcap_fig2_contrasts$Target <- gene
  total_contrast_df <- rbind(total_contrast_df,
                             pcap_fig2_contrasts)
  
}

total_contrast_df <- filter(total_contrast_df, !is.na(p.value))
sfig3_pvals <- filter(total_contrast_df, !is.na(p.value))

sampl_num <- as.data.frame(group_by(fig2_data, Treatment) %>%
				dplyr::summarize(Freq = length(unique(Sample))))
				
sfig3_pvals$no_plants_left <- sampl_num[match(sfig3_pvals$contrast, sampl_num$Treatment), "Freq"]
sfig3_pvals$no_plants_right <- sampl_num[which(sampl_num$Treatment=="TRV_PVY"), "Freq"]

fig2 <- ggplot(data = total_contrast_df,
               aes(x=contrast,
                   y=estimate,
                   fill=Target)) +
  geom_bar(stat="identity",
           colour="black") +
  geom_errorbar(aes(ymin=estimate-2*SE,
                    ymax=estimate+2*SE),
                width=0.3) +
  geom_hline(yintercept = 0) +
  facet_grid(~Target, scales = "free_x", space="free") +
  labs(x = "Treatment",
       y = expression(Log[2]*"FC")) +
  theme_bw() +
  scale_fill_manual(values=c("white", "gray80", "gray38")) +
  theme(legend.position = "None",
        text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.text.x = element_text(angle=60, hjust=1)) +
  scale_y_continuous(limits=c(-5, 7),
                     breaks = seq(-5, 7, by = 2))

ggsave("resubmission/figures/SupFigure3.pdf")


#########################
##### Sup. Figure 5 #####
#########################


VIGS_rep_data <- read.csv("resubmission/data/SupFigure5_data.csv")

total_contrast_df <- data.frame()

for (gene in c("PcaP","TRV"))
{
  df <- select(VIGS_rep_data,
               Sample,
               Treatment,
               UBI,
               gene) %>%
    rename(Ct_ref="UBI", Ct_target=gene) %>%
    filter(!is.na(Ct_target)) %>%
    gather(gene,
           "Ct",
           "Ct_ref",
           "Ct_target")
  
  
  df$Treatment <- factor(df$Treatment, 
                         levels = unique(df$Treatment),
                         ordered=TRUE)
  
  df$gene <- factor(df$gene,
                    levels = c("Ct_ref", "Ct_target"),
                    ordered=TRUE)
  
  df_lmm <- lmer(Ct ~ Treatment*gene + (1|Sample),
                 data=df)
  
  lsmod <- lsmeans(df_lmm,
                   pairwise ~ Treatment*gene,
                   lmer.df = lmer_df_method,
                   adjust="tukey")
  
  combs <- paste(summary(lsmod$lsmeans)[,1],
                 summary(lsmod$lsmeans)[,2],
                 sep=",")
  
                                   
    contrastVIGS <- confirm_contrasts(combs,
                                    c("TRV/PVY,Ct_target", "VIGS,Ct_ref"),
                                    c("TRV/PVY,Ct_ref", "VIGS,Ct_target"))
  contrastTRV <- confirm_contrasts(combs,
                                    c("TRV/PVY,Ct_target", "TRV,Ct_ref"),
                                    c("TRV/PVY,Ct_ref", "TRV,Ct_target"))                                     
  contrastVIGS_PVY <- confirm_contrasts(combs,
                                        c("TRV/PVY,Ct_target", "VIGS/PVY,Ct_ref"),
                                        c("TRV/PVY,Ct_ref", "VIGS/PVY,Ct_target"))
  contrastH <- confirm_contrasts(combs,
                                 c("TRV/PVY,Ct_target", "H,Ct_ref"),
                                 c("TRV/PVY,Ct_ref", "H,Ct_target"))
  contrastPVY <- confirm_contrasts(combs,
                                   c("TRV/PVY,Ct_target", "PVY,Ct_ref"),
                                   c("TRV/PVY,Ct_ref", "PVY,Ct_target"))                                   
  
  if (gene=="PcaP")
  {
    contrast_list <- list(VIGS = contrastVIGS,
                          VIGS_PVY = contrastVIGS_PVY,
                          H = contrastH,
                          PVY = contrastPVY,
                          TRV = contrastTRV)
   
    pcap_fig2_contrasts <- as.data.frame(emmeans::contrast(lsmod[[1]],
                                                           contrast_list,
                                                           by=NULL,
                                                           adjust="tukey"))
    
    pcap_fig2_contrasts$contrast <- factor(pcap_fig2_contrasts$contrast,
                                           levels=c("H", "PVY", "TRV", "VIGS", "VIGS_PVY"),
                                           ordered = TRUE)
  }
  else if (gene=="TRV")
  {
    contrast_list <- list(TRV = contrastTRV,
      VIGS = contrastVIGS,
      VIGS_PVY = contrastVIGS_PVY)
      
    pcap_fig2_contrasts <- as.data.frame(emmeans::contrast(lsmod[[1]],
                                                           contrast_list,
                                                           by=NULL,
                                                           adjust="tukey"))
    
    pcap_fig2_contrasts$contrast <- factor(pcap_fig2_contrasts$contrast,
                                           levels=c("TRV", "VIGS", "VIGS_PVY"),
                                           ordered = TRUE)
  }
  
  
  pcap_fig2_contrasts$Target <- gene
  total_contrast_df <- rbind(total_contrast_df,
                             pcap_fig2_contrasts)
  
}

total_contrast_df <- filter(total_contrast_df, !is.na(p.value))
sfig5_pvals <- filter(total_contrast_df, !is.na(p.value))

sampl_num <- as.data.frame(group_by(VIGS_rep_data, Treatment) %>%
				dplyr::summarize(Freq = length(unique(Sample))))
				
sfig5_pvals$no_plants_left <- sampl_num[match(sfig5_pvals$contrast, sampl_num$Treatment), "Freq"]
sfig5_pvals$no_plants_right <- sampl_num[which(sampl_num$Treatment=="TRV"), "Freq"]

#add no_plants_left in two rows where contrast has / instead of _
sfig5_pvals[2,"no_plants_left"] <- 3
sfig5_pvals[8,"no_plants_left"] <- 3

vigs_rep_fig <- ggplot(data = total_contrast_df,
                       aes(x=contrast,
                           y=estimate,
                           fill=Target)) +
  geom_bar(stat="identity",
           colour="black") +
  geom_errorbar(aes(ymin=estimate-2*SE,
                    ymax=estimate+2*SE),
                width=0.3) +
  geom_hline(yintercept = 0) +
  facet_grid(~Target,
             scales="free_x",
             space = "free") +
  labs(x = "Treatment",
       y = expression(Log[2]*"FC")) +
  scale_fill_manual(values=c("white","gray80")) +
  theme_bw() +
  theme(text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "None")
  
ggsave("resubmission/figures/SupFigure5.pdf")


#######################
##### Sup. Fig 1  #####
#######################

tc_data <- read.csv("resubmission/data/SupFigure1_data.csv")

total_contrast_df <- data.frame()

for (tm in c(1, 12, 19))
{
  df <-filter(tc_data,
              timepoint==tm) %>% 
    select(Sample,
           Treatment,
           UBI,
           PcaP) %>%
    rename(Ct_ref="UBI", Ct_target="PcaP") %>%
    filter(!is.na(Ct_target)) %>%
    gather(gene,
           "Ct",
           "Ct_ref",
           "Ct_target")
  
  
  df$Treatment <- factor(df$Treatment, 
                         levels = unique(df$Treatment),
                         ordered=TRUE)
  
  df$gene <- factor(df$gene,
                    levels = c("Ct_ref", "Ct_target"),
                    ordered=TRUE)
  
  df_lmm <- lmer(Ct ~ Treatment*gene + (1|Sample),
                 data=df)
  
  lsmod <- lsmeans(df_lmm,
                   pairwise ~ Treatment*gene,
                   lmer.df = lmer_df_method,
                   adjust="tukey")
  
  combs <- paste(summary(lsmod$lsmeans)[,1],
                 summary(lsmod$lsmeans)[,2],
                 sep=",")
  
  contrastVIGS <- confirm_contrasts(combs,
                                    c("TRV,Ct_target", "VIGS,Ct_ref"),
                                    c("TRV,Ct_ref", "VIGS,Ct_target"))
  contrastH <- confirm_contrasts(combs,
                                 c("TRV,Ct_target", "H,Ct_ref"),
                                 c("TRV,Ct_ref", "H,Ct_target"))
  
  
  contrast_list <- list(VIGS = contrastVIGS,
                        H = contrastH)
  
  pcap_fig2_contrasts <- as.data.frame(emmeans::contrast(lsmod[[1]],
                                                         contrast_list,
                                                         by=NULL,
                                                         adjust="tukey"))
  
  pcap_fig2_contrasts$contrast <- factor(pcap_fig2_contrasts$contrast,
                                         levels=c("H", "VIGS"),
                                         ordered = TRUE)
  
  pcap_fig2_contrasts$Target <- tm
  total_contrast_df <- rbind(total_contrast_df,
                             pcap_fig2_contrasts)                             
  
}

total_contrast_df <- filter(total_contrast_df, !is.na(p.value))
total_contrast_df$Target <- as.character(total_contrast_df$Target)

sfig1_pvals <- total_contrast_df

sampl_num <- as.data.frame(group_by(tc_data, Treatment, timepoint) %>%
				dplyr::summarize(Freq = length(unique(Sample))))
				
sfig1_pvals$no_plants_left <- sampl_num[match(paste(sfig1_pvals$contrast, sfig1_pvals$Target), paste(sampl_num$Treatment, sampl_num$timepoint)), "Freq"]
sfig1_pvals$no_plants_right <- sampl_num[which(sampl_num$Treatment=="TRV"), "Freq"]

timecourse_figure <- ggplot(data = total_contrast_df,
                            aes(x=contrast,
                                y=estimate,
                                fill=Target)) +
  geom_bar(stat="identity",
           colour="black") +
  geom_errorbar(aes(ymin=estimate-2*SE,
                    ymax=estimate+2*SE),
                width=0.3) +
  geom_hline(yintercept = 0) +
  facet_grid(~Target,
             scales="free_x",
             space = "free") +
  labs(x = "Treatment",
       y = expression(Log[2]*"FC")) +
  scale_fill_manual(values=c("white", "gray80", "gray38")) +
  theme_bw() +
  theme(text = element_text(size=14),
        axis.title = element_text(size=16))

ggsave("resubmission/figures/SupFigure1.pdf")

###### pvals ########

Fig1_pvals$dataset <- "Figure 1"
Fig3_pvals$dataset <- "Figure 3"
Fig4_pvals$dataset <- "Figure 4"
sfig1_pvals$dataset <- "Supp. Fig. 1"
sfig2_pvals$dataset <- "Supp. Fig. 2"
add_sfig2_pvals$dataset <- "Supp. Fig. 2"
sfig3_pvals$dataset <- "Supp. Fig. 3"
sfig5_pvals$dataset <- "Supp. Fig. 5"

all_pvals <- rbind(Fig1_pvals,
                   Fig3_pvals,
                   Fig4_pvals,
                   sfig1_pvals,
                   sfig2_pvals,
                   add_sfig2_pvals,
                   sfig3_pvals,
                   sfig5_pvals) %>%
			 mutate(estimate=round(estimate, 3),
				    SE=round(SE, 3),
				    df=round(df, 3),
				    t.ratio=round(t.ratio, 3),
				    p.value=signif(p.value, digits=3))
                   

write.csv(all_pvals, 
            "resubmission/Sup_Table_3.csv",
          row.names = FALSE)


#######################
##### Sup. Fig 4  #####
#######################

df <- read.csv("resubmission/data/SupFigure4_data.csv")

df_s <- gather(df, 
			type,
		    measurement,
           VIGS_PVY_GFP,
           VIGS_PVY_RFP,
           TRV_PVY_GFP,
           TRV_PVY_RFP,
           PVY_GFP,
           PVY_RFP)
           
wilcox.test(df$VIGS_PVY_GFP, df$TRV_PVY_GFP)
t.test(df$VIGS_PVY_GFP, df$TRV_PVY_GFP, paired=FALSE)
ks.test(df$VIGS_PVY_GFP, df$TRV_PVY_GFP)
           

df_s$type <- factor(df_s$type, levels = c("VIGS_PVY_GFP", "TRV_PVY_GFP", "PVY_GFP", "VIGS_PVY_RFP", "TRV_PVY_RFP", "PVY_RFP"), ordered=TRUE)

pl <- ggplot() +
		geom_boxplot(data = df_s, aes(y=measurement, x=type, fill=type, colour=type), width=0.8, alpha=0.3, outlier.shape=NA) +
		geom_jitter(data = df_s, aes(y=measurement, x=type, colour=type), width=0.3) +
		scale_fill_manual(values=c("darkgreen","darkgreen","darkgreen","red","red","red")) +
		scale_colour_manual(values=c("darkgreen","darkgreen","darkgreen","red","red","red")) +
		theme_bw() + 
		theme(legend.position="None") +
		xlab("")
		
ggsave("resubmission/figures/SupFigure4.pdf")		
