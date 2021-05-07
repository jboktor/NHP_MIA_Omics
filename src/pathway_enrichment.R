
source("src/_load_packages.R")
source("src/_misc_functions.R")

d1 <- read_xlsx('input_files/Pathway Enrichment.xlsx', sheet = 'TGC')
d2 <- read_xlsx('input_files/Pathway Enrichment.xlsx', sheet = 'MC')
d3 <- read_xlsx('input_files/Pathway Enrichment.xlsx', sheet = 'FC')

Order <- c("CSF", "Plasma", "Jejunum", "Ileum", "Colon", "Feces")
Sample <- factor(d1$Origin, level = Order)

################ COMBINED  MALE & FEMALE SUPERPATHWAY ANALYSIS ################ 
#Pathway on X
p1 <- ggplot(data=d1) + 
  geom_col(mapping=aes(x=Superpathway, y=P_S, fill=Origin)) + 
  theme_classic() + 
  labs(x="Superpathways", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", 
       title = "Rhesus Monkey Metabolic Superpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "RdYlBu")
p1 

#Tissue on X
p2 <- ggplot(data=d1) + 
  geom_col(mapping=aes(x=factor(Origin, level = Order), y=P_S, fill=Superpathway), width = 0.7) + 
  theme_classic() + 
  labs(x="Tissue Origin", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", title = "Superpathway Enrichment in Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  coord_flip() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu")
p2 
ggsave(p2, filename = "figures/pathway_enrichment/CombinedSuperPathwayByTissue.png", 
       height = 4 , width = 10)


######## MALES  ########

#Pathway on X
Mp1 <- ggplot(data=d2) + 
  geom_col(mapping=aes(x=Superpathway, y=P_S, fill=factor(Origin, level = Order)), width = 0.7) + 
  theme_classic() + 
  labs(x="Superpathways", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", 
       title = "Male Rhesus Monkey Metabolic Superpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme(axis.title =element_text(vjust = 0) ) +
  theme(legend.title = element_blank()) + 
  coord_flip() +
  scale_fill_brewer(palette = "RdYlBu") +
  theme(axis.title.y =element_blank())
Mp1 
ggsave("figures/pathway_enrichment/males_SuperPathwayByPATH.png",  height = 4 , width =10)


#Tissue on X
Mp2 <- ggplot(data=d2) + 
  geom_col(mapping=aes(x=factor(Origin, level = Order), y=P_S, fill=Superpathway), width = 0.7) + 
  theme_classic() + 
  labs(x="Tissue Origin", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", title = "Superpathway Enrichment in Male Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  coord_flip() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu")
Mp2 
ggsave("figures/pathway_enrichment/males_SuperPathwayByTissue.png",  height = 4 , width =10)


#### Individual Data Points (with n) Tissue on X
co2 <- ggplot(data=d2) + 
  geom_count(mapping=aes(x=factor(Origin, level = Order), y=P_S, color=Superpathway)) + 
  theme_classic() + 
  labs(x="Tissue Origin", y="Pathway Enrichment Score (P/S)", title = "Male Rhesus Monkey Metabolic Superpathway Enrichment by Tissue") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3")
co2


######## FEMALES  ########

#Pathway on X
Fp1 <- ggplot(data=d3) + 
  geom_col(mapping=aes(x=Superpathway, y=P_S, fill=factor(Origin, level = Order)), width = 0.7) + 
  theme_classic() + 
  labs(x="Superpathways", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", 
       y = NULL,
       title = "Female Rhesus Monkey Metabolic Superpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95), 
        legend.title = "Origin") +
  scale_fill_brewer(palette = "RdYlBu") + 
  coord_flip()
Fp1 
ggsave("figures/pathway_enrichment/females_SuperPathwayByPATH.png",  height = 4 , width =10)

#Tissue on X
Fp2 <- ggplot(data=d3) + 
  geom_col(mapping=aes(x=factor(Origin, level = Order), y=P_S, fill=Superpathway), width = 0.7) + 
  theme_classic() + 
  labs(x="Tissue Origin", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", title = "Superpathway Enrichment in Female Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  coord_flip() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu")
Fp2
ggsave("figures/pathway_enrichment/females_SuperPathwayByTissue.png",  height = 4 , width =10)


#### Individual Data Points (with n) Tissue on X
co3 <- ggplot(data=d3) + 
  geom_count(mapping=aes(x=factor(Origin, level = Order), y=P_S, color=Superpathway)) + 
  theme_classic() + 
  labs(x="Tissue Origin", y="Pathway Enrichment Score (P/S)", title = "Female Rhesus Monkey Metabolic Superpathway Enrichment by Tissue") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3")
co3


###################################### Subpathway Analysis ###################################### 

######################### AMINO ACIDS #########################

AminoAllonly <- filter(d1, Superpathway=="Amino Acid")
AminoAllonly
 ######## MALES & females_  ########
# X is subpathway
S1 <- ggplot(data=AminoAllonly) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=factor(Origin, level = Order)), width = 0.7) + 
  theme_classic() + 
  labs(x="Superpathways", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", title = "Rhesus Monkey Metabolic Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3")
S1 


# X is Tissue
S2 <- ggplot(data=AminoAllonly) + 
  geom_col(mapping=aes(x=factor(Origin, level = Order), y=P_S, fill=Subpathway), width = 0.7) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", title = "Amino Acid Subpathway Enrichment in Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  coord_flip() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") 
S2 
ggsave("figures/pathway_enrichment/CombinedAminoAcidSubPathwayByTissue.png",  height = 4 , width = 10)

#### Individual Data Points (with n) Tissue on X
co4 <- ggplot(data=AminoAllonly) + 
  geom_point(mapping=aes(x=factor(Origin, level = Order), y=P_S, color=Subpathway)) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Individual Pathway Enrichment Score (P/S)", title = "Rhesus Monkey Metabolic Amino Acid Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3") 
co4

######## MALES  ########
AminoMaleonly <- filter(d2, Superpathway=="Amino Acid")

# X is Tissue
colourCount = length(unique(AminoMaleonly$Subpathway))
getPalette = colorRampPalette(brewer.pal(12, "RdYlBu"))

S3 <- ggplot(data=AminoMaleonly) + 
  geom_col(mapping=aes(x=factor(Origin, level = Order), y=P_S, fill=Subpathway), width = 0.7) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", 
       title = "Amino Acid Subpathway Enrichment in Male Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme(axis.title.y =element_blank()) +
  coord_flip() +
  scale_fill_manual(values = getPalette(colourCount)) 
S3 
ggsave("figures/pathway_enrichment/males_AminoAcidSubPathwayByTissue.png", height = 4 , width = 10)

# X is Subpathway
ggplot(data=AminoMaleonly) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=factor(Origin, level = Order)), width = 0.7) + 
  theme_classic() + 
  labs(y="Summated Pathway Enrichment \nPoly(I:C) / Saline", 
       title = "Male Rhesus Monkey Metabolic Amino Acid Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  coord_flip() +
  scale_fill_brewer(palette = "RdYlBu") +
  theme(axis.title.y =element_blank())
ggsave("figures/pathway_enrichment/males_AminoAcidSubPathwayByTissue.png", height = 4 , width = 10)



#### Individual Data Points (with n) Tissue on X
co5 <- ggplot(data=AminoMaleonly) + 
  geom_point(mapping=aes(x=factor(Origin, level = Order), y=P_S, color=Subpathway)) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Individual Pathway Enrichment Score (P/S)", title = "Male Rhesus Monkey Metabolic Amino Acid Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3") 
co5

######## FEMALES  ########

AminoAFemaleOnly <- filter(d3, Superpathway=="Amino Acid")

# X is Tissue
colourCount = length(unique(AminoAFemaleOnly$Subpathway))
getPalette = colorRampPalette(brewer.pal(12, "RdYlBu"))

S4 <- ggplot(data=AminoAFemaleOnly) + 
  geom_col(mapping=aes(x=factor(Origin, level = Order), y=P_S, fill=Subpathway), width = 0.7) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", title = "Amino Acid Subpathway Enrichment in Female Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme(axis.title.y =element_blank()) +
  coord_flip() +
  scale_fill_manual(values = getPalette(colourCount)) 
S4 
ggsave("figures/pathway_enrichment/females_AminoAcidSubPathwayByTissue.png", height = 4 , width = 10)



ggplot(data=AminoAFemaleOnly) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=factor(Origin, level = Order)), width = 0.7) + 
  theme_classic() + 
  labs(y="Summated Pathway Enrichment \nPoly(I:C) / Saline", title = "Female Rhesus Monkey Metabolic Amino Acid Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  coord_flip() +
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") 

ggsave("figures/pathway_enrichment/females_SubPathwayByPATH.png",  height = 4 , width =10)







#### Individual Data Points (with n) Tissue on X
co6 <- ggplot(data=AminoAFemaleOnly) + 
  geom_point(mapping=aes(x=factor(Origin, level = Order), y=P_S, color=Subpathway)) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Individual Pathway Enrichment Score (P/S)", title = "Female Rhesus Monkey Metabolic Amino Acid Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3") 
co6


######################### LIPID ######################### 

######## MALES & FEMALES  ########

LipidAllonly <- filter(d1, Superpathway=="Lipid")


# X is Tissue
l1 <- ggplot(data=LipidAllonly) + 
  geom_col(mapping=aes(x=factor(Origin, level = Order), y=P_S, fill=Subpathway), width = 0.7) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", title = "Lipid Subpathway Enrichment in Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  coord_flip() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") 
l1 
ggsave("figures/pathway_enrichment/CombinedLipidSubPathwayByTissue.png",  height = 4 , width = 10)


#### Individual Data Points (with n) Tissue on X
co7 <- ggplot(data=LipidAllonly) + 
  geom_point(mapping=aes(x=factor(Origin, level = Order), y=P_S, color=Subpathway)) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Individual Pathway Enrichment Score (P/S)", title = "Rhesus Monkey Metabolic Amino Acid Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3") 
co7


######## MALES  ########

LipidonlyMale <- filter(d2, Superpathway=="Lipid")

# Tissue on X
colourCount = length(unique(AminoMaleonly$Subpathway))
getPalette = colorRampPalette(brewer.pal(12, "RdYlBu"))

l2 <- ggplot(data=LipidonlyMale) + 
  geom_col(mapping=aes(x=factor(Origin, level = Order), y=P_S, fill=Subpathway), width = 0.7) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", 
       title = "Lipid Subpathway Enrichment in Male Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme(axis.title.y =element_blank()) +
  coord_flip() +
  scale_fill_manual(values = getPalette(colourCount)) 
l2 
ggsave("figures/pathway_enrichment/males_LipidSubPathwayByTissue.png", height = 4 , width = 10)


#### Individual Data Points (with n) Tissue on X
co8 <- ggplot(data=LipidonlyMale) + 
  geom_point(mapping=aes(x=factor(Origin, level = Order), y=P_S, color=Subpathway)) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Individual Pathway Enrichment Score (P/S)", title = "Male Rhesus Monkey Metabolic Lipid Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3") 
co8

######## FEMALES ########

LipidonlyFemale <- filter(d3, Superpathway=="Lipid")

# X is Tissue
colourCount = length(unique(LipidonlyFemale$Subpathway))
getPalette = colorRampPalette(brewer.pal(15, "RdYlBu"))

l3 <- ggplot(data=LipidonlyFemale) + 
  geom_col(mapping=aes(x=factor(Origin, level = Order), y=P_S, fill=Subpathway), width = 0.7) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", title = "Lipid Subpathway Enrichment in Female Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme(axis.title.y =element_blank()) +
  coord_flip() +
  scale_fill_manual(values = getPalette(colourCount)) 
l3 
ggsave("figures/pathway_enrichment/females_LipidSubPathwayByTissue.png", height = 4 , width = 10)

#### Individual Data Points (with n) Tissue on X
co9 <- ggplot(data=LipidonlyFemale) + 
  geom_point(mapping=aes(x=factor(Origin, level = Order), y=P_S, color=Subpathway)) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Individual Pathway Enrichment Score (P/S)", title = "Female Rhesus Monkey Metabolic Lipid Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3") 
co9

######################### CoFactors #########################

######## MALES & FEMALES  ########

COFAllonly <- filter(d1, Superpathway=="Cofactors and Vitamins")

cf1 <- ggplot(data=COFAllonly) + 
  geom_col(mapping=aes(x=factor(Origin, level = Order), y=P_S, fill=Subpathway), width = 0.7) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", title = "Rhesus Monkey Metabolic Cofactors and Vitamins Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3")
cf1 

#### Individual Data Points (with n) Tissue on X
co10 <- ggplot(data=COFAllonly) + 
  geom_point(mapping=aes(x=factor(Origin, level = Order), y=P_S, color=Subpathway)) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Individual Pathway Enrichment Score (P/S)", title = "Rhesus Monkey Metabolic Cofactors and Vitamins Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3") 
co10

######## MALES  ########

COFAllonlyMale <- filter(d2, Superpathway=="Cofactors and Vitamins")

cf2 <- ggplot(data=COFAllonlyMale) + 
  geom_col(mapping=aes(x=factor(Origin, level = Order), y=P_S, fill=Subpathway), width = 0.7) + 
  theme_classic() + 
  labs(x="Sample Origine", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", title = "Male Rhesus Monkey Metabolic Cofactors and Vitamins Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3")
cf2

#### Individual Data Points (with n) Tissue on X
co11 <- ggplot(data=COFAllonlyMale) + 
  geom_point(mapping=aes(x=factor(Origin, level = Order), y=P_S, color=Subpathway)) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Individual Pathway Enrichment Score (P/S)", title = "Male Rhesus Monkey Metabolic Cofactors and Vitamins Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3") 
co11

######## FEMALES ########

COFAllonlyFemale <- filter(d3, Superpathway=="Cofactors and Vitamins")

cf3 <- ggplot(data=COFAllonlyFemale) + 
  geom_col(mapping=aes(x=factor(Origin, level = Order), y=P_S, fill=Subpathway), width = 0.7) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", title = "Female Rhesus Monkey Metabolic Cofactors and Vitamins Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3")
cf3 

#### Individual Data Points (with n) Tissue on X
co12 <- ggplot(data=COFAllonlyFemale) + 
  geom_point(mapping=aes(x=factor(Origin, level = Order), y=P_S, color=Subpathway)) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Individual Pathway Enrichment Score (P/S)", title = "Female Rhesus Monkey Metabolic Cofactors and Vitamins Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3") 
co12


######################### Carbohydrates  #########################

######## MALES & FEMALES  ########

CarbAllonly <- filter(d1, Superpathway=="Carbohydrate")


# X is Tissue
cb1 <- ggplot(data=CarbAllonly) + 
  geom_col(mapping=aes(x=factor(Origin, level = Order), y=P_S, fill=Subpathway), width = 0.7) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", title = "Carbohydrate Subpathway Enrichment in Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  coord_flip() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") 
cb1 
ggsave("figures/pathway_enrichment/CombinedCarbohydratesSubPathwayByTissue.png",  height = 4 , width = 10)



#### Individual Data Points (with n) Tissue on X
co13 <- ggplot(data=CarbAllonly) + 
  geom_point(mapping=aes(x=factor(Origin, level = Order), y=P_S, color=Subpathway)) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Individual Pathway Enrichment Score (P/S)", title = "Rhesus Monkey Carbohydrate Metabolic Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3") 
co13

######## MALES  ########

# X is Tissue
CarbOnlyMale <- filter(d2, Superpathway=="Carbohydrate")

cb2 <- ggplot(data=CarbOnlyMale) + 
  geom_col(mapping=aes(x=factor(Origin, level = Order), y=P_S, fill=Subpathway), width = 0.7) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", title = "Carbohydrate Metabolism Subpathway Enrichment in Male Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "RdYlBu")+
  theme(axis.title.y =element_blank()) +
  coord_flip() 
cb2
ggsave("figures/pathway_enrichment/males_CarbohydrateSubPathwayByTissue.png", height = 4 , width = 10)


#### Individual Data Points (with n) Tissue on X
co14 <- ggplot(data=CarbOnlyMale) + 
  geom_point(mapping=aes(x=factor(Origin, level = Order), y=P_S, color=Subpathway)) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Individual Pathway Enrichment Score (P/S)", title = "Male Rhesus Monkey Carbohydrate Metabolic Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3") 
co14

######## FEMALES ########

CarbOnlyFemale <- filter(d3, Superpathway=="Carbohydrate")

cb3 <- ggplot(data=CarbOnlyFemale) + 
  geom_col(mapping=aes(x=factor(Origin, level = Order), y=P_S, fill=Subpathway), width = 0.7) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", title = "Carbohydrate Metabolism Subpathway Enrichment in Female Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "RdYlBu")+
  theme(axis.title.y =element_blank()) +
  coord_flip() 
cb3 
ggsave("figures/pathway_enrichment/females_CarbohydrateSubPathwayByTissue.png", height = 4 , width = 10)

#### Individual Data Points (with n) Tissue on X
co15 <- ggplot(data=CarbOnlyFemale) + 
  geom_point(mapping=aes(x=factor(Origin, level = Order), y=P_S, color=Subpathway)) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Individual Pathway Enrichment Score (P/S)", title = "Female Rhesus Monkey Carbohydrate Metabolism Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3") 
co15


######################### Energy  #########################

######## MALES & FEMALES  ########

EnergyAll <- filter(d1, Superpathway=="Energy")

e1 <- ggplot(data=EnergyAll) + 
  geom_col(mapping=aes(x=factor(Origin, level = Order), y=P_S, fill=Subpathway), width = 0.7) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", title = "Rhesus Monkey Energy Metabolism Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3")
e1 

#### Individual Data Points (with n) Tissue on X
co16 <- ggplot(data=EnergyAll) + 
  geom_point(mapping=aes(x=factor(Origin, level = Order), y=P_S, color=Subpathway)) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Individual Pathway Enrichment Score (P/S)", title = "Rhesus Monkey Energy Metabolic Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3") 
co16

######## MALES  ########

EnergyOnlyMale <- filter(d2, Superpathway=="Energy")

e2 <- ggplot(data=EnergyOnlyMale) + 
  geom_col(mapping=aes(x=factor(Origin, level = Order), y=P_S, fill=Subpathway), width = 0.7) + 
  theme_classic() + 
  labs(x="Sample Origine", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", title = "Male Rhesus Monkey Energy Metabolism Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3")
e2

#### Individual Data Points (with n) Tissue on X
co17 <- ggplot(data=EnergyOnlyMale) + 
  geom_point(mapping=aes(x=factor(Origin, level = Order), y=P_S, color=Subpathway)) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Individual Pathway Enrichment Score (P/S)", title = "Male Rhesus Monkey Energy Metabolism Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3") 
co17

######## FEMALES ########

EnergyOnlyFemale <- filter(d3, Superpathway=="Energy")

e3 <- ggplot(data=EnergyOnlyFemale) + 
  geom_col(mapping=aes(x=factor(Origin, level = Order), y=P_S, fill=Subpathway), width = 0.7) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", title = "Female Rhesus Monkey Energy Metabolism Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3")
e3 

#### Individual Data Points (with n) Tissue on X
co18 <- ggplot(data=EnergyOnlyFemale) + 
  geom_point(mapping=aes(x=factor(Origin, level = Order), y=P_S, color=Subpathway)) + 
  theme_classic() + 
  labs(x="Sample Origin", y="Individual Pathway Enrichment Score (P/S)", title = "Female Rhesus Monkey Energy Metabolism Subpathway Enrichment") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  scale_fill_brewer(palette = "Set3") 
co18




########################   ########################   ########################   


# Bar Charts for Top Subpathways within each tissue (for each Treatment group)
## ALL Tissues
# Combined
d1_Total <- read_xlsx('input_files/Pathway Enrichment.xlsx', sheet = 'TGC_forTotalEnrichPlot')
d1_Total$Subpathway <- as.character(d1_Total$Subpathway)
d1_Total$Subpathway <- factor(d1_Total$Subpathway, levels = unique(d1_Total$Subpathway))

# Male
d2_Total <- read_xlsx('input_files/Pathway Enrichment.xlsx', sheet = 'MC_forTotalEnrichPlot')
d2_Total$Subpathway <- as.character(d2_Total$Subpathway)
d2_Total$Subpathway <- factor(d2_Total$Subpathway, levels = unique(d2_Total$Subpathway))
MTG <- ggplot(data=d2_Total) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=Superpathway), width = 0.7) + 
  labs(x="Superpathways", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", 
       title = "Metabolic Pathway Enrichment in Male Treatment Condition") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme_classic() + 
  theme(axis.title.y =element_blank()) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(d2_Total$Subpathway))) +
  scale_fill_brewer(palette = "RdYlBu")
MTG 
ggsave("figures/pathway_enrichment/all_tissues_male_Treatment_Enrichment.png", height = 7, width = 10)

#Female
d3_Total <- read_xlsx('input_files/Pathway Enrichment.xlsx', sheet = 'FC_forTotalEnrichPlot')
d3_Total$Subpathway <- as.character(d3_Total$Subpathway)
d3_Total$Subpathway <- factor(d3_Total$Subpathway, levels = unique(d3_Total$Subpathway))
FTG <- ggplot(data=d3_Total) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=Superpathway), width = 0.7) + 
  labs(x="Superpathways", y="Summated Pathway Enrichment \nPoly(I:C) / Saline", title = "Metabolic Pathway Enrichment in Female Treatment Condition") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme_classic() + 
  theme(axis.title.y =element_blank()) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(d3_Total$Subpathway))) +
  scale_fill_brewer(palette = "RdYlBu")
FTG 
ggsave("figures/pathway_enrichment/all_tissues_female_Treatment_Enrichment.png", height = 7, width = 10)






########################################################################

## PLASMA
# Combined

PlasmaTreatment <- filter(d1, Origin=="Plasma")
PlasmaTreatment$Subpathway <- as.character(PlasmaTreatment$Subpathway)
PlasmaTreatment$Subpathway <- factor(PlasmaTreatment$Subpathway, levels = unique(PlasmaTreatment$Subpathway))

Pl1 <- ggplot(data=PlasmaTreatment) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=Superpathway), width = 0.7) + 
  labs(x="Tissue Origin", y=" Enrichment Score", title = "Plasma Pathway Enrichment within Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme_classic() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") +
  scale_x_discrete(limits = rev(levels(PlasmaTreatment$Subpathway))) +
  coord_flip()
Pl1 
ggsave("figures/pathway_enrichment/Plasma_combinedTreatment_Enrichment.png", height = 4, width = 9)


# Male

PlasmaTreatment <- filter(d2, Origin=="Plasma")
PlasmaTreatment$Subpathway <- as.character(PlasmaTreatment$Subpathway)
PlasmaTreatment$Subpathway <- factor(PlasmaTreatment$Subpathway, levels = unique(PlasmaTreatment$Subpathway))

mPl1 <- ggplot(data=PlasmaTreatment) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=Superpathway), width = 0.7) + 
  labs(x="Tissue Origin", y=" Enrichment Score", title = "Plasma Pathway Enrichment within Male Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme_classic() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") +
  scale_x_discrete(limits = rev(levels(PlasmaTreatment$Subpathway))) +
  coord_flip()
mPl1 
ggsave("figures/pathway_enrichment/Plasma_male_Treatment_Enrichment.png", height = 4, width = 9)


# Female

PlasmaTreatment <- filter(d3, Origin=="Plasma")
PlasmaTreatment$Subpathway <- as.character(PlasmaTreatment$Subpathway)
PlasmaTreatment$Subpathway <- factor(PlasmaTreatment$Subpathway, levels = unique(PlasmaTreatment$Subpathway))

mPl1 <- ggplot(data=PlasmaTreatment) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=Superpathway), width = 0.7) + 
  labs(x="Tissue Origin", y=" Enrichment Score", title = "Plasma Pathway Enrichment within Female Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme_classic() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") +
  scale_x_discrete(limits = rev(levels(PlasmaTreatment$Subpathway))) +
  coord_flip()
mPl1 
ggsave("figures/pathway_enrichment/Plasma_female_Treatment_Enrichment.png", height = 4, width = 9)



## CSF
# Combined

CSFTreatment <- filter(d1, Origin=="CSF")
CSFTreatment$Subpathway <- as.character(CSFTreatment$Subpathway)
CSFTreatment$Subpathway <- factor(CSFTreatment$Subpathway, levels = unique(CSFTreatment$Subpathway))

CS1 <- ggplot(data=CSFTreatment) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=Superpathway), width = 0.7) + 
  labs(x="Tissue Origin", y=" Enrichment Score", title = "CSF Pathway Enrichment within Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme_classic() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") +
  scale_x_discrete(limits = rev(levels(CSFTreatment$Subpathway))) +
  coord_flip()
CS1 
ggsave("figures/pathway_enrichment/CSF_combinedTreatment_Enrichment.png", height = 4, width = 9)


# Male

CSFTreatment <- filter(d2, Origin=="CSF")
CSFTreatment$Subpathway <- as.character(CSFTreatment$Subpathway)
CSFTreatment$Subpathway <- factor(CSFTreatment$Subpathway, levels = unique(CSFTreatment$Subpathway))

mCS1 <- ggplot(data=CSFTreatment) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=Superpathway), width = 0.7) + 
  labs(x="Tissue Origin", y=" Enrichment Score", title = "CSF Pathway Enrichment within Male Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme_classic() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") +
  scale_x_discrete(limits = rev(levels(CSFTreatment$Subpathway))) +
  coord_flip()
mCS1 
ggsave("figures/pathway_enrichment/CSF_male_Treatment_Enrichment.png", height = 4, width = 9)


# Female

CSFTreatment <- filter(d3, Origin=="CSF")
CSFTreatment$Subpathway <- as.character(CSFTreatment$Subpathway)
CSFTreatment$Subpathway <- factor(CSFTreatment$Subpathway, levels = unique(CSFTreatment$Subpathway))

fCS1 <- ggplot(data=CSFTreatment) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=Superpathway), width = 0.7) + 
  labs(x="Tissue Origin", y=" Enrichment Score", title = "CSF Pathway Enrichment within Female Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme_classic() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") +
  scale_x_discrete(limits = rev(levels(CSFTreatment$Subpathway))) +
  coord_flip()
fCS1 
ggsave("figures/pathway_enrichment/CSF_female_Treatment_Enrichment.png", height = 4, width = 9)



## JEJUNUM
# Combined

JejunumTreatment <- filter(d1, Origin=="Jejunum")
JejunumTreatment$Subpathway <- as.character(JejunumTreatment$Subpathway)
JejunumTreatment$Subpathway <- factor(JejunumTreatment$Subpathway, levels = unique(JejunumTreatment$Subpathway))

Jej1 <- ggplot(data=JejunumTreatment) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=Superpathway), width = 0.7) + 
  labs(x="Tissue Origin", y=" Enrichment Score", title = "Jejunum Pathway Enrichment within Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme_classic() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") +
  scale_x_discrete(limits = rev(levels(JejunumTreatment$Subpathway))) +
  coord_flip()
Jej1 
ggsave("figures/pathway_enrichment/Jejunum_combinedTreatment_Enrichment.png", height = 4, width = 9)


# Male

JejunumTreatment <- filter(d2, Origin=="Jejunum")
JejunumTreatment$Subpathway <- as.character(JejunumTreatment$Subpathway)
JejunumTreatment$Subpathway <- factor(JejunumTreatment$Subpathway, levels = unique(JejunumTreatment$Subpathway))

mJej1 <- ggplot(data=JejunumTreatment) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=Superpathway), width = 0.7) + 
  labs(x="Tissue Origin", y=" Enrichment Score", title = "Jejunum Pathway Enrichment within Male Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme_classic() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") +
  scale_x_discrete(limits = rev(levels(JejunumTreatment$Subpathway))) +
  coord_flip()
mJej1 
ggsave("figures/pathway_enrichment/Jejunum_male_Treatment_Enrichment.png", height = 4, width = 9)


# Female

JejunumTreatment <- filter(d3, Origin=="Jejunum")
JejunumTreatment$Subpathway <- as.character(JejunumTreatment$Subpathway)
JejunumTreatment$Subpathway <- factor(JejunumTreatment$Subpathway, levels = unique(JejunumTreatment$Subpathway))

fJej1 <- ggplot(data=JejunumTreatment) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=Superpathway), width = 0.7) + 
  labs(x="Tissue Origin", y=" Enrichment Score", title = "Jejunum Pathway Enrichment within Female Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme_classic() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") +
  scale_x_discrete(limits = rev(levels(JejunumTreatment$Subpathway))) +
  coord_flip()
fJej1 
ggsave("figures/pathway_enrichment/Jejunum_female_Treatment_Enrichment.png", height = 4, width = 9)



## ILEUM
# Combined

IleumTreatment <- filter(d1, Origin=="Ileum")
IleumTreatment$Subpathway <- as.character(IleumTreatment$Subpathway)
IleumTreatment$Subpathway <- factor(IleumTreatment$Subpathway, levels = unique(IleumTreatment$Subpathway))

Jej1 <- ggplot(data=IleumTreatment) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=Superpathway), width = 0.7) + 
  labs(x="Tissue Origin", y=" Enrichment Score", title = "Ileum Pathway Enrichment within Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme_classic() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") +
  scale_x_discrete(limits = rev(levels(IleumTreatment$Subpathway))) +
  coord_flip()
Jej1 
ggsave("figures/pathway_enrichment/Ileum_combinedTreatment_Enrichment.png", height = 4, width = 9)


# Male

IleumTreatment <- filter(d2, Origin=="Ileum")
IleumTreatment$Subpathway <- as.character(IleumTreatment$Subpathway)
IleumTreatment$Subpathway <- factor(IleumTreatment$Subpathway, levels = unique(IleumTreatment$Subpathway))

mJej1 <- ggplot(data=IleumTreatment) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=Superpathway), width = 0.7) + 
  labs(x="Tissue Origin", y=" Enrichment Score", title = "Ileum Pathway Enrichment within Male Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme_classic() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") +
  scale_x_discrete(limits = rev(levels(IleumTreatment$Subpathway))) +
  coord_flip()
mJej1 
ggsave("figures/pathway_enrichment/Ileum_male_Treatment_Enrichment.png", height = 4, width = 9)


# Female

IleumTreatment <- filter(d3, Origin=="Ileum")
IleumTreatment$Subpathway <- as.character(IleumTreatment$Subpathway)
IleumTreatment$Subpathway <- factor(IleumTreatment$Subpathway, levels = unique(IleumTreatment$Subpathway))

fJej1 <- ggplot(data=IleumTreatment) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=Superpathway), width = 0.7) + 
  labs(x="Tissue Origin", y=" Enrichment Score", title = "Ileum Pathway Enrichment within Female Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme_classic() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") +
  scale_x_discrete(limits = rev(levels(IleumTreatment$Subpathway))) +
  coord_flip()
fJej1 
ggsave("figures/pathway_enrichment/Ileum_female_Treatment_Enrichment.png", height = 4, width = 9)



## COLON
# Combined

ColonTreatment <- filter(d1, Origin=="Colon")
ColonTreatment$Subpathway <- as.character(ColonTreatment$Subpathway)
ColonTreatment$Subpathway <- factor(ColonTreatment$Subpathway, levels = unique(ColonTreatment$Subpathway))

Jej1 <- ggplot(data=ColonTreatment) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=Superpathway), width = 0.7) + 
  labs(x="Tissue Origin", y=" Enrichment Score", title = "Colon Pathway Enrichment within Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme_classic() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") +
  scale_x_discrete(limits = rev(levels(ColonTreatment$Subpathway))) +
  coord_flip()
Jej1 
ggsave("figures/pathway_enrichment/Colon_combinedTreatment_Enrichment.png", height = 4, width = 9)


# Male

ColonTreatment <- filter(d2, Origin=="Colon")
ColonTreatment$Subpathway <- as.character(ColonTreatment$Subpathway)
ColonTreatment$Subpathway <- factor(ColonTreatment$Subpathway, levels = unique(ColonTreatment$Subpathway))

mJej1 <- ggplot(data=ColonTreatment) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=Superpathway), width = 0.7) + 
  labs(x="Tissue Origin", y=" Enrichment Score", title = "Colon Pathway Enrichment within Male Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme_classic() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") +
  scale_x_discrete(limits = rev(levels(ColonTreatment$Subpathway))) +
  coord_flip()
mJej1 
ggsave("figures/pathway_enrichment/Colon_male_Treatment_Enrichment.png", height = 4, width = 9)


# Female

ColonTreatment <- filter(d3, Origin=="Colon")
ColonTreatment$Subpathway <- as.character(ColonTreatment$Subpathway)
ColonTreatment$Subpathway <- factor(ColonTreatment$Subpathway, levels = unique(ColonTreatment$Subpathway))

fJej1 <- ggplot(data=ColonTreatment) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=Superpathway), width = 0.7) + 
  labs(x="Tissue Origin", y=" Enrichment Score", title = "Colon Pathway Enrichment within Female Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme_classic() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") +
  scale_x_discrete(limits = rev(levels(ColonTreatment$Subpathway))) +
  coord_flip()
fJej1 
ggsave("figures/pathway_enrichment/Colon_female_Treatment_Enrichment.png", height = 4, width = 9)



## FECES
# Combined

FecesTreatment <- filter(d1, Origin=="Feces")
FecesTreatment$Subpathway <- as.character(FecesTreatment$Subpathway)
FecesTreatment$Subpathway <- factor(FecesTreatment$Subpathway, levels = unique(FecesTreatment$Subpathway))

Jej1 <- ggplot(data=FecesTreatment) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=Superpathway), width = 0.7) + 
  labs(x="Tissue Origin", y=" Enrichment Score", title = "Feces Pathway Enrichment within Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme_classic() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") +
  scale_x_discrete(limits = rev(levels(FecesTreatment$Subpathway))) +
  coord_flip()
Jej1 
ggsave("figures/pathway_enrichment/Feces_combinedTreatment_Enrichment.png", height = 4, width = 9)


# Male

FecesTreatment <- filter(d2, Origin=="Feces")
FecesTreatment$Subpathway <- as.character(FecesTreatment$Subpathway)
FecesTreatment$Subpathway <- factor(FecesTreatment$Subpathway, levels = unique(FecesTreatment$Subpathway))

mJej1 <- ggplot(data=FecesTreatment) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=Superpathway), width = 0.7) + 
  labs(x="Tissue Origin", y=" Enrichment Score", title = "Feces Pathway Enrichment within Male Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme_classic() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") +
  scale_x_discrete(limits = rev(levels(FecesTreatment$Subpathway))) +
  coord_flip()
mJej1 
ggsave("figures/pathway_enrichment/Feces_male_Treatment_Enrichment.png", height = 4, width = 9)


# Female

FecesTreatment <- filter(d3, Origin=="Feces")
FecesTreatment$Subpathway <- as.character(FecesTreatment$Subpathway)
FecesTreatment$Subpathway <- factor(FecesTreatment$Subpathway, levels = unique(FecesTreatment$Subpathway))

fJej1 <- ggplot(data=FecesTreatment) + 
  geom_col(mapping=aes(x=Subpathway, y=P_S, fill=Superpathway), width = 0.7) + 
  labs(x="Tissue Origin", y=" Enrichment Score", title = "Feces Pathway Enrichment within Female Treatment Group") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.95)) +
  theme_classic() + 
  theme(axis.title.y =element_blank()) +
  scale_fill_brewer(palette = "RdYlBu") +
  scale_x_discrete(limits = rev(levels(FecesTreatment$Subpathway))) +
  coord_flip()
fJej1 
ggsave("figures/pathway_enrichment/Feces_female_Treatment_Enrichment.png", height = 4, width = 9)



