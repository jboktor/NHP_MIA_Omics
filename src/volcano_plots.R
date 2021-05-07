source("src/_load_packages.R")
source("src/_misc_functions.R")


######################## PLASMA MALE Poly IC/Saline ######################## 


df <- read_xlsx('input_files/NHP Plasma HD4 Metabolon.xlsx', sheet = 'LogData')
Ind <- read_xlsx('input_files/NHP Plasma HD4 Metabolon.xlsx', sheet = 'ScaledImpData')


df <- mutate(df, thresh = if_else(`Negative Log(PM/SM p-value)` >= 1.3 & `Log2 PM/SM` > 0, "A",
                                  if_else( `Negative Log(PM/SM p-value)` >= 1.3 & `Log2 PM/SM` < 0, "B", "C" )))
df <- mutate(df, FEMthresh = if_else(`Negative Log(PF/SF p-value)` >= 1.3 & `Log2 PF/SF` > 0, "A",
                                     if_else(`Negative Log(PF/SF p-value)` >= 1.3 & `Log2 PF/SF` < 0, "B", "C" )))
df <- mutate(df, combthresh = if_else(`Negative Log(P/S p-value)` >= 1.3 & `Log2 P/S` > 0, "A",
                                  if_else(`Negative Log(P/S p-value)` >= 1.3 & `Log2 P/S` < 0, "B", "C" )))

# Create list of significant molecules 
MaleSigonly <- filter(df, thresh=="A" | thresh=="B" )
FEMSigonly <- filter(df, FEMthresh=="A" | FEMthresh=="B" )
CombSigonly<- filter(df, combthresh=="A" | combthresh=="B" )

############################ PLASMA Volcano - ALL Significant MALE   ############################ 
ggplot(data=df, aes(x=`Log2 PM/SM` , y=`Negative Log(PM/SM p-value)`, color=thresh)) + geom_point() + 
  xlab("log2 fold change (Male PolyIC/Saline)") + ylab(expression(paste(-log[10], "[ P-value ]"))) +
  scale_colour_manual(values = c("A"= "red", "B"= "blue", "C"= "black")) +
  geom_label_repel(data=MaleSigonly, aes(label = `Biochemical Name`)) +
  geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = 1.00, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -1.00, linetype = 2, alpha = 0.5) +
  theme_classic() + 
  theme(legend.position="none") # Hide the legend

ggsave("figures/volcano_plots/MalePlasmaVolcanoPlot.svg", height = 10, width = 10)

############################  PLASMA Volcano - ALL Significant FEMALE   ############################ 
ggplot(data=df, aes(x=`Log2 PF/SF` , y=`Negative Log(PF/SF p-value)`, color=FEMthresh)) + geom_point() + 
  xlab("log2 fold change (Female PolyIC/Saline)") + ylab(expression(paste(-log[10], "[ P-value ]"))) + 
  scale_colour_manual(values = c("A"= "red", "B"= "blue", "C"= "black")) +
  geom_label_repel(data=FEMSigonly, aes(label = `Biochemical Name`)) +
  geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = 1.00, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -1.00, linetype = 2, alpha = 0.5) +
  theme_classic() + 
  theme(legend.position="none") # Hide the legend

ggsave("figures/volcano_plots/FemalePlasmaVolcanoPlot.svg", height = 10, width =10)

############################  PLASMA Volcano - ALL Significant Male + Female   ############################ 
ggplot(data=df, aes(x=`Log2 P/S`, y=`Negative Log(P/S p-value)`, color=combthresh)) + geom_point() + 
  xlab(expression(paste(log[2], "[ Poly(I:C)/Saline ]"))) + ylab(expression(paste(-log[10], "[ P-value ]"))) + 
  scale_colour_manual(values = c("A"= "red", "B"= "blue", "C"= "black")) +
  geom_label_repel(data=CombSigonly, aes(label = `Biochemical Name`)) +
  geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = 1.00, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -1.00, linetype = 2, alpha = 0.5) +
  theme_classic() + 
  theme(legend.position="none") # Hide the legend

ggsave("figures/volcano_plots/AllPlasmaVolcanoPlot.svg", height = 10, width =10)


################################################  CSF ################################################  


dfcs <- read_xlsx('input_files/NHP CSF Metabolon.xlsx', sheet = 'LogData')
Indc <- read_xlsx('input_files/NHP CSF Metabolon.xlsx', sheet = 'ScaledImpData')

dfcs <- mutate(dfcs, mthresh = if_else(`Negative Log(PM/SM p-value)` >= 1.3 & `Log2 PM/SM` > 0, "A",
                                  if_else( `Negative Log(PM/SM p-value)` >= 1.3 & `Log2 PM/SM` < 0, "B", "C" )))
dfcs <- mutate(dfcs, fthresh = if_else(`Negative Log(PF/SF p-value)` >= 1.3 & `Log2 PF/SF` > 0, "A",
                                     if_else(`Negative Log(PF/SF p-value)` >= 1.3 & `Log2 PF/SF` < 0, "B", "C" )))
dfcs <- mutate(dfcs, combthresh = if_else(`Negative Log(P/S p-value)` >= 1.3 & `Log2 P/S` > 0, "A",
                                      if_else(`Negative Log(P/S p-value)` >= 1.3 & `Log2 P/S` < 0, "B", "C" )))

# Create list of significant molecules
CSFmsig <- filter(dfcs, mthresh=="A" | mthresh=="B" )
CSFfsig <- filter(dfcs, fthresh=="A" | fthresh=="B" )
CSFCombSigonly<- filter(dfcs, combthresh=="A" | combthresh=="B" )


############################ CSF Volcano - ALL Significant MALE   ############################ 


ggplot(data=dfcs, aes(x=`Log2 PM/SM`, y=`Negative Log(PM/SM p-value)`, color=mthresh)) + geom_point() + 
  xlab(expression(paste(log[2], "[ Male Poly(I:C)/Saline ]"))) + ylab(expression(paste(-log[10], "[ P-value ]"))) +
  scale_colour_manual(values = c("A"= "red", "B"= "blue", "C"= "black")) +
  geom_label_repel(data=CSFmsig, aes(label=`Biochemical Name`)) +
  geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = 1.00, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -1.00, linetype = 2, alpha = 0.5) +
  theme_classic() + 
  theme(legend.position="none") # Hide the legend
ggsave("figures/volcano_plots/MaleCSFVolcanoPlot.svg", height = 10, width =10)

############################  CSF Volcano - ALL Significant FEMALE   ############################ 
ggplot(data=dfcs, aes(x=`Log2 PF/SF` , y=`Negative Log(PF/SF p-value)`, color=fthresh)) + geom_point() + 
  xlab(expression(paste(log[2], "[ Female Poly(I:C)/Saline ]"))) + ylab(expression(paste(-log[10], "[ P-value ]"))) + 
  scale_colour_manual(values = c("A"= "red", "B"= "blue", "C"= "black")) +
  geom_label_repel(data=CSFfsig, aes(label = `Biochemical Name`)) +
  geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = 1.00, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -1.00, linetype = 2, alpha = 0.5) +
  theme_classic() + 
  theme(legend.position="none") # Hide the legend
ggsave("figures/volcano_plots/FemaleCSFVolcanoPlot.svg", height = 10, width =10)


############################  CSF Volcano - ALL Significant Male + Female   ############################ 
ggplot(data=dfcs, aes(x=`Log2 P/S`, y=`Negative Log(P/S p-value)`, color=combthresh)) + geom_point() + 
  xlab(expression(paste(log[2], "[ Poly(I:C)/Saline ]"))) + ylab(expression(paste(-log[10], "[ P-value ]"))) + 
  scale_colour_manual(values = c("A"= "red", "B"= "blue", "C"= "black")) +
  geom_label_repel(data=CSFCombSigonly, aes(label = `Biochemical Name`)) +
  geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = 1.00, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -1.00, linetype = 2, alpha = 0.5) +
  theme_classic() + 
  theme(legend.position="none") # Hide the legend

ggsave("figures/volcano_plots/AllCSFVolcanoPlot.svg", height = 10, width =10)



################################################  Jejunum  ################################################  


dfj <- read_xlsx('input_files/NHP Jejunum Metabolon.xlsx', sheet = 'LogData')
Indj <- read_xlsx('input_files/NHP CSF Metabolon.xlsx', sheet = 'ScaledImpData')


dfj <- mutate(dfj, mthresh = if_else(`Negative Log(PM/SM p-value)` >= 1.3 & `Log2 PM/SM` > 0, "A",
                                     if_else( `Negative Log(PM/SM p-value)` >= 1.3 & `Log2 PM/SM` < 0, "B", "C" )))
dfj <- mutate(dfj, fthresh = if_else(`Negative Log(PF/SF p-value)` >= 1.3 & `Log2 PF/SF` > 0, "A",
                                     if_else(`Negative Log(PF/SF p-value)` >= 1.3 & `Log2 PF/SF` < 0, "B", "C" )))
dfj <- mutate(dfj, combthresh = if_else(`Negative Log(P/S p-value)` >= 1.3 & `Log2 P/S` > 0, "A",
                                          if_else(`Negative Log(P/S p-value)` >= 1.3 & `Log2 P/S` < 0, "B", "C" )))

# Create list of significant molecules
Jejmsig <- filter(dfj, mthresh=="A" | mthresh=="B" )
Jejfsig <- filter(dfj, fthresh=="A" | fthresh=="B" )
JejCombSigonly<- filter(dfj, combthresh=="A" | combthresh=="B" )

############################  Jejunum - ALL Significant MALE   ############################ 

ggplot(data=dfj, aes(x=`Log2 PM/SM`, y=`Negative Log(PM/SM p-value)`, color=mthresh)) + geom_point() + 
  xlab(expression(paste(log[2], "[ Male Poly(I:C)/Saline ]"))) + ylab(expression(paste(-log[10], "[ P-value ]"))) +
  scale_colour_manual(values = c("A"= "red", "B"= "blue", "C"= "black")) +
  geom_label_repel(data=Jejmsig, aes(label=`Biochemical Name`)) +
  geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = 1.00, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -1.00, linetype = 2, alpha = 0.5) +
  theme_classic() + 
  theme(legend.position="none") # Hide the legend
ggsave("figures/volcano_plots/MaleJejunumVolcanoPlot.svg", height = 10, width =10)

############################   Jejunum - ALL Significant FEMALE   ############################ 
ggplot(data=dfj, aes(x=`Log2 PF/SF` , y=`Negative Log(PF/SF p-value)`, color=fthresh)) + geom_point() + 
  xlab(expression(paste(log[2], "[ Female Poly(I:C)/Saline ]"))) + ylab(expression(paste(-log[10], "[ P-value ]"))) + 
  scale_colour_manual(values = c("A"= "red", "B"= "blue", "C"= "black")) +
  geom_label_repel(data=Jejfsig, aes(label = `Biochemical Name`)) +
  geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = 1.00, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -1.00, linetype = 2, alpha = 0.5) +
  theme_classic() + 
  theme(legend.position="none") # Hide the legend
ggsave("figures/volcano_plots/FemaleJejunumVolcanoPlot.svg", height = 10, width =10)


############################  Jejunum Volcano - ALL Significant Male + Female   ############################ 
ggplot(data=dfj, aes(x=`Log2 P/S`, y=`Negative Log(P/S p-value)`, color=combthresh)) + geom_point() + 
  xlab(expression(paste(log[2], "[ Poly(I:C)/Saline ]"))) + ylab(expression(paste(-log[10], "[ P-value ]"))) + 
  scale_colour_manual(values = c("A"= "red", "B"= "blue", "C"= "black")) +
  geom_label_repel(data=JejCombSigonly, aes(label = `Biochemical Name`)) +
  geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = 1.00, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -1.00, linetype = 2, alpha = 0.5) +
  theme_classic() + 
  theme(legend.position="none") # Hide the legend

ggsave("figures/volcano_plots/AllJejVolcanoPlot.svg", height = 10, width =10)




################################################  Ileum  ################################################  


dfi <- read_xlsx('input_files/NHP Ileum Metabolon.xlsx', sheet = 'LogData')
Indi <- read_xlsx('input_files/NHP Ileum Metabolon.xlsx', sheet = 'ScaledImpData')



dfi <- mutate(dfi, mthresh = if_else(`Negative Log(PM/SM p-value)` >= 1.3 & `Log2 PM/SM` > 0, "A",
                                     if_else( `Negative Log(PM/SM p-value)` >= 1.3 & `Log2 PM/SM` < 0, "B", "C" )))
dfi <- mutate(dfi, fthresh = if_else(`Negative Log(PF/SF p-value)` >= 1.3 & `Log2 PF/SF` > 0, "A",
                                     if_else(`Negative Log(PF/SF p-value)` >= 1.3 & `Log2 PF/SF` < 0, "B", "C" )))
dfi <- mutate(dfi, combthresh = if_else(`Negative Log(P/S p-value)` >= 1.3 & `Log2 P/S` > 0, "A",
                                        if_else(`Negative Log(P/S p-value)` >= 1.3 & `Log2 P/S` < 0, "B", "C" )))

# Create list of significant molecules
Ilmsig <- filter(dfi, mthresh=="A" | mthresh=="B" )
Ilfsig <- filter(dfi, fthresh=="A" | fthresh=="B" )
IlCombSigonly<- filter(dfi, combthresh=="A" | combthresh=="B" )

############################  Ileum - ALL Significant MALE   ############################ 

ggplot(data=dfi, aes(x=`Log2 PM/SM`, y=`Negative Log(PM/SM p-value)`, color=mthresh)) + geom_point() + 
  xlab(expression(paste(log[2], "[ Male Poly(I:C)/Saline ]"))) + ylab(expression(paste(-log[10], "[ P-value ]"))) +
  scale_colour_manual(values = c("A"= "red", "B"= "blue", "C"= "black")) +
  geom_label_repel(data=Ilmsig, aes(label=`Biochemical Name`)) +
  geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = 1.00, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -1.00, linetype = 2, alpha = 0.5) +
  theme_classic() + 
  theme(legend.position="none") # Hide the legend
ggsave("figures/volcano_plots/MaleIleumVolcanoPlot.svg", height = 10, width =10)

############################   Ileum - ALL Significant FEMALE   ############################ 
ggplot(data=dfi, aes(x=`Log2 PF/SF` , y=`Negative Log(PF/SF p-value)`, color=fthresh)) + geom_point() + 
  xlab(expression(paste(log[2], "[ Female Poly(I:C)/Saline ]"))) + ylab(expression(paste(-log[10], "[ P-value ]"))) + 
  scale_colour_manual(values = c("A"= "red", "B"= "blue", "C"= "black")) +
  geom_label_repel(data=Ilfsig, aes(label = `Biochemical Name`)) +
  geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = 1.00, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -1.00, linetype = 2, alpha = 0.5) +
  theme_classic() + 
  theme(legend.position="none") # Hide the legend
ggsave("figures/volcano_plots/FemaleIleumVolcanoPlot.svg", height = 10, width =10)

############################  Ileum Volcano - ALL Significant Male + Female   ############################ 
ggplot(data=dfi, aes(x=`Log2 P/S`, y=`Negative Log(P/S p-value)`, color=combthresh)) + geom_point() + 
  xlab(expression(paste(log[2], "[ Poly(I:C)/Saline ]"))) + ylab(expression(paste(-log[10], "[ P-value ]"))) + 
  scale_colour_manual(values = c("A"= "red", "B"= "blue", "C"= "black")) +
  geom_label_repel(data=IlCombSigonly, aes(label = `Biochemical Name`)) +
  geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = 1.00, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -1.00, linetype = 2, alpha = 0.5) +
  theme_classic() + 
  theme(legend.position="none") # Hide the legend
ggsave("figures/volcano_plots/AllIleumVolcanoPlot.svg", height = 10, width =10)

################################################  Colon  ################################################  


dfc <- read_xlsx('input_files/NHP Colon Metabolon.xlsx', sheet = 'LogData')
Indc <- read_xlsx('input_files/NHP Colon Metabolon.xlsx', sheet = 'ScaledImpData')


dfc <- mutate(dfc, mthresh = if_else(`Negative Log(PM/SM p-value)` >= 1.3 & `Log2 PM/SM` > 0, "A",
                                     if_else( `Negative Log(PM/SM p-value)` >= 1.3 & `Log2 PM/SM` < 0, "B", "C" )))
dfc <- mutate(dfc, fthresh = if_else(`Negative Log(PF/SF p-value)` >= 1.3 & `Log2 PF/SF` > 0, "A",
                                     if_else(`Negative Log(PF/SF p-value)` >= 1.3 & `Log2 PF/SF` < 0, "B", "C" )))
dfc <- mutate(dfc, combthresh = if_else(`Negative Log(P/S p-value)` >= 1.3 & `Log2 P/S` > 0, "A",
                                        if_else(`Negative Log(P/S p-value)` >= 1.3 & `Log2 P/S` < 0, "B", "C" )))

# Create list of significant molecules
Cmsig <- filter(dfc, mthresh=="A" | mthresh=="B" )
Cfsig <- filter(dfc, fthresh=="A" | fthresh=="B" )
ColCombSigonly<- filter(dfc, combthresh=="A" | combthresh=="B" )

############################  Colon - ALL Significant MALE   ############################ 

ggplot(data=dfc, aes(x=`Log2 PM/SM`, y=`Negative Log(PM/SM p-value)`, color=mthresh)) + geom_point() + 
  xlab(expression(paste(log[2], "[ Male Poly(I:C)/Saline ]"))) + ylab(expression(paste(-log[10], "[ P-value ]"))) +
  scale_colour_manual(values = c("A"= "red", "B"= "blue", "C"= "black")) +
  geom_label_repel(data=Cmsig, aes(label=`Biochemical Name`)) +
  geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = 1.00, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -1.00, linetype = 2, alpha = 0.5) +
  theme_classic() + 
  theme(legend.position="none") # Hide the legend
ggsave("figures/volcano_plots/MaleColonVolcanoPlot.svg", height = 10, width =10)

############################   Colon - ALL Significant FEMALE   ############################ 
ggplot(data=dfc, aes(x=`Log2 PF/SF` , y=`Negative Log(PF/SF p-value)`, color=fthresh)) + geom_point() + 
  xlab(expression(paste(log[2], "[ Female Poly(I:C)/Saline ]"))) + ylab(expression(paste(-log[10], "[ P-value ]"))) + 
  scale_colour_manual(values = c("A"= "red", "B"= "blue", "C"= "black")) +
  geom_label_repel(data=Cfsig, aes(label = `Biochemical Name`)) +
  geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = 1.00, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -1.00, linetype = 2, alpha = 0.5) +
  theme_classic() + 
  theme(legend.position="none") # Hide the legend
ggsave("figures/volcano_plots/FemaleColonVolcanoPlot.svg", height = 10, width =10)

############################  Colon Volcano - ALL Significant Male + Female   ############################ 
ggplot(data=dfc, aes(x=`Log2 P/S`, y=`Negative Log(P/S p-value)`, color=combthresh)) + geom_point() + 
  xlab(expression(paste(log[2], "[ Poly(I:C)/Saline ]"))) + ylab(expression(paste(-log[10], "[ P-value ]"))) + 
  scale_colour_manual(values = c("A"= "red", "B"= "blue", "C"= "black")) +
  geom_label_repel(data=ColCombSigonly, aes(label = `Biochemical Name`)) +
  geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = 1.00, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -1.00, linetype = 2, alpha = 0.5) +
  theme_classic() + 
  theme(legend.position="none") # Hide the legend
ggsave("figures/volcano_plots/AllColonVolcanoPlot.svg", height = 10, width =10)



################################################  Feces  ################################################  


dffc <- read_xlsx('input_files/NHP Feces Metabolon.xlsx', sheet = 'LogData')
Indfc <- read_xlsx('input_files/NHP Feces Metabolon.xlsx', sheet = 'ScaledImpData')


dffc <- mutate(dffc, mthresh = if_else(`Negative Log(PM/SM p-value)` >= 1.3 & `Log2 PM/SM` > 0, "A",
                                     if_else( `Negative Log(PM/SM p-value)` >= 1.3 & `Log2 PM/SM` < 0, "B", "C" )))
dffc <- mutate(dffc, fthresh = if_else(`Negative Log(PF/SF p-value)` >= 1.3 & `Log2 PF/SF` > 0, "A",
                                     if_else(`Negative Log(PF/SF p-value)` >= 1.3 & `Log2 PF/SF` < 0, "B", "C" )))
dffc <- mutate(dffc, combthresh = if_else(`Negative Log(P/S p-value)` >= 1.3 & `Log2 P/S` > 0, "A",
                                        if_else(`Negative Log(P/S p-value)` >= 1.3 & `Log2 P/S` < 0, "B", "C" )))

# Create list of significant molecules
Fcmsig <- filter(dffc, mthresh=="A" | mthresh=="B" )
Fcfsig <- filter(dffc, fthresh=="A" | fthresh=="B" )
FcCombSigonly<- filter(dffc, combthresh=="A" | combthresh=="B" )

############################  Feces - ALL Significant MALE   ############################ 

ggplot(data=dffc, aes(x=`Log2 PM/SM`, y=`Negative Log(PM/SM p-value)`, color=mthresh)) + geom_point() + 
  xlab(expression(paste(log[2], "[ Male Poly(I:C)/Saline ]"))) + ylab(expression(paste(-log[10], "[ P-value ]"))) +
  scale_colour_manual(values = c("A"= "red", "B"= "blue", "C"= "black")) +
  geom_label_repel(data=Fcmsig, aes(label=`Biochemical Name`)) +
  geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = 1.00, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -1.00, linetype = 2, alpha = 0.5) +
  theme_classic() + 
  theme(legend.position="none") # Hide the legend
ggsave("figures/volcano_plots/MaleFecesVolcanoPlot.svg", height = 10, width =10)

############################   Feces - ALL Significant FEMALE   ############################ 
ggplot(data=dffc, aes(x=`Log2 PF/SF` , y=`Negative Log(PF/SF p-value)`, color=fthresh)) + geom_point() + 
  xlab("log2 fold change (Female PolyIC/Saline)") + ylab(expression(paste(-log[10], "[ P-value ]"))) + 
  scale_colour_manual(values = c("A"= "red", "B"= "blue", "C"= "black")) +
  geom_label_repel(data=Fcfsig, aes(label = `Biochemical Name`), ) +
  geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = 1.00, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -1.00, linetype = 2, alpha = 0.5) +
  theme_classic() + 
  theme(legend.position="none") # Hide the legend
ggsave("figures/volcano_plots/FemaleFecesVolcanoPlot.svg", height = 10, width =10)

############################  Feces Volcano - ALL Significant Male + Female   ############################ 
ggplot(data=dffc, aes(x=`Log2 P/S`, y=`Negative Log(P/S p-value)`, color=combthresh)) + geom_point() + 
  xlab(expression(paste(log[2], "[ Poly(I:C)/Saline ]"))) + ylab(expression(paste(-log[10], "[ P-value ]"))) + 
  scale_colour_manual(values = c("A"= "red", "B"= "blue", "C"= "black")) +
  geom_label_repel(data=FcCombSigonly, aes(label = `Biochemical Name`)) +
  geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = 1.00, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -1.00, linetype = 2, alpha = 0.5) +
  theme_classic() + 
  theme(legend.position="none") # Hide the legend

ggsave("figures/volcano_plots/AllFecesVolcanoPlot.svg", height = 10, width =10)


