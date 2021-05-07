source("src/_load_packages.R")
source("src/_misc_functions.R")


df <- read_xlsx('/Users/joeboktor/Documents/Rhesus\ Monkey\ MIA\ Metabolomics/NHP Plasma HD4 Metabolon.xlsx', sheet = 'LogData')
Ind <- read_xlsx('/Users/joeboktor/Documents/Rhesus\ Monkey\ MIA\ Metabolomics/NHP Plasma HD4 Metabolon.xlsx', sheet = 'ScaledImpData')

df <- mutate(df, thresh = if_else(`Negative Log(PM/SM p-value)` >= 1.3 & `Log2 PM/SM` > 0, "A",
                                  if_else( `Negative Log(PM/SM p-value)` >= 1.3 & `Log2 PM/SM` < 0, "B", "C" )))
df <- mutate(df, FEMthresh = if_else(`Negative Log(PF/SF p-value)` >= 1.3 & `Log2 PF/SF` > 0, "A",
                                     if_else(`Negative Log(PF/SF p-value)` >= 1.3 & `Log2 PF/SF` < 0, "B", "C" )))
df <- mutate(df, combthresh = if_else(`Negative Log(P/S p-value)` >= 1.3 & `Log2 P/S` > 0, "A",
                                      if_else(`Negative Log(P/S p-value)` >= 1.3 & `Log2 P/S` < 0, "B", "C" )))

# Create lists of significant molecules
MaleSigonly <- filter(df, thresh=="A" | thresh=="B" )
FEMSigonly <- filter(df, FEMthresh=="A" | FEMthresh=="B" )
CombSigonly<- filter(df, combthresh=="A" | combthresh=="B" )

AllPlasmaSig <- filter(df, thresh=="A" | thresh=="B" | FEMthresh=="A" | FEMthresh=="B" | combthresh=="A" | combthresh=="B") 
# write.csv(AllPlasmaSig, file = "TEMP_Plasma.csv")

################################################  CSF ################################################  


dfcs <- read_xlsx('/Users/joeboktor/Documents/Rhesus\ Monkey\ MIA\ Metabolomics/NHP CSF Metabolon.xlsx', sheet = 'LogData')
Indc <- read_xlsx('/Users/joeboktor/Documents/Rhesus\ Monkey\ MIA\ Metabolomics/NHP CSF Metabolon.xlsx', sheet = 'ScaledImpData')

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

AllCSFSig <- filter(dfcs, mthresh=="A" | mthresh=="B" | fthresh=="A" | fthresh=="B" | combthresh=="A" | combthresh=="B")

# write.csv(AllCSFSig, file = "TEMP_CSF.csv")






################################################  Jejunum  ################################################  

dfj <- read_xlsx('/Users/joeboktor/Documents/Rhesus\ Monkey\ MIA\ Metabolomics/NHP Jejunum Metabolon.xlsx', sheet = 'LogData')
Indj <- read_xlsx('/Users/joeboktor/Documents/Rhesus\ Monkey\ MIA\ Metabolomics/NHP CSF Metabolon.xlsx', sheet = 'ScaledImpData')


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

AlljejSig <- filter(dfj, mthresh=="A" | mthresh=="B" | fthresh=="A" | fthresh=="B" | combthresh=="A" | combthresh=="B")

# write.csv(AlljejSig, file = "TEMP_jej.csv")


################################################  Ileum  ################################################  


dfi <- read_xlsx('/Users/joeboktor/Documents/Rhesus\ Monkey\ MIA\ Metabolomics/NHP Ileum Metabolon.xlsx', sheet = 'LogData')
Indi <- read_xlsx('/Users/joeboktor/Documents/Rhesus\ Monkey\ MIA\ Metabolomics/NHP Ileum Metabolon.xlsx', sheet = 'ScaledImpData')

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

AllIlSig <- filter(dfi, mthresh=="A" | mthresh=="B" | fthresh=="A" | fthresh=="B" | combthresh=="A" | combthresh=="B")

# write.csv(AllIlSig, file = "TEMP_ILE.csv")

################################################  Colon  ################################################  

dfc <- read_xlsx('/Users/joeboktor/Documents/Rhesus\ Monkey\ MIA\ Metabolomics/NHP Colon Metabolon.xlsx', sheet = 'LogData')
Indc <- read_xlsx('/Users/joeboktor/Documents/Rhesus\ Monkey\ MIA\ Metabolomics/NHP Colon Metabolon.xlsx', sheet = 'ScaledImpData')


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

AllColSig <- filter(dfc, mthresh=="A" | mthresh=="B" | fthresh=="A" | fthresh=="B" | combthresh=="A" | combthresh=="B" )


# write.csv(AllColSig, file = "TEMP_Col.csv")

################################################  Feces  ################################################  


dffc <- read_xlsx('/Users/joeboktor/Documents/Rhesus\ Monkey\ MIA\ Metabolomics/NHP Feces Metabolon.xlsx', sheet = 'LogData')
Indfc <- read_xlsx('/Users/joeboktor/Documents/Rhesus\ Monkey\ MIA\ Metabolomics/NHP Feces Metabolon.xlsx', sheet = 'ScaledImpData')


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

AllFecSig <- filter(dffc, mthresh=="A" | mthresh=="B" | fthresh=="A" | fthresh=="B" | combthresh=="A" | combthresh=="B")

# write.csv(AllFecSig, file = "TEMP_FECES.csv")

################################################  VENN DIAGRAMS  ################################################  


Met <- read_xlsx('input_files/TissueSharedMetabolites.xlsx', sheet = 'ALL')

nrow(subset(Met, Plasma == 1))
nrow(subset(Met, CSF == 1))
nrow(subset(Met, Jejunum == 1))
nrow(subset(Met, Ileum == 1))
nrow(subset(Met, Colon == 1))
nrow(subset(Met, Feces == 1))

# Number of Metabolites Shared between Plasma and CSF
nrow(subset(Met, Plasma == 1 & CSF == 1))

# Checking metabolites shared with CSF
nrow(subset(Met, Feces == 1 & CSF == 1))

#Venn Diagram with three component
grid.newpage()
draw.triple.venn(area1 = nrow(subset(Met, Plasma == 1)), nrow(subset(Met, CSF == 1)), 
                 nrow(subset(Met, Ileum == 1)), n12 = nrow(subset(Met, Plasma == 1 & CSF == 1)), 
                 nrow(subset(Met, CSF == 1 & Ileum == 1)), n13 = nrow(subset(Met, Plasma == 1 & Ileum == 1)), 
                 n123 = nrow(subset(Met, Plasma == 1 & CSF == 1 & Ileum == 1)), 
                 category = c("Plasma Metabolites", "CSF Metabolites", "Ileum Metabolites"), 
                 lty = "blank", fill = c("skyblue", "pink1", "mediumorchid"))


fMet <- Met[, 2:7]

##### USING EULERR
set.seed(1)

NoCSF<- fMet[-2]
plot(venn(NoCSF))
plot(euler(NoCSF), quantities = TRUE)
plot(euler(NoCSF, shape = "ellipse"), quantities = TRUE)


# JUST GI 
jGI <- fMet[3:5]
plot(venn(jGI))
plot(euler(jGI, shape = "ellipse"), quantities = TRUE)

#  GI + Feces
fGI <- fMet[3:6]
plot(venn(fGI))

#  GI + Plasma
pGI <- NoCSF[1:4]
plot(venn(pGI))

#  CSF + Plasma
PC <- fMet[1:2]
plot(venn(PC))
plot(euler(PC), quantities = TRUE)

#  No jejunum
NoJej <- fMet[-3]
plot(venn(NoJej))
plot(euler(NoJej), quantities = TRUE)

#  No Plasma
NoP <- fMet[-1]
plot(venn(NoP))
plot(euler(NoP), quantities = TRUE)






