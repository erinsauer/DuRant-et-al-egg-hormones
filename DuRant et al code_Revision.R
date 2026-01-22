#############################################################################
#############   Code for DuRant et al. Diet-MG egg yolks  ###################
#############      developed by Erin L. Sauer             ###################
#############################################################################

############## load packages ################################
library(tidyverse)
library(paletteer)
library(lme4)
library(nlme)
library(car)
library(emmeans)
library(gridExtra)
library(DescTools)
library(factoextra)
library(ggpubr)
library(car)
library(emmeans)

############# tidy data #####################################
horm <- read.csv("hormones.csv")
str(horm)
colnames(horm)

# viz of all hormones and create summary table#
hormbp <- horm
horm.hist <- horm[c(15:24,29,36:37,39:41)]
for (column in horm.hist) {
  dev.new()
  hist(column)
}
sum.table <- as.data.frame(summary(horm.hist))
View(sum.table)
write.csv(sum.table, "hormone summary table.csv")

# remove hormones where <5 samples had concentrations >0 and scale
horm.list <- scale(horm[c(15:18,20:24,36:37,39:41)])
summary(horm.list)
horm.pred <- horm[1:12]
horm2 <- cbind(horm.pred,horm.list)
colnames(horm2) #remaining hormones are [13:26]
horm2 <- horm2 %>% mutate(MGD = paste(Treatment, Diet, sep="_"))

#more conservative hormone list without zero inflated groups
colnames(horm2)
horm2 <- (horm2[,c(1:14,16:18,20,25:27)])
#dimentions should be 122, 21
dim(horm2)

############## creation of BCI and FI   ####################
#mom data
moms <- read.csv("moms2.csv")
str(moms)

#body condition index
BCIm <- lm(Mass.y ~ Fat.y, data=moms)
summary(BCIm)
moms$BCI <- residuals(BCIm)
summary(moms$BCI)

#fecundity index
FIm <- lm(Mass.y ~ csize, data=moms)
summary(FIm)
moms$FI <- residuals(FIm)
summary(moms$FI)

################ merge moms and horm 2 ###########################
colnames(moms)
moms2 <- (moms[,c(7,8,11,12,13,17, 20, 21)])
nrow(horm2);nrow(moms2)
unique(moms2$Bird.ID)
unique(horm2$Bird.ID)

horm3 <- merge(moms2, horm2, by="Bird.ID")
unique(horm3$Bird.ID);nrow(horm3)
colnames(horm3)

#make total clutch mass
horm3$clutch.mass <- (horm3$csize * horm3$Egg.Mass.Day.Laid)

############## PCA of egg yolk hormones #####################
PCdb <- horm2[c(13:20)]
colnames(PCdb)
colnames(PCdb)[7] <- "Pregnanedione"  
colnames(PCdb)[6] <- "17a.hydroxypregnenolone"
res.pca <- prcomp(PCdb, scale = TRUE)
str(res.pca)
summary(res.pca)
eigenvec <- as.data.frame(res.pca$rotation)
write.csv(eigenvec,file=
            "PCA_eigenvec.csv")
fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             geom="point",
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

PCA.dim <-fviz_pca_var(res.pca,
                       col.var = "contrib", # Color by contributions to the PC
                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                       labelsize=6,
                       repel = TRUE# Avoid text overlapping
                       
)+labs(title ="", x = "PC1 (41%)", y = "PC2 (24%)")+
  font("xlab",size=25)+
  font("ylab",size=25)+
  font("xy.text",size=20)+
  font("legend.title",size=15)+
  font("legend.text",size=15)
#extract PCs
PCs <- res.pca[["x"]]
hormPC <- cbind(horm2,PCs)
head(hormPC)
hormPC2 <- merge(moms2, hormPC, by="Bird.ID")
hormPC2$clutch.mass <- (hormPC2$csize * hormPC2$Egg.Mass.Day.Laid)

######## multivariate PC model ####################
part1 <- cbind(hormPC2, response = (hormPC2$PC1), whichpart = factor(1))
part2 <- cbind(hormPC2, response = (hormPC2$PC2), whichpart = factor(2))
part3 <- cbind(hormPC2, response = (hormPC2$PC3), whichpart = factor(3))

PC_disaggregated <- rbind(part1, part2, part3)
PC_disaggregated <-
  PC_disaggregated[order(PC_disaggregated$Bird.ID),]
names(PC_disaggregated)
rownames(PC_disaggregated) <- 1:dim(PC_disaggregated)[1]
dim(PC_disaggregated) 
#remove nas if there are any
PC_disaggregated <- PC_disaggregated[!is.na(PC_disaggregated$response),]
dim(PC_disaggregated) #there are none

PC.model <- lme(response ~ whichpart - 1 +
                    whichpart:(Diet*Treatment+Yolk.Mass+Egg.+Mass.y 
                                 ),
                  random = ~ whichpart -1 |Bird.ID,
                  weights=varIdent(form=~1|whichpart),
                  corr=corSymm(form=~1|Bird.ID/Tube.number),
                  data=PC_disaggregated, method="ML",
                  control = list(maxIter=500, msMaxIter=500, tolerance=1e-8,
                                 niterEM=250))
AIC(PC.model) 
summary(PC.model)
VarCorr(PC.model)
#extract model tables
horm.diet <- emmeans(PC.model, ~ Diet , at = list(whichpart = c("1","2","3")))
diet.table <- as.data.frame(horm.diet)
diet.table
horm.mg <- emmeans(PC.model, ~ Treatment , at = list(whichpart = c("1","2","3")))
mg.table <- as.data.frame(horm.mg)
mg.table
horm.interaction <- emmeans(PC.model, ~ Diet|Treatment, at = list(whichpart = c("1","2","3")))
interaction.table <- as.data.frame(horm.interaction)
interaction.table
interaction.table <- interaction.table %>% mutate(MGD = paste(Treatment, Diet, sep="_"))
interaction.table.PC1 <- subset(interaction.table, whichpart=="1")
interaction.table.PC2 <- subset(interaction.table, whichpart=="2")
interaction.table.PC3 <- subset(interaction.table, whichpart=="3")

PC.model.PC1 <- emmeans(PC.model, ~ Diet*Treatment, 
                       at = list(whichpart = c("1")))
PC.model.PC1.tukey <- as.data.frame(pairs(PC.model.PC1, adjust="tukey"))
PC.model.PC2 <- emmeans(PC.model, ~ Diet*Treatment, 
                        at = list(whichpart = c("2")))
PC.model.PC2.tukey <- as.data.frame(pairs(PC.model.PC2, adjust="tukey"))
PC.model.PC3 <- emmeans(PC.model, ~ Diet*Treatment, 
                        at = list(whichpart = c("3")))
PC.model.PC3.tukey <- as.data.frame(pairs(PC.model.PC3, adjust="tukey"))

PC1 <- ggplot(interaction.table.PC1, aes(x=MGD, y=emmean, color=Diet)) +
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL)) +
  theme_bw()+
  theme(axis.text=element_text(size=20, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=25),
        legend.position="right",
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        plot.title=element_text(size=30, hjust=0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab("Effect size")+xlab("Control     MG")+ggtitle("PC1")
PC2 <- ggplot(interaction.table.PC2, aes(x=MGD, y=emmean, color=Diet)) +
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL)) +
  theme_bw()+
  theme(axis.text=element_text(size=20, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=25),
        legend.position="right",
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        plot.title=element_text(size=30, hjust=0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab("Effect size")+xlab("Control     MG")+ggtitle("PC2")
PC3 <- ggplot(interaction.table.PC3, aes(x=MGD, y=emmean, color=Diet)) +
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL)) +
  theme_bw()+
  theme(axis.text=element_text(size=20, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=25),
        legend.position="right",
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        plot.title=element_text(size=30, hjust=0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab("Effect size")+xlab("Control     MG")+ggtitle("PC3")
PC.panel <- ggarrange(PC1,PC2,PC3,PCA.dim, ncol=2, nrow=2,
                      labels=c("A","B","C","D"), font.label=list(size=25))
PC.panel 

############# Androgens ######################
part1 <- cbind(horm3, response = (horm3$Androstenedione), whichpart = factor(1))
part2 <- cbind(horm3, response = (horm3$Testosterone), whichpart = factor(2))
part3 <- cbind(horm3, response = (horm3$Etiocholanolone), whichpart = factor(3))

andro_disaggregated <- rbind(part1, part2, part3)
andro_disaggregated <- 
  andro_disaggregated[order(andro_disaggregated$Bird.ID),]
names(andro_disaggregated)
rownames(andro_disaggregated) <- 1:dim(andro_disaggregated)[1]
# dimensions should be 336, 36
dim(andro_disaggregated)

andro.model <- lme(response ~ whichpart - 1 + 
                   whichpart:(Diet*Treatment+Yolk.Mass+Egg.+
                                +Fat.y+csize),  
                   random = ~ whichpart -1 |Bird.ID, 
                   weights=varIdent(form=~1|whichpart),
                   corr=corSymm(form=~as.numeric(whichpart)|Bird.ID/Tube.number),
                   data=andro_disaggregated, method="ML",
                   control = list(maxIter=500, msMaxIter=500, tolerance=1e-8,
                                  niterEM=250))
AIC(andro.model) 
summary(andro.model)
VarCorr(andro.model)
#extract model tables
andro.diet <- emmeans(andro.model, ~ Diet*Treatment, 
                      at = list(whichpart = c("1","2" ,"3")))
andro.diet1 <- emmeans(andro.model, ~ Diet*Treatment, 
                      at = list(whichpart = c("1")))
andro.diet2 <- emmeans(andro.model, ~ Diet*Treatment, 
                      at = list(whichpart = c("2")))
andro.diet3 <- emmeans(andro.model, ~ Diet*Treatment, 
                      at = list(whichpart = c("3")))
andro.diet1.tukey <- as.data.frame(pairs(andro.diet1, adjust="tukey"))
andro.diet2.tukey <- as.data.frame(pairs(andro.diet2, adjust="tukey"))
andro.diet3.tukey <- as.data.frame(pairs(andro.diet3, adjust="tukey"))
write.csv(andro.diet1.tukey, file="andro_tukey_table.csv")
write.csv(andro.diet2.tukey, file="testo_tukey_table.csv")
write.csv(andro.diet3.tukey, file="ei_tukey_table.csv")
#save tukey

#save emmeans
summary(andro.diet)
andro.table <- as.data.frame(andro.diet)
andro.table
write.csv(andro.table, file="andro_table.csv")
andro.table <- andro.table %>% mutate(MGD = paste(Treatment, Diet, sep="_"))
andro.table.andro <- subset(andro.table, whichpart=="1")
andro.table.testo <- subset(andro.table, whichpart=="2")
andro.table.ei <- subset(andro.table, whichpart=="3")

andro.fig <- ggplot(andro.table.andro, aes(x=MGD, y=emmean, color=Diet)) +
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL)) +
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none",
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        plot.title=element_text(size=20, hjust=0.5))+
  ylab("")+xlab("")+ggtitle("Androstenedione")+
  scale_x_discrete(labels=c("Control_Lipid" = "Control \n lipid", 
                            "Control_Protein" = "Control \n protein", 
                            "MG_Lipid" = "MG \n lipid", 
                            "MG_Protein" = "MG \n protein"))
testo.fig <- ggplot(andro.table.testo, aes(x=MGD, y=emmean, color=Diet)) +
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL)) +
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none",
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        plot.title=element_text(size=20, hjust=0.5))+
  ylab("Effect size")+xlab("")+ggtitle("Testosterone")+
  scale_x_discrete(labels=c("Control_Lipid" = "Control \n lipid", 
                            "Control_Protein" = "Control \n protein", 
                            "MG_Lipid" = "MG \n lipid", 
                            "MG_Protein" = "MG \n protein"))
etio.fig <- ggplot(andro.table.ei, aes(x=MGD, y=emmean, color=Diet)) +
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL)) +
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none",
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        plot.title=element_text(size=20, hjust=0.5))+
  ylab("")+xlab("")+ggtitle("Etiocholanolone") +
  scale_x_discrete(labels=c("Control_Lipid" = "Control \n lipid", 
                            "Control_Protein" = "Control \n protein", 
                            "MG_Lipid" = "MG \n lipid", 
                            "MG_Protein" = "MG \n protein"))
andro.panel <- ggarrange(andro.fig,testo.fig,etio.fig, ncol=3, nrow=1, 
                         labels=c("A","B","C"), font.label=list(size=25))
andro.panel

########## andro model without controls ##############
andro_disaggregated.MG <- subset(andro_disaggregated, Treatment=="MG")
andro.model.MG <- lme(response ~ whichpart - 1 + 
                        whichpart:(Yolk.Mass+Egg.+
                                     Fat.y+csize+load.AUC*Diet),  
                      random = ~ whichpart -1 |Bird.ID, 
                      weights=varIdent(form=~1|whichpart),
                      corr=corSymm(form=~1|Bird.ID/Tube.number),
                      data=andro_disaggregated.MG, method="ML",
                      control = list(maxIter=500, msMaxIter=500, tolerance=1e-8,
                                     niterEM=250))
AIC(andro.model.MG)
summary(andro.model.MG)
VarCorr(andro.model.MG)

#plotting
colnames(andro_disaggregated.MG)
levels(andro_disaggregated.MG$whichpart) <- c("Androstenedione", 
                                           "Testosterone", 
                                           "Etiocholanolone")
protein.andro.mg <- subset(andro_disaggregated.MG, Diet=="Protein")
lipid.andro.mg <- subset(andro_disaggregated.MG, Diet=="Lipid")

PAML <- ggplot(protein.andro.mg, aes(x=load.AUC, y=response, color=whichpart)) +
  geom_smooth(method=lm) +
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none",
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        plot.title=element_text(size=30, hjust=0.5))+
  scale_color_paletteer_d("nationalparkcolors::CraterLake")+
  ylab("Scaled concentration in eggs")+xlab("")+ggtitle("Protein-MG")
LAML <- ggplot(lipid.andro.mg, aes(x=load.AUC, y=response, color=whichpart)) +
  geom_smooth(method=lm) +
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position=c(0.3,0.8),
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        plot.title=element_text(size=30, hjust=0.5))+
  scale_color_paletteer_d("nationalparkcolors::CraterLake")+
  ylab("")+xlab("")+ggtitle("Lipid-MG")
PAME <- ggplot(protein.andro.mg, aes(x=ES.AUC, y=response, color=whichpart)) +
  geom_smooth(method=lm) +
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="right",
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        plot.title=element_text(size=30, hjust=0.5))+
  ylab("Scaled concentration in eggs")+xlab("Eye score AUC")+ggtitle("Protein-MG")
LAME <- ggplot(lipid.andro.mg, aes(x=ES.AUC, y=response, color=whichpart)) +
  geom_smooth(method=lm) +
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="right",
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        plot.title=element_text(size=30, hjust=0.5))+
  ylab("Scaled concentration in eggs")+xlab("Eye score AUC")+ggtitle("Lipid-MG")

andro.panel2 <- ggarrange(PAML,LAML,PAME,LAME, ncol=2, nrow=2, 
                          labels=c("A","B","C","D"), font.label=list(size=25)
)
andro.panel2

andro.panel3 <- ggarrange(PAML,LAML, ncol=2, nrow=1, 
                          labels=c("A","B"), font.label=list(size=25)
)
andro.panel3

########## progestagens #########################
part1 <- cbind(horm3, response = (horm3$Progesterone), whichpart = factor(1))
part2 <- cbind(horm3, response = (horm3$Pregnenolone), whichpart = factor(2))
part3 <- cbind(horm3, response = (horm3$X17a.hydroxypregnenolone), whichpart = factor(3))
part4 <- cbind(horm3, response = (horm3$X5_.dihydroproggesterone), whichpart = factor(4))
part5 <- cbind(horm3, response = (horm3$Pregnanolone), whichpart = factor(5))

prog_disaggregated <- rbind(part1, part2, part3, part4, part5)
prog_disaggregated <- 
  prog_disaggregated[order(prog_disaggregated$Bird.ID),]
names(prog_disaggregated)
rownames(prog_disaggregated) <- 1:dim(prog_disaggregated)[1]
# A short look at what has been produced:
prog_disaggregated[1:10,]
# Omit rows with missings in response variable: there are none
prog_disaggregated <- prog_disaggregated[!is.na(prog_disaggregated$response),]
# How much is left:
dim(prog_disaggregated)

prog.model <- lme(response ~ whichpart - 1 + 
                    whichpart:(Diet*Treatment+Egg.+BCI+Yolk.Mass
                               +csize),  
                  random = ~ whichpart -1 |Bird.ID, 
                  weights=varIdent(form=~1|whichpart),
                  corr=corSymm(form=~1|Bird.ID/Tube.number),
                  data=prog_disaggregated, method="ML",
                  control = list(maxIter=800, msMaxIter=800, tolerance=1e-8,
                                 niterEM=250))
AIC(prog.model)
summary(prog.model)
VarCorr(prog.model)
#extract model tables
prog.diet <- emmeans(prog.model, ~ Diet*Treatment , 
                     at = list(whichpart = c("1","2","3", "4","5")))
prog.table <- as.data.frame(prog.diet)
prog.table
write.csv(prog.table, file="prog_table.csv")
#tukey
prog.diet1 <- as.data.frame(pairs(
  emmeans(prog.model, ~ Diet*Treatment,
          at = list(whichpart = c("1"))),adjust="tukey"))
prog.diet2 <- as.data.frame(pairs(
  emmeans(prog.model, ~ Diet*Treatment,
          at = list(whichpart = c("2"))),adjust="tukey"))
prog.diet3 <- as.data.frame(pairs(
  emmeans(prog.model, ~ Diet*Treatment,
          at = list(whichpart = c("3"))),adjust="tukey"))
prog.diet4 <- as.data.frame(pairs(
  emmeans(prog.model, ~ Diet*Treatment,
          at = list(whichpart = c("4"))),adjust="tukey"))
prog.diet5 <- as.data.frame(pairs(
  emmeans(prog.model, ~ Diet*Treatment,
          at = list(whichpart = c("5"))),adjust="tukey"))

write.csv(prog.diet1, file="prog.diet1.tukey.csv")
write.csv(prog.diet2, file="prog.diet2.tukey.csv")
write.csv(prog.diet3, file="prog.diet3.tukey.csv")
write.csv(prog.diet4, file="prog.diet4.tukey.csv")
write.csv(prog.diet5, file="prog.diet5.tukey.csv")
#plot
prog.table <- prog.table %>% mutate(MGD = paste(Treatment, Diet, sep="_"))
prog.table.proges <- subset(prog.table, whichpart=="1")
prog.table.pregne <- subset(prog.table, whichpart=="2")
prog.table.hp <- subset(prog.table, whichpart=="3")
prog.table.dhp <- subset(prog.table, whichpart=="4")
prog.table.pregna <- subset(prog.table, whichpart=="5")


proges.fig <- ggplot(prog.table.proges, aes(x=MGD, y=emmean, color=Diet)) +
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL)) +
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none",
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        plot.title=element_text(size=20, hjust=0.5))+
  ylab("Effect size")+xlab("")+ggtitle("Progesterone")+
  scale_x_discrete(labels=c("Control_Lipid" = "Control \n lipid", 
                            "Control_Protein" = "Control \n protein", 
                            "MG_Lipid" = "MG \n lipid", 
                            "MG_Protein" = "MG \n protein"))
pregne.fig <- ggplot(prog.table.pregne, aes(x=MGD, y=emmean, color=Diet)) +
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL)) +
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none",
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        plot.title=element_text(size=20, hjust=0.5))+
  ylab("")+xlab("")+ggtitle("Pregnenolone")+
  scale_x_discrete(labels=c("Control_Lipid" = "Control \n lipid", 
                            "Control_Protein" = "Control \n protein", 
                            "MG_Lipid" = "MG \n lipid", 
                            "MG_Protein" = "MG \n protein"))
hp.fig <- ggplot(prog.table.hp, aes(x=MGD, y=emmean, color=Diet)) +
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL)) +
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none",
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        plot.title=element_text(size=17, hjust=0.5))+
  ylab("")+xlab("")+ggtitle("17a-hydroxypregnenolone")+
  scale_x_discrete(labels=c("Control_Lipid" = "Control \n lipid", 
                            "Control_Protein" = "Control \n protein", 
                            "MG_Lipid" = "MG \n lipid", 
                            "MG_Protein" = "MG \n protein"))
dhp.fig <- ggplot(prog.table.dhp, aes(x=MGD, y=emmean, color=Diet)) +
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL)) +
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none",
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        plot.title=element_text(size=20, hjust=0.5))+
  ylab("Effect size")+xlab("")+ggtitle("Pregnanedione")+
  scale_x_discrete(labels=c("Control_Lipid" = "Control \n lipid", 
                            "Control_Protein" = "Control \n protein", 
                            "MG_Lipid" = "MG \n lipid", 
                            "MG_Protein" = "MG \n protein"))
pregna.fig <- ggplot(prog.table.pregna, aes(x=MGD, y=emmean, color=Diet)) +
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL)) +
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none",
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        plot.title=element_text(size=20, hjust=0.5))+
  ylab("")+xlab("")+ggtitle("Pregnanolone")+
  scale_x_discrete(labels=c("Control_Lipid" = "Control \n lipid", 
                            "Control_Protein" = "Control \n protein", 
                            "MG_Lipid" = "MG \n lipid", 
                            "MG_Protein" = "MG \n protein"))
prog.panel <- ggarrange(proges.fig,pregne.fig,hp.fig,dhp.fig,pregna.fig,
                        ncol=3, nrow=2, 
                        labels=c("","","","",""), font.label=list(size=25))
prog.panel

########### Make Figure 3 ######################
Figure3 <- ggarrange(proges.fig,pregne.fig,hp.fig,dhp.fig,pregna.fig,
                     andro.fig,testo.fig,etio.fig,
                     ncol=3, nrow=3, 
                     labels=c("A","B","C","D","E","F","G","H"), 
                     font.label=list(size=25))
Figure3 

############ Figure 4 #####################
colnames(andro_disaggregated)
levels(andro_disaggregated$whichpart) <- c("Androstenedione", 
                                           "Testosterone", 
                                           "Etiocholanolone")
fat.fig <- ggplot(andro_disaggregated, aes(x=Fat.y, y=response, color=whichpart)) +
  #geom_point() +
  geom_smooth(method="lm")+
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position=c(0.7,0.9),
        legend.text=element_text(size=15),
        legend.title=element_blank(),
        plot.title=element_text(size=20, hjust=0.5))+
  scale_color_paletteer_d("nationalparkcolors::CraterLake")+
  ylab("Scaled concentration")+xlab("Maternal fat score")+ggtitle("Androgens")
fat.fig

colnames(prog_disaggregated)
levels(prog_disaggregated$whichpart) <- c("Progesterone", 
                                          "Pregnenolone",
                                          "17a.hydroxypregnenolon",
                                          "Pregnanedione",
                                          "Pregnanolone")

mass.fig <- ggplot(prog_disaggregated, aes(x=BCI, y=response, color=whichpart)) +
  #geom_point() +
  geom_smooth(method="lm")+
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position=c(0.31,0.87),
        legend.text=element_text(size=15),
        legend.title=element_blank(),
        plot.title=element_text(size=20, hjust=0.5))+
  scale_color_paletteer_d("nationalparkcolors::Saguaro")+
  ylab("Scaled concentration")+xlab("Body condition index")+ggtitle("Progestagens")
mass.fig

Figure4 <- ggarrange(mass.fig,fat.fig, ncol=2, nrow=1, 
                     labels=c("A","B"), font.label=list(size=20))
Figure4 
############ MG only prog models ###########################

prog_disaggregated.MG <- subset(prog_disaggregated, Treatment == "MG")
prog.model.MG <- lme(response ~ whichpart - 1 +
                       whichpart:(Yolk.Mass+Diet*load.AUC+csize+scale(BCI)) ,  
                     random = ~ whichpart -1 |Bird.ID, 
                     weights=varIdent(form=~1|whichpart),
                     corr=corSymm(form=~1|Bird.ID/Tube.number),
                     data=prog_disaggregated.MG, method="ML",
                     control = list(maxIter=800, msMaxIter=800, tolerance=1e-8,
                                    niterEM=250))
AIC(prog.model.MG)
summary(prog.model.MG)
VarCorr(prog.model.MG)

levels(prog_disaggregated.MG$whichpart) <- c("Progesterone", 
                                          "Pregnenolone",
                                          "17a.hydroxypregnenolon",
                                          "Pregnanedione",
                                          "Pregnanolone")
protein.prog.mg <- subset(prog_disaggregated.MG, Diet=="Protein")
lipid.prog.mg <- subset(prog_disaggregated.MG, Diet=="Lipid")
PPML <- ggplot(protein.prog.mg, aes(x=load.AUC, y=response, color=whichpart)) +
  geom_smooth(method=lm) +
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none",
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        plot.title=element_text(size=30, hjust=0.5))+
  scale_color_paletteer_d("nationalparkcolors::Saguaro")+
  ylab("")+xlab("")+ggtitle("Protein-MG")
LPML <- ggplot(lipid.prog.mg, aes(x=load.AUC, y=response, color=whichpart)) +
  geom_smooth(method=lm) +
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position=c(0.32,0.87),
        legend.text=element_text(size=15),
        legend.title=element_blank(),
        plot.title=element_text(size=30, hjust=0.5))+
  scale_color_paletteer_d("nationalparkcolors::Saguaro")+
  ylab("Scaled concentration in eggs")+xlab("")+ggtitle("Lipid-MG")
prog.panel.MG <- ggarrange(PPML,LPML,
                           ncol=2, nrow=1, 
                           labels=c("",""), font.label=list(size=25))
prog.panel.MG

############ Figure 5 ########################
PAML <- ggplot(protein.andro.mg, aes(x=load.AUC, y=response, color=whichpart)) +
  geom_smooth(method=lm) +
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none",
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        plot.title=element_text(size=30, hjust=0.5))+
  scale_color_paletteer_d("nationalparkcolors::CraterLake")+
  ylab("")+xlab("MG load AUC")+ggtitle("Protein-MG")
LAML <- ggplot(lipid.andro.mg, aes(x=load.AUC, y=response, color=whichpart)) +
  geom_smooth(method=lm) +
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position=c(0.22,0.87),
        legend.text=element_text(size=15),
        legend.title=element_blank(),
        plot.title=element_text(size=30, hjust=0.5))+
  scale_color_paletteer_d("nationalparkcolors::CraterLake")+
  ylab("Scaled concentration in eggs")+xlab("MG load AUC")+ggtitle("Lipid-MG")

Figure5 <- ggarrange(LPML,PPML,LAML,PAML, ncol=2, nrow=2, 
                           labels=c("A","B","C","D"), 
                     font.label=list(size=25))
Figure5 
############# terminal investment analysis #############
hist(horm3$Yolk.Mass)
hist(horm3$Egg.Mass.Day.Laid)

massm1 <- lmer(Yolk.Mass ~ Diet * Treatment + 
                 csize+ (1|Bird.ID), data=horm3)
AIC(massm1)
summary(massm1)
Anova(massm1)
#save tables
massm1S <- summary(massm1)
str(massm1S$coefficients)
massm1S <- (massm1S$coefficients)
massm1S<- as.matrix(massm1S)
write.csv(massm1S, file = "massm1S_lmer.csv")
massm1A<- as.matrix(Anova(massm1))
colnames(massm1A)<-c("Chi squared", "df", "p value")
write.csv(massm1A, file = "massm1A_ANOVA.csv")

#egg mass model
massm2 <- lmer(Egg.Mass.Day.Laid ~ Diet * Treatment + 
                 csize + (1|Bird.ID), data=horm3)
AIC(massm2)
summary(massm2)
Anova(massm2)
#save tables
massm2S <- summary(massm2)
str(massm2S$coefficients)
massm2S <- (massm2S$coefficients)
massm2S<- as.matrix(massm2S)
write.csv(massm2S, file = "massm2S_lmer.csv")
massm2A<- as.matrix(Anova(massm2))
colnames(massm2A)<-c("Chi squared", "df", "p value")
write.csv(massm2A, file = "massm2A_ANOVA.csv")

#models using clutch-level data
str(moms)
colnames(horm3)
cs <- horm3[,c(1,11, 12,29)]
CLD <- merge(moms, cs, by="Bird.ID")


#clutch size model - 
csizem <- glm(csize ~ Diet * Treatment+Fat.y, data=CLD)
AIC(csizem)
summary(csizem)
Anova(csizem)
#save tables
csizemS <- summary(csizem)
str(csizemS$coefficients)
csizemS <- (csizemS$coefficients)
csizemS<- as.matrix(csizemS)
write.csv(csizemS, file = "csizemS_lmer.csv")
csizemA<- as.matrix(Anova(csizem))
colnames(csizemA)<-c("Chi squared", "df", "p value")
write.csv(csizemA, file = "csizemA_ANOVA.csv")

csizep <- ggplot(CLD, aes(x=Diet, y=csize, color=Treatment)) +
  geom_violin()+
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position=c(0.7,0.9),
        legend.text=element_text(size=15),
        legend.title=element_blank(),
        plot.title=element_text(size=20, hjust=0.5))+
  ylab("# of eggs")+xlab("Treatment")+ggtitle("Clutch size")
csizep


#total clutch mass
tmassm <- glm(clutch.mass ~ Diet * Treatment +Fat.y, data=CLD)
AIC(tmassm)
summary(tmassm)
Anova(tmassm)
#save tables
tmassmS <- summary(tmassm)
str(tmassmS$coefficients)
tmassmS <- (tmassmS$coefficients)
tmassmS<- as.matrix(tmassmS)
write.csv(tmassmS, file = "tmassmS_lmer.csv")
tmassmA<- as.matrix(Anova(tmassm))
colnames(tmassmA)<-c("Chi squared", "df", "p value")
write.csv(tmassmA, file = "tmassmA_ANOVA.csv")

cmassp <- ggplot(CLD, aes(x=Diet, y=clutch.mass, color=Treatment)) +
  geom_violin()+
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position=c(0.7,0.9),
        legend.text=element_text(size=15),
        legend.title=element_blank(),
        plot.title=element_text(size=20, hjust=0.5))+
  ylab("mass (g)")+xlab("Treatment")+ggtitle("Clutch mass")
cmassp


#total clutch mass
FIm <- glm(FI ~ Diet * Treatment +Fat.y, data=CLD)
AIC(FIm)
summary(FIm)
Anova(FIm)
#save tables
tmassmS <- summary(tmassm)
str(tmassmS$coefficients)
tmassmS <- (tmassmS$coefficients)
tmassmS<- as.matrix(tmassmS)
write.csv(tmassmS, file = "tmassmS_lmer.csv")
tmassmA<- as.matrix(Anova(tmassm))
colnames(tmassmA)<-c("Chi squared", "df", "p value")
write.csv(tmassmA, file = "tmassmA_ANOVA.csv")



############# egg antibody anaylsis ######################
eggAB <- read.csv("eggAB.csv")
str(eggAB)
hist(eggAB$ng_ml)

ABm1 <- lmer(ng_ml ~ Treatment*Diet+Egg.Number+
                  (1|Female.ID) + (1|plate), data=eggAB)
summary(ABm1)
Anova(ABm1)


IgY1 <- ggplot(eggAB, aes(x=TD, y=ng_ml, color=Diet)) +
  geom_boxplot() +
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none",
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        plot.title=element_text(size=30, hjust=0.5))+
  ylab("Antibody concentration (ng/mL)")+xlab("")
IgY2 <- ggplot(eggAB, aes(x=Egg.Number, y=ng_ml)) +
  geom_smooth(method="lm") +
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="none",
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        plot.title=element_text(size=30, hjust=0.5))+
  ylab("")+xlab("Egg lay order")

Figure1 <- ggarrange(IgY1,IgY2, ncol=2, nrow=1, 
                     labels=c("A","B"), 
                     font.label=list(size=25))
Figure1

