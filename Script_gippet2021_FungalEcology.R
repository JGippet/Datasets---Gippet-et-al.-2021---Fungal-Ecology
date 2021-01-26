###################################################################################################################
##############################○○○○                                               ○○○○##############################
##################○○○○                          Gippet et al. - 2020                         ○○○○##################
###########○○○○                    Laboulbenia formicarum in the upper Rhone Valley                 ○○○○###########
##################○○○○                                                                       ○○○○##################
##############################○○○○                                               ○○○○##############################
###################################################################################################################

# R version 3.6.2 (2019-12-12)
library(Rphylip)
library(pegas)
library(adegenet)
library(poppr)
library(ecodist)
library(vegan)
library(ape)
library(dartR)
library(pvclust)
library(mpmcorrelogram)
library(hierfstat)
library(PopGenReport)
library(resample)

library(ade4)
library(performance)
library(ggplot2)
library(sjPlot)
library(DHARMa)

library(lme4)
library(ggeffects)
library(spatialreg)
library(rlist)
library(spdep)


### Documentation used for analyses

## mantel tests
# from https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecs2.1558:
# In all cases, Mantel tests were conducted by constructing a “predictor matrix,” a pairwise Euclidean distance matrix between the capture locations of all individuals in the data set. A “dependent distance matrix” was created with the binary exposure status of all pairs of individuals in the data set for each pathogen. A pair of individuals that were either both seropositive or both seronegative had a distance of “0,” whereas a pair of individuals that had one seropositive individual and one seronegative individual had a distance of “1.” Distance matrices were calculated and Mantel tests conducted using the stats and ade4 packages in R (R Core Team 2014). Mantel test results include P‐values (P), ±SE, and the observed correlation (ρ).


################### ################### ################### 
################# Regional scale analysis ################### 
################### ################### ################### 

#####  load dataset - epidemio + land cover  #####
data1<-read.table("dataset1_landscape_gippet2020.txt", h=T, row.names=1)
head(data1)
dim(data1)
##################################################

# dataset with colonies with at least 10 workers screened
data2 <- data1[data1$worker_see>9,]
dim(data2)
head(data2)

# dataset with genotyped colonies only
data3 <- data2[data2$genet=="Yes",]
dim(data3)
head(data3)

#####  Load dataset - microsatellite data  #######
gen <- read.gtx("dataset3_microsat_gippet2020.gtx")
##################################################

gen1 <- loci2genind(gen)
unique(gen1@pop)
gen2 <- gen1[gen1@pop %in% data3$nom_genet]
unique(gen2@pop)


## analyses

# Calculate geographic distance between colonies
geoD <- dist(data2[,1:2]) # euclidean geographic distance in meters
geoV <- as.vector(geoD)

# Calculate infection distance between colonies
infD <- dist(data2[,9])
infD <- as.matrix(infD)
row.names(infD)  <- row.names(data2)
colnames(infD)  <- row.names(data2)
infD <- as.dist(infD)
infV <- as.vector(infD)


df1 <- as.data.frame(cbind(c(0,1), 
                           c(mean(geoV[which(infV==0)]), mean(geoV[which(infV==1)])), 
                           c(sd(geoV[which(infV==0)]), sd(geoV[which(infV==1)]))))
colnames(df1) <- c("Infection_status", "mean_Geodist", "sd_Geodist")
df1
r1 <- mantel.rtest(geoD,infD, nrepet=10000)
r1 # obs: 0.04003157, p-value: 0.05819418 
plot(r1)

# p <- ggplot(df1, aes(x=Infection_status, y=mean_Geodist)) + 
#   geom_bar(stat="identity", position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean_Geodist-sd_Geodist, ymax=mean_Geodist+sd_Geodist), width=.2,
#                 position=position_dodge(.9))
# 
# p + scale_fill_brewer(palette="Paired") + theme_minimal()


#### 
# Calculate infection distance
infD3 <- dist(data3[,9])
infD3 <- as.matrix(infD3)
row.names(infD3)  <- data3$nom_genet
colnames(infD3)  <- data3$nom_genet
infD3 <- as.dist(infD3)
infV3 <- as.vector(infD3)


# Calculate genetic distance (FST*) between colonies
FST <- pairwise.fst(gen2, pop = gen2$pop)
FST <- as.matrix(FST)
row.names(FST) <- unique(gen2$pop)
colnames(FST) <- unique(gen2$pop)
FST1 <- FST/(1-FST)
dim(FST1)
colnames(FST1)
FST1 <- as.dist(FST1)
FST1V <- as.vector(FST1)

df2 <- as.data.frame(cbind(c(0,1), 
                           c(mean(FST1V[which(infV3==0)]), mean(FST1V[which(infV3==1)])), 
                           c(sd(FST1V[which(infV3==0)]), sd(FST1V[which(infV3==1)]))))
colnames(df2) <- c("Infection_status", "mean_Genetdist", "sd_Genetdist")
df2
r2 <- mantel.rtest(FST1, infD3, nrepet=10000)
r2 # obs: 0.04260753, p-value: 0.1393861 
plot(r2, nclass=20) # no genetic signal in infection status



# Calculate landcover distance between colonies
head(data2)
envD <- dist(scale(data2[,c(10:16)]), method="euclidean")
envV <- as.vector(envD)

df3 <- as.data.frame(cbind(c(0,1), 
                           c(mean(envV[which(infV==0)]), mean(envV[which(infV==1)])), 
                           c(sd(envV[which(infV==0)]), sd(envV[which(infV==1)]))))
colnames(df3) <- c("Infection_status", "mean_Envdist", "sd_Envdist")
df3
r3 <- mantel.rtest(envD,infD, nrepet=10000)
r3 # obs: 0.08335536, p-value: 0.0059994 
plot(r3)


##############################################
##### Environment-infection relationship #####
##############################################

head(data2)
dim(data2)
range(data2$ESM_vegetation)
mean(data2$ESM_vegetation)
sqrt(var(data2$ESM_vegetation))

pca0 <- dudi.pca(data2[,c(10:16)], scannf = FALSE, nf = 3)
s.corcircle(pca0$co, 1, 2)
data2$ax1 <- pca0$li[,1]
data2$ax2 <- pca0$li[,2]
# data2$ax3 <- pca0$li[,3]
(pca0$eig)/sum(pca0$eig)*100

m4 <- glm(inf1_noinf0  ~ 1 +
            ax1 +
            ax2,
          family=binomial, data=data2)
summary(m4)
drop1(m4, test="Chisq")
m4_2 <- update(m4, .~.- ax2 )
drop1(m4_2, test="LRT")
summary(m4_2)
# no predictors of infection status

r2_nagelkerke(m4_2) # 0.11
plot(simulateResiduals(m4_2))

theme_set(theme_classic())
effects4 <- plot_model(m4_2, type="eff")
testplot4 <- effects4$ax1 + 
  geom_point(data = data2, mapping = 
               aes(x =ax1, y = inf1_noinf0, size = log(worker_see) ) )
testplot4

# pdf("fig3a_v2.pdf")
# testplot4
# dev.off()
# 
# pdf("s.corcircle__fig3a_v2.pdf")
# s.corcircle(pca0$co, 1, 2, grid=F)
# dev.off()


#### only parasited colonies #### 
data4 <- data2[data2$inf1_noinf0==1,]
dim(data4)
head(data4)

pca0b <- dudi.pca(data4[,c(10:16)], scannf = FALSE, nf = 3)
(pca0b$eig)/sum(pca0b$eig)*100
s.corcircle(pca0b$co, 1, 2, grid=F)
data4$ax1 <- pca0b$li[,1]
data4$ax2 <- pca0b$li[,2]
# data4$ax3 <- pca0b$li[,3]


m5 <- glm(prevalence ~ 1 + # cbind(worker_inf, worker_see-worker_inf)
            ax1 +
            ax2,
          weights=log(worker_see),
          family=quasibinomial, 
          data=data4) #
summary(m5)
drop1(m5, test="LRT")
m5_2 <- update(m5, .~. -ax1)
drop1(m5_2, test="Chisq")
summary(m5_2)

par(mfrow=c(2,2))
plot(m5)
shapiro.test(m5$residuals) # W = 0.96186, p-value = 0.2181
hist(m5$residuals)

r2_nagelkerke(m5) # 0.55

theme_set(theme_classic())
effects1 <- plot_model(m5, type="eff")
testplot0 <- effects1$ax2 + 
  geom_point(data = data4, mapping = 
               aes(x = (ax2), y = prevalence, size = log(worker_see) ) )
testplot0


# pdf("fig3b_v2.pdf")
# testplot0
# dev.off()
# 
# pdf("s.corcircle__fig3b_v2.pdf")
# s.corcircle(pca0b$co, 1, 2, grid=F)
# dev.off()



################### ################### ################### 
################### Local scale analysis ################### 
################### ################### ################### 

#####  Load dataset  #####
dataLocal1 <-read.table("dataset2_local_gippet2020.txt", h=T, row.names=1)
head(dataLocal1)
dim(dataLocal1)
###################################### 

dataLocal1$inf1_noninf0 <- rep(0, dim(dataLocal1)[1])
dataLocal1$inf1_noninf0[which(dataLocal1$prevalence>0)] <- 1

# dataLocal1$regional <- rep(NA, dim(dataLocal1)[1])
# for (i in dataLocal1$colony){
#   if(dim(data4[which(row.names(data4)==i),])[1]>0){
#     dataLocal1$regional[dataLocal1$colony==i] <- data4$ax1[which(row.names(data4)==i)]
#   }
#   
# }


# keeping only colonies with at least one infected worker
dataLocal2 <- dataLocal1[which(dataLocal1$colony %in% row.names(data4)),]
# keeping only infected colonies with at least 5 samples
dataLocal3 <- dataLocal2[dataLocal2$colony %in%  names(which(table(dataLocal2$colony)>=5)),]
dim(dataLocal3) # n= 219
unique(dataLocal3$colony)
sum(dataLocal3$workers_see)
var(sort(table(dataLocal3$colony))[18:33])

range(dataLocal3$unsealed_ways)
mean(dataLocal3$unsealed_ways)
sqrt(var(dataLocal3$unsealed_ways))

range(dataLocal3$impervious)
mean(dataLocal3$impervious)
sqrt(var(dataLocal3$impervious))

head(dataLocal3)

pca1b <- dudi.pca(dataLocal3[,c(6,7,8,9)], scannf = FALSE, nf = 3)
s.corcircle((-1) * pca1b$co, 1, 2, grid=F)
dataLocal3$ax1 <- (-1) * pca1b$li[,1]
dataLocal3$ax2 <- (-1) * pca1b$li[,2]
pca1b$eig/sum(pca1b$eig)*100


sum(table(dataLocal1$colony))
sum(table(dataLocal3$colony))


hist(dataLocal3$prevalence)

###############################################################################################
###############################################################################################
glm1 <- glmer(cbind(workers_inf, workers_see-workers_inf) ~ # 
                ax1 + 
                ax2 + 
                (1|colony),
              family=binomial, 
              data=dataLocal3)

summary(glm1)
drop1(glm1, test="Chisq")

plot(glm1)
plot(simulateResiduals(glm1))
r2(glm1)

eff2 <- ggpredict(glm1, c("ax2 [all]"), type = "fe")
plot(eff2, colors="set1",line.size=1.25) # , limits=c(0,11)

theme_set(theme_classic())
effects5 <- plot_model(glm1, type="eff")
testplot5 <- effects5$ax2 + 
  geom_point(data = dataLocal3, mapping = 
               aes(x = (ax2), y = prevalence, size = log(workers_see) ) )
testplot5

###############################################################################################
###############################################################################################




