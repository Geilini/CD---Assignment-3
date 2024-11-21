#######################
# Exercise 1
#######################
library(lattice)
library(lmerTest)
library(MASS)
library(xtable)
library(stringr) # works like f_strings in python
library(emmeans)
library(multcomp)

getwd()
remove(list=ls())
par(mfrow=c(1,1))
options(digits=4)

#Data:
(Cr <-rep(c(1,2,3),each=6))
(Fam<- seq(1,18))
(Fa<-c(5,11,3,3,10,9, 11,10,5,1,5,7, 6,9,2,0,3,6))
(Mo<-c(7,8,12,19,9,0, 7,5,4,9,5,3, 3,6,2,2,2,2))
(Ch1<-c(6,11,19,12,15,6, 7,8,3,4,10,13, 5,6,6,10,0,4))
(Ch2<-c(25,33,6,17,11,9, 15,13,18,16,16,17, 7,14,15,16,3,7))
(Ch3<-c(19,35,21,17,17,5, 13,17,10,8,20,18, 3,10,8,21,14,20))

df_wide<-data.frame(crowding=factor(Cr), family=factor(Fam),
            father=Fa, mother=Mo, ch1=Ch1, ch2=Ch2, ch3=Ch3)
str(df_wide)

df_wide

levels(df_wide$crowding)<-c("Over", "Crowded", "Under")
summary(df_wide)
apply(df_wide[,3:7], 2, sum)
apply(df_wide[,3:7], 1, sum)

matplot(df_wide[,3:7], type="b", xlab="Families")

matplot(t(df_wide[,3:7]), type="b",xlab="Family status")


?reshape

(aux<-reshape(df_wide, direction="long", 
              idvar="family",varying = list(3:7),
              v.names = "swabs",
              timevar="status"))

(df<- aux[order(aux$family),1:4])

df$family<-as.factor(df$family)
df$status<-as.factor(df$status)
levels(df$status)
df$status<- fct_recode(df$status, "father"="1", "mother"="2",
                            "ch1"="3", "ch2"="4","ch3"="5")
str(df)



################ Exploratory Analysis ################ 
# baseline anova
lm_init = lm(swabs ~ crowding*status,data=df)
anova(lm_init)


# boxplots

boxplot(swabs ~ crowding, data = df, 
        xlab = " ", ylab = "Positive Swabs", cex.lab = 1.2)
means <- tapply(df$swabs, df$crowding, mean, na.rm = TRUE)
points(x = 1:length(means), y = means, col = "red", pch = 19, cex = 1.2)

boxplot(swabs ~ status, data = df, 
        xlab = " ", ylab = "Positive Swabs", cex.lab = 1.2)
means <- tapply(df$swabs, df$status, mean, na.rm = TRUE)
points(x = 1:length(means), y = means, col = "red", pch = 19, cex = 1.2)


# Interaction plots
matplot(df_wide[,3:7], type="b", xlab="Families", ylab = "Positive Swabs", 
        main = "",  cex.lab = 1.2)
#matplot(t(df_wide[,3:7]), type="b",xlab="Family status", ylab = "Positive Swabs")

interaction.plot(df$status, df$crowding, df$swabs,
                 col = c("red", "green", "blue"), lty = 1,
                 xlab = " ", ylab = "Mean Positive Swabs", trace.label = "Crowding",
                 cex.lab = 1.2)

interaction.plot(df$crowding, df$status, df$swabs,
                 col = c("red", "green", "blue", "orange", "black"), lty = 1,
                 xlab = "", ylab = "Mean Positive Swabs", trace.label = "Family Status",
                 cex.lab = 1.2)



################ LMM's ################

# Initial model
m0 = lmer(swabs ~ crowding*status + (1|family),data=df,REML=F)
ranova(m0)
drop1(m0)

m1 = update(m0,~. - crowding:status)
drop1(m1) # cannot reduce further
summary(m1)

# Final model fitted using REML
final_model = lmer(swabs ~crowding + status + (1|family),data=df,REML=T)
tab1 = xtable(cbind(coefficients(summary(final_model)),confint(final_model)[3:9,1:2]))
print(tab1, type = "latex", digits = 3)

tab2= xtable(confint(final_model))
print(tab2, type = "latex", digits = 3)

################ Post-hoc Analysis ################
# Post-hoc
print(xtable(emmeans(final_model,pairwise~status)$contrasts),type='latex',digits=3)
print(xtable(emmeans(final_model,pairwise~crowding)$contrasts),type='latex',digits=3)


################ Model Diagnostics ################

# Cooks distance
png(str_glue('diagnostics_plot_cooks.png'), width = 800, height = 600)
par(mfrow=c(1,1),mai=c(1,1,1,1))
infl <- influence(final_model, obs = TRUE)
plot(cooks.distance(infl),ylab="Cook's distance",xlab='Observations',lty=1,type="h",cex=1.5,cex.main=1.5,cex.axis=1.5,lwd=3,cex.lab=1.5,main="Cook's Distance")
dev.off()

# Residuals vs. fitted
png(str_glue('diagnostics_plot_residual_fitted.png'), width = 800, height = 600)
par(mfrow=c(1,1),mai=c(1,1,1,1))
plot(fitted(final_model),residuals(final_model),cex=1.5,cex.main=1.5,cex.axis=1.5,lwd=3,cex.lab=1.5,xlab='Fitted Values',ylab='Residuals',main='Residuals vs Fitted')
abline(h=0,col='red')
dev.off()

png(str_glue('diagnostics_plot_qq_plot_residuals.png'))
par(mfrow=c(1,1),mai=c(1,1,1,1))
standardized_residuals = (residuals(final_model) - mean(residuals(final_model)))/(sd(residuals(final_model)))
qqnorm(standardized_residuals,cex=1.5,cex.main=1.5,cex.axis=1.5,lwd=3,cex.lab=1.5)
qqline(standardized_residuals)
dev.off()


# Random qq plots
temp = names(ranef(final_model))
ranefs = ranef(final_model)
for(i in 1:length(names)){
  png(str_glue('diagnostics_ranef_{temp[i]}.png'), width = 800, height = 600)
  par(mfrow=c(1,1),mai=c(1.1,1,1,1))
  qqnorm(unlist(ranefs[[i]]),main=paste(temp[i]),cex=1.5,cex.main=1.5,cex.axis=1.5,lwd=3,cex.lab=1.5)
  lines((-3):3,sd(unlist(ranefs[[i]]))*((-3):3),col="red",lwd=3)
  dev.off()
}
par(mfrow=c(1,1))




