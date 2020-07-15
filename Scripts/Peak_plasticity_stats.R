#This is the master code for the paper 
#"Refining the timing of recombination rate plasticity in response to temperature in Drosophila pseudoobscura"
#written by Ulku Huma Altindag

#It is split into 5 code sections which corresponds to the each sequential experiment in this study
#These sections can be expanded/collapsed in RStudio 
#Expand: either click on the arrow in the gutter or on the icon that overlays the folded code 
#Collapse: click on the arrow in the gutter

#necessary packages are listed below:

library(ggplot2)
library(ggthemes)
library(emmeans)
library(lme4)
library(lmerTest)
library(doBy)
library(reshape2)
library(car)

#Only pre-requisite is to have all the datasets listed in the github page for the paper
#in the same folder as the code in.


# Experiment1 -------------------------------------------------------------

#load the datasets
exp1=read.csv("Exp2_rawdata.csv", header= TRUE)
exp1_backcross=read.csv("exp2_backcross.csv", header=T)


#The 24 hour transfers are aggregated into the same with experiment 3 due to insufficient sample sizes.
exp1$new_day=ifelse(exp1$Day=="A","A",
                  ifelse(exp1$Day..letter.of.vial.=="B","A",
                         ifelse(exp1$Day..letter.of.vial.=="C","B",
                                ifelse(exp1$Day..letter.of.vial.=="D","B",
                                       ifelse(exp1$Day..letter.of.vial.=="E","C",
                                              ifelse(exp1$Day..letter.of.vial.=="F","C",
                                                     ifelse(exp1$Day..letter.of.vial.=="G","D",
                                                            ifelse(exp1$Day..letter.of.vial.=="H","D",
                                                                   ifelse(exp1$Day..letter.of.vial.=="I","D",
                                                                          ifelse(exp1$Day..letter.of.vial.=="J","E",
                                                                                 ifelse(exp1$Day..letter.of.vial.=="K","E",
                                                                                        ifelse(exp1$Day..letter.of.vial.=="L","E",
                                                                                               ifelse(exp1$Day..letter.of.vial.=="M","F",
                                                                                                      ifelse(exp1$Day..letter.of.vial.=="N","F",
                                                                                                             ifelse(exp1$Day..letter.of.vial.=="O","F",NA)))))))))))))))
#backcross$Treatment=as.numeric(as.character(backcross$Treatment)) #to make sure R reads it as a character rather than a number.

exp1$Wildtype=as.numeric(as.character(exp1$Wildtype))
exp1$Vellow.vermillion=as.numeric(as.character(exp1$Vellow.vermillion))
exp1$Yellow.only=as.numeric(as.character(exp1$Yellow.only))
exp1$Vermillion.only=as.numeric(as.character(exp1$Vermillion.only))
exp1$Day..letter.of.vial.=as.character(exp1$Day..letter.of.vial.)

#CO groups were defined

exp1$SCO1=exp1$Yellow.only
exp1$SCO2=exp1$Vermillion.only
exp1$NCO1=exp1$Wildtype
exp1$NCO2=exp1$Vellow.vermillion

sco_count=sum(exp1$SCO1, na.rm = TRUE)+sum(exp1$SCO2, na.rm = TRUE)
nco_count=sum(exp1$NCO1+exp1$NCO2, na.rm = TRUE)
num_samples=sum(nco_count+sco_count, na.rm = TRUE)

#rate between yellow and vermillion

(sum(exp1$SCO1, na.rm = TRUE)+sum(exp1$SCO2, na.rm = TRUE))/num_samples


#merge with treatment data
colnames(exp1_backcross)
exp1_merged <- merge(exp1, exp1_backcross, by.x = "Vial..", by.y = "?..Vial.Number", all=T)
exp1_merged = na.omit(exp1_merged)
exp1_merged=na.omit(exp1_merged)

dataset=summaryBy(SCO1+SCO2+NCO1+NCO2~F1.Vial+new_day+Treatment,data=exp1_merged, FUN=sum,na.rm=T)

#add in a column for total offspring

dataset$total_offspring=dataset$SCO1.sum + dataset$SCO2.sum + dataset$NCO1.sum + dataset$NCO2.sum

num_moms=vector(mode="numeric",length=length(dataset$F1.Vial))

#loop through dataset to get the count
for (h in 1:length(dataset$F1.Vial)) { 
  f1_vial=dataset$F1.Vial[h]
  #f1_vial=as.numeric(gsub("[^0-9\\.]","",f1_vial))
  #print(f1_vial)
  mom_ct=length(unique(sort(subset(exp1_merged,exp1_merged$F1.Vial==f1_vial)$Vial..)))
  #store result in the vector
  num_moms[h]=mom_ct
}
dataset$Num_moms=num_moms

dataset$fecundity=dataset$total_offspring/dataset$Num_moms

write.csv(dataset,file="Experiment1_data.csv")


#Statistics of the vy interval
#dataset2=read.csv("Experiment2_data.csv",header=T,stringsAsFactors = F)
dataset$Treatment=as.character(dataset$Treatment)

#get mean fecundity
tapply(dataset$fecundity,dataset$Treatment,mean)
tapply(dataset$fecundity,dataset$new_day,mean)

#poisson regression, similar to a t-test for count data
fit=glm(fecundity~Treatment*new_day,data=dataset,family=quasipoisson)
summary(fit)
#Anova
Fecundity_anova=anova(fit,test="Chisq")
write.csv(Fecundity_anova,"Experiment1_fecundity_anova.csv")
#post-hoc
fit_contrast <- emmeans::emmeans(fit, "Treatment", by="new_day", mode="kenward-roger")
fit_contr <- contrast(fit_contrast, method="trt.vs.ctrl")

pheno_contr <- as.data.frame(summary(fit_contr))
pheno_contr

write.csv(pheno_contr,"Experiment1_fecundity_posthoc_table.csv")

#convert p-values to stars for plot
sig=ifelse(pheno_contr$p.value<0.001,"***",ifelse(pheno_contr$p.value<0.01,"**",ifelse(pheno_contr$p.value<0.05,"*","")))
pdf("Experiment1_fecundity.pdf")

Fecund_figure=ggplot(aes(y=fecundity,x=new_day, col=Treatment,label=Num_moms),data=dataset)+ylab("# Progeny per mom")+ggtitle("Total Fecundity vs. Days post-mating")+theme_base()

Fecund_figure=Fecund_figure+scale_colour_manual(values=c("blue", "red"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("1-2","3-4","5-6","7-9","10-12","13-15"))
Fecund_figure=Fecund_figure+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)+
  annotate(geom="text", x=1, y=25, label=sig[1],size=10)+annotate(geom="text", x=2, y=20, label=sig[2],size=10)+annotate(geom="text", x=3, y=25, label=sig[3],size=10)+annotate(geom="text", x=4, y=25, label=sig[4],size=10)
Fecund_figure

dev.off()

#Now do the same, but with recombination rate

#The progeny in vials lower than 10 is removed from the recombination analysis
dataset=dataset[dataset$total_offspring>=10,]

#sum of crossovers in intervals 1
dataset$num_CO=dataset$SCO1.sum+dataset$SCO2.sum
#sum of non-crossovers in intervals 1
dataset$num_NCO=dataset$NCO1.sum+dataset$NCO2.sum
#total recombination rate
dataset$rec_rate_total=dataset$num_CO/dataset$total_offspring

#get mean recombination
tapply(dataset$rec_rate_total,dataset$Treatment,mean)
tapply(dataset$rec_rate_total,dataset$new_day,mean)
dataset$F1.Vial=as.numeric(gsub("[^0-9\\.]","",dataset$F1.Vial))
#print(dataset$F1.Vial)
fit_recrate=glmer(cbind(num_CO,num_NCO)~(1|F1.Vial)+Treatment*new_day,data=dataset,family=binomial(link="logit"))
summary(fit_recrate)
#Anova
Recratemodel=Anova(fit_recrate,test="Chisq")
write.csv(Recratemodel,"Experiment1_recrate_anova.csv")
#posthoc
fit_contrast_rec <- emmeans::emmeans(fit_recrate, "Treatment", by="new_day", mode="kenward-roger")
fit_contr_rec <- contrast(fit_contrast_rec, method="trt.vs.ctrl")

pheno_contr_rec <- as.data.frame(summary(fit_contr_rec))
pheno_contr_rec
write.csv(pheno_contr_rec,"Experiment1_recrate_posthoc.csv")

#add in odds ratio and standard error
vy_or=exp(pheno_contr_rec$estimate)
vy_error=(pheno_contr_rec$SE)
x=cbind(vy_or,vy_error)
day=c("a","b","c","d","e","f")
x=cbind(day,x)

#convert p-values to stars for plot
vy_sig=ifelse(pheno_contr_rec$p.value<0.001,"***",ifelse(pheno_contr_rec$p.value<0.01,"**",ifelse(pheno_contr_rec$p.value<0.05,"*","")))
y=cbind(x,vy_sig)
odds_ratios=y
write.csv(y,"Experiment1_odds.csv")

pdf("Experiment1_odds.pdf")
odds_ratios=read.csv("Experiment1_odds.csv", header=TRUE)

odds_figure=ggplot(aes(y=vy_or,x=day,group=1),data=odds_ratios)+scale_colour_manual(values=c("black"))+
  geom_point(alpha=0.6,size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.5,2.6)+
  scale_x_discrete(name="Day",labels=c("1-2","3-4","5-6","7-9","10-12","13-15"))+geom_line()+geom_errorbar(aes(ymin=vy_or-vy_error,ymax=vy_or+vy_error))+
  annotate(geom="text", x=1, y=1.75, label=vy_sig[1],color="black",size=10)+annotate(geom="text", x=2, y=1.75, label=vy_sig[2],color="black",size=10)+annotate(geom="text", x=3, y=1.75, label=vy_sig[3],color="black",size=10)+annotate(geom="text", x=4, y=1.75, label=vy_sig[4],color="black",size=10)+annotate(geom="text", x=5, y=1.75, label=vy_sig[5],color="black",size=10)+annotate(geom="text", x=6, y=1.75, label=vy_sig[6],color="black",size=10)
odds_figure
dev.off()

#Summarization of the data

length(unique(dataset$F1.Vial[dataset$Treatment=="18"]))
length(unique(dataset$F1.Vial[dataset$Treatment=="24"]))

median(dataset$Num_moms[dataset$Treatment=="18"])
median(dataset$Num_moms[dataset$Treatment=="24"])

sum(dataset$num_CO[dataset$Treatment=="18"])+sum(dataset$num_NCO[dataset$Treatment=="18"])
sum(dataset$num_NCO[dataset$Treatment=="24"])+sum(dataset$num_CO[dataset$Treatment=="24"])

# Experiment2 -------------------------------------------------------------

#load the datasets
exp2=read.csv("Exp1_rawdata.csv", header= TRUE, na.strings="-")
exp2_backcross=read.csv("Exp1_bc_setup.csv", header=T)
exp2=na.omit(exp2)

exp2_backcross$Treatment=as.numeric(as.character(exp2_backcross$Treatment))

#recode phenotypes as characters; mutant screen used 1 s and 0 s for scoring data
#conversion to character tells R that our treatment is not a number and to treat it as a character.

#If "na.strings" is used correctly for reading in the data file, these lines are code become unnecessary
#yv$wildtype=as.numeric(as.character(yv$wildtype))
#yv$yellow.vermillion=as.numeric(as.character(yv$yellow.vermillion))
#yv$yellow.only=as.numeric(as.character(yv$yellow.only))
#yv$vermillion.only=as.numeric(as.character(yv$vermillion.only))
#yv=na.omit(yv)

#CO groups were defined

exp2$SCO1=exp2$yellow.only
exp2$SCO2=exp2$vermillion.only
exp2$NCO1=exp2$wildtype
exp2$NCO2=exp2$yellow.vermillion

sco_count=sum(exp2$SCO1+exp2$SCO2, na.rm = TRUE)
nco_count=sum(exp2$NCO1+exp2$NCO2, na.rm = TRUE)
num_samples=sum(nco_count+sco_count, na.rm = TRUE)

#merge with treatment data
exp2_merged <- merge(exp2, exp2_backcross, by.x = "Vial.number", by.y = "Vial.number", all=T)
exp2_merged = na.omit(exp2_merged)

dataset=summaryBy(SCO1+SCO2+NCO1+NCO2~F1.Vial+Day..letter.of.vial.+Treatment,data=exp2_merged, FUN=sum,na.rm=T)

#add in a column for total offspring
dataset$total_offspring=dataset$SCO1.sum+dataset$NCO1.sum+dataset$SCO2.sum+dataset$NCO2.sum

num_moms=vector(mode="numeric",length=length(dataset$F1.Vial))

#loop through dataset for each replicate to get the count
for (h in 1:length(dataset$F1.Vial)) { 
  f1_vial=dataset$F1.Vial[h]
  F1.vial=as.numeric(gsub("[^0-9\\.]","",f1_vial))
  #print(f1_vial)
  mom_ct=length(unique(sort(subset(exp2_merged,exp2_merged$F1.Vial==f1_vial)$Vial.number)))
  #store result in the vector
  num_moms[h]=mom_ct
}
dataset$Num_moms=num_moms

#use data to get fecundity calculation per replicate
dataset$fecundity=dataset$total_offspring/dataset$Num_moms

#We will write our data into a file so we can read it in for later analysis. 
write.csv(dataset,file="Experiment2_cleanedup.csv")

dataset2=read.csv("Experiment2_cleanedup.csv")
dataset2$Treatment=as.character(dataset2$Treatment)

#get mean fecundity
tapply(dataset2$fecundity,dataset2$Treatment,mean)
tapply(dataset2$fecundity,dataset2$Day..letter.of.vial.,mean)

#poisson regression, similar to a t-test for count data
fit=glm(fecundity~Treatment*Day..letter.of.vial.,data=dataset2,family=quasipoisson)
summary(fit)

#Anova
Fecundity_model=anova(fit,test="Chisq")
write.csv(Fecundity_model,"Experiment2_fecudity_model.csv")
#posthoc
fit_contrast <- emmeans::emmeans(fit, "Treatment", by="Day..letter.of.vial.", mode="kenward-roger")
fit_contr <- contrast(fit_contrast, method="trt.vs.ctrl")

pheno_contr <- as.data.frame(summary(fit_contr))
pheno_contr
write.csv(pheno_contr,"Experiment2_fecundity_posthoc.csv")

#convert p-values to stars for plot
sig=ifelse(pheno_contr$p.value<0.001,"***",ifelse(pheno_contr$p.value<0.01,"**",ifelse(pheno_contr$p.value<0.05,"*","")))

#Fecundity figure
pdf("Experiment2_fecundity.pdf")
Fecund_figure=ggplot(aes(y=fecundity,x=Day..letter.of.vial., col=Treatment,label=Num_moms),data=dataset2)+ylab("# Progeny per mom")+theme_base()

Fecund_figure=Fecund_figure+scale_colour_manual(values=c("blue", "red"))+geom_point(size=1.5)+scale_x_discrete(name="Day",labels=c("1-5","6-10","11-15","16-20"))+theme(legend.title = element_text(size = 10),legend.text = element_text(size = 8))
Fecund_figure=Fecund_figure+stat_summary(fun = median, geom="line",aes(group=Treatment),size=1.5)+geom_text(check_overlap = F,hjust = 0, nudge_x = 0.1,angle=45,size=2.5)+
  annotate(geom="text", x=1, y=40, label=sig[1],size=5)+annotate(geom="text", x=2, y=40, label=sig[2],size=5)+annotate(geom="text", x=3, y=40, label=sig[3],size=5)+annotate(geom="text", x=4, y=40, label=sig[4],size=5)
Fecund_figure

dev.off()

#Recombination rate analysis

#The progeny in vials lower than 10 is removed from the recombination analysis
dataset=dataset[dataset$total_offspring>=10,]

#sum of crossovers in intervals 1 
dataset2$num_CO=dataset2$SCO1.sum+dataset2$SCO2.sum
#sum of non-crossovers in intervals 1
dataset2$num_NCO_1=dataset2$NCO1.sum+dataset2$NCO2.sum
#total recombination rate
dataset2$rec_rate_total=dataset2$SCO1.sum/dataset2$total_offspring

#get mean recombination
tapply(dataset2$rec_rate_total,dataset2$Treatment,mean)
tapply(dataset2$rec_rate_total,dataset2$Day..letter.of.vial.,mean)

#Recombination rate model
#logistic regression, similar to a t-test for count data
dataset2$F1.Vial=as.numeric(gsub("[^0-9\\.]","",dataset2$F1.Vial))
dataset2$F1.Vial=as.numeric(dataset2$F1.Vial)

fit2=glmer(cbind(num_CO,num_NCO_1)~(1|F1.Vial)+Treatment*Day..letter.of.vial.,data=dataset2,
           family=binomial(link="logit"), control=glmerControl(optimizer="bobyqa"))
summary(fit2)
#coefs=coef(fit2)
#coefs

#Anova
Recrate_model=Anova(fit2,test="Chisq")
write.csv(Recrate_model,"Experiment2_recrate_model.csv")

#Posthoc
fit_contrast_rec <- emmeans::emmeans(fit2, "Treatment", by="Day..letter.of.vial.", mode="kenward-roger")
fit_contr_rec <- contrast(fit_contrast_rec, method="trt.vs.ctrl")

pheno_contr_rec <- as.data.frame(summary(fit_contr_rec))
pheno_contr_rec
write.csv(pheno_contr_rec,"Experiment2_recrate_posthoc.csv")

#repeat, but only for day a
#need odds ratio and standard error
#can extract from model for each time point
fit3a=glm(cbind(num_CO,num_NCO_1)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day..letter.of.vial.=="A"))
exp(coef(fit3a)[3]) #extract odds ratio for sd-y at time point 1

#Sanity Check: this is the same as extracting the exp of the estimate of the posthoc table!
exp(pheno_contr_rec$estimate[1])

#the next few lines of code don't seem necessary. 
#repeat, but only for day b
fit2b=glm(cbind(num_CO,num_NCO_1)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day..letter.of.vial.=="B"))
summary(fit2b)
exp(coef(fit2b)[3]) 

#repeat, but only for day c
fit2c=glm(cbind(num_CO,num_NCO_1)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day..letter.of.vial.=="C"))
summary(fit2c)
exp(coef(fit2c)[3]) 

#repeat, but only for day d
fit2d=glm(cbind(num_CO,num_NCO_1)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day..letter.of.vial.=="D"))
summary(fit2d)
exp(coef(fit2d)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day d

#SE and CI are related: the confidence limits are usually = estimate +/- 1.96*se (at least approximately 1.96; it actually depends on sample size depending on function)

#Now we can extract the odd ratios and SE from the posthoc table, which is a LOT cleaner and does not require additional model fits for each time point!
exp2_or=exp(pheno_contr_rec$estimate)
exp2_error=(pheno_contr_rec$SE)

#add in standard error
x=cbind(exp2_or,exp2_error)

#convert p-values to stars for plot
exp2_sig=ifelse(pheno_contr_rec$p.value<0.001,"***",ifelse(pheno_contr_rec$p.value<0.01,"**",ifelse(pheno_contr_rec$p.value<0.05,"*","")))
y=cbind(x,exp2_sig)
Day=c("a","b","c","d")
z=cbind(y,Day)
odds_ratios=z
write.csv(odds_ratios,"Experiment2_odds.csv")

odds_ratios=read.csv("Experiment2_odds.csv", header=TRUE)

#Odds ratio graph
pdf("Experiment2_odds.pdf")
odds_figure_exp2=ggplot(aes(y=exp2_or,x=Day,group=1),data=odds_ratios)+scale_colour_manual(values=c("black"))+
  geom_line(size=1.5)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.5,2.0)+
  scale_x_discrete(name="Day",labels=c("1-5","6-10","11-15","16-20"))+geom_line()+geom_errorbar(aes(ymin=exp2_or-exp2_error,ymax=exp2_or+exp2_error))+ggtitle("Recombination Rate Model")+theme(plot.title = element_text(size = 15))+
  annotate(geom="text", x=1, y=1.4, label=exp2_sig[1],color="black",size=5)+annotate(geom="text", x=2, y=1.4, label=exp2_sig[2],color="black",size=5)+annotate(geom="text", x=3, y=1.4, label=exp2_sig[3],color="black",size=5)+annotate(geom="text", x=4, y=1.4, label=exp2_sig[4],color="black",size=5)
odds_figure_exp2
dev.off()
#Summarizing the data

length(unique(dataset2$F1.Vial[dataset2$Treatment=="20"]))
length(unique(dataset2$F1.Vial[dataset2$Treatment=="25"]))

median(dataset2$Num_moms[dataset2$Treatment=="20"])
median(dataset2$Num_moms[dataset2$Treatment=="25"])

sum(dataset2$total_offspring[dataset2$Treatment=="20"])
sum(dataset2$total_offspring[dataset2$Treatment=="25"])


# Experiment3 -------------------------------------------------------------

##load the datasets

exp3=read.csv("Exp3_rawdata.csv", header= TRUE, stringsAsFactors = TRUE, na.strings="-")
exp3_backcross=read.csv("Exp3bc_setup.csv", header=T,na.strings="-")


##to make sure R reads it as a character rather than a number.

exp3_backcross$Treatment=as.character(exp3_backcross$Treatment)
exp3$wildtype=as.numeric(as.character(exp3$wildtype))
exp3$yellow.vermillion=as.numeric(as.character(exp3$yellow.vermillion))
exp3$yellow.only=as.numeric(as.character(exp3$yellow.only))
exp3$vermillion.only=as.numeric(as.character(exp3$vermillion.only))
exp3$wildtype.1=as.numeric(as.character(exp3$wildtype.1))
exp3$yellow.vermillion.1=as.numeric(as.character(exp3$yellow.vermillion.1))
exp3$yellow.only.1=as.numeric(as.character(exp3$yellow.only.1))
exp3$vermillion.only.1=as.numeric(as.character(exp3$vermillion.only.1))

##CO groups were defined

exp3$SCO1=exp3$yellow.only+exp3$yellow.only.1
exp3$SCO2=exp3$vermillion.only+exp3$vermillion.only.1
exp3$NCO1=exp3$wildtype+exp3$wildtype.1
exp3$NCO2=exp3$yellow.vermillion+exp3$vermillion.only.1

##Using the sanity check, we can see whether R interprets the results accurately.

sco_count=sum(exp3$SCO1, na.rm = TRUE)+sum(exp3$SCO2, na.rm = TRUE)
nco_count=sum(exp3$NCO1+exp3$NCO2, na.rm = TRUE)
num_samples=sum(nco_count+sco_count, na.rm = TRUE)

(sum(exp3$SCO1, na.rm = TRUE)+sum(exp3$SCO2, na.rm = TRUE))/num_samples

##merge with treatment data (This is important for summarization of the data)
exp3_merged <- merge(exp3, exp3_backcross, by.x = "Vial..", by.y = "Vial.number", all=T)
exp3_merged = na.omit(exp3_merged)
dataset=summaryBy(SCO1+SCO2+NCO1+NCO2~F1.Vial+Day..letter.of.vial.+Treatment,data=exp3_merged, FUN=sum,na.rm=T)

##calculating the fecundity
##add in a column for total offspring
dataset$total_offspring=dataset$SCO1.sum + dataset$SCO2.sum + dataset$NCO1.sum + dataset$NCO2.sum

##add a column for mothers from the F1 generation
num_moms=vector(mode="numeric",length=length(dataset$F1.Vial))


##loop through dataset to get the count
for (h in 1:length(dataset$F1.Vial)) { 
  f1_vial=dataset$F1.Vial[h]
  #  f1_vial=as.numeric(gsub("[^0-9\\.]","",f1_vial))
  #  print(f1_vial)
  mom_ct=length(unique(sort(subset(exp3_merged,exp3_merged$F1.Vial==f1_vial)$Vial..)))
  #store result in the vector
  num_moms[h]=mom_ct
}
dataset$Num_moms=num_moms

##fecundity equals to total # of offsprings divided by # of mothers
dataset$fecundity=dataset$total_offspring/dataset$Num_moms

##producing a cleanup version of the datasets
write.csv(dataset,file="Experiment3_data.csv")

##Fecundity_figure
pdf("Experiment3_fecundity.pdf")
Fecund_figure=ggplot(aes(y=fecundity,x=Day..letter.of.vial., col=Treatment,label=Num_moms),data=dataset)+ylab("# Progeny per mom")+ggtitle("Total Fecundity vs. Days post-mating")+theme_base()

Fecund_figure=Fecund_figure+scale_colour_manual(values=c("blue", "red"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("1-2","3-4","5-6","7-9","10-12","13-15"))
Fecund_figure=Fecund_figure+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
Fecund_figure
dev.off()
##remove vials with few progeny
dataset=dataset[dataset$total_offspring>=10,]


##Calculate the recombination rate for each replicate
dataset$rec_rate_total=(dataset$SCO1.sum+dataset$SCO2.sum)/(dataset$NCO1.sum+dataset$NCO2.sum+dataset$SCO1.sum+dataset$SCO2.sum)

##Recomb_figure_total
pdf("Experiment3_recombination.pdf")
Recomb_figure=ggplot(aes(y=rec_rate_total,x=Day..letter.of.vial., col=Treatment,label=total_offspring),data=dataset)+ylab("% recombination")+ggtitle("Total Recombination rate vs. Days post-mating")+theme_base()
Recomb_figure=Recomb_figure+scale_colour_manual(values=c("blue", "red"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("1-2","3-4","5-6","7-9","10-12","13-15"))
Recomb_figure=Recomb_figure+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)
Recomb_figure=Recomb_figure+ggtitle("Total Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3) +ylim(0,1.2)
Recomb_figure
dev.off()

##Statistical analysis
###the dataset necessary for the statistics is the cleaned up version of the summary data
dataset2=read.csv("Experiment3_data.csv",header=T,stringsAsFactors = F)
dataset2$Treatment=as.character(dataset2$Treatment)

###get mean fecundity
tapply(dataset2$fecundity,dataset2$Treatment,mean)
tapply(dataset2$fecundity,dataset2$Day,mean)

###poisson regression, similar to a t-test for count data
fit=glm(fecundity~Treatment*Day..letter.of.vial.,data=dataset2,family=quasipoisson)
summary(fit)

###anova results for fecundity
Experiment3_fecundity_stats=anova(fit, test="Chisq")
write.csv(Experiment3_fecundity_stats,"Experiment3_fecundity_anova.csv")

###post-hoc test for fecundity
fit_contrast <- emmeans::emmeans(fit, "Treatment", by="Day..letter.of.vial.", mode="kenward-roger")
fit_contr <- contrast(fit_contrast, method="trt.vs.ctrl")
pheno_contr <- as.data.frame(summary(fit_contr))
pheno_contr
write.csv(pheno_contr,"Experiment3_fecundity_posthoc.csv")


###Now do the same, but with recombination rate

###remove vials with few progeny
dataset2=dataset2[dataset2$total_offspring>=10,]

###sum of crossovers in exp3
dataset2$num_CO_1=dataset2$SCO1.sum+dataset2$SCO2.sum


###sum of non-crossovers in exp3
dataset2$num_NCO_1=dataset2$total_offspring-(dataset2$SCO1.sum+dataset2$SCO2.sum)

###total recombination rate
dataset2$rec_rate_exp3=(dataset2$SCO1.sum+dataset2$SCO2.sum)/dataset2$total_offspring

###get mean recombination
tapply(dataset2$rec_rate_exp3,dataset2$Treatment,mean)
tapply(dataset2$rec_rate_exp3,dataset2$Day,mean)

###to clean the dataset we get rid of the characters on the left and right of the F1 vial numbers
dataset2$F1.Vial=as.numeric(gsub("[^0-9\\.]","",dataset2$F1.Vial))
#print(dataset2$F1.Vial)

###logistic regression, similar to a t-test for count data
fit2=glmer(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+Treatment*Day..letter.of.vial.,data=dataset2,
           family=binomial(link="logit"))

###anova for the recombination rates in vy region
anova_vy=Anova(fit2,test="Chisq")
write.csv(anova_vy,"experiment3_recrate_anova.csv")

fit_contrast_rec <- emmeans::emmeans(fit2, "Treatment", by="Day..letter.of.vial.", mode="kenward-roger")
fit_contr_rec <- contrast(fit_contrast_rec, method="trt.vs.ctrl")

pheno_contr_rec <- as.data.frame(summary(fit_contr_rec))
pheno_contr_rec
write.csv(pheno_contr_rec,"Experiment3_recrate_vy_posthoc_table.csv")

#add in odds ratio and standard error
vy_or=exp(pheno_contr_rec$estimate)
vy_error=(pheno_contr_rec$SE)
x=cbind(vy_or,vy_error)
day=c("a","b","c","d","e","f")
x=cbind(day,x)
#convert p-values to stars for plot
vy_sig=ifelse(pheno_contr_rec$p.value<0.001,"***",ifelse(pheno_contr_rec$p.value<0.01,"**",ifelse(pheno_contr_rec$p.value<0.05,"*","")))
y=cbind(x,vy_sig)
odds_ratios=y
write.csv(y,"Experiment3_odds.csv")

odds_ratios=read.csv("Experiment3_odds.csv", header=TRUE)

pdf("Experiment3.odds.pdf")
odds_figure_exp3=ggplot(aes(y=vy_or,x=day,group=1),data=odds_ratios)+scale_colour_manual(values=c("black"))+
  geom_point(alpha=0.6,size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.5,2.6)+
  scale_x_discrete(name="Day",labels=c("1-2","3-4","5-6","7-9","10-12","13-15"))+geom_line()+geom_errorbar(aes(ymin=vy_or-vy_error,ymax=vy_or+vy_error))+
  annotate(geom="text", x=1, y=1.75, label=vy_sig[1],color="black",size=10)+annotate(geom="text", x=2, y=1.75, label=vy_sig[2],color="black",size=10)+annotate(geom="text", x=3, y=1.75, label=vy_sig[3],color="black",size=10)+annotate(geom="text", x=4, y=1.75, label=vy_sig[4],color="black",size=10)+annotate(geom="text", x=5, y=1.75, label=vy_sig[5],color="black",size=10)+annotate(geom="text", x=6, y=2.55, label=vy_sig[6],color="black",size=10)
odds_figure_exp3
dev.off()

sum(dataset2$total_offspring[dataset2$Treatment=="20°"])
sum(dataset2$total_offspring[dataset2$Treatment=="25°"])

# Experiment4 -------------------------------------------------------------

#read in cross data with treatment information
exp4=read.csv("Exp4_rawdata.csv",header=T)
exp4_bc_worksheet=read.csv("Exp4_backcrosses.csv",header=T,stringsAsFactors = F) 
exp4_female_counts=read.csv(file="Exp4_female.csv",header=T,stringsAsFactors = F)
exp4_bc_worksheet$Treatment=as.character(exp4_bc_worksheet$Treatment)

#recode phenotypes as characters we told R they were 1 s and 0 s so they were defined
#as characters R can understand.
#that our treatment is not a number and to treat it as a character.
exp4$sd=as.character(exp4$sd)
exp4$y=as.character(exp4$y)
exp4$se=as.character(exp4$se)

#Define Crossovers
exp4$co_class=ifelse(exp4$sd==exp4$y & exp4$y==exp4$se,"non_CO", 
                   ifelse(exp4$sd!=exp4$y & exp4$y==exp4$se,"single_CO_1",
                          ifelse(exp4$sd==exp4$y & exp4$y!=exp4$se,"single_CO_2",
                                 ifelse(exp4$sd!=exp4$y & exp4$y!=exp4$se,"double_CO",
                                        "error"))))

#Sanity check, there should not be any error!
exp4[exp4$co_class=="error",]

#add columns for counting
exp4$NCO=ifelse(exp4$co_class=="non_CO",exp4$numbMales,0)
exp4$SCO_1=ifelse(exp4$co_class=="single_CO_1",exp4$numbMales,0)
exp4$SCO_2=ifelse(exp4$co_class=="single_CO_2",exp4$numbMales,0)
exp4$DCO=ifelse(exp4$co_class=="double_CO",exp4$numbMales,0)

#get rough crossover rate as a sanity check the numbers of crossovers at the intervals 
#should equal to the total crossover rate.
nco_count=sum(exp4$numbMales*exp4$NCO)
sco_count=sum(exp4$numbMales*exp4$SCO_1, na.rm = TRUE)+sum(exp4$numbMales*exp4$SCO_2, na.rm = TRUE)
dco_count=sum(exp4$numbMales*exp4$DCO, na.rm = TRUE)
num_samples=sum(nco_count+sco_count+dco_count, na.rm = TRUE) 

#sanity check the rate of recombination for the intervals should be equal to Total.
(sco_count+(2*dco_count))/num_samples

#rate between scalloped and yellow
(sum(exp4$numbMales*exp4$SCO_1,na.rm = TRUE)+sum(exp4$numbMales*exp4$DCO,na.rm = TRUE))/num_samples

#rate between sepia and yellow
(sum(exp4$numbMales*exp4$SCO_2, na.rm = TRUE)+sum(exp4$numbMales*exp4$DCO,na.rm = TRUE))/num_samples

#we defined the crossover groups as if non crossovers will be identified as 0s; this is same as co_class column, but numerical
#single crossovers 1s
#double crossovers 2s
exp4$num_co=ifelse(exp4$y==exp4$sd & exp4$y==exp4$se,0, 
                 ifelse(exp4$sd==exp4$y & exp4$y!=exp4$se,1*exp4$numbMales, 
                        ifelse(exp4$sd!=exp4$y & exp4$y==exp4$se,1*exp4$numbMales,  
                               ifelse(exp4$sd!=exp4$y & exp4$y!=exp4$se,2*exp4$numbMales, 
                                      NA))))

#This is for summarizing our data
exp4$male=c(exp4$numbMales)

#merge with treatment data
exp4_merged <- merge(exp4, exp4_bc_worksheet, by.x = "ViaNumber", by.y = "Vial.Number", all=T)

#merge female data
exp4_female_merged=merge(exp4_female_counts, exp4_bc_worksheet, by.x = "Vial", by.y = "Vial.Number", all=T)

#summarize long form data
dataset=summaryBy(male+num_co+SCO_1+SCO_2+DCO~F1.Vial+Day+Treatment,data=exp4_merged, FUN=sum,na.rm=T)

#add in female data
exp4_female_short=summaryBy(Numbfemales~F1.Vial+Day+Treatment,data=exp4_female_merged, FUN=sum,na.rm=T)
exp4_female_short=na.omit(exp4_female_short)

#one more merge
dataset2=merge(exp4_female_short,dataset, by=c("F1.Vial","Day","Treatment"))

#add in a column for total offspring
dataset2$total_offspring=dataset2$Numbfemales.sum + dataset2$male.sum

#make a vector to store the data
num_moms=vector(mode="numeric",length=length(dataset2$F1.Vial))

#loop through dataset to get the count
for (h in 1:length(dataset2$F1.Vial)) { 
  f1_vial=dataset2$F1.Vial[h]
  mom_ct=length(unique(sort(subset(exp4_merged,exp4_merged$F1.Vial==f1_vial)$ViaNumber)))
  #store result in the vector
  num_moms[h]=mom_ct
}

#add vector as a column in dataset2
dataset2$Num_moms=num_moms

#use data to get fecundity calculation
dataset2$fecundity=dataset2$total_offspring/dataset2$Num_moms
dataset2$Treatment=as.factor(dataset2$Treatment)
#We will write our data into a file so we can read it in for later analysis. 
write.csv(dataset2,file="Experiment4_cleanedup.csv")

#Fecundity_figure
#dataset2=read.csv("Experiment4_cleanedup.csv")

#poisson regression, similar to a t-test for count data
fit=glm(fecundity~Treatment*Day,data=dataset2,family=quasipoisson)
summary(fit)
anova_fec=anova(fit, test="Chisq")
anova_fec
write.csv(anova_fec,"Exp4_fecundity_model_table.csv")


fit_contrast <- emmeans::emmeans(fit, "Treatment", by="Day", mode="kenward-roger")
fit_contr <- contrast(fit_contrast, method="trt.vs.ctrl")

pheno_contr <- as.data.frame(summary(fit_contr))
pheno_contr
write.csv(pheno_contr,"Exp4_fecundity_posthoc_table.csv")

#convert p-values to stars for plot
sig=ifelse(pheno_contr$p.value<0.001,"***",ifelse(pheno_contr$p.value<0.01,"**",ifelse(pheno_contr$p.value<0.05,"*","")))

#Fecundity figure for the paper
pdf("Exp4_Fecundity.pdf")
Fecund_figure=ggplot(aes(y=fecundity,x=Day, col=Treatment,label=Num_moms),data=dataset2)+ylab("# Progeny per mom")+ggtitle("Total Fecundity vs. Days post-mating")+theme_base()

Fecund_figure=Fecund_figure+scale_colour_manual(values=c("blue", "red"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))
Fecund_figure=Fecund_figure+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)+
  annotate(geom="text", x=1, y=120, label=sig[1],size=10)+annotate(geom="text", x=2, y=120, label=sig[2],size=10)+annotate(geom="text", x=3, y=120, label=sig[3],size=10)+annotate(geom="text", x=4, y=120, label=sig[4],size=10)
Fecund_figure

dev.off()

#get averages for text; These are the numbers we put in the paper!
tapply(dataset2$fecundity,dataset2$Treatment,mean,na.rm=T)
tapply(dataset2$fecundity,dataset2$Day,mean,na.rm=T)

#remove vials with few progeny
dataset2=dataset2[dataset2$male.sum>=10,]
dataset2=dataset2[dataset2$Numbfemales.sum>=10,]

#sum of crossovers in intervals 1 & 2 
dataset2$num_CO_1=dataset2$SCO_1.sum+dataset2$DCO.sum
dataset2$num_CO_2=dataset2$SCO_2.sum+dataset2$DCO.sum

#sum of non-crossovers in intervals 1 & 2
dataset2$num_NCO_1=dataset2$male.sum-(dataset2$SCO_1.sum+dataset2$DCO.sum)
dataset2$num_NCO_2=dataset2$male.sum-(dataset2$SCO_2.sum+dataset2$DCO.sum)

#total recombination rate
dataset2$rec_rate_total=dataset2$num_co.sum/dataset2$male.sum
dataset2$rec_rate_ysd=(dataset2$SCO_1.sum+dataset2$DCO.sum)/dataset2$male.sum
dataset2$rec_rate_yse=(dataset2$SCO_2.sum+dataset2$DCO.sum)/dataset2$male.sum

pdf("Experiment4_vialsremoved.pdf")
#Recomb_figure_total
Recomb_figure=ggplot(aes(y=rec_rate_total,x=Day, col=as.factor(Treatment),label=male.sum),data=dataset2)+ylab("% recombination")+ggtitle("Total Recombination rate vs. Days post-mating")+theme_base()

#Now add points and change labels for x-axis
Recomb_figure=Recomb_figure+scale_colour_manual(values=c("blue", "red"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))
Recomb_figure
#Add a line through the median of the points
Recomb_figure=Recomb_figure+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)

#Add a label for sample size
Recomb_figure=Recomb_figure+ggtitle("Total Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3) +ylim(0.2,1)

#print figure
Recomb_figure

#Recomb_figure_yellow_scalloped
Recomb_figure=ggplot(aes(y=rec_rate_ysd,x=Day, col=as.factor(Treatment),label=male.sum),data=dataset2)+scale_colour_manual(values=c("red", "blue"))+geom_point(alpha=0.6,size=3)+ylab("% recombination")+theme_base()+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))
Recomb_figure=Recomb_figure+ggtitle("y-sd Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
Recomb_figure

#Recomb_figure_yellow_sepia
Recomb_figure=ggplot(aes(y=rec_rate_yse,x=Day, col=as.factor(Treatment),label=male.sum),data=dataset2)+scale_colour_manual(values=c("red", "blue"))+geom_point(alpha=0.6,size=3)+ylab("% recombination")+theme_base()+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))
Recomb_figure=Recomb_figure+ggtitle("y-se Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
Recomb_figure

dev.off()

#get averages for text; These are the numbers we put in the paper!
tapply(dataset2$rec_rate_total,dataset2$Treatment,mean,na.rm=T)
tapply(dataset2$rec_rate_total,dataset2$Day,mean,na.rm=T)

tapply(dataset2$rec_rate_ysd,dataset2$Treatment,mean,na.rm=T)
tapply(dataset2$rec_rate_ysd,dataset2$Day,mean,na.rm=T)

tapply(dataset2$rec_rate_yse,dataset2$Treatment,mean,na.rm=T)
tapply(dataset2$rec_rate_yse,dataset2$Day,mean,na.rm=T)

#Odds ratios
#SD-Y REGION

#logistic regression, similar to a t-test for count data
fit3=glmer(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+Treatment*Day,data=dataset2,
           family=binomial(link="logit"))
#coefs=coef(fit3)
#coefs
anova_sdy=Anova(fit3,test="Chisq")
write.csv(anova_sdy,"Exp4_recrate_sd-y_model_table.csv")

fit_contrast3 <- emmeans::emmeans(fit3, "Treatment", by="Day", mode="kenward-roger")
fit_contr3 <- contrast(fit_contrast3, method="trt.vs.ctrl")

pheno_contr3 <- as.data.frame(summary(fit_contr3))
pheno_contr3
write.csv(pheno_contr3,"Exp4_recrate_sd-y_posthoc_table.csv")

#need odds ratio and standard error
#can extract from model for each time point
#repeat, but only for day A

fit3a=glm(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day=="A"))
exp(coef(fit3a)[3]) #extract odds ratio for sd-y at time point 1

#Sanity Check: this is the same as extracting the exp of the estimate of the posthoc table!
exp(pheno_contr3$estimate[1])
#At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day v

#SE and CI are related: the confidence limits are usually = estimate +/- 1.96*se (at least approximately 1.96; it actually depends on sample size depending on function)
#We can use this as a sanity check to make sure the SE in the posthoc table is similar to the CI we could extract from the model

#extract lower and upper 95% CI from model above
exp(confint(fit3a)[3])
exp(confint(fit3a)[6])

#compare to SE in posthoc table
exp((pheno_contr3$estimate[1])+(pheno_contr3$SE[1]*1.96))
exp((pheno_contr3$estimate[1])-(pheno_contr3$SE[1]*1.96))

#great, they are VERY similar. Now we can extract the odd ratios and SE from the posthoc table, which is a LOT cleaner and does not require additional model fits for each time point!
y_sd_or=exp(pheno_contr3$estimate)
y_sd_error=(pheno_contr3$SE)


#Y-SE REGION

#logistic regression, similar to a t-test for count data
fit4=glmer(cbind(num_CO_2,num_NCO_2)~(1|F1.Vial)+Treatment*Day,data=dataset2,
         family=binomial(link="logit"))
#summary(fit4)
#coefs=coef(fit4)
#coefs
anova_yse=Anova(fit4,test="Chisq")
write.csv(anova_yse,"Exp4_recrate_y-se_model_table.csv")

fit_contrast4 <- emmeans::emmeans(fit4, "Treatment", by="Day", mode="kenward-roger")
fit_contr4 <- contrast(fit_contrast4, method="trt.vs.ctrl")

pheno_contr4 <- as.data.frame(summary(fit_contr4))
pheno_contr4
write.csv(pheno_contr4,"Exp4_recrate_y-se_posthoc_table.csv")

y_se_or=exp(pheno_contr4$estimate)
y_se_error=(pheno_contr4$SE)

#setup data frame
odds_ratios=cbind(y_sd_or,y_se_or)
colnames(odds_ratios)=c("sd-y","y-se")
rownames(odds_ratios)=c("1-3","4-6","7-9","10-12")

#melt dataframe into long form
odds_ratios=melt(odds_ratios)
colnames(odds_ratios)=c("Day","Interval","OR")
odds_ratios

#add in standard error
SE=c(y_sd_error,y_se_error)
x=cbind(odds_ratios,SE)
#x

#convert p-values to stars for plot
ysd_sig=ifelse(pheno_contr3$p.value<0.001,"***",ifelse(pheno_contr3$p.value<0.01,"**",ifelse(pheno_contr3$p.value<0.05,"*","")))
yse_sig=ifelse(pheno_contr4$p.value<0.001,"***",ifelse(pheno_contr4$p.value<0.01,"**",ifelse(pheno_contr4$p.value<0.05,"*","")))

#add significance to table
sig=c(ysd_sig,yse_sig)
y=cbind(x,sig)
odds_ratios=y
odds_ratios
write.csv(odds_ratios,"Experiment4_odds.csv")

#Odds ratio plot
pdf("Exp4_odds_ratio.pdf")
odds_figure_exp4=ggplot(aes(y=OR,x=Day, col=Interval,group=Interval),data=odds_ratios)+scale_colour_manual(values=c("#f1a340", "#998ec3"))+
  geom_point(alpha=0.6,size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.5,2.0)+
  scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))+geom_line()+geom_errorbar(aes(ymin=OR-SE,ymax=OR+SE))+
  annotate(geom="text", x=1, y=1.75, label=ysd_sig[1],color="#f1a340",size=10)+annotate(geom="text", x=2, y=1.75, label=ysd_sig[2],color="#f1a340",size=10)+annotate(geom="text", x=3, y=1.75, label=ysd_sig[3],color="#f1a340",size=10)+annotate(geom="text", x=4, y=1.75, label=ysd_sig[4],color="#f1a340",size=10)+
  annotate(geom="text", x=1, y=1.7, label=yse_sig[1],color="#998ec3",size=10)+annotate(geom="text", x=2, y=1.7, label=yse_sig[2],color="#998ec3",size=10)+annotate(geom="text", x=3, y=1.7, label=yse_sig[3],color="#998ec3",size=10)+annotate(geom="text", x=4, y=1.7, label=yse_sig[4],color="#998ec3",size=10)
odds_figure_exp4
dev.off()


# Experiment5 -------------------------------------------------------------

#Data files are uploaded
exp5=read.csv("Exp5_rawdata.csv",header=T)
exp5_bc_worksheet=read.csv("Exp5_backcrosses.csv",header=T,stringsAsFactors = F) 
exp5_female_counts=read.csv(file="Exp5_females.csv",header=T,stringsAsFactors = F)

#recode phenotypes as characters we told R they were 1 s and 0 s so they were defined
#as characters R can understand.
#that our treatment is not a number and to treat it as a character.

exp5_bc_worksheet$Treatment=as.character(exp5_bc_worksheet$Treatment)

exp5$sd=as.character(exp5$sd)
exp5$y=as.character(exp5$y)
exp5$se=as.character(exp5$se)

#Define Crossovers
exp5$co_class=ifelse(exp5$sd==exp5$y & exp5$y==exp5$se,"non_CO", 
                   ifelse(exp5$sd!=exp5$y & exp5$y==exp5$se,"single_CO_1",
                          ifelse(exp5$sd==exp5$y & exp5$y!=exp5$se,"single_CO_2",
                                 ifelse(exp5$sd!=exp5$y & exp5$y!=exp5$se,"double_CO",
                                        "error"))))

#check for errors, which are the removed cut phenotypes.
exp5[exp5$co_class=="error",]

#add columns for counting
exp5$NCO=ifelse(exp5$co_class=="non_CO",exp5$numbMales,0)
exp5$SCO_1=ifelse(exp5$co_class=="single_CO_1",exp5$numbMales,0)
exp5$SCO_2=ifelse(exp5$co_class=="single_CO_2",exp5$numbMales,0)
exp5$DCO=ifelse(exp5$co_class=="double_CO",exp5$numbMales,0)

#get rough crossover rate as a sanity check the numbers of crossovers at the intervals 
#should equal to the total crossover rate.
nco_count=sum(exp5$numbMales*exp5$NCO)
sco_count=sum(exp5$numbMales*exp5$SCO_1, na.rm = TRUE)+sum(exp5$numbMales*exp5$SCO_2, na.rm = TRUE)
dco_count=sum(exp5$numbMales*exp5$DCO, na.rm = TRUE)
num_samples=sum(nco_count+sco_count+dco_count, na.rm = TRUE) 

#total rate
(sco_count+(2*dco_count))/num_samples

#rate between scalloped and yellow
(sum(exp5$numbMales*exp5$SCO_1,na.rm = TRUE)+sum(exp5$numbMales*exp5$DCO,na.rm = TRUE))/num_samples

#rate between sepia and yellow
(sum(exp5$numbMales*exp5$SCO_2, na.rm = TRUE)+sum(exp5$numbMales*exp5$DCO,na.rm = TRUE))/num_samples

#we defined the crossover groups as if non crossovers will be identified as 0s; this is same as co_class column, but numerical
#single crossovers 1s
#double crossovers 2s
exp5$num_co=ifelse(exp5$y==exp5$sd & exp5$y==exp5$se,0, 
                 ifelse(exp5$sd==exp5$y & exp5$y!=exp5$se,1*exp5$numbMales, 
                        ifelse(exp5$sd!=exp5$y & exp5$y==exp5$se,1*exp5$numbMales,  
                               ifelse(exp5$sd!=exp5$y & exp5$y!=exp5$se,2*exp5$numbMales, 
                                      NA))))

#This is for summarizing our data
exp5$male=c(exp5$numbMales)

#merge with treatment data
exp5_merged <- merge(exp5, exp5_bc_worksheet, by.x = "ViaNumber", by.y = "Vial.Number", all=T)

#merge female data
exp5_female_merged=merge(exp5_female_counts, exp5_bc_worksheet, by.x = "Vial", by.y = "Vial.Number", all=T)

#summarize long form data
dataset=summaryBy(male+num_co+SCO_1+SCO_2+DCO~F1.Vial+Day+Treatment,data=exp5_merged, FUN=sum,na.rm=T)

#add in female data
exp5_female_short=summaryBy(Numbfemales~F1.Vial+Day+Treatment,data=exp5_female_merged, FUN=sum,na.rm=T)
exp5_female_short=na.omit(exp5_female_short)

#one more merge
dataset2=merge(exp5_female_short,dataset, by=c("F1.Vial","Day","Treatment"))

#add in a column for total offspring
dataset2$total_offspring=dataset2$Numbfemales.sum + dataset2$male.sum

#make a vector to store the data
num_moms=vector(mode="numeric",length=length(dataset2$F1.Vial))

#loop through dataset to get the count
for (h in 1:length(dataset2$F1.Vial)) { 
  f1_vial=dataset2$F1.Vial[h]
  mom_ct=length(unique(sort(subset(exp5_merged,exp5_merged$F1.Vial==f1_vial)$ViaNumber)))
  #store result in the vector
  num_moms[h]=mom_ct
}

#add vector as a column in dataset2
dataset2$Num_moms=num_moms

#use data to get fecundity calculation
dataset2$fecundity=dataset2$total_offspring/dataset2$Num_moms

#We will write our data into a file so we can read it in for later analysis. 
write.csv(dataset2,file="Experiment5_cleanedup.csv")

#poisson regression, similar to a t-test for count data
fit=glm(fecundity~Treatment*Day,data=dataset2,family=quasipoisson)
#summary(fit)
anova_fec=anova(fit, test="Chisq")
write.csv(anova_fec,"Exp5_fecundity_model_table.csv")

fit_contrast <- emmeans::emmeans(fit, "Treatment", by="Day", mode="kenward-roger")
fit_contr <- contrast(fit_contrast, method="trt.vs.ctrl")

pheno_contr <- as.data.frame(summary(fit_contr))
pheno_contr
write.csv(pheno_contr,"Exp5_fecundity_posthoc_table.csv")

#convert p-values to stars for plot
sig=ifelse(pheno_contr$p.value<0.001,"***",ifelse(pheno_contr$p.value<0.01,"**",ifelse(pheno_contr$p.value<0.05,"*","")))

#Fecundity figure for the paper
pdf("Exp5_Fecundity.pdf")
Fecund_figure=ggplot(aes(y=fecundity,x=Day, col=Treatment,label=Num_moms),data=dataset2)+ylab("# Progeny per mom")+ggtitle("Total Fecundity vs. Days post-mating")+theme_base()

Fecund_figure=Fecund_figure+scale_colour_manual(values=c("blue", "red"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("6","7","8","9","10"))+ylim(0,45)
Fecund_figure=Fecund_figure+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)+
  annotate(geom="text", x=1, y=40, label=sig[1],size=10)+annotate(geom="text", x=2, y=40, label=sig[2],size=10)+annotate(geom="text", x=3, y=40, label=sig[3],size=10)+annotate(geom="text", x=4, y=40, label=sig[4],size=10)
Fecund_figure

dev.off()

#get averages for text; These are the numbers we put in the paper!
tapply(dataset2$fecundity,dataset2$Treatment,mean,na.rm=T)
tapply(dataset2$fecundity,dataset2$Day,mean,na.rm=T)

#remove vials with few progeny
dataset2=dataset2[dataset2$male.sum>=10,]
dataset2=dataset2[dataset2$Numbfemales.sum>=10,]

#sum of crossovers in intervals 1 & 2 
dataset2$num_CO_1=dataset2$SCO_1.sum+dataset2$DCO.sum
dataset2$num_CO_2=dataset2$SCO_2.sum+dataset2$DCO.sum

#sum of non-crossovers in intervals 1 & 2
dataset2$num_NCO_1=dataset2$male.sum-(dataset2$SCO_1.sum+dataset2$DCO.sum)
dataset2$num_NCO_2=dataset2$male.sum-(dataset2$SCO_2.sum+dataset2$DCO.sum)

#total recombination rate
dataset2$rec_rate_total=dataset2$num_co.sum/dataset2$male.sum
dataset2$rec_rate_ysd=(dataset2$SCO_1.sum+dataset2$DCO.sum)/dataset2$male.sum
dataset2$rec_rate_yse=(dataset2$SCO_2.sum+dataset2$DCO.sum)/dataset2$male.sum
mean(dataset2$rec_rate_total)
mean(dataset2$rec_rate_ysd)
mean(dataset2$rec_rate_yse)

###Recombination figures

pdf("Experiment5_recombination.pdf")
#Recomb_figure_total
Recomb_figure=ggplot(aes(y=rec_rate_total,x=Day, col=as.factor(Treatment),label=male.sum),data=dataset2)+ylab("% recombination")+ggtitle("Total Recombination rate vs. Days post-mating")+theme_base()
Recomb_figure=Recomb_figure+scale_colour_manual(values=c("blue", "red"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("6","7","8","9","10"))
Recomb_figure=Recomb_figure+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)
Recomb_figure=Recomb_figure+ggtitle("Total Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3) +ylim(0.2,1)
Recomb_figure

#Recomb_figure_yellow_scalloped
Recomb_figure=ggplot(aes(y=rec_rate_ysd,x=Day, col=as.factor(Treatment),label=male.sum),data=dataset2)+scale_colour_manual(values=c("blue", "red"))+geom_point(alpha=0.6,size=3)+ylab("% recombination")+theme_base()+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))
Recomb_figure=Recomb_figure+ggtitle("y-sd Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
Recomb_figure

#Recomb_figure_yellow_sepia
Recomb_figure=ggplot(aes(y=rec_rate_yse,x=Day, col=as.factor(Treatment),label=male.sum),data=dataset2)+scale_colour_manual(values=c("blue","red"))+geom_point(alpha=0.6,size=3)+ylab("% recombination")+theme_base()+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+scale_x_discrete(name="Day",labels=c("6","7","8","9","10"))
Recomb_figure=Recomb_figure+ggtitle("y-se Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
Recomb_figure
dev.off()

#get averages for text; These are the numbers we put in the paper!
tapply(dataset2$rec_rate_total,dataset2$Treatment,mean,na.rm=T)
tapply(dataset2$rec_rate_total,dataset2$Day,mean,na.rm=T)

d=subset(dataset2, dataset2$Day=="Y")
tapply(d$rec_rate_total,d$Treatment,mean,na.rm=T)
tapply(d$rec_rate_ysd,d$Treatment,mean,na.rm=T)
tapply(d$rec_rate_yse,d$Treatment,mean,na.rm=T)
e=subset(dataset2,dataset2$Day=="V")
tapply(e$rec_rate_yse,e$Treatment,mean,na.rm=T)
tapply(e$rec_rate_yse,e$Day,mean,na.rm=T)

tapply(dataset2$rec_rate_ysd,dataset2$Treatment,mean,na.rm=T)
tapply(dataset2$rec_rate_ysd,dataset2$Day,mean,na.rm=T)

tapply(dataset2$rec_rate_yse,dataset2$Treatment,mean,na.rm=T)
tapply(dataset2$rec_rate_yse,dataset2$Day,mean,na.rm=T)

###Odds ratios
#SD-Y REGION

#logistic regression, similar to a t-test for count data
fit3=glmer(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+Treatment*Day,data=dataset2,
         family=binomial(link="logit"))
#coefs=coef(fit3)
anova_sdy=Anova(fit3,test="Chisq")
write.csv(anova_sdy,"Exp5_recrate_sd-y_model_table.csv")

fit_contrast3 <- emmeans::emmeans(fit3, "Treatment", by="Day", mode="kenward-roger")
fit_contr3 <- contrast(fit_contrast3, method="trt.vs.ctrl")

pheno_contr3 <- as.data.frame(summary(fit_contr3))
pheno_contr3
write.csv(pheno_contr3,"Exp5_recrate_sd-y_posthoc_table.csv")

#we can extract the odd ratios and SE from the posthoc table, which is a LOT cleaner and does not require additional model fits for each time point!
y_sd_or=exp(pheno_contr3$estimate)
y_sd_error=(pheno_contr3$SE)

#Y-SE REGION

#logistic regression, similar to a t-test for count data
fit4=glmer(cbind(num_CO_2,num_NCO_2)~(1|F1.Vial)+Treatment*Day,data=dataset2,
         family=binomial(link="logit"))
#summary(fit4)
#coefs=coef(fit4)
#coefs
anova_yse=Anova(fit4,test="Chisq")
write.csv(anova_yse,"Exp5_recrate_y-se_model_table.csv")

fit_contrast4 <- emmeans::emmeans(fit4, "Treatment", by="Day", mode="kenward-roger")
fit_contr4 <- contrast(fit_contrast4, method="trt.vs.ctrl")

pheno_contr4 <- as.data.frame(summary(fit_contr4))
pheno_contr4
write.csv(pheno_contr4,"Exp5_recrate_y-se_posthoc_table.csv")

y_se_or=exp(pheno_contr4$estimate)
y_se_error=(pheno_contr4$SE)

#setup data frame
odds_ratios=cbind(y_sd_or,y_se_or)
colnames(odds_ratios)=c("sd-y","y-se")
rownames(odds_ratios)=c("6","7","8","9","10")

#melt dataframe into long form
odds_ratios=melt(odds_ratios)
colnames(odds_ratios)=c("Day","Interval","OR")
odds_ratios

#add in standard error
SE=c(y_sd_error,y_se_error)
x=cbind(odds_ratios,SE)

#convert p-values to stars for plot
ysd_sig=ifelse(pheno_contr3$p.value<0.001,"***",ifelse(pheno_contr3$p.value<0.01,"**",ifelse(pheno_contr3$p.value<0.05,"*","")))
yse_sig=ifelse(pheno_contr4$p.value<0.001,"***",ifelse(pheno_contr4$p.value<0.01,"**",ifelse(pheno_contr4$p.value<0.05,"*","")))

#add significance to table
sig=c(ysd_sig,yse_sig)
y=cbind(x,sig)
odds_ratios=y

#Odds ratio plot
pdf("Exp5_odds_ratio.pdf")
odds_figure_exp5=ggplot(aes(y=OR,x=Day, col=Interval,group=Interval),data=odds_ratios)+scale_colour_manual(values=c("#f1a340", "#998ec3"))+
  geom_point(alpha=0.6,size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.5,1.8)+
  scale_x_discrete(name="Days post-mating",labels = c("6","7","8","9","10"))+geom_line()+geom_errorbar(aes(ymin=OR-SE,ymax=OR+SE))+
  annotate(geom="text", x=6, y=1.75, label=ysd_sig[1],color="#f1a340",size=10)+annotate(geom="text", x=7, y=1.75, label=ysd_sig[2],color="#f1a340",size=10)+annotate(geom="text", x=8, y=1.75, label=ysd_sig[3],color="#f1a340",size=10)+annotate(geom="text", x=9, y=1.75, label=ysd_sig[4],color="#f1a340",size=10)+
  annotate(geom="text", x=6, y=1.7, label=yse_sig[1],color="#998ec3",size=10)+annotate(geom="text", x=7, y=1.7, label=yse_sig[2],color="#998ec3",size=10)+annotate(geom="text", x=8, y=1.7, label=yse_sig[3],color="#998ec3",size=10)+annotate(geom="text", x=9, y=1.7, label=yse_sig[4],color="#998ec3",size=10)
odds_figure_exp5
dev.off()

#boxplot for total recombination rate
pdf("Exp5_sdyboxplotonly.pdf")

Recomb_figure=ggplot(aes(y=rec_rate_ysd,x=Day, col=Treatment,label=male.sum),data=dataset2)+ylab("% recombination")+ggtitle("Total Recombination rate vs. Days post-mating")+theme_update(text = element_text(size=30))
Recomb_figure= Recomb_figure+ theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
Recomb_figure=Recomb_figure+scale_colour_manual(values=c("blue", "red"))+annotate(geom="text", x=1, y=0.4, label=ysd_sig[1],color="black",size=10)+annotate(geom="text", x=2, y=0.4, label=ysd_sig[2],color="black",size=7)+annotate(geom="text", x=3, y=0.4, label=ysd_sig[3],color="black",size=10)+annotate(geom="text", x=4, y=0.4, label=ysd_sig[4],color="black",size=10)
Recomb_figure <- Recomb_figure + geom_boxplot() # Adds color
Recomb_figure <- Recomb_figure + scale_x_discrete(name="Days", labels=c("6","7","8","9","10")) # Adds kaveks
Recomb_figure # displays the boxplots

dev.off()


pdf("Exp5_yseboxplotonly.pdf")

Recomb_figure=ggplot(aes(y=rec_rate_yse,x=Day, col=Treatment,label=male.sum),data=dataset2)+ylab("% recombination")+ggtitle("Total Recombination rate vs. Days post-mating")+theme_update(text = element_text(size=30))
Recomb_figure= Recomb_figure+ theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))+ylim(0.1,0.65)
Recomb_figure=Recomb_figure+scale_colour_manual(values=c("blue", "red"))+annotate(geom="text", x=1, y=0.4, label=yse_sig[1],color="black",size=10)+annotate(geom="text", x=2, y=0.4, label=yse_sig[2],color="black",size=7)+annotate(geom="text", x=3, y=0.4, label=yse_sig[3],color="black",size=10)+annotate(geom="text", x=4, y=0.4, label=yse_sig[4],color="black",size=10)
Recomb_figure <- Recomb_figure + geom_boxplot() # Adds color
Recomb_figure <- Recomb_figure + scale_x_discrete(name="Days", labels=c("6","7","8","9","10")) # Adds kaveks
Recomb_figure # displays the boxplots

dev.off()



# Reproducibility ---------------------------------------------------------

#Fecundity and recombination rate comparison
#Subset is for the flies who are collected post mating day 9.
x=subset(dataset2,dataset2$Day=="Y")

#So we can compare against the control and stress temperetures.
control=subset(dataset2,dataset2$Day=="Y" & dataset2$Treatment=="21")
hitemp=subset(dataset2,dataset2$Day=="Y" & dataset2$Treatment=="26")

pdf("recrate_vs_fecundity.pdf")
plot(x$fecundity,x$rec_rate_ysd,cex.axis=1.5)
points(control$fecundity,control$rec_rate_ysd, col="blue",pch=19,cex=1.5)
points(hitemp$fecundity,hitemp$rec_rate_ysd, col="red",pch=19,cex=1.5)
dev.off()

#Statistics that reports the p values as well as the R^2.
a=cor.test(control$fecundity,control$rec_rate_ysd)
b=cor.test(hitemp$fecundity,hitemp$rec_rate_ysd)

#Next step is to investigate the producibility between our experiments. 
#As Experiment 2-3 and 4-5 are comperable based on time scales.

#Odds for 7-9 merged in dataset2
# this is done to compare directly between the experiments for sdy and yse
#First we merge the days 7, 8, and 9
dataset2$new_day=ifelse(dataset2$Day=="W","C",
                        ifelse(dataset2$Day=="X","C",
                               ifelse(dataset2$Day=="Y","C","NA")))
only_C=subset(dataset2,dataset2$new_day=="C")
write.csv(only_C, file= "Experiment_5_onlyc.csv")

#model for the new day arrangement and odds ratio at the days 7-9 for Experiment5
#Odds for the timepoint c and the rest. (Remaining will be days 6 and 10 and they are annotated with "NA")
fit3c=glmer(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+Treatment*new_day,data=dataset2,
          family=binomial(link="logit"))
fit_contrast3c <- emmeans::emmeans(fit3c, "Treatment", by="new_day", mode="kenward-roger")
fit_contr3c <- contrast(fit_contrast3c, method="trt.vs.ctrl")

pheno_contr3c <- as.data.frame(summary(fit_contr3c))
pheno_contr3c
write.csv(pheno_contr3,"Exp5_recrate_sd-y_c_posthoc_table.csv")

#Let's repeat the same for y-se interval
fit4c=glmer(cbind(num_CO_2,num_NCO_2)~(1|F1.Vial)+Treatment*new_day,data=dataset2,
          family=binomial(link="logit"))

fit_contrast4c <- emmeans::emmeans(fit4c, "Treatment", by="new_day", mode="kenward-roger")
fit_contr4c <- contrast(fit_contrast4c, method="trt.vs.ctrl")

pheno_contr4c <- as.data.frame(summary(fit_contr4c))
pheno_contr4c
write.csv(pheno_contr4c,"Exp5_recrate_y-se_c_posthoc_table.csv")

#Extract the "c" odds and errors
y_sd_or_c=exp(pheno_contr3c$estimate)
y_sd_error_c=(pheno_contr3c$SE)
y_se_or_c=exp(pheno_contr4c$estimate)
y_se_error_c=(pheno_contr4c$SE)

odds_ratios=cbind(y_sd_or_c,y_se_or_c)
colnames(odds_ratios)=c("sd-y","y-se")
rownames(odds_ratios)=c("7-9","NA")

odds_ratios=melt(odds_ratios)
colnames(odds_ratios)=c("Day","Interval","OR")
odds_ratios

#add in standard error
SE=c(y_sd_error_c,y_se_error_c)
Exp5_odds_c=cbind(odds_ratios,SE)

write.csv(Exp5_odds_c,"Experiment5_odds_c.csv")

#The odds from previous experiments for the reproducibility.
#We will be comparing between Experiments 2-3 for v-y, Experiment4-5 for sd-y, Experiment4-5 for y-se.
Exp5_odds=read.csv("Experiment5_odds_c.csv",header = TRUE)
Exp4_odds=read.csv("Experiment4_odds.csv",header = TRUE)
Exp3_odds=read.csv("Experiment3_odds.csv",header = TRUE)
Exp2_odds=read.csv("Experiment2_odds.csv",header = TRUE)

#Data frame includes the points for 7-9days from each experiment as well as the errors. 
Reproducibility <- data.frame(interval=c("yst","yst","sdy","yse","sdy","yse"),
                   odds=c(Exp2_odds[4,3],Exp3_odds[4,3],Exp4_odds[3,4],Exp4_odds[7,4],Exp5_odds[1,4],Exp5_odds[3,4]),
                   error=c(Exp2_odds[4,4],Exp3_odds[4,4],Exp4_odds[3,5],Exp4_odds[7,5],Exp5_odds[1,5],Exp5_odds[3,5]),
                   experiment=c("Exp1","Exp3","Exp4","Exp4","Exp5","Exp5"))

#Graph for reproducibility

pdf("Reproducibility.pdf")
Reproducibility_figure=ggplot(Reproducibility, aes(x=interval, y=odds, group=experiment, col=experiment))+
  geom_point(size=3)+theme_bw()+
  geom_errorbar(aes(ymin=odds-error, ymax=odds+error), width=.2,
                position=position_dodge(0.05))
Reproducibility_figure
dev.off()


# SNP analysis ------------------------------------------------------------

#Cleanedup dataset from the Preliminary genotyping analysis
genotyping=read.csv("SNP_genotyping.csv",header = TRUE)
genotyping$Treatment=as.character(genotyping$Treatment)
genotyping$F1.Vial=as.numeric(genotyping$F1.Vial)

genotyping$repID=as.numeric(paste(genotyping$F1.Vial,genotyping$Treatment,sep="."))

#Summarizing the data

length(unique(genotyping$F1.Vial[genotyping$Treatment=="18"]))
length(unique(genotyping$F1.Vial[genotyping$Treatment=="23"]))

median(genotyping$Num_moms[genotyping$Treatment=="18"])
median(genotyping$Num_moms[genotyping$Treatment=="23"])

sum(genotyping$count.sum[genotyping$Treatment=="18"])
sum(genotyping$count.sum[genotyping$Treatment=="23"])

#get mean recombination
tapply(genotyping$rec_rate_total,genotyping$Treatment,mean)
tapply(genotyping$rec_rate_total,genotyping$Day,mean)

#i1 REGION
#logistic regression, similar to a t-test for count data
fit1=glmer(cbind(num_CO_1,num_NCO_1)~(1|repID)+Treatment*Day,data=genotyping,
         family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa"))
#summary(fit1)

I1_model=Anova(fit1,test="Chisq")
write.csv(I1_model,"SNP_i1_model.csv")
fit_contrast_i1 <- emmeans::emmeans(fit1, "Treatment", by="Day", mode="kenward-roger")
fit_contr_i1 <- contrast(fit_contrast_i1, method="trt.vs.ctrl")

pheno_contr_i1 <- as.data.frame(summary(fit_contr_i1))
pheno_contr_i1
write.csv(pheno_contr_i1,"SNP_interval1_posthoc.csv")

#extract the odd ratios and SE from the posthoc table, which is a LOT cleaner and does not require additional model fits for each time point!
i1_or=exp(pheno_contr_i1$estimate)
i1_error=(pheno_contr_i1$SE)

#i2 region
#logistic regression, similar to a t-test for count data
fit3=glmer(cbind(num_CO_2,num_NCO_2)~(1|repID)+Treatment*Day,data=genotyping,
         family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa"))
summary(fit3)

I2_model=Anova(fit3,test="Chisq")
write.csv(I2_model,"SNP_i2_model.csv")
fit_contrast_i2 <- emmeans::emmeans(fit3, "Treatment", by="Day", mode="kenward-roger")
fit_contr_i2 <- contrast(fit_contrast_i2, method="trt.vs.ctrl")

pheno_contr_i2 <- as.data.frame(summary(fit_contr_i2))
pheno_contr_i2
write.csv(pheno_contr_i2,"SNP_interval2_posthoc.csv")
i2_or=exp(pheno_contr_i2$estimate)
i2_error=(pheno_contr_i2$SE)

#i3 region

#geno2=subset(genotyping, genotyping$F1.Vial!="5")
#brief=subset(genotyping[,c(2:5,34,38,43,49,52)])

#logistic regression, similar to a t-test for count data
fit5=glmer(cbind(num_CO_3,num_NCO_3)~(1|repID)+Treatment*Day,data=genotyping,
         family=binomial(link="logit"), control=glmerControl(optimizer="bobyqa"))
summary(fit5)

I3_model=Anova(fit5,test="Chisq")
write.csv(I3_model,"SNP_i3_model.csv")
fit_contrast_i3 <- emmeans::emmeans(fit5, "Treatment", by="Day", mode="kenward-roger")
fit_contr_i3 <- contrast(fit_contrast_i3, method="trt.vs.ctrl")

pheno_contr_i3 <- as.data.frame(summary(fit_contr_i3))
pheno_contr_i3
write.csv(pheno_contr_i3,"SNP_interval3_posthoc.csv")

i3_or=exp(pheno_contr_i3$estimate)
i3_error=(pheno_contr_i3$SE)

#i4 region
#logistic regression, similar to a t-test for count data
fit7=glmer(cbind(num_CO_4,num_NCO_4)~(1|repID)+Treatment*Day,data=genotyping,
         family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa"))

I4_model=Anova(fit7,test="Chisq")
write.csv(I4_model,"SNP_i4_model.csv")
fit_contrast_i4 <- emmeans::emmeans(fit7, "Treatment", by="Day", mode="kenward-roger")
fit_contr_i4 <- contrast(fit_contrast_i4, method="trt.vs.ctrl")

pheno_contr_i4 <- as.data.frame(summary(fit_contr_i4))
pheno_contr_i4
write.csv(pheno_contr_i4,"SNP_interval4_posthoc.csv")

i4_or=exp(pheno_contr_i4$estimate)
i4_error=(pheno_contr_i4$SE)


#i5 region
#logistic regression, similar to a t-test for count data
fit9=glmer(cbind(num_CO_5,num_NCO_5)~(1|repID)+Treatment*Day,data=genotyping,
         family=binomial(link="logit"))
summary(fit9)

I5_model=Anova(fit9,test="Chisq")
write.csv(I5_model,"SNP_i5_model.csv")
fit_contrast_i5 <- emmeans::emmeans(fit9, "Treatment", by="Day", mode="kenward-roger")
fit_contr_i5 <- contrast(fit_contrast_i5, method="trt.vs.ctrl")

pheno_contr_i5 <- as.data.frame(summary(fit_contr_i5))
pheno_contr_i5
write.csv(pheno_contr_i5,"SNP_interval5_posthoc.csv")

i5_or=exp(pheno_contr_i5$estimate)
i5_error=(pheno_contr_i5$SE)

#setup data frame
odds_ratios=cbind(i1_or,i2_or,i3_or,i4_or,i5_or)
colnames(odds_ratios)=c("i1","i2","i3","i4","i5")
rownames(odds_ratios)=c("1-2","3-4","5-6","7-8","9-10")

#melt dataframe into long form
odds_ratios=melt(odds_ratios)
colnames(odds_ratios)=c("Day","Interval","OR")
odds_ratios

#add in standard error
SE=c(i1_error,i2_error,i3_error,i4_error,i5_error)
x=cbind(odds_ratios,SE)
x

#convert p-values to stars for plot
i1_sig=ifelse(pheno_contr_i1$p.value<0.001,"***",ifelse(pheno_contr_i1$p.value<0.01,"**",ifelse(pheno_contr_i1$p.value<0.05,"*","")))
i2_sig=ifelse(pheno_contr_i2$p.value<0.001,"***",ifelse(pheno_contr_i2$p.value<0.01,"**",ifelse(pheno_contr_i2$p.value<0.05,"*","")))
i3_sig=ifelse(pheno_contr_i3$p.value<0.001,"***",ifelse(pheno_contr_i3$p.value<0.01,"**",ifelse(pheno_contr_i3$p.value<0.05,"*","")))
i4_sig=ifelse(pheno_contr_i4$p.value<0.001,"***",ifelse(pheno_contr_i4$p.value<0.01,"**",ifelse(pheno_contr_i4$p.value<0.05,"*","")))
i5_sig=ifelse(pheno_contr_i5$p.value<0.001,"***",ifelse(pheno_contr_i5$p.value<0.01,"**",ifelse(pheno_contr_i5$p.value<0.05,"*","")))

#add significance to table
sig=c(i1_sig,i2_sig,i3_sig,i4_sig,i5_sig)
y=cbind(x,sig)
odds_ratios=y
odds_ratios

#Odds ratio plot
pdf("SNP_genotyping_odds.pdf")

odds_figure_SNP=ggplot(aes(y=OR,x=Day, col=Interval,group=Interval,shape=Interval),data=odds_ratios)+scale_colour_manual(values=c("black","black","black","#f1a340","#f1a340"))+ 
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(-1.25,8)+
  scale_x_discrete(name="Days post-mating",labels = c("1-2","3-4","5-6","7-8","9-10"))+geom_errorbar(aes(ymin=OR-SE,ymax=OR+SE))+ geom_line(aes(linetype=Interval))+ #scale_linetype_manual("", values=c(1,2,1,2,3))+ 
  annotate(geom="text", x=1, y=7.6, label=i1_sig[1],color="black",size=10)+annotate(geom="text", x=2, y=7.6, label=i1_sig[2],color="black",size=10)+annotate(geom="text", x=3, y=7.6, label=i1_sig[3],color="black",size=10)+annotate(geom="text", x=4, y=7.6, label=i1_sig[4],color="black",size=10)+annotate(geom="text", x=5, y=7.6, label=i1_sig[5],color="black",size=10)+
  annotate(geom="text", x=1, y=7.7, label=i2_sig[1],color="black",size=10)+annotate(geom="text", x=2, y=7.7, label=i2_sig[2],color="black",size=10)+annotate(geom="text", x=3, y=7.7, label=i2_sig[3],color="black",size=10)+annotate(geom="text", x=4, y=7.7, label=i2_sig[4],color="black",size=10)+annotate(geom="text", x=5, y=7.7, label=i2_sig[5],color="black",size=10)+
  annotate(geom="text", x=1, y=7.8, label=i3_sig[1],color="black",size=10)+annotate(geom="text", x=2, y=7.8, label=i3_sig[2],color="black",size=10)+annotate(geom="text", x=3, y=7.8, label=i3_sig[3],color="black",size=10)+annotate(geom="text", x=4, y=7.8, label=i3_sig[4],color="black",size=10)+annotate(geom="text", x=5, y=7.8, label=i3_sig[5],color="black",size=10)+
annotate(geom="text", x=1, y=7.9, label=i4_sig[1],color="#f1a340",size=10)+annotate(geom="text", x=2, y=7.9, label=i4_sig[2],color="#f1a340",size=10)+annotate(geom="text", x=3, y=7.9, label=i4_sig[3],color="#f1a340",size=10)+annotate(geom="text", x=4, y=7.9, label=i4_sig[4],color="#f1a340",size=10)+annotate(geom="text", x=5, y=7.9, label=i4_sig[5],color="#f1a340",size=10)+
  annotate(geom="text", x=1, y=8, label=i5_sig[1],color="#f1a340",size=10)+annotate(geom="text", x=2, y=8, label=i5_sig[2],color="#f1a340",size=10)+annotate(geom="text", x=3, y=8, label=i5_sig[3],color="#f1a340",size=10)+annotate(geom="text", x=4, y=8, label=i5_sig[4],color="#f1a340",size=10)+annotate(geom="text", x=5, y=8, label=i5_sig[5],color="#f1a340",size=10)
  
  odds_figure_SNP
dev.off()
