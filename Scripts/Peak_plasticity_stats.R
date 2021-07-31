#This is the master code for the paper 
#"Refining the timing of recombination rate plasticity in response to temperature in Drosophila pseudoobscura"
#written by Ulku Huma Altindag

#It is split into 18 code sections which corresponds to the each sequential experiment in this study
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


# Pilot Experiment 1 (exp1) -------------------------------------------------------------

#load the datasets
exp1=read.csv("Exp1_rawdata.csv", header= TRUE)
exp1_backcross=read.csv("exp1_backcross.csv", header=T)


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

yonly=sco_count=sum(exp1$SCO1, na.rm = TRUE)
stonly=sum(exp1$SCO2, na.rm = TRUE)
wild=sum(exp1$NCO1,na.rm=T)
mutant=sum(exp1$NCO2, na.rm = TRUE)

exp1_haplotypes=rbind(yonly,stonly,wild,mutant)
write.csv(exp1_haplotypes,"exp1_haplotypes.csv")

rownames(exp1_haplotypes)=c("y+", "+st", "++","yst")
colnames(exp1_haplotypes)=("number of progeny")

exp1_pvalues_4_haplotype_analysistotal=cbind(binom.test(c(exp1_haplotypes[1,1],exp1_haplotypes[2,1]),p=0.5)[3],
                                        binom.test(c(exp1_haplotypes[3,1],exp1_haplotypes[4,1]),p=0.5)[3])
                                        
colnames(exp1_pvalues_4_haplotype_analysistotal)=c("y+ and +st","++ and yst")

write.csv(exp1_pvalues_4_haplotype_analysistotal,"exp1_pvalues_4_haplotype_analysistotal.csv")


sco_count=sum(exp1$SCO1, na.rm = TRUE)+sum(exp1$SCO2, na.rm = TRUE)
nco_count=sum(exp1$NCO1+exp1$NCO2, na.rm = TRUE)
num_samples=sum(nco_count+sco_count, na.rm = TRUE)

#rate between yellow and vermillion

(sum(exp1$SCO1, na.rm = TRUE)+sum(exp1$SCO2, na.rm = TRUE))/num_samples


#merge with treatment data
colnames(exp1_backcross)
colnames(exp1_backcross)[1] <- "?..Vial.Number"
exp1_merged <- merge(exp1, exp1_backcross, by.x = "Vial..", by.y = "?..Vial.Number", all=T)
exp1_merged = na.omit(exp1_merged)
exp1_merged=na.omit(exp1_merged)

dataset=summaryBy(SCO1+SCO2+NCO1+NCO2~F1.Vial+new_day+Treatment,data=exp1_merged, FUN=sum,na.rm=T)

#continuing the haplotype analysis in order to calculate them with repsect to temperature and day
yonly=tapply(dataset$SCO1.sum,list(dataset$Treatment,dataset$new_day),sum,na.rm=T)
stonly=tapply(dataset$SCO2.sum,list(dataset$Treatment,dataset$new_day),sum,na.rm=T)
wild=tapply(dataset$NCO1.sum,list(dataset$Treatment,dataset$new_day),sum,na.rm=T)
mutant=tapply(dataset$NCO2.sum,list(dataset$Treatment,dataset$new_day),sum,na.rm=T)

exp1_haplotype_byday_and_treatment=rbind(yonly,stonly,wild,mutant)
#20C y+ and +st
pvaluesfor_18C_SCO_yonly_stonly=cbind(binom.test(c(yonly[1,1],stonly[1,1]),p=0.5)[3], 
                                      binom.test(c(yonly[1,2],stonly[1,2]),p=0.5)[3],
                                      binom.test(c(yonly[1,3],stonly[1,3]),p=0.5)[3], 
                                      binom.test(c(yonly[1,4],stonly[1,4]),p=0.5)[3], 
                                      binom.test(c(yonly[1,5],stonly[1,5]),p=0.5)[3], 
                                      binom.test(c(yonly[1,6],stonly[1,6]),p=0.5)[3]) 
colnames(pvaluesfor_18C_SCO_yonly_stonly)=c("1-2","3-4","5-6","7-9","10-12","13-15")
#20C ++ and yst
pvaluesfor_18C_NCO_wild_mutant=cbind(binom.test(c(wild[1,1],mutant[1,1]),p=0.5)[3], # p-value = 0.474
                                     binom.test(c(wild[1,2],mutant[1,2]),p=0.5)[3], # p-value = 0.9419
                                     binom.test(c(wild[1,3],mutant[1,3]),p=0.5)[3], # p-value = 0.7315
                                     binom.test(c(wild[1,4],mutant[1,4]),p=0.5)[3], # p-value = 0.00579 **
                                     binom.test(c(wild[1,5],mutant[1,5]),p=0.5)[3], # p-value = 0.03797 *
                                     binom.test(c(wild[1,6],mutant[1,6]),p=0.5)[3]) # p-value = 2.942e-08 ***
colnames(pvaluesfor_18C_NCO_wild_mutant)=c("1-2","3-4","5-6","7-9","10-12","13-15")
#25C y+ and +st
pvaluesfor_24C_SCO_yonly_stonly=cbind(binom.test(c(yonly[2,1],stonly[2,1]),p=0.5)[3] ,# p-value = 0.19
                                      binom.test(c(yonly[2,2],stonly[2,2]),p=0.5)[3] ,# p-value = 0.05176
                                      binom.test(c(yonly[2,3],stonly[2,3]),p=0.5)[3] ,# p-value = 0.5485
                                      binom.test(c(yonly[2,4],stonly[2,4]),p=0.5)[3] ,# p-value = 2.49e-10 ***
                                      binom.test(c(yonly[2,5],stonly[2,5]),p=0.5)[3] ,# p-value = 0.1055 
                                      binom.test(c(yonly[2,6],stonly[2,6]),p=0.5)[3] )# p-value = 3.761e-05 ***
colnames(pvaluesfor_24C_SCO_yonly_stonly)=c("1-2","3-4","5-6","7-9","10-12","13-15")
#25C ++ and yst
pvaluesfor_24C_NCO_wild_mutant=cbind(binom.test(c(wild[2,1],mutant[2,1]),p=0.5)[3], # p-value = 0.8626
                                     binom.test(c(wild[2,2],mutant[2,2]),p=0.5)[3], # p-value = 1
                                     binom.test(c(wild[2,3],mutant[2,3]),p=0.5)[3],# p-value = 0.02208 *
                                     binom.test(c(wild[2,4],mutant[2,4]),p=0.5)[3], # p-value = 0.6691
                                     binom.test(c(wild[2,5],mutant[2,5]),p=0.5)[3], # p-value = 0.0006782 ***
                                     binom.test(c(wild[2,6],mutant[2,6]),p=0.5)[3]) # p-value = 0.0357 *
colnames(pvaluesfor_24C_NCO_wild_mutant)=c("1-2","3-4","5-6","7-9","10-12","13-15")
pvalues_4_haplotype_analysis=rbind(pvaluesfor_18C_SCO_yonly_stonly,pvaluesfor_18C_NCO_wild_mutant,pvaluesfor_24C_SCO_yonly_stonly,pvaluesfor_24C_NCO_wild_mutant)
rownames(pvalues_4_haplotype_analysis)=c("pvaluesfor_18C_SCO_yonly_stonly","pvaluesfor_18C_NCO_wild_mutant","pvaluesfor_24C_SCO_yonly_stonly","pvaluesfor_24C_NCO_wild_mutant")

write.csv(exp1_haplotype_byday_and_treatment,"exp1_haplotype_byday_and_treatment.csv")
write.csv(pvalues_4_haplotype_analysis,"exp1_pvalues_4_haplotype_analysis.csv")

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
#line plot for fecundity analysis
#Fecund_figure_exp1=ggplot(aes(y=fecundity,x=new_day, col=Treatment,label=Num_moms),data=dataset)+ylab("# Progeny per mom")+theme_base()+

#Fecund_figure_exp1=Fecund_figure_exp1+scale_colour_manual(values=c("blue", "red"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("1-2","3-4","5-6","7-9","10-12","13-15"))+ylim(0,50)
#Fecund_figure_exp1=Fecund_figure_exp1+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)+
 # annotate(geom="text", x=1, y=25, label=sig[1],size=10)+annotate(geom="text", x=2, y=20, label=sig[2],size=10)+annotate(geom="text", x=3, y=25, label=sig[3],size=10)+annotate(geom="text", x=4, y=25, label=sig[4],size=10)
#Fecund_figure_exp1

#boxplot for fecundity analysis
Fecund_figure_exp1=ggplot(dataset, aes(x=new_day, y=fecundity, col=Treatment)) + theme_base()+ylab("# Progeny per mom")+ggtitle("Experiment 1")+
  geom_boxplot()+scale_x_discrete(name="Days", labels=c("1-2","3-4","5-6","7-9","10-12","13-15"))+theme(axis.text.x = element_text(angle = 45))+
   annotate(geom="text", x=2, y=20, label=sig[2],size=10)+annotate(geom="text", x=3, y=25, label=sig[3],size=10)+annotate(geom="text", x=4, y=25, label=sig[4],size=10)+
    scale_color_manual(values = c("blue","red"))
Fecund_figure_exp1

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

write.csv(dataset,"experiment1_combined dataset.csv")
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

odds_figure_exp1=ggplot(aes(y=vy_or,x=day,group=1),data=odds_ratios)+scale_colour_manual(values=c("black"))+ggtitle("Experiment 1")+
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.5,2.6)+theme(axis.text.x = element_text(angle = 45))+
  scale_x_discrete(name="Days post-mating",labels=c("1-2","3-4","5-6","7-9","10-12","13-15"))+geom_line()+geom_errorbar(aes(ymin=vy_or-vy_error,ymax=vy_or+vy_error))+
  annotate(geom="text", x=1, y=1.75, label=vy_sig[1],color="black",size=10)+annotate(geom="text", x=2, y=1.75, label=vy_sig[2],color="black",size=10)+annotate(geom="text", x=3, y=1.75, label=vy_sig[3],color="black",size=10)+annotate(geom="text", x=4, y=1.75, label=vy_sig[4],color="black",size=10)+annotate(geom="text", x=5, y=1.75, label=vy_sig[5],color="black",size=10)+annotate(geom="text", x=6, y=1.75, label=vy_sig[6],color="black",size=10)
odds_figure_exp1
dev.off()

#Summarization of the data

length(unique(dataset$F1.Vial[dataset$Treatment=="18"]))
length(unique(dataset$F1.Vial[dataset$Treatment=="24"]))

median(dataset$Num_moms[dataset$Treatment=="18"])
median(dataset$Num_moms[dataset$Treatment=="24"])

sum(dataset$num_CO[dataset$Treatment=="18"])+sum(dataset$num_NCO[dataset$Treatment=="18"])
sum(dataset$num_NCO[dataset$Treatment=="24"])+sum(dataset$num_CO[dataset$Treatment=="24"])

# Pilot Experiment2  (exp2)-------------------------------------------------------------

#load the datasets
exp2=read.csv("Exp2_rawdata.csv", header= TRUE, na.strings="-")
exp2_backcross=read.csv("Exp2_bc_setup.csv", header=T)
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

yonly=sco_count=sum(exp2$SCO1, na.rm = TRUE)
stonly=sum(exp2$SCO2, na.rm = TRUE)
wild=sum(exp2$NCO1,na.rm=T)
mutant=sum(exp2$NCO2, na.rm = TRUE)

exp2_haplotypes=rbind(yonly,stonly,wild,mutant)


rownames(exp2_haplotypes)=c("y+", "+st", "++","yst")
colnames(exp2_haplotypes)=("number of progeny")
write.csv(exp2_haplotypes,"exp2_haplotypes.csv")

exp2_pvalues_4_haplotype_analysistotal=cbind(binom.test(c(exp2_haplotypes[1,1],exp2_haplotypes[2,1]),p=0.5)[3],
                                             binom.test(c(exp2_haplotypes[3,1],exp2_haplotypes[4,1]),p=0.5)[3])

colnames(exp2_pvalues_4_haplotype_analysistotal)=c("y+ and +st","++ and yst")

write.csv(exp2_pvalues_4_haplotype_analysistotal,"exp2_pvalues_4_haplotype_analysistotal.csv")


sco_count=sum(exp2$SCO1+exp2$SCO2, na.rm = TRUE)
nco_count=sum(exp2$NCO1+exp2$NCO2, na.rm = TRUE)
num_samples=sum(nco_count+sco_count, na.rm = TRUE)

#merge with treatment data
exp2_merged <- merge(exp2, exp2_backcross, by.x = "Vial.number", by.y = "Vial.number", all=T)
exp2_merged = na.omit(exp2_merged)

dataset=summaryBy(SCO1+SCO2+NCO1+NCO2~F1.Vial+Day..letter.of.vial.+Treatment,data=exp2_merged, FUN=sum,na.rm=T)

#continuing the haplotype analysis in order to calculate them with repsect to temperature and day
yonly=tapply(dataset$SCO1.sum,list(dataset$Treatment,dataset$Day..letter.of.vial.),sum,na.rm=T)
stonly=tapply(dataset$SCO2.sum,list(dataset$Treatment,dataset$Day..letter.of.vial.),sum,na.rm=T)
wild=tapply(dataset$NCO1.sum,list(dataset$Treatment,dataset$Day..letter.of.vial.),sum,na.rm=T)
mutant=tapply(dataset$NCO2.sum,list(dataset$Treatment,dataset$Day..letter.of.vial.),sum,na.rm=T)

exp2_haplotype_byday_and_treatment=rbind(yonly,stonly,wild,mutant)
#20C y+ and +st
pvaluesfor_20C_SCO_yonly_stonly=cbind(binom.test(c(yonly[1,1],stonly[1,1]),p=0.5)[3], 
                                      binom.test(c(yonly[1,2],stonly[1,2]),p=0.5)[3],
                                      binom.test(c(yonly[1,3],stonly[1,3]),p=0.5)[3], 
                                      binom.test(c(yonly[1,4],stonly[1,4]),p=0.5)[3])
colnames(pvaluesfor_20C_SCO_yonly_stonly)=c("1-5","6-10","11-15","16-20")
#20C ++ and yst
pvaluesfor_20C_NCO_wild_mutant=cbind(binom.test(c(wild[1,1],mutant[1,1]),p=0.5)[3], # p-value = 0.474
                                     binom.test(c(wild[1,2],mutant[1,2]),p=0.5)[3], # p-value = 0.9419
                                     binom.test(c(wild[1,3],mutant[1,3]),p=0.5)[3], # p-value = 0.7315
                                     binom.test(c(wild[1,4],mutant[1,4]),p=0.5)[3])
colnames(pvaluesfor_20C_NCO_wild_mutant)=c("1-5","6-10","11-15","16-20")
#25C y+ and +st
pvaluesfor_25C_SCO_yonly_stonly=cbind(binom.test(c(yonly[2,1],stonly[2,1]),p=0.5)[3] ,# p-value = 0.19
                                      binom.test(c(yonly[2,2],stonly[2,2]),p=0.5)[3] ,# p-value = 0.05176
                                      binom.test(c(yonly[2,3],stonly[2,3]),p=0.5)[3] ,# p-value = 0.5485
                                      binom.test(c(yonly[2,4],stonly[2,4]),p=0.5)[3])
colnames(pvaluesfor_25C_SCO_yonly_stonly)=c("1-5","6-10","11-15","16-20")
#25C ++ and yst
pvaluesfor_25C_NCO_wild_mutant=cbind(binom.test(c(wild[2,1],mutant[2,1]),p=0.5)[3], # p-value = 0.8626
                                     binom.test(c(wild[2,2],mutant[2,2]),p=0.5)[3], # p-value = 1
                                     binom.test(c(wild[2,3],mutant[2,3]),p=0.5)[3],# p-value = 0.02208 *
                                     binom.test(c(wild[2,4],mutant[2,4]),p=0.5)[3])
colnames(pvaluesfor_25C_NCO_wild_mutant)=c("1-5","6-10","11-15","16-20")
pvalues_4_haplotype_analysis=rbind(pvaluesfor_20C_SCO_yonly_stonly,pvaluesfor_20C_NCO_wild_mutant,pvaluesfor_25C_SCO_yonly_stonly,pvaluesfor_25C_NCO_wild_mutant)
rownames(pvalues_4_haplotype_analysis)=c("pvaluesfor_20C_SCO_yonly_stonly","pvaluesfor_20C_NCO_wild_mutant","pvaluesfor_25C_SCO_yonly_stonly","pvaluesfor_25C_NCO_wild_mutant")

write.csv(exp2_haplotype_byday_and_treatment,"exp2_haplotype_byday_and_treatment.csv")
write.csv(pvalues_4_haplotype_analysis,"exp2_pvalues_4_haplotype_analysis.csv")


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
#Fecund_figure_exp2=ggplot(aes(y=fecundity,x=Day..letter.of.vial., col=Treatment,label=Num_moms),data=dataset2)+ylab("# Progeny per mom")+theme_base()

#Fecund_figure_exp2=Fecund_figure_exp2+scale_colour_manual(values=c("blue", "red"))+geom_point(size=1.5)+scale_x_discrete(name="Day",labels=c("1-5","6-10","11-15","16-20"))+ylim(0,50)+theme(legend.title = element_text(size = 10),legend.text = element_text(size = 8))
#Fecund_figure_exp2=Fecund_figure_exp2+stat_summary(fun = median, geom="line",aes(group=Treatment),size=1.5)+geom_text(check_overlap = F,hjust = 0, nudge_x = 0.1,angle=45,size=2.5)+
#  annotate(geom="text", x=1, y=40, label=sig[1],size=5)+annotate(geom="text", x=2, y=40, label=sig[2],size=5)+annotate(geom="text", x=3, y=40, label=sig[3],size=5)+annotate(geom="text", x=4, y=40, label=sig[4],size=5)
#Fecund_figure_exp2

#boxplot for fecundity analysis
Fecund_figure_exp2=ggplot(dataset2, aes(x=Day..letter.of.vial., y=fecundity, col=Treatment)) + theme_base()+ylab("# Progeny per mom")+ggtitle("Experiment 2")+theme(axis.text.x = element_text(angle = 45))+
  geom_boxplot()+scale_x_discrete(name="Days", labels=c("1-5","6-10","11-15","16-20"))+
  annotate(geom="text", x=1, y=40, label=sig[1],size=5)+annotate(geom="text", x=2, y=40, label=sig[2],size=5)+annotate(geom="text", x=3, y=40, label=sig[3],size=5)+annotate(geom="text", x=4, y=40, label=sig[4],size=5)+
  scale_color_manual(values = c("blue","red"))
Fecund_figure_exp2

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

write.csv(dataset2,"experiment2_combined dataset.csv")

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
odds_figure_exp2=ggplot(aes(y=exp2_or,x=Day,group=1),data=odds_ratios)+scale_colour_manual(values=c("black"))+ggtitle("Experiment 2")+
  geom_line()+geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.5,2.5)+
  scale_x_discrete(name="Days post-mating",labels=c("1-5","6-10","11-15","16-20"))+geom_line()+geom_errorbar(aes(ymin=exp2_or-exp2_error,ymax=exp2_or+exp2_error))+theme(axis.text.x = element_text(angle = 45))+
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


# Experiment 1 (Exp3) -------------------------------------------------------------

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

##CO groups were defined for total progeny combined

exp3$SCO1=exp3$yellow.only+exp3$yellow.only.1
exp3$SCO2=exp3$vermillion.only+exp3$vermillion.only.1
exp3$NCO1=exp3$wildtype+exp3$wildtype.1
exp3$NCO2=exp3$yellow.vermillion+exp3$yellow.vermillion.1

sum(exp3$SCO1,na.rm = T)
sum(exp3$SCO2,na.rm = T)
sum(exp3$NCO1,na.rm = T)
sum(exp3$NCO2,na.rm = T)

haplotypes_exp3=rbind(sum(exp3$SCO1,na.rm = T),
                      sum(exp3$SCO2,na.rm = T),
                      sum(exp3$NCO1,na.rm = T),
                      sum(exp3$NCO2,na.rm = T))

rownames(haplotypes_exp3)=c("y+", "+st", "++","yst")
colnames(haplotypes_exp3)=("number of progeny")

exp3$SCO1male=exp3$yellow.only
exp3$SCO2male=exp3$vermillion.only
exp3$NCO1male=exp3$wildtype
exp3$NCO2male=exp3$yellow.vermillion

sum(exp3$SCO1male,na.rm = T)
sum(exp3$SCO2male,na.rm = T)
sum(exp3$NCO1male,na.rm = T)
sum(exp3$NCO2male,na.rm = T)

#male haplotypes

haplotypes_exp3_male=rbind(sum(exp3$SCO1male,na.rm = T),
                           sum(exp3$SCO2male,na.rm = T),
                           sum(exp3$NCO1male,na.rm = T),
                           sum(exp3$NCO2male,na.rm = T))
rownames(haplotypes_exp3_male)=c("y+", "+st", "++","yst")
colnames(haplotypes_exp3_male)=("number of progeny")

write.csv(haplotypes_exp3_male,"haplotypes_exp3_male.csv")

#female haplotypes

exp3$SCO1female=exp3$yellow.only.1
exp3$SCO2female=exp3$vermillion.only.1
exp3$NCO1female=exp3$wildtype.1
exp3$NCO2female=exp3$yellow.vermillion.1
sum(exp3$SCO1female,na.rm = T)
sum(exp3$SCO2female,na.rm = T)
sum(exp3$NCO1female,na.rm = T)
sum(exp3$NCO2female,na.rm = T)

haplotypes_exp3_female=rbind(sum(exp3$SCO1female,na.rm = T),
                             sum(exp3$SCO2female,na.rm = T),
                             sum(exp3$NCO1female,na.rm = T),
                             sum(exp3$NCO2female,na.rm = T))
rownames(haplotypes_exp3_female)=c("y+", "+st", "++","yst")
colnames(haplotypes_exp3_female)=("number of progeny")

write.csv(haplotypes_exp3_female,"haplotypes_exp3_female.csv")

pvalues_4_haplotype_analysistotal=cbind(binom.test(c(haplotypes_exp3[1,1],haplotypes_exp3[2,1]),p=0.5)[3],#p-value=0.01388 *
                                   binom.test(c(haplotypes_exp3[3,1],haplotypes_exp3[4,1]),p=0.5)[3],#p-value = 0.002021 **
                                   binom.test(c(haplotypes_exp3_female[1,1],haplotypes_exp3_female[2,1]),p=0.5)[3],#p-value = 1.397e-06
                                   binom.test(c(haplotypes_exp3_female[3,1],haplotypes_exp3_female[4,1]),p=0.5)[3],#p-value = 2.2e-16
                                   binom.test(c(haplotypes_exp3_male[1,1],haplotypes_exp3_male[2,1]),p=0.5)[3],#p-value=0.1405
                                   binom.test(c(haplotypes_exp3_male[3,1],haplotypes_exp3_male[4,1]),p=0.5)[3])#p-value = 3.089e-09
                                   
colnames(pvalues_4_haplotype_analysistotal)=c("y+ and +st","++ and yst","female y+ and +st","female ++ and yst","male y+ and +st","male ++ and yst")

write.csv(haplotypes_exp3,"haplotypes_exp3.csv")
write.csv(pvalues_4_haplotype_analysistotal,"pvalues_4_haplotype_analysistotal.csv")


binom.test(2388,4608, p=0.5)
binom.test(1178,2371, p=0.5)
binom.test(2452,5126, p=0.5)
##Using the sanity check, we can see whether R interprets the results accurately.

sco_count=sum(exp3$SCO1, na.rm = TRUE)+sum(exp3$SCO2, na.rm = TRUE)
nco_count=sum(exp3$NCO1+exp3$NCO2, na.rm = TRUE)
num_samples=sum(nco_count+sco_count, na.rm = TRUE)

(sum(exp3$SCO1, na.rm = TRUE)+sum(exp3$SCO2, na.rm = TRUE))/num_samples

##merge with treatment data (This is important for summarization of the data)
exp3_merged <- merge(exp3, exp3_backcross, by.x = "Vial..", by.y = "Vial.number", all=T)
exp3_merged = na.omit(exp3_merged)
dataset=summaryBy(SCO1+SCO2+NCO1+NCO2~F1.Vial+Day..letter.of.vial.+Treatment,data=exp3_merged, FUN=sum,na.rm=T)

#continuing the haplotype analysis in order to calculate them with repsect to temperature and day
yonly=tapply(dataset$SCO1.sum,list(dataset$Treatment,dataset$Day..letter.of.vial.),sum,na.rm=T)
stonly=tapply(dataset$SCO2.sum,list(dataset$Treatment,dataset$Day..letter.of.vial.),sum,na.rm=T)
wild=tapply(dataset$NCO1.sum,list(dataset$Treatment,dataset$Day..letter.of.vial.),sum,na.rm=T)
mutant=tapply(dataset$NCO2.sum,list(dataset$Treatment,dataset$Day..letter.of.vial.),sum,na.rm=T)

exp3_haplotype_byday_and_treatment=rbind(yonly,stonly,wild,mutant)
#20C y+ and +st
pvaluesfor_20C_SCO_yonly_stonly=cbind(binom.test(c(yonly[1,1],stonly[1,1]),p=0.5)[3], # p-value = 0.7713
                                      binom.test(c(yonly[1,2],stonly[1,2]),p=0.5)[3], # p-value = 0.1769
                                      binom.test(c(yonly[1,3],stonly[1,3]),p=0.5)[3], # p-value = 0.8831
                                      binom.test(c(yonly[1,4],stonly[1,4]),p=0.5)[3], # p-value = 2.701e-05 ***
                                      binom.test(c(yonly[1,5],stonly[1,5]),p=0.5)[3], # p-value = 0.01135 *
                                      binom.test(c(yonly[1,6],stonly[1,6]),p=0.5)[3]) # p-value = 2.939e-06 ***
colnames(pvaluesfor_20C_SCO_yonly_stonly)=c("1-2","3-4","5-6","7-9","10-12","13-15")
#20C ++ and yst
pvaluesfor_20C_NCO_wild_mutant=cbind(binom.test(c(wild[1,1],mutant[1,1]),p=0.5)[3], # p-value = 0.474
                                      binom.test(c(wild[1,2],mutant[1,2]),p=0.5)[3], # p-value = 0.9419
                                      binom.test(c(wild[1,3],mutant[1,3]),p=0.5)[3], # p-value = 0.7315
                                      binom.test(c(wild[1,4],mutant[1,4]),p=0.5)[3], # p-value = 0.00579 **
                                      binom.test(c(wild[1,5],mutant[1,5]),p=0.5)[3], # p-value = 0.03797 *
                                      binom.test(c(wild[1,6],mutant[1,6]),p=0.5)[3]) # p-value = 2.942e-08 ***
colnames(pvaluesfor_20C_NCO_wild_mutant)=c("1-2","3-4","5-6","7-9","10-12","13-15")
#25C y+ and +st
pvaluesfor_25C_SCO_yonly_stonly=cbind(binom.test(c(yonly[2,1],stonly[2,1]),p=0.5)[3] ,# p-value = 0.19
                                     binom.test(c(yonly[2,2],stonly[2,2]),p=0.5)[3] ,# p-value = 0.05176
                                     binom.test(c(yonly[2,3],stonly[2,3]),p=0.5)[3] ,# p-value = 0.5485
                                     binom.test(c(yonly[2,4],stonly[2,4]),p=0.5)[3] ,# p-value = 2.49e-10 ***
                                     binom.test(c(yonly[2,5],stonly[2,5]),p=0.5)[3] ,# p-value = 0.1055 
                                     binom.test(c(yonly[2,6],stonly[2,6]),p=0.5)[3] )# p-value = 3.761e-05 ***
colnames(pvaluesfor_25C_SCO_yonly_stonly)=c("1-2","3-4","5-6","7-9","10-12","13-15")
#25C ++ and yst
pvaluesfor_25C_NCO_wild_mutant=cbind(binom.test(c(wild[2,1],mutant[2,1]),p=0.5)[3], # p-value = 0.8626
                                      binom.test(c(wild[2,2],mutant[2,2]),p=0.5)[3], # p-value = 1
                                      binom.test(c(wild[2,3],mutant[2,3]),p=0.5)[3],# p-value = 0.02208 *
                                      binom.test(c(wild[2,4],mutant[2,4]),p=0.5)[3], # p-value = 0.6691
                                      binom.test(c(wild[2,5],mutant[2,5]),p=0.5)[3], # p-value = 0.0006782 ***
                                      binom.test(c(wild[2,6],mutant[2,6]),p=0.5)[3]) # p-value = 0.0357 *
colnames(pvaluesfor_25C_NCO_wild_mutant)=c("1-2","3-4","5-6","7-9","10-12","13-15")
pvalues_4_haplotype_analysis=rbind(pvaluesfor_20C_SCO_yonly_stonly,pvaluesfor_20C_NCO_wild_mutant,pvaluesfor_25C_SCO_yonly_stonly,pvaluesfor_25C_NCO_wild_mutant)
rownames(pvalues_4_haplotype_analysis)=c("pvaluesfor_20C_SCO_yonly_stonly","pvaluesfor_20C_NCO_wild_mutant","pvaluesfor_25C_SCO_yonly_stonly","pvaluesfor_25C_NCO_wild_mutant")

write.csv(exp3_haplotype_byday_and_treatment,"exp3_haplotype_byday_and_treatment.csv")
write.csv(pvalues_4_haplotype_analysis,"exp3_pvalues_4_haplotype_analysis.csv")
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
#Fecund_figure_exp3=ggplot(aes(y=fecundity,x=Day..letter.of.vial., col=Treatment,label=Num_moms),data=dataset)+ylab("# Progeny per mom")+theme_base()

#Fecund_figure_exp3=Fecund_figure_exp3+scale_colour_manual(values=c("blue", "red"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("1-2","3-4","5-6","7-9","10-12","13-15"))+ylim(0,50)
#Fecund_figure_exp3=Fecund_figure_exp3+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
#Fecund_figure_exp3

##Statistical analysis
###the dataset necessary for the statistics is the cleaned up version of the summary data
dataset2=read.csv("Experiment3_combined dataset.csv",header=T,stringsAsFactors = F)
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



#boxplot for fecundity analysis
Fecund_figure_exp3=ggplot(dataset, aes(x=Day..letter.of.vial., y=fecundity, col=Treatment)) + theme_base()+ylab("# Progeny per mom")+ggtitle("Experiment 3")+theme(axis.text.x = element_text(angle = 45))+
  geom_boxplot()+scale_x_discrete(name="Days", labels=c("1-2","3-4","5-6","7-9","10-12","13-15"))+
  #annotate(geom="text", x=1, y=40, label=sig[1],size=5)+annotate(geom="text", x=2, y=40, label=sig[2],size=5)+annotate(geom="text", x=3, y=40, label=sig[3],size=5)+annotate(geom="text", x=4, y=40, label=sig[4],size=5)+
  scale_color_manual(values = c("blue","red"))
Fecund_figure_exp3

dev.off()
##remove vials with few progeny
dataset=dataset[dataset$total_offspring>=10,]


##Calculate the recombination rate for each replicate
dataset$rec_rate_total=(dataset$SCO1.sum+dataset$SCO2.sum)/(dataset$NCO1.sum+dataset$NCO2.sum+dataset$SCO1.sum+dataset$SCO2.sum)
dataset$kosambi_rec_rate=((2.71^(4*dataset$rec_rate_total)-1)/2*(2.71^(-4*dataset$rec_rate_total)+1))/10



write.csv(dataset,"experiment3_combined dataset.csv")
for (n in 1:length(dataset$F1.Vial)) {
  dataset$NCO_skew[n]=min(dataset$NCO1.sum[n],dataset$NCO2.sum[n])/max(dataset$NCO1.sum[n],dataset$NCO2.sum[n])
}

for (n in 1:length(dataset$F1.Vial)) {
  dataset$SCO1_skew[n]=min(dataset$SCO1.sum[n],dataset$SCO2.sum[n])/max(dataset$SCO1.sum[n],dataset$SCO2.sum[n])
}



skew=dataset[, c(3,13,14)]

skew2=melt(skew, id="Treatment")


pdf("skew_exp3.pdf")
skew_figure_exp3=ggplot(skew2, aes(x = Treatment, y = value, color = variable)) +  # ggplot function
  geom_boxplot()+theme_base()+ylab("haplotype bias")+xlab("Treatment")

skew_figure_exp3
dev.off()

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

g=subset(dataset2, dataset2$Day..letter.of.vial.=="D")
g
###get mean recombination
tapply(g$rec_rate_exp3,g$Treatment,mean)
tapply(g$rec_rate_exp3,g$Day,mean)
tapply(g$total_offspring,g$Treatment,sum)

h=subset(dataset2, dataset2$Day..letter.of.vial.=="F")
h

tapply(h$rec_rate_exp3,h$Treatment,mean)
tapply(h$rec_rate_exp3,h$Day,mean)
tapply(h$total_offspring,h$Treatment,sum)

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


##Recomb_figure_total
pdf("Experiment3_recombination.pdf")
Recomb_figure=ggplot(aes(y=rec_rate_total,x=Day..letter.of.vial., col=Treatment,label=total_offspring),data=dataset)+ylab("% recombination")+ggtitle("Total Recombination rate vs. Days post-mating")+theme_base()
Recomb_figure=Recomb_figure+scale_colour_manual(values=c("blue", "red"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("1-2","3-4","5-6","7-9","10-12","13-15"))
Recomb_figure=Recomb_figure+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)
Recomb_figure=Recomb_figure+ggtitle("Total Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3) +ylim(0,1.2)
Recomb_figure
dev.off()


pdf("Exp3_boxplot.pdf")

Recomb_figure_yst=ggplot(aes(y=rec_rate_total,x=Day..letter.of.vial., col=Treatment),data = dataset)+scale_colour_manual(values=c("blue","red"))+geom_boxplot()+ylab("% recombination")+theme_base()+scale_x_discrete(name="Days post-mating",labels=c("1-2","3-4","5-6","7-9","10-12","13-15"))+
  ggtitle("sd-y Interval")+ylim(0,1)+
  annotate(geom="text", x=1, y=0.58, label=vy_sig[1],size=10)+annotate(geom="text", x=2, y=0.58, label=vy_sig[2],size=10)+annotate(geom="text", x=3, y=0.58, label=vy_sig[3],size=10)+annotate(geom="text", x=4, y=0.78, label=vy_sig[4],size=10)+annotate(geom="text", x=5, y=0.78,label=vy_sig[5],size=10)+annotate(geom="text",x=6,y=0.78,label=vy_sig[6],size=10)
Recomb_figure_yst

dev.off()


pdf("Experiment3.odds_new.pdf")
odds_figure_exp3=ggplot(aes(y=vy_or,x=day,group=1),data=odds_ratios)+scale_colour_manual(values=c("black"))+ggtitle("Experiment 3")+
  geom_point(size=3)+geom_line(size=2)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.5,2.6)+theme(axis.text.x = element_text(angle = 45))+
  scale_x_discrete(name="Days post-mating",labels=c("1-2","3-4","5-6","7-9","10-12","13-15"))+geom_line(show.legend = TRUE)+geom_errorbar(aes(ymin=vy_or-vy_error,ymax=vy_or+vy_error))+
  annotate(geom="text", x=1, y=1.75, label=vy_sig[1],color="black",size=10)+annotate(geom="text", x=2, y=1.75, label=vy_sig[2],color="black",size=10)+annotate(geom="text", x=3, y=1.75, label=vy_sig[3],color="black",size=8)+annotate(geom="text", x=4, y=1.75, label=vy_sig[4],color="black",size=8)+annotate(geom="text", x=5, y=1.75, label=vy_sig[5],color="black",size=10)+annotate(geom="text", x=5.75, y=1.75, label=vy_sig[6],color="black",size=8)
odds_figure_exp3
dev.off()

sum(dataset2$total_offspring[dataset2$Treatment=="20°"])
sum(dataset2$total_offspring[dataset2$Treatment=="25°"])

# Experiment 2 (Exp4) -------------------------------------------------------------

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

#Define Crossovers
exp4$gamete_class=ifelse(exp4$sd==exp4$y & exp4$y==exp4$se & exp4$y==1,"gt1a",
                         ifelse(exp4$sd==exp4$y & exp4$y==exp4$se & exp4$y==0,"gt1b",
                                ifelse(exp4$sd!=exp4$y & exp4$y==exp4$se & exp4$y==1,"gt2a",
                                       ifelse(exp4$sd!=exp4$y & exp4$y==exp4$se & exp4$y==0,"gt2b",
                                              ifelse(exp4$sd==exp4$y & exp4$y!=exp4$se & exp4$y==1,"gt3a",
                                                     ifelse(exp4$sd==exp4$y & exp4$y!=exp4$se & exp4$y==0,"gt3b",
                                                            ifelse(exp4$sd!=exp4$y & exp4$y!=exp4$se & exp4$y==1,"gt4a",
                                                                   ifelse(exp4$sd!=exp4$y & exp4$y!=exp4$se & exp4$y==0,"gt4b",
                                                                          "error"))))))))
#Sanity check, there should not be any error!
exp4[exp4$gamete_class=="error",]

#add columns for counting
exp4$gt1a_111=ifelse(exp4$gamete_class=="gt1a",exp4$numbMales,0)
exp4$gt1b_000=ifelse(exp4$gamete_class=="gt1b",exp4$numbMales,0)
exp4$gt2a_011=ifelse(exp4$gamete_class=="gt2a",exp4$numbMales,0)
exp4$gt2b_100=ifelse(exp4$gamete_class=="gt2b",exp4$numbMales,0)
exp4$gt3a_110=ifelse(exp4$gamete_class=="gt3a",exp4$numbMales,0)
exp4$gt3b_001=ifelse(exp4$gamete_class=="gt3b",exp4$numbMales,0)
exp4$gt4a_010=ifelse(exp4$gamete_class=="gt4a",exp4$numbMales,0)
exp4$gt4b_101=ifelse(exp4$gamete_class=="gt4b",exp4$numbMales,0)

gt1a_111=sum(exp4$numbMales*exp4$gt1a_111)
gt1b_000=sum(exp4$numbMales*exp4$gt1b_000)
gt2a_011=sum(exp4$numbMales*exp4$gt2a_011)
gt2b_100=sum(exp4$numbMales*exp4$gt2b_100)
gt3a_110=sum(exp4$numbMales*exp4$gt3a_110)
gt3b_001=sum(exp4$numbMales*exp4$gt3b_001)
gt4a_010=sum(exp4$numbMales*exp4$gt4a_010)
gt4b_101=sum(exp4$numbMales*exp4$gt4b_101)

exp4_haplotype_types=cbind(gt1a_111,gt1b_000,gt2a_011,gt2b_100,gt3a_110,gt3b_001,gt4a_010,gt4b_101)

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

#summarize long form data, include the gamete classes in the summary
dataset=summaryBy(male+gt1a_111+gt1b_000+gt2a_011+gt2b_100+gt3a_110+gt3b_001+gt4a_010+gt4b_101+num_co+SCO_1+SCO_2+DCO~F1.Vial+Day+Treatment,data=exp4_merged, FUN=sum,na.rm=T)

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

#gamete analysis by treatment
gt111=tapply(dataset2$gt1a_111.sum,dataset2$Treatment,sum,na.rm=T)
gt000=tapply(dataset2$gt1b_000.sum,dataset2$Treatment,sum,na.rm=T)
gt011=tapply(dataset2$gt2a_011.sum,dataset2$Treatment,sum,na.rm=T)
gt100=tapply(dataset2$gt2b_100.sum,dataset2$Treatment,sum,na.rm=T)
gt110=tapply(dataset2$gt3a_110.sum,dataset2$Treatment,sum,na.rm=T)
gt001=tapply(dataset2$gt3b_001.sum,dataset2$Treatment,sum,na.rm=T)
gt010=tapply(dataset2$gt4a_010.sum,dataset2$Treatment,sum,na.rm=T)
gt101=tapply(dataset2$gt4b_101.sum,dataset2$Treatment,sum,na.rm=T)
i=rbind(gt111, gt000, gt011, gt100, gt110, gt001, gt010, gt101)
write.csv(i,"Exp4_gamete_class_by_treatment.csv")
#gamete analysis by day
gt111=tapply(dataset2$gt1a_111.sum,dataset2$Day,sum,na.rm=T)
gt000=tapply(dataset2$gt1b_000.sum,dataset2$Day,sum,na.rm=T)
gt011=tapply(dataset2$gt2a_011.sum,dataset2$Day,sum,na.rm=T)
gt100=tapply(dataset2$gt2b_100.sum,dataset2$Day,sum,na.rm=T)
gt110=tapply(dataset2$gt3a_110.sum,dataset2$Day,sum,na.rm=T)
gt001=tapply(dataset2$gt3b_001.sum,dataset2$Day,sum,na.rm=T)
gt010=tapply(dataset2$gt4a_010.sum,dataset2$Day,sum,na.rm=T)
gt101=tapply(dataset2$gt4b_101.sum,dataset2$Day,sum,na.rm=T)
s=rbind(gt111, gt000, gt011, gt100, gt110, gt001, gt010, gt101)
write.csv(s,"exp4_gamete_class_byday.csv")
gt111=tapply(dataset2$gt1a_111.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt000=tapply(dataset2$gt1b_000.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt011=tapply(dataset2$gt2a_011.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt100=tapply(dataset2$gt2b_100.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt110=tapply(dataset2$gt3a_110.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt001=tapply(dataset2$gt3b_001.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt010=tapply(dataset2$gt4a_010.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt101=tapply(dataset2$gt4b_101.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
t=rbind(gt111, gt000, gt011, gt100, gt110, gt001, gt010, gt101)
write.csv(t,"exp4_gamete_class_byday_and_treatment.csv")
exp4_pvalues_111_000_21=cbind(binom.test(c(gt111[1,1],gt000[1,1]),p=0.5)[3],
                              binom.test(c(gt111[1,2],gt000[1,2]),p=0.5)[3],
                              binom.test(c(gt111[1,3],gt000[1,3]),p=0.5)[3],
                              binom.test(c(gt111[1,4],gt000[1,4]),p=0.5)[3])

#20 +yse and sd++
exp4_pvalues_011_100_21=cbind(binom.test(c(gt011[1,1],gt100[1,1]),p=0.5)[3],
                              binom.test(c(gt011[1,2],gt100[1,2]),p=0.5)[3],
                              binom.test(c(gt011[1,3],gt100[1,3]),p=0.5)[3],
                              binom.test(c(gt011[1,4],gt100[1,4]),p=0.5)[3])

#20 sdy+ and ++se
exp4_pvalues_110_001_21=cbind(binom.test(c(gt110[1,1],gt001[1,1]),p=0.5)[3],
                              binom.test(c(gt110[1,2],gt001[1,2]),p=0.5)[3],
                              binom.test(c(gt110[1,3],gt001[1,3]),p=0.5)[3],
                              binom.test(c(gt110[1,4],gt001[1,4]),p=0.5)[3])

#20 +y+ and sd+se
exp4_pvalues_010_101_21=cbind(binom.test(c(gt010[1,1],gt101[1,1]),p=0.5)[3],
                              binom.test(c(gt010[1,2],gt101[1,2]),p=0.5)[3],
                              binom.test(c(gt010[1,3],gt101[1,3]),p=0.5)[3],
                              binom.test(c(gt010[1,4],gt101[1,4]),p=0.5)[3])


#26 sdyse and +++
exp4_pvalues_111_000_26=cbind(binom.test(c(gt111[2,1],gt000[2,1]),p=0.5)[3],
                              binom.test(c(gt111[2,2],gt000[2,2]),p=0.5)[3],
                              binom.test(c(gt111[2,3],gt000[2,3]),p=0.5)[3],
                              binom.test(c(gt111[2,4],gt000[2,4]),p=0.5)[3])


#26 +yse and sd++
exp4_pvalues_011_100_26=cbind(binom.test(c(gt011[2,1],gt100[2,1]),p=0.5)[3],
                              binom.test(c(gt011[2,2],gt100[2,2]),p=0.5)[3],
                              binom.test(c(gt011[2,3],gt100[2,3]),p=0.5)[3],
                              binom.test(c(gt011[2,4],gt100[2,4]),p=0.5)[3])


#26 sdy+ and ++se
exp4_pvalues_110_001_26=cbind(binom.test(c(gt110[2,1],gt001[2,1]),p=0.5)[3],
                              binom.test(c(gt110[2,2],gt001[2,2]),p=0.5)[3],
                              binom.test(c(gt110[2,3],gt001[2,3]),p=0.5)[3],
                              binom.test(c(gt110[2,4],gt001[2,4]),p=0.5)[3])


#26 +y+ and sd+se
exp4_pvalues_010_101_26=cbind(binom.test(c(gt010[2,1],gt101[2,1]),p=0.5)[3],
                              binom.test(c(gt010[2,2],gt101[2,2]),p=0.5)[3],
                              binom.test(c(gt010[2,3],gt101[2,3]),p=0.5)[3],
                              binom.test(c(gt010[2,4],gt101[2,4]),p=0.5)[3])
exp4_pvalues_haplotype_types_bydayandtreatment=rbind(exp4_pvalues_111_000_21,exp4_pvalues_011_100_21,exp4_pvalues_110_001_21,exp4_pvalues_010_101_21,exp4_pvalues_111_000_26,exp4_pvalues_011_100_26,exp4_pvalues_110_001_26,exp4_pvalues_010_101_26)
colnames(exp4_pvalues_haplotype_types_bydayandtreatment)=c("1-3","4-6","7-9","10-12")
rownames(exp4_pvalues_haplotype_types_bydayandtreatment)=c("sdyse_+++_21","+yse_sd++_21","sdy+_++se_21","+y+_sd+se_21","sdyse_+++_26","+yse_sd++_26","sdy+_++se_26","+y+_sd+se_26")
write.csv(exp4_pvalues_haplotype_types_bydayandtreatment,"exp4_pvalues_haplotype_types_bydayandtreatment.csv")


#use data to get fecundity calculation
dataset2$fecundity=dataset2$total_offspring/dataset2$Num_moms
dataset2$Treatment=as.factor(dataset2$Treatment)

#divide total_offspring females and males
dataset2$fecundity_female=dataset2$Numbfemales.sum/dataset2$Num_moms
dataset2$fecundity_male=dataset2$male.sum/dataset2$Num_moms

#We will write our data into a file so we can read it in for later analysis. 
write.csv(dataset2,file="Experiment4_cleanedup.csv")

for (n in 1:length(dataset2$F1.Vial)) {
dataset2$NCO_skew[n]=min(dataset2$gt1a_111.sum[n],dataset2$gt1b_000.sum[n])/max(dataset2$gt1a_111.sum[n],dataset2$gt1b_000.sum[n])
}
summary(dataset2$NCO_skew)


for (n in 1:length(dataset2$F1.Vial)) {
  dataset2$SCO1_skew[n]=min(dataset2$gt2a_011.sum[n],dataset2$gt2b_100.sum[n])/max(dataset2$gt2a_011.sum[n],dataset2$gt2b_100.sum[n])
}


for (n in 1:length(dataset2$F1.Vial)) {
  dataset2$SCO2_skew[n]=min(dataset2$gt3a_110.sum[n],dataset2$gt3b_001.sum[n])/max(dataset2$gt3a_110.sum[n],dataset2$gt3b_001.sum[n])
}

for (n in 1:length(dataset2$F1.Vial)) {
  dataset2$DCO_skew[n]=min(dataset2$gt4a_010.sum[n],dataset2$gt4b_101.sum[n])/max(dataset2$gt4a_010.sum[n],dataset2$gt4b_101.sum[n])
}


skew=dataset2[, c(3,21,22,23,24)]

skew2=melt(skew, id="Treatment")


pdf("skew_exp4.pdf")
skew_figure_exp4=ggplot(skew2, aes(x = Treatment, y = value, color = variable)) +  # ggplot function
  geom_boxplot()+theme_base()+ylab("haplotype bias")+xlab("Treatment")

skew_figure_exp4
dev.off()
#Fecundity_figure
#dataset2=read.csv("Experiment4_cleanedup.csv")

#poisson regression, similar to a t-test for count data
fit=glm(fecundity_male~Treatment*Day,data=dataset2,family=quasipoisson)
summary(fit)
anova_fec=anova(fit, test="Chisq")
anova_fec
write.csv(anova_fec,"Exp4_fecundity_model_table.csv")


fit_contrast <- emmeans::emmeans(fit, "Treatment", by="Day", mode="kenward-roger")
fit_contr <- contrast(fit_contrast, method="trt.vs.ctrl")

pheno_contr <- as.data.frame(summary(fit_contr))
pheno_contr
write.csv(pheno_contr,"Exp4_fecundity_male_posthoc_table.csv")

#convert p-values to stars for plot
sig=ifelse(pheno_contr$p.value<0.001,"***",ifelse(pheno_contr$p.value<0.01,"**",ifelse(pheno_contr$p.value<0.05,"*","")))

#Fecundity figure for the paper
pdf("Exp4_Fecundity.pdf")
#Fecund_figure_exp4=ggplot(aes(y=fecundity,x=Day, col=Treatment,label=Num_moms),data=dataset2)+ylab("# Progeny per mom")+theme_base()

#Fecund_figure_exp4=Fecund_figure_exp4+scale_colour_manual(values=c("blue", "red"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))
#Fecund_figure_exp4=Fecund_figure_exp4+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)+
#  annotate(geom="text", x=1, y=135, label=sig[1],size=10)+annotate(geom="text", x=2, y=135, label=sig[2],size=10)+annotate(geom="text", x=3, y=135, label=sig[3],size=10)+annotate(geom="text", x=4, y=135, label=sig[4],size=10)
#Fecund_figure_exp4


#Boxplots for fecundity analysis
Fecund_figure_exp4=ggplot(dataset2, aes(x=Day, y=fecundity, col=Treatment)) + theme_base()+ylab("# Progeny per mom")+ggtitle("Experiment 4")+theme(axis.text.x = element_text(angle = 45))+
  geom_boxplot()+scale_x_discrete(name="Days", labels=c("1-3","4-6","7-9","10-12"))+
  annotate(geom="text", x=1, y=135, label=sig[1],size=10)+annotate(geom="text", x=2, y=135, label=sig[2],size=10)+annotate(geom="text", x=3, y=135, label=sig[3],size=10)+annotate(geom="text", x=4, y=135, label=sig[4],size=10)+
scale_color_manual(values = c("blue","red"))
Fecund_figure_exp4

fecundity_female_male=dataset2[,c(2,3,21,22)]
fecundity_female_male$treatment_day=paste(fecundity_female_male$Treatment,fecundity_female_male$Day,sep = "_")
fecundity_female_male=fecundity_female_male[,c(3,4,5)]
fecundity_female_male2=melt(fecundity_female_male, id="treatment_day")


pdf("fecundity_male_female_exp4.pdf")
Fecund_figure_exp5_maleandfemale=ggplot(fecundity_female_male2, aes(x=treatment_day, y=value, col=variable)) + theme_base()+ylab("# Progeny per mom")+ggtitle("Experiment 5")+theme(axis.text.x = element_text(angle = 45))+
  geom_boxplot()+scale_x_discrete(name="Days", labels=c("21_day1-3","21_day4-6","21_day7-9","21_day10-12","26_day1-3","26_day4-6","26_day7-9","26_day10-12"))+ylim(0,50)+
  #annotate(geom="text", x=1, y=38, label=sig[1],size=10)+annotate(geom="text", x=2, y=38, label=sig[2],size=10)+annotate(geom="text", x=3, y=38, label=sig[3],size=10)+annotate(geom="text", x=4, y=38, label=sig[4],size=10)+
  scale_color_manual(values = c("blue","red","green","orange","purple","yellow","pink","black","grey","dark blue"))
Fecund_figure_exp5_maleandfemale
dev.off()


pdf("fecundity male and female_exp4.pdf")

#poisson regression, similar to a t-test for count data
fit=glm(fecundity_male~Treatment*Day,data=dataset2,family=quasipoisson)
summary(fit)
anova_fec=anova(fit, test="Chisq")
anova_fec
write.csv(anova_fec,"Exp4_fecundity_model_table.csv")


fit_contrast <- emmeans::emmeans(fit, "Treatment", by="Day", mode="kenward-roger")
fit_contr <- contrast(fit_contrast, method="trt.vs.ctrl")

pheno_contr <- as.data.frame(summary(fit_contr))
pheno_contr
write.csv(pheno_contr,"Exp4_fecundity_male_posthoc_table.csv")

#convert p-values to stars for plot
sig=ifelse(pheno_contr$p.value<0.001,"***",ifelse(pheno_contr$p.value<0.01,"**",ifelse(pheno_contr$p.value<0.05,"*","")))


Fecund_figure_exp4=ggplot(dataset2, aes(x=Day, y=fecundity_male, col=Treatment)) + theme_base()+ylab("# Progeny per mom")+ggtitle("Experiment 4")+theme(axis.text.x = element_text(angle = 45))+
  geom_boxplot()+scale_x_discrete(name="Days", labels=c("1-3","4-6","7-9","10-12"))+
  annotate(geom="text", x=1, y=135, label=sig[1],size=10)+annotate(geom="text", x=2, y=135, label=sig[2],size=10)+annotate(geom="text", x=3, y=135, label=sig[3],size=10)+annotate(geom="text", x=4, y=135, label=sig[4],size=10)+
  scale_color_manual(values = c("blue","red"))
Fecund_figure_exp4


#poisson regression, similar to a t-test for count data
fit=glm(fecundity_female~Treatment*Day,data=dataset2,family=quasipoisson)
summary(fit)
anova_fec=anova(fit, test="Chisq")
anova_fec
write.csv(anova_fec,"Exp4_fecundity_model_table.csv")


fit_contrast <- emmeans::emmeans(fit, "Treatment", by="Day", mode="kenward-roger")
fit_contr <- contrast(fit_contrast, method="trt.vs.ctrl")

pheno_contr <- as.data.frame(summary(fit_contr))
pheno_contr
write.csv(pheno_contr,"Exp4_fecundity_male_posthoc_table.csv")

#convert p-values to stars for plot
sig=ifelse(pheno_contr$p.value<0.001,"***",ifelse(pheno_contr$p.value<0.01,"**",ifelse(pheno_contr$p.value<0.05,"*","")))

Fecund_figure_exp4=ggplot(dataset2, aes(x=Day, y=fecundity_female, col=Treatment)) + theme_base()+ylab("# Progeny per mom")+ggtitle("Experiment 4")+theme(axis.text.x = element_text(angle = 45))+
  geom_boxplot()+scale_x_discrete(name="Days", labels=c("1-3","4-6","7-9","10-12"))+
  annotate(geom="text", x=1, y=135, label=sig[1],size=10)+annotate(geom="text", x=2, y=135, label=sig[2],size=10)+annotate(geom="text", x=3, y=135, label=sig[3],size=10)+annotate(geom="text", x=4, y=135, label=sig[4],size=10)+
  scale_color_manual(values = c("blue","red"))
Fecund_figure_exp4

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


write.csv(dataset2,"experiment4_combined dataset.csv")

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
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.75,1.75)+ggtitle("Experiment 4")+
  scale_x_discrete(name="Days post-mating",labels=c("1-3","4-6","7-9","10-12"))+geom_line(size=2,show.legend = FALSE)+geom_errorbar(aes(ymin=OR-SE,ymax=OR+SE))+theme(axis.text.x = element_text(angle = 45))+
  annotate(geom="text", x=1, y=1.45, label=ysd_sig[1],color="#f1a340",size=10)+annotate(geom="text", x=2, y=1.45, label=ysd_sig[2],color="#f1a340",size=10)+annotate(geom="text", x=3, y=1.45, label=ysd_sig[3],color="#f1a340",size=10)+annotate(geom="text", x=4, y=1.45, label=ysd_sig[4],color="#f1a340",size=10)+
  annotate(geom="text", x=1, y=1.4, label=yse_sig[1],color="#998ec3",size=10)+annotate(geom="text", x=2, y=1.4, label=yse_sig[2],color="#998ec3",size=10)+annotate(geom="text", x=3, y=1.4, label=yse_sig[3],color="#998ec3",size=10)+annotate(geom="text", x=4, y=1.4, label=yse_sig[4],color="#998ec3",size=10)
odds_figure_exp4
dev.off()

dataset2$Exp_DCO=(((dataset2$SCO_1.sum+dataset2$DCO.sum)/dataset2$male.sum)*((dataset2$SCO_2.sum+dataset2$DCO.sum)/dataset2$male.sum))
dataset2$Obs_DCO=dataset2$DCO.sum/dataset2$male.sum

dataset2$COC=dataset2$Obs_DCO/dataset2$Exp_DCO
dataset2$Interference=1-dataset2$COC

#model
fit5=glm(Interference~Treatment*Day,data=dataset2)
summary(fit5)
anova_coi=anova(fit5, test="Chisq")
anova_coi
write.csv(anova_coi,"exp4_interference_model_table.csv")


fit_contrast5 <- emmeans::emmeans(fit5, "Treatment", by="Day", mode="kenward-roger")
fit_contr5 <- contrast(fit_contrast5, method="trt.vs.ctrl")

pheno_contr5 <- as.data.frame(summary(fit_contr5))
pheno_contr5
write.csv(pheno_contr5,"exp4_interference_posthoc_table.csv")

or=exp(pheno_contr5$estimate)
or

#convert p-values to stars for plot
sig2=ifelse(pheno_contr5$p.value<0.001,"***",ifelse(pheno_contr5$p.value<0.01,"**",ifelse(pheno_contr5$p.value<0.05,"*","")))

#Interference figure for the paper
#postscript ("../Figures/Figure.4B.eps", width=4, height=3, horizontal=FALSE, pointsize=5)
#png("../Figures/Figure4B.png")
COI_figure=ggplot(aes(y=Interference,x=Day, col=Treatment,label=Num_moms),data=dataset2)+ylab("COI")+ggtitle("Interference")+theme_base()+ylim(-0.25,0.4)
COI_figure=COI_figure+scale_colour_manual(values=c("#7b3294","#008837"))+geom_point(size=1.5)+scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))
COI_figure=COI_figure+stat_summary(fun = median, geom="line",aes(group=Treatment),size=1)+
  annotate(geom="text", x=1, y=0.25, label=sig2[1],size=5)+annotate(geom="text", x=2, y=0.25, label=sig2[2],size=5)+annotate(geom="text", x=3, y=0.25, label=sig2[3],size=5)+annotate(geom="text", x=4, y=0.25, label=sig2[4],size=5)
COI_figure


# Experiment 3 (Exp5)-------------------------------------------------------------

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
#Define Crossovers
exp5$gamete_class=ifelse(exp5$sd==exp5$y & exp5$y==exp5$se & exp5$y==1,"gt1a",
                         ifelse(exp5$sd==exp5$y & exp5$y==exp5$se & exp5$y==0,"gt1b",
                                ifelse(exp5$sd!=exp5$y & exp5$y==exp5$se & exp5$y==1,"gt2a",
                                       ifelse(exp5$sd!=exp5$y & exp5$y==exp5$se & exp5$y==0,"gt2b",
                                              ifelse(exp5$sd==exp5$y & exp5$y!=exp5$se & exp5$y==1,"gt3a",
                                                     ifelse(exp5$sd==exp5$y & exp5$y!=exp5$se & exp5$y==0,"gt3b",
                                                            ifelse(exp5$sd!=exp5$y & exp5$y!=exp5$se & exp5$y==1,"gt4a",
                                                                   ifelse(exp5$sd!=exp5$y & exp5$y!=exp5$se & exp5$y==0,"gt4b",
                                                                          "error"))))))))
#Sanity check, there should not be any error!
exp5[exp5$gamete_class=="error",]

#add columns for counting
exp5$gt1a_111=ifelse(exp5$gamete_class=="gt1a",exp5$numbMales,0)
exp5$gt1b_000=ifelse(exp5$gamete_class=="gt1b",exp5$numbMales,0)
exp5$gt2a_011=ifelse(exp5$gamete_class=="gt2a",exp5$numbMales,0)
exp5$gt2b_100=ifelse(exp5$gamete_class=="gt2b",exp5$numbMales,0)
exp5$gt3a_110=ifelse(exp5$gamete_class=="gt3a",exp5$numbMales,0)
exp5$gt3b_001=ifelse(exp5$gamete_class=="gt3b",exp5$numbMales,0)
exp5$gt4a_010=ifelse(exp5$gamete_class=="gt4a",exp5$numbMales,0)
exp5$gt4b_101=ifelse(exp5$gamete_class=="gt4b",exp5$numbMales,0)

gt1a_111=sum(exp5$numbMales*exp5$gt1a_111)
gt1b_000=sum(exp5$numbMales*exp5$gt1b_000) #p-value < 2.2e-16
gt2a_011=sum(exp5$numbMales*exp5$gt2a_011)
gt2b_100=sum(exp5$numbMales*exp5$gt2b_100)
gt3a_110=sum(exp5$numbMales*exp5$gt3a_110)
gt3b_001=sum(exp5$numbMales*exp5$gt3b_001) #p-value = 1.306e-15
gt4a_010=sum(exp5$numbMales*exp5$gt4a_010)
gt4b_101=sum(exp5$numbMales*exp5$gt4b_101)

cbind(gt1a_111,gt1b_000,gt2a_011,gt2b_100,gt3a_110,gt3b_001,gt4a_010,gt4b_101)
binom.test(gt1a_111,gt1b_000,p=0.5)
binom.test(gt2a_011,gt2b_100,p=0.5)
binom.test(gt3a_110,gt3b_001,p=0.5)
binom.test(gt4a_010,gt4b_101,p=0.5)

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
dataset=summaryBy(male+gt1a_111+gt1b_000+gt2a_011+gt2b_100+gt3a_110+gt3b_001+gt4a_010+gt4b_101+num_co+SCO_1+SCO_2+DCO~F1.Vial+Day+Treatment,data=exp5_merged, FUN=sum,na.rm=T)

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
dataset2$fecundity_male=dataset2$Numbfemales.sum/dataset2$Num_moms
dataset2$fecundity_female=dataset2$male.sum/dataset2$Num_moms
#to calculate the skew based on the inviability of mutations
dataset2$skew=dataset$

#We will write our data into a file so we can read it in for later analysis. 
write.csv(dataset2,file="Experiment5_cleanedup.csv")


for (n in 1:length(dataset2$F1.Vial)) {
  dataset2$NCO_skew[n]=min(dataset2$gt1a_111.sum[n],dataset2$gt1b_000.sum[n])/max(dataset2$gt1a_111.sum[n],dataset2$gt1b_000.sum[n])
}
summary(dataset2$NCO_skew)


for (n in 1:length(dataset2$F1.Vial)) {
  dataset2$SCO1_skew[n]=min(dataset2$gt2a_011.sum[n],dataset2$gt2b_100.sum[n])/max(dataset2$gt2a_011.sum[n],dataset2$gt2b_100.sum[n])
}


for (n in 1:length(dataset2$F1.Vial)) {
  dataset2$SCO2_skew[n]=min(dataset2$gt3a_110.sum[n],dataset2$gt3b_001.sum[n])/max(dataset2$gt3a_110.sum[n],dataset2$gt3b_001.sum[n])
}

for (n in 1:length(dataset2$F1.Vial)) {
  dataset2$DCO_skew[n]=min(dataset2$gt4a_010.sum[n],dataset2$gt4b_101.sum[n])/max(dataset2$gt4a_010.sum[n],dataset2$gt4b_101.sum[n])
}


skew=dataset2[, c(3,21,22,23,24)]

skew2=melt(skew, id="Treatment")


pdf("skew_exp5.pdf")
skew_figure_exp5=ggplot(skew2, aes(x = Treatment, y = value, color = variable)) +  # ggplot function
  geom_boxplot()+theme_base()+ylab("haplotype bias")+xlab("Treatment")

skew_figure_exp5
dev.off()

#gamete analysis by treatment

gt111=tapply(dataset2$gt1a_111.sum,dataset2$Treatment,sum,na.rm=T)

gt000=tapply(dataset2$gt1b_000.sum,dataset2$Treatment,sum,na.rm=T)
gt011=tapply(dataset2$gt2a_011.sum,dataset2$Treatment,sum,na.rm=T)
gt100=tapply(dataset2$gt2b_100.sum,dataset2$Treatment,sum,na.rm=T)
gt110=tapply(dataset2$gt3a_110.sum,dataset2$Treatment,sum,na.rm=T)
gt001=tapply(dataset2$gt3b_001.sum,dataset2$Treatment,sum,na.rm=T)
gt010=tapply(dataset2$gt4a_010.sum,dataset2$Treatment,sum,na.rm=T)
gt101=tapply(dataset2$gt4b_101.sum,dataset2$Treatment,sum,na.rm=T)
i=rbind(gt111, gt000, gt011, gt100, gt110, gt001, gt010, gt101)
write.csv(i,"gamete_class_by_treatment.csv")
#gamete analysis by day
gt111=tapply(dataset2$gt1a_111.sum,dataset2$Day,sum,na.rm=T)
gt000=tapply(dataset2$gt1b_000.sum,dataset2$Day,sum,na.rm=T)
gt011=tapply(dataset2$gt2a_011.sum,dataset2$Day,sum,na.rm=T)
gt100=tapply(dataset2$gt2b_100.sum,dataset2$Day,sum,na.rm=T)
gt110=tapply(dataset2$gt3a_110.sum,dataset2$Day,sum,na.rm=T)
gt001=tapply(dataset2$gt3b_001.sum,dataset2$Day,sum,na.rm=T)
gt010=tapply(dataset2$gt4a_010.sum,dataset2$Day,sum,na.rm=T)
gt101=tapply(dataset2$gt4b_101.sum,dataset2$Day,sum,na.rm=T)
s=rbind(gt111, gt000, gt011, gt100, gt110, gt001, gt010, gt101)
write.csv(s,"gamete_class_byday.csv")
gt111=tapply(dataset2$gt1a_111.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt000=tapply(dataset2$gt1b_000.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt011=tapply(dataset2$gt2a_011.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt100=tapply(dataset2$gt2b_100.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt110=tapply(dataset2$gt3a_110.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt001=tapply(dataset2$gt3b_001.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt010=tapply(dataset2$gt4a_010.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt101=tapply(dataset2$gt4b_101.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
t=rbind(gt111, gt000, gt011, gt100, gt110, gt001, gt010, gt101)
write.csv(t,"gamete_class_byday_and_treatment.csv")
#21 sdyse and +++
exp5_pvalues_111_000_21=cbind(binom.test(c(gt111[1,1],gt000[1,1]),p=0.5)[3],
                              binom.test(c(gt111[1,2],gt000[1,2]),p=0.5)[3],
                              binom.test(c(gt111[1,3],gt000[1,3]),p=0.5)[3],
                              binom.test(c(gt111[1,4],gt000[1,4]),p=0.5)[3],
                              binom.test(c(gt111[1,5],gt000[1,5]),p=0.5)[3])

#20 +yse and sd++
exp5_pvalues_011_100_21=cbind(binom.test(c(gt011[1,1],gt100[1,1]),p=0.5)[3],
                              binom.test(c(gt011[1,2],gt100[1,2]),p=0.5)[3],
                              binom.test(c(gt011[1,3],gt100[1,3]),p=0.5)[3],
                              binom.test(c(gt011[1,4],gt100[1,4]),p=0.5)[3],
                              binom.test(c(gt011[1,5],gt100[1,5]),p=0.5)[3])

#20 sdy+ and ++se
exp5_pvalues_110_001_21=cbind(binom.test(c(gt110[1,1],gt001[1,1]),p=0.5)[3],
                              binom.test(c(gt110[1,2],gt001[1,2]),p=0.5)[3],
                              binom.test(c(gt110[1,3],gt001[1,3]),p=0.5)[3],
                              binom.test(c(gt110[1,4],gt001[1,4]),p=0.5)[3],
                              binom.test(c(gt110[1,5],gt001[1,5]),p=0.5)[3])

#20 +y+ and sd+se
exp5_pvalues_010_101_21=cbind(binom.test(c(gt010[1,1],gt101[1,1]),p=0.5)[3],
                              binom.test(c(gt010[1,2],gt101[1,2]),p=0.5)[3],
                              binom.test(c(gt010[1,3],gt101[1,3]),p=0.5)[3],
                              binom.test(c(gt010[1,4],gt101[1,4]),p=0.5)[3],
                              binom.test(c(gt010[1,5],gt101[1,5]),p=0.5)[3])


#26 sdyse and +++
exp5_pvalues_111_000_26=cbind(binom.test(c(gt111[2,1],gt000[2,1]),p=0.5)[3],
                              binom.test(c(gt111[2,2],gt000[2,2]),p=0.5)[3],
                              binom.test(c(gt111[2,3],gt000[2,3]),p=0.5)[3],
                              binom.test(c(gt111[2,4],gt000[2,4]),p=0.5)[3],
                              binom.test(c(gt111[2,5],gt000[2,5]),p=0.5)[3])


#26 +yse and sd++
exp5_pvalues_011_100_26=cbind(binom.test(c(gt011[2,1],gt100[2,1]),p=0.5)[3],
                              binom.test(c(gt011[2,2],gt100[2,2]),p=0.5)[3],
                              binom.test(c(gt011[2,3],gt100[2,3]),p=0.5)[3],
                              binom.test(c(gt011[2,4],gt100[2,4]),p=0.5)[3],
                              binom.test(c(gt011[2,5],gt100[2,5]),p=0.5)[3])


#26 sdy+ and ++se
exp5_pvalues_110_001_26=cbind(binom.test(c(gt110[2,1],gt001[2,1]),p=0.5)[3],
                              binom.test(c(gt110[2,2],gt001[2,2]),p=0.5)[3],
                              binom.test(c(gt110[2,3],gt001[2,3]),p=0.5)[3],
                              binom.test(c(gt110[2,4],gt001[2,4]),p=0.5)[3],
                              binom.test(c(gt110[2,5],gt001[2,5]),p=0.5)[3])


#26 +y+ and sd+se
exp5_pvalues_010_101_26=cbind(binom.test(c(gt010[2,1],gt101[2,1]),p=0.5)[3],
                              binom.test(c(gt010[2,2],gt101[2,2]),p=0.5)[3],
                              binom.test(c(gt010[2,3],gt101[2,3]),p=0.5)[3],
                              binom.test(c(gt010[2,4],gt101[2,4]),p=0.5)[3],
                              binom.test(c(gt010[2,5],gt101[2,5]),p=0.5)[3])
exp5_pvalues_haplotype_types_bydayandtreatment=rbind(exp5_pvalues_111_000_21,exp5_pvalues_011_100_21,exp5_pvalues_110_001_21,exp5_pvalues_010_101_21,exp5_pvalues_111_000_26,exp5_pvalues_011_100_26,exp5_pvalues_110_001_26,exp5_pvalues_010_101_26)
colnames(exp5_pvalues_haplotype_types_bydayandtreatment)=c("6","7","8","9","10")
rownames(exp5_pvalues_haplotype_types_bydayandtreatment)=c("sdyse_+++_21","+yse_sd++_21","sdy+_++se_21","+y+_sd+se_21","sdyse_+++_26","+yse_sd++_26","sdy+_++se_26","+y+_sd+se_26")
write.csv(exp5_pvalues_haplotype_types_bydayandtreatment,"exp5_pvalues_haplotype_types_bydayandtreatment.csv")


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
#Fecund_figure_exp5=ggplot(aes(y=fecundity,x=Day, col=Treatment,label=Num_moms),data=dataset2)+ylab("# Progeny per mom")+theme_base()

#Fecund_figure_exp5=Fecund_figure_exp5+scale_colour_manual(values=c("blue", "red"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("6","7","8","9","10"))+ylim(0,50)
#Fecund_figure_exp5=Fecund_figure_exp5+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)+
 # annotate(geom="text", x=1, y=40, label=sig[1],size=10)+annotate(geom="text", x=2, y=40, label=sig[2],size=10)+annotate(geom="text", x=3, y=40, label=sig[3],size=10)+annotate(geom="text", x=4, y=40, label=sig[4],size=10)
#Fecund_figure_exp5


#boxplots for fecundity analysis
Fecund_figure_exp5=ggplot(dataset2, aes(x=Day, y=fecundity, col=Treatment)) + theme_base()+ylab("# Progeny per mom")+ggtitle("Experiment 5")+theme(axis.text.x = element_text(angle = 45))+
  geom_boxplot()+scale_x_discrete(name="Days", labels=c("6","7","8","9","10"))+ylim(0,50)+
  annotate(geom="text", x=1, y=38, label=sig[1],size=10)+annotate(geom="text", x=2, y=38, label=sig[2],size=10)+annotate(geom="text", x=3, y=38, label=sig[3],size=10)+annotate(geom="text", x=4, y=38, label=sig[4],size=10)+
scale_color_manual(values = c("blue","red"))
Fecund_figure_exp5

fecundity_female_male=dataset2[,c(2,3,21,22)]
fecundity_female_male$treatment_day=paste(fecundity_female_male7$Treatment,fecundity_female_male$Day,sep = "_")
fecundity_female_male=fecundity_female_male[,c(3,4,5)]
fecundity_female_male2=melt(fecundity_female_male, id="treatment_day")

pdf("fecundity_male_female.pdf")
Fecund_figure_exp5_maleandfemale=ggplot(fecundity_female_male2, aes(x=treatment_day, y=value, col=variable)) + theme_base()+ylab("# Progeny per mom")+ggtitle("Experiment 5")+theme(axis.text.x = element_text(angle = 45))+
  geom_boxplot()+scale_x_discrete(name="Days", labels=c("21_day6","_21_day7","_21_day8","21_day9","21_day10","26_day6","26_day7","26_day8","26_day9","26_day10"))+ylim(0,50)+
  annotate(geom="text", x=1, y=38, label=sig[1],size=10)+annotate(geom="text", x=2, y=38, label=sig[2],size=10)+annotate(geom="text", x=3, y=38, label=sig[3],size=10)+annotate(geom="text", x=4, y=38, label=sig[4],size=10)+
  scale_color_manual(values = c("blue","red","green","orange","purple","yellow","pink","black","grey","dark blue"))
Fecund_figure_exp5_maleandfemale
dev.off()


pdf("fecundity for males and females_exp5.pdf")
#fecundity graph for males

fit=glm(fecundity_male~Treatment*Day,data=dataset2,family=quasipoisson)
#summary(fit)
anova_fec=anova(fit, test="Chisq")
write.csv(anova_fec,"Exp5_fecundity_model_table.csv")

fit_contrast <- emmeans::emmeans(fit, "Treatment", by="Day", mode="kenward-roger")
fit_contr <- contrast(fit_contrast, method="trt.vs.ctrl")

pheno_contr <- as.data.frame(summary(fit_contr))
pheno_contr
write.csv(pheno_contr,"Exp5_fecundity_male_posthoc_table.csv")

#convert p-values to stars for plot
sig=ifelse(pheno_contr$p.value<0.001,"***",ifelse(pheno_contr$p.value<0.01,"**",ifelse(pheno_contr$p.value<0.05,"*","")))

Fecund_figure_exp5_male=ggplot(dataset2, aes(x=Day, y=fecundity_male, col=Treatment)) + theme_base()+ylab("# Progeny per mom")+ggtitle("Experiment 5")+theme(axis.text.x = element_text(angle = 45))+
  geom_boxplot()+scale_x_discrete(name="Days", labels=c("6","7","8","9","10"))+ylim(0,50)+
  annotate(geom="text", x=1, y=38, label=sig[1],size=10)+annotate(geom="text", x=2, y=38, label=sig[2],size=10)+annotate(geom="text", x=3, y=38, label=sig[3],size=10)+annotate(geom="text", x=4, y=38, label=sig[4],size=10)+
  scale_color_manual(values = c("blue","red"))
Fecund_figure_exp5_male

#fecundity graph for females

fit=glm(fecundity_female~Treatment*Day,data=dataset2,family=quasipoisson)
#summary(fit)
anova_fec=anova(fit, test="Chisq")
write.csv(anova_fec,"Exp5_fecundity_model_table.csv")

fit_contrast <- emmeans::emmeans(fit, "Treatment", by="Day", mode="kenward-roger")
fit_contr <- contrast(fit_contrast, method="trt.vs.ctrl")

pheno_contr <- as.data.frame(summary(fit_contr))
pheno_contr

write.csv(pheno_contr,"Exp5_fecundity_female_posthoc_table.csv")

#convert p-values to stars for plot
sig=ifelse(pheno_contr$p.value<0.001,"***",ifelse(pheno_contr$p.value<0.01,"**",ifelse(pheno_contr$p.value<0.05,"*","")))

Fecund_figure_exp5_male=ggplot(dataset2, aes(x=Day, y=fecundity_female, col=Treatment)) + theme_base()+ylab("# Progeny per mom")+ggtitle("Experiment 5")+theme(axis.text.x = element_text(angle = 45))+
  geom_boxplot()+scale_x_discrete(name="Days", labels=c("6","7","8","9","10"))+ylim(0,50)+
  annotate(geom="text", x=1, y=38, label=sig[1],size=10)+annotate(geom="text", x=2, y=38, label=sig[2],size=10)+annotate(geom="text", x=3, y=38, label=sig[3],size=10)+annotate(geom="text", x=4, y=38, label=sig[4],size=10)+
  scale_color_manual(values = c("blue","red"))
Fecund_figure_exp5_male

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

#Kosambi correction
dataset2$kosambi_rec_rate_ysd=((2.71^(4*dataset2$rec_rate_ysd)-1)/2*(2.71^(-4*dataset2$rec_rate_ysd)+1))/10
dataset2$kosambi_rec_rate_yse=((2.71^(4*dataset2$rec_rate_yse)-1)/2*(2.71^(-4*dataset2$rec_rate_yse)+1))/10



write.csv(dataset2,"experiment5_combined dataset.csv")

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
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.75,1.76)+
  geom_line(size=2)+geom_errorbar(aes(ymin=OR-SE,ymax=OR+SE))+ggtitle("Experiment 5")+theme(axis.text.x = element_text(angle = 45))+
  annotate(geom="text", x=6, y=1.65, label=ysd_sig[1],color="#f1a340",size=10)+annotate(geom="text", x=7, y=1.65, label=ysd_sig[2],color="#f1a340",size=10)+annotate(geom="text", x=8, y=1.65, label=ysd_sig[3],color="#f1a340",size=10)+annotate(geom="text", x=9, y=1.65, label=ysd_sig[4],color="#f1a340",size=10)+annotate(geom="text", x=9, y=1.65, label=yse_sig[4],color="#f1a340",size=10)+
  annotate(geom="text", x=6, y=1.6, label=yse_sig[1],color="#998ec3",size=10)+annotate(geom="text", x=7, y=1.6, label=yse_sig[2],color="#998ec3",size=10)+annotate(geom="text", x=8, y=1.6, label=yse_sig[3],color="#998ec3",size=10)+annotate(geom="text", x=9, y=1.6, label=yse_sig[4],color="#998ec3",size=10)+annotate(geom="text", x=9, y=1.6, label=yse_sig[4],color="#998ec3",size=10)
odds_figure_exp5
dev.off()


#boxplot for total recombination rate
pdf("Exp5_sdyboxplotonly.pdf")

Recomb_figure_ysd=ggplot(aes(y=kosambi_rec_rate_ysd,x=Day, col=Treatment,label=male.sum),data=dataset2)+scale_colour_manual(values=c("blue","red"))+geom_boxplot()+ylab("% recombination")+theme_base()+scale_x_discrete(name="Days post-mating",labels=c("6","7","8","9","10"))
Recomb_figure_ysd=Recomb_figure_ysd+ggtitle("sd-y Interval")+ylim(0,0.6)
Recomb_figure_ysd=Recomb_figure_ysd+annotate(geom="text", x=1, y=0.58, label=ysd_sig[1],size=10)+annotate(geom="text", x=2, y=0.58, label=ysd_sig[2],size=10)+annotate(geom="text", x=3, y=0.58, label=ysd_sig[3],size=10)+annotate(geom="text", x=4, y=0.58, label=ysd_sig[4],size=10)
Recomb_figure_ysd

dev.off()


pdf("Exp5_yseboxplotonly.pdf")


Recomb_figure_yse=ggplot(aes(y=kosambi_rec_rate_yse,x=Day, col=Treatment,label=male.sum),data=dataset2)+scale_colour_manual(values=c("blue","red"))+geom_boxplot()+ylab("% recombination")+theme_base()+scale_x_discrete(name="Days post-mating",labels=c("6","7","8","9","10"))
Recomb_figure_yse=Recomb_figure_yse+ggtitle("y-se Interval")+ylim(0.1,0.6)
Recomb_figure_yse=Recomb_figure_yse+annotate(geom="text", x=1, y=0.58, label=yse_sig[1],size=10)+annotate(geom="text", x=2, y=0.58, label=yse_sig[2],size=10)+annotate(geom="text", x=3, y=0.58, label=yse_sig[3],size=10)+annotate(geom="text", x=4, y=0.58, label=yse_sig[4],size=10)
Recomb_figure_yse

dev.off()

dataset2$Exp_DCO=(((dataset2$SCO_1.sum+dataset2$DCO.sum)/dataset2$male.sum)*((dataset2$SCO_2.sum+dataset2$DCO.sum)/dataset2$male.sum))
dataset2$Obs_DCO=dataset2$DCO.sum/dataset2$male.sum

dataset2$COC=dataset2$Obs_DCO/dataset2$Exp_DCO
dataset2$Interference=1-dataset2$COC

#model
fit5=glm(Interference~Treatment*Day,data=dataset2)
summary(fit5)
anova_coi=anova(fit5, test="Chisq")
anova_coi
write.csv(anova_coi,"exp4_interference_model_table.csv")


fit_contrast5 <- emmeans::emmeans(fit5, "Treatment", by="Day", mode="kenward-roger")
fit_contr5 <- contrast(fit_contrast5, method="trt.vs.ctrl")

pheno_contr5 <- as.data.frame(summary(fit_contr5))
pheno_contr5
write.csv(pheno_contr5,"exp4_interference_posthoc_table.csv")

or=exp(pheno_contr5$estimate)
or

#convert p-values to stars for plot
sig2=ifelse(pheno_contr5$p.value<0.001,"***",ifelse(pheno_contr5$p.value<0.01,"**",ifelse(pheno_contr5$p.value<0.05,"*","")))

#Interference figure for the paper
#postscript ("../Figures/Figure.4B.eps", width=4, height=3, horizontal=FALSE, pointsize=5)
#png("../Figures/Figure4B.png")
COI_figure=ggplot(aes(y=Interference,x=Day, col=Treatment,label=Num_moms),data=dataset2)+ylab("COI")+ggtitle("Interference")+theme_base()+ylim(-1,1)
COI_figure=COI_figure+scale_colour_manual(values=c("#7b3294","#008837"))+geom_point(size=1.5)+scale_x_discrete(name="Day",labels=c("6","7","8","9","10"))
COI_figure=COI_figure+stat_summary(fun = median, geom="line",aes(group=Treatment),size=1)+
  annotate(geom="text", x=1, y=0.25, label=sig2[1],size=5)+annotate(geom="text", x=2, y=0.25, label=sig2[2],size=5)+annotate(geom="text", x=3, y=0.25, label=sig2[3],size=5)+annotate(geom="text", x=4, y=0.25, label=sig2[4],size=5)
COI_figure


# experiment 3 least biased -----------------------------------------------

yv1=subset(exp5,exp5$gamete_class=="gt2a")
yv2=subset(exp5,exp5$gamete_class=="gt4a")
yv3=subset(exp5,exp5$gamete_class=="gt1b")
yv4=subset(exp5,exp5$gamete_class=="gt3b")
yv=rbind(yv1,yv2,yv3,yv4)

#Define Crossovers
yv$co_class=ifelse(yv$sd==yv$y & yv$y==yv$se,"non_CO", 
                   ifelse(yv$sd!=yv$y & yv$y==yv$se,"single_CO_1",
                          ifelse(yv$sd==yv$y & yv$y!=yv$se,"single_CO_2",
                                 ifelse(yv$sd!=yv$y & yv$y!=yv$se,"double_CO",
                                        "error"))))

#add columns for counting
yv$NCO=ifelse(yv$co_class=="non_CO",yv$numbMales,0)
yv$SCO_1=ifelse(yv$co_class=="single_CO_1",yv$numbMales,0)
yv$SCO_2=ifelse(yv$co_class=="single_CO_2",yv$numbMales,0)
yv$DCO=ifelse(yv$co_class=="double_CO",yv$numbMales,0)

#get rough crossover rate as a sanity check the numbers of crossovers at the intervals 
#should equal to the total crossover rate.
nco_count=sum(yv$numbMales*yv$NCO)
sco_count=sum(yv$numbMales*yv$SCO_1, na.rm = TRUE)+sum(yv$numbMales*yv$SCO_2, na.rm = TRUE)
dco_count=sum(yv$numbMales*yv$DCO, na.rm = TRUE)
num_samples=sum(nco_count+sco_count+dco_count, na.rm = TRUE) 


#rate between scalloped and yellow
(sum(yv$numbMales*yv$SCO_1,na.rm = TRUE)+sum(yv$numbMales*yv$DCO,na.rm = TRUE))/num_samples

#rate between sepia and yellow
(sum(yv$numbMales*yv$SCO_2, na.rm = TRUE)+sum(yv$numbMales*yv$DCO,na.rm = TRUE))/num_samples

#we defined the crossover groups as if non crossovers will be identified as 0s; this is same as co_class column, but numerical
#single crossovers 1s
#double crossovers 2s
yv$num_co=ifelse(yv$y==yv$sd & yv$y==yv$se,0, 
                 ifelse(yv$sd==yv$y & yv$y!=yv$se,1*yv$numbMales, 
                        ifelse(yv$sd!=yv$y & yv$y==yv$se,1*yv$numbMales,  
                               ifelse(yv$sd!=yv$y & yv$y!=yv$se,2*yv$numbMales, 
                                      NA))))

#This is for summarizing our data
yv$male=c(yv$numbMales)

#merge with treatment data
yv_merged <- merge(yv, exp5_bc_worksheet, by.x = "ViaNumber", by.y = "Vial.Number", all=T)

#merge female data
yv_female_merged=merge(exp5_female_counts, exp5_bc_worksheet, by.x = "Vial", by.y = "Vial.Number", all=T)

#summarize long form data
dataset=summaryBy(male+num_co+SCO_1+SCO_2+DCO~F1.Vial+Day+Treatment,data=yv_merged, FUN=sum,na.rm=T)

#add in female data
yv_female_short=summaryBy(Numbfemales~F1.Vial+Day+Treatment,data=yv_female_merged, FUN=sum,na.rm=T)
yv_female_short=na.omit(yv_female_short)

#one more merge
dataset2=merge(yv_female_short,dataset, by=c("F1.Vial","Day","Treatment"))

#add in a column for total offspring
dataset2$total_offspring=dataset2$Numbfemales.sum + dataset2$male.sum

#make a vector to store the data
num_moms=vector(mode="numeric",length=length(dataset2$F1.Vial))

#loop through dataset to get the count
for (h in 1:length(dataset2$F1.Vial)) { 
  f1_vial=dataset2$F1.Vial[h]
  mom_ct=length(unique(sort(subset(yv_merged,yv_merged$F1.Vial==f1_vial)$ViaNumber)))
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
write.csv(anova_fec,"yv_fecundity_model_table.csv")

fit_contrast <- emmeans::emmeans(fit, "Treatment", by="Day", mode="kenward-roger")
fit_contr <- contrast(fit_contrast, method="trt.vs.ctrl")

pheno_contr <- as.data.frame(summary(fit_contr))
pheno_contr
write.csv(pheno_contr,"yv_fecundity_posthoc_table.csv")

#convert p-values to stars for plot
sig=ifelse(pheno_contr$p.value<0.001,"***",ifelse(pheno_contr$p.value<0.01,"**",ifelse(pheno_contr$p.value<0.05,"*","")))


#boxplots for fecundity analysis
Fecund_figure_yv=ggplot(dataset2, aes(x=Day, y=fecundity, col=Treatment)) + theme_base()+ylab("# Progeny per mom")+ggtitle("Experiment 5")+theme(axis.text.x = element_text(angle = 45))+
  geom_boxplot()+scale_x_discrete(name="Days", labels=c("6","7","8","9","10"))+ylim(0,50)+
  annotate(geom="text", x=1, y=38, label=sig[1],size=10)+annotate(geom="text", x=2, y=38, label=sig[2],size=10)+annotate(geom="text", x=3, y=38, label=sig[3],size=10)+annotate(geom="text", x=4, y=38, label=sig[4],size=10)+
  scale_color_manual(values = c("blue","red"))
Fecund_figure_yv


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

#Kosambi correction
dataset2$kosambi_rec_rate_ysd=((2.71^(4*dataset2$rec_rate_ysd)-1)/2*(2.71^(-4*dataset2$rec_rate_ysd)+1))/10
dataset2$kosambi_rec_rate_yse=((2.71^(4*dataset2$rec_rate_yse)-1)/2*(2.71^(-4*dataset2$rec_rate_yse)+1))/10



###Recombination figures

#Recomb_figure_yellow_scalloped
Recomb_figure=ggplot(aes(y=rec_rate_ysd,x=Day, col=as.factor(Treatment),label=male.sum),data=dataset2)+scale_colour_manual(values=c("blue", "red"))+geom_point(alpha=0.6,size=3)+ylab("% recombination")+theme_base()+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+scale_x_discrete(name="Day",labels=c("6","7","8","9","10"))
Recomb_figure=Recomb_figure+ggtitle("y-sd Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
Recomb_figure

#Recomb_figure_yellow_sepia
Recomb_figure=ggplot(aes(y=rec_rate_yse,x=Day, col=as.factor(Treatment),label=male.sum),data=dataset2)+scale_colour_manual(values=c("blue","red"))+geom_point(alpha=0.6,size=3)+ylab("% recombination")+theme_base()+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+scale_x_discrete(name="Day",labels=c("6","7","8","9","10"))
Recomb_figure=Recomb_figure+ggtitle("y-se Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
Recomb_figure


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

tapply(dataset2$kosambi_rec_rate_ysd,dataset2$Treatment,mean,na.rm=T)
tapply(dataset2$rec_rate_ysd,dataset2$Day,mean,na.rm=T)

tapply(dataset2$kosambi_rec_rate_yse,dataset2$Treatment,mean,na.rm=T)
tapply(dataset2$rec_rate_yse,dataset2$Day,mean,na.rm=T)

###Odds ratios
#SD-Y REGION

#logistic regression, similar to a t-test for count data
fit3=glmer(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+Treatment*Day,data=dataset2,
           family=binomial(link="logit"))
#coefs=coef(fit3)
anova_sdy=Anova(fit3,test="Chisq")
write.csv(anova_sdy,"yv_recrate_sd-y_model_table.csv")

fit_contrast3 <- emmeans::emmeans(fit3, "Treatment", by="Day", mode="kenward-roger")
fit_contr3 <- contrast(fit_contrast3, method="trt.vs.ctrl")

pheno_contr3 <- as.data.frame(summary(fit_contr3))
pheno_contr3
write.csv(pheno_contr3,"yv_recrate_sd-y_posthoc_table.csv")

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
write.csv(anova_yse,"yv_recrate_y-se_model_table.csv")

fit_contrast4 <- emmeans::emmeans(fit4, "Treatment", by="Day", mode="kenward-roger")
fit_contr4 <- contrast(fit_contrast4, method="trt.vs.ctrl")

pheno_contr4 <- as.data.frame(summary(fit_contr4))
pheno_contr4
write.csv(pheno_contr4,"yv_recrate_y-se_posthoc_table.csv")

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
pdf("yv_odds_ratio.pdf")
odds_figure_yv=ggplot(aes(y=OR,x=Day, col=Interval,group=Interval),data=odds_ratios)+scale_colour_manual(values=c("#f1a340", "#998ec3"))+
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.50,1.76)+
  geom_line(size=2)+geom_errorbar(aes(ymin=OR-SE,ymax=OR+SE))+ggtitle("Experiment 5")+theme(axis.text.x = element_text(angle = 45))+
  annotate(geom="text", x=6, y=1.65, label=ysd_sig[1],color="#f1a340",size=10)+annotate(geom="text", x=7, y=1.65, label=ysd_sig[2],color="#f1a340",size=10)+annotate(geom="text", x=8, y=1.65, label=ysd_sig[3],color="#f1a340",size=10)+annotate(geom="text", x=9, y=1.65, label=ysd_sig[4],color="#f1a340",size=10)+annotate(geom="text", x=9, y=1.65, label=yse_sig[4],color="#f1a340",size=10)+
  annotate(geom="text", x=6, y=1.6, label=yse_sig[1],color="#998ec3",size=10)+annotate(geom="text", x=7, y=1.6, label=yse_sig[2],color="#998ec3",size=10)+annotate(geom="text", x=8, y=1.6, label=yse_sig[3],color="#998ec3",size=10)+annotate(geom="text", x=9, y=1.6, label=yse_sig[4],color="#998ec3",size=10)+annotate(geom="text", x=9, y=1.6, label=yse_sig[4],color="#998ec3",size=10)
odds_figure_yv
dev.off()


#boxplot for total recombination rate
pdf("yv_sdyboxplotonly.pdf")
Recomb_figure=ggplot(aes(y=kosambi_rec_rate_ysd,x=Day, col=Treatment,label=male.sum),data=dataset2)+scale_colour_manual(values=c("blue","red"))+geom_boxplot()+ylab("% recombination")+theme_base()+scale_x_discrete(name="Days post-mating",labels=c("6","7","8","9","10"))
Recomb_figure=Recomb_figure+ggtitle("sd-y Interval")+ylim(0,0.6)
Recomb_figure=Recomb_figure+annotate(geom="text", x=1, y=0.58, label=ysd_sig[1],size=10)+annotate(geom="text", x=2, y=0.58, label=ysd_sig[2],size=10)+annotate(geom="text", x=3, y=0.58, label=ysd_sig[3],size=10)+annotate(geom="text", x=4, y=0.58, label=ysd_sig[4],size=10)
Recomb_figure

dev.off()


pdf("yv_yseboxplotonly.pdf")
Recomb_figure=ggplot(aes(y=kosambi_rec_rate_yse,x=Day, col=Treatment,label=male.sum),data=dataset2)+scale_colour_manual(values=c("blue","red"))+geom_boxplot()+ylab("% recombination")+theme_base()+scale_x_discrete(name="Days post-mating",labels=c("6","7","8","9","10"))
Recomb_figure=Recomb_figure+ggtitle("y-se Interval")+ylim(0.1,0.7)
Recomb_figure=Recomb_figure+annotate(geom="text", x=1, y=0.58, label=yse_sig[1],size=10)+annotate(geom="text", x=2, y=0.58, label=yse_sig[2],size=10)+annotate(geom="text", x=3, y=0.58, label=yse_sig[3],size=10)+annotate(geom="text", x=4, y=0.58, label=yse_sig[4],size=10)
Recomb_figure

dev.off()



# Reproducibility Experiments 2-3 ---------------------------------------------------------

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
Reproducibility <- data.frame(interval=c("sdy","yse","sdy","yse"),
                              odds=c(Exp4_odds[3,4],Exp4_odds[7,4],Exp5_odds[1,4],Exp5_odds[3,4]),
                              error=c(Exp4_odds[3,5],Exp4_odds[7,5],Exp5_odds[1,5],Exp5_odds[3,5]),
                              experiment=c("Exp4","Exp4","Exp5","Exp5"))

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

genotyping_haplotypes=rbind(sum(genotyping$SCO_1.sum,na.rm = T),
                      sum(genotyping$SCO_2.sum,na.rm = T),
                      sum(genotyping$SCO_3.sum,na.rm = T),
                      sum(genotyping$SCO_4.sum,na.rm = T),
                      sum(genotyping$SCO_5.sum,na.rm = T),
                      sum(genotyping$NCO_1.sum,na.rm = T),
                      sum(genotyping$NCO_2.sum,na.rm = T),
                      sum(genotyping$NCO_3.sum,na.rm = T),
                      sum(genotyping$NCO_4.sum,na.rm = T),
                      sum(genotyping$NCO_5.sum,na.rm = T))

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
pdf("SNP_genotyping_odds_with5colors.pdf")

odds_figure_SNP=ggplot(aes(y=OR,x=Day, col=Interval,group=Interval,shape=Interval),data=odds_ratios)+scale_colour_manual(values=c("#000000","#56B4E9","#009E73","#D55E00","#CC79A7"))+ 
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(-1.25,8)+theme(axis.text.x = element_text(angle = 45))+ggtitle("SNP genotyping")+
  scale_x_discrete(name="Days post-mating",labels = c("1-2","3-4","5-6","7-8","9-10"))+geom_errorbar(aes(ymin=OR-SE,ymax=OR+SE))+ geom_line(aes(linetype=Interval), size=1.5)+ #scale_linetype_manual("", values=c(1,2,1,2,3))+ 
  annotate(geom="text", x=1, y=6.6, label=i1_sig[1],color="#000000",size=10)+annotate(geom="text", x=2, y=6.6, label=i1_sig[2],color="#000000",size=10)+annotate(geom="text", x=3, y=6.6, label=i1_sig[3],color="#000000",size=10)+annotate(geom="text", x=4, y=6.6, label=i1_sig[4],color="#000000",size=10)+annotate(geom="text", x=5, y=6.6, label=i1_sig[5],color="#000000",size=10)+
  annotate(geom="text", x=1, y=6.7, label=i2_sig[1],color="#56B4E9",size=10)+annotate(geom="text", x=2, y=6.7, label=i2_sig[2],color="#56B4E9",size=10)+annotate(geom="text", x=3, y=6.7, label=i2_sig[3],color="#56B4E9",size=10)+annotate(geom="text", x=4, y=6.7, label=i2_sig[4],color="#56B4E9",size=10)+annotate(geom="text", x=5, y=6.7, label=i2_sig[5],color="#56B4E9",size=10)+
  annotate(geom="text", x=1, y=6.8, label=i3_sig[1],color="#009E73",size=10)+annotate(geom="text", x=2, y=6.8, label=i3_sig[2],color="#009E73",size=10)+annotate(geom="text", x=3, y=6.8, label=i3_sig[3],color="#009E73",size=10)+annotate(geom="text", x=4, y=6.8, label=i3_sig[4],color="#009E73",size=10)+annotate(geom="text", x=5, y=6.8, label=i3_sig[5],color="#009E73",size=10)+
  annotate(geom="text", x=1, y=6.9, label=i4_sig[1],color="#D55E00",size=10)+annotate(geom="text", x=2, y=6.9, label=i4_sig[2],color="#D55E00",size=10)+annotate(geom="text", x=3, y=6.9, label=i4_sig[3],color="#D55E00",size=10)+annotate(geom="text", x=4, y=6.9, label=i4_sig[4],color="#D55E00",size=10)+annotate(geom="text", x=5, y=6.9, label=i4_sig[5],color="#D55E00",size=10)+
  annotate(geom="text", x=1, y=7, label=i5_sig[1],color="#CC79A7",size=10)+annotate(geom="text", x=2, y=7, label=i5_sig[2],color="#CC79A7",size=10)+annotate(geom="text", x=3, y=7, label=i5_sig[3],color="#CC79A7",size=10)+annotate(geom="text", x=4, y=7, label=i5_sig[4],color="#CC79A7",size=10)+annotate(geom="text", x=5, y=7, label=i5_sig[5],color="#CC79A7",size=10)
  
  odds_figure_SNP
dev.off()

#Figure 3 for all the combined odds ratios
library(ggpubr)
postscript("Figure3.eps",width=4, height=3, horizontal=FALSE, pointsize=10)

#ggarrange(odds_figure_SNP,odds_figure_exp1,odds_figure_exp2,odds_figure_exp3,odds_figure_exp4,odds_figure_exp5, labels = c("a. SNP Genotyping","b. Experiment 1","c. Experiment 2","d. Experiment 3","e. Experiment 4","f. Experiment 5"), ncol = 1,nrow = 3)
pdf("figure3_snp.pdf")
ggarrange(odds_figure_SNP,                                                 # First row with scatter plot
          ggarrange(odds_figure_exp1, odds_figure_exp2,odds_figure_exp3, ncol = 3, labels = c("B", "C", "D")), # Second row with box and dot plots
          ggarrange(odds_figure_exp4,odds_figure_exp5, ncol = 2, labels = c("E","F")),
          nrow = 3, 
          labels = "A"                                        # Labels of the scatter plot
) 

dev.off()

#supplementary figure 1 for fecundity measurements
pdf("Supplementary_figure_1.pdf")
ggarrange(Fecund_figure_exp1,Fecund_figure_exp2,Fecund_figure_exp3,Fecund_figure_exp4,Fecund_figure_exp5, Fecund_figure_age,labels = c("A","B","C","D","E","F"), ncol = 2,nrow = 3)
dev.off()

pdf("SNP_i3_boxplot.pdf")
SNP_i3_boxplot=ggplot(aes(y=rec_rate_i3,x=Day, col=Treatment,label=count.sum),data=genotyping)+scale_colour_manual(values=c("blue","red"))+geom_boxplot()+ylab("% recombination")+theme_base()+scale_x_discrete(name="Days post-mating",labels=c("1-2","3-4","5-6","7-8","9-10"))
SNP_i3_boxplot=SNP_i3_boxplot+ggtitle("interval 3")+ylim(0,0.6)
SNP_i3_boxplot=SNP_i3_boxplot+annotate(geom="text", x=1, y=0.58, label=i3_sig[1],size=10)+annotate(geom="text", x=2, y=0.58, label=i3_sig[2],size=10)+annotate(geom="text", x=3, y=0.58, label=i3_sig[3],size=10)+annotate(geom="text", x=4, y=0.58, label=i3_sig[4],size=10 +annotate(geom="text",x=5,y=0.58,label=i3_sig[5],size=10))
SNP_i3_boxplot
dev.off()

pdf("SNP_i4_boxplot.pdf")
SNP_i4_boxplot=ggplot(aes(y=rec_rate_i4,x=Day, col=Treatment,label=count.sum),data=genotyping)+scale_colour_manual(values=c("blue","red"))+geom_boxplot()+ylab("% recombination")+theme_base()+scale_x_discrete(name="Days post-mating",labels=c("1-2","3-4","5-6","7-8","9-10"))
SNP_i4_boxplot=SNP_i4_boxplot+ggtitle("interval 4")+ylim(0,0.6)
SNP_i4_boxplot=SNP_i4_boxplot+annotate(geom="text", x=1, y=0.58, label=i4_sig[1],size=10)+annotate(geom="text", x=2, y=0.58, label=i4_sig[2],size=10)+annotate(geom="text", x=3, y=0.58, label=i4_sig[3],size=10)+annotate(geom="text", x=4, y=0.58, label=i4_sig[4],size=10+annotate(geom="text",x=5,y=0.58,label=i4_sig[5],size=10))
SNP_i4_boxplot
dev.off()

# Experiment 4 (Age Experiment) ----------------------------------------------------------

#This script includes analysis for the manuscript titled:
#"Maternal age alters recombination rate in Drosophila pseudoobscura"

#It is split into 7 code sections which can be expanded/collapsed in RStudio 
#Expand: either click on the arrow in the gutter or on the icon that overlays the folded code 
#Collapse: click on the arrow in the gutter

#The third section of code shows the steps to convert the raw data file "Mutant_screen_data_raw.csv" 
#to a cleaned-up version that is used in the following sections called "Mutant_screen_data_cleanedup.csv".

# Experiment 4 Section 1: Fecundity Pilot Data Analysis -------------------------------------------

fec=read.csv("FecundityPilotExperiment_updated.csv", header=T,na.strings = "na",stringsAsFactors = T)

fec$Treatment=as.factor(fec$Treatment)

#postscript ("../Figures/Figure.2A.eps", width=4, height=3, horizontal=FALSE, pointsize=10) 
#png("../Figures/Figure2A.png",res=300,height = 3.5,width = 4.5,units="in",pointsize = 10)
#plot(fec$Treatment,fec$Fecundity,xlab="Maternal Age (Days)", ylab="Fecundity")
#dev.off()

#boxplot(Fecundity ~ Treatment, data=fec )
pdf("fecunditypilot.pdf")
fec_boxplot = ggplot(data=fec, aes(x=Treatment,y=Fecundity))+geom_boxplot(show.legend = FALSE)+theme_base()+ylim(0,181)+
  labs(y="Fecundity", x="Maternal Age (Days)")+theme(legend.title = element_blank(),legend.position = "none")
fec_boxplot
dev.off()
#add colors: scale_fill_manual(values=c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8","#253494"))

hist(fec$Fecundity)
summary(fec$Fecundity)

model_full=aov(fec$Fecundity~fec$Treatment)
summary(model_full)
posthoc <- TukeyHSD(x=model_full, 'fec$Treatment', conf.level=0.95)
posthoc

tapply(fec$Fecundity,fec$Treatment,mean,na.rm=T)

# Experiment 4_Section 2: Survivorship Analysis ---------------------------------------------------

#read in data
sv=read.csv("../datasets/survivorship_per_rep_updated.csv",header=T)
colnames(sv)=c("F1 Vial","0","7","14","28","35","42","49","56","63","77","84","91")

#compute median survivorship
col0=median(sv$`0`)
col1=median(sv$`7`)
col2=median(sv$`14`)
col3=median(sv$`28`)
col4=median(sv$`35`)
col5=median(sv$`42`)
col6=median(sv$`49`)
col7=median(sv$`56`)
col8=median(sv$`63`)
col9=median(sv$`77`)
col10=median(sv$`84`)
col11=median(sv$`91`)

sv2=data.frame(x=c(0,7,14,28,35,42,49,56,63,77,84,91),y=c(col0,col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11))

#building the survivorship curve
postscript ("../Figures/Figure.3.eps", paper="special", width=4, height=3, horizontal=FALSE, pointsize=10)
#png("../Figures/Figure3.png",res=300,height = 3.5,width = 4.5,units="in",pointsize = 10)
plot(as.numeric(sv2$x),sv2$y,xlab="Age in Days",ylab="Median survival",type="l",lwd=3)

for (v in 1:18) {
  lines(as.numeric(sv2$x),sv[v,c(2:13)],lwd=1,col="grey")
}

lines(as.numeric(sv2$x),sv2$y,lwd=3,col="#018571")

abline(v=35, col="#a6611a")

text(70,0.9,paste("Median Survival = 42 days"),col="#018571",cex=0.75)
text(70,0.8,paste("Selected Treatment Age = 35 days"),col="#a6611a",cex=0.75)

dev.off()

# Experiment 4_Section 3: Mutant Screen Data Clean-up ---------------------------------------------

#read in raw data, combine with treatment data, and female counts
yv=read.csv("Mutant_screen_data_raw.csv",header=T)
unique(yv$Initials)
#read in cross data with treatment information
bc_worksheet=read.csv("Mutant_screen_treatment_data_raw.csv",header=T,stringsAsFactors = F) 
yv=na.omit(yv)
bc_worksheet$Treatment=as.character(bc_worksheet$Treatment)

#read in female count data
female_counts=read.csv(file="Mutant_screen_female_count_data.csv",header=T,stringsAsFactors = F)

#recode phenotypes as characters; mutant screen used 1 s and 0 s for scoring data
#conversion to character tells R that our treatment is not a number and to treat it as a character.
yv$sd=as.character(yv$sd)
yv$y=as.character(yv$y)
yv$se=as.character(yv$se)

#Define Crossovers
yv$co_class=ifelse(yv$sd==yv$y & yv$y==yv$se,"non_CO", 
                   ifelse(yv$sd!=yv$y & yv$y==yv$se,"single_CO_1",
                          ifelse(yv$sd==yv$y & yv$y!=yv$se,"single_CO_2",
                                 ifelse(yv$sd!=yv$y & yv$y!=yv$se,"double_CO",
                                        "error"))))

#check for errors, which are the removed cut phenotypes. This should ALWAYS be empty!
yv[yv$co_class=="error",]

#Define Crossovers
yv$gamete_class=ifelse(yv$sd==yv$y & yv$y==yv$se & yv$y==1,"gt1a",
                         ifelse(yv$sd==yv$y & yv$y==yv$se & yv$y==0,"gt1b",
                                ifelse(yv$sd!=yv$y & yv$y==yv$se & yv$y==1,"gt2a",
                                       ifelse(yv$sd!=yv$y & yv$y==yv$se & yv$y==0,"gt2b",
                                              ifelse(yv$sd==yv$y & yv$y!=yv$se & yv$y==1,"gt3a",
                                                     ifelse(yv$sd==yv$y & yv$y!=yv$se & yv$y==0,"gt3b",
                                                            ifelse(yv$sd!=yv$y & yv$y!=yv$se & yv$y==1,"gt4a",
                                                                   ifelse(yv$sd!=yv$y & yv$y!=yv$se & yv$y==0,"gt4b",
                                                                          "error"))))))))
#Sanity check, there should not be any error!
yv[yv$gamete_class=="error",]

#add columns for counting
yv$gt1a_111=ifelse(yv$gamete_class=="gt1a",yv$numbMales,0)
yv$gt1b_000=ifelse(yv$gamete_class=="gt1b",yv$numbMales,0)
yv$gt2a_011=ifelse(yv$gamete_class=="gt2a",yv$numbMales,0)
yv$gt2b_100=ifelse(yv$gamete_class=="gt2b",yv$numbMales,0)
yv$gt3a_110=ifelse(yv$gamete_class=="gt3a",yv$numbMales,0)
yv$gt3b_001=ifelse(yv$gamete_class=="gt3b",yv$numbMales,0)
yv$gt4a_010=ifelse(yv$gamete_class=="gt4a",yv$numbMales,0)
yv$gt4b_101=ifelse(yv$gamete_class=="gt4b",yv$numbMales,0)

gt1a_111=sum(yv$numbMales*yv$gt1a_111)
gt1b_000=sum(yv$numbMales*yv$gt1b_000)
gt2a_011=sum(yv$numbMales*yv$gt2a_011)
gt2b_100=sum(yv$numbMales*yv$gt2b_100)
gt3a_110=sum(yv$numbMales*yv$gt3a_110)
gt3b_001=sum(yv$numbMales*yv$gt3b_001)
gt4a_010=sum(yv$numbMales*yv$gt4a_010)
gt4b_101=sum(yv$numbMales*yv$gt4b_101)

age_haplotype_types=cbind(gt1a_111,gt1b_000,gt2a_011,gt2b_100,gt3a_110,gt3b_001,gt4a_010,gt4b_101)
write.csv(age_haplotype_types,"age_haplotype_types.csv")

#add columns for counting
yv$NCO=ifelse(yv$co_class=="non_CO",yv$numbMales,0)
yv$SCO_1=ifelse(yv$co_class=="single_CO_1",yv$numbMales,0)
yv$SCO_2=ifelse(yv$co_class=="single_CO_2",yv$numbMales,0)
yv$DCO=ifelse(yv$co_class=="double_CO",yv$numbMales,0)

#get rough crossover rate as a sanity check the numbers of crossovers at the intervals 
#should equal to the total crossover rate.
nco_count=sum(yv$NCO)
sco_count=sum(yv$SCO_1, na.rm = TRUE)+sum(yv$SCO_2, na.rm = TRUE)
dco_count=sum(yv$DCO, na.rm = TRUE)
num_samples=sum(nco_count+sco_count+dco_count, na.rm = TRUE) 
num_samples

#total rate
(sco_count+(2*dco_count))/num_samples

#rate between scalloped and yellow
(sum(yv$SCO_1,na.rm = TRUE)+sum(yv$DCO,na.rm = TRUE))/num_samples

#rate between sepia and yellow
(sum(yv$SCO_2, na.rm = TRUE)+sum(yv$DCO,na.rm = TRUE))/num_samples

#we defined the crossover groups as if non crossovers will be identified as 0s; this is same as co_class column, but numerical
#single crossovers 1s
#double crossovers 2s
yv$num_co=ifelse(yv$y==yv$sd & yv$y==yv$se,0, 
                 ifelse(yv$sd==yv$y & yv$y!=yv$se,1*yv$numbMales, 
                        ifelse(yv$sd!=yv$y & yv$y==yv$se,1*yv$numbMales,  
                               ifelse(yv$sd!=yv$y & yv$y!=yv$se,2*yv$numbMales, 
                                      NA))))

#This is for summarizing our data
yv$male=c(yv$numbMales)

#merge with treatment data
yv_merged <- merge(yv, bc_worksheet, by.x = "ViaNumber", by.y = "Vial.Number", all=T)

#summarize long form data, include the gamete classes in the summary
dataset=summaryBy(male+gt1a_111+gt1b_000+gt2a_011+gt2b_100+gt3a_110+gt3b_001+gt4a_010+gt4b_101+num_co+SCO_1+SCO_2+DCO~F1.Vial+Day+Treatment,data=yv_merged, FUN=sum,na.rm=T)


#merge female data with treatment data
female_merged=merge(female_counts, bc_worksheet, by.x = "Vial", by.y = "Vial.Number", all=T)

#summarize female data by day, replicate and treatment
female_short=summaryBy(Numbfemales~F1.Vial+Day+Treatment,data=female_merged, FUN=sum,na.rm=T)

#get rid of weird NAs
female_short=na.omit(female_short)

#final merge combines male and female datasets
dataset2=merge(female_short,dataset, by=c("F1.Vial","Day","Treatment"))

#add in a column for total offspring
dataset2$total_offspring=dataset2$Numbfemales.sum + dataset2$male.sum

#total sample size
sum(dataset2$total_offspring)

#now to calculate fecundity, we need to add a column to report how many individual backcrosses were conducted per replicate vial
#first make a vector to store the data
num_moms=vector(mode="numeric",length=length(dataset2$F1.Vial))

#loop through dataset to get the cross count per rep
for (h in 1:length(dataset2$F1.Vial)) { 
  f1_vial=dataset2$F1.Vial[h]
  mom_ct=length(unique(sort(subset(yv_merged,yv_merged$F1.Vial==f1_vial)$ViaNumber)))
  #store result in the vector
  num_moms[h]=mom_ct
}

#add vector as a column in dataset2
dataset2$Num_moms=num_moms

#use data to get fecundity calculation per replicate
dataset2$fecundity=dataset2$total_offspring/dataset2$Num_moms

#We will write our data into a file so we can read it in for later analysis. 
write.csv(dataset2,file="Mutant_screen_data_cleanedup.csv")


for (n in 1:length(dataset2$F1.Vial)) {
  dataset2$NCO_skew[n]=min(dataset2$gt1a_111.sum[n],dataset2$gt1b_000.sum[n])/max(dataset2$gt1a_111.sum[n],dataset2$gt1b_000.sum[n])
}
summary(dataset2$NCO_skew)


for (n in 1:length(dataset2$F1.Vial)) {
  dataset2$SCO1_skew[n]=min(dataset2$gt2a_011.sum[n],dataset2$gt2b_100.sum[n])/max(dataset2$gt2a_011.sum[n],dataset2$gt2b_100.sum[n])
}


for (n in 1:length(dataset2$F1.Vial)) {
  dataset2$SCO2_skew[n]=min(dataset2$gt3a_110.sum[n],dataset2$gt3b_001.sum[n])/max(dataset2$gt3a_110.sum[n],dataset2$gt3b_001.sum[n])
}

for (n in 1:length(dataset2$F1.Vial)) {
  dataset2$DCO_skew[n]=min(dataset2$gt4a_010.sum[n],dataset2$gt4b_101.sum[n])/max(dataset2$gt4a_010.sum[n],dataset2$gt4b_101.sum[n])
}

skew=dataset2[, c(3,21,22,23,24)]

skew2=melt(skew, id="Treatment")


pdf("skew_expage.pdf")
skew_figure_age=ggplot(skew2, aes(x = Treatment, y = value, color = variable)) +  # ggplot function
  geom_boxplot()+theme_base()+ylab("haplotype bias")+xlab("Treatment")

skew_figure_age
dev.off()

#gamete analysis by treatment
gt111=tapply(dataset2$gt1a_111.sum,dataset2$Treatment,sum,na.rm=T)
gt000=tapply(dataset2$gt1b_000.sum,dataset2$Treatment,sum,na.rm=T)
gt011=tapply(dataset2$gt2a_011.sum,dataset2$Treatment,sum,na.rm=T)
gt100=tapply(dataset2$gt2b_100.sum,dataset2$Treatment,sum,na.rm=T)
gt110=tapply(dataset2$gt3a_110.sum,dataset2$Treatment,sum,na.rm=T)
gt001=tapply(dataset2$gt3b_001.sum,dataset2$Treatment,sum,na.rm=T)
gt010=tapply(dataset2$gt4a_010.sum,dataset2$Treatment,sum,na.rm=T)
gt101=tapply(dataset2$gt4b_101.sum,dataset2$Treatment,sum,na.rm=T)
i=rbind(gt111, gt000, gt011, gt100, gt110, gt001, gt010, gt101)
write.csv(i,"age_gamete_class_by_treatment.csv")
#gamete analysis by day
gt111=tapply(dataset2$gt1a_111.sum,dataset2$Day,sum,na.rm=T)
gt000=tapply(dataset2$gt1b_000.sum,dataset2$Day,sum,na.rm=T)
gt011=tapply(dataset2$gt2a_011.sum,dataset2$Day,sum,na.rm=T)
gt100=tapply(dataset2$gt2b_100.sum,dataset2$Day,sum,na.rm=T)
gt110=tapply(dataset2$gt3a_110.sum,dataset2$Day,sum,na.rm=T)
gt001=tapply(dataset2$gt3b_001.sum,dataset2$Day,sum,na.rm=T)
gt010=tapply(dataset2$gt4a_010.sum,dataset2$Day,sum,na.rm=T)
gt101=tapply(dataset2$gt4b_101.sum,dataset2$Day,sum,na.rm=T)
s=rbind(gt111, gt000, gt011, gt100, gt110, gt001, gt010, gt101)
write.csv(s,"age_gamete_class_byday.csv")
gt111=tapply(dataset2$gt1a_111.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt000=tapply(dataset2$gt1b_000.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt011=tapply(dataset2$gt2a_011.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt100=tapply(dataset2$gt2b_100.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt110=tapply(dataset2$gt3a_110.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt001=tapply(dataset2$gt3b_001.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt010=tapply(dataset2$gt4a_010.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt101=tapply(dataset2$gt4b_101.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
t=rbind(gt111, gt000, gt011, gt100, gt110, gt001, gt010, gt101)
write.csv(t,"age_gamete_class_byday_and_treatment.csv")
age_pvalues_111_000_age=cbind(binom.test(c(gt111[1,1],gt000[1,1]),p=0.5)[3],
                              binom.test(c(gt111[1,2],gt000[1,2]),p=0.5)[3],
                              binom.test(c(gt111[1,3],gt000[1,3]),p=0.5)[3],
                              binom.test(c(gt111[1,4],gt000[1,4]),p=0.5)[3])

#20 +yse and sd++
age_pvalues_011_100_age=cbind(binom.test(c(gt011[1,1],gt100[1,1]),p=0.5)[3],
                              binom.test(c(gt011[1,2],gt100[1,2]),p=0.5)[3],
                              binom.test(c(gt011[1,3],gt100[1,3]),p=0.5)[3],
                              binom.test(c(gt011[1,4],gt100[1,4]),p=0.5)[3])

#20 sdy+ and ++se
age_pvalues_110_001_age=cbind(binom.test(c(gt110[1,1],gt001[1,1]),p=0.5)[3],
                              binom.test(c(gt110[1,2],gt001[1,2]),p=0.5)[3],
                              binom.test(c(gt110[1,3],gt001[1,3]),p=0.5)[3],
                              binom.test(c(gt110[1,4],gt001[1,4]),p=0.5)[3])

#20 +y+ and sd+se
age_pvalues_010_101_age=cbind(binom.test(c(gt010[1,1],gt101[1,1]),p=0.5)[3],
                              binom.test(c(gt010[1,2],gt101[1,2]),p=0.5)[3],
                              binom.test(c(gt010[1,3],gt101[1,3]),p=0.5)[3],
                              binom.test(c(gt010[1,4],gt101[1,4]),p=0.5)[3])


#26 sdyse and +++
age_pvalues_111_000_control=cbind(binom.test(c(gt111[2,1],gt000[2,1]),p=0.5)[3],
                              binom.test(c(gt111[2,2],gt000[2,2]),p=0.5)[3],
                              binom.test(c(gt111[2,3],gt000[2,3]),p=0.5)[3],
                              binom.test(c(gt111[2,4],gt000[2,4]),p=0.5)[3])


#26 +yse and sd++
age_pvalues_011_100_control=cbind(binom.test(c(gt011[2,1],gt100[2,1]),p=0.5)[3],
                              binom.test(c(gt011[2,2],gt100[2,2]),p=0.5)[3],
                              binom.test(c(gt011[2,3],gt100[2,3]),p=0.5)[3],
                              binom.test(c(gt011[2,4],gt100[2,4]),p=0.5)[3])


#26 sdy+ and ++se
age_pvalues_110_001_control=cbind(binom.test(c(gt110[2,1],gt001[2,1]),p=0.5)[3],
                              binom.test(c(gt110[2,2],gt001[2,2]),p=0.5)[3],
                              binom.test(c(gt110[2,3],gt001[2,3]),p=0.5)[3],
                              binom.test(c(gt110[2,4],gt001[2,4]),p=0.5)[3])


#26 +y+ and sd+se
age_pvalues_010_101_control=cbind(binom.test(c(gt010[2,1],gt101[2,1]),p=0.5)[3],
                              binom.test(c(gt010[2,2],gt101[2,2]),p=0.5)[3],
                              binom.test(c(gt010[2,3],gt101[2,3]),p=0.5)[3],
                              binom.test(c(gt010[2,4],gt101[2,4]),p=0.5)[3])
age_pvalues_haplotype_types_bydayandtreatment=rbind(age_pvalues_111_000_age,age_pvalues_011_100_age,age_pvalues_110_001_age,age_pvalues_010_101_age,age_pvalues_111_000_control,age_pvalues_011_100_control,age_pvalues_110_001_control,age_pvalues_010_101_control)
colnames(age_pvalues_haplotype_types_bydayandtreatment)=c("1-3","4-6","7-9","10-12")
rownames(age_pvalues_haplotype_types_bydayandtreatment)=c("sdyse_+++_age","+yse_sd++_age","sdy+_++se_age","+y+_sd+se_age","sdyse_+++_control","+yse_sd++_control","sdy+_++se_control","+y+_sd+se_control")
write.csv(age_pvalues_haplotype_types_bydayandtreatment,"age_pvalues_haplotype_types_bydayandtreatment.csv")



#Summarizing the data

length(unique(dataset2$F1.Vial[dataset2$Treatment=="Control"]))
length(unique(dataset2$F1.Vial[dataset2$Treatment=="Age"]))

median(dataset2$Num_moms[dataset2$Treatment=="Control"])
median(dataset2$Num_moms[dataset2$Treatment=="Age"])

sum(dataset2$male.sum[dataset2$Treatment=="Control"])
sum(dataset2$male.sum[dataset2$Treatment=="Age"])

# Experiment 4_Section 4: Fecundity Analysis of Mutant Screen -------------------------------------

#Read in cleaned up data
dataset2=read.csv("Mutant_screen_data_cleanedup.csv")

###Fecundity Model 

#poisson regression, similar to a t-test for count data
fit=glm(fecundity~Treatment*Day,data=dataset2,family=quasipoisson)
summary(fit)
anova_fec=anova(fit, test="Chisq")
anova_fec
write.csv(anova_fec,"age_fecundity_model_table.csv")

fit_contrast <- emmeans::emmeans(fit, "Treatment", by="Day", mode="kenward-roger")
fit_contr <- contrast(fit_contrast, method="trt.vs.ctrl")

pheno_contr <- as.data.frame(summary(fit_contr))
pheno_contr
write.csv(pheno_contr,"fecundity_posthoc_table.csv")

#convert p-values to stars for plot
sig=ifelse(pheno_contr$p.value<0.001,"***",ifelse(pheno_contr$p.value<0.01,"**",ifelse(pheno_contr$p.value<0.05,"*","")))

#draw line plot
Fecund_figure=ggplot(aes(y=fecundity,x=Day, col=Treatment,label=Num_moms),data=dataset2)+ylab("Fecundity")+theme_base()
Fecund_figure=Fecund_figure+scale_colour_manual(values=c("#7b3294","#008837"))+geom_point(size=2)+scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))+theme(legend.title = element_text(size = 10),legend.text = element_text(size = 8))
Fecund_figure=Fecund_figure+stat_summary(fun = median, geom="line",aes(group=Treatment),size=1.5)+geom_text(check_overlap = F,hjust = 0, nudge_x = 0.1,angle=45,size=3)+ylim(0,181)+
  annotate(geom="text", x=1, y=150, label=sig[1],size=7)+annotate(geom="text", x=2, y=150, label=sig[2],size=7)+annotate(geom="text", x=3, y=150, label=sig[3],size=7)+annotate(geom="text", x=4, y=150, label=sig[4],size=7)
Fecund_figure

#Fecundity figure for the paper
#postscript ("../Figures/Figure.2B.eps", width=4, height=3, horizontal=FALSE, pointsize=4)
#png("../Figures/Figure2B.png")
pdf("age_fecundity.pdf")
Fecund_figure_age=ggplot(aes(y=fecundity,x=Day, col=Treatment,label=Num_moms),data=dataset2)+ylab("Fecundity")+theme_base()
Fecund_figure_age=Fecund_figure_age+scale_colour_manual(values=c("#7b3294","#008837"))+geom_boxplot()+scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))+theme(legend.title = element_text(size = 10),legend.text = element_text(size = 8))
Fecund_figure_age=Fecund_figure_age+ylim(0,181)+
  annotate(geom="text", x=1, y=150, label=sig[1],size=7)+annotate(geom="text", x=2, y=150, label=sig[2],size=7)+annotate(geom="text", x=3, y=150, label=sig[3],size=7)+annotate(geom="text", x=4, y=150, label=sig[4],size=7)
Fecund_figure_age

dev.off()


# Experiment 4_Make multi-panel fecundity figure ---------------------------------------


#now print combined figure 2
#postscript ("../Figures/Figure.2.eps", width=8, height=4, horizontal=FALSE, pointsize=5)
jpeg("../Figures/Figure.2.jpg", width=8, height=4, pointsize=5,units="in",res=300)
ggarrange(fec_boxplot,Fecund_figure,nrow=1,ncol=2,labels=c("A","B"),widths=c(1,1.75))
dev.off()

#get averages for text; These are the numbers we put in the paper!
tapply(dataset2$fecundity,dataset2$Treatment,mean,na.rm=T)
tapply(dataset2$fecundity,dataset2$Day,mean,na.rm=T)


# Experiment 4_Section 5: Recombination Analysis for Mutant Screen --------------------------------

#Read in cleaned up data
dataset2=read.csv("Mutant_screen_data_cleanedup.csv")

#remove vials with too few progeny
dataset2=dataset2[dataset2$male.sum>=10,]
dataset2=dataset2[dataset2$Numbfemales.sum>=10,]

#sum of crossovers in intervals 1 (sd-y) & 2 (y-se)
dataset2$num_CO_1=dataset2$SCO_1.sum+dataset2$DCO.sum
dataset2$num_CO_2=dataset2$SCO_2.sum+dataset2$DCO.sum

#sum of non-crossovers in intervals 1 & 2
dataset2$num_NCO_1=dataset2$male.sum-(dataset2$SCO_1.sum+dataset2$DCO.sum)
dataset2$num_NCO_2=dataset2$male.sum-(dataset2$SCO_2.sum+dataset2$DCO.sum)

#total recombination rate
dataset2$rec_rate_total=dataset2$num_co.sum/dataset2$male.sum
dataset2$rec_rate_ysd=(dataset2$SCO_1.sum+dataset2$DCO.sum)/dataset2$male.sum
dataset2$rec_rate_yse=(dataset2$SCO_2.sum+dataset2$DCO.sum)/dataset2$male.sum

#Recomb_figure_total
Recomb_figure_total=ggplot(aes(y=rec_rate_total,x=Day, col=Treatment,label=male.sum),data=dataset2)+ylab("% recombination")+ggtitle("Total Recombination rate vs. Days post-mating")+theme_base()

#Now add points and change labels for x-axis
Recomb_figure_total=Recomb_figure_total+scale_colour_manual(values=c("#7b3294","#008837"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))

#Add a line through the median of the points
Recomb_figure_total=Recomb_figure_total+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)

#Add a label for sample size
Recomb_figure_total=Recomb_figure_total+ggtitle("Total Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3) +ylim(0.2,1)

#print figure
Recomb_figure_total

#get averages for text; These are the numbers we put in the paper!
tapply(dataset2$rec_rate_total,dataset2$Treatment,mean,na.rm=T)
#tapply(dataset2$rec_rate_total,dataset2$Day,mean,na.rm=T)

tapply(dataset2$rec_rate_ysd,dataset2$Treatment,mean,na.rm=T)
#tapply(dataset2$rec_rate_ysd,dataset2$Day,mean,na.rm=T)

tapply(dataset2$rec_rate_yse,dataset2$Treatment,mean,na.rm=T)
#tapply(dataset2$rec_rate_yse,dataset2$Day,mean,na.rm=T)



write.csv(dataset2,"experimentage_combined dataset.csv")

#Recomb_figure_yellow_scalloped
Recomb_figure_ysd=ggplot(aes(y=rec_rate_ysd,x=Day, col=Treatment,label=male.sum),data=dataset2)+scale_colour_manual(values=c("#7b3294","#008837"))+geom_point(alpha=0.6,size=3)+ylab("% recombination")+theme_base()+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))
Recomb_figure_ysd=Recomb_figure_ysd+ggtitle("y-sd Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
Recomb_figure_ysd

#Recomb_figure_yellow_sepia
Recomb_figure_yse=ggplot(aes(y=rec_rate_yse,x=Day, col=Treatment,label=male.sum),data=dataset2)+scale_colour_manual(values=c("#7b3294","#008837"))+geom_point(alpha=0.6,size=3)+ylab("% recombination")+theme_base()+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))
Recomb_figure_yse=Recomb_figure_yse+ggtitle("y-se Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
Recomb_figure_yse



# Experiment 4_Section 6: Recombination Rate Model & Odds Ratio Analysis ---------------------------------

#Recombination Rate Models per regions
#NOTE: YOU MUST RUN SECTION 5 CODE before this code section for it to work!

#Add new column to make treatment vs control sort properly. Odd ratios are alphabetical and thus age sorts first. 
#Change to "Zage" and it will sort properly to get odd ratios for Figure 4A

dataset2$NewTreatment=ifelse(dataset2$Treatment=="Age","Zage","Control")

#SD-Y REGION

#logistic regression, similar to a t-test for count data
fit3=glmer(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+NewTreatment*Day,data=dataset2,
           family=binomial(link="logit"))
#coefs=coef(fit3)
#coefs
anova_sdy=Anova(fit3,test="Chisq")
write.csv(anova_sdy,"age_recrate_sd-y_model_table.csv")

fit_contrast3 <- emmeans::emmeans(fit3, "NewTreatment", by="Day", mode="kenward-roger")
fit_contr3 <- contrast(fit_contrast3, method="trt.vs.ctrl")

pheno_contr3 <- as.data.frame(summary(fit_contr3))
pheno_contr3
write.csv(pheno_contr3,"age_recrate_sd-y_posthoc_table.csv")

#need odds ratio and standard error
#can extract from model for each time point
fit3a=glm(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+NewTreatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day=="A"))
exp(coef(fit3a)[3]) #extract odds ratio for sd-y at time point 1

#Sanity Check: this is the same as extracting the exp of the estimate of the posthoc table!
exp(pheno_contr3$estimate[1])

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
fit4=glmer(cbind(num_CO_2,num_NCO_2)~(1|F1.Vial)+NewTreatment*Day,data=dataset2,
           family=binomial(link="logit"))
summary(fit4)
#coefs=coef(fit4)
#coefs
anova_yse=Anova(fit4,test="Chisq")
write.csv(anova_yse,"age_recrate_y-se_model_table.csv")

fit_contrast4 <- emmeans::emmeans(fit4, "NewTreatment", by="Day", mode="kenward-roger")
fit_contr4 <- contrast(fit_contrast4, method="trt.vs.ctrl")

pheno_contr4 <- as.data.frame(summary(fit_contr4))
pheno_contr4
write.csv(pheno_contr4,"age_recrate_y-se_posthoc_table.csv")

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

#convert p-values to stars for plot
ysd_sig=ifelse(pheno_contr3$p.value<0.001,"***",ifelse(pheno_contr3$p.value<0.01,"**",ifelse(pheno_contr3$p.value<0.05,"*","")))
yse_sig=ifelse(pheno_contr4$p.value<0.001,"***",ifelse(pheno_contr4$p.value<0.01,"**",ifelse(pheno_contr4$p.value<0.05,"*","")))

#add significance to table
sig=c(ysd_sig,yse_sig)
y=cbind(x,sig)
odds_ratios=y
odds_ratios

#Odds ratio plot
#postscript ("../Figures/Figure.5.eps", width=4, height=3, horizontal=FALSE, pointsize=10)
#jpeg("../Figures/Figure.5.jpg", width=4, height=3, pointsize=10,units="in",res=300)
pdf("age_odds_ratio.pdf")
odds_figure=ggplot(aes(y=OR,x=Day, col=Interval,group=Interval),data=odds_ratios)+scale_colour_manual(values=c("#f1a340", "#998ec3"))+
  geom_point(size=1.5)+geom_line(size=2)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.75,1.3)+
  scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))+geom_line()+geom_errorbar(aes(ymin=OR-SE,ymax=OR+SE))+ggtitle("Recombination Rate Model")+theme(plot.title = element_text(size = 15))+
  annotate(geom="text", x=1, y=1.25, label=ysd_sig[1],color="#f1a340",size=5)+annotate(geom="text", x=2, y=1.25, label=ysd_sig[2],color="#f1a340",size=5)+annotate(geom="text", x=3, y=1.25, label=ysd_sig[3],color="#f1a340",size=5)+annotate(geom="text", x=4, y=1.25, label=ysd_sig[4],color="#f1a340",size=5)+
  annotate(geom="text", x=1, y=1.3, label=yse_sig[1],color="#998ec3",size=5)+annotate(geom="text", x=2, y=1.3, label=yse_sig[2],color="#998ec3",size=5)+annotate(geom="text", x=3, y=1.3, label=yse_sig[3],color="#998ec3",size=5)+annotate(geom="text", x=4, y=1.3, label=yse_sig[4],color="#998ec3",size=5)
odds_figure
dev.off()

#get percent difference for time point 1-3
tp1=subset(dataset2,dataset2$Day=="A")
s=tapply(tp1$rec_rate_ysd,tp1$Treatment,mean,na.rm=T)
100*(s[1]-s[2])

se=tapply(tp1$rec_rate_yse,tp1$Treatment,mean,na.rm=T)
100*(se[1]-se[2])



# Experiment 4_Section 7: Crossover interference analysis -----------------------------------------

#Read in cleaned up data
dataset2=read.csv("Mutant_screen_data_cleanedup.csv")

#remove vials with too few progeny
dataset2=dataset2[dataset2$male.sum>=10,]
dataset2=dataset2[dataset2$Numbfemales.sum>=10,]

dataset2$Exp_DCO=(((dataset2$SCO_1.sum+dataset2$DCO.sum)/dataset2$male.sum)*((dataset2$SCO_2.sum+dataset2$DCO.sum)/dataset2$male.sum))
dataset2$Obs_DCO=dataset2$DCO.sum/dataset2$male.sum

dataset2$COC=dataset2$Obs_DCO/dataset2$Exp_DCO
dataset2$Interference=1-dataset2$COC

#model
fit9=glm(Interference~Treatment*Day,data=dataset2)
summary(fit9)
anova_coi=anova(fit9, test="Chisq")
anova_coi
write.csv(anova_coi,"interference_model_table.csv")


fit_contrast9 <- emmeans::emmeans(fit9, "Treatment", by="Day", mode="kenward-roger")
fit_contr9 <- contrast(fit_contrast9, method="trt.vs.ctrl")

pheno_contr9 <- as.data.frame(summary(fit_contr9))
pheno_contr9
write.csv(pheno_contr9,"../stats/interference_posthoc_table.csv")

or=exp(pheno_contr9$estimate)
or

#convert p-values to stars for plot
sig2=ifelse(pheno_contr9$p.value<0.001,"***",ifelse(pheno_contr9$p.value<0.01,"**",ifelse(pheno_contr9$p.value<0.05,"*","")))

#Interference figure for the paper
#postscript ("../Figures/Figure.4B.eps", width=4, height=3, horizontal=FALSE, pointsize=5)
#png("../Figures/Figure4B.png")
COI_figure=ggplot(aes(y=Interference,x=Day, col=Treatment,label=Num_moms),data=dataset2)+ylab("COI")+ggtitle("Interference")+theme_base()+ylim(-0.25,0.4)
COI_figure=COI_figure+scale_colour_manual(values=c("#7b3294","#008837"))+geom_point(size=1.5)+scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))
COI_figure=COI_figure+stat_summary(fun = median, geom="line",aes(group=Treatment),size=1)+
  annotate(geom="text", x=1, y=0.25, label=sig2[1],size=5)+annotate(geom="text", x=2, y=0.25, label=sig2[2],size=5)+annotate(geom="text", x=3, y=0.25, label=sig2[3],size=5)+annotate(geom="text", x=4, y=0.25, label=sig2[4],size=5)
COI_figure

#add sample sizes to figure with the code below:
#geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=1)+

#dev.off()



# Experiment 4_Make multi-panel Recombination Rate and COI Figure -------------------------------------------------


#Combine Recombination Figure and COI Figure into multipanel figure:

#add columns for rec rate again
#sum of crossovers in intervals 1 (sd-y) & 2 (y-se)
dataset2$num_CO_1=dataset2$SCO_1.sum+dataset2$DCO.sum
dataset2$num_CO_2=dataset2$SCO_2.sum+dataset2$DCO.sum

#sum of non-crossovers in intervals 1 & 2
dataset2$num_NCO_1=dataset2$male.sum-(dataset2$SCO_1.sum+dataset2$DCO.sum)
dataset2$num_NCO_2=dataset2$male.sum-(dataset2$SCO_2.sum+dataset2$DCO.sum)

#total recombination rate
dataset2$rec_rate_total=dataset2$num_co.sum/dataset2$male.sum
dataset2$rec_rate_ysd=(dataset2$SCO_1.sum+dataset2$DCO.sum)/dataset2$male.sum
dataset2$rec_rate_yse=(dataset2$SCO_2.sum+dataset2$DCO.sum)/dataset2$male.sum


#Kosambi correction
dataset2$kosambi_rec_rate_ysd=((2.71^(4*dataset2$rec_rate_ysd)-1)/2*(2.71^(-4*dataset2$rec_rate_ysd)+1))/10
dataset2$kosambi_rec_rate_yse=((2.71^(4*dataset2$rec_rate_yse)-1)/2*(2.71^(-4*dataset2$rec_rate_yse)+1))/10



#Redraw recombination plots:

#Recomb_figure_yellow_scalloped
pdf("age_sdyboxplot.pdf")
Recomb_figure_ysd=ggplot(aes(y=kosambi_rec_rate_ysd,x=Day, col=Treatment,label=male.sum),data=dataset2)+scale_colour_manual(values=c("#7b3294","#008837"))+geom_boxplot()+ylab("% recombination")+theme_base()+scale_x_discrete(name="Days post-mating",labels=c("1-3","4-6","7-9","10-12"))
Recomb_figure_ysd=Recomb_figure_ysd+ggtitle("sd-y Interval")+ylim(0,0.6)
Recomb_figure_ysd=Recomb_figure_ysd+annotate(geom="text", x=1, y=0.58, label=ysd_sig[1],size=10)+annotate(geom="text", x=2, y=0.58, label=ysd_sig[2],size=10)+annotate(geom="text", x=3, y=0.58, label=ysd_sig[3],size=10)+annotate(geom="text", x=4, y=0.58, label=ysd_sig[4],size=10)
Recomb_figure_ysd
dev.off()
#Recomb_figure_yellow_sepia
pdf("age_yseboxplot.pdf")
Recomb_figure_yse=ggplot(aes(y=kosambi_rec_rate_yse,x=Day, col=Treatment,label=male.sum),data=dataset2)+scale_colour_manual(values=c("#7b3294","#008837"))+geom_boxplot()+ylab("% recombination")+theme_base()+scale_x_discrete(name="Days post-mating",labels=c("1-3","4-6","7-9","10-12"))
Recomb_figure_yse=Recomb_figure_yse+ggtitle("y-se Interval")+ylim(0,0.6)
Recomb_figure_yse=Recomb_figure_yse+annotate(geom="text", x=1, y=0.58, label=yse_sig[1],size=10)+annotate(geom="text", x=2, y=0.58, label=yse_sig[2],size=10)+annotate(geom="text", x=3, y=0.58, label=yse_sig[3],size=10)+annotate(geom="text", x=4, y=0.58, label=yse_sig[4],size=10)
Recomb_figure_yse
dev.off()
#Redraw COI plot:

#COI Figure
COI_figure=ggplot(aes(y=Interference,x=Day, col=Treatment,label=Num_moms),data=dataset2)+ylab("COI")+ggtitle("Interference")+theme_base()+ylim(-0.25,0.4)
COI_figure=COI_figure+scale_colour_manual(values=c("#7b3294","#008837"))+geom_boxplot()+scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))
COI_figure=COI_figure+
  annotate(geom="text", x=1, y=0.25, label=sig2[1],size=10)+annotate(geom="text", x=2, y=0.25, label=sig2[2],size=10)+annotate(geom="text", x=3, y=0.25, label=sig2[3],size=10)+annotate(geom="text", x=4, y=0.25, label=sig2[4],size=10)
COI_figure


#make multipanel plot
#postscript ("../Figures/Figure.4.eps", width=8, height=3, horizontal=FALSE, pointsize=2, paper="special")
jpeg("../Figures/Figure.4.jpg", width=7, height=3.5, pointsize=1,units="in",res=300)
ggarrange(Recomb_figure_ysd,Recomb_figure_yse,COI_figure,nrow=1,ncol=3,labels=c("A","B","C"),common.legend = T,legend="bottom")
dev.off()


