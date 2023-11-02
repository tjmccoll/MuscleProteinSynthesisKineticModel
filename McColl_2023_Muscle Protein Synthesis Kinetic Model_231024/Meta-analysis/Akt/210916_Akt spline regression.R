# clear data
rm(list=ls())

library(caret)
library(mgcv)
library(dplyr)
library(ggplot2)

## Load data
setwd("/Users/taylormccoll/Documents/SFU/1. Ph.D./1. Research Projects/1. Network Feedback/6. Manuscript/Manuscript Versions/230912_revision/Manuscript/MuscleProteinSynthesisKineticModel/McColl_2023_Muscle Protein Synthesis Kinetic Model_230919/Meta-analysis/Akt")
data.set <- read.csv("210916_Akt analysis, young adults.csv",header=T,sep=',')

# data.set = data.set[!(data.set$Lab.Code==1),] # removing data from Rasmussen's lab
data.set = data.set[!(data.set$Exercise.code==1),]

data.set$time = data.set$Measurement.Time..min.
data.set$foldChange = data.set$Fold.Change.from.Baseline
data.set$dose = data.set$Total..Leucine..g.
exercised = data.set$Exercised.

# Visualize
ggplot(data.set, aes(time, foldChange) ) +
  geom_point() +
  stat_smooth()

#Generalized additive model - no rasmussen, no exercise
predict.time = seq(0, 300, 30)
predict.dose = 3.5
new.data = data.frame(time = predict.time, dose = rep(predict.dose, each=length(predict.time)))

# Spline regression (5 knots)
knots=5
model = gam(foldChange ~ s(time, k=knots) + s(dose, k=knots), data = data.set)
summary(model)
#---

prd <- data.frame(time = seq(from = range(data.set$time)[1], to = range(data.set$time)[2], length.out = 301),
                  dose = rep(3.5, each=301))
err <- predict.gam(model, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

Akt_analysis = ggplot(prd, aes(x = time, y = fit)) +
  theme_bw() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity", color = "black") +
  xlab("Time (minutes)") + 
  ylab("Fold Change") +
  scale_y_continuous(minor_breaks=NULL, breaks=seq(-0.5, 3.5, 0.5)) +
  scale_x_continuous(minor_breaks=NULL, breaks=seq(0, 300, 60)) +
  geom_point(data=data.set, mapping = aes(x=Measurement.Time..min., y=foldChange, size=dose, color=Lab), 
             position=position_jitter(height=.1, width=10), alpha=0.5) +
  scale_size_continuous(range=c(2.5, 10)) +
  geom_hline(yintercept=1, linetype="dashed", color="gray50") +
  labs(size = "Leucine Dose (g)") +
  labs(color = "Lab Group") 
 Akt_analysis

Akt_analysis_theme = Akt_analysis + theme(legend.position = c(0.98, 0.98), 
                                                legend.text = element_text(size=8), 
                                                legend.title = element_text(size=8, face="bold"), 
                                                legend.box = "horizontal",
                                                legend.justification = c("right","top"),
                                                legend.background = element_blank(),
                                                legend.box.background = element_rect(fill ="white", colour = "black")) +
                                    guides(size = guide_legend(order=2), col = guide_legend(order=1))

FileName = paste("Akt Analysis - Spline Regression, no exercise - ", format(Sys.time(), "%m-%d-%y %X"), ".pdf")
ggsave(FileName,  Akt_analysis_theme, 
       device = "pdf",
       width = 16.5*0.75, height = 16.5*0.75, unit = "cm",
       path = "/Users/taylormccoll/Documents/SFU/1. Ph.D./1. Research Projects/1. Network Feedback/6. Manuscript/Manuscript Versions/MATLAB Scripts/220623_manuscript code/2. Meta-analyzed Data/Meta-analysis/Akt")

ggsave("Akt Analysis - Spline Regression, no exercise.pdf",  Akt_analysis, width = 15/2, height = 10/2, device = "pdf",
       path = "/Users/tjmccoll/Documents/SFU/1. Ph.D./1. Research Projects/1. Network Feedback/6. Manuscript/MATLAB Scripts/15. 210713 manuscript code/2. Meta-analyzed Data/Meta-analysis/Akt")

# Data predicted for the p-Akt time-course (fold change)
predict.change = predict.gam(model, new.data, se.fit=TRUE)
predict.data = data.frame(time = predict.time, change = predict.change$fit, se = err$se.fit[predict.time+1], lci = prd$lci[predict.time+1], uci = prd$uci[predict.time+1]) # +1 because time 0 = row 1 in data.frame
round(predict.data,2)
