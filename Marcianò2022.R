#Beta-regression analysis
#Load required libraries
library(readxl)
library(betareg)
library(emmeans)

#Load .xlsx file containing data
df <- read_excel("Marciano2022.xlsx", sheet = "Mandipropamid")
df <- df[df$Dose != 0, ]

#Transform GIPs (as proportions) for beta-regression (Smithson and Verkuilen, 2006)
SandV2006 <- function(y) { (y*(n-1)+0.5)/n }
n <- nrow(df)
df$Response_GIP <- SandV2006(df$Response_GIP)

#Set factor variables
df$Assay <- as.factor(df$Assay)
df$Dose <- as.factor(df$Dose)

#Set Beta-regression model (BRM)
BRM <- with(df, betareg(Response_GIP ~ Assay + as.factor(Dose), link = "logit"))
summary(BRM)
joint_tests(BRM)
EMM = emmeans(BRM, ~ Dose + Assay)
pairs(EMM, adjust="sidak")
Sum = cld(EMM, alpha   = 0.05, Letters = letters, adjust  = "sidak")
Sum

#Dose-Response (DRA) analysis
#Load required libraries
library(drc)
library(multcomp)

#Load .xlsx file containing data
df <- read_excel("Marciano2022.xlsx", sheet = "Mandipropamid")

#DRA
DRA <- drm(Response_GIP*100 ~ Dose, curveid = Assay, data = df, fct= L.4(fixed = c(NA, 0, NA, NA), names = c("b", "c", "d", "e")))
summary(DRA)
mf <- modelFit(DRA, test = NULL, method = c("cum"))
mf
ED <-ED(DRA, c(50), interval = "delta", type = c("absolute"), multcomp = TRUE)
plot(DRA, broken = FALSE, type = "all", main = "", xlab = "mandipropamid dose (ng/l)", ylab = "Percentage growth inhibition (%)")
axis(1, at = c(0, 0.1, 1, 10, 100, 1000), las=1, labels=c("0","0.1","1","10", "100", "1000"))


#Multiple linear regression (MLR) analysis
#Load required libraries
library(readxl)
library(multcomp)
library(emmeans)

#Load .xlsx file containing data
df <- read_excel("Marciano2022.xlsx", sheet = "Metalaxyl")
df$Dose <- as.factor(df$Dose)

m <- lm(log(OD620) ~ DPI + Dose, data = df)
summary(m)
EMM <- emmeans(m, specs= "Dose", by="DPI")
contrast <- pairs(EMM)
plot(EMM, comparisons = TRUE, ylab = "metalaxyl M dose (mg/l)", xlab = "OD620") + theme_classic(base_size = 20)
marginal = emmeans(m, ~ Dose + DPI)
pairs(marginal, adjust="sidak")
Sum = cld(marginal, alpha   = 0.05, Letters = letters)
Sum

#Kruskal wallis and post-hoc
#Load required libraries
library(agricolae)
library(readxl)

#Load .xlsx file containing data
df <- read_excel("Marciano2022.xlsx", sheet = "Metalaxyl")

KW_1000 <- subset(df, df$Dose == 1000)
PH_1000 <-with(KW_1000,kruskal(OD620,DPI,p.adj="BH",group=TRUE, main=" "))
PH_1000

KW_0 <- subset(df, df$Dose == 0)
PH_0 <-with(KW_0,kruskal(OD620,DPI,p.adj="BH",group=TRUE, main=" "))
PH_0
