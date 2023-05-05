######################
#### Clear memory ####
######################

rm(list = ls())
set.seed(55)

#################################################
#### Set WD and import packages and datasets ####
#################################################
library(ggplot2)
library(R2jags)

data <- read.delim("~/Documents/MEME/S4 Uppsala/MS/data_submission/data_adult_experiment.txt")


######################
#### prepare data ####
######################

data$nestbox <- as.factor(data$nestbox)
data$treatment <- as.factor(data$treatment)
data$graded_resp <- factor(data$graded_resp, ordered=T)

###############
#### model ####
###############

# this model will take about 10-15 minutes to run

# Specify model in BUGS language
sink("graded_resp")
cat("
    model {
    
    # Priors
    mu_treatment[1] <- 0
    mu_treatment[2] ~ dnorm(0, 0.001)
    
    ## priors over thresholds
    for(r in 1:3){
    alpha0[r] ~ dnorm(0, 0.001)
    }
    alpha <- sort(alpha0)
    
    for (i in 1:q) {                                                                
    mu_nest[i] ~ dnorm(0, tau)
    }
    
    sigma ~ dunif(0,10)                                                               
    tau <- 1 / (sigma * sigma)                                                        
    
    # Likelihood
    for (i in 1:n) { 
    lambda[i] <- mu_treatment[treatment[i]] + mu_nest[nestbox[i]]
    logit(Q[i,1]) <- alpha[1]-lambda[i]
    p[i,1] <- Q[i,1]
    for(j in 2:3){
    logit(Q[i,j]) <- alpha[j]-lambda[i]
    p[i,j] <- Q[i,j] - Q[i,j-1]
    }
    p[i,4] <- 1 - Q[i,3] 
    
    y[i] ~ dcat(p[i,])
    
    }
    
    for (i in 1:n) {
    predicted[i] <- lambda[i]
    residual[i] <- y[i]-predicted[i]
    sq[i] <- pow(residual[i], 2)
    
    y.new[i] ~ dcat(p[i,])       
    sq.new[i] <- pow(y.new[i]-predicted[i], 2)  
    }
    
    fit <- sum(sq[])              
    fit.new <- sum(sq.new[])      
    test <- step(fit.new-fit) 		
    bpvalue <- mean(test) 		  	
    
    # derived values
    treat_1_3 <- 1 - (exp(alpha0[3])/(1+exp(alpha0[3])))
    treat_2_3 <- 1 - (exp(alpha0[3]- mu_treatment[2])/(1+exp(alpha0[3]- mu_treatment[2])))
    
    treat_1_0 <- (exp(alpha0[1])/(1+exp(alpha0[1])))
    treat_2_0 <- (exp(alpha0[1]- mu_treatment[2])/(1+exp(alpha0[1]- mu_treatment[2])))
    
    treat_1_1 <- (exp(alpha0[2])/(1+exp(alpha0[2]))) - treat_1_0 
    treat_2_1 <- (exp(alpha0[2]- mu_treatment[2])/(1+exp(alpha0[2]- mu_treatment[2]))) - treat_2_0 
    
    treat_1_2 <- (exp(alpha0[3])/(1+exp(alpha0[3]))) - treat_1_0 - treat_1_1 
    treat_2_2 <- (exp(alpha0[3]- mu_treatment[2])/(1+exp(alpha0[3]- mu_treatment[2]))) - treat_2_0 - treat_2_1 
    
    var_nest <- sigma*sigma
    
    }
    ",fill = TRUE)
sink()

jags.data <- list(y = data$graded_resp, n = dim(data)[1], q = length(unique(data$nestbox)),
                  nestbox = data$nestbox, treatment = data$treatment)

inits <- function(){list(mu_treatment = c(NA, rnorm(1)),
                         alpha0 = c(-0.15, 0.25, 0.55))}  

# Parameters monitored
parameters <- c("mu_treatment", "alpha0", "sigma", "var_nest", "fit", "fit.new", "bpvalue",
                "treat_1_3", "treat_2_3", 
                "treat_1_0", "treat_2_0", 
                "treat_1_1", "treat_2_1", 
                "treat_1_2", "treat_2_2")


# Call JAGS from R
mod_graded_resp <- jags(jags.data, inits, parameters, "graded_resp", n.chains = 3,
                           n.thin = 100, n.iter = 200000, n.burnin = 100000, working.directory = getwd())


print(mod_graded_resp, digits = 3, intervals=c(0.025, 0.975))
#traceplot(mod_graded_resp)

#######################
#### Create figure ####
#######################

mod_alarm2 <- as.mcmc(mod_graded_resp)
int.mcmc.mat <- as.matrix(mod_alarm2)
int.mcmc.dat <- as.data.frame(int.mcmc.mat)   

int.mcmc.dat$treat_1 <- int.mcmc.dat$treat_1_1 + 2*int.mcmc.dat$treat_1_2 + 3*int.mcmc.dat$treat_1_3     # graded response for treatment 1
int.mcmc.dat$treat_2 <- int.mcmc.dat$treat_2_1 + 2*int.mcmc.dat$treat_2_2 + 3*int.mcmc.dat$treat_2_3     # graded response for treatment 2

int.sim <- int.mcmc.dat$treat_1
int.sim2 <- int.mcmc.dat$treat_2

bayes.c.eff.mean <- mean(int.sim)
bayes.c.eff.lower <- quantile(int.sim, probs = c(0.025))
bayes.c.eff.upper <- quantile(int.sim, probs = c(0.975))
bayes.c.eff.lower.bis <- quantile(int.sim, probs = c(0.25))
bayes.c.eff.upper.bis <- quantile(int.sim, probs = c(0.75)) 
plot.dat <- data.frame(bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper, bayes.c.eff.lower.bis, bayes.c.eff.upper.bis)

bayes.c.eff.mean2 <- mean(int.sim2)
bayes.c.eff.lower2 <- quantile(int.sim2, probs = c(0.025))
bayes.c.eff.upper2 <- quantile(int.sim2, probs = c(0.975))
bayes.c.eff.lower2.bis <- quantile(int.sim2, probs = c(0.25))
bayes.c.eff.upper2.bis <- quantile(int.sim2, probs = c(0.75)) 
plot.dat2 <- data.frame(bayes.c.eff.mean2, bayes.c.eff.lower2, bayes.c.eff.upper2, bayes.c.eff.lower2.bis, bayes.c.eff.upper2.bis)
colnames(plot.dat2) <- c("bayes.c.eff.mean", "bayes.c.eff.lower", "bayes.c.eff.upper", "bayes.c.eff.lower.bis", "bayes.c.eff.upper.bis")

plot.dat <- rbind(plot.dat, plot.dat2)

sp <- as.data.frame(matrix(c("Collared", "Pied"), nrow = 2, ncol = 1))
plot.dat <- cbind(plot.dat, sp)

data$graded_resp <- as.numeric(data$graded_resp)-1

p <- ggplot(data, aes(x = treatment, y = graded_resp, colour=treatment))+
  geom_point(position = position_jitter(w = 0.3, h = 0.2),alpha=0.2, size=4)+
  scale_color_manual(values=c("purple", "orange", "seagreen4"))+
  ggtitle("")+
  xlab("")+
  ylab("Graded behavioral response")+
  scale_x_discrete(labels = c("Collared call - Collared song", "Collared call - Pied song"))+
  theme_bw() +
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))
p <- p + annotate("segment", x = 0.7, xend = 1.3, y = plot.dat$bayes.c.eff.mean[1], yend = plot.dat$bayes.c.eff.mean[1], size=2)
p <- p + annotate("segment", x = 1.7, xend = 2.3, y = plot.dat$bayes.c.eff.mean[2], yend = plot.dat$bayes.c.eff.mean[2], size=2)
p <- p + annotate("segment", x = 0.7, xend = 1.3, y = plot.dat$bayes.c.eff.lower[1], yend = plot.dat$bayes.c.eff.lower[1], size=1.1)
p <- p + annotate("segment", x = 1.7, xend = 2.3, y = plot.dat$bayes.c.eff.lower[2], yend = plot.dat$bayes.c.eff.lower[2], size=1.1)
p <- p + annotate("segment", x = 0.7, xend = 1.3, y = plot.dat$bayes.c.eff.upper[1], yend = plot.dat$bayes.c.eff.upper[1], size=1.1)
p <- p + annotate("segment", x = 1.7, xend = 2.3, y = plot.dat$bayes.c.eff.upper[2], yend = plot.dat$bayes.c.eff.upper[2], size=1.1)
p <- p + annotate("segment", x = 1, xend = 1, y = plot.dat$bayes.c.eff.lower[1], yend = plot.dat$bayes.c.eff.lower.bis[1], size=1.1)
p <- p + annotate("segment", x = 1, xend = 1, y = plot.dat$bayes.c.eff.upper.bis[1], yend = plot.dat$bayes.c.eff.upper[1], size=1.1)
p <- p + annotate("segment", x = 2, xend = 2, y = plot.dat$bayes.c.eff.lower[2], yend = plot.dat$bayes.c.eff.lower.bis[2], size=1.1)
p <- p + annotate("segment", x = 2, xend = 2, y = plot.dat$bayes.c.eff.upper.bis[2], yend = plot.dat$bayes.c.eff.upper[2], size=1.1)
p <- p + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower.bis[1], ymax = plot.dat$bayes.c.eff.upper.bis[1], alpha=0.15, color="black", size=1.1, fill="purple")
p <- p + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower.bis[2], ymax = plot.dat$bayes.c.eff.upper.bis[2], alpha=0.15, color="black", size=1.1, fill="orange")
p


############################
#### Plot for the table ####
############################

library(magrittr)
library(dplyr)
library(purrr)
library(forcats)
library(tidyr)
library(modelr)
library(tidybayes)
library(ggplot2)
library(ggstance)
library(ggridges)
library(cowplot)
library(rstan)
library(brms)
library(ggrepel)
library(RColorBrewer)
library(gganimate)

theme_set(theme_tidybayes() + panel_border())

#### fixed effects

vec_post <- c(int.mcmc.dat$`mu_treatment[2]`)
vec_key <- c(rep("Collared call - Pied song", 3000))
data_post <- tibble(vec_key, vec_post)


ggplot(data_post, aes(y = vec_key, x = vec_post)) +
  geom_halfeyeh(scale = "width") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = "gray80")+
  coord_cartesian(xlim = c(-5, 2.3), ylim = c(1.3,1.5))+
  xlab("Effect size")+
  ylab("")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#### random effects

vec_post <- c(int.mcmc.dat$var_nest)
vec_key <- c(rep("Within individual variance", 3000))
data_post <- tibble(vec_key, vec_post)


ggplot(data_post, aes(y = fct_rev(vec_key), x = vec_post)) +
  geom_halfeyeh(scale = "width") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = "gray80")+
  coord_cartesian(xlim = c(0, 20), ylim = c(1.3,1.5))+
  xlab("")+
  ylab("")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#### Intercept ####

vec_post <- c(int.mcmc.dat$`alpha0[1]`, int.mcmc.dat$`alpha0[2]`, int.mcmc.dat$`alpha0[3]`)
vec_key <- c(rep("Intercept 0|1", 3000), rep("Intercept 1|2", 3000), rep("Intercept 2|3", 3000))
data_post <- tibble(vec_key, vec_post)


ggplot(data_post, aes(y = fct_rev(vec_key), x = vec_post)) +
  geom_halfeyeh(scale = "width") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  coord_cartesian(xlim = c(-3, 5), ylim = c(1.3,3.5))+
  xlab("")+
  ylab("")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())



