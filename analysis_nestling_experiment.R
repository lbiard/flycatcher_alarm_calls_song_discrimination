######################
#### Clear memory ####
######################

rm(list = ls())

#################################################
#### Set WD and import packages and datasets ####
#################################################
library(R2jags)
library(ggplot2)

data <- read.delim("~/data_nestling_experiment.txt")

#### check zero inflation index ####
1 + log(length(data$Beg_play[data$Beg_play==0])/length(data$Beg_play))/mean(data$Beg_play)  #more than 0 so it is ZI

#### check overdispersion ####
var(data$Beg_play)/mean(data$Beg_play) # it is overdispersed


######################
#### prepare data ####
######################

nest <- as.numeric(as.factor(data$Box))
chick <- as.numeric(as.factor(data$Chick_ID))

y <- as.numeric(data$Beg_play)
rest <- as.numeric(scale(data$Beg_rest))
condition <- as.numeric(scale(data$Mass))
treatment <- as.factor(data$Treatment)
playback <- as.factor(data$Playback_ID)

data_re <- data.frame(cbind(chick, nest))
data_re <- data_re[!duplicated(data_re),]
chickbis <- data_re$chick
nestbis <- data_re$nest


###############
#### model ####
###############

# this model will take about 2h-3h to run

# Specify model in BUGS language
sink("zip.jags")
cat("
    model {
    
    # Priors
    psi ~ dunif(0,1)                                                                
    alpha ~ dnorm(0,0.001)
    beta_rest ~ dnorm(0,0.001)
    beta_cond ~ dnorm(0,0.001)
    mu_treatment[1] <- 0
    mu_treatment[2] ~ dnorm(0, 0.001)
    z ~ dunif(0,50)
    
    for (i in 1:q) {                                                                
    mu_nest[i] ~ dnorm(0, tau)
    }
    
    for (i in 1:r) {                                                                
    mu_chick[i] ~ dnorm(mu_nest[nest_vector[i]], tau_bis)
    }
    
    for (i in 1:s) {                                                               
    mu_pb[i] ~ dnorm(0, tau_ter)
    }    
    
    sigma ~ dunif(0,3)                                                             
    tau <- 1 / (sigma * sigma)                                                     
    var_nest <- (sigma * sigma)
    
    sigma_bis ~ dunif(0,3)                                                         
    tau_bis <- 1 / (sigma_bis * sigma_bis)                                         
    var_id <- (sigma_bis * sigma_bis)
    
    sigma_ter ~ dunif(0,3)                                                        
    tau_ter <- 1 / (sigma_ter * sigma_ter)
    var_pb <- (sigma_ter * sigma_ter)
    
    
    # Likelihood Zero Inflated negative binomial
    for (i in 1:n) {
    C[i] ~ dnegbin(p[i], z)
    p[i] <- z/(z+eff.lambda[i])
    eff.lambda[i] <- W[i]*lambda[i] + 0.00001  
    W[i] ~ dbern(psi) 
    log(lambda[i]) <- alpha + beta_rest*rest[i] + beta_cond*condition[i] + mu_treatment[treatment[i]] + mu_chick[chick_vector[i]] + mu_pb[playback[i]]
    }
    
    for (i in 1:n) {
    predicted[i] <- eff.lambda[i]
    residual[i] <- C[i]-predicted[i]
    sq[i] <- pow(residual[i], 2)
    
    C.new[i]~dnegbin(p[i], z)      
    sq.new[i] <- pow(C.new[i]-predicted[i], 2)  
    }
    
    fit <- sum(sq[])              
    fit.new <- sum(sq.new[])     
    test <- step(fit.new-fit) 	
    bpvalue <- mean(test) 		  	
    
    }
    ",fill = TRUE)
sink()

jags.data <- list(C = y, n = dim(data)[1], q = length(unique(nestbis)), r = length(unique(chick)), s = length(unique(playback)),
                  nest_vector = nestbis, chick_vector = chick, rest = rest, condition = condition,
                  treatment = treatment, playback = playback)

inits <- function(){list(psi = runif(0,1), beta_rest = rnorm(1), beta_cond = rnorm(1), alpha = rnorm(1),
                         mu_treatment = c(NA, rnorm(1)))}  

# Parameters monitored
parameters <- c("beta_rest", "beta_cond", "mu_treatment", "z", "psi", "alpha", "var_id", "var_pb", "var_nest", "fit", "fit.new", "bpvalue")


# Call JAGS from R
zip.trial <- jags(jags.data, inits, parameters, "zip.jags", n.chains = 3,
                           n.thin = 250, n.iter = 500000, n.burnin = 250000, working.directory = getwd())

options(scipen=999)
print(zip.trial, digits = 3, intervals=c(0.025, 0.975))
traceplot(zip.trial)


#######################
#### Create figure ####
#######################

data_no_zero <- subset(data, Beg_play>0)
data_zero <- subset(data, Beg_play==0)

p <- ggplot(data, aes(x = Treatment, y = log(Beg_play), colour=Treatment))+
  annotate("segment", x = 0.7, xend = 1.3, y = log(mean(data$Beg_play[data$Treatment=="Collared"])), 
           yend = log(mean(data$Beg_play[data$Treatment=="Collared"])), size=2)+
  annotate("segment", x = 1.7, xend = 2.3, y = log(mean(data$Beg_play[data$Treatment=="Pied"])), 
           yend = log(mean(data$Beg_play[data$Treatment=="Pied"])), size=2)+
  ggtitle("")+
  xlab("")+
  ylab("Log begging calls")+
  scale_x_discrete(labels = c("Collared call - Collared song", "Collared call - Pied song"))+
  coord_cartesian(ylim = c(-0.5, 4.4))+
  theme_bw() +
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))
p <- p + geom_point(data = data_no_zero, aes(x = Treatment, y = log(Beg_play), colour=Treatment), inherit.aes = F,
                    position = position_jitter(w = 0.3, h = 0), alpha = 0.2, size=4)+
  scale_color_manual(values=c("purple", "orange"))
p <- p + geom_point(data = data_zero, aes(x = Treatment, y = log(Beg_play), colour=Treatment), inherit.aes = F,
                    position = position_jitter(w = 0.3, h = 0), alpha = 0.05, size=4)
p



############################
#### Plot for the table ####
############################

zip.trial2 <- as.mcmc(zip.trial)
int.mcmc.mat <- as.matrix(zip.trial2)
int.mcmc.dat <- as.data.frame(int.mcmc.mat)

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

vec_post <- c(int.mcmc.dat$beta_cond, int.mcmc.dat$beta_rest, int.mcmc.dat$`mu_treatment[2]`)
vec_key <- c(rep("Condition", 3000), rep("Rest", 3000), rep("Collared call - Pied song", 3000))
data_post <- tibble(vec_key, vec_post)

ggplot(data_post, aes(y = vec_key, x = vec_post)) +
  geom_halfeyeh(scale = "width") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = "gray80")+
  coord_cartesian(xlim = c(-1.2, 2.1), ylim = c(1.3,3.5))+
  xlab("Effect size")+
  ylab("")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#### random effects

vec_post <- c(int.mcmc.dat$var_id, int.mcmc.dat$var_nest, int.mcmc.dat$var_pb)
vec_key <- c(rep("Within nestling variance", 3000), rep("Within nest variance", 3000), rep("Within playback variance", 3000))
data_post <- tibble(vec_key, vec_post)


ggplot(data_post, aes(y = fct_rev(vec_key), x = vec_post)) +
  geom_halfeyeh(scale = "width") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = "gray80")+
  coord_cartesian(xlim = c(0, 2.5), ylim = c(1.3,3.5))+
  xlab("")+
  ylab("")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#### Intercept ####

vec_post <- c(int.mcmc.dat$alpha)
vec_key <- c(rep("Intercept", 3000))
data_post <- tibble(vec_key, vec_post)


ggplot(data_post, aes(y = fct_rev(vec_key), x = vec_post)) +
  geom_halfeyeh(scale = "width") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  coord_cartesian(xlim = c(-1.5, 3.5), ylim = c(1.3,1.5))+
  xlab("")+
  ylab("")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

