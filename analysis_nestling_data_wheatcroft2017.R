######################
#### Clear memory ####
######################

rm(list = ls())
set.seed(5)
#################################################
#### Set WD and import packages and datasets ####
#################################################

library(R2jags)
library(ggplot2)

data <- read.delim("~/data_nestling_wheatcroft2017.txt")



#### check zero inflation index ####
1 + log(length(data$beg.play[data$beg.play==0])/length(data$beg.play))/mean(data$beg.play)  #more than 0 so it is ZI

#### check overdispersion ####
var(data$beg.play)/mean(data$beg.play) # it is overdispersed


######################
#### prepare data ####
######################

nest <- as.numeric(as.factor(data$nestbox))

y <- as.numeric(data$beg.play)
rest <- as.numeric(scale(data$beg.rest))
condition <- as.numeric(scale(data$mass))
prop <- as.numeric(scale(data$Call.prop))
playback <- as.factor(data$recording)
playback <- droplevels(playback)


###############
#### model ####
###############

# Specify model in BUGS language
sink("zip_prop.jags")
cat("
    model {
    
    # Priors
    psi ~ dunif(0,1)                                                               
    alpha ~ dnorm(0,0.001)
    beta_rest ~ dnorm(0,0.001)
    beta_cond ~ dnorm(0,0.001)
    beta_prop ~ dnorm(0,0.001)
    z ~ dunif(0,50)    
    
    
    for (i in 1:q) {                                                                
    mu_nest[i] ~ dnorm(0, tau)
    }
    
    for (i in 1:s) {                                                                
    mu_pb[i] ~ dnorm(0, tau_bis)
    }    
    
    sigma_bis ~ dunif(0,5)
    tau_bis <- 1/(sigma_bis * sigma_bis)
    var_pb <- sigma_bis * sigma_bis    
    
    sigma ~ dunif(0,5)                                                             
    tau <- 1 / (sigma * sigma)                                                     
    var_nest <- sigma*sigma
    
    
    # Likelihood Zero Inflated negative binomial
    for (i in 1:n) {
    C[i] ~ dnegbin(p[i], z)
    p[i] <- z/(z+eff.lambda[i])
    eff.lambda[i] <- W[i]*lambda[i] + 0.00001  
    W[i] ~ dbern(psi)   
    log(lambda[i]) <- alpha + beta_rest*rest[i] + beta_cond*condition[i] + beta_prop*prop[i] + mu_nest[nest_vector[i]] + mu_pb[playback[i]]          
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
    bpvalue <- mean(test) 		  	# Bayesian p-value
    
    }
    ",fill = TRUE)
sink()

jags.data <- list(C = y, n = dim(data)[1], q = length(unique(nest)),
                  nest_vector = nest, rest = rest, condition = condition, prop = prop, 
                  s = length(unique(playback)), playback = playback)

inits <- function(){list(psi = runif(0,1), beta_rest = rnorm(1), beta_cond = rnorm(1), alpha = rnorm(1))}  

# Parameters monitored
parameters <- c("beta_rest", "beta_cond", "beta_prop", "psi", 
                "alpha", "sigma", "var_nest", "sigma_bis", "var_pb",
                "fit", "fit.new", "bpvalue")


# Call JAGS from R
zip_prop <- jags(jags.data, inits, parameters, "zip_prop.jags", n.chains = 3,
                           n.thin = 250, n.iter = 500000, n.burnin = 250000, working.directory = getwd())

print(zip_prop, digits = 3, intervals=c(0.025, 0.975))
traceplot(zip_prop)


############################
#### Plot for the table ####
############################

model_prop <- as.mcmc(zip_prop)
int.mcmc.mat <- as.matrix(model_prop)
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

vec_post <- c(int.mcmc.dat$beta_cond, int.mcmc.dat$beta_rest, int.mcmc.dat$beta_prop)
vec_key <- c(rep("Condition", 3000), rep("Rest", 3000), rep("Aproportion", 3000))
data_post <- tibble(vec_key, vec_post)


ggplot(data_post, aes(y = vec_key, x = vec_post)) +
  geom_halfeyeh(scale = "width") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = "gray80")+
  coord_cartesian(xlim = c(-2.8, 3), ylim = c(1.3,3.5))+
  xlab("Effect size")+
  ylab("")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#### random effects

vec_post <- c(int.mcmc.dat$var_nest, int.mcmc.dat$var_pb)
vec_key <- c(rep("Within nest variance", 3000), rep("Within playback variance", 3000))
data_post <- tibble(vec_key, vec_post)

ggplot(data_post, aes(y = fct_rev(vec_key), x = vec_post)) +
  geom_halfeyeh(scale = "width") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = "gray80")+
  coord_cartesian(xlim = c(0, 17), ylim = c(1.3,2.5))+
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
  coord_cartesian(xlim = c(-3.2, 1.2), ylim = c(1.3,1.5))+
  xlab("")+
  ylab("")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())



data$new <- NA
data$new[data$Call.prop < 0.5] <- "a_p<0.5"
data$new[data$Call.prop >= 0.5] <- "b_p>0.5"
data$new <- as.factor(data$new)
data_no_zero <- subset(data, beg.play>0)
data_zero <- subset(data, beg.play==0)

p <- ggplot(data, aes(x=new, y=log(beg.play)))+
  geom_point(alpha=0)+
  annotate("segment", x = 0.7, xend = 1.3, y = log(mean(data$beg.play[data$new=="a_p<0.5"])), 
           yend = log(mean(data$beg.play[data$new=="a_p<0.5"])), size=2)+
  annotate("segment", x = 1.7, xend = 2.3, y = log(mean(data$beg.play[data$new=="b_p>0.5"])), 
           yend = log(mean(data$beg.play[data$new=="b_p>0.5"])), size=2)+
  ggtitle("")+
  scale_x_discrete(labels = c("P<0.5", "P>0.5"))+
  xlab("Proportion of song phrases starting with a call")+
  ylab("Log begging calls")+
  theme_bw() +
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))
p <- p + geom_point(data = data_no_zero, aes(x=new, y=log(beg.play)), inherit.aes = F,
                    position = position_jitter(w = 0.3, h = 0), alpha = 0.2, size=4)
p <- p + geom_point(data = data_zero, aes(x=new, y=log(beg.play)), inherit.aes = F,
                    position = position_jitter(w = 0.3, h = 0), alpha = 0.05, size=4)
p
