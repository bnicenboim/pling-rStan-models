rm(list=ls())
set.seed(123)

library(MASS)

#############
### Simulation of typical dataset of psycholinguistics
### For simplicity of the code not in latin-square
N_subj <- 40
N_item <- 20
N_coef <- 4  #number of predictors

(means <- c(-1000/550, -1000/500,-1000/540, -1000/480))
#-1.818182 -2.000000 -1.851852 -2.083333

# Correlation matrix by subj
rho_subj<- matrix(rep(.2,N_coef*N_coef),nrow=N_coef)
diag(rho_subj) <-1

# Correlation matrix by item
rho_item<- matrix(rep(.02,N_coef*N_coef),nrow=N_coef)
diag(rho_item) <-1


sdev_subj <- c(.7,.05,.05,.05) 
sdev_item <- c(.1,.01,.01,.01) # 
sdev <- .05                    # sigma in bayesian model

b_subj <- sdev_subj %*% t(sdev_subj)  
Sigma_subj <- b_subj * rho_subj  #variance covariance matrix for subj (Sigma_u in bayesian model)
raneff_subj <- mvrnorm(n = N_subj, rep(0,N_coef), Sigma_subj)


b_item <- sdev_item %*% t(sdev_item)  
Sigma_item <- b_item * rho_item  #variance covariance matrix for item (Sigma_w in bayesian model)
raneff_item <- mvrnorm(n = N_item, rep(0,N_coef), Sigma_subj)

datawide <- data.frame(subj =factor(rep(seq(1:N_subj),each=N_item)),
                    item =factor(rep(seq(1:N_item),times=N_subj)),
                    trialN=factor(seq(1:N_subj*N_item)),
                    conditions=
matrix(rep(means,times=N_subj*N_item ),nrow=N_subj*N_item,byrow=T)+
raneff_subj[rep(1:nrow(raneff_subj),each=N_item),]+
raneff_item[rep(1:nrow(raneff_item),times=N_subj),])

head(datawide,100)
library(reshape2)
datalong <- melt(datawide, variable.name="condition",
    value.name="rrt")

head(datalong)
dcast(datalong, condition~.,mean)

######################################3
###
### Analysis with lmer
library(lme4)

# model with 4 levels
contrasts(datalong$condition) <- contr.sum(4)
summary(m1 <- lmer(rrt ~ condition + (1|subj) +(1|item), datalong))

# 2x2 model
datalong$c1 <- factor(ifelse(datalong$condition %in% c('conditions.1','conditions.2'),"A","B"))
datalong$c2 <- factor(ifelse(datalong$condition %in% c('conditions.1','conditions.3'),"X","Y"))

contrasts(datalong$c1) <- contr.sum(2)
contrasts(datalong$c2) <- contr.sum(2)

summary(m2 <- lmer(rrt ~ c1*c2 + (1|subj) +(1|item), datalong))


# 2x2 model with full random-effects structures

#summary(m2_full <- lmer(rrt ~ c1*c2 + (c1*c2|subj) +(c1*c2|item), datalong))
#It didn't converge with the simulated data from above.

#################################3
## Bayesian LMM with RStan

library(rstan)
stan_lmm <- stan_model(file="lmm.stan") 


niter <- 500  #should be more than 2000
nchains <- 4

subj <- as.numeric(as.character(datalong$subj))
N_subj <- length(unique(subj))
item <-as.numeric(as.character(datalong$item))
N_item <- length(unique(item))

rt <-datalong$rrt

# Model equivalent to m1 lmer(rrt ~ condition + (1|subj) +(1|item), datalong)
x_1 <- model.matrix(~ 1+ condition , data=datalong)
x_u_1 <- model.matrix(~ 1, data=datalong) #for subj
x_w_1 <-  model.matrix(~ 1, data=datalong) #for item


lsdata <- list(rt=rt, 
        subj=subj,
        item=item,
        N_obs=nrow(datalong),
        N_coef=ncol(x_1),
        N_coef_u=ncol(x_u_1),
        N_coef_w=ncol(x_w_1),
        x =x_1,
        x_u=x_u_1,
        x_w=x_w_1,
        N_subj=N_subj,
        N_item=N_item
) 


samples_lmm_1 <- sampling(stan_lmm,    
                data=lsdata, 
                iter=niter,
                chains=nchains
                )



print(samples_lmm_1,pars=c("beta","sigma_u","sigma_w","sigma","Cor_u","Cor_w"))




# Model equivalent to m2 lmer(rrt ~ c1*c2 + (1|subj) +(1|item), datalong)
x_2 <- model.matrix(~ 1+ c1*c2 , data=datalong)
x_u_2 <- model.matrix(~ 1, data=datalong) #for subj
x_w_2 <-  model.matrix(~ 1, data=datalong) #for item


lsdata <- list(rt=rt, 
        subj=subj,
        item=item,
        N_obs=nrow(datalong),
        N_coef=ncol(x_2),
        N_coef_u=ncol(x_u_2),
        N_coef_w=ncol(x_w_2),
        x =x_2,
        x_u=x_u_2,
        x_w=x_w_2,
        N_subj=N_subj,
        N_item=N_item
) 




samples_lmm_2 <- sampling(stan_lmm,    
                data=lsdata, 
                iter=niter,
                chains=nchains
                )



print(samples_lmm_2,pars=c("beta","sigma_u","sigma_w","sigma","Cor_u","Cor_w"))



# Model equivalent to full random effects model:
# lmer(rrt ~ c1*c2 + (c1*c2|subj) +(c1*c2|item), datalong)
x_full <- model.matrix(~ 1+ c1*c2 , data=datalong)
x_u_full <- model.matrix(~ 1+c1*c2, data=datalong) #for subj
x_w_full <-  model.matrix(~ 1+c1*c2, data=datalong) #for item


lsdata <- list(rt=rt, 
        subj=subj,
        item=item,
        N_obs=nrow(datalong),
        N_coef=ncol(x_full),
        N_coef_u=ncol(x_u_full),
        N_coef_w=ncol(x_w_full),
        x =x_full,
        x_u=x_u_full,
        x_w=x_w_full,
        N_subj=N_subj,
        N_item=N_item
) 




samples_lmm_full <- sampling(stan_lmm,    
                data=lsdata, 
                iter=niter,
                chains=nchains
                )



print(samples_lmm_full,pars=c("beta","sigma_u","sigma_w","sigma","Cor_u","Cor_w"))
