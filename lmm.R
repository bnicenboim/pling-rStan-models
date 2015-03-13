rm(list=ls())
set.seed(123)

library(MASS)

#############
### Simulation of typical dataset of psycholinguistics
### in latin-square
N_subj <- 80
N_item <- 40
N_coef <- 4  
            
#we assume that the -1000/DV is normally distributed (as typically happens with RT)
#2 conditions with 2 levels: 1a,1b x 2a,2b

(coefs <- c(-1000/500, -1000/540+1000/500, -1000/540+1000/500, -1000/550+1000/500))

sdev_subj <-c(.9,.03,.03,.03) 
sdev_item <- c(.2,.005,.005,.005) 
sdev <- .3                   # sigma in bayesian model
rho_subjs <- .25
rho_items <- .02

#latin square:

conditions <- data.frame(
  c1=rep(c("a","b"),(N_item+1)*N_subj/4),
    c2=rep(c("a","b"),(N_item+1),each=N_item/2),r=rep(1:(N_item+1),N_subj) )

conditions<-conditions[conditions$r!=(N_item+1),c("c1","c2")]

contrasts(conditions$c1) <- contr.sum(2)
contrasts(conditions$c2) <- contr.sum(2)

X <- model.matrix(~1+c1*c2,data=conditions)

# Correlation matrix by subj
rho_subj<- matrix(rep(rho_subjs,N_coef*N_coef),nrow=N_coef)
diag(rho_subj) <-1

# Correlation matrix by item
rho_item<- matrix(rep(rho_items,N_coef*N_coef),nrow=N_coef)
diag(rho_item) <-1




b_subj <- sdev_subj %*% t(sdev_subj)  
Sigma_subj <- b_subj * rho_subj  #variance covariance matrix for subj (Sigma_u in bayesian model)
raneff_subj <- mvrnorm(n = N_subj, rep(0,N_coef), Sigma_subj)


b_item <- sdev_item %*% t(sdev_item)  
Sigma_item <- b_item * rho_item  #variance covariance matrix for item (Sigma_w in bayesian model)
raneff_item <- mvrnorm(n = N_item, rep(0,N_coef), Sigma_item)


#auxiliary matrixes of 1 for the random effects:
subjones <- matrix(rep( c(rep(1,N_item),rep(0,N_item*N_subj)),N_subj),
    nrow=N_subj *N_item,ncol=N_subj )

tempones<-rep(c(rep( c(1,rep(0,N_subj-1)),N_item),0),N_item)
itemones <- matrix(tempones[1:(length(tempones)-N_item)],
    nrow=N_item *N_subj,ncol=N_item )


datalong <- data.frame(subj =factor(rep(seq(1:N_subj),each=N_item)),
                    item =factor(rep(seq(1:N_item),times=N_subj)),
                    trialN=factor(seq(1:N_subj*N_item)),
                    conditions,
                    rt= round(-1000/ rnorm(N_subj*N_item,
                        X %*% matrix(coefs)+
                    rowSums(X %*% t(raneff_subj) * subjones)+
                    rowSums(X %*% t(raneff_item) * itemones)
                    , sdev) ,0)
                    )
head(datalong,10)
library(reshape2)

datalong[datalong$rt <0,]$rt <- 150

dcast(datalong, c1+c2~.,function(x) mean(-1000/x))
dcast(datalong, c1+c2~.,mean)

######################################3
###
### Analysis with lmer
library(lme4)



summary(m1 <- lmer(I(-1000/rt) ~ c1*c2 + (1|subj) +(1|item), datalong))




# 2x2 model with full random-effects structures

summary(m2_full <- lmer(I(-1000/rt) ~ c1*c2 + (c1*c2|subj) +(c1*c2|item), datalong))

#It fails to converge



#################################3
## Bayesian LMM with RStan
library(parallel)

parallelizeStan <- function(fit, data,iter=2000,cores=4){
      print("Starting...")

    temp0 <- sampling(fit,data=data,chains=0)
    print("Starting first chain...")
    sflist <- 
    mclapply(1:cores, mc.cores = cores, 
               function(i) #{sink(sprintf("chain_progress_%d.txt", i))
                stan(fit = temp0, data = data, chains = 1, chain_id = i, iter=iter,
                                refresh = -1)
              #}

               )
    return(sflist2stanfit(sflist))
}


library(rstan)
stan_lmm <- stan_model(file="lmm.stan") 


niter <- 500
nchains <- 4

subj <- as.numeric(as.character(datalong$subj))
N_subj <- length(unique(subj))
item <-as.numeric(as.character(datalong$item))
N_item <- length(unique(item))

rt <- -1000/datalong$rt

# Model equivalent to m1 lmer(rrt ~ c1*c2 + (1|subj) +(1|item), datalong)
x_1 <- model.matrix(~ 1+ c1*c2 , data=datalong)
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


# samples_lmm_1 <- sampling(stan_lmm,    
#                 data=lsdata, 
#                 iter=niter,
#                 chains=nchains
#                 )

samples_lmm_1 <- parallelizeStan(stan_lmm,    
                data=lsdata, 
                iter=niter,
                cores=nchains
                )




#print(samples_lmm_1,pars=c("beta","sigma_u","sigma_w","sigma","Cor_u","Cor_w"))

sum_lmm1 <- summary(samples_lmm_1,pars=c("beta","sigma_u","sigma_w","sigma","Cor_u","Cor_w"),probs = c(0.025,  0.975) ,digits_summary = 3)

#summary in table:
round(sum_lmm1$summary,2)




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




# samples_lmm_full <- sampling(stan_lmm,    
#                 data=lsdata, 
#                 iter=niter,
#                 chains=nchains
#                 )


 samples_lmm_full <- parallelizeStan(stan_lmm,    
                data=lsdata, 
                iter=niter,
                cores=nchains
                )



# print(samples_lmm_full,pars=c("beta","sigma_u","sigma_w","sigma","Cor_u","Cor_w"))

sum_lmm_full<- summary(samples_lmm_full,pars=c("beta","sigma_u","sigma_w","sigma","Cor_u","Cor_w"),probs = c(0.025,  0.975) ,digits_summary = 3)

#summary in table:
round(sum_lmm_full$summary,2)

