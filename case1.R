
rm(list = ls())

source("/home/grad/syi/geninf13/source_inf.R")

library(hdi)

n = 100 ; p = 500 ; type = 1 ; s0 = 3 ; J = 1:p ; reptim = 100 ; tau <- 2
m <- length(J) ; thres <- sqrt(tau*log(p))  ; ngrid = 20 ; ngrid2 = 100
B = 200 ; num.cores <- 10
sigsq <- sigsq2 <- sigsq3 <- c()
surr.list <- list()

chosen.tune <- matrix(NA, reptim, m)
mod.chosen.tune <- matrix(NA, reptim, m)

# T1 : our original / T2 : modified / T3 : DB
# T4 : JM  / T5 : our-boots / T6 : DB-boots 

test.stat <- T1 <- T2 <- T3 <- T4 <- T44 <- T5 <- T6 <- B1 <- B2 <- B3 <- B4 <- matrix(NA, reptim, m)
L1 <- L2 <- L3 <- L4 <- L44 <- L5 <- L6 <- matrix(NA, reptim, m)
pval1 <- pval2 <- pval3 <- pval4 <- pval5 <- pval6 <- matrix(NA, reptim, m)

library(doMC) ; registerDoMC(cores = num.cores) 

for(jj in 1:reptim){
  
  itt <- indc <- 0 
  
  while(indc<1){
    
    itt <- itt +1
    set.seed(jj + itt) 
    
    beta0 <- rep(0, p)  
    beta0[1:s0] <- runif(s0, 0, 4)
    
    X <- design(type, n, p) ; Y <- X%*%beta0 + rnorm(n)
    glmnetfit <- cv.glmnet(X, Y) ;  beta.hat <- as.vector(coef(glmnetfit,s=glmnetfit$lambda.1se))[-1]
    err.vec <- (Y-X%*%beta.hat) ; sigsq[jj] <- (sum(err.vec**2))/n
    
    # Use Our tuning-free procedure to get a surrogate set for the strong signals
    vv1 <- Proj.L(X, J, parallel = T)
    for(k in 1:m){
      v <- vv1[k,]
      eta <- Y - X[,-J[k]]%*%beta.hat[-J[k]]
      beta.u <- sum(v*eta)/n
      test.stat[jj,k] <- abs( sqrt(n)*(beta.u-0)/sqrt(sigsq[jj]*sum(v**2)/n))
    }
    surr.list[[jj]] <- est.act.ind <- J[test.stat[jj,] > thres]
    if(length(est.act.ind)>0){
      indc <- 1
    } else{
      indc <- 0
    }
  }
  
  # VBRD's approach
  fit.VBRD <- lasso.proj(X, Y, parallel = T, ncores = num.cores, betainit = beta.hat, sigma = sqrt(sigsq[jj]), return.Z = TRUE, suppress.grouptesting = FALSE)
  T3[jj,] <- abs(((fit.VBRD$bhat-beta0)/fit.VBRD$se)[J])
  L3[jj,] <- (2*qnorm(1-0.05/2)*fit.VBRD$se)[J]
  noise.factor <- c()
  
  for(k in 1:m){
    noise.factor[k] <- (sqrt(sum(fit.VBRD$Z[,J[k]]^2))/abs(sum(fit.VBRD$Z[,J[k]]*X[,J[k]])))
    B3[jj,k] <- abs(((noise.factor[k]^-1)*sum(apply(fit.VBRD$Z[,J[k]]*X[,-J[k]], 2, sum)*(beta0[-J[k]]-beta.hat[-J[k]]))/sum(fit.VBRD$Z[,J[k]]*X[,J[k]]))/sqrt(sigsq[jj]))
  }
  pval3[jj,] = fit.VBRD$pval[J]
  
  # JM's approach
  
  JMfit <- SSLasso(X, Y, alpha=0.05, lambda = glmnetfit$lambda.1se, mu = NULL, intercept = FALSE, 
                                resol=1.3, maxiter=50, threshold=1e-2, verbose = FALSE)
  
  sigsq3[jj] <- (JMfit$noise.sd)^2
  
  for(k in 1:m){
    
    JM.lb <- JMfit$low.lim[J[k]] ; JM.ub <- JMfit$up.lim[J[k]]
    if( JM.lb<=beta0[J[k]] & beta0[J[k]]<= JM.ub){
      T4[jj,k] <- 1
    } else{
      T4[jj,k] <- 0
    }
    L4[jj,k] <- JM.ub-JM.lb
    
    JM.lb2 <- JMfit$low.lim2[J[k]] ; JM.ub2 <- JMfit$up.lim2[J[k]]
    if( JM.lb2<=beta0[J[k]] & beta0[J[k]]<= JM.ub2){
      T44[jj,k] <- 1
    } else{
      T44[jj,k] <- 0
    }
    L44[jj,k] <- JM.ub2-JM.lb2
    
  }

  B4[jj,] <- abs((sqrt(n)*(JMfit$M%*%(t(X)%*%X)/n-diag(1,p))%*%(beta0-JMfit$coef))/(JMfit$noise.sd*sqrt(JMfit$std.errsq)))
  pval4[jj,] <- JMfit$pvals
    
  # For the original procedure
  C.grid <- (seq(0, 15, len = ngrid+1)[-1])/sigsq[jj]
  tune.parms <- expand.grid(C1=C.grid,C2=C.grid) ;  ntune <- dim(tune.parms)[1]
  proj.array <- array(NA, dim = c(ntune, m, n))
  for(i in 1:ntune){
    proj.array[i,,] <- Proj.gen(X, J, tune.parms[i,], est.act.ind, parallel=T)
  }
  
  # refitted OLS for the modified procedure
  
  mod.beta.hat <- rep(0, p) ;  mod.fit <- lm(Y~-1+X[,est.act.ind])
  mod.beta.hat[est.act.ind] <- coef(mod.fit)
  sigsq2[jj] <- sum((Y-X%*%mod.beta.hat)^2)/n
    
  # For the modified procedure

  C.grid2 <- (seq(0, 15, len = ngrid2+1)[-1])/sigsq[jj]
  proj.array2 <- array(NA, dim = c(ngrid2, m, n))
  
  for(i in 1:ngrid2){
    proj.array2[i,,] = Proj.mod(X, J, C.grid2[i], est.act.ind, parallel = T)
  }

  # For Bootstrap
  cent.resid <- err.vec-mean(err.vec)
  boots.cover <- boots.len <- array(NA, dim=c(ntune, B, m))
  boots.cover2 <- boots.len2 <- array(NA, dim=c(ngrid2, B, m))
  
  fun1 <- function(b){
    
    set.seed(b+2017+jj)
    new.Y <- X%*%beta.hat + sample(cent.resid, n, replace = T)
    
    glmnetfit2 <- cv.glmnet(X, new.Y)
    boots.beta.hat <- as.vector(coef(glmnetfit2,s=glmnetfit2$lambda.1se))[-1]
    boots.resid <- new.Y-X%*%boots.beta.hat ; boots.sigmasq <- (sum(boots.resid**2))/n
    
    mod.boots.beta.hat <- rep(0, p) ;  mod.boots.fit <- lm(new.Y~-1+X[,est.act.ind])
    mod.boots.beta.hat[est.act.ind] <- coef(mod.boots.fit)
    mod.boots.sigmasq <- sum((new.Y-X%*%mod.boots.beta.hat)^2)/n
    
    cv1 <- len1 <- matrix(NA, ntune, m) 
    cv2 <- len2 <- matrix(NA, ngrid2, m) 
    
    for(k in 1:m){
      part.resid <- (new.Y - X[,-J[k]]%*%boots.beta.hat[-J[k]])
      for(i in 1:ntune){
        vk <- proj.array[i,k,]
        beta.uk <- sum(vk*part.resid)/n
        Tk <- abs( sqrt(n)*( beta.uk-beta.hat[J[k]] )/sqrt(boots.sigmasq*sum(vk**2)/n) )
        cv1[i,k] <- Tk<qnorm(1-0.05/2)
        len1[i,k] <- 2*qnorm(1-0.05/2)*sqrt(boots.sigmasq*sum(vk**2)/n)/sqrt(n)
      }
      part.resid2 <- (new.Y - X[,-J[k]]%*%mod.boots.beta.hat[-J[k]])
      for(ii in 1:ngrid2){
        vk2 <- proj.array2[ii,k,]
        beta.uk22 <- sum(vk2*part.resid2)/n
        Tk22 <- abs( sqrt(n)*( beta.uk22-beta.hat[J[k]] )/sqrt(mod.boots.sigmasq*sum(vk2**2)/n) )
        cv2[ii,k] <- Tk22<qnorm(1-0.05/2)
        len2[ii,k] <- 2*qnorm(1-0.05/2)*sqrt(mod.boots.sigmasq*sum(vk2**2)/n)/sqrt(n)
      }
    }
    if(b%%50==0) print(paste("Boots",sep=":",b))
    return(list(cov1=cv1,len1=len1,cov2=cv2,len2=len2))
  }
  
  bb2 <- foreach(b=1:B) %dopar% fun1(b)

  for(b in 1:B){
    boots.cover[,b,] <- bb2[[b]]$cov1
    boots.len[,b,] <- bb2[[b]]$len1
    boots.cover2[,b,] <- bb2[[b]]$cov2
    boots.len2[,b,] <- bb2[[b]]$len2
  }
  rm(bb2)
  
  for(k in 1:m){
    cover.k <- apply(boots.cover[,,k], 1, mean) ; cover.k2 <- apply(boots.cover2[,,k], 1, mean)
    len.k <- apply(boots.len[,,k], 1, mean) ; len.k2 <- apply(boots.len2[,,k], 1, mean)
    sd.k <- sqrt(cover.k*(1-cover.k)/B) ; sd.k2 <- sqrt(cover.k2*(1-cover.k2)/B)
    ww1 <- which(cover.k+sd.k-0.95>0) ; ww2 <- which(cover.k2+sd.k2-0.95>0)
    
    if(length(ww1)>0){
      chosen.tune[jj,k] <- ww1[which.min(len.k[ww1])]  
    } else{
      chosen.tune[jj,k] <- which.min(abs(cover.k-0.95))
    }
    
    if(length(ww2)>0){
      mod.chosen.tune[jj,k] <- ww2[which.min(len.k2[ww2])]  
    } else{
      mod.chosen.tune[jj,k] <- which.min(abs(cover.k2-0.95))
    }
  }
  
  fun2 <- function(b){
    set.seed(b)
    new.Y <- X%*%beta.hat + sample(cent.resid, n, replace = T)
    
    glmnetfit2 <- cv.glmnet(X, new.Y)
    boots.beta.hat <- as.vector(coef(glmnetfit2,s=glmnetfit2$lambda.1se))[-1]
    boots.resid <- new.Y-X%*%boots.beta.hat
    boots.sigmasq <- sum(boots.resid^2)/n
    t.star <- t.star.dez <- c()
    
    for(k in 1:m){
      
      debiased <- boots.beta.hat[J[k]] + sum(fit.VBRD$Z[,J[k]]*boots.resid)/sum(fit.VBRD$Z[,J[k]]*X[,J[k]])
      debiased.se <- ((1/n)*sqrt(boots.sigmasq))*sqrt(sum(fit.VBRD$Z[,J[k]]^2))/(abs(sum(fit.VBRD$Z[,J[k]]*X[,J[k]])/n)) 
      t.star.dez[k] <- (debiased - beta.hat[J[k]])/(debiased.se)
      
      vk <- proj.array[chosen.tune[jj,k],k,]
      eta <- (new.Y - X[,-J[k]]%*%boots.beta.hat[-J[k]])
      beta.tilde <- sum(vk*eta)/n
      Tk <- sqrt(n)*( beta.tilde-beta.hat[J[k]] )/sqrt(boots.sigmasq*sum(vk**2)/n) 
      t.star[k] <- Tk
    }
    return(cbind(t.star,t.star.dez))
  }
  
  bb22 <- foreach(b=1:B) %dopar% fun2(b)
  
  
  boots.ours <- boots.dez <- matrix(NA, B, m)
  
  for(b in 1:B){
    boots.ours[b,] <- bb22[[b]][,1]
    boots.dez[b,] <- bb22[[b]][,2]
  }
  rm(bb22)
  
  for(k in 1:m){
    
    v1 <- proj.array[chosen.tune[jj,k],k,]
    ss1 <- sqrt(sum(v1**2)/n)
    eta <- Y - X[,-J[k]]%*%beta.hat[-J[k]]
    beta.u1 <- sum(v1*eta)/n
    
    T1[jj,k] <- abs( sqrt(n)*( beta.u1-beta0[J[k]] )/sqrt(sigsq[jj]*sum(v1**2)/n) )
    L1[jj,k] <- 2*qnorm(1-0.05/2)*sqrt(sigsq[jj]*sum(v1**2)/n)/sqrt(n)
    pval1[jj,k] <- 2*pnorm(abs( sqrt(n)*(beta.u1-0)/sqrt(sigsq[jj]*sum(v1**2)/n)), lower.tail = F)
    B1[jj,k] <- (abs(sum(v1*(X[,-J[k]]%*%(beta0[-J[k]]-beta.hat[-J[k]]))))/sqrt(n))/ss1
    
    lb <- beta.u1-quantile(boots.ours[,k],0.975)*sqrt(sigsq[jj]*sum(v1**2))/n
    ub <- beta.u1-quantile(boots.ours[,k],0.025)*sqrt(sigsq[jj]*sum(v1**2))/n
    if( lb<=beta0[J[k]] & beta0[J[k]]<= ub){
      T5[jj,k] <- 1
    } else{
      T5[jj,k] <- 0
    }
    L5[jj,k] <- (ub-lb)
    pval5[jj, k] <- sum(boots.ours[,k]>T1[jj,k] | boots.ours[,k]<(-T1[jj,k]))/B
    
    v2 <- proj.array2[mod.chosen.tune[jj,k],k,]
    ss2 <- sqrt(sum(v2**2)/n)
    eta2 <- Y - X[,-J[k]]%*%mod.beta.hat[-J[k]]
    beta.u2 <- sum(v2*eta2)/n
    T2[jj,k] <- abs( sqrt(n)*( beta.u2-beta0[J[k]] )/sqrt(sigsq2[jj]*sum(v2**2)/n) )
    L2[jj,k] <- 2*qnorm(1-0.05/2)*sqrt(sigsq2[jj]*sum(v2**2)/n)/sqrt(n)
    pval2[jj,k] <- 2*pnorm(abs( sqrt(n)*(beta.u2-0)/sqrt(sigsq2[jj]*sum(v2**2)/n)), lower.tail = F)
    B2[jj,k] <- (abs(sum(v2*(X[,-J[k]]%*%(beta0[-J[k]]-mod.beta.hat[-J[k]]))))/sqrt(n))/ss2
    
    
    TS = fit.VBRD$bhat[J[k]]/fit.VBRD$se[J[k]]
    lb2 <- fit.VBRD$bhat[J[k]]-quantile(boots.dez[,k],0.975)*fit.VBRD$se[J[k]]
    ub2 <- fit.VBRD$bhat[J[k]]-quantile(boots.dez[,k],0.025)*fit.VBRD$se[J[k]]
    if( lb2<=beta0[J[k]] & beta0[J[k]]<= ub2){
      T6[jj,k] <- 1
    } else{
      T6[jj,k] <- 0
    }
    L6[jj,k] <- (ub2-lb2)
    pval6[jj, k] <- sum(boots.dez[,k]>abs(TS) | boots.dez[,k]<(-abs(TS)))/B
    
  }
  
  print(paste("Iteration",sep=":",jj))
  
}

rm(boots.cover,boots.cover2,boots.len,boots.len2,proj.array,proj.array2)

save.image(paste("/home/grad/syi/geninf13/set1/",sep="",paste("sim",sep="_",type,s0),".RData"))
