########## categorical data analysis chap2

# aspirin 
library(magrittr)
library(vcd) #for threeway, likelihood ratio
library(gmodels) #crosstable
library(fmsb) # for relative risk

############# data
table2.3 = matrix(c(189,10845,104,10933),2,2,byrow = T)
dimnames(table2.3) = list(group=c("Placebo","Aspirin"),MI=c("yes","no"))

### risk and difference
#test for difference of proportions
(table2.3.prop.test = prop.test(table2.3))

#relative risk
fmsb::riskratio(189,104,(189+10845),(104+10933))

#odds.ratio

#1.selfmade function
odds.ratio =function(m,correct=F,print=T){
  if(!is.array(m)) return("must be array or matrix")
  if(length(dim(m))>2) return("must be 2x2 matrix")
  if(dim(m)[1]!=2 | dim(m)[2]!=2) return("must be 2x2 matrix")
  if(correct) m+0.5
  
  theta = (m[1,1]*m[2,2])/(m[1,2]*m[2,1])
  
  log.theta = log(theta)
  ase = (1/m) %>% sum() %>% sqrt()
  log.ci = c(log.theta-qnorm(0.975)*ase,log.theta+qnorm(0.975)*ase)
  ci = exp(log.ci)
  
  ##prepare result
  if(dimnames(m) %>% length()>0){
    categ = sprintf("(%s x %s)",names(dimnames(m)[1]),names(dimnames(m)[2]))
    names(theta) = sprintf("odds.ratio %s",categ)
  }else{
    names(theta) = sprintf("odds.ratio")    
  } 
  names(ci) = c("lower","upper")
  if(print) print(theta)
  invisible(list(statistic=theta,confidence.interval=ci))
}

#2.vcd
vcd::oddsratio(table2.3)

##2.4 chi square test
#godness of fit 
godness.example =  rmultinom(1, 300, p<-c(0.2,0.3,0.1,0.4)) %>% as.numeric() 
chisq.test(godness.example,p=p)
#chi square test
table2.3.chi = chisq.test(table2.3,cor=F)

##SAS-like Output
table2.3.chi2 = gmodels::CrossTable(m,prop.c = T,prop.r = T,expected = T,prop.chisq = F,prop.t = T,
                             resid=T,sresid = T,chisq=T,format = "SAS") 

## likelihood test
llh.test = function(m){
  total = m %>% sum()
  mu = replicate(dim(m)[2],rowSums(m)) * t(replicate(dim(m)[1],colSums(m)))/total
  g.square = (m*log(m/mu)) %>% sum()*2  
  return(g.square)
}

llh.test(table2.3)
##vcd
vcd::assocstats(table2.3)



# example : Gender Gap in Political Affiliation

table2.5 = matrix(c(279,73,225,165,47,191),2,3,byrow = T)
dimnames(table2.5) = list(gender=c("F","M"),partyIden=c("demo","inde","repub"))
table2.5.chi = chisq.test(table2.5)

#standardized pearson residual:
table2.3.chi$stdres

######## three-way table ############
table2.10 = c(53,11,414,37,0,4,16,139)
dim(table2.10) = c(2,2,2)
dimnames(table2.10) = list(defendant=c("W","B"),death.penalty=c("Y","N"),victim=c("W","B"))
#view 3way table
vcd::structable(death.penalty~victim+defendant,data=table2.10) 



threeway.odds.ratio = function(m,print=T){
  if(!is.array(m)) return("must be three dim array")
  #len = length(dim(m))
  if(length(dim(m))==2) return(odds.ratio(m))
  if(sum(dim(m))!=6) return("must be 2x2x2 matrix")
  
  #conditional odds.ratio for 3way table
  if(print) cat("\n##conditional odds ratio:\n")
  con =  sapply(3:1,function(x){apply(m,x,odds.ratio,print=print)})
  con.vec= c()
  for(i in 1:6){
    element = con[[i]]$statistic
    names(element) = strsplit(names(con[[i]]$statistic),split = "odds.ratio ")[[1]][2]
    con.vec = c(con.vec,element)
  }
  #marginal
  if(print) cat("\n##marginal odds ratio:\n")
  mar = sapply(3:1,function(x){(margin.table(m,-x) %>% odds.ratio(.,print=print))$statistic})
  
  res = list(conditional = con.vec, marginal = mar)
  invisible(res)
}

########## Example: con.independent vs mar.independent
table.2.11 = c(18,12,12,8,2,8,8,32) 
dim(table.2.11)=c(2,2,2)
dimnames(table.2.11) = list(Treatment=c("A","B"),Response=c("Success","Failure"),Clinic=c(1,2))
vcd::structable(Response~Clinic+Treatment,table.2.11)
threeway.odds.ratio(table.2.11)





########### 4.3.4 COCHRAN-MANTEL-HAENSZEL METHODS
tab.mantel = c(126,35,100,61,
               908,497,688,807,
               913,336,747,598,
               235,58,172,121,
               402,121,308,215,
               182,72,156,98,
               60,11,99,43,
               104,21,89,36)
dim(tab.mantel) = c(2,2,8)
dimnames(tab.mantel) = list(smoking=c("Y","N"),lung.cancer = c("Y","N"),
                            city=c("Beijing","Shanghai","Shenyang","Nanjing","Harbin","Zhengzhou","Taiyuan","Nanchang"))

structable(lung.cancer~city+smoking,data=tab.mantel)

tab.mantel.test =  mantelhaen.test(tab.mantel,correct = F)


### common odds ratio
numer = sapply(1:8,function(x){
  sum = tab.mantel[,,x] %>% sum()
  product = tab.mantel[1,1,x]*tab.mantel[2,2,x]/sum
}) %>% sum()
denom = sapply(1:8,function(x){
  sum = tab.mantel[,,x] %>% sum()
  product = tab.mantel[1,2,x]*tab.mantel[2,1,x]/sum
}) %>% sum()
numer/denom

tab.mantel.test$estimate

##Breslow-Day statistic:
# source : https://onlinecourses.science.psu.edu/stat504/node/114
breslowday.test <- function(x) {
  #Find the common OR based on Mantel-Haenszel
  or.hat.mh <- mantelhaen.test(x)$estimate
  #Number of strata
  K <- dim(x)[3]
  #Value of the Statistic
  X2.HBD <- 0
  #Value of aj, tildeaj and Var.aj
  a <- tildea <- Var.a <- numeric(K)
  
  for (j in 1:K) {
    #Find marginals of table j
    mj <- apply(x[,,j], MARGIN=1, sum)
    nj <- apply(x[,,j], MARGIN=2, sum)
    
    #Solve for tilde(a)_j
    coef <- c(-mj[1]*nj[1] * or.hat.mh, nj[2]-mj[1]+or.hat.mh*(nj[1]+mj[1]),
              1-or.hat.mh)
    sols <- Re(polyroot(coef))
    #Take the root, which fulfills 0 < tilde(a)_j <= min(n1_j, m1_j)
    tildeaj <- sols[(0 < sols) &  (sols <= min(nj[1],mj[1]))]
    #Observed value
    aj <- x[1,1,j]
    
    #Determine other expected cell entries
    tildebj <- mj[1] - tildeaj
    tildecj <- nj[1] - tildeaj
    tildedj <- mj[2] - tildecj
    
    #Compute \hat{\Var}(a_j | \widehat{\OR}_MH)
    Var.aj <- (1/tildeaj + 1/tildebj + 1/tildecj + 1/tildedj)^(-1)
    
    #Compute contribution
    X2.HBD <- X2.HBD + as.numeric((aj - tildeaj)^2 / Var.aj)
    
    #Assign found value for later computations
    a[j] <- aj ;  tildea[j] <- tildeaj ; Var.a[j] <- Var.aj
  }
  
  #Compute Tarone corrected test
  X2.HBDT <-as.numeric( X2.HBD -  (sum(a) - sum(tildea))^2/sum(Var.aj) )
  
  #Compute p-value based on the Tarone corrected test
  p <- 1-pchisq(X2.HBDT, df=K-1)
  
  res <- list(X2.HBD=X2.HBD,X2.HBDT=X2.HBDT,p=p)
  class(res) <- "bdtest"
  return(res)
}

print.bdtest <- function(x) {
  cat("Breslow and Day test (with Tarone correction):\n")
  cat("Breslow-Day X-squared         =",x$X2.HBD,"\n")
  cat("Breslow-Day-Tarone X-squared  =",x$X2.HBDT,"\n\n")
  cat("Test for test of a common OR: p-value = ",x$p,"\n\n")
}


