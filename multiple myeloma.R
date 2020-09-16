# libraries -------------------------------------------------------------
library(rjags)
library(stringr)
library(netmeta)
library(taRifx)
library(ggplot2)
library(xlsx)
library(MCMCvis)
library(igraph)
options(digits=2)
rm(list=ls()) 

# load auxiliary functions
source("funcs.R")    


# set up data -------------------------------------------------------------
## Data by Susanne Schmitz (BMC Med Res Meth 2018, 18:66)
## Hazard ratios such that small values are good


## White Network
##
alive <- rbind(
  c(63, 186), c(161, 162), c(26, 24), c(160, 158), 
  c(168, 166), c(6, 8), c(88, 88), c(88, 88), 
  c(198, 198), c(56, 54), c(76, 151), c(57, 55), 
  c(162, 160), c(182, 180), c(29, 28), c(142, 143))
n <- rbind(
  c(126, 373), c(322, 324), c(53, 49), c(320, 317), 
  c(336, 333), c(12, 16), c(175, 176), c(176, 177), 
  c(396, 396), c(113, 108), c(153, 302), c(114, 110), 
  c(325, 321), c(363, 360), c(58, 57), c(283, 286))
treat.nr <- rbind(
  c(1, 4), c(2, 6), c(2, 7), c(2, 8),
  c(1, 2), c(4, 9), c(1, 3), c(1, 3),
  c(3, 10), c(5, 11), c(1, 5), c(1, 12),
  c(3, 13), c(3, 14), c(3, 14), c(3, 15))
time <- rbind(
  c(6.00, 7.40), c(6.50, 9.00), c(5.10, 6.20), c(6.83, 7.63),
  c(3.49, 6.22), c(7.90, 1.50), c(4.70, 11.30), c(4.70, 11.10),
  c(17.60, 26.30), c(4.20, 2.70), c(1.90, 4.00), c(3.55, 3.09),
  c(14.90, 19.40), c(14.70, 20.60), c(4.00, 6.70), c(18.40, 54.10))
treats.white <- c("dex", "bor", "dex + len", "thal", "pom + dex", "bor + PLD",
                  "bor + bev", "bor + vor", "thal + IFN", "carf + dex + len",
                  "pom", "ob + dex", "elo + dex + len", "ixa + dex + len",
                  "dara + dex + len")
treat <- matrix(treats.white[treat.nr], nrow = 16, ncol = 2)
##
white <- data.frame(alive, n, treat.nr, treat, time, source = "white")
colnames(white) <- c("alive1", "alive2", "n1", "n2",
                     "treat1.nr", "treat2.nr", "treat1", "treat2",
                     "time1", "time2", "source")


## Calculate hazard ratios and their standard errors from the given
## data under the exponential distribution assumption:
## - alive / n = S(t) = exp(-lambda*t)
## - lambda = (log(n) - log(alive)) / t
##
white$lambda1 <- (log(white$n1) - log(white$alive1)) / white$time1
white$lambda2 <- (log(white$n2) - log(white$alive2)) / white$time2
white$HR <- white$lambda1 / white$lambda2
white$logHR <- log(white$HR)


## Derivation of variances see text ("PoissonExponentialEstimation.pdf")
##
white$var1 <-
  (1 / white$alive1 - 1 / white$n1) / (white$lambda1 * white$time1)^2
white$var2 <-
  (1 / white$alive2 - 1 / white$n2) / (white$lambda2 * white$time2)^2
white$se.logHR <- sqrt(white$var1 + white$var2)


## Black Network
##
alive <- rbind(
  c(32, 34), c(72, 71), c(33, 34), 
  c(190, 194), c(232, 232), c(22, 24), 
  c(38, 38), c(67, 68), c(124, 126)) 
n <- rbind(
  c(64, 67), c(144, 142), c(66, 69), 
  c(381, 387), c(465, 464), c(43, 47), 
  c(75, 77), c(134, 135), c(247, 251)) 
treat.nr <- rbind(
  c(1, 2), c(1, 3), c(1, 4), 
  c(1, 5), c(1, 6), c(1, 7), 
  c(1, 8), c(2, 9), c(1, 10)) 
time <- rbind(
  c(7.20, 9.00), c(7.60, 8.10), c(8.29, 5.23), 
  c(8.10, 12.00), c(9.40, 18.70), c(12.60, 9.90), 
  c(6.90, 9.70), c(13.60, 18.30), c(7.20, 18.50))
treats.black <- c("bor + dex", "thal + dex", "bor + dex + sil",
                  "bor + dex + peri", "bor + dex + pan", "carf + dex",
                  "bor + dex + cyc", "elo + bor + dex", "thal + bor + dex",
                  "dara + bor + dex")
treat <- matrix(treats.black[treat.nr], nrow = 9, ncol = 2)
##
black <- data.frame(alive, n, treat.nr, treat, time, source = "black")
colnames(black) <- c("alive1", "alive2", "n1", "n2",
                     "treat1.nr", "treat2.nr", "treat1", "treat2",
                     "time1", "time2", "source")


## Calculate hazard ratios and their standard errors from the given data
## under the exponential distribution assumption
##
black$lambda1 <- (log(black$n1) - log(black$alive1)) / black$time1
black$lambda2 <- (log(black$n2) - log(black$alive2)) / black$time2
black$HR <- black$lambda1 / black$lambda2
black$logHR <- log(black$HR)


## Derivation of variances see Text ("PoissonExponentialEstimation.pdf")
##
black$var1 <- (1 / black$alive1 -
                 1 / black$n1) / (black$lambda1 * black$time1)^2
black$var2 <- (1 / black$alive2 -
                 1 / black$n2) / (black$lambda2 * black$time2)^2
black$se.logHR <- sqrt(black$var1 + black$var2)


## Both networks together
##
alldata <- rbind(white, black)
alldata$id <- seq_len(nrow(alldata))
alldata2=alldata[, c("id", "treat1", "treat2", "logHR", "se.logHR")]

dput(unique(c(alldata2$treat1, alldata2$treat2)))
dput(unique(c("dex","bor","thal","dex","len","pom","dex","bor","dex",
"thal","dex","bor","PLD","bor","bev","bor","vor","thal","IFN",
"carf","dex","len","pom","ob","dex","elo","dex","len","ixa","dex","len",
"dara","dex","len","bor","dex","sil","bor","dex","peri","bor","dex","pan",
"carf","dex","bor","dex","cyc","elo","bor","dex","thal","bor","dex",
"dara","bor","dex")))
components=c("dex", "bor", "thal", "len", "pom", "PLD", "bev", "vor", "IFN", 
             "carf", "ob", "elo", "ixa", "dara", "sil", "peri", "pan", "cyc"
)
Nc=length(components)
for (i in 1:Nc){
  eval(parse(text=paste("alldata2$", components[i], "=0")))
}

for (i in 1:Nc){
  for (j in 1:length(alldata2$id)){
   alldata2[j, which(colnames(alldata2)==components[i])]=
     1* grepl(components[i], alldata2$treat1[j], fixed = TRUE)-1*grepl(components[i], alldata2$treat2[j], fixed = TRUE)
  }
}


for (i in 1:Nc){eval(parse(text=paste("c", i, "=alldata2$", components[i], sep="")))}
c=alldata2[, c(which(colnames(alldata2)==components[1]):which(colnames(alldata2)==components[Nc]))]
ns=length(alldata2$id)

#  with discomb -----------------------------------------------------------
net.additive <- discomb(logHR, se.logHR, treat1, treat2, id,
                        data = alldata,
                        sm = "MD", 
                        comb.fixed = FALSE)

observed.inter=interactions(net.additive)



# create plots ------------------------------------------------------------

## plot at the component level
tt <- graph.data.frame(rbind(alldata[,c(7,8)]), directed = F)
plot(tt, vertex.size=5,vertex.label.cex=0.7, edge.color="gray60",vertex.color="gray85",
     vertex.label.font=3, vertex.label.color="black", vertex.label.dist=2,layout=layout.circle(tt))

## plot the White and black subnetworks
white$id <- seq_len(nrow(white))
black$id <- seq_len(nrow(black))

## 
net.white <- netmeta(logHR, se.logHR, treat1, treat2, id,
                     data = alldata, subset = source == "white",
                     sm = "HR", comb.fixed = FALSE)
netgraph(net.white, seq = "optimal", col = "black", plastic = FALSE,
         points = TRUE, pch = 21, cex.points = 3, col.points = "black",
         bg.points = "gray", thickness ="equal",lwd=2,
         multiarm = FALSE, number.of.studies = TRUE)
##
net.black <- netmeta(logHR, se.logHR, treat1, treat2, id,
                     data = alldata, subset = source == "black",
                     sm = "HR", comb.fixed = FALSE)

netgraph(net.black, seq = "optimal", col = "black", plastic = FALSE,
         points = TRUE, pch = 21, cex.points = 3, col.points = "black",
         bg.points = "gray", thickness ="equal",lwd=2,
         multiarm = FALSE, number.of.studies = TRUE)


# no interactions model ---------------------------------------------------  -------------------------------------------------------
model1.string <-  "
model {
for(i in 1:Ns) { 
prec[i]<-1/(se[i]*se[i])
  y[i]~dnorm(phi[i],prec[i])
}

for(i in 1:Ns) { 
phi[i]~dnorm(mean[i], prec.tau)
mean[i]<- inprod(d[], c[i,])
}

##prior distribution for heterogeneity	
tau ~ dnorm(0,0.1)I(0,)                                      
prec.tau<- 1/tau.sq
tau.sq<- pow(tau,2)
for(i in 1:Nc){d[i]~dnorm(0,0.01)}

example1<- d[1]+d[3]-d[1]-d[2]-d[3] ## thal + dex VS thal + bor + dex
example2<- d[1]+d[2]+d[16]-d[3]-d[9] ## bor + dex + peri VS thal + IFN
example3<- d[1]+d[11]-d[14]-d[1]-d[4] ## ob + dex VS dara + dex + len
example4<- d[1]+d[2]+d[18]-d[3]-d[9]  ## bor + dex + cyc VS thal + IFN
}
"


model1.spec<-textConnection(model1.string) 
data <- list(y=alldata2$logHR,se=alldata2$se.logHR, Nc=Nc, Ns=ns,c=c)
jags.m=0
jags.m <- jags.model(model1.spec, data = data, n.chains =4, n.adapt = 10000)

params <- c("tau", "d", "example1" , "example2" , "example3", "example4") 
closeAllConnections()
samps<- coda.samples(jags.m, params, n.iter =30000)
MCMCtrace(samps,pdf = FALSE, params = "d") 
A1= MCMCsummary(samps)
A1["HR.1",]=exp(A1["example1",])
A1["HR.2",]=exp(A1["example2",])
A1["HR.3",]=exp(A1["example3",])
A1["HR.4",]=exp(A1["example4",])
rownames(A1)[1:Nc]=components
nointer=round(A1[,c(4,3,5)], digits=3)
A1$results=paste(format(round(A1$`50%`,2), nsmall=2), "[", 
                 format(round(A1$`2.5%`,2), nsmall=2), ";", 
                 format(round(A1$`97.5%`,2), nsmall=2), "]", sep="")




# usefull functions ---------------------------------------------------  -------------------------------------------------------

place=function(i, j){
  if(i<j){pl=(2*Nc-i)*(i-1)/2+j-i}
  if(j<i){
    j1=j; i1=i;i=j1; j=i1 
    pl=(2*Nc-i)*(i-1)/2+j-i}
  return(pl)}

which.place=function(n){
  for (k in 1:(Nc-1)){
    for (l in (k+1):Nc){
      if(place(k,l)==n){pl=(paste(k,l))
      compo=paste(components[k],"-", components[l])
      }
      
    }}
  return(c(pl, compo))}

all.interactions=function(combination){
  ## eg combination=c(2,4,5,11,12,15,16)
  nc=length(combination)
  combination=sort(combination)
  interactions=c()
  for ( i in 1:(nc-1)){
    for (j in (i+1):nc){
      interactions=c(interactions,place(combination[i],combination[j]))      }    }
  return(interactions)}




# define interactions -----------------------------------------------------
alldata3=alldata2[,c(1:(which(colnames(alldata2)==components[1])-1))]


for (i in 1:Nc){  eval(parse(text=paste("alldata3$c1.", i, "=0", sep="")))}
for (i in 1:Nc){  eval(parse(text=paste("alldata3$c2.", i, "=0", sep="")))}
#grepl(value, chars, fixed = TRUE)
for (i in 1:Nc){
  for (j in 1:length(alldata3$id)){
    alldata3[j, which(colnames(alldata3)==paste("c1.", i, sep=""))]=
      1* grepl(components[i], alldata3$treat1[j], fixed = TRUE)
    
    alldata3[j, which(colnames(alldata3)==paste("c2.", i, sep=""))]=
      1* grepl(components[i], alldata3$treat2[j], fixed = TRUE)
  }
}

Ninter=Nc*(Nc-1)/2
for (i in 1:Ninter){eval(parse(text=paste("alldata3$int", i, "=0", sep="")))}

  for (j in 1:length(alldata3$id)){
    for (i in 1:Ninter){
      k1=as.numeric(word(which.place(i),1,sep = " ")[1])
      k2=as.numeric(word(which.place(i),2,sep = " ")[1])
    alldata3[j, which(colnames(alldata3)==paste("int", i, sep=""))]=
      alldata3[j, which(colnames(alldata3)==paste("c1.",k1, sep=""))]*alldata3[j, which(colnames(alldata3)==paste("c1.",k2, sep=""))]-
      alldata3[j, which(colnames(alldata3)==paste("c2.",k1, sep=""))]*alldata3[j, which(colnames(alldata3)==paste("c2.",k2, sep=""))]
  }
}


interactions=alldata3[, which(colnames(alldata3)=="int1"):which(colnames(alldata3)==paste("int", Ninter, sep=""))]

# SSVS interactions model ---------------------------------------------------  -------------------------------------------------------
model1.string <-  "
model {
for(i in 1:Ns) { 
prec[i]<-1/(se[i]*se[i])
y[i]~dnorm(phi[i],prec[i])
}

for(i in 1:Ns) { 
phi[i]~dnorm(mean[i], prec.tau)
mean[i]<- inprod(d[], c[i,])+inprod(gamma[], interactions[i,])
}

##prior distribution for heterogeneity	
tau ~ dnorm(0,0.1)I(0,)                                      
prec.tau<- 1/tau.sq
tau.sq<- pow(tau,2)
for(i in 1:Nc){
d[i]~dnorm(0,0.01)
}

## SSVS
for(k in 1:Ninter){
IndA[k] ~ dcat(Pind[])
Ind[k] <- IndA[k] - 1
gamma[k] ~ dnorm(0, tauCov[IndA[k]])  }

zeta <- pow(eta, -2)
eta ~ dnorm(0,1000)I(0,)
tauCov[1] <- zeta
tauCov[2] <- zeta * 0.01# g = 100

## all interactions equiprobable
Pind[1] <- 0.5 
Pind[2] <- 0.5 


### Example comparisons
## use the function all.interactions to find the which gamma to include
## e.g. run all.interactions(c(2,3,4,5,8))

## thal + dex VS thal + bor + dex
example1<- d[1]+d[3]+gamma[2]-d[1]-d[2]-d[3] -gamma[1]-gamma[2]-gamma[18]           

 ## bor + dex + peri VS thal + IFN
example2<- d[1]+d[2]+d[16]+gamma[1]+
gamma[15]+gamma[31]-d[3]-d[9]-gamma[39]        

## ob + dex VS dara + dex + len
example3<- d[1]+d[11]+gamma[10]-d[14]-d[1]-d[4] -
gamma[3]-gamma[13]-gamma[58]    

## bor + dex + cyc VS thal + IFN
example4<- d[1]+d[2]+d[18]+gamma[1]+gamma[17]+
gamma[33]-d[3]-d[9] -gamma[39]      
}
"


model1.spec<-textConnection(model1.string) 
data <- list(y=alldata2$logHR,se=alldata2$se.logHR, Nc=Nc, Ns=ns,c=c, 
             interactions=interactions, Ninter=Ninter)
jags.m=0
jags.m <- jags.model(model1.spec, data = data, n.chains =4, n.adapt =10000)



params <- c("tau", "d" , "gamma", "Ind", "example1", "example2", "example3", "example4", "eta") 
closeAllConnections()
samps2<- coda.samples(jags.m, params, n.iter =30000)
MCMCtrace(samps2,pdf = FALSE, params = "d") 


A2= MCMCsummary(samps2)
A2["HR.1",]=exp(A2["example1",])
A2["HR.2",]=exp(A2["example2",])
A2["HR.3",]=exp(A2["example3",])
A2["HR.4",]=exp(A2["example4",])
rownames(A2)[which(rownames(A2)=="d[1]"):which(rownames(A2)=="d[18]")]=components
inter=round(A2[, c(4,3,5)], digits=3)

meanInd=c() 
for(i in 1:Ninter){
  meanInd=c(meanInd, mean(rbind(samps2[[1]],samps2[[2]], samps2[[3]], samps2[[4]]) [, paste("Ind[", i,"]", sep="")]))}

medianG=c() 
for(i in 1:Ninter){  medianG=c(medianG, median(rbind(samps2[[1]],samps2[[2]], samps2[[3]], samps2[[4]])[, paste("gamma[", i,"]", sep="")]))}
nmax=which(medianG==max(medianG))
nmin=which(medianG==min(medianG))
order(medianG)
medianG[nmax]
medianG[nmin]
meanInd[nmax]
meanInd[order(-medianG)[1]] #   "len - carf" 
meanInd[order(-medianG)[2]] #  "dex - bor"
meanInd[nmin] #"dex - pom"

which.place(nmax)
which.place(nmin)

dev.off()
ggplot(data.frame("gamma"=medianG, "Ind"=meanInd), aes(x=medianG, y=meanInd)) + 
  geom_point() + labs(x = "Estimated coefficent of interaction terms (SSVS)")+
  labs(y = "Variable selection frequency for interaction terms (SSVS)")

which.place(nmin) # ""dex - pom" weak negative interaction
which.place(nmax) # "len - carf"  strong positive interaction
which.place(1) # "dex - bor"  strong positive interaction
A2$results=paste(format(round(A2$`50%`,2), nsmall=2), "[", 
                 format(round(A2$`2.5%`,2), nsmall=2), ";", 
                 format(round(A2$`97.5%`,2), nsmall=2), "]", sep="")




# LASSO interactions model ---------------------------------------------------  -------------------------------------------------------
model1.string <-  "
model {
for(i in 1:Ns) { 
prec[i]<-1/(se[i]*se[i])
y[i]~dnorm(phi[i],prec[i])
}

for(i in 1:Ns) { 
phi[i]~dnorm(mean[i], prec.tau)
mean[i]<- inprod(d[], c[i,])+inprod(gamma[], interactions[i,])
}

##prior distribution for heterogeneity	
tau ~ dnorm(0,0.1)I(0,)                                      
prec.tau<- 1/tau.sq
tau.sq<- pow(tau,2)
for(i in 1:Nc){
d[i]~dnorm(0,0.01)
}

## Bayesian LASSO
tauGamma <- pow(sdGamma,-1)
sdGamma ~ dunif(0, 5)
for(k in 1:Ninter){gamma[k] ~ ddexp(0, tauGamma) }


### Example comparisons
## use the function all.interactions to find the which gamma to include
## e.g. run all.interactions(c(2,3,4,5,8))

## thal + dex VS thal + bor + dex
example1<- d[1]+d[3]+gamma[2]-d[1]-d[2]-d[3] -gamma[1]-gamma[2]-gamma[18]           

 ## bor + dex + peri VS thal + IFN
example2<- d[1]+d[2]+d[16]+gamma[1]+
gamma[15]+gamma[31]-d[3]-d[9]-gamma[39]        

## ob + dex VS dara + dex + len
example3<- d[1]+d[11]+gamma[10]-d[14]-d[1]-d[4] -
gamma[3]-gamma[13]-gamma[58]    

## bor + dex + cyc VS thal + IFN
example4<- d[1]+d[2]+d[18]+gamma[1]+gamma[17]+
gamma[33]-d[3]-d[9] -gamma[39]          
}
"


model1.spec<-textConnection(model1.string) 
data <- list(y=alldata2$logHR,se=alldata2$se.logHR, Nc=Nc, Ns=ns,c=c, interactions=interactions, 
             Ninter=Ninter)
jags.m=0
jags.m <- jags.model(model1.spec, data = data, n.chains =4, n.adapt =10000)



params <- c("tau", "d" , "gamma",  "example1", "example2", "example3", "example4") 
closeAllConnections()
samps3<- coda.samples(jags.m, params, n.iter =30000)
MCMCtrace(samps3,pdf = FALSE, params = "d") 


A3= MCMCsummary(samps3)
A3["HR.1",]=exp(A3["example1",])
A3["HR.2",]=exp(A3["example2",])
A3["HR.3",]=exp(A3["example3",])
A3["HR.4",]=exp(A3["example4",])
rownames(A3)[which(rownames(A3)=="d[1]"):which(rownames(A3)=="d[18]")]=components
interLASSO=round(A3[, c(4,3,5)], digits=3)

medianG=c() 
for(i in 1:Ninter){  medianG=c(medianG, median(rbind(samps3[[1]],samps3[[2]], samps3[[3]], samps3[[4]])[, paste("gamma[", i,"]", sep="")]))}
nmax=which(medianG==max(medianG))
nmin=which(medianG==min(medianG))
order(medianG)
medianG[nmax]
medianG[nmin]
meanInd[nmax]

which.place(nmax)
which.place(nmin)
A3$results=paste(format(round(A3$`50%`,2), nsmall=2), "[", 
                 format(round(A3$`2.5%`,2), nsmall=2), ";", 
                 format(round(A3$`97.5%`,2), nsmall=2), "]", sep="")


# summarize results -------------------------------------------------------
summ=cbind(A1[c(which(rownames(A1)=="dex"):which(rownames(A1)=="cyc"), which(rownames(A1)=="tau"),
                which(rownames(A1)=="HR.1"),which(rownames(A1)=="HR.2"),
                which(rownames(A1)=="HR.3"),which(rownames(A1)=="HR.4")),c("results")], 
           A2[c(which(rownames(A2)=="dex"):which(rownames(A2)=="cyc"), which(rownames(A2)=="tau"),
                which(rownames(A2)=="HR.1"),which(rownames(A2)=="HR.2"),
                which(rownames(A2)=="HR.3"),which(rownames(A2)=="HR.4")),c("results")]
           , 
           A3[c(which(rownames(A3)=="dex"):which(rownames(A3)=="cyc"), which(rownames(A3)=="tau"),
                which(rownames(A3)=="HR.1"),which(rownames(A3)=="HR.2"),
                which(rownames(A3)=="HR.3"),which(rownames(A3)=="HR.4")),c("results")]
           )
rownames(summ)=c(components, "tau", "comparison I", "comparison II", "comparison III", "comparison IV")


#write.xlsx2(summ, file="results myeloma.xlsx")
