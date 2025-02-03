library(glmmLasso)
set.seed(1909)
data(soccer)
N<-dim(soccer)[1]
ind<-sample(N,N)
lambda <- seq(500,0,by=-5)

kk<-5
nk <- floor(N/kk)

Devianz_ma<-matrix(Inf,ncol=kk,nrow=length(lambda))
family <- poisson(link = log)

## first fit good starting model
library(MASS);library(nlme)
PQL <- glmmPQL(points~1,random = ~1|team,family=family,data=soccer)

Delta.start <- as.matrix(t(rep(0,7+23)))
Q.start<-as.numeric(VarCorr(PQL)[1,1])


BIC_vec<-rep(Inf,length(lambda))

for (i in 1:kk)
{
  print(paste("CV Loop ", i,sep=""))

  if (i < kk)
  {
    indi <- ind[(i-1)*nk+(1:nk)]
  }else{
    indi <- ind[((i-1)*nk+1):N]
  }

  soccer.train<-soccer[-indi,]
  soccer.test<-soccer[indi,]

  Delta.temp <- Delta.start
  Q.temp <- Q.start
  ## loop over lambda grid
  for(j in 1:length(lambda))
  {
    #print(paste("Lambda Iteration ", j,sep=""))

    glm4 <- try(glmmLasso(points~transfer.spendings
                          + ave.unfair.score  + ball.possession
                          + tackles + ave.attend + sold.out, rnd = list(team=~1),
                          family = family, data = soccer.train, lambda=lambda[j],switch.NR=FALSE,final.re=FALSE,
                          control=list(start=Delta.temp[j,],q_start=Q.temp[j]))
                ,silent=TRUE)

    if(!inherits(glm4, "try-error"))
    {
      y.hat<-predict(glm4,soccer.test)
      Delta.temp<-rbind(Delta.temp,glm4$Deltamatrix[glm4$conv.step,])
      Q.temp<-c(Q.temp,glm4$Q_long[[glm4$conv.step+1]])

      Devianz_ma[j,i]<-sum(family$dev.resids(soccer.test$points,y.hat,wt=rep(1,length(y.hat))))
    }
  }
}

Devianz_vec<-apply(Devianz_ma,1,sum)
opt4<-which.min(Devianz_vec)

## now fit full model until optimal lambda (which is at opt4)
for(j in 1:opt4)
{
  glm4.big <- glmmLasso(points~transfer.spendings
                        + ave.unfair.score + ball.possession
                        + tackles + ave.attend + sold.out, rnd = list(team=~1),
                        family = family, data = soccer, lambda=lambda[j],switch.NR=FALSE,final.re=FALSE,
                        control=list(start=Delta.start[j,],q_start=Q.start[j]))
  Delta.start<-rbind(Delta.start,glm4.big$Deltamatrix[glm4.big$conv.step,])
  Q.start<-c(Q.start,glm4.big$Q_long[[glm4.big$conv.step+1]])
}

glm4_final <- glm4.big

summary(glm4_final)
