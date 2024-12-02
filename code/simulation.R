geneData=function(N=40,Nag=3,Nab=15,sd,type=1){
  n=N*Nag*Nab
  subject=paste("sub",factor(paste(formatC(seq(1, N), width=2, flag=0),sep = "")), sep = "")
  data=data.frame()
  for (i in 1:N) {
    antigen=rep(factor(paste(formatC(seq(1, Nag), width=2, flag=0),sep = "")),rep(Nab,Nag))
    antibody=rep(factor(paste(formatC(seq(1, Nab), width=2, flag=0),sep = "")),Nag)
    datai=data.frame(subjectid=subject[i],antigen=antigen,antibody=antibody)
    data=rbind(data,datai)
  }

  data$treatment=factor(c(rep(1,nrow(data)/2),rep(0,nrow(data)/2)))
  ff=as.formula("~treatment+antigen+antibody+treatment:antigen+
                treatment:antibody+antigen:antibody+treatment:antigen:antibody")
n1=Nag
n2=Nab
nn=n1+n2
nnn=nn+1*(n1-1)+1*(n2-1)+(n1-1)*(n2-1)
  sele=list(antigen=c(2,(nn+1):(nn+n1-1)),antibody=c(2,(nn+n1):(nn+n1-1+n2-1)),
            three=c(2,(nnn+1):(nnn+(n1-1)*(n2-1))))
  names(sele$antigen)=levels(data$antigen)
  names(sele$antibody)=levels(data$antibody)
## homo error, no random effects
  if (type==1){
errorVector=sd["Error"]*rnorm(n=n)
coef=c(0,rnorm(1),rnorm(n=Nag-1,sd=sd["Ag"]),rnorm(n=Nab-1,sd=sd["Ab"]),
       rnorm(n=Nag-1,sd=sd["Trag"]),rnorm(n=Nab-1,sd=sd["Trab"]), rnorm(n=(Nag-1)*(Nab-1),sd=sd["Agab"]),
       rnorm(n=(Nag-1)*(Nab-1),sd=sd["Tragab"]))
names(coef)=c("intercept","treatment",
              paste0("Ag",1:(Nag-1)),
              paste0("Ab",1:(Nab-1)),
              paste0("Trag",1:(Nag-1)),
              paste0("Trab",1:(Nab-1)),
              paste0("Agab",1:((Nab-1)*(Nag-1))),
              paste0("Tragab",1:((Nab-1)*(Nag-1)))
              )
MM=model.matrix(ff,data = data)
data$response= MM%*%coef+errorVector
#sele=list(antigen=c(2,19,20),detect_reagent=c(2,21:34),three=c(2,63:90))
bb=contrast(sele=sele,coeVec = coef,varMatrix = diag(1,length(coef)))
  }

## homo error, random intercept

  if (type==2){
    errorVector=sd["Error"]*rnorm(n=n)
    coef=c(0,rnorm(1),rnorm(n=Nag-1,sd=sd["Ag"]),rnorm(n=Nab-1,sd=sd["Ab"]),
           rnorm(n=Nag-1,sd=sd["Trag"]),rnorm(n=Nab-1,sd=sd["Trab"]), rnorm(n=(Nag-1)*(Nab-1),sd=sd["Agab"]),
           rnorm(n=(Nag-1)*(Nab-1),sd=sd["Tragab"]))
    names(coef)=c("intercept","treatment",
                  paste0("Ag",1:(Nag-1)),
                  paste0("Ab",1:(Nab-1)),
                  paste0("Trag",1:(Nag-1)),
                  paste0("Trab",1:(Nab-1)),
                  paste0("Agab",1:((Nab-1)*(Nag-1))),
                  paste0("Tragab",1:((Nab-1)*(Nag-1)))
    )
    MM=model.matrix(ff,data = data)
    rdata=data.frame(subjectid=subject,rintercept=rnorm(N,sd=sd["Rdin"]))
    data=merge(data,rdata)
    data$response= MM%*%coef+errorVector+data$rintercept
    #sele=list(antigen=c(2,19,20),detect_reagent=c(2,21:34),three=c(2,63:90))
    bb=contrast(sele=sele,coeVec = coef,varMatrix = diag(1,length(coef)))
  }

  ## homo error, random intercept, random antigen effects and random antibody effects
  if (type==3){

    errorVector=sd["Error"]*rnorm(n=n)
    coef=c(0,rnorm(1),rnorm(n=Nag-1,sd=sd["Ag"]),rnorm(n=Nab-1,sd=sd["Ab"]),
           rnorm(n=Nag-1,sd=sd["Trag"]),rnorm(n=Nab-1,sd=sd["Trab"]), rnorm(n=(Nag-1)*(Nab-1),sd=sd["Agab"]),
           rnorm(n=(Nag-1)*(Nab-1),sd=sd["Tragab"]))
    names(coef)=c("intercept","treatment",
                  paste0("Ag",1:(Nag-1)),
                  paste0("Ab",1:(Nab-1)),
                  paste0("Trag",1:(Nag-1)),
                  paste0("Trab",1:(Nab-1)),
                  paste0("Agab",1:((Nab-1)*(Nag-1))),
                  paste0("Tragab",1:((Nab-1)*(Nag-1)))
    )
    MM=model.matrix(ff,data = data)
    rdata1=data.frame(subjectid=subject,rintercept=rnorm(N,sd=sd["Rdin"]))
    rdata2=data.frame(matrix(rnorm(N*(Nag-1),sd=sd["Rdag"]),nrow = N))
    colnames(rdata2)=paste0("Rdag",2:Nag)
    rdata3=data.frame(matrix(rnorm(N*(Nab-1),sd=sd["Rdab"]),nrow = N))
    colnames(rdata3)=paste0("Rdab",2:Nab)
    rdata=cbind(rdata1,rdata2)
    rrdata=data.frame(subjectid=rdata$subjectid,rsum=rowSums(rdata[,-1]))
    data=merge(data,rrdata)
    data$response= MM%*%coef+errorVector+data$rsum
    #sele=list(antigen=c(2,19,20),detect_reagent=c(2,21:34),three=c(2,63:90))
    bb=contrast(sele=sele,coeVec = coef,varMatrix = diag(1,length(coef)))
  }

  return(list(data=data,coef=coef,modelMatrix=MM,errorVector=errorVector,effects=bb$effect))
}


