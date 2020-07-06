library(glmnet)
library(ISLR)
library(MASS)
data =
###GENE4
GENE = function(k=5,n=100,imp=0.5,evo.prob=0.05,sex.prob=0.5,
                choose.type="unrank",data=NULL,model.type="lda",glm_family=binomial){
  ldaac = vector( length = k)
  pac = vector("list", k)
  same = c(rep(F,n))
  tfflag = matrix(nrow=length(data)-1,ncol=k)
  
  m1 = names(data)[length(data)]
  m2 = names(data)[-length(data)]

  for(i in 1:k){
    start = T
    while(start){
    tf = rbinom(length(m2),1,prob = 0.5)
    if(sum(tf)>0){
      start=F
    }
    tf2<-ifelse(tf==1,T,F)
    tfflag[,i] = tf2
    }
  }
  for(iteration in 1:n){
    train = sample(1:nrow(data),round(nrow(data)*4/5))
    trainWine = data[train,]
    testWine = data[-train,]
    for(i in 1:k){
      choose = m2[tfflag[,i]]
      a = paste(choose,sep = "",collapse = "+")
      b = paste(m1,"~",a,sep = "")
      switch (model.type,
              lda = {
                model = lda(as.formula(b),data= trainWine)
                model_pred = predict(model,testWine)
                ldaac[i]  = sum(diag(table(model_pred[[1]],testWine[,length(data)])))/length(testWine[,length(data)])
                
              },
                glm ={
                levels(trainWine[,length(data)]) = c(0,1)
                levels(testWine[,length(data)]) = c(0,1)
                
                model = glm(as.formula(b),family = glm_family,data= trainWine,control = list(maxit = 50))
                model_prob = predict(model,testWine,type="response")
                model_pred = rep(levels(trainWine[,length(data)])[1],length(model_prob))
                model_pred[model_prob>0.5] = levels(trainWine[,length(data)])[2]
                ldaac[i]  = sum(diag(table(model_pred,testWine[,length(data)])))/length(testWine[,length(data)])
                
              },
              stop("Wrong model")
      )
      
      pac[[i]]=choose
      names(ldaac)[i] = b
    }

    ##選擇
    tempFlag =switch (choose.type,
                      rank = sample(1:k,size=k,replace=T,prob=(imp*rank(round(1/apply(tfflag,2,sum),digit=4))
                      +(1-imp)*rank(ldaac))/(sum(imp*rank(round(1/apply(tfflag,2,sum),digit=4)))
                      +(1-imp)*sum(rank(ldaac))))
                      ,unrank = sample(1:k,size=k,replace=T,prob=(imp*round(1/apply(tfflag,2,sum),digit=4)
                       +(1-imp)*ldaac)/(sum(imp*round(1/apply(tfflag,2,sum),digit=4))+(1-imp)*sum(ldaac)))
                      ,{
                        stop("Wrong choose.type") }
    )  
    tfflagChoose = tfflag[,tempFlag]
    ##交配
    
    sexCUM = combn(1:k,2)[,sample(1:ncol(combn(1:k,2)),replace = F,prob = rep(1/ncol(combn(1:k,2))
                                                                              ,ncol(combn(1:k,2))) )]
    sexTF = sample(c(T,F),size =ncol(combn(1:k,2)),replace=T,prob=c(sex.prob,1-sex.prob) )
    sexFlag = sexCUM[,sexTF]
    sexPOS = sample(1:length(m2),size=sum(sexTF),replace=T,prob=rep(1/length(m2),length(m2)))
    if(sum(sexTF>0)){
      for(i in 1:length(sexPOS)){
        if(is.matrix(sexFlag)==T){
          for(j in 0:(length(m2)-sexPOS[i])){
          if(tfflagChoose[sexPOS[i]+j,sexFlag[1,i]]!=tfflagChoose[sexPOS[i]+j,sexFlag[2,i]]){
            tfflagChoose[sexPOS[i]+j,sexFlag[1,i]] = !tfflagChoose[sexPOS[i]+j,sexFlag[1,i]]
            tfflagChoose[sexPOS[i]+j,sexFlag[2,i]] = !tfflagChoose[sexPOS[i]+j,sexFlag[2,i]]
          }
          }
        }else{
          for(j in 0:(length(m2)-sexPOS[i])){
          if(tfflagChoose[sexPOS[i]+j,sexFlag[1]]!=tfflagChoose[sexPOS[i],sexFlag[2]]){
            tfflagChoose[sexPOS[i]+j,sexFlag[1]] = !tfflagChoose[sexPOS[i],sexFlag[2]]
            tfflagChoose[sexPOS[i]+j,sexFlag[1]] = !tfflagChoose[sexPOS[i],sexFlag[2]]
          }
          }
        }
      }
    }
    ##相同指標
    if(length(table(ldaac))==1){
      same[iteration]=T
    }
    ##終止
    if(iteration>4){
      if(n>100&same[iteration]==T&same[iteration-1]==T&same[iteration-2]&same[iteration-3]&same[iteration-4]){
        break
      }
    }
    

    ##突變
    for(i in 1:ncol(tfflagChoose)){
      evoFlag = rbinom(1,1,prob=evo.prob)
      if(evoFlag==1){
        evoPOS=sample(1:nrow(tfflagChoose),1,prob=rep(1/nrow(tfflagChoose),nrow(tfflagChoose)))
        tfflagChoose[evoPOS,i] = !tfflagChoose[evoPOS,i]
      }
    }
    ##防BUG
    for (i in ncol(tfflag) ) {
    if(sum(tfflagChoose[,i])>0)
      tfflag[,i] = tfflagChoose[,i]
    }
    
  }
  monkey = list(mod = ldaac[which.max(imp*rank(round(1/apply(tfflag,2,sum),digit=4))+(1-imp)*rank(ldaac))],para=pac[which.max(imp*rank(round(1/apply(tfflag,2,sum),digit=4))+(1-imp)*rank(ldaac))])
  return(monkey)
}

GENE4InBox = function(k=5,n=100,imp=0.5,evo.prob=0.05,sex.prob=0.5,choose.type="unrank",data=NULL,model.type="lda",glm_family="binomial",r=4){
  re = c()
  pac = vector("list", r)
  for(i in 1:r){
    model =GENE(k=k,n=n,imp = imp,evo.prob = evo.prob,sex.prob=sex.prob,choose.type = choose.type,data = data,model.type=model.type,glm_family=glm_family)
    re=c(re,model$mod)
    pac[[i]] = model$para
  } 
  ans = list(mod=re,para=pac,data=data)
  class(ans)="GENE"
  return(ans)
}


##CV

GENE.cv = function(GENE=NULL,k=5){
  a = vector("list",length(GENE$mod))
  names(a)= names(GENE$mod)
  ldaac = vector(length = k)
  for (i in 1:length(GENE$mod)) {
    data2 =GENE$data[sample(1:nrow(GENE$data),nrow(GENE$data)),]
    dataFlag = vector("list",k)
    for(x in 0:(k-1)){
      dataFlag[[x+1]] = ((round(nrow(GENE$data)/k))*x+1):(round((nrow(GENE$data)/k)*(x+1)))
    }
    

    for(j in 1:k){
      trainWine = data2[-dataFlag[[j]],]
      testWine = data2[dataFlag[[j]],]
      model = lda(as.formula(names(GENE$mod)[i]),data= trainWine)
      model_pred = predict(model,testWine)
      ldaac[j]  =sum(diag(table(model_pred[[1]],testWine[,length(testWine)])))/length(testWine[,length(testWine)])
      
      
    }
    a[[i]] = ldaac
    
  }
  return(a)
}


print.GENE = function(dick){
  print(dick$mod)
  
}