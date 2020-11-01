# ---------------------------------------------------------------------------#
# -                            function to run HyDaP                       - #
# -                                                                        - #
# -                                  10/28/2020                            - #
# -                                                                        - #
# ---------------------------------------------------------------------------#

library(cluster)

library(DataCombine)

library(klaR)

hydap<- function(data,conti.pos=NULL,cate.pos=NULL,k){
  
  # data is the data used for clustering contains clustering-related variables only, without missingness
  
  # conti.pos is a vector of the column number of continuous variables
  
  # cate.pos is a vector of the column number of categorical variables
  
  # k is pre-specified number of clusters
  
  
  n<-nrow(data)
  p<-ncol(data)
  
  if(any(apply(data, 2, function(x) any(is.na(x) | is.infinite(x))))){
    
    print('Missingness or Inf/-Inf exists in the data!')
    
    break
  }
  
  if(is.null(conti.pos)){
    
    output<- kmodes(data,modes=k)
    
    warning('As no continuous variable exists in the data, kmodes is used for clustering')
    
    
  }else if(is.null(cate.pos)){
    
    output<- kmeans(data,centers=k,nstart=20)
    
    warning('As no categorical variable exists in the data, kmeans is used for clustering')
    
  }else{
    
    # check whether conti.pos and cate.pos exceed data dimensionality
    
    if(max(conti.pos)>p | max(cate.pos)>p){
      
      print('Some number in conti.pos or cate.pos exceeds number of columns in the data!')
      
      break
    }
    
    # check whether conti.pos corresponds to continuous variables in data
    
    check.class<- apply(data[,conti.pos],2,class)
    
    if(!(length(unique(check.class))==1 & 'numeric' %in% check.class)){
      
      print('Some number in conti.pos does not correspond to numerical variable in the data!')
      
      break
    }
    
    pair<-data.frame(combn(n,2))
    
    gower<-matrix(rep(NA,ncol(pair)*p),ncol=p)
    
    for(i in 1:p){
      
      if(i %in% conti.pos){
        
        x<-data[,i]
        
        range_x<-max(x)-min(x)
        
        gower[,i]<-unlist(lapply(pair, function(y) abs(x[y[1]]-x[y[2]])/range_x)) # dissimilarity
        
      }else if(i%in% cate.pos){
        
        a<-data[,i]
        
        score<-lapply(pair, function(y){
          
          if(a[y[1]]!=a[y[2]]){
            1
          }else{0}
          
        }) # dissimilarity
        
        gower[,i]<-unlist(score)
      }
      
    }
    
    for(i in 1:p){
      
      gower[,i]<- gower[,i]/sum(gower[,i])
      
    }
    
    # change distance matrix to dist format
    
    matrixsum<- unlist(apply(gower,1,mean))
    
    ids<- t(pair)
    
    matrixsum2<- data.frame(ids,matrixsum)
    
    colnames(matrixsum2)<- c('v1','v2','gower_sum')
    
    for(j in 1:nrow(data_s)){
      
      new<- c(j,j,0)
      
      matrixsum2<- InsertRow(matrixsum2,NewRow = new)
      
    }
    
    matrixsum3<- matrixsum2[order(matrixsum2[,1],matrixsum2[,2]),]
    
    finallist1<- as.dist(xtabs(matrixsum3$gower_sum~matrixsum3$v2+matrixsum3$v1))
    
    pam1<- pam(finallist1,k,diss=T)
    
    output<- list(clustering=pam1$clustering,
                  dissMatrix=gower)
    
  }
 
  return(output)
  
}
