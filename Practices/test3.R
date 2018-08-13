data1<-seq(1:100)
data2<-1+2*data1+rnorm(100)
data3<-data1+100*rnorm(100)
data4<-data1+10*rnorm(100)



routine<-function(lm_func,input){
  output<-lm_func(input)
  return(output)
}


majorfunc<-function(){
  lm_temp<-function(x){
    return(lm(x ~ data1))
  }
  a<-routine(lm_temp,data2)
  return(a)
}


b<-majorfunc()
summary(b)


a<-lm(data2~data1+data3+data4)

c<-vif(a)
