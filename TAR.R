data(lynx);
loglynx <- log10(lynx)

pacf(loglynx)
plot(loglynx,type="l")
m0=ar(loglynx)
m0=arima(loglynx,order=c(11,0,0))
acf(m0$res)
pacf(m0$res)
x=loglynx

########### arranged autoregression
n=length(x)
p = 3
#Perform the threshold F-test to decide whether it is necessary to use a TAR model. Perform the test for d = 1, 2, 3. If a TAR model is necessary, identify d.

# for (d in 1:3) {
d=2 
m=n%/%10+p
X=array(0,c(n-p,(p+1)))
xd=x[(p+1-d):(n-d)]
index=order(xd)
xdo=xd[index]
for (i in 1:(n-p)){
  ind=index[i]+p-d+d
  X[i,]=x[ind:(ind-p)]
}

res=array(0,c(n-p,2))
AR=array(0,c(n-p,p+1,3))
for (i in m:(n-p-1)){
  fit=lm(X[1:i,1]~X[1:i,2:(p+1)])
  xnew=as.vector(X[i+1,2:(p+1)])
  xnew=c(1,xnew)
  pred=sum(fit$coef[1:(p+1)]*xnew)
  fit.s=summary(fit)
  res[i+1,1]=X[i+1]-pred
  temp=1+t(xnew)%*%fit.s$cov%*%(xnew)
  res[i+1,2]=res[i+1,1]/sqrt(temp)
  #res[i+1,2]=res[i+1,2]/fit.s$sigma
  AR[i+1,,1]=fit.s$coef[,1]
  AR[i+1,,2]=fit.s$coef[,2]
  AR[i+1,,3]=fit.s$coef[,3]
}

res.m=lm(res[(m+1):(n-p),2]~X[(m+1):(n-p),2:(p+1)])
res.s=summary(res.m)
summary(res.m)
F=(sum(res[(m+1):(n-p),2]^2)-res.s$sigma^2*(n-p-m-p-1))/(p+1)/res.s$sigma^2
F
 print( 1-pf(F,p+1,(n-p-m-p-1)) ) # reject null hypothesis 

#there is a tar 
#d = 1   0.001591441
#d=2 0.0002000959
#d=3  0.002315181

plot(xdo[(m+1):(n-p-m)],res[(m+1):(n-p-m),1],main="Ordinary Predictive Residuals")
plot(xdo[(m+1):(n-p-m)],res[(m+1):(n-p-m),2],main="Standardized Predictive Residuals")

plot(xdo[(m+1):(n-p-m)],AR[(m+1):(n-p-m),3,3])


xdo[(m+1):(n-p-m)]

abline(v=c( 2.475671,2.685742,2.835056 ))

threshold1= c ( 2.475671, 2.475671 ,2.506505)
threshold2= c (  2.671173 ,2.674861, 2.685742 )
threshold3= c (2.835056 ,2.836957, 2.863917)
aicforth1 <- rep(0.0, 3)
aicforth2 <- rep(0.0, 3)
aicforth3 <- rep(0.0, 3)
th1<-2.506505
th2 <- 2.685742
th3<-2.835056
for(th3 in threshold3)
{
num1=length(which( xd< th1))
num2=length(which((xd>= th1) & (xd< th2)))
num3=length(which((xd>=th2) & (xd< th3 )))
num4=length(which(xd>=th3))
ind1=1:num1
ind2=(num1+1):(num1+num2)
ind3=(num1+num2+1):(num1+num2+num3)
ind4=(num1+num2+num3+1):(num1+num2+num3+num4)
## Regime 1
for (i in 2:4) {
  fit1=lm(X[ind1,1]~X[ind1,2:(i)])
  fit1.s=summary(fit1)
  print("Regime 1 ")
  print(AIC(fit1))
}
p1=3

fit1=lm(X[ind1,1]~X[ind1,2:(p1+1)])
fit1.s=summary(fit1)
## num1*log(fit1.s$sigma^2*(num1-p1-1)/num1) + 2*(p1+1)
AIC(fit1)
fit1.s
## Regime 2

for (i in 2:4) {
  fit2=lm(X[ind2,1]~X[ind2,2:i]) ## order identified by aic
  fit2.s=summary(fit2)
  print("Regime 2 ")
  print( AIC(fit2) )
}

fit2=lm(X[ind2,1]~X[ind2,2:4]) ## order identified by aic
fit2.s=summary(fit2)
AIC(fit2)
fit2.s

## Regime 3
for (i in 2:4) {
  fit3=lm(X[ind3,1]~X[ind3,2:i]) ## order identified by aic
  fit3.s=summary(fit3)
  print("Regime 3 ")
  print(AIC(fit3))
  
}
fit3=lm(X[ind3,1]~X[ind3,2:2]) ## order identified by aic
fit3.s=summary(fit3)
AIC(fit3)
fit3.s

## Regime 4
for (i in 2:4) {
  fit4=lm(X[ind4,1]~X[ind4,2:i]) ## order identified by aic
  fit4.s=summary(fit4)
  print("Regime 4 ")
  print(AIC(fit4))
}
fit4=lm(X[ind4,1]~X[ind4,2:2]) ## order identified by aic
fit4.s=summary(fit4)
AIC(fit4)
fit4.s
}
AIC(fit1)+AIC(fit2)+AIC(fit3)+AIC(fit4)

coef=array(0,c(4,4))
coef[1,1:4]=fit1$coef
coef[2,1:4]=fit2$coef
coef[3,1:2]=fit3$coef
coef[4,1:2]=fit4$coef
sigma=c(fit1.s$sigma*sqrt((num1-4)/num1),fit2.s$sigma*sqrt((num2-4)/num2),fit3.s$sigma*sqrt((num3-2)/num3),fit4.s$sigma*sqrt((num4-2)/num4))
tar=function(x){
  #2.475671,2.685742,2.835056 
  xx=c(1,x)
  if (x[1]<2.475671){y=coef[1,]%*%xx+rnorm(1,0,sigma[1])}
  else if (x[1]<2.685742)   {y=coef[2,]%*%xx+rnorm(1,0,sigma[2])}
  else if (x[1]< 2.835056 )  {y=coef[3,]%*%xx+rnorm(1,0,sigma[3])}
  else  {y=coef[4,]%*%xx+rnorm(1,0,sigma[4])}
}
coef[1,]
temp=x[n:(n-2)]

#Perform the parametric bootstrapping to make forecasts of horizon 1, 2, 3.
#Take T = 114 as the forecast origin.
h=3
B=5000
BS=array(0,c(h,B))
for (i in 1:B){
  temp=x[n:(n-2)]
  BS[1,i]=tar(temp)
  temp=c(BS[1,i],temp)
  temp=temp[1:3]
  BS[2,i]=tar(temp)
  temp=c(BS[2,i],temp)
  temp=temp[1:3]
  BS[3,i]=tar(temp)
}
mean(BS[1,])
mean(BS[2,])
mean(BS[3,])
###################### behaviro of t-ratio in a AR model

x=arima.sim(300,model=list(ar=c(.7)))
########### arranged autoregression
n=length(x)
p=3
d=2
m=n%/%10+p
X=array(0,c(n-p,(p+1)))
xd=x[(p+1-d):(n-d)]
index=order(xd)
xdo=xd[index]
for (i in 1:(n-p)){
  ind=index[i]+p-d+d
  X[i,]=x[ind:(ind-p)]
}

res=array(0,c(n-p,2))
AR=array(0,c(n-p,p+1,4))
for (i in m:(n-p-1)){
  fit=lm(X[1:i,1]~X[1:i,2:(p+1)])
  xnew=as.vector(X[i+1,2:(p+1)])
  xnew=c(1,xnew)
  pred=sum(fit$coef[1:(p+1)]*xnew)
  fit.s=summary(fit)
  res[i+1,1]=X[i+1]-pred
  temp=1+t(xnew)%*%fit.s$cov%*%(xnew)
  res[i+1,2]=res[i+1,1]/sqrt(temp)
  #res[i+1,2]=res[i+1,2]/fit.s$sigma
  AR[i+1,,1]=fit.s$coef[,1]
  AR[i+1,,2]=fit.s$coef[,2]
  AR[i+1,,3]=fit.s$coef[,3]
}


res.m=lm(res[(m+1):(n-p),2]~X[(m+1):(n-p),2:(p+1)])
res.s=summary(res.m)
summary(res.m)
F=(sum(res[(m+1):(n-p),2]^2)-res.s$sigma^2*(n-p-m-p-1))/(p+1)/res.s$sigma^2
F
1-pf(F,p+1,(n-p-m-p-1))

plot(xdo[(m+1):(n-p-m)],res[(m+1):(n-p-m),1],main="Ordinary Predictive Residuals")
plot(xdo[(m+1):(n-p-m)],res[(m+1):(n-p-m),2],main="Standardized Predictive Residuals")

plot(xdo[(m+1):(n-p-m)],AR[(m+1):(n-p-m),2,3])


