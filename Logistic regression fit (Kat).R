# -------------------------------------          Bernoulli data:
AZT = c(rep(1,170),rep(0,168))
AIDS = c(rep(1,25),rep(0,145),rep(1,44),rep(0,124))
#Creating a dataframe with data
df = data.frame(AZT,AIDS)
print(df)

xy <- scan('surgery.dat',what=list(x=0,y=0))
  x<- AZT
  y<- AIDS
  
  (x)
  (y)

par(mfrow=c(2,2))
plot(x,y,xlab='treatment',ylab='death status',type='n',ylim=c(0,1))
  points(x,y,cex=.5)
  title(expression('Binary data'))

# glm
x0 <- x
xglm <- glm(y~x0,family=binomial)
  b<- summary(xglm)$coef
  print(b)
  lines(x,xglm$fit)

# 2d likelihood

b0hat<-xglm$coef[1]; se0<- b[1,2]
b1hat<-xglm$coef[2]; se1<- b[2,2]

x0 <- x
fun1 <- function(bo,b1){
	sum(y*(bo+b1*x0)- log(1+exp(bo+b1*x0)))
        }
  lik<- NULL
  bbo<- seq(b0hat-3*se0, b0hat+3*se0 ,len=20)
  bb1 <- seq(b1hat-3*se1, b1hat+3*se1,len=20)
  for (bo in bbo)
  {
    for (b1 in bb1){
       lik <- c(lik,fun1(bo,b1))
    }
  }
  maxlik <- max(lik)
  lik <- lik- maxlik
  lik<- exp(lik)

contour(bbo,bb1,matrix(lik,20,byrow=T),level=seq(.1,1,by=.2),
        xlab=expression(beta[0]),
        ylab=expression(beta[1]))
  title(expression('Contour of the likelihood'))
  abline(h=xglm$coef[2],v=xglm$coef[1],lwd=.3)

lik<- matrix(lik,20,byrow=T)
  lik1 <- apply(lik,2,max)


plot(bb1,lik1,type='n',xlab='beta<-1',ylab='likelihood')
  lines(bb1,lik1,lwd=.4)
  abline(h=.15,lwd=.3)
  source('li.r')
  cat('95% CI for beta1=',round(li(bb1,lik1,.15),3),'\n')
  title(expression('(c) Profile likelihood'))

# normal approx
b<- summary(xglm)$coef
  se<- b[2,2]
ll<- log(lik1)
  ll2<- -1/2/se^2*(bb1-b[2,1])^2
plot(bb1,ll,type='n',xlab=expression(beta[1]),ylim=c(-6,0),
     ylab='log-likelihood')
  lines(bb1,ll,lwd=.5)
  lines(bb1,ll2,lty=2)
  title(expression('(d) Quadratic approximation'))
  
  (b0hat)
  (b1hat)

