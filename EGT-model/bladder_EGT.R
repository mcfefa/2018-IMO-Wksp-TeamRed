library(deSolve)
X0<-10^8; E0<-10^8; IL0<-100
yini<-c(X0,E0,IL0)
r<- 0.01; d1<- 10^(-10); e3<- 0.1; g3<- 100; a<- e3/g3; d2<- 3*10^(-2); e4<- 10^(-4); g4<- 10^(10); b<-e4/g4; d3<-10; m<-10^7; lambda1<-7*10^(-20); ILmin<-1 ; lambda2<-10*e4 ; lambda3<-200*d3;
m<-0; lambda1<-0; lambda2<-0; ILmin<-0; lambda3<-0; 
cancer <- function(t, y, parms) {
dy1 <- y[1] * (0.01 - 10^(-10)*y[2])
dy2 <- y[2] *( (0.1/100) * y[3] - 3*10^(-2))
dy3 <- 10^(-14) * y[1] * y[2] - 10*y[3]
list(c(dy1, dy2, dy3))
}
cancerpar <- function(t, y, parms) {
  dy1 <- y[1] * (r - d1*y[2])
  dy2 <- y[2] *(a * y[3] - d2)
  dy3 <- b * y[1] * y[2] - d3*y[3]
  list(c(dy1, dy2, dy3))
}
cancer2 <- function(t, y, parms) {
 #print(t)
  dy1 <- y[1] * (r - d1*y[2]/(1+y[1]/10^9))
  dy2 <- y[2] *(a * y[3] - d2) + m*d2-lambda1*((y[2])^3)
  dy3 <- b * y[1] * y[2] - d3*y[3] - lambda2*(y[3])^3 + d3 *ILmin + lambda3 * y[1] / (10^10 + y[1])
  list(c(dy1, dy2, dy3))}
cancer2rescaled <- function(t, y, parms) {
  #print(t)
  dy1 <- y[1] * (r - d1*y[2]/(1+y[1]/10^9))
  dy2 <- y[2] *(a * y[3]/10^6 - d2) + m*d2-lambda1*((y[2])^3)
  dy3 <- 10^6 * b * y[1] * y[2] - d3*y[3] - lambda2*((y[3])^3)/10^(12) + 10^6*d3*ILmin + 10^6*lambda3 * y[1] / (10^10 + y[1])
  list(c(dy1, dy2, dy3))}
cancer2fullyrescaled <- function(t, y, parms) {
  #print(t)
  dy1 <- y[1] * (r - d1*10^6*y[2]/(1+10^6*y[1]/10^9))
  dy2 <- y[2] *(a * y[3] - d2) + 10^(-6)*m*d2-lambda1*(10^(12)(y[2])^3)
  dy3 <- 10^(12)*b * y[1] * y[2] - d3*y[3] - lambda2*(y[3])^3 + d3 *ILmin + lambda3 * 10^6*y[1] / (10^10 + 10^6*y[1])
  list(c(dy1, dy2, dy3))}
times <- seq(from = 0, to = 10, by = 0.01)
times <- seq(from = 0, to = 200, by = 0.01)
times <- seq(from = 0, to = 400, by = 0.01)
out <- ode (times = times, y = yini, func = cancer, parms = NULL, method = rkMethod("rk45ck"))
out <- ode (times = times, y = yini, func = cancer2rescaled, parms = NULL, method = rkMethod("rk45ck"))
out <- ode (times = times, y = yini, func = cancer2, parms = NULL, method = rkMethod("rk45ck"))
head(out, n=5)
matplot(x = out[,1], y = out[,-1], type = "l", lwd = 2,
        lty = "solid", col = c("red", "blue", "black"),
        xlab = "time", ylab = "y", main = "Minimal resection")
legend("topleft", col = c("red", "blue", "black"),
       legend = c("Tumor", "T-cells", "IL"), lwd = 2)


