library(ggplot2)
library(xts)
library(quantmod)
# Read in the data from CSV files
GC1 <- read.csv("~/Desktop/MDP/DATA/GC1RB.csv",stringsAsFactors = FALSE)
SP1 <- read.csv("~/Desktop/MDP/DATA/SP1RB.csv",stringsAsFactors = FALSE)
TY1 <- read.csv("~/Desktop/MDP/DATA/TY1RB.csv",stringsAsFactors = FALSE)

# Calculate the monthly returns
GC1$Return <- ROC(as.numeric(GC1$PX_LAST),type='discrete',na.pad=TRUE)
SP1$Return <- ROC(as.numeric(SP1$PX_LAST),type='discrete',na.pad=TRUE)
TY1$Return <- ROC(as.numeric(TY1$PX_LAST),type='discrete',na.pad=TRUE)
# View(GC1)
Returns <- cbind(GC1$Return,SP1$Return,TY1$Return)
Returns[is.na(Returns)==TRUE] <- 0
colnames(Returns) <- c('GC1','SP1','TY1')

# Lookback period to be set here
lookback <- 120

# Calculate the equal weighted return
equal.weight.return <- apply(Returns,MARGIN=1,FUN = mean)
equal.weight.return[1:lookback] <- 0
equal.weight.portfolio <- cumprod(1+equal.weight.return)

# Calculate the equal risk parity
risk.parity.weights <- matrix(data = 0, nrow = nrow(Returns), ncol = ncol(Returns))
colnames(risk.parity.weights) <- c('GC1','SP01','TY01')
for(i in (lookback+1):nrow(Returns)){
  std.dev <- sqrt(diag(cov(Returns[(i-lookback):(i-1),])))
  sigma.inv <- std.dev^-1
  risk.parity.weights[i,] = sigma.inv/sum(sigma.inv)
}
risk.parity.return <- apply(Returns*risk.parity.weights,MARGIN=1,FUN = sum)
risk.parity.portfolio <- cumprod(1+risk.parity.return)

# Calculate the maximum diversification
max.divers.weights <- matrix(data = 0, nrow = nrow(Returns), ncol = ncol(Returns))
colnames(max.divers.weights) <- c('GC1','SP1','TY1')
for(i in (lookback+1):nrow(Returns)){
  V <- cov(Returns[(i-30):(i-1),])
  std.dev <- sqrt(diag(V))
  sigma.inv <- std.dev^-1
  C <- cov2cor(V)
  max.divers.weights[i,] <- (sigma.inv %*% (solve(C)))/sum((sigma.inv %*% (solve(C))))
}
max.divers.return <- apply(Returns*max.divers.weights,MARGIN=1,FUN = sum)
max.divers.portfolio <- cumprod(1+max.divers.return)

# Create report
port.returns <- cbind(equal.weight.portfolio,risk.parity.portfolio,max.divers.portfolio)
PDF <- data.frame(cbind(GC1$Index,port.returns))
pg <- ggplot(PDF, aes(x=1:nrow(Returns))) + 
  geom_line(aes(y = equal.weight.portfolio, colour = "EW")) + 
  geom_line(aes(y = risk.parity.portfolio, colour = "RP")) +
  geom_line(aes(y = max.divers.portfolio, colour = "MDP"))
print(pg)