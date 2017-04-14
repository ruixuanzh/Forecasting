library(forecast)
library(tseries)
library(lawstat)
library(vars)

# read data in
df <- read.csv("train.csv")
train <- df[1:228,]
test <- df[229:288,]

# train
unemployment_train <- ts(train$Unemployment_Rate, start = 1987, frequency = 12)
pop_train <- ts(train$Population, start = 1987, frequency = 12)
bankruptcy_train <- ts(train$Bankruptcy_Rate, start = 1987, frequency = 12)
hprice_train <- ts(train$House_Price_Index, start = 1987, frequency = 12)

# test
unemployment_test <- ts(test$Unemployment_Rate, start = 2006, frequency = 12)
pop_test <- ts(test$Population, start = 2006, frequency = 12)
bankruptcy_test <- ts(test$Bankruptcy_Rate, start = 2006, frequency = 12)
hprice_test <- ts(test$House_Price_Index, start = 2006, frequency = 12)


# plot 
par(mfrow=c(4,1))
plot(bankruptcy)
acf(bankruptcy)
pacf(bankruptcy)

plot(unemployment)
plot(pop)
plot(bankruptcy)
plot(hprice)

adf.test(bankruptcy)
ndiffs(bankruptcy)
nsdiffs(diff(bankruptcy))


acf(diff(diff(bankruptcy),12),36)
pacf(diff(diff(bankruptcy),12),36)


# unit root test: augmented dicky-fuller test
adf.test(diff(unemployment))   # unit root of order 1
adf.test(diff(diff(pop)))   # unit root of order 2
adf.test(diff(bankruptcy))   # unit root of order 1
adf.test(diff(hprice))   # unit root of order 1

# univariate time series 
new <- diff(log(bankruptcy))

##################  SARIMA  ###################################

# implemented grid search to find the model with smallest test rmse and with residual has contant variance
# there's none model's leven test p_value larger than 0.001
vector = c()
p_vector = c()
q_vector = c()
P_vector = c()
Q_vector = c()
for (p in 1:3){
  for (q in 1:3){
    for (P in 1:3){
      for (Q in 1:3){
  m2x4 <- arima(bankruptcy_train, order = c(p,1,q), seasonal = list(order=c(P,1,Q), period = 12), method = 'CSS')
  fc1 <- forecast(m2x4, h = 60, level = 95)
  rmse <- sqrt(mean((fc1$mean - bankruptcy_test)^2))
  vector <- c(vector, rmse)
  p_vector <- c(p_vector, p)
  q_vector <- c(q_vector, q)
  P_vector <- c(P_vector, P)
  Q_vector <- c(Q_vector, Q)
  group <- c(rep(1,57),rep(2,57),rep(3,57),rep(4,57))
  variance_test <-levene.test(m2x4$residuals, group)
  print(c(p, q, P, Q, rmse))
  if (variance_test$p.value > 0.001) print("there's one")
     }
   }
  }
}


################## log  SARIMA  ###################################

# implemented grid search to find the model with smallest test rmse
# there's none model's leven test p_value larger than 0.001
vector <- c()
p_vector <- c()
q_vector <- c()
P_vector <- c()
Q_vector <- c()
whether_has_constant_variance <- c()
P_value_leven_test <- c()
for (p in 1:3){
  for (q in 1:3){
    for (P in 1:3){
      for (Q in 1:3){
        m2x4 <- arima(log(bankruptcy_train), order = c(p,1,q), seasonal = list(order=c(P,1,Q), period = 12), method = 'CSS')
        fc1 <- forecast(m2x4, h = 60, level = 95)
        rmse <- sqrt(mean((exp(fc1$mean) - bankruptcy_test)^2)) 
        vector <- c(vector, rmse)
        p_vector <- c(p_vector, p)
        q_vector <- c(q_vector, q)
        P_vector <- c(P_vector, P)
        Q_vector <- c(Q_vector, Q)
        group <- c(rep(1,57),rep(2,57),rep(3,57),rep(4,57))
        variance_test <-levene.test(m2x4$residuals, group)
        print(c(p, q, P, Q, rmse))
        result = 0
        if (variance_test$p.value > 0.01){
          result = 'yes'} else {
            result = 'no'
          }
        P_value_leven_test <- c(P_value_leven_test, variance_test$p.value)
        whether_has_constant_variance <- c(whether_has_constant_variance,result)
      }
    }
  }
}

log_SARIMA_result <- data.frame(p_vector, q_vector, P_vector, Q_vector, vector, whether_has_constant_variance, P_value_leven_test)
write.table(log_SARIMA_result, "project_log_SARIMA_result.csv", sep = ",")

# rank the models that gives us low test rmse and has constant variance, I also checked if the residual are serial uncorrelated or not
# c(2,1,3),c(2,1,3)
m2x4 <- arima(log(bankruptcy_train), 
              order = c(2,1,3), 
              seasonal = list(order=c(1,1,3), period = 12), method = 'CSS')
acf(m2x4$residuals,36)
pacf(m2x4$residuals,36)
tsdiag(m2x4)

group <- c(rep(1,57),rep(2,57),rep(3,57),rep(4,57))
variance_test <-levene.test(m2x4$residuals, group)
plot(m2x4$residuals)
fc1 <- forecast(m2x4, h = 60, level = 95)
plot(fc1)
sqrt(mean((exp(fc1$mean) - bankruptcy_test)^2))

################## log ARIMAX model ~ X2, X3 ###################################

# model selection : choose the best model that has smallest test mset, meanwhile satisfy all the underlying assumptions
vector <- c()
p_vector <- c()
q_vector <- c()
P_vector <- c()
Q_vector <- c()
whether_has_constant_variance <- c()
P_value_leven_test <- c()
for (p in 1:4){
  for (q in 1:4){
    for (P in 1:4){
      for (Q in 1:4){
        m2x4 <- arima(log(bankruptcy_train), order = c(p,1,q), seasonal = list(order=c(P,1,Q), period = 12), xreg = data.frame(log(hprice_train), log(unemployment_train)), method = 'CSS')
        fc1 <- forecast(m2x4, h = 60, xreg = data.frame(log(hprice_test), log(unemployment_test)), level = 95)
        rmse <- sqrt(mean((exp(fc1$mean) - bankruptcy_test)^2)) 
        vector <- c(vector, rmse)
        p_vector <- c(p_vector, p)
        q_vector <- c(q_vector, q)
        P_vector <- c(P_vector, P)
        Q_vector <- c(Q_vector, Q)
        group <- c(rep(1,57),rep(2,57),rep(3,57),rep(4,57))
        variance_test <-levene.test(m2x4$residuals, group)
        print(c(p, q, P, Q, rmse))
        result = 0
        if (variance_test$p.value > 0.01){
          result = 'yes'} else {
            result = 'no'
          }
        P_value_leven_test <- c(P_value_leven_test, variance_test$p.value)
        whether_has_constant_variance <- c(whether_has_constant_variance,result)
      }
    }
  }
}
log_ARIMAX_result <- data.frame(p_vector, q_vector, P_vector, Q_vector, vector, whether_has_constant_variance, P_value_levene_test)
write.table(log_ARIMAX_result, "project_log_log_ARIMAX_result.csv", sep = ",")



m2x4 <- arima(log(bankruptcy_train), order = c(4,1,3), seasonal = list(order=c(1,1,4), period = 12), xreg = data.frame(log(hprice_train), log(unemployment_train)), method = 'CSS')
fc1 <- forecast(m2x4, h = 60, xreg = data.frame(log(hprice_test), log(unemployment_test)), level = 95)
sqrt(mean((exp(fc1$mean) - bankruptcy_test)^2))

par(mfrow=c(2,1))

acf(m2x4$residuals,36)
pacf(m2x4$residuals,36)
tsdiag(m2x4)
group <- c(rep(1,57),rep(2,57),rep(3,57),rep(4,57))
levene.test(m2x4$residuals, group)

#################### VAR model #######################################
# winning model --> serial correlation
VARselect(y = data.frame(log(bankruptcy_train), log(unemployment_train), log(hprice_train)),season = 12)  # SIC: choose 3

# model1
var3 <- VAR(y = data.frame(log(bankruptcy_train), unemployment_train, hprice_train), p = 3, season =12 ) # test mse : 0.006717102

# model2
var5 <- VAR(y = data.frame(log(bankruptcy_train), unemployment_train, hprice_train, pop_train), p = 5, season =12 ) # test mse : 0.00699

# best model: model3
var3 <- VAR(y = data.frame(log(bankruptcy_train), log(unemployment_train), log(hprice_train)), p = 3, season =12 ) # test mse : 0.00535

# plot
plot(var3)
res <- var3$varresult$log.bankruptcy_train.$residuals
plot(res)
par(mfrow = c(2,1))
acf(res)
pacf(res)
Box.test(res, lag = 12, type = "Ljung-Box")

# heteroskedacity
group <- c(rep(1,56),rep(2,56),rep(3,56),rep(4,57))
levene.test(var3$varresult$log.bankruptcy_train.$residuals, group = group)

# Normal
shapiro.test(var3$varresult$log.bankruptcy_train.$residuals)

# predict on the test set
pred <- predict(var3, n.ahead = 60, ci = 0.95)
fc <- exp(pred$fcst$log.bankruptcy_train.[,1])
rmse <- sqrt(mean((bankruptcy_test- fc)^2))  # 0.00671
rmse












jpeg("qqplot.jpeg")
qqPlot(model2)
dev.off()


install.packages("apsrtable")
library(apsrtable)
s <- summary(m2x4)
apsrtable(s)


