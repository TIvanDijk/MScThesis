###############################################################################
#                             FUNCTIONS                                       #
# Contains all functions used in main.R & sensitivity.R                       #
############################################################################### 

# functions to work with the dataset
splitData <- function(data, order){
  y.train = data[(order + 1):length(data)]
  x.train = matrix(NA, nrow = (length(data)-order), ncol = order)
  
  for (row in 1:nrow(x.train)){
    x.train[row, ] = data[row:(row + order - 1)]
  }
  
  return( list(x = x.train, y = y.train))
}

doRidge <- function(param, lambda, order, newx, seed = 42){
  param.split = splitData(param, order)
  x.param = param.split[[1]]; y.param = param.split[[2]]

  set.seed(42)
  if (length(lambda > 1)){
    ridge_cv <- cv.glmnet(x.param, y.param, alpha = 0, lambda = lambda)
    lambda.opt = ridge_cv$lambda.min
    #cat('\n Optimal lambda (for given range) of: ', lambda.opt, '\n')
  } else {lambda.opt = lambda}
  
  param.ridge = glmnet(x.param, y.param, alpha = 0, lambda = lambda.opt)
  #print(coef(param.ridge))
  
  param.hat = c()
  all.params = y.param
  for (i in 1:newx){
    p.i <- as.numeric(predict(param.ridge, s = lambda.opt, newx = matrix(tail(all.params, order), nrow = 1)))
    param.hat = c(param.hat, p.i)
    all.params = c(all.params, p.i)
  }
  
  return( list(param.hat = param.hat, param = y.param))
}

# calculate MSE 
mse <- function(actual, pred){
  mean( (actual - pred)^2 )
}

# functions for the SIR model and their extensions 
classical.SIR <- function(time, state, parameters){
  with(as.list( c(state,parameters) ),{
    N=S+R+I
    dS=-beta*S*I/N
    dI=beta*S*I/N - gamma*I
    dR=gamma*I
    
    return(list(c(dS,dI,dR)))
  }
  )
}

timedependent.SIR <- function(time, state, parameters){
  with(as.list( c(state, parameters) ),{
    beta = eval(parse(text = paste0('beta', round(time, 0))))
    gamma = eval(parse(text = paste0('gamma', round(time, 0))))
    
    N=S+R+I
    dS=-beta*S*I/N
    dI=beta*S*I/N - gamma*I
    dR=gamma*I
    
    return(list(c(dS,dI,dR)))
  }
  )
}

asymptomatic.SIR <- function(time, state, parameters){
  with(as.list( c(state, parameters) ),{
    beta.s = eval(parse(text = paste0('beta.s', round(time, 0))))
    gamma = eval(parse(text = paste0('gamma', round(time, 0))))
    
    N=S+R+Is+A
    dS = -(beta.s*Is + beta.a*A)*S/N
    dIs = (beta.s*Is + beta.a*A)*S*w.s/N - gamma * Is
    dA = (beta.s*Is + beta.a*A)*S*w.a/N - gamma * A
    dR=gamma*(Is + A)
    
    return(list(c(dS,dIs,dA,dR)))
  }
  )
}

asymptomatic2.SIR <- function(time, state, parameters){
  with(as.list( c(state, parameters) ),{
    beta.a = eval(parse(text = paste0('beta.a', round(time, 0))))
    beta.s = eval(parse(text = paste0('beta.s', round(time, 0))))
    gamma = eval(parse(text = paste0('gamma', round(time, 0))))
    
    N=S+R+Is+A
    dS = -(beta.s*Is + beta.a*A)*S/N
    dIs = (beta.s*Is + beta.a*A)*S*w.s/N - gamma * Is
    dA = (beta.s*Is + beta.a*A)*S*w.a/N - gamma * A
    dR=gamma*(Is + A)
    
    return(list(c(dS,dIs,dA,dR)))
  }
  )
}

vaccinated.SIR <- function(time, state, parameters){
  with(as.list( c(state, parameters)), {
    alpha = eval(parse(text = paste0('alpha', round(time, 0))))
    beta = eval(parse(text = paste0('beta', round(time, 0))))
    gamma = eval(parse(text = paste0('gamma', round(time, 0))))
    
    N = S + I + R + V
    dS=-(beta*S*I)/N + alpha*S
    dI=(beta*S*I + sigma*beta*V*I)/N - gamma*I
    dR=gamma*I
    dV=alpha*S - sigma*beta*V*I/N
    
    return(list(c(dS, dI, dR, dV)))
  } )
}




