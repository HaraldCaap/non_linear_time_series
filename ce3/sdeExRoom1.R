sdeExRoom1 <- function(data, yTi,Ph){
  data$yTi <- yTi
  data$Ph <- Ph
  
  # Generate a new object of class ctsm
  model = ctsm()
  
  # Add a system equation and thereby also a state
  model$addSystem(dTi ~  1/Ci*(1/Ria*(Ta-Ti) + 1/Rim*(Tm-Ti)  + H + (a1*bs1+a2*bs2+a3*bs3+a4*bs4+a5*bs5)*Gv)*dt + exp(p11)*dw1)
  
  model$addSystem(dH ~ 1/tau *(K*Ph-H) *dt + exp(phh)*dwh)
  model$addSystem(dTm ~  1/Cm*(1/Rim*(Ti-Tm))*dt + exp(p22)*dw2)
  # Set the names of the inputs
  model$addInput(Ta,Gv,Ph,bs1,bs2,bs3,bs4,bs5)
  # Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(yTi ~ Ti)
  # Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  
  ##----------------------------------------------------------------
  # Set the initial value (for the optimization) of the value of the state at the starting time point
  model$setParameter(Ti = c(init = 15, lb = 0, ub = 40))
  model$setParameter(Tm = c(init = 15, lb = 0, ub = 40))
  model$setParameter(H = c(init = Ph[1], lb = 0, ub = 40))
  
  ##----------------------------------------------------------------
  # Set the initial value for the optimization
  model$setParameter(Ci = c(init = 1, lb = 1E-5, ub = 1E5))
  model$setParameter(Cm = c(init = 1000, lb = 1E-5, ub = 1E5))
  model$setParameter(Ria = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Rim = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(p11 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(p22 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(e11 = c(init = -1, lb = -50, ub = 10))
  model$setParameter(a1 = c(init = 1, lb = -100, ub = 100))
  model$setParameter(a2 = c(init = 1, lb = -100, ub = 100))
  model$setParameter(a3 = c(init = 1, lb = -100, ub = 100))
  model$setParameter(a4 = c(init = 1, lb = -100, ub = 100))
  model$setParameter(a5 = c(init = 1, lb = -100, ub = 100))
  model$setParameter(tau = c(init = 1, lb = 0, ub = 20))
  model$setParameter(phh = c(init = 1, lb = -30, ub = 10))
  model$setParameter(K = c(init=1, lb=0, ub=20))
  ##----------------------------------------------------------------   
  # Run the parameter optimization
  
  fit = model$estimate(data,firstorder = TRUE)
  return(fit)
}