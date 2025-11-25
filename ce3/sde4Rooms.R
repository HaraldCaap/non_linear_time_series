sde4Rooms <- function(data){
  # Create model object
  model = ctsm()
  
  # Add states for rooms
  model$addSystem(dT1 ~  
                    1/(flArea1*C)*(1/(Ra/exArea1)*(Ta-T1) + 1/(Rm/flArea1)*(Tm1-T1)  
                                           + Ph1 
                                           + Aw1*Gv 
                                           + 1/(R12)*(T2-T1))*dt 
                  + exp(p11)*dw1
                  )
  model$addSystem(dTm1 ~  1/(flArea1*Cm)*(1/(Rm/flArea1)*(T1-Tm1))*dt + exp(pm11)*dwm1)
  
  model$addSystem(dT2 ~  
                    1/(flArea2*C)*(1/(Ra/exArea2)*(Ta-T2) + 1/(Rm/flArea2)*(Tm2-T2)  
                                      + Ph1 
                                      + Aw2*Gv 
                                      + 1/(R12)*(T1-T2) + 1/(R23)*(T3-T2))*dt 
                  + exp(p22)*dw2
  )
  model$addSystem(dTm2 ~  1/(flArea2*Cm)*(1/(Rm/flArea2)*(T2-Tm2))*dt + exp(pm22)*dwm2)
  
  model$addSystem(dT3 ~  
                    1/(flArea3*C)*(1/(Ra/exArea3)*(Ta-T3) + 1/(Rm/flArea3)*(Tm3-T3)  
                                      + Ph2 
                                      + Aw3*Gv 
                                      + 1/(R34)*(T4-T3) + 1/R23*(T2-T3))*dt 
                  + exp(p33)*dw3
  )
  model$addSystem(dTm3 ~  1/(flArea3*Cm)*(1/(Rm/flArea3)*(T3-Tm3))*dt + exp(pm33)*dwm3)
  
  model$addSystem(dT4 ~  
                    1/(flArea4*C)*(1/(Ra/exArea4)*(Ta-T4) + 1/(Rm/flArea4)*(Tm4-T4)  
                                      + Ph2 
                                      + Aw4*Gv 
                                      + 1/R34*(T3-T4))*dt 
                  + exp(p44)*dw4
  )
  model$addSystem(dTm4 ~  1/(flArea4*Cm)*(1/(Rm/flArea4)*(T4-Tm4))*dt + exp(pm44)*dwm4)
  
  # Add inputs
  model$addInput(Ta,Gv,Ph1,Ph2)
  
  # Add measurements
  model$addObs(yTi1 ~ T1)
  model$setVariance(yTi1 ~ exp(e11))
  
  model$addObs(yTi2 ~ T2)
  model$setVariance(yTi2 ~ exp(e22))
  
  model$addObs(yTi3 ~ T3)
  model$setVariance(yTi3 ~ exp(e33))
  
  model$addObs(yTi4 ~ T4)
  model$setVariance(yTi4 ~ exp(e44))
  
  ## Set parameters
  # Set the initial value (for the optimization) of the value of the state at the starting time point
  model$setParameter(T1 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(Tm1 = c(init = 15, lb = 0, ub = 40))
  #model$setParameter(H = c(init = Ph[1], lb = 0, ub = 40))
  model$setParameter(T2 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(Tm2 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(T3 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(Tm3 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(T4 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(Tm4 = c(init = 15, lb = 0, ub = 40))
  
  # Set the initial value for the optimization
  model$setParameter(C = c(init = 1, lb = 1E-5, ub = 1E5))
  model$setParameter(Cm = c(init = 1000, lb = 1E-5, ub = 1E5))
  model$setParameter(R12 = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(R23 = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(R34 = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Ra = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Rm = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(p11 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(p22 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(p33 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(p44 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(pm11 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(pm22 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(pm33 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(pm44 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(e11 = c(init = -1, lb = -50, ub = 10))
  model$setParameter(e22 = c(init = -1, lb = -50, ub = 10))
  model$setParameter(e33 = c(init = -1, lb = -50, ub = 10))
  model$setParameter(e44 = c(init = -1, lb = -50, ub = 10))
  model$setParameter(Aw1 = c(init = 6, lb = 1E-2, ub = 7.5+4.8+5))
  model$setParameter(Aw2 = c(init = 6, lb = 1E-2, ub = 7.5+4.8+5))
  model$setParameter(Aw3 = c(init = 6, lb = 1E-2, ub = 7.5+4.8+5))
  model$setParameter(Aw4 = c(init = 6, lb = 1E-2, ub = 7.5+4.8+5))
  
  # Set fixed parameters
  model$setParameter(flArea1 =c(init = 3))
  model$setParameter(flArea2 =c(init = 21))
  model$setParameter(flArea3 =c(init = 17))
  model$setParameter(flArea4 =c(init = 3))
  model$setParameter(exArea1 =c(init = 1.5))
  model$setParameter(exArea2 =c(init = 3))
  model$setParameter(exArea3 =c(init = 2.5))
  model$setParameter(exArea4 =c(init = 1.5))
  
  ##----------------------------------------------------------------   
  # Run the parameter optimization
  
  fit = model$estimate(data,firstorder = TRUE)
  return(fit)
}