original_model <- function(data, yTi,Ph){
  data$yTi <- yTi
  data$Ph <- Ph
  # Generate a new object of class ctsm
  model = ctsm()
  # Add a system equation and thereby also a state
  # Gv in Ti: Aw/Ci*Gv or Tm: Aw/Cm*Gv
  model$addSystem(dTi ~  1/Ci*(1/Ria*(Ta-Ti) + 1/Rim*(Tm-Ti)  + Ph + Aw*Gv)*dt + exp(p11)*dw1)
  model$addSystem(dTm ~  1/Cm*(1/Rim*(Ti-Tm))*dt + exp(p22)*dw2)
  # Set the names of the inputs
  model$addInput(Ta,Gv,Ph)
  # Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(yTi ~ Ti)
  # Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ##----------------------------------------------------------------
  # Set the initial value (for the optimization) of the value of the state at the starting time point
  model$setParameter(Ti = c(init = 15, lb = 0, ub = 40))
  model$setParameter(Tm = c(init = 15, lb = 0, ub = 40))
  ##----------------------------------------------------------------
  # Set the initial value for the optimization
  model$setParameter(Ci = c(init = 1, lb = 1E-5, ub = 1E5))
  model$setParameter(Cm = c(init = 1000, lb = 1E-5, ub = 1E5))
  model$setParameter(Ria = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Rim = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Aw = c(init = 6, lb = 1E-2, ub = 7.5+4.8+5))
  model$setParameter(p11 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(p22 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(e11 = c(init = -1, lb = -50, ub = 10))
  ##----------------------------------------------------------------    
  
  # Run the parameter optimization
  
  fit = model$estimate(data,firstorder = TRUE)
  return(fit)
}


splines_model <- function(data, yTi, Ph) {
  data$yTi <- yTi
  data$Ph <- Ph
  
  model <- ctsm()
  
  model$addSystem(
    dTi ~ 1 / Ci * (
      1 / Ria * (Ta - Ti) +
        1 / Rim * (Tm - Ti) +
        Ph +
        Gv * (Aw1 * bs1 + Aw2 * bs2 + Aw3 * bs3 + Aw4 * bs4 + Aw5 * bs5)
    ) * dt + exp(p11) * dw1
  )
  model$addSystem(dTm ~ 1 / Cm * (1 / Rim * (Ti - Tm)) * dt + exp(p22) * dw2)
  
  model$addInput(Ta, Gv, Ph, bs1, bs2, bs3, bs4, bs5)
  
  model$addObs(yTi ~ Ti)
  model$setVariance(yTi ~ exp(e11))
  
  model$setParameter(Ti = c(init = 15, lb = 0, ub = 40))
  model$setParameter(Tm = c(init = 15, lb = 0, ub = 40))
  
  model$setParameter(Ci = c(init = 1, lb = 1e-5, ub = 1e5))
  model$setParameter(Cm = c(init = 1000, lb = 1e-5, ub = 1e5))
  model$setParameter(Ria = c(init = 20, lb = 1e-4, ub = 1e5))
  model$setParameter(Rim = c(init = 20, lb = 1e-4, ub = 1e5))
  
  ub_aw <- 7.5 + 4.8 + 5
  model$setParameter(Aw1 = c(init = 4, lb = 1e-2, ub = ub_aw))
  model$setParameter(Aw2 = c(init = 4, lb = 1e-2, ub = ub_aw))
  model$setParameter(Aw3 = c(init = 4, lb = 1e-2, ub = ub_aw))
  model$setParameter(Aw4 = c(init = 4, lb = 1e-2, ub = ub_aw))
  model$setParameter(Aw5 = c(init = 4, lb = 1e-2, ub = ub_aw))
  
  model$setParameter(p11 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(p22 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(e11 = c(init = -1, lb = -50, ub = 10))
  
  fit <- model$estimate(data, firstorder = TRUE)
  return(fit)
}

splines_neighbor_model <- function(data, yTi, Ph, yTj) {
  data$yTi <- yTi
  data$Ph  <- Ph
  data$yTj <- yTj  # temperature in room j
  
  model <- ctsm()
  
  # ---------------------------------------------------
  # Room 1 Air Temperature with Room j Coupling
  # ---------------------------------------------------
  model$addSystem(
    dTi ~ 1 / Ci * (
      1 / Ria * (Ta - Ti) +
        1 / Rim * (Tm - Ti) +
        1 / Rij * (yTj - Ti) +                  # NEW inter-room term
        Ph +
        Gv * (Aw1 * bs1 + Aw2 * bs2 + Aw3 * bs3 + Aw4 * bs4 + Aw5 * bs5)
    ) * dt + exp(p11) * dw1
  )
  
  # ---------------------------------------------------
  # Room 1 Thermal Mass
  # ---------------------------------------------------
  model$addSystem(
    dTm ~ 1 / Cm * (
      1 / Rim * (Ti - Tm)
    ) * dt + exp(p22) * dw2
  )
  
  # ---------------------------------------------------
  # Inputs
  # ---------------------------------------------------
  model$addInput(Ta, Gv, Ph, yTj, bs1, bs2, bs3, bs4, bs5)
  
  # Observation equation
  model$addObs(yTi ~ Ti)
  model$setVariance(yTi ~ exp(e11))
  
  # ---------------------------------------------------
  # Parameters
  # ---------------------------------------------------
  model$setParameter(Ti = c(init = 24, lb = 0, ub = 40))
  model$setParameter(Tm = c(init = 21, lb = 0, ub = 40))
  
  model$setParameter(Ci = c(init = 8.4, lb = 1e-5, ub = 1e5))
  model$setParameter(Cm = c(init = 50000, lb = 1e-5, ub = 1e5))
  model$setParameter(Ria = c(init = 12, lb = 1e-4, ub = 1e5))
  model$setParameter(Rim = c(init = 0.5, lb = 1e-4, ub = 1e5))
  
  # NEW parameter for inter-room heat flow
  model$setParameter(Rij = c(init = 50, lb = 1e-4, ub = 1e5))
  
  ub_aw <- 7.5 + 4.8 + 5
  model$setParameter(Aw1 = c(init = 17, lb = 1e-2, ub = ub_aw))
  model$setParameter(Aw2 = c(init = 4, lb = 1e-2, ub = ub_aw))
  model$setParameter(Aw3 = c(init = 3, lb = 1e-2, ub = ub_aw))
  model$setParameter(Aw4 = c(init = 3, lb = 1e-2, ub = ub_aw))
  model$setParameter(Aw5 = c(init = 5, lb = 1e-2, ub = ub_aw))
  
  model$setParameter(p11 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(p22 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(e11 = c(init = -1, lb = -50, ub = 10))
  
  # Fit the model
  fit <- model$estimate(data, firstorder = TRUE)
  return(fit)
}

splines_heating_delay_model <- function(data, yTi, Ph) {
  
  data$yTi <- yTi
  data$Ph <- Ph
  
  model <- ctsm()
  
  # --- Radiator dynamic state (delay) ---
  model$addSystem(
    dQh ~ (1/tau_h) * (Ph1 - Qh) * dt + exp(p33) * dw3
  )
  
  # --- Room 1 temperature ---
  model$addSystem(
    dTi ~ 1/Ci * (
      1/Ria*(Ta - Ti) +
        1/Rim*(Tm - Ti)+
        Qh + 
        Gv*(Aw1*bs1 + Aw2*bs2 + Aw3*bs3 + Aw4*bs4 + Aw5*bs5)
    ) * dt + exp(p11)*dw1
  )
  
  # --- Thermal mass Tm ---
  model$addSystem(
    dTm ~ 1/Cm * (
      1/Rim*(Ti - Tm)
    ) * dt + exp(p22)*dw2
  )
  
  # --- Inputs ---
  model$addInput(Ta, Gv, Ph, bs1, bs2, bs3, bs4, bs5)
  
  # --- Observation ---
  model$addObs(yTi ~ Ti)
  model$setVariance(yTi ~ exp(e11))
  
  # --- Parameters ---
  model$setParameter(Ti = c(init = 24, lb = 0, ub = 40))
  model$setParameter(Tm = c(init = 21, lb = 0, ub = 40))
  model$setParameter(Qh = c(init=Ph[1], lb=-1000, ub=1000))
  
  model$setParameter(Ci = c(init = 8.4, lb = 1e-5, ub = 1e5))
  model$setParameter(Cm = c(init = 50000, lb = 1e-5, ub = 1e5))
  model$setParameter(Ria = c(init = 12, lb = 1e-4, ub = 1e5))
  model$setParameter(Rim = c(init = 0.5, lb = 1e-4, ub = 1e5))
  
  # Radiator time constant: 1–4 hours typical
  model$setParameter(tau_h = c(init=2, lb=0.2, ub=24))
  
  # Solar parameters
  ub_aw <- 7.5 + 4.8 + 5
  model$setParameter(Aw1 = c(init = 17, lb = 1e-2, ub = ub_aw))
  model$setParameter(Aw2 = c(init = 4, lb = 1e-2, ub = ub_aw))
  model$setParameter(Aw3 = c(init = 3, lb = 1e-2, ub = ub_aw))
  model$setParameter(Aw4 = c(init = 3, lb = 1e-2, ub = ub_aw))
  model$setParameter(Aw5 = c(init = 5, lb = 1e-2, ub = ub_aw))
  
  # Noise parameters
  model$setParameter(p11 = c(init=1, lb=-20, ub=10))
  model$setParameter(p22 = c(init=1, lb=-20, ub=10))
  model$setParameter(p33 = c(init=-5, lb=-20, ub=5))
  model$setParameter(e11 = c(init=-1, lb=-50, ub=10))
  
  fit <- model$estimate(data, firstorder=TRUE)
  return(fit)
}


splines_solar_delay_model <- function(data, yTi, Ph) {
  
  data$yTi <- yTi
  data$Ph <- Ph
  
  model <- ctsm()
  
  # --- Heating dynamic state (delay) ---
  model$addSystem(
    dSh ~ (1/tau_s) * (Gv - Sh) * dt + exp(p44) * dw4
  )
  
  # --- Room 1 temperature ---
  model$addSystem(
    dTi ~ 1/Ci * (
      1/Ria*(Ta - Ti) +
        1/Rim*(Tm - Ti)+
        Ph + 
        Sh*(Aw1*bs1 + Aw2*bs2 + Aw3*bs3 + Aw4*bs4 + Aw5*bs5)
    ) * dt + exp(p11)*dw1
  )
  
  # --- Thermal mass Tm ---
  model$addSystem(
    dTm ~ 1/Cm * (
      1/Rim*(Ti - Tm)
    ) * dt + exp(p22)*dw2
  )
  
  # --- Inputs ---
  model$addInput(Ta, Gv, Ph, bs1, bs2, bs3, bs4, bs5)
  
  # --- Observation ---
  model$addObs(yTi ~ Ti)
  model$setVariance(yTi ~ exp(e11))
  
  # --- Parameters ---
  model$setParameter(Ti = c(init = 24, lb = 0, ub = 40))
  model$setParameter(Tm = c(init = 21, lb = 0, ub = 40))
  model$setParameter(Sh = c(init=data$Gv[1], lb=0, ub=10))
  
  model$setParameter(Ci = c(init = 8.4, lb = 1e-5, ub = 1e5))
  model$setParameter(Cm = c(init = 50000, lb = 1e-5, ub = 1e5))
  model$setParameter(Ria = c(init = 12, lb = 1e-4, ub = 1e5))
  model$setParameter(Rim = c(init = 0.5, lb = 1e-4, ub = 1e5))
  
  # Radiator time constant: 1–4 hours typical
  model$setParameter(tau_s = c(init=2, lb=0.2, ub=24))
  
  # Solar parameters
  ub_aw <- 7.5 + 4.8 + 5
  model$setParameter(Aw1 = c(init = 17, lb = 1e-2, ub = ub_aw))
  model$setParameter(Aw2 = c(init = 4, lb = 1e-2, ub = ub_aw))
  model$setParameter(Aw3 = c(init = 3, lb = 1e-2, ub = ub_aw))
  model$setParameter(Aw4 = c(init = 3, lb = 1e-2, ub = ub_aw))
  model$setParameter(Aw5 = c(init = 5, lb = 1e-2, ub = ub_aw))
  
  # Noise parameters
  model$setParameter(p11 = c(init=1, lb=-20, ub=10))
  model$setParameter(p22 = c(init=1, lb=-20, ub=10))
  model$setParameter(p44 = c(init=-5, lb=-20, ub=5))
  model$setParameter(e11 = c(init=-1, lb=-50, ub=10))
  
  fit <- model$estimate(data, firstorder=TRUE)
  return(fit)
}

splines_radiatorTemp_model <- function(data, yTi,Ph){
  data$yTi <- yTi
  data$Ph <- Ph

  # Generate a new object of class ctsm
  model = ctsm()
  
  # Add a system equation and thereby also a state
  
  model$addSystem(dTi ~  1/Ci*(1/Ria*(Ta-Ti) + 1/Rim*(Tm-Ti) + 1/Rir*(Tr-Ti) + (a1*bs1+a2*bs2+a3*bs3+a4*bs4+a5*bs5)*Gv)*dt + exp(p11)*dw1)
  model$addSystem(dTr ~ 1/Cr*(1/Rir*(Ti-Tr) + Ph)*dt + exp(prr)*dwr)
  
  model$addSystem(dTm ~  1/Cm*(1/Rim*(Ti-Tm))*dt + exp(p22)*dw2)
  # Set the names of the inputs
  model$addInput(Ta,Gv,Ph,bs1,bs2,bs3,bs4,bs5)
  
  # Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(yTi ~ Ti)
  # Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  
  ##----------------------------------------------------------------
  # Set the initial value (for the optimization) of the value of the state at the starting time point
  model$setParameter(Ti = c(init = 24, lb = 0, ub = 40))
  model$setParameter(Tm = c(init = 21, lb = 0, ub = 40))
  
  model$setParameter(Tr = c(init = 20, lb = 0, ub = 40))
  
  ##----------------------------------------------------------------
  # Set the initial value for the optimization
  
  ub_aw <- 20
  
  model$setParameter(Ci = c(init = 8.4, lb = 1e-5, ub = 1e5))
  model$setParameter(Cm = c(init = 50000, lb = 1e-5, ub = 1e5))
  model$setParameter(Ria = c(init = 12, lb = 1e-4, ub = 1e5))
  model$setParameter(Rim = c(init = 0.5, lb = 1e-4, ub = 1e5))
  
  model$setParameter(p11 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(p22 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(e11 = c(init = -1, lb = -50, ub = 10))
  
  model$setParameter(a1 = c(init = 17, lb = 1e-2, ub = ub_aw))
  model$setParameter(a2 = c(init = 4, lb = 1e-2, ub = ub_aw))
  model$setParameter(a3 = c(init = 3, lb = 1e-2, ub = ub_aw))
  model$setParameter(a4 = c(init = 3, lb = 1e-2, ub = ub_aw))
  model$setParameter(a5 = c(init = 5, lb = 1e-2, ub = ub_aw))
  
  model$setParameter(Cr = c(init = 10, lb = 1E-5, ub = 1E5))
  model$setParameter(Rir = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(prr = c(init = 1, lb = -30, ub = 10))
  
  ##----------------------------------------------------------------   
  # Run the parameter optimization
  
  fit = model$estimate(data,firstorder = TRUE)
  return(fit)
}

combined_spline_delay_neighbor_model <- function(data, yTi, Ph, yTj) {
  
  data$yTi <- yTi
  data$Ph  <- Ph
  data$yTj <- yTj
  
  model <- ctsm()
  
  # =====================================================
  # 1) Heating delay state: Qh(t)
  # =====================================================
  model$addSystem(
    dQh ~ (1/tau_h) * (Ph - Qh) * dt + exp(p33) * dw3
  )
  
  # =====================================================
  # 2) Solar delay state: Sh(t)
  # =====================================================
  model$addSystem(
    dSh ~ (1/tau_s) * (Gv - Sh) * dt + exp(p44) * dw4
  )
  
  # =====================================================
  # 3) Room 1 Air Temperature Ti(t) with all effects
  # =====================================================
  model$addSystem(
    dTi ~ 1 / Ci * (
      1/Ria * (Ta - Ti) +
        1/Rim * (Tm - Ti) +
        1/Rij * (yTj - Ti) +                 # inter-room heat connection
        Qh +                                 # radiator delayed power
        Sh * (Aw1*bs1 + Aw2*bs2 + Aw3*bs3 +
                Aw4*bs4 + Aw5*bs5)             # delayed solar + spline shape
    ) * dt + exp(p11) * dw1
  )
  
  # =====================================================
  # 4) Thermal Mass Tm(t)
  # =====================================================
  model$addSystem(
    dTm ~ 1/Cm * (
      1/Rim * (Ti - Tm)
    ) * dt + exp(p22) * dw2
  )
  
  # =====================================================
  # Inputs
  # =====================================================
  model$addInput(Ta, Gv, Ph, yTj, bs1, bs2, bs3, bs4, bs5)
  
  # =====================================================
  # Observation equation
  # =====================================================
  model$addObs(yTi ~ Ti)
  model$setVariance(yTi ~ exp(e11))
  
  # =====================================================
  # Parameters: states
  # =====================================================
  model$setParameter(Ti  = c(init = 23.6, lb = 0, ub = 40))
  model$setParameter(Tm  = c(init = 21, lb = 0, ub = 40))
  model$setParameter(Qh  = c(init = 12.7,    lb = -500, ub = 500))
  model$setParameter(Sh  = c(init = 1.297, lb = -10, ub = 50))
  
  # =====================================================
  # Physical parameters
  # =====================================================
  model$setParameter(Ci = c(init = 8.4,    lb = 1e-6, ub = 1e3))
  model$setParameter(Cm = c(init = 400,  lb = 1e-6, ub = 1e4))
  model$setParameter(Ria = c(init = 12,    lb = 1e-4, ub = 1e3))
  model$setParameter(Rim = c(init = 0.5,   lb = 1e-4, ub = 1e3))
  
  # Inter-room resistance
  model$setParameter(Rij = c(init = 0.1, lb = 1e-5, ub = 1e3))
  
  # Delay time constants
  model$setParameter(tau_h = c(init = 1.18, lb = 1e-4, ub = 48))
  model$setParameter(tau_s = c(init = 1.34, lb = 1e-4, ub = 48))
  
  # =====================================================
  # Solar spline parameters
  # =====================================================
  ub_aw <- 7.5 + 4.8 + 5
  model$setParameter(Aw1 = c(init = 17, lb = 1e-3, ub = ub_aw))
  model$setParameter(Aw2 = c(init = 4,  lb = 1e-3, ub = ub_aw))
  model$setParameter(Aw3 = c(init = 3,  lb = 1e-3, ub = ub_aw))
  model$setParameter(Aw4 = c(init = 3,  lb = 1e-3, ub = ub_aw))
  model$setParameter(Aw5 = c(init = 5,  lb = 1e-3, ub = ub_aw))
  
  # =====================================================
  # Noise parameters
  # =====================================================
  model$setParameter(p11 = c(init = 1,   lb = -20, ub = 10))
  model$setParameter(p22 = c(init = 1,   lb = -20, ub = 10))
  model$setParameter(p33 = c(init = -5,  lb = -20, ub = 10))
  model$setParameter(p44 = c(init = -5,  lb = -20, ub = 10))
  model$setParameter(e11 = c(init = -1,  lb = -50, ub = 10))
  
  # Fit
  fit <- model$estimate(data, firstorder = TRUE)
  return(fit)
}