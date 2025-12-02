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
  model$setParameter(tau_h = c(init = 1.18, lb = 1e-4, ub = 8))
  model$setParameter(tau_s = c(init = 1.34, lb = 1e-4, ub = 5))
  
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

combined_spline_temp_delay_neighbor_model <- function(data, yTi, Ph, yTj) {
  
  data$yTi <- yTi
  data$Ph  <- Ph
  data$yTj <- yTj
  
  model <- ctsm()
  
  # =====================================================
  # 1) Radiator temperature state: Tr(t)
  # =====================================================
  model$addSystem(
    dTr ~ 1/Cr*(1/Rir*(Ti-Tr) + 1/Rjr*(yTj-Tr) + Ph)*dt + exp(p33) * dw3
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
        1/Rir * (Tr - Ti) +                     # heating from radiator
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
  model$setParameter(Tr  = c(init = 20,    lb = -500, ub = 500))
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
  
  # Room-radiator resistance
  model$setParameter(Rir = c(init = 0.1, lb = 1e-5, ub = 1e3))
  model$setParameter(Rjr = c(init = 0.1, lb = 1e-5, ub = 1e3))
  
  # Radiator heat capacity
  model$setParameter(Cr = c(init = 8.4,    lb = 1e-6, ub = 1e3))
  
  # Delay time constant
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

sde4Rooms_orig <- function(data){
  
  model = ctsm()
  
  ## ------------------------------
  ## Room 1
  ## ------------------------------
  
  model$addSystem(
    dT1 ~ ( 1/C1 * ( 1/Ra1*(Ta - T1)
                     + 1/Rm1*(Tm1 - T1)
                     + Ph1
                     + Aw1*Gv
                     + 1/R12*(T2 - T1)
    ) )*dt
    + exp(p11)*dw1
  )
  
  model$addSystem(
    dTm1 ~ ( 1/Cm1 * ( 1/Rm1*(T1 - Tm1) ) )*dt
    + exp(pm11)*dwm1
  )
  
  ## ------------------------------
  ## Room 2
  ## ------------------------------
  
  model$addSystem(
    dT2 ~ ( 1/C2 * ( 1/Ra2*(Ta - T2)
                     + 1/Rm2*(Tm2 - T2)
                     + Ph1
                     + Aw2*Gv
                     + 1/R12*(T1 - T2)
                     + 1/R23*(T3 - T2)
    ) )*dt
    + exp(p22)*dw2
  )
  
  model$addSystem(
    dTm2 ~ ( 1/Cm2 * ( 1/Rm2*(T2 - Tm2) ) )*dt
    + exp(pm22)*dwm2
  )
  
  ## ------------------------------
  ## Room 3
  ## ------------------------------
  
  model$addSystem(
    dT3 ~ ( 1/C3 * ( 1/Ra3*(Ta - T3)
                     + 1/Rm3*(Tm3 - T3)
                     + Ph2
                     + Aw3*Gv
                     + 1/R23*(T2 - T3)
                     + 1/R34*(T4 - T3)
    ) )*dt
    + exp(p33)*dw3
  )
  
  model$addSystem(
    dTm3 ~ ( 1/Cm3 * ( 1/Rm3*(T3 - Tm3) ) )*dt
    + exp(pm33)*dwm3
  )
  
  ## ------------------------------
  ## Room 4
  ## ------------------------------
  
  model$addSystem(
    dT4 ~ ( 1/C4 * ( 1/Ra4*(Ta - T4)
                     + 1/Rm4*(Tm4 - T4)
                     + Ph2
                     + Aw4*Gv
                     + 1/R34*(T3 - T4)
    ) )*dt
    + exp(p44)*dw4
  )
  
  model$addSystem(
    dTm4 ~ ( 1/Cm4 * ( 1/Rm4*(T4 - Tm4) ) )*dt
    + exp(pm44)*dwm4
  )
  
  ## ------------------------------
  ## Inputs
  ## ------------------------------
  
  model$addInput(Ta, Gv, Ph1, Ph2)
  
  ## ------------------------------
  ## Observations (one per room)
  ## ------------------------------
  
  model$addObs(yTi1 ~ T1); model$setVariance(yTi1 ~ exp(e11))
  model$addObs(yTi2 ~ T2); model$setVariance(yTi2 ~ exp(e22))
  model$addObs(yTi3 ~ T3); model$setVariance(yTi3 ~ exp(e33))
  model$addObs(yTi4 ~ T4); model$setVariance(yTi4 ~ exp(e44))
  
  ## ------------------------------
  ## Parameters
  ## ------------------------------
  
  # initial states
  for(v in c("T1","Tm1","T2","Tm2","T3","Tm3","T4","Tm4")) {
    args <- list()
    args[[v]] <- c(init = 20, lb = -10, ub = 40)
    do.call(model$setParameter, args)
  }
  
  # heat capacities
  model$setParameter(C1  = c(init=1, lb=1e-4, ub=1e5))
  model$setParameter(C2  = c(init=5, lb=1e-4, ub=1e5))
  model$setParameter(C3  = c(init=5, lb=1e-4, ub=1e5))
  model$setParameter(C4  = c(init=1, lb=1e-4, ub=1e5))
  
  model$setParameter(Cm1 = c(init=1, lb=1e-4, ub=1e5))
  model$setParameter(Cm2 = c(init=5, lb=1e-4, ub=1e5))
  model$setParameter(Cm3 = c(init=5, lb=1e-4, ub=1e5))
  model$setParameter(Cm4 = c(init=1, lb=1e-4, ub=1e5))
  
  # resistances
  model$setParameter(Ra1 = c(init=10, lb=1e-2, ub=1e3)); model$setParameter(Ra2 = c(init=10, lb=1e-2, ub=1e3))
  model$setParameter(Ra3 = c(init=10, lb=1e-2, ub=1e3)); model$setParameter(Ra4 = c(init=10, lb=1e-2, ub=1e3))
  
  model$setParameter(Rm1 = c(init=100, lb=1e-2, ub=1e4));  model$setParameter(Rm2 = c(init=100, lb=1e-2, ub=1e4))
  model$setParameter(Rm3 = c(init=100, lb=1e-2, ub=1e4));  model$setParameter(Rm4 = c(init=100, lb=1e-2, ub=1e4))
  
  model$setParameter(R12 = c(init=1, lb=1e-2, ub=1e4))
  model$setParameter(R23 = c(init=10, lb=1e-2, ub=1e4))
  model$setParameter(R34 = c(init=1, lb=1e-2, ub=1e4))
  
  # wind gains
  model$setParameter(Aw1 = c(init = 6, lb = 1E-2, ub = 7.5+4.8+5))
  model$setParameter(Aw2 = c(init = 6, lb = 1E-2, ub = 7.5+4.8+5))
  model$setParameter(Aw3 = c(init = 6, lb = 1E-2, ub = 7.5+4.8+5))
  model$setParameter(Aw4 = c(init = 6, lb = 1E-2, ub = 7.5+4.8+5))
  
  # diffusion parameters
  model$setParameter(p11 = c(init = 1,   lb = -20, ub = 10))
  model$setParameter(p22 = c(init = 1,   lb = -20, ub = 10))
  model$setParameter(p33 = c(init = 1,  lb = -20, ub = 10))
  model$setParameter(p44 = c(init = 1,  lb = -20, ub = 10))
  
  model$setParameter(pm11 = c(init=-1,lb = -20, ub = 10))
  model$setParameter(pm22 = c(init=-1,lb = -20, ub = 10))
  model$setParameter(pm33 = c(init=-1,lb = -20, ub = 10))
  model$setParameter(pm44 = c(init=-1,lb = -20, ub = 10))
  
  # measurement noise
  model$setParameter(e11 = c(init = -1,  lb = -50, ub = 10))
  model$setParameter(e22 = c(init = -1,  lb = -50, ub = 10))
  model$setParameter(e33 = c(init = -1,  lb = -50, ub = 10))
  model$setParameter(e44 = c(init = -1,  lb = -50, ub = 10))
  
  ## ------------------------------
  ## Estimate
  ## ------------------------------
  
  fit <- model$estimate(data, firstorder=TRUE)
  return(fit)
}

sde4Rooms_floorplan <- function(data){
  
  model = ctsm()
  
  # --------------------------------------------------------
  # Room 1
  # --------------------------------------------------------
  model$addSystem(
    dT1 ~  
      1/(flArea1*C) * ( 1/(Ra/exArea1)*(Ta - T1)
                        + 1/(Rm/flArea1)*(Tm1 - T1)
                        + Ph1
                        + Aw1*Gv 
                        + 1/R12*(T2 - T1) ) * dt 
    + exp(p11)*dw1
  )
  model$addSystem(
    dTm1 ~ 1/(flArea1*Cm)*(1/(Rm/flArea1)*(T1 - Tm1))*dt + exp(pm11)*dwm1
  )
  
  # --------------------------------------------------------
  # Room 2
  # --------------------------------------------------------
  model$addSystem(
    dT2 ~  
      1/(flArea2*C) * ( 1/(Ra/exArea2)*(Ta - T2)
                        + 1/(Rm/flArea2)*(Tm2 - T2)
                        + Ph1
                        + Aw2*Gv 
                        + 1/R12*(T1 - T2)
                        + 1/R23*(T3 - T2) ) * dt 
    + exp(p22)*dw2
  )
  model$addSystem(
    dTm2 ~ 1/(flArea2*Cm)*(1/(Rm/flArea2)*(T2 - Tm2))*dt + exp(pm22)*dwm2
  )
  
  # --------------------------------------------------------
  # Room 3
  # --------------------------------------------------------
  model$addSystem(
    dT3 ~  
      1/(flArea3*C) * ( 1/(Ra/exArea3)*(Ta - T3)
                        + 1/(Rm/flArea3)*(Tm3 - T3)
                        + Ph2
                        + Aw3*Gv 
                        + 1/R34*(T4 - T3)
                        + 1/R23*(T2 - T3) ) * dt 
    + exp(p33)*dw3
  )
  model$addSystem(
    dTm3 ~ 1/(flArea3*Cm)*(1/(Rm/flArea3)*(T3 - Tm3))*dt + exp(pm33)*dwm3
  )
  
  # --------------------------------------------------------
  # Room 4
  # --------------------------------------------------------
  model$addSystem(
    dT4 ~  
      1/(flArea4*C) * ( 1/(Ra/exArea4)*(Ta - T4)
                        + 1/(Rm/flArea4)*(Tm4 - T4)
                        + Ph2
                        + Aw4*Gv 
                        + 1/R34*(T3 - T4) ) * dt 
    + exp(p44)*dw4
  )
  model$addSystem(
    dTm4 ~ 1/(flArea4*Cm)*(1/(Rm/flArea4)*(T4 - Tm4))*dt + exp(pm44)*dwm4
  )
  
  # --------------------------------------------------------
  # Inputs
  # --------------------------------------------------------
  model$addInput(Ta, Gv, Ph1, Ph2)
  
  # --------------------------------------------------------
  # Observations
  # --------------------------------------------------------
  model$addObs(yTi1 ~ T1); model$setVariance(yTi1 ~ exp(e11))
  model$addObs(yTi2 ~ T2); model$setVariance(yTi2 ~ exp(e22))
  model$addObs(yTi3 ~ T3); model$setVariance(yTi3 ~ exp(e33))
  model$addObs(yTi4 ~ T4); model$setVariance(yTi4 ~ exp(e44))
  
  # --------------------------------------------------------
  # Initial states
  # --------------------------------------------------------
  for(v in c("T1","Tm1","T2","Tm2","T3","Tm3","T4","Tm4")) {
    args <- list()
    args[[v]] <- c(init = 20, lb = -10, ub = 40)
    do.call(model$setParameter, args)
  }
  
  # --------------------------------------------------------
  # Thermal parameters
  # --------------------------------------------------------
  model$setParameter(C  = c(init=1, lb=1e-4, ub=1e5))
  model$setParameter(Cm = c(init=1, lb=1e-4, ub=1e5))
  
  model$setParameter(Ra = c(init=10, lb=1e-2, ub=1e3))
  model$setParameter(Rm = c(init=100, lb=1e-2, ub=1e4))
  
  model$setParameter(R12 = c(init=1, lb=1e-2, ub=1e4))
  model$setParameter(R23 = c(init=10, lb=1e-2, ub=1e4))
  model$setParameter(R34 = c(init=1, lb=1e-2, ub=1e4))
  
  # --------------------------------------------------------
  # Wind gains
  # --------------------------------------------------------
  model$setParameter(Aw1 = c(init=6, lb=1e-2, ub=7.5+4.8+5))
  model$setParameter(Aw2 = c(init=6, lb=1e-2, ub=7.5+4.8+5))
  model$setParameter(Aw3 = c(init=6, lb=1e-2, ub=7.5+4.8+5))
  model$setParameter(Aw4 = c(init=6, lb=1e-2, ub=7.5+4.8+5))
  
  # --------------------------------------------------------
  # Diffusion parameters
  # --------------------------------------------------------
  model$setParameter(p11 = c(init=1, lb=-20, ub=10))
  model$setParameter(p22 = c(init=1, lb=-20, ub=10))
  model$setParameter(p33 = c(init=1, lb=-20, ub=10))
  model$setParameter(p44 = c(init=1, lb=-20, ub=10))
  
  model$setParameter(pm11 = c(init=-1, lb=-20, ub=10))
  model$setParameter(pm22 = c(init=-1, lb=-20, ub=10))
  model$setParameter(pm33 = c(init=-1, lb=-20, ub=10))
  model$setParameter(pm44 = c(init=-1, lb=-20, ub=10))
  
  # --------------------------------------------------------
  # Measurement noise
  # --------------------------------------------------------
  model$setParameter(e11 = c(init=-1, lb=-50, ub=10))
  model$setParameter(e22 = c(init=-1, lb=-50, ub=10))
  model$setParameter(e33 = c(init=-1, lb=-50, ub=10))
  model$setParameter(e44 = c(init=-1, lb=-50, ub=10))
  
  # --------------------------------------------------------
  # Fixed geometry parameters
  # --------------------------------------------------------
  model$setParameter(flArea1 =c(init=3))
  model$setParameter(flArea2 =c(init=21))
  model$setParameter(flArea3 =c(init=17))
  model$setParameter(flArea4 =c(init=3))
  
  model$setParameter(exArea1 =c(init=1.5))
  model$setParameter(exArea2 =c(init=3))
  model$setParameter(exArea3 =c(init=2.5))
  model$setParameter(exArea4 =c(init=1.5))
  
  # --------------------------------------------------------
  # Estimate
  # --------------------------------------------------------
  fit <- model$estimate(data, firstorder=TRUE)
  return(fit)
}

sde4Rooms_floorplan_heatdelay <- function(data){
  model <- ctsm()
  
  # Heating delay states
  model$addSystem(
    dH1 ~ (1/tau1 * (Ph1 - H1))*dt + exp(ph1_sd)*dwH1
  )
  model$addSystem(
    dH2 ~ (1/tau2 * (Ph2 - H2))*dt + exp(ph2_sd)*dwH2
  )
  
  # Room 1
  model$addSystem(
    dT1 ~ 1/(flArea1*C)*(1/(Ra/exArea1)*(Ta-T1) + 1/(Rm/flArea1)*(Tm1-T1)  
                         + H1
                         + Aw1*Gv 
                         + 1/R12*(T2-T1))*dt 
    + exp(p11)*dw1
  )
  model$addSystem(
    dTm1 ~ 1/(flArea1*Cm)*(1/(Rm/flArea1)*(T1-Tm1))*dt + exp(pm11)*dwm1
  )
  
  # Room 2
  model$addSystem(
    dT2 ~ 1/(flArea2*C)*(1/(Ra/exArea2)*(Ta-T2) + 1/(Rm/flArea2)*(Tm2-T2)  
                         + H1
                         + Aw2*Gv 
                         + 1/R12*(T1-T2) + 1/R23*(T3-T2))*dt 
    + exp(p22)*dw2
  )
  model$addSystem(
    dTm2 ~ 1/(flArea2*Cm)*(1/(Rm/flArea2)*(T2-Tm2))*dt + exp(pm22)*dwm2
  )
  
  # Room 3
  model$addSystem(
    dT3 ~ 1/(flArea3*C)*(1/(Ra/exArea3)*(Ta-T3) + 1/(Rm/flArea3)*(Tm3-T3)  
                         + H2
                         + Aw3*Gv 
                         + 1/R34*(T4-T3) + 1/R23*(T2-T3))*dt 
    + exp(p33)*dw3
  )
  model$addSystem(
    dTm3 ~ 1/(flArea3*Cm)*(1/(Rm/flArea3)*(T3-Tm3))*dt + exp(pm33)*dwm3
  )
  
  # Room 4
  model$addSystem(
    dT4 ~ 1/(flArea4*C)*(1/(Ra/exArea4)*(Ta-T4) + 1/(Rm/flArea4)*(Tm4-T4)  
                         + H2
                         + Aw4*Gv 
                         + 1/R34*(T3-T4))*dt 
    + exp(p44)*dw4
  )
  model$addSystem(
    dTm4 ~ 1/(flArea4*Cm)*(1/(Rm/flArea4)*(T4-Tm4))*dt + exp(pm44)*dwm4
  )
  
  # Inputs
  model$addInput(Ta, Gv, Ph1, Ph2)
  
  # Observations
  model$addObs(yTi1 ~ T1); model$setVariance(yTi1 ~ exp(e11))
  model$addObs(yTi2 ~ T2); model$setVariance(yTi2 ~ exp(e22))
  model$addObs(yTi3 ~ T3); model$setVariance(yTi3 ~ exp(e33))
  model$addObs(yTi4 ~ T4); model$setVariance(yTi4 ~ exp(e44))
  
  # Initial states
  for (v in c("T1","Tm1","T2","Tm2","T3","Tm3","T4","Tm4","H1","H2")) {
    args <- list()
    args[[v]] <- c(init = 20, lb = -10, ub = 40)
    do.call(model$setParameter, args)
  }
  
  # Heat capacities
  model$setParameter(C  = c(init = 1, lb = 1e-4, ub = 1e5))
  model$setParameter(Cm = c(init = 1, lb = 1e-4, ub = 1e5))
  
  # Resistances
  model$setParameter(Ra = c(init = 10,  lb = 1e-2, ub = 1e3))
  model$setParameter(Rm = c(init = 100, lb = 1e-2, ub = 1e4))
  model$setParameter(R12 = c(init = 1,  lb = 1e-2, ub = 1e4))
  model$setParameter(R23 = c(init = 10, lb = 1e-2, ub = 1e4))
  model$setParameter(R34 = c(init = 1,  lb = 1e-2, ub = 1e4))
  
  # Wind gains
  model$setParameter(Aw1 = c(init = 6, lb = 1e-2, ub = 7.5+4.8+5))
  model$setParameter(Aw2 = c(init = 6, lb = 1e-2, ub = 7.5+4.8+5))
  model$setParameter(Aw3 = c(init = 6, lb = 1e-2, ub = 7.5+4.8+5))
  model$setParameter(Aw4 = c(init = 6, lb = 1e-2, ub = 7.5+4.8+5))
  
  # Diffusion parameters (room temps)
  model$setParameter(p11 = c(init = 1, lb = -20, ub = 10))
  model$setParameter(p22 = c(init = 1, lb = -20, ub = 10))
  model$setParameter(p33 = c(init = 1, lb = -20, ub = 10))
  model$setParameter(p44 = c(init = 1, lb = -20, ub = 10))
  
  # Diffusion parameters (thermal masses)
  model$setParameter(pm11 = c(init = -1, lb = -20, ub = 10))
  model$setParameter(pm22 = c(init = -1, lb = -20, ub = 10))
  model$setParameter(pm33 = c(init = -1, lb = -20, ub = 10))
  model$setParameter(pm44 = c(init = -1, lb = -20, ub = 10))
  
  # Heating delay parameters
  model$setParameter(tau1   = c(init = 1, lb = 0, ub = 5))
  model$setParameter(tau2   = c(init = 1, lb = 0, ub = 5))
  model$setParameter(ph1_sd = c(init = -10,   lb = -20, ub = 10))
  model$setParameter(ph2_sd = c(init = -10,   lb = -20, ub = 10))
  
  # Measurement noise
  model$setParameter(e11 = c(init = -1, lb = -50, ub = 10))
  model$setParameter(e22 = c(init = -1, lb = -50, ub = 10))
  model$setParameter(e33 = c(init = -1, lb = -50, ub = 10))
  model$setParameter(e44 = c(init = -1, lb = -50, ub = 10))
  
  # Fixed parameters
  model$setParameter(flArea1 = c(init = 3))
  model$setParameter(flArea2 = c(init = 21))
  model$setParameter(flArea3 = c(init = 17))
  model$setParameter(flArea4 = c(init = 3))
  model$setParameter(exArea1 = c(init = 1.5))
  model$setParameter(exArea2 = c(init = 3))
  model$setParameter(exArea3 = c(init = 2.5))
  model$setParameter(exArea4 = c(init = 1.5))
  
  fit <- model$estimate(data, firstorder = TRUE)
  return(fit)
}

sde4Rooms_orig_heatdelay <- function(data){
  
  model = ctsm()
  
  ## ---------------------------------------------------------
  ## Heating delay states (north & south circuits)
  ## ---------------------------------------------------------
  
  model$addSystem(
    dH1 ~ (1/tau1 * (Ph1 - H1)) * dt + exp(pH1) * dwH1
  )
  
  model$addSystem(
    dH2 ~ (1/tau2 * (Ph2 - H2)) * dt + exp(pH2) * dwH2
  )
  
  ## ---------------------------------------------------------
  ## Room 1 (north circuit → H1)
  ## ---------------------------------------------------------
  
  model$addSystem(
    dT1 ~ ( 1/C1 * ( 1/Ra1*(Ta - T1)
                     + 1/Rm1*(Tm1 - T1)
                     + H1
                     + Aw1*Gv
                     + 1/R12*(T2 - T1)
    ) )*dt
    + exp(p11)*dw1
  )
  
  model$addSystem(
    dTm1 ~ ( 1/Cm1 * ( 1/Rm1*(T1 - Tm1) ) )*dt
    + exp(pm11)*dwm1
  )
  
  ## ---------------------------------------------------------
  ## Room 2 (north circuit → H1)
  ## ---------------------------------------------------------
  
  model$addSystem(
    dT2 ~ ( 1/C2 * ( 1/Ra2*(Ta - T2)
                     + 1/Rm2*(Tm2 - T2)
                     + H1
                     + Aw2*Gv
                     + 1/R12*(T1 - T2)
                     + 1/R23*(T3 - T2)
    ) )*dt
    + exp(p22)*dw2
  )
  
  model$addSystem(
    dTm2 ~ ( 1/Cm2 * ( 1/Rm2*(T2 - Tm2) ) )*dt
    + exp(pm22)*dwm2
  )
  
  ## ---------------------------------------------------------
  ## Room 3 (south circuit → H2)
  ## ---------------------------------------------------------
  
  model$addSystem(
    dT3 ~ ( 1/C3 * ( 1/Ra3*(Ta - T3)
                     + 1/Rm3*(Tm3 - T3)
                     + H2
                     + Aw3*Gv
                     + 1/R23*(T2 - T3)
                     + 1/R34*(T4 - T3)
    ) )*dt
    + exp(p33)*dw3
  )
  
  model$addSystem(
    dTm3 ~ ( 1/Cm3 * ( 1/Rm3*(T3 - Tm3) ) )*dt
    + exp(pm33)*dwm3
  )
  
  ## ---------------------------------------------------------
  ## Room 4 (south circuit → H2)
  ## ---------------------------------------------------------
  
  model$addSystem(
    dT4 ~ ( 1/C4 * ( 1/Ra4*(Ta - T4)
                     + 1/Rm4*(Tm4 - T4)
                     + H2
                     + Aw4*Gv
                     + 1/R34*(T3 - T4)
    ) )*dt
    + exp(p44)*dw4
  )
  
  model$addSystem(
    dTm4 ~ ( 1/Cm4 * ( 1/Rm4*(T4 - Tm4) ) )*dt
    + exp(pm44)*dwm4
  )
  
  ## ---------------------------------------------------------
  ## Inputs
  ## ---------------------------------------------------------
  
  model$addInput(Ta, Gv, Ph1, Ph2)
  
  ## ---------------------------------------------------------
  ## Observations
  ## ---------------------------------------------------------
  
  model$addObs(yTi1 ~ T1); model$setVariance(yTi1 ~ exp(e11))
  model$addObs(yTi2 ~ T2); model$setVariance(yTi2 ~ exp(e22))
  model$addObs(yTi3 ~ T3); model$setVariance(yTi3 ~ exp(e33))
  model$addObs(yTi4 ~ T4); model$setVariance(yTi4 ~ exp(e44))
  
  ## ---------------------------------------------------------
  ## Parameters
  ## ---------------------------------------------------------
  
  # initial states
  for(v in c("T1","Tm1","T2","Tm2","T3","Tm3","T4","Tm4","H1","H2")) {
    args <- list()
    args[[v]] <- c(init = 20, lb = -20, ub = 40)
    do.call(model$setParameter, args)
  }
  
  # heat capacities
  model$setParameter(C1  = c(init=1, lb=1e-4, ub=1e5))
  model$setParameter(C2  = c(init=5, lb=1e-4, ub=1e5))
  model$setParameter(C3  = c(init=5, lb=1e-4, ub=1e5))
  model$setParameter(C4  = c(init=1, lb=1e-4, ub=1e5))
  
  model$setParameter(Cm1 = c(init=1, lb=1e-4, ub=1e5))
  model$setParameter(Cm2 = c(init=5, lb=1e-4, ub=1e5))
  model$setParameter(Cm3 = c(init=5, lb=1e-4, ub=1e5))
  model$setParameter(Cm4 = c(init=1, lb=1e-4, ub=1e5))
  
  # resistances
  model$setParameter(Ra1 = c(init=10, lb=1e-2, ub=1e3)); model$setParameter(Ra2 = c(init=10, lb=1e-2, ub=1e3))
  model$setParameter(Ra3 = c(init=10, lb=1e-2, ub=1e3)); model$setParameter(Ra4 = c(init=10, lb=1e-2, ub=1e3))
  
  model$setParameter(Rm1 = c(init=100, lb=1e-2, ub=1e4)); model$setParameter(Rm2 = c(init=100, lb=1e-2, ub=1e4))
  model$setParameter(Rm3 = c(init=100, lb=1e-2, ub=1e4)); model$setParameter(Rm4 = c(init=100, lb=1e-2, ub=1e4))
  
  model$setParameter(R12 = c(init=1, lb=1e-2, ub=1e4))
  model$setParameter(R23 = c(init=10, lb=1e-2, ub=1e4))
  model$setParameter(R34 = c(init=1, lb=1e-2, ub=1e4))
  
  # wind gains
  model$setParameter(Aw1 = c(init=6, lb=1e-2, ub=7.5+4.8+5))
  model$setParameter(Aw2 = c(init=6, lb=1e-2, ub=7.5+4.8+5))
  model$setParameter(Aw3 = c(init=6, lb=1e-2, ub=7.5+4.8+5))
  model$setParameter(Aw4 = c(init=6, lb=1e-2, ub=7.5+4.8+5))
  
  # diffusion parameters (rooms)
  model$setParameter(p11 = c(init = 1, lb = -20, ub = 10))
  model$setParameter(p22 = c(init = 1, lb = -20, ub = 10))
  model$setParameter(p33 = c(init = 1, lb = -20, ub = 10))
  model$setParameter(p44 = c(init = 1, lb = -20, ub = 10))
  
  # diffusion parameters (middle masses)
  model$setParameter(pm11 = c(init=-1, lb=-20, ub=10))
  model$setParameter(pm22 = c(init=-1, lb=-20, ub=10))
  model$setParameter(pm33 = c(init=-1, lb=-20, ub=10))
  model$setParameter(pm44 = c(init=-1, lb=-20, ub=10))
  
  # diffusion parameters (heating delay)
  model$setParameter(pH1 = c(init=-10, lb=-20, ub=10))
  model$setParameter(pH2 = c(init=-10, lb=-20, ub=10))
  
  # heating delay time constants
  model$setParameter(tau1 = c(init = 1, lb = 0, ub = 10))
  model$setParameter(tau2 = c(init = 1, lb = 0, ub = 10))
  
  # measurement noise
  model$setParameter(e11 = c(init=-1, lb=-50, ub=10))
  model$setParameter(e22 = c(init=-1, lb=-50, ub=10))
  model$setParameter(e33 = c(init=-1, lb=-50, ub=10))
  model$setParameter(e44 = c(init=-1, lb=-50, ub=10))
  
  ## ---------------------------------------------------------
  ## Fit model
  ## ---------------------------------------------------------
  
  fit <- model$estimate(data, firstorder=TRUE)
  return(fit)
}