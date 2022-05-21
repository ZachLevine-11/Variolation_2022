##Compute R0 based on the calculation in simplemaskmodel.pdf.
compute_R0 <- function(params = make_params(), use_unequal_transmission = TRUE){
  if (!use_unequal_transmission){
    R0 <- as.numeric(params["Beta"]*(params["gamma1"]*(1 - params["rho"]) + params["gamma2"] * params["rho"] + params["mu"])/((params["gamma1"] + params["mu"])*(params["gamma2"] + params["mu"])))
  }
  else{
    ##Allow to solve for R0 using a reduction in transmission rate.
    R0 <- as.numeric(params["Beta"]*((params["mu"] + params["gamma2"])*params["rho"]*params["bm"] + (1- params["rho"])*(params["mu"] + params["gamma1"]))/((params["gamma1"] + params["mu"])*(params["gamma2"] + params["mu"])))
  }
  return(R0)
}

##Compute the expression for the EE based on the params and the model selection.
##Set use_unequal_transsmission to TRUE to access that model
##To use the same transmission rates for Im and IS, just set gamma2 = gamma1 in the params
compute_EE <- function(params = make_params(), use_unequal_transmission = TRUE, withP = FALSE){
  ##S is the same regardless of the rate of vaccination, p
  S <- 1/compute_R0(params = params, use_unequal_transmission = use_unequal_transmission)
  ##Easier to unpack into the environment for all the other ones
  with(as.list(c(params)), {
    if (!use_unequal_transmission){
      bm <- 1
    }
    else{
    }
    if (!withP){
      ##Without vaccination rates
      Im <- ((gamma2 + mu)* (delta + mu)* rho*((-1)*((gamma1 + mu)*(gamma2 + mu)) + N*Beta*(gamma1 + mu - (gamma1 + mu - bm * (gamma2 + mu)) * rho)))/(Beta*((gamma1 + mu)* (gamma2 + delta + mu) + (-gamma1 + gamma2)*delta*rho)*(gamma1 + mu - (gamma1 + mu - bm *(gamma2 + mu)) * rho))
      Is <- ((gamma1 + mu)* (delta + mu)*(rho - 1)*(-N*Beta*(gamma1 + mu) + (gamma1 + mu)*(gamma2 + mu) + N*Beta*(gamma1 + mu - bm*(gamma2 + mu)) * rho))/(Beta*((gamma1 + mu)* (gamma2 + delta + mu) + (-gamma1 + gamma2)*delta*rho)*(gamma1 + mu - (gamma1 + mu - bm *(gamma2 + mu)) * rho))
    }
    else{
      ##With vaccination rates
      Im <- -(((gamma2 + mu)*rho*((gamma1 + mu)*(gamma2 + mu)*(delta + mu) -N*Beta*(delta + mu - p*mu)*(gamma1 + mu - (gamma1 + mu - bm*(gamma2 + mu)) * rho)))/(Beta * (( gamma1 + mu)*(gamma2 + delta + mu) + (-gamma1 + gamma2)* delta* rho)*(gamma1 + mu - (gamma1 + mu - bm*(gamma2 + mu))*rho)))
      Is <- ((gamma1 + mu)*(rho - 1)*((gamma1 + mu)*(gamma2 + mu)*(delta + mu) - N * Beta*(delta + mu - p * mu)*(gamma1 + mu -(gamma1 + mu - bm*(gamma2 + mu))*rho))) /(Beta * (( gamma1 + mu)*(gamma2 + delta + mu) + (-gamma1 + gamma2)* delta* rho)*(gamma1 + mu - (gamma1 + mu - bm*(gamma2 + mu))*rho))
    }
    ##R is the same either way
    R <- as.numeric(params["N"] - S - Im - Is)
    EE <- c("S" = S,
            "Im" = Im,
            "Is" = Is,
            "R" = R)
    return(EE)
  })
}

##Compute gamma1 based on a choice of parameters and R0.
compute_gamma1 <- function(params = make_params(), R0){
  gamma1 <- as.numeric((params["gamma2"]*params["rho"])/((R0*(params["gamma2"] + params["mu"])/params["Beta"]) + params["rho"] - 1))
  return(as.numeric(gamma1))
}
##Derive gamma2 based on a choice of parameters and R0.
compute_gamma2 <- function(params = make_params(), R0){
  gamma2 <- as.numeric((params["gamma1"]*(1-params["rho"]))/((R0*(params["gamma1"] + params["mu"])/params["Beta"]) - params["rho"]))
  return(as.numeric(gamma2))
}
##Derive rho from R0 and other parameters.
##Now updated to use unequal transmimssion
compute_rho <- function(R0, params = make_params(), use_unequal_transmission = TRUE){
  if (use_unequal_transmission == TRUE){
    with(as.list(c(params)), {
      rho <- (R0 - Beta/(gamma2 + mu))/(Beta*bm/(gamma1 + mu) - Beta/(gamma2 + mu))
      return(rho)
    })
  }
  else{
    print("This has only been used for the case with unequal transmission. Please try setting bm = 1 and try again")
    return("NA")
  }
}

##Derive Beta rho from R0 and other parameters.
compute_beta <- function(R0, params, use_unequal_transmission = TRUE){
  if (!use_unequal_transmission){
    Beta <- R0*((params["gamma1"] + params["mu"])*(params["gamma2"] + params["mu"]))/((params["gamma1"]*(1 - params["rho"])) + params["gamma2"]*params["rho"] + params["mu"])
  }
  else{
    ##Use with a reduction in transmission rate, grabbing the tranmsission reduction from variolation from the params.
    Beta <- R0*((params["gamma1"] + params["mu"])*(params["gamma2"] + params["mu"]))/((params["mu"] + params["gamma2"])*params["rho"]*params["bm"] + (1- params["rho"])*(params["mu"] + params["gamma1"]))
  }
  return(as.numeric(Beta))
}

###Package params to pass to the simulation.
##Include defaults.
##You can override any default selection by passing the new value on its own. 
##Because we're trying to investigate final size, vital dynamics are by default set to zero.
make_params <- function(Beta = 0.3287973, ##Sets R0 to 3 given the rest of the parameters, set using compute_beta
                        rho = 0.6,
                        ##Define based on the generation interval
                        gamma1 = 1/(3.7+9),
                        gamma2 = 1/(3.7+14),
                        delta = 0,
                        bm = 0.5,
                        nu = 0.0105,
                        mu = nu,
                        N = 1,
                        p = 0.5,
                        ##For EE things
                        ##Pass this to overwrite any passed values with those from a built in parameter set, except for those in change_pars.
                        base_params = NULL,
                        change_pars = NULL
){
  
  if (is.null(base_params)){
    return(c("Beta" = Beta,
             "rho" = rho,
             "gamma1" = gamma1,
             "gamma2" = gamma2,
             "delta" = delta,
             ##We don't need bm unless we're considering unequal transmission rates.
             "bm" = bm,
             "nu" = nu,
             "mu" = mu,
             "N" = N,
             "p" = p))
  }
  else {
    library("dplyr")
    i <- 1
    while (i <= length(change_pars)){
      change_par <- change_pars[i]
      base_params[change_par] <- eval(parse(text = paste0(change_par)))
      i <- i + 1
    }
    return(base_params)
  }
}

##Make the reduced parameter set for the simple SIR model.
make_simple_SIR_params <- function(Beta =  0.2,
                                   gamma = (make_params()[["gamma1"]] +  make_params()[["gamma2"]])/2,
                                   delta = 0,
                                   N = 0,
                                   mu = 0,
                                   nu = 0){
  return(c("Beta" = Beta,
           "gamma" = gamma,
           "delta" = delta,
           "nu" = nu,
           "mu" = mu,
           "N" = N))
  
}

make_sir_inits <- function(S = 0.9999,
                           I = 0.0001,
                           R = 0){
  return(c(
    "S" = S,
    "I" = I,
    "R" = R)
  )
  
}

###Package initial values to pass to the simulation, including defaults and allow for overriding in the same way as above.
##Set useExpandedModel to TRUE to pass initial values for the state variables of the expanded system too.
make_inits <- function(S = 0.9999,
                       I1 = 0.00008,
                       I2 = 0.00002,
                       R = 0,
                       ##Allow for
                       useExpandedModel = FALSE,
                       R1 = 0,
                       R2 = 0){
  
  if (useExpandedModel){
    ##Return initial values for the state variables of the expanded system for final size explorations.
    return(c(
      "S" = S,
      "I1" = I1,
      "I2" = I2,
      "R1" = R1,
      "R2" = R2)
    )
  }
  else{
    return(c(
      "S" = S,
      "I1" = I1,
      "I2" = I2,
      "R" = R)
    )
  }
}
##Solve the differential equation model.
##Set compute_final_size to TRUE to simulate over a long time period and return the final size of the epidemic (R(infinity)) only.
##Set useExpandedModel = TRUE to use the expanded differential equation model instead of the standard one.
simulate_solve <- function(params = make_params(),
                           useExpandedModel = FALSE,
                           init = make_inits(useExpandedModel = useExpandedModel),
                           length = 360,
                           ##for any reopening or time varying Beta, the time step needs to be set to one.
                           by = 1,
                           compute_final_size = FALSE,
                           use_unequal_transmission = TRUE,
                           use_simple_SIR = FALSE,
                           ##General option
                           use_timevarying_Beta = FALSE,
                           ##Specific for opening plot, pass if you want to override
                           ##what date to simulate reopening?
                           opening_time_index = 100,
                           ##double contact on that date.
                           opening_beta_multiplier = 2,
                           ##Beta_vector should be a vector with the ith value being the ith value in the time series for beta
                           ##Needs to be off by one ( - 2 so that the Beta changes on the opening index - 1th date because we have to add one to avoid indexing by one)
                           Beta_vector = c(rep(params[["Beta"]], opening_time_index - 2), opening_beta_multiplier*params[["Beta"]], rep(opening_beta_multiplier*params[["Beta"]], length - opening_time_index + 2))){
  library("deSolve")
  library("tibble")
  ## Create an SIR function
  if (useExpandedModel & !use_simple_SIR & !use_unequal_transmission){
    ##Use the model with expanded S and R compartments S1, S2, R1, and R2.
    theModel <- function(time, state, parameters) {
      with(as.list(c(state, parameters)), {
        dS <-  nu * N - Beta*(I1 + I2) * S - mu * S + delta * R1
        dI1 <- rho * Beta*(I1 + I2) * S - gamma1*I1 - mu * I1
        dI2 <- (1 - rho) * Beta*(I1 + I2) * S - gamma2 * I2 - mu * I2
        dR1 <- gamma1 * I1 - mu * R1 - delta * R1
        dR2 <- gamma2 * I2 - mu * R2 - delta * R2
        return(list(c(dS, dI1, dI2, dR1, dR2)))
      })
    }
  }
  else if (!useExpandedModel & !use_simple_SIR & !use_unequal_transmission){
    ##Use the basic differential equations model.
    theModel <- function(time, state, parameters) {
      with(as.list(c(state, parameters)), {
        dS <-  nu * N - Beta*(I1 + I2) * S - mu * S + delta * R
        dI1 <- rho * Beta*(I1 + I2) * S - gamma1*I1 - mu * I1
        dI2 <- (1 - rho) * Beta*(I1 + I2) * S - gamma2 * I2 - mu * I2
        dR <- gamma1 * I1 + gamma2 * I2 - mu * R - delta * R
        return(list(c(dS, dI1, dI2, dR)))
      })
    }
  }
  else if (!use_unequal_transmission){
    ##Just solve the equations for the simple SIR model.
    ##Use the basic differential equations model.
    theModel <- function(time, state, parameters) {
      with(as.list(c(state, parameters)), {
        dS <-  nu*N - Beta*(I)*S - mu*S + delta*R
        dI <- Beta*I*S - gamma*I - mu*I
        dR <- gamma * I - delta * R
        return(list(c(dS, dI, dR)))
      })
    }
  }
  else if (use_unequal_transmission & !useExpandedModel){
    ##Use unequal transmission rates.
    ##use_unequal_transmission = TRUE in this case
    ##But don't use the expanded model.
    theModel <- function(time, state, parameters) {
      with(as.list(c(state, parameters)), {
        if (use_timevarying_Beta){
          ##There is no 0th element in the Beta_vector, so index by one
          Beta <- Beta_vector[[time  + 1]]
        }
        else{
          ##Beta is as defined in the params
        }
        dS <-  nu * N - Beta*(I1 * bm + I2) * S - mu * S + delta * R
        dI1 <- rho * Beta*(I1 * bm + I2) * S - gamma1*I1 - mu * I1
        dI2 <- (1 - rho) * Beta*(I1 * bm + I2) * S - gamma2 * I2 - mu * I2
        dR <- gamma1 * I1 + gamma2 * I2 - mu * R - delta * R
        return(list(c(dS, dI1, dI2, dR)))
      })
    }
  }
  else{
    ##Now use the expanded model.
    ##Use the model with expanded S and R compartments S1, S2, R1, and R2.
    theModel <- function(time, state, parameters) {
      with(as.list(c(state, parameters)), {
        dS <-  nu * N - Beta*(I1*bm + I2) * S - mu * S + delta * R1
        dI1 <- rho * Beta*(I1*bm + I2) * S - gamma1*I1 - mu * I1
        dI2 <- (1 - rho) * Beta*(I1*bm + I2) * S - gamma2 * I2 - mu * I2
        dR1 <- gamma1 * I1 - mu * R1 - delta * R1
        dR2 <- gamma2 * I2 - mu * R2 - delta * R2
        return(list(c(dS, dI1, dI2, dR1, dR2)))
      })
    }
  }
  if (compute_final_size){
    time <- seq(0, length*5, by = 1)
    if (useExpandedModel & !use_simple_SIR){
      ##Return Z1 and Z2 separately.
      Z <- c("Z1" = as.data.frame(ode(y = init , times = time, func = theModel, parms = params))[[length*5, "R1"]],
             "Z2" =  as.data.frame(ode(y = init , times = time,func = theModel, parms = params))[[length*5, "R2"]])
    }
    else if (!useExpandedModel & !use_simple_SIR){
      ##There is no Z1 or Z2, so just return Z.
      Z <- as.data.frame(ode(y = init,
                             times = time,
                             func = theModel,
                             parms = params))[[length*5, "R"]]
    }
    else{
      ##Just solve the simple SIR model, updating the Beta to match the R0 of the params.
      Z <- as.data.frame(ode(y = make_sir_inits(),
                             times = time,
                             func = theModel,
                             parms = make_simple_SIR_params(Beta = params[["Beta"]])))[[length*5, "R"]]
    }
    return(Z)
  }
  else{
    ##Return the whole simulation.
    time <- seq(0, length, by = by)
    if (!use_simple_SIR){
      solved <- as.data.frame(ode(y = init ,
                                  times = time,
                                  func = theModel,
                                  parms = params))
    }
    else{
      ##The simple SIR model is the only one with a different call.
      solved <- as.data.frame(ode(y = make_sir_inits() ,
                                  times = time,
                                  func = theModel,
                                  parms = make_simple_SIR_params()))
    }
    
    return(solved)
  }
}

all_curve_sim_palette <- c("#F28C8E", "#377EB8", "#E41A1C" ,"#4DAF4A")
three_curve_sim_palette <- c("#F28C8E", "#E41A1C","#000000")
recovered_palette <- c("#9BBFDC", "#377EB8")
plot_sim <- function(use_unequal_transmission = TRUE,
                     params = make_params(),
                     twoPanel = FALSE,
                     textSize = 12,
                     end =  200,
                     use_timevarying_Beta = FALSE,
                     ##if using two panels, we'll need the expanded model.
                     theSol = simulate_solve(params = params, use_unequal_transmission = use_unequal_transmission, useExpandedModel = twoPanel, use_timevarying_Beta = use_timevarying_Beta, opening_time_index = opening_time_index, opening_beta_multiplier = opening_beta_multiplier),
                     ##Drop only matters for the single (not two) panel plot
                     drop = c(),
                     show_sum = FALSE,
                     opening_time_index = 100,
                     ##double contact on that date.
                     opening_beta_multiplier = 2){
  library("dplyr")
  library("ggpubr")
  library("ggplot2")
  library("tidyr")
  library("directlabels")
  library("tibble")
  if (!twoPanel){
    theSol <- as_tibble(theSol)
    theSol$"Mildly infectious (variolated)" <- theSol$I1
    theSol$"Severely infectious" <- theSol$I2
    theSol$Susceptible <- theSol$S
    theSol$Recovered <- theSol$R
    theSol <- theSol[,!colnames(theSol) %in% c("I1", "I2", "S", "R")]
    theSol <- theSol[,!colnames(theSol) %in% drop]
    ##clip the simulations to fit the labels in.
    ##Assumes a by of 1
    # theSol <- theSol[1:(nrow(theSol) - 50),]
    if (length(drop) == 2){
      ##Not using the expended model
      pal <- three_curve_sim_palette
    }
    else{
      pal <- all_curve_sim_palette
    }
    if (show_sum){
      ##Add a curve for the sum of I1 and I2
      theSol$Total <- theSol$"Severely infectious" + theSol$"Mildly infectious (variolated)"
    }
    else{
    }
    ncurves <- length(theSol)
    theSol <- tidyr::pivot_longer(theSol, col = !time)
    p <- ggline(size = 1.5,plot_type = "l",data = theSol, x= "time", y = "value", color = "name", palette = pal) +
      labs(title = "Prevalence",
           x = "Time (days)",
           y = "Proportion of population") + ggplot2::coord_cartesian(clip  = "off")
    # p <- direct.label(p, method = list(cex = 1 ,"top.points"))
    p <- p + ggplot2::scale_x_discrete(expand = c(0,0), breaks = seq(from = 0, to = end, by = 25)) 
    return(p)
  }
  else{
    ##Make the two panel plot
    ##Since we're playing around with final size, set all the vital dynamics to zero
    params["nu"] <- params["mu"] <- 0
    theSol <- as_tibble(theSol)
    theSol$"Mildly infectious (variolated)" <- theSol$I1
    theSol$"Severely infectious" <- theSol$I2
    theSol$"$\\Rmild$" <- theSol$R1
    theSol$"$\\Rsevere$" <- theSol$R2
    theSol_prevalence <- theSol[,colnames(theSol) %in% c("Mildly infectious (variolated)", "Severely infectious", "time")]
    theSol_prevalence <- theSol_prevalence[1:end,]
    theSol_cumulative <- theSol[1:(end - 12),]
    theSol_cumulative <-  theSol_cumulative[,colnames(theSol) %in% c("$\\Rmild$", "$\\Rsevere$", "time")]
    if (show_sum){
      ##Add a curve for the sum of I1 and I2
      theSol$Total <- theSol$"Severely infectious" + theSol$"Mildly infectious (variolated)"
    }
    else{
    }
    ncurves <- length(theSol)
    theSol_prevalence <- tidyr::pivot_longer(theSol_prevalence, col = !time)
    ##so we don't pass the labels
    theSol_cumulative <- tidyr::pivot_longer(theSol_cumulative, col = !time)[0:(2*end - 10),]
    library("patchwork")
    Z <- sum(simulate_solve(params = params, use_unequal_transmission = use_unequal_transmission, useExpandedModel = twoPanel, compute_final_size = TRUE))
    p1 <- ggplot(data = theSol_prevalence, mapping = aes(x= time, y = value, color = name)) + geom_line(lwd = 1.5) + scale_color_manual(values = three_curve_sim_palette) +
      labs(title = "prevalence time series",
           x = "time (days)",
           y = "proportion of population") +
      theme_pubr()
    #   p1 <- direct.label(p1, method = list(cex = 1.3 ,"top.points"))
    p1 <- p1 + geom_text(aes(100, 0.16, label = "Mildly infectious (variolated)", vjust = 1), color =  three_curve_sim_palette[1], size = textSize/3) +
      ##The labels are centered at the x coordinate
      geom_text(aes(125, 0.05, label = "Severely Infectious", vjust = 0.5), color = three_curve_sim_palette[2], size  = textSize/3) +  theme(legend.position = "none")
    p2 <- ggplot(data = theSol_cumulative, mapping = aes(x = time, y = value, color = name)) + geom_line(lwd = 1.5) +
      # geom_segment(x = 0, y = params["rho"]*Z, xend = end, yend = params["rho"]*Z, color = recovered_palette[1], lwd = 1.5) +
      ##Subtract 2 from the xend of this line so that the label for (1-m)Z, which is longer than the one above it, doesn't hit the line
      #  geom_segment(x = 0, y = (1-params["rho"])*Z, xend = end -2, yend = (1-params["rho"])*Z,  color = recovered_palette[2], lwd = 1.5) +
      ##Labels
      geom_text(aes(end - 5, params["rho"]*Z, label = "$\\mildprop Z$"), color =  recovered_palette[1], size = textSize/3) +
      geom_text(aes(end - 5, (1-params["rho"])*Z, label = "$(1-\\mildprop) Z$"), color = recovered_palette[2], size  = textSize/3) +
      scale_color_manual(values = recovered_palette) +
      labs(title = "cumulative incidence",
           x = "time (days)",
           y = "proportion of population") +
      theme_pubr() +
      theme(legend.position = "none") + 
      geom_text(aes(88, 0.45, label = "$\\Rmild$", vjust = 1), color =  recovered_palette[1], size = textSize/3) +
      geom_text(aes(80, 0.2, label = "$\\Rsevere$", vjust = 0.5), color = recovered_palette[2], size  = textSize/3)
    #    p2 <- direct.label(p2, method = list(cex = 1.3 ,"smart.grid"))
    p1 <- p1 + ggplot2::scale_x_continuous(expand = c(0,0), breaks = seq(from = 0, to = end, by = 20), limits = c(0, end)) + ggplot2::theme(text = element_text(size=textSize)) + ggplot2::coord_cartesian(clip="off") 
    ##Do this first in order to not obliterate the scale command underneath
    p2 <- p2 + ggplot2::scale_x_continuous(expand = c(0,0), limits = c(0, end), breaks = seq(from = 0, to = end, by = 20)) + ggplot2::theme(text = element_text(size=textSize))  + ggplot2::coord_cartesian(clip="off") 
    return(p1 + p2 + plot_layout(ncol = 1))
  }
}

reopening_plot <- function(use_unequal_transmission = TRUE,
                           params = make_params(),
                           ##Assumes use_timevarying_Beta = TRUE
                           opening_time_index = 100,
                           use_directlabels = TRUE,
                           ##double contact on that date.
                           opening_beta_multiplier = 2,
                           time = seq(from = 0, to = 360, by = 1),
                           rhoSeq = seq(from = 0, to = 1, by = 0.25)){
  library("ggplot2")
  library("ggpubr")
  library("tibble")
  df <- tibble("$\\Isevere$" = c(), "Time" = c(), "$\\mildprop$" = c())
  for (rhoVal in rhoSeq){
    params["rho"] <- rhoVal
    severe_prevalence_time_series <- simulate_solve(params = params, use_unequal_transmission = use_unequal_transmission, useExpandedModel = FALSE, use_timevarying_Beta = TRUE, opening_time_index = opening_time_index, opening_beta_multiplier = opening_beta_multiplier)$I2
    tempdf <- tibble("$\\Isevere$" = severe_prevalence_time_series, "Time" = time, "$\\mildprop$" = c(rhoVal))
    df <- dplyr::bind_rows(df, tempdf)
  }
  p <- ggplot2::ggplot(data = df, mapping = aes(x = Time, y = `$\\Isevere$`, color = as.factor(`$\\mildprop$`))) +
    geom_line(lwd = 1.5) +
    ##Reopening date line
    geom_vline(xintercept = opening_time_index, lwd = 1.2, linetype = "dotted") + 
    geom_text(mapping = aes(x = opening_time_index + 40, y = 0.4), label = "Reopening date", color =  "black")
  if (use_directlabels){
    p <- p +    directlabels::geom_dl(mapping = aes(label = `$\\mildprop$`), method = "last.bumpup") 
  }
  p <- p + theme_pubr() +
    theme(legend.position = "none")
  p <-  p + labs(title = "Reopening and masks", y = "Severe prevalence", x = "Time (days)")
  p
}


##Estimate the final size relationship of a pandemic.
##Only works for the original, not expanded set of differential equations
##For a sequence of R0 values, back compute Beta fixing the other parameters, solve the differential equations of the model numerically, and run the simulation for a long period of time to estimate the final size.
##For each set offin parameters, find R0 and then run a simulation and find the final size.
##Currently only updating Beta to match R0 is supported.
##Set checkFormula to test different final size formulas.
long_run_simulator <- function(R0by = 0.1,
                               R0range = seq("from" = 0, "to" = 4, by = R0by),
                               change_param = "Beta",
                               base_params = make_params(delta = 0, nu = 0, mu = 0),
                               useExpandedModel = FALSE,
                               use_equal_removal = FALSE,
                               use_unequal_transmission = FALSE,
                               fixed_params = setdiff(names(base_params), change_param),
                               inits = make_inits(useExpandedModel = useExpandedModel),
                               use_simple_SIR = FALSE){
  finalSizes <- Z1s <- Z2s <- c()
  if (use_unequal_transmission){
    base_params["gamma2"] <- base_params["gamma1"]
  }
  else{
    
  }
  for (R0 in R0range){
    ##Update the value of base_params for the change_param to match the new R0 value.
    ##Currently only Beta is supported.
    if (change_param == "Beta"){
      base_params[[change_param]] <- compute_beta(R0, base_params, use_unequal_transmission = use_unequal_transmission)
    }
    else{
      return("Could not compute")
    }
    if (useExpandedModel){
      ##First solve the expanded model.
      finalSize <- simulate_solve(params = base_params,
                                  use_simple_SIR = use_simple_SIR,
                                  init = inits, 
                                  use_unequal_transmission = use_unequal_transmission,
                                  useExpandedModel = useExpandedModel,
                                  length = 500,
                                  compute_final_size = TRUE)
      
      ##Store Z1 and Z2 separately
      Z1 <- finalSize["Z1"]
      Z2 <- finalSize["Z2"]
      Z1s <- c(Z1s, Z1)
      Z2s <- c(Z2s, Z2)
    }
    else{
      ##Not using the expanded model.
      finalSize <- simulate_solve(params = base_params,
                                  init = inits, 
                                  use_simple_SIR = use_simple_SIR,
                                  use_unequal_transmission = use_unequal_transmission,
                                  useExpandedModel = useExpandedModel,
                                  length = 500,
                                  compute_final_size = TRUE)
      ##Just return the final size from the original model.
      finalSizes <- c(finalSizes, finalSize)
    }
  }
  if (useExpandedModel){
    finalSizedf <- data.frame("Z1" = Z1s,
                              "Z2" = Z2s,
                              "R0" = R0range)
  }
  else{
    ##Either use the expanded model with a single Z based on the formula for Z1 and Z2 or use the basic model with only one Z anyway.
    finalSizedf <- data.frame("Z" = finalSizes,
                              "R0" = R0range)
  }
  
  class(finalSizedf) <- c("finalSize", "data.frame")
  return(finalSizedf) 
}

##Plot Z vs R0.
plot_Z_R0 <- function(finalSizedf){
  library("ggpubr")
  p <- ggpubr::ggline(size = 1.5,plot_type = "l",finalSizedf, x = "R0", y = "Z") +
    labs(title = "Final size vs. R0", y = "Final size of epidemic (proportion)")
  return(p)
}

##Do the same plot as above but solve for Ro(rho)
##We want to take those values and from them solve for rho if we fix the gammas.
finalsize_rho <- function(finalsizedf, gamma1 = 0.08, gamma2 = 0.02, mu = 0, Beta = 0.1){
  finalsizedf$rho <- compute_rho(finalsizedf$R0, gamma1, gamma2, mu = mu, Beta = Beta)
  ##Only consider viable values of rho.
  finalsizedf <- finalsizedf[finalsizedf$rho <= 1,]
  library("ggpubr")
  ggpubr::ggline(size = 1.5,plot_type = "l",finalsizedf, x = "rho", y = "Z") +
    labs(title = "Total proportion of infections", x = "Probability of mild infection $(\\mildprop)$", y = "Final size of epidemic (proportion)", y = "Probability of mild infection")
  
}

##Use the Lambert W function based on the 2006 paper by of Dr. Earn and Dr. Junling Ma.
##Use the standard formula to compute the final size of the epidemic based on an input of R0s
standard_sir_final_size <- function(R0vec, inits = c("S" = 0.99, "I" = 0.01, "R" = 0)){
  library(emdbook)
  finalSizevec <- 1 + (1/R0vec)*emdbook::lambertW(-R0vec*exp(-R0vec))
  return(finalSizevec)
}

##Plot the SIR standard final size plot on top of the final size plot from the original model.
##SIR final size available via Lambert W function or numerical integration of the differential equations.
final_size_plot_show_SIR <- function(useExpandedModel = FALSE,
                                     use_unequal_transmission = FALSE,
                                     use_equal_removal = FALSE,
                                     finalSize = long_run_simulator(useExpandedModel = useExpandedModel,
                                                                    inits = make_inits(0.999999,I1 = 0.0000008,I2 = 0.0000002),
                                                                    use_unequal_transmission = use_unequal_transmission,
                                                                    use_equal_removal = use_equal_removal),
                                     numerically_compute_SIR = FALSE){
  library("ggpubr")
  if (!useExpandedModel){
    ##Plot for just the basic model.
    ##Use the analytical solution for the simple SIR model.
    if (!numerically_compute_SIR){
      df <- dplyr::bind_rows("R0" = rep(finalSize$R0, 2), "Z" = c(finalSize$Z, standard_sir_final_size(finalSize$R0)), type = c(rep("Model with variolation", nrow(finalSize)), rep("Analytical standard SIR model", nrow(finalSize))))
    }
    else{
      simpleSIR_sim <- long_run_simulator(use_simple_SIR = TRUE)
      df <- dplyr::bind_rows("R0" = c(finalSize$R0, simpleSIR_sim$R0), "Z" = c(finalSize$Z, simpleSIR_sim$Z), type = c(rep("Model with variolation", nrow(finalSize)), rep("Numerical standard SIR model", nrow(simpleSIR_sim))))
    }
  }
  else{
    ##Display both Z1 and Z2 separately on the graph.
    df <- dplyr::bind_rows("R0" = rep(finalSize$R0, 3), "Z" = c(finalSize$Z1, finalSize$Z2, standard_sir_final_size(finalSize$R0)), type = c(rep("model with variolation Z1", nrow(finalSize)), rep("model with variolation Z2", nrow(finalSize)), rep("standard SIR model", nrow(finalSize))))
  }
  p <- ggplot2::ggplot(df, mapping = aes(x = R0, y = Z, color = type)) + geom_line(lwd = 1.5,data = subset(df, type == "Analytical standard SIR model" | type == "Numerical standard SIR model"), lwd = 1.5)
  p <- p + geom_line(lwd = 1.5, data = subset(df, type == "Model with variolation"), size = 3) + labs(title = "Z(R0)")
  ##Apply ggppubr publication style theme
  p <- p + ggpubr::theme_pubr()
  p <- p + directlabels::geom_dl(mapping = aes(label = type), method = "extreme.grid")
  ##Remove the legend
  p <- p + theme(legend.position = "none") +  scale_color_manual(values=c("#000000","#FF0000"))
  p
}

##For peak sizes, rename I1 and I2 to mildly infectious and severely infectious
##df should be the intermediary output from compare_peak_size. 
fix_names <- function(df, vary){
  library("tibble")
  library("tidyr")
  df <- as_tibble(df)
  df <- pivot_wider(df, names_from = "infectious_class", values_from = c(peak))
  ##Rename columns.
  df$"Mildly infectious (variolated)" <- df$I1
  df$"Severely infectious" <- df$I2
  df$"Total" <- df$sum
  ##Drop original columns
  df <- df[,!colnames(df) %in% c("I1", "I2", "sum")]
  if (vary == "rho"){
    df <- pivot_longer(df, cols = !rho, names_to = "infectious_class")
  }
  else if (vary == "gammaratio") {
    df <- pivot_longer(df, cols = !gamma_ratio, names_to = "infectious_class")
  }
  if (vary == "bm"){
    df <- pivot_longer(df, cols = !bm, names_to = "infectious_class")
  }
  df$peak <- df$value
  df <- df[,colnames(df) != "value"]
  return(df)
}
##Vary can be "gammaratio", "rho" or "bm", to look at unequal transmission rates.
##Vary gamma 2 if we're varying the gammaratio.
##To make the plot a heatmap, set plot = TRUE, do_heatmap = TRUE, and select the heatmap variable from R0 (heatmap_var = "R0"), or Z.
##Set fixed_R0 to any (not NULL) to fix R0 at that value and see whether the changes in peak height are driven by changing R0 or not.
compare_peak_size <- function(vary = "rho",
                              ##For varying rho and the gamma ratio with unequal transmission rates.
                              use_unequal_transmission = FALSE,
                              params = make_params(),
                              gamma_ratio_from = 0.1,
                              gamma_ratio_to = 0.99,
                              gamma_ratio_by = 0.01,
                              gammaseq = seq(from = gamma_ratio_from,
                                             to = gamma_ratio_to,
                                             by = gamma_ratio_by),
                              rho_to = 1,
                              rho_by = 0.01,
                              rhoseq = seq(from = 0.01,
                                           to = rho_to,
                                           by = rho_by),
                              plot = TRUE,
                              fixed_R0 = NULL,
                              ##This makes sense.
                              bmseq = rhoseq
){
  library(directlabels)
  R0vec <- Zvec <- c()
  if (vary == "gammaratio"){
    df <- tibble::tibble("peak" = c(),
                         "infectious_class" = c(),
                         "gamma_ratio" = c())
    for (gammaratio in gammaseq){
      if (!is.null(fixed_R0)){
        ##If we want to hold R0 constant, find a beta value that cancels out the change in I1 or I2 so that R0 is held constant.
        ##First, fix gamma1 and set gamma2 based on the ratio.
        correction_params <- make_params(gamma1 = params[["gamma1"]],
                                         gamma2 = params[["gamma1"]]*gammaratio,
                                         bm = as.numeric(params["bm"]))
        ##Then calculate beta based on the I1 and I2 changed params and the R0 target value
        correction_beta <- compute_beta(R0 = fixed_R0, params = correction_params, use_unequal_transmission = use_unequal_transmission)
        ##Then put everything into the parameter set
        varied_gamma_params <- make_params(gamma1 = params[["gamma1"]],
                                           gamma2 = params[["gamma1"]]*gammaratio,
                                           Beta = correction_beta,
                                           bm = as.numeric(params["bm"]))
      }
      else{
        ##Don't fix R0 with a correction Beta.
        varied_gamma_params <-  make_params(gamma1 = params[["gamma1"]],
                                            gamma2 = params[["gamma1"]]*gammaratio,
                                            bm = as.numeric(params["bm"]))
      }
      theSim <- simulate_solve(params = varied_gamma_params, use_unequal_transmission = use_unequal_transmission)
      ##Find the peak I1 and I2 height based on the gamma ratio, and also the peak of the sum.
      varied_gamma_peak_I1 <- max(theSim$I1)
      varied_gamma2_peak_I2 <- max(theSim$I2)
      sum_peak <- max(theSim$I1 + theSim$I2)
      ##Compute and store R0 and Z for this combination of parameters.
      R0 <- compute_R0(params = varied_gamma_params)
      Z <- simulate_solve(params = varied_gamma_params, compute_final_size = TRUE, use_unequal_transmission = use_unequal_transmission)
      R0vec <- c(R0vec, R0)
      Zvec <- c(Zvec, Z)
      ##Put these results into a mini tibble.
      holdertibble <- tibble::tibble("peak" = c(varied_gamma_peak_I1,
                                                varied_gamma2_peak_I2,
                                                sum_peak),
                                     "infectious_class" = c("I1",
                                                            "I2",
                                                            "sum"),
                                     "gamma_ratio" = rep(gammaratio, 3))
      df <- dplyr::bind_rows(df, holdertibble)
    }
    df <- fix_names(df, vary = vary)
    if (plot){
      library("ggpubr")
      p <- ggpubr::ggline(size = 1.5, plot_type = "l",df, x = "gamma_ratio", y = "peak", color = "infectious_class", palette = three_curve_sim_palette) +
        labs(title = paste0("Peak Height"), x = "Relative length of mild infections", y = "Peak height (proportion of population)")
      p <- direct.label(p, method = list(cex = 1 ,"smart.grid"))
      p <- p + ggplot2::scale_x_discrete(expand = c(0,0), breaks = seq(0, gamma_ratio_to, 0.2)) 
      p
    }
  }
  else if (vary == "rho"){
    ##The maximum investigated value of rho should be one.
    to <- 1
    ##Store for later.
    df <- tibble::tibble("peak" = c(),
                         "infectious_class" = c(),
                         "rho" = c())
    for (therho in rhoseq){
      if (!is.null(fixed_R0)){
        ##If we want to hold R0 constant, find a beta value that cancels out the change in rho so that R0 is held constant.
        ##First create a set of params with changed rho.
        correction_params <- make_params(rho = therho, bm = as.numeric(params["bm"]))
        ##Then calculate beta based on the rho changed params and the R0 target value
        correction_beta <- compute_beta(R0 = fixed_R0, params = correction_params, use_unequal_transmission = use_unequal_transmission)
        ##Then put everything into the parameter set
        params <- make_params(rho = therho, Beta = correction_beta, bm = as.numeric(params["bm"]))
      }
      else{
        ##Don't fix R0 with a correction Beta.
        params <- make_params(rho = therho, bm = as.numeric(params["bm"]))
      }
      ##Find the peak I1 and I2 height based on the gamma ratio.
      theSim <- simulate_solve(params = params, use_unequal_transmission = use_unequal_transmission)
      ##Find the peak I1 and I2 height based on the gamma ratio, and also the peak of the sum.
      peak_I1 <- max(theSim$I1)
      peak_I2 <- max(theSim$I2)
      sum_peak <- max(theSim$I1 + theSim$I2)
      ##Compute and store R0 and Z for this combination of parameters.
      R0 <- compute_R0(params = params, use_unequal_transmission = use_unequal_transmission)
      Z <- simulate_solve(params = params, compute_final_size = TRUE, use_unequal_transmission = use_unequal_transmission)
      R0vec <- c(R0vec, R0)
      Zvec <- c(Zvec, Z)
      ##Put these results into a mini tibble.
      holdertibble <- tibble::tibble("peak" = c(peak_I1,
                                                peak_I2,
                                                sum_peak),
                                     "infectious_class" = c("I1",
                                                            "I2",
                                                            "sum"),
                                     "rho" = rep(therho, 3))
      df <- dplyr::bind_rows(df, holdertibble)
    }
    if (plot){
      df <- fix_names(df, vary = vary)
      library("ggpubr")
      p <- ggpubr::ggline(size = 1.5, plot_type = "l",df, x = "rho", y = "peak", color = "infectious_class", palette = three_curve_sim_palette) +
        labs(title = "Peak height", x = "Probability of mild infection $(\\mildprop)$", y = "Peak height (proportion)")
      p <- direct.label(p)
      p <- p + ggplot2::scale_x_discrete(expand = c(0,0), breaks = seq(0, 1, 0.1)) 
      p
    }
  }
  else if (vary == "bm"){
    ##Assumes the model with unequal transmission rates, as is reasonable.
    ##The maximum investigated value of bm, like rho, should be one.
    ##Store for later.
    df <- tibble::tibble("peak" = c(),
                         "infectious_class" = c(),
                         "bm" = c())
    for (thebm in bmseq){
      if (!is.null(fixed_R0)){
        ##If we want to hold R0 constant, find a beta value that cancels out the change in rho so that R0 is held constant.
        ##First create a set of params with changed rho.
        correction_params <- make_params(bm = thebm)
        ##Then calculate beta based on the rho changed params and the R0 target value
        correction_beta <- compute_beta(use_unequal_transmission = TRUE, R0 = fixed_R0, params = correction_params)
        ##Then put everything into the parameter set
        params <- make_params(bm = thebm, Beta = correction_beta)
      }
      else{
        ##Don't fix R0 with a correction Beta.
        params <- make_params(bm = thebm)
      }
      ##Find the peak I1 and I2 height based on the gamma ratio.
      theSim <- simulate_solve(params = params, use_unequal_transmission = TRUE)
      ##Find the peak I1 and I2 height based on the gamma ratio, and also the peak of the sum.
      peak_I1 <- max(theSim$I1)
      peak_I2 <- max(theSim$I2)
      sum_peak <- max(theSim$I1 + theSim$I2)
      ##Compute and store R0 and Z for this combination of parameters.
      R0 <- compute_R0(params = params, use_unequal_transmission = TRUE)
      Z <- simulate_solve(params = params, compute_final_size = TRUE, use_unequal_transmission = TRUE)
      R0vec <- c(R0vec, R0)
      Zvec <- c(Zvec, Z)
      ##Put these results into a mini tibble.
      holdertibble <- tibble::tibble("peak" = c(peak_I1,
                                                peak_I2,
                                                sum_peak),
                                     "infectious_class" = c("I1",
                                                            "I2",
                                                            "sum"),
                                     "bm" = rep(thebm, 3))
      df <- dplyr::bind_rows(df, holdertibble)
    }
    if (plot){
      df <- fix_names(df, vary = vary)
      library("ggpubr")
      p <- ggpubr::ggline(size = 1.5, plot_type = "l",df, x = "bm", y = "peak", color = "infectious_class", palette = three_curve_sim_palette) +
        labs(title = "Peak height", x = "Relative infectiousness of mild cases $(b_ \\mathrm{m})$", y = "Peak height (proportion)")
      p <- direct.label(p, method = list(cex = 1 ,"smart.grid"))
      p <- p + ggplot2::scale_x_discrete(expand = c(0,0), breaks = seq(0, 1, 0.1)) 
      p
    }
  }
  else{
    
  }
}

##Rho and gamma ratio and heat map showing the final size.
##Vary gamma 2 to obtain the ratio of the gammas.
##In other words, gammaratio = gamma2/gamma1
##Return the data.
##Set colour_var to chose whether the heatmap shows Z or R0.
##Set fixed_R0 to any (not NULL) to fix R0 at that value and see whether the changes in Z are affected.
do_heatmap <- function(params = make_params(),
                       from = 0,
                       to = 0.5,
                       by = 0.05,
                       theSeq = seq(from = from,
                                    to = to,
                                    by = by),
                       colour_var = "Z",
                       fixed_R0 = NULL,
                       useExpandedModel = FALSE){
  tempdf <- sapply(theSeq, function(gammaratio){
    sapply(theSeq, function(theRho){
      adjust_params <- make_params(gamma1 = params[["gamma1"]],
                                   gamma2 = params[["gamma1"]]*gammaratio,
                                   ##Also set rho accordingly.
                                   rho = theRho)
      ##Compute Z for this combination of parameters.
      if (colour_var == "Z"){
        ##If we're fixing R0.
        if (!is.null(fixed_R0)){
          ##If we want to hold R0 constant, find a beta value that cancels out the change in I1 or I2 so that R0 is held constant.
          ##Then calculate beta based on the I1 and I2 changed params and the R0 target value
          correction_beta <- compute_beta(R0 = fixed_R0, params = adjust_params)
          ##Then put everything into the parameter set
          varied_gamma_params <- make_params(gamma1 = params[["gamma1"]],
                                             gamma2 = params[["gamma1"]]*gammaratio,
                                             rho = theRho,
                                             Beta = correction_beta)
          if (useExpandedModel){
            Z <- simulate_solve(params = varied_gamma_params, useExpandedModel = useExpandedModel, compute_final_size = TRUE)
          }
          else{
            Z <- simulate_solve(params = varied_gamma_params, useExpandedModel = useExpandedModel, compute_final_size = TRUE)
          }
        }
        else{
          ##We're not fixing R0
          if (useExpandedModel){
            Z <- simulate_solve(params = adjust_params, useExpandedModel = useExpandedModel, compute_final_size = TRUE)
          }
          else{
            Z <- simulate_solve(params = adjust_params, useExpandedModel = useExpandedModel, compute_final_size = TRUE)
          }
        }
      }
      else{
        ##We just want R0 for this combination of parameters.
        Z <- compute_R0(params = adjust_params)
      }
      ##If useExpandedModel = TRUE, then Z will have values for both Z1 and Z2 in a named vector. Otherwise, Z will just be a number.
      Z
    })
  })
  
  if (useExpandedModel){
    ##tempdf is double the normal length, every odd entry is Z1, and every even entry is Z2.
    ##If we do this after, the rows won't necessarily be Z2 even, Z1 odd odd.
    z1tempdf <- tempdf[c(TRUE,FALSE),]
    z2tempdf <- tempdf[c(FALSE,TRUE),]
  }
  else{
  }
  if (colour_var == "Z"){
    ##as.vector proceeds down columns of the data frame, so make sure everything else matches.
    if (useExpandedModel){
      z1df <- data.frame("Z" = as.vector(z1tempdf), "rho" = rep(theSeq, length(z1tempdf)), "gammaratio" = rep(theSeq, each = nrow(z1tempdf)))
      z2df <- data.frame("Z" = as.vector(z2tempdf), "rho" = rep(theSeq, length(z2tempdf)), "gammaratio" = rep(theSeq, each = nrow(z2tempdf)))
    }
    else{
      df <- data.frame("Z" = as.vector(tempdf), "rho" = rep(theSeq, length(tempdf)), "gammaratio" = rep(theSeq, each = nrow(tempdf)))
    }
    library(ggplot2)
    ##Use a print statement to show the plot but still return the data.
    if (!is.null(fixed_R0)){
      if (useExpandedModel){
        ##tempdf is double the normal length, every odd entry is Z1, and every even entry is Z2.
        ##Grab those corresponding entries.
        mildPlot <- ggplot(z1df, mapping = aes(x = rho, y = gammaratio, fill = Z)) + geom_tile() + labs(x = "Probability of mild infection $(\\mildprop)$", y = "Relative length of mild infections", title = paste0("Total mild infections")) + scale_fill_distiller(palette = "Reds", direction = 1) + ggpubr::theme_pubr()
        severePlot <- ggplot(z2df, mapping = aes(x = rho, y = gammaratio, fill = Z)) + geom_tile() + labs(x = "Probability of mild infection $(\\mildprop)$", y = "Relative length of mild infections", title = paste0("Total severe infections")) + scale_fill_distiller(palette = "Reds", direction = 1) + ggpubr::theme_pubr()
        ##Display both plot.
        return(ggarrange(mildPlot, severePlot, ncol = 1))
      }
      else{
        return(ggplot(df, mapping = aes(x = rho, y = gammaratio, fill = Z)) + geom_tile() + labs(x = "Probability of mild infection $(\\mildprop)$", y = "Relative length of mild infections", title = paste0("Final Size")) + scale_fill_distiller(palette = "Reds", direction = 1) + ggpubr::theme_pubr())
      }
    }
    else{
      if (useExpandedModel){
        ##Grab those corresponding entries.
        mildPlot <- ggplot(z1df, mapping = aes(x = rho, y = gammaratio, fill = Z)) + geom_tile() + labs(y = "Relative length of mild infections", title = "Total mild infections") + scale_fill_distiller(palette = "Reds", direction = 1) + ggpubr::theme_pubr() + labs(x = "Probability of mild infection", y = "Relative length of mild infections", title = paste0("Total mild infections")) + scale_fill_distiller(palette = "Reds", direction = 1) + ggpubr::theme_pubr()
        severePlot <- ggplot(z2df, mapping = aes(x = rho, y = gammaratio, fill = Z)) + geom_tile() + labs(y = "Relative length of mild infections", title = "Total severe infections") + scale_fill_distiller(palette = "Reds", direction = 1) + ggpubr::theme_pubr()
        ##Display both plot.
        return(ggarrange(mildPlot, severePlot, ncol = 1))
      }
      else{
        return(ggplot(df, mapping = aes(x = rho, y = gammaratio, fill = Z)) + geom_tile() + labs(x = "Probability of mild infection.", y = "Relative length of mild infections", title = "Final Size") + scale_fill_distiller(palette = "Reds", direction = 1) + ggpubr::theme_pubr())
      }
    }
  }
  ##Colour variable is R0.
  else{
    ##as.vector proceeds down columns of the data frame, so make sure everything else matches
    df <- data.frame("R0" = as.vector(tempdf), "rho" = rep(theSeq, length(tempdf)), "gammaratio" = rep(theSeq, each = nrow(tempdf)))
    library(ggplot2)
    ##Use a print statement to show the plot but still return the data.
    return(ggplot(df, mapping = aes(x = rho, y = gammaratio, fill = R0)) + geom_tile() + labs(y = "Relative Length of Mild Infections", title = "R0 as a function of rho and gamma2/gamma1") + scale_fill_distiller(palette = "Reds", direction = 1) + ggpubr::theme_pubr())
  }
  return(df)
}

##Alternative to make_inits.
###Package initial values to pass to the simulation to test for the herd immunity threshold by declaring how much of the population is susceptible when time = 0.
##In other words, S = 1-P.
immunity_inits <- function(P = 0.5,
                           I1 = 0.008,
                           I2 = 0.002){
  return(
    c(
      "S" = 1 - P,
      "I1" = I1,
      "I2" = I2,
      "R" = 0)
  )
}
##Given a set of parameters for a simulation, find the smallest starting immunization value for the simulation so that I1 and I2 don't increase above their initial values.
find_herd_immunity_threshold <- function(pars = make_params(),
                                         length = 360,
                                         from = 0.01,
                                         to = 1,
                                         by = 0.1,
                                         theSeq = seq(from = from,
                                                      to = to,
                                                      by = by),
                                         I1_init = 0.008,
                                         I2_init = 0.002){
  theP <- 1
  for (P in theSeq){
    ##For given levels of population level immunity, simulate a pandemic and find the herd immunity threhold.
    sim <- simulate_solve(params = pars, init = immunity_inits(P = P, I1 = I1_init, I2 = I2_init), length = length)
    ##If I_1 and I_2 don't increase at all (or very little) for the course of the epidemic (they stay at or below their initial values), we are at below the herd immunity threshold.
    if (sum(sim$I1 > I1_init) == 0 & sum(sim$I2 > I2_init) == 0){
      ##We want to be at the herd immunity threshold, not above it, so find the smallest of these values.
      if (P < theP)
        theP <- P
    }
  }
  return(theP)
}

##For a range of initial values for I1 and I2, find the herd immunity threshold (the level of initial population level immunity such that the number of infectives does not grow past initial values (only decreases)) and return that in a dataset.
herd_immunity_relationship_old <- function(from = 0.00001,
                                           ##Only consider a sequence of realistic initial values for I1 and I2.
                                           to = 0.01,
                                           by = 0.001,
                                           theSeq = seq(from = from,
                                                        to = to,
                                                        by = by),
                                           plot = TRUE){
  Ps <- I1s <- I2s <- c()
  for (I1 in theSeq){
    for (I2 in theSeq){
      ##For a range of initial values, find the herd immunity threshold for a simulation with those stating values.
      P <- find_herd_immunity_threshold(pars = make_params(), I1_init = I1, I2_init = I2)
      Ps <- c(Ps, P)
      I1s <- c(I1s, I1)
      I2s <- c(I2s, I2)
    }
  }
  df <- data.frame("P" = Ps,
                   "I1_init" = I1s,
                   "I2_init" = I2s)
  if (plot){
    library(ggplot2)
    p <- ggplot(df, mapping = aes(x = I1_init, y = I2_init, fill = P)) + geom_tile() + labs(x = "initial infections in Im",                                                                                          y = "initial infections in I2",
                                                                                            title = "Initial infections in Im, Is, herd immunity threshold")
    p <- direct.label(p)
    return(p)
  }
  return(df)                             
}

##Plot the herd immunity threshold for a given set of parameters.
##Include a switch to grab compute parameters from R0.
p_func <- function(params, inits, R0 = NULL){
  if (is.null(R0)){
    p <- 1 - (params[["gamma1"]]*(inits[["I1"]]/inits[["I2"]]) + params[["gamma2"]])/(params["Beta"]*((inits[["I1"]]/inits[["I2"]]) + 1))
    ##Remove names
    p <- as.numeric(p)
    return(p)
  }
  else{
    correction_beta <- compute_beta(R0 = R0, params = params)
    ##Ger rid of the original beta and back compute it to hit the right R0 value.
    p <- 1 - (params[["gamma1"]]*(inits[["I1"]]/inits[["I2"]]) + params[["gamma2"]])/(correction_beta*((inits[["I1"]]/inits[["I2"]]) + 1))
  }
}

herd_immunity_relationship <- function(from = 0.1,
                                       to = 1,
                                       by = 0.05,
                                       theSeq = seq(from = from,
                                                    to = to,
                                                    by = by),
                                       params = make_params(),
                                       ##To grab the starting values from.
                                       base_inits = make_inits(),
                                       R0_from = 1,
                                       R0_to = 3,
                                       R0_by = 0.5,
                                       R0seq = seq(from = R0_from,
                                                   to = R0_to,
                                                   by = R0_by)){
  library(dplyr)
  pData <- data.frame("p" = c(),"I1_to_I2" = c())
  ##Repeat for a range of values of R0.
  theData <- lapply(R0seq, function(R0){
    pvec <- sapply(theSeq, function(ratio){
      ##Fix I2 and vary I1 to vary the ratio.
      I_2 <-  base_inits[["I1"]]
      I_1 <-  I_2 * ratio
      ##Balance out the initial values with S0 so that we don't start above one.
      inits <- make_inits(S = 1- (I_1 + I_2),
                          I1 = I_1,
                          I2 = I_2)
      p <- p_func(params, inits = inits, R0 = R0)
      return(p)
    })
    pData <- dplyr::bind_rows(pData, data.frame("p" = pvec ,"I1_to_I2" = theSeq, "R0" = R0))
    pData
  })   
  ##Combine it all into one data frame.
  theData <- bind_rows(theData)
  library(ggpubr)
  p <- ggline(size = 1.5,plot_type = "l",theData, x = "I1_to_I2", y = "p", color = "R0") + labs(x = "Initial Im/Is", y = "Herd Immunity Threshold", title = "Herd Immunity Threshold and initial Im/Is")
  p <- directlabels::direct.label(p)
  p
}

final_size_line <- function(
  gammaratio = NULL,
  ##Can set rho in here.
  params = make_params(),
  from = 0,
  to = 0.5,
  by = 0.05,
  ##Allow to fix value of bm instead of gammaratio, which only makes sense when varying rho.
  bm = NULL,
  theSeq = seq(from = from,
               to = to,
               by = by),
  fixed_R0 = NULL,
  ##Only varying rho, bm, and the gammaratio are supported right now.
  vary = "rho",
  ##Just relevant for varying rho, changing varying bm assumes the model with differential transmission between severe and variolated cases
  use_unequal_transmission = FALSE,
  ##Can also vary bm
  useExpandedModel = TRUE){
  library(directlabels)
  if (vary == "rho"){
    tempdf <- sapply(theSeq, function(theRho){
      if (!is.null(gammaratio)){
        ##Fix the gammaratio
        adjust_params <- make_params(gamma1 = params[["gamma1"]],
                                     gamma2 = params[["gamma1"]]*gammaratio,
                                     ##Also set rho accordingly.
                                     rho = theRho)
      }
      else{
        ##Fix bm
        adjust_params <- make_params(bm = bm,
                                     rho = theRho)
      }
      
      ##If we're fixing R0.
      if (!is.null(fixed_R0)){
        ##If we want to hold R0 constant, find a beta value that cancels out the change in I1 or I2 so that R0 is held constant.
        ##Then calculate beta based on the I1 and I2 changed params and the R0 target value
        correction_beta <- compute_beta(R0 = fixed_R0, params = adjust_params, use_unequal_transmission = use_unequal_transmission)
        if (!is.null(gammaratio)){
          ##Then put everything into the parameter set
          varied_gamma_params <- make_params(gamma1 = params[["gamma1"]],
                                             gamma2 = params[["gamma1"]]*gammaratio,
                                             rho = theRho,
                                             Beta = correction_beta)
        }
        else{
          # we're fixing bm'
          varied_gamma_params <- make_params(bm = bm,
                                             rho = theRho,
                                             Beta = correction_beta)
        }
        
        if (useExpandedModel){
          Z <- simulate_solve(params = varied_gamma_params, useExpandedModel = useExpandedModel, compute_final_size = TRUE, use_unequal_transmission =  use_unequal_transmission)
        }
        ##We're not fixing R0
        else{
          Z <- simulate_solve(params = varied_gamma_params, useExpandedModel = useExpandedModel, compute_final_size = TRUE, use_unequal_transmission =  use_unequal_transmission)
        }
      }
      else{
        Z <- simulate_solve(params = adjust_params, useExpandedModel = useExpandedModel, compute_final_size = TRUE, use_unequal_transmission =  use_unequal_transmission)
      }
      ##If useExpandedModel = TRUE, then Z will have values for both Z1 and Z2 in a named vector. Otherwise, Z will just be a number.
      Z
    })
    
    if (useExpandedModel){
      ##tempdf is double the normal length, every odd entry is Z1, and every even entry is Z2.
      ##If we do this after, the rows won't necessarily be Z2 even, Z1 odd.
      z1tempdf <- tempdf[c(TRUE,FALSE),]
      z2tempdf <- tempdf[c(FALSE,TRUE),]
    }
    else{
    }
    if (useExpandedModel){
      z1df <- data.frame("Z" = as.vector(z1tempdf), "rho" = rep(theSeq, length(z1tempdf)), type = "Mildly infectious (variolated)")
      z2df <- data.frame("Z" = as.vector(z2tempdf), "rho" = rep(theSeq, length(z2tempdf)),  type = "Severely infectious")
      dftotal <- data.frame("Z" = z1df$Z + z2df$Z, "rho" = z1df$rho,  type = "Total")
      df <- dplyr::bind_rows(z1df, z2df, dftotal)
    }
    else{
      df <- data.frame("Z" = as.vector(tempdf), "rho" = rep(theSeq, length(tempdf)))
    }
    library(ggpubr)
    p <- ggline(plot_type = "l",data = df, x = "rho", y  = "Z", color = "type", size  = 3, palette = three_curve_sim_palette)
    p <- direct.label(p, method = list(cex = 1 ,"smart.grid"))
    p <- p + labs(x = "Probability of mild infection $(\\mildprop)$", y = "Total proportion infected")
  }
  
  else if (vary == "bm"){
    ##Vary the reduction in transmission for variolated cases.
    tempdf <- sapply(theSeq, function(thebm){
      adjust_params <- make_params(gamma1 = params[["gamma1"]],
                                   gamma2 = params[["gamma1"]]*gammaratio,
                                   ##Also set bm accordingly.
                                   bm = thebm)
      ##If we're fixing R0.
      if (!is.null(fixed_R0)){
        ##If we want to hold R0 constant, find a beta value that cancels out the change in I1 or I2 so that R0 is held constant.
        ##Then calculate beta based on the I1 and I2 changed params and the R0 target value
        correction_beta <- compute_beta(R0 = fixed_R0,  use_unequal_transmission = TRUE, params = adjust_params)
        ##Then put everything into the parameter set
        varied_gamma_params <- make_params(gamma1 = params[["gamma1"]],
                                           gamma2 = params[["gamma1"]]*gammaratio,
                                           bm = thebm,
                                           Beta = correction_beta)
        if (useExpandedModel){
          Z <- simulate_solve(use_unequal_transmission = TRUE, params = varied_gamma_params, useExpandedModel = useExpandedModel, compute_final_size = TRUE)
        }
        ##We're not fixing R0
        else{
          Z <- simulate_solve(use_unequal_transmission = TRUE, params = varied_gamma_params, useExpandedModel = useExpandedModel, compute_final_size = TRUE)
        }
      }
      else{
        Z <- simulate_solve(use_unequal_transmission = TRUE, params = adjust_params, useExpandedModel = useExpandedModel, compute_final_size = TRUE)
      }
      ##If useExpandedModel = TRUE, then Z will have values for both Z1 and Z2 in a named vector. Otherwise, Z will just be a number.
      Z
    })
    
    if (useExpandedModel){
      ##tempdf is double the normal length, every odd entry is Z1, and every even entry is Z2.
      ##If we do this after, the rows won't necessarily be Z2 even, Z1 odd.
      z1tempdf <- tempdf[c(TRUE,FALSE),]
      z2tempdf <- tempdf[c(FALSE,TRUE),]
    }
    else{
    }
    if (useExpandedModel){
      z1df <- data.frame("Z" = as.vector(z1tempdf), "bm" = rep(theSeq, length(z1tempdf)), type = "Mildly infectious (variolated)")
      z2df <- data.frame("Z" = as.vector(z2tempdf), "bm" = rep(theSeq, length(z2tempdf)),  type = "Severely infectious")
      dftotal <- data.frame("Z" = z1df$Z + z2df$Z, "bm" = z1df$bm,  type = "Total")
      df <- dplyr::bind_rows(z1df, z2df, dftotal)
    }
    else{
      df <- data.frame("Z" = as.vector(tempdf), "bm" = rep(theSeq, length(tempdf)))
    }
    library(ggpubr)
    p <- ggline(plot_type = "l", data = df, x = "bm", y  = "Z", color = "type", size  = 3, palette = three_curve_sim_palette)
    p <- direct.label(p, method = list(cex = 1 ,"smart.grid"))
    p <- p + labs(x = "Relative infectivity of variolated cases $(b_ \\mathrm{m})$", y = "Total proportion infected")
  }
  else if (vary == "gammaratio"){
    ##Vary the reduction in transmission for variolated cases.
    tempdf <- sapply(theSeq, function(thegammaratio){
      adjust_params <- make_params(gamma1 = params[["gamma1"]],
                                   gamma2 = params[["gamma1"]]*thegammaratio,
                                   ##Also set bm accordingly.
                                   bm = bm,
                                   ##Let us set rho too.
                                   rho = as.numeric(params["rho"]))
      ##If we're fixing R0.
      if (!is.null(fixed_R0)){
        ##If we want to hold R0 constant, find a beta value that cancels out the change in I1 or I2 so that R0 is held constant.
        ##Then calculate beta based on the I1 and I2 changed params and the R0 target value
        correction_beta <- compute_beta(R0 = fixed_R0,  use_unequal_transmission = use_unequal_transmission, params = adjust_params)
        ##Then put everything into the parameter set
        varied_gamma_params <- make_params(gamma1 = params[["gamma1"]],
                                           gamma2 = params[["gamma1"]]*thegammaratio,
                                           bm = bm,
                                           rho = as.numeric(params["rho"]),
                                           Beta = correction_beta)
        if (useExpandedModel){
          Z <- simulate_solve(use_unequal_transmission = use_unequal_transmission, params = varied_gamma_params, useExpandedModel = useExpandedModel, compute_final_size = TRUE)
        }
        ##We're not fixing R0
        else{
          Z <- simulate_solve(use_unequal_transmission = use_unequal_transmission, params = varied_gamma_params, useExpandedModel = useExpandedModel, compute_final_size = TRUE)
        }
      }
      else{
        Z <- simulate_solve(use_unequal_transmission = use_unequal_transmission, params = adjust_params, useExpandedModel = useExpandedModel, compute_final_size = TRUE)
      }
      ##If useExpandedModel = TRUE, then Z will have values for both Z1 and Z2 in a named vector. Otherwise, Z will just be a number.
      Z
    })
    
    if (useExpandedModel){
      ##tempdf is double the normal length, every odd entry is Z1, and every even entry is Z2.
      ##If we do this after, the rows won't necessarily be Z2 even, Z1 odd.
      z1tempdf <- tempdf[c(TRUE,FALSE),]
      z2tempdf <- tempdf[c(FALSE,TRUE),]
    }
    else{
    }
    if (useExpandedModel){
      z1df <- data.frame("Z" = as.vector(z1tempdf), "gammaratio" = rep(theSeq, length(z1tempdf)), type = "Mildly infectious (variolated)")
      z2df <- data.frame("Z" = as.vector(z2tempdf), "gammaratio" = rep(theSeq, length(z2tempdf)),  type = "Severely infectious")
      dftotal <- data.frame("Z" = z1df$Z + z2df$Z, "gammaratio" = z1df$gammaratio,  type = "Total")
      df <- dplyr::bind_rows(z1df, z2df, dftotal)
    }
    else{
      df <- data.frame("Z" = as.vector(tempdf), "gammaratio" = rep(theSeq, length(tempdf)))
    }
    library(ggpubr)
    p <- ggline(plot_type = "l",data = df, x = "gammaratio", y  = "Z", color = "type", size  = 3, palette = three_curve_sim_palette)
    p <- direct.label(p, method = list(cex = 1 ,"smart.grid"))
    p <- p + labs(x = "Relative length of mild infections $\\frac{\\gammasevere}{\\gammamild}$", y = "Total proportion infected")
  }
  else{
  }
  ##Add a title with a properly typset R0 to the plot  if it's fixed
  if (!is.null(fixed_R0)){
    #  p <- p + labs(title = paste0("R0 = " , fixed_R0))
  }
  else{
  }
  return(p)
}

##Used for EE plot
internal_renamer <- function(tempdf,
                             theSeq,
                             theNames,
                             vary){
  library("tibble")
  library("tidyr")
  tempdf <- as_tibble(tempdf)
  ##We have columns for different values of the parameter we are varying 
  colnames(tempdf) <- theSeq
  ##Properly name the states
  tempdf$Severity <- theNames
  ##Rename columns and reformat into long form.
  tempdf <- pivot_longer(
    tempdf, 
    cols = !Severity,
    names_to = vary,
    values_to = "value")
  return(tempdf)
}

##Pretty axis labels
pretty_axis_labels <- function(p, to = 0.0004, at = seq(from = 0, to = to, by = 0.0001)){
  library("scales")
  return(p +
           scale_y_continuous(breaks = at,
                              labels=latex2exp::TeX(sprintf(at, fmt = "$%g$"))))
}

##vary can be "rho", "bm", "gammaratio", "R0", or "p"
EE_plot <- function(##Parameter to vary on the x axis
  vary = "rho",
  use_unequal_transmission = FALSE,
  from = 0,
  to = 1,
  by = 0.05,
  ##Allow to fix anything not being varied.
  ##Supports all params including gammaratio (which uses the base gamma1 value)
  fixed_params = c(NULL),
  ##Fix R0 if you want.
  fixed_R0 = NULL,
  ##Compute and show total prevalence too.
  compute_show_total = FALSE,
  withP = FALSE,
  ##Show multiple R0s on the plot for p
  p_multiple_R0 = FALSE,
  p_R0_seq = c(1.5, 2, 4, 8),
  ##Sequence of values to interate over for whatever parameter we are changing on the x axis
  theSeq = seq(from = from,
               to = to,
               by = by),
  ##Equilibrium values to hide from the final plot
  drop = c("S", "R"),
  ##Only used if vary = "R0"
  ##Main values
  R0_from = 1,
  R0_to = 6,
  R0_by = 0.5,
  R0Seq =  seq(from = R0_from,
               to = R0_to,
               by = R0_by)){
  ##Avoid errors at startup
  library("directlabels")
  if (vary %in% names(fixed_params)){
    print("can't both fix and vary the same parameter")
    return(NA)
  }
  else {
  }
  ##If we're fixing the gammaratio, set gamma2 in the fixed params beind the scenes and use that from now on
  if ("gammaratio" %in% names(fixed_params)){
    fixed_params <- c(fixed_params, "gamma2" = as.numeric(make_params()["gamma1"]*fixed_params["gammaratio"]))
    ##Get rid of the gammaratio from the fixed params
    fixed_params <- fixed_params[names(fixed_params) != "gammaratio"]
  }
  else{
  }
  if (vary == "R0"){
    ##Switch the sequence of values to iterate over to be ones that make sense for R0
    theSeq <- R0Seq
  }
  else if (vary == "p"){
    if (is.null(fixed_R0)){
      thePars <- make_params()
      thePars[names(fixed_params)] <- c(fixed_params)
      p_crit <- 1-1/compute_R0(params = thePars, use_unequal_transmission = use_unequal_transmission)
    }
    else{
      ##We can do this here because changing p has no effect on R0. So we can either compute it from the values or fix it.
      p_crit <- 1 - 1/fixed_R0
    }
    ##If we're varying the proportion vaccinated, only iterate up to the critical vaccination proportion (pcrit)
    theSeq <- seq(from = 0,
                  to = p_crit,
                  by = by)
  }
  else{
  }
  if (!p_multiple_R0){
    tempdf <- sapply(theSeq, function(paramval){
      ##Start with a fresh set of params and add the adjust parameters and the varying parameters to it.
      adjust_params <- make_params()
      if (vary == "gammaratio"){
        if (!is.null(fixed_params)){
          adjust_params[names(fixed_params)] <- fixed_params
        }
        else{
        }
        adjust_params["gamma2"] <- as.numeric(make_params()["gamma1"]*paramval)
      }
      else{
        ## vary is rho, bm, R0, or p
        if (!is.null(fixed_params)){
          adjust_params[names(fixed_params)] <- c(fixed_params)
        }
        else{
        }
        if (vary != "R0"){
          ##This doesn't make sense to do if we're varying R0
          adjust_params[vary] <- paramval
        }
        
        else{
          ##Vary = R0
          ##Vary R0, named as R0, by varying Beta
          correction_beta <- compute_beta(R0 = paramval,  use_unequal_transmission = use_unequal_transmission, params = adjust_params)
          adjust_params <- make_params(base_params = adjust_params,
                                       change_pars = c("Beta"),
                                       Beta = correction_beta)
        }
      }
      if (vary != "R0"){
        ##Again, this doesn't make sense to do if we're varying R0
        ##Let us fix R0
        if (!is.null(fixed_R0) ){
          correction_beta <- compute_beta(R0 = fixed_R0,  use_unequal_transmission = use_unequal_transmission, params = adjust_params)
          adjust_params <- make_params(base_params = adjust_params,
                                       change_pars = c("Beta"),
                                       Beta = correction_beta)
        }
        else{
          #We're not fixing R0
        }
      }
      else{
        ##Do nothing.
      }
      EE <- compute_EE(params = adjust_params, use_unequal_transmission = use_unequal_transmission, withP = withP)
      if (compute_show_total){
        EE["total"] <- EE["Im"] + EE["Is"]
      }
      else{
      }
      return(EE)
    })
  }
  else{
    ##Vary p for multiple values of R0
    ##Sort the R0 values in increasing order.
    p_R0_seq <- p_R0_seq[order(p_R0_seq)]
    full_p_seq <-  seq(from = from,
                       to = 1-1/p_R0_seq[length(p_R0_seq)],
                       by = by)
    tempdf <- lapply(p_R0_seq, function(R0val){
      theSeq <- seq(from = from,
                    to = 1-1/R0val,
                    by = by)
      sapply(theSeq, function(paramval){
        ##Start with a fresh set of params and add the adjust parameters and the varying parameters to it.
        adjust_params <- make_params()
        if (!is.null(fixed_params)){
          adjust_params[names(fixed_params)] <- c(fixed_params)
        }
        else{
        }
        adjust_params["p"] <- paramval
        ##Again, this doesn't make sense to do if we're varying R0
        ##Let us fix R0
        correction_beta <- compute_beta(R0 = R0val,  use_unequal_transmission = use_unequal_transmission, params = adjust_params)
        adjust_params <- make_params(base_params = adjust_params,
                                     change_pars = c("Beta"),
                                     Beta = correction_beta)
        EE <- compute_EE(params = adjust_params, use_unequal_transmission = use_unequal_transmission, withP = withP)
        if (compute_show_total){
          EE["total"] <- EE["Im"] + EE["Is"]
        }
        else{
        }
        return(EE)
      })
    })
  }
  ##Grab the names of the states
  theNames <- c("S", "Mildly Infectious", "Severely Infectious", "R")
  if (compute_show_total == TRUE){
    ##Add a column name for total prevalence
    theNames <- c(theNames, "Total")
  }
  else{
  }
  ##Use a tibble so we can keep the spaces in the names
  library("tidyr")
  library("tibble")
  if (!p_multiple_R0){
    tempdf <- internal_renamer(tempdf, theSeq = theSeq, theNames = theNames, vary = vary)
  }
  else{
    ##We have instead columns for different values of R0
    i <- 1
    while (i <= length(p_R0_seq)){
      ##The theSeq here is the longest sequence, corresponding to the largest R0 value
      tempdf[[i]] <-  internal_renamer(tempdf = tempdf[[i]], theSeq = full_p_seq[1:length(seq(from = from, to = 1-1/p_R0_seq[i], by = by))], theNames = theNames, vary = vary)
      tempdf[[i]]$R0 <- paste0(p_R0_seq[i])
      ##Differentiate between both R0 and State on the plot
      tempdf[[i]]$id <- paste0(tempdf[[i]]$Severity, ", R0 = ", tempdf[[i]]$R0)
      i <- i + 1
    }
    ##Smash them all together
    tempdf <- dplyr::bind_rows(tempdf)
    tempdf$R0 <- factor(tempdf$"R0", levels = unique(tempdf$"R0"[order(tempdf$R0)]))
  }
  ##Drop the curves we don't want
  ##Drop the entries we don't want
  ##We have to do it this way because the data is in long form
  tempdf <- tempdf[!tempdf$Severity %in% drop,]
  ##Proper colours.
  pal <- three_curve_sim_palette
  ##Plotting time.
  library("ggpubr")
  if (vary == "rho"){
    ##Properly colour the I2 curve if that's the only curve we're plotting
    if ("S" %in% drop & "Mildly Infectious" %in% drop & "R" %in% drop){
      pal <- pal[2]
    }
    else{
    }
    p <- ggline(size = 2, plot_type = "l",data = tempdf, x= "rho", y = "value", color = "Severity", palette = pal) +
      labs(title = "Prevalence of Cases at Endemic Equilibrium as a function of $\\mildprop$",
           x = "Probability of mild infection $(\\mildprop)$",
           y = "Equilibrium Prevalence")
    
  }
  else if (vary == "bm"){
    p <- ggline(size = 2, plot_type = "l",data = tempdf, x= "bm", y = "value", color = "Severity", palette = pal) +
      labs(title = "Prevalence of Cases at Endemic Equilibrium",
           x = "Relative transmission of mild cases ($\\frac{\\betamild}{\\betasevere}$)",
           y = "Equilibrium Prevalence")
  }
  else if (vary == "gammaratio"){
    p <- ggline(size = 2, plot_type = "l",data = tempdf, x= "gammaratio", y = "value", color = "Severity", palette = pal) +
      labs(title = "Prevalence of Cases at Endemic Equilibrium",
           x = "Relative length of mild cases ($\\frac{\\gammasevere}{\\gammamild}$)",
           y = "Equilibrium Prevalence")
  }
  else if (vary == "p"){
    if (!p_multiple_R0){
      p <- ggline(size = 2 ,plot_type = "l",data = tempdf, x= "p", y = "value", color = "Severity", palette = pal) +
        labs(title = "Prevalence of Cases at Endemic Equilibrium",
             x = "Proportion vaccinated",
             y = "Equilibrium Prevalence")
    }
    else{
      library("ggplot2")
      ##group = id draws a line for each value of R0 and severity
      p <- ggplot(data = tempdf, mapping = aes(x = p, y = value, linetype = Severity, colour = R0, group = id)) + geom_line(lwd = 1.5, size = 1.2) +
        labs(title = "Prevalence of Cases at Endemic Equilibrium",
             x = "Proportion vaccinated",
             y = "Equilibrium Prevalence") +
        ##Because we used a ggplot call rather than a ggpubr call. 
        theme_pubr()
      p <- pretty_axis_labels(p)
      return(p)
    }
  }
  else{
    ##vary == R0
    if (!is.null(fixed_params["rho"])){
      if (!is.na(fixed_params["rho"])){
        theRho <- fixed_params["rho"]
        ##Put rho in the title if it's been set
        theTitle <- as.expression(bquote(rho~"="~.(theRho)))
      }
    }
    else{
      ##Title as normal
      theTitle <- "Prevalence of Cases at Endemic Equilibrium"
    }
    p <- ggline(plot_type = "l",data = tempdf, x= "R0", y = "value", color = "Severity", palette = pal, size = 2) +
      labs(title = theTitle,
           x = "R0",
           y = "Proportion of population")
  }
  ##Direct label all plots
  p <- direct.label(p, method = "smart.grid")
  ##For the gammaratio, the plot isn't on a log scale, so skip this step.
  if (vary != "gammaratio"){
    ##The + 0.001 lets us get the top axis label
    p <- pretty_axis_labels(p, to = max(tempdf$value) + 0.0001)
  }
  else{
  }
  return(p)
}

##To work in main_results_plot, order of params needs to be R0Seq,  "vary", params, use_unequal_transmission
peak_R0_bm_plot <- function(R0Seq = c(3, 3*1.5, 28.4),
                            vary = "rho",
                            use_directlabels = FALSE,
                            params = make_params(),
                            use_unequal_transmission = TRUE,
                            let_R0_vary = FALSE,
                            bmseq = c(1),
                            ## bmseq = c(0.10, 0.50),
                            ##Set to "severe" to maximize and only plot severe peak prevalence
                            maximize = "severe",
                            rho_to = 1,
                            rho_by = 0.1,
                            ##Use the same sequence for rho and bm
                            rhoseq = seq(from = 0.001,
                                         to = rho_to,
                                         by = rho_by),
                            gamma_to = 2,
                            gamma_by = 0.05,
                            ##Add some new values for the gammasequence
                            gammaseq = c(seq(from = 0.01,
                                             to = 0.2,
                                             by = 0.1),
                                         seq(from = 0.2,
                                             to = gamma_to,
                                             by = gamma_by))){
  if (vary != "gammaratio"){
    fix_R0_at <- params[[vary]]
  }
  else{
    fix_R0_at <- params[["gamma2"]]/params[["gamma1"]]
  }
  if (maximize == "total"){
    df <- tibble::tibble("peak" = c(),
                         "rho" = c(),
                         "$\\R_0$" = c(),
                         "bm" = c())
  }
  else{
    ##The varying bm is too much information for this plot
    df <- tibble::tibble("peak" = c(),
                         "rho" = c(),
                         "$\\R_0$" = c())
  }
  if (vary == "rho"){
    for (therho in rhoseq){
      for (fixed_R0 in R0Seq){
        for (bm in bmseq){
          if (maximize == "severe"){
            ##If we're maximizing severe prevalence, grab bm from the params and stick with it, overwriting the valye grabbed from the bmseq.
            bm <- as.numeric(params["bm"])
          }
          else{
          }
          if (!let_R0_vary){
            correction_params <- make_params(base_params = params, change_pars = c("rho", "bm"), rho = therho, bm = bm)
            correction_beta <- compute_beta(R0 = fixed_R0, params = correction_params, use_unequal_transmission = use_unequal_transmission)
            params <- make_params(base_params = correction_params,
                                  change_pars = c("Beta"),
                                  Beta = correction_beta)
          }
          else{
            ##compute the Beta so that at vary = fix_R0_at, we hit R0 = R0val
            fix_R0_params <- make_params(base_params = params,
                                         change_pars = c("rho"),
                                         rho = fix_R0_at)
            correction_beta <- compute_beta(R0 = fixed_R0,  use_unequal_transmission = use_unequal_transmission, params = fix_R0_params)
            adjust_params <- make_params(base_params = params,
                                         change_pars = c("Beta"),
                                         Beta = correction_beta)
            adjust_params["rho"] <- therho
            params <- adjust_params
          }
          
          theSim <- simulate_solve(params = params, use_unequal_transmission = use_unequal_transmission)
          
          if (maximize == "total"){
            sum_peak <- max(theSim$I1 + theSim$I2)
            holdertibble <- tibble::tibble("peak" = c(sum_peak),
                                           "rho" = c(therho),
                                           "$\\R_0$" = c(toString(fixed_R0)),
                                           "bm" = c(paste0(bm*100, "\\%")))
          }
          else{
            ##Maximize just severe prevalence
            sum_peak <- max(theSim$I2)
            holdertibble <- tibble::tibble("peak" = c(sum_peak),
                                           "rho" = c(therho),
                                           "$\\R_0$" = c(toString(fixed_R0)))
          }
          df <- dplyr::bind_rows(df, holdertibble)
        }
      }
    }
    library("ggplot2")
    library("ggpubr")
    library("latex2exp")
    ##Rename bm column
    if (maximize == "total"){
      df <- tibble::tibble(df)
      df$"Relative Infectivity of Mild Cases" <- df$bm
      df <- df[, colnames(df) != "bm"]
      p <- ggplot(df, mapping = aes(rho, peak, color = factor(`$\\R_0$`, levels = R0Seq), linetype = `Relative Infectivity of Mild Cases`)) +
        geom_line(lwd = 1.5,size = 1.2) +
        theme_bw()
      p <- p +  labs(x = "probability of mild infection ($\\mildprop$)", y = "peak height", title = "peak height vs $\\mildprop$")
    }
    else{
      ##Maximize == severe
      df <- tibble::tibble(df)
      p <- ggplot(df, mapping = aes(rho, peak, color = factor(`$\\R_0$`, levels = R0Seq))) +
        geom_line(lwd = 1.5, size = 1.2)
      if (use_directlabels){
        p <- p +   directlabels::geom_dl(mapping = aes(label = `$\\R_0$`) , method = "last.bumpup")
      }
      else{
        ##no direct labels, but still no legend
      }
      p <- p  + theme_bw() +  theme(legend.position = "none")
      p <- p +  labs(x = "probability of mild infection ($\\mildprop$)", y = "severe peak prevalence", title = "peak prevalence (severe cases)")
    }
  }
  else if (vary == "bm" | vary == "gammaratio") {
    ## we're varying bm or the gammaratio
    if (vary == "gammaratio"){
      holdertibble <- tibble::tibble("peak" = c(),
                                     "$\\frac{\\gammasevere}{\\gammamild}$" = c(),
                                     "$\\R_0$" = c())
    }
    else{
      df <- tibble::tibble("peak" = c(),
                           "$\\reducemild$" = c(),
                           "$\\R_0$" = c())
    }
    ##Use the rho sequence for the other two varying parameters, but if we're varying the gammaratio, consider values above 1
    if (vary == "gammaratio"){
      rhoseq <- gammaseq
    }
    else{
    }
    for (theVal in rhoseq){
      for (fixed_R0 in R0Seq){
        ##If we're maximizing severe prevalence, grab bm from the params and stick with it.
        bm <- as.numeric(params["bm"])
        if (!let_R0_vary){
          if (vary == "bm"){
            correction_params <- make_params(base_params = params, change_pars = c("bm"), bm = theVal)
          }
          else{
            ##vary is the gammaratio
            correction_params <- make_params(base_params = params, change_pars = c("gamma1", "bm"), gamma1 = as.numeric(make_params()[["gamma2"]])/theVal, bm = bm)
          }
          correction_beta <- compute_beta(R0 = fixed_R0, params = correction_params, use_unequal_transmission = use_unequal_transmission)
          params <- make_params(base_params = correction_params,
                                change_pars = c("Beta"),
                                Beta = correction_beta)
        }
        else{
          if (vary == "bm"){
            ##compute the Beta so that at vary = fix_R0_at, we hit R0 = R0val
            fix_R0_params <- make_params(base_params = params,
                                         change_pars = c("bm"),
                                         bm = fix_R0_at)
            correction_beta <- compute_beta(R0 = fixed_R0,  use_unequal_transmission = use_unequal_transmission, params = fix_R0_params)
            adjust_params <- make_params(base_params = params,
                                         change_pars = c("Beta"),
                                         Beta = correction_beta)
            adjust_params["bm"] <- theVal
            params <- adjust_params
          }
          else{
            ##compute the Beta so that at vary = fix_R0_at, we hit R0 = R0val
            fix_R0_params <- make_params(base_params = params,
                                         change_pars = c("gamma1"),
                                         gamma1 = params[["gamma2"]]/fix_R0_at)
            correction_beta <- compute_beta(R0 = fixed_R0,  use_unequal_transmission = use_unequal_transmission, params = fix_R0_params)
            adjust_params <- make_params(base_params = params,
                                         change_pars = c("Beta"),
                                         Beta = correction_beta)
            adjust_params["gamma1"] <- adjust_params[["gamma2"]]/theVal
            params <- adjust_params            
          }
        }
        
        theSim <- simulate_solve(params = params, use_unequal_transmission = use_unequal_transmission)
        
        ##Maximize just severe prevalence
        sum_peak <- max(theSim$I2)
        if (vary == "bm"){
          holdertibble <- tibble::tibble("peak" = c(sum_peak),
                                         "$\\reducemild$" = c(theVal),
                                         "$\\R_0$" = c(toString(fixed_R0)))
        }
        else{
          ##vary == gammaratio
          holdertibble <- tibble::tibble("peak" = c(sum_peak),
                                         "$\\frac{\\gammasevere}{\\gammamild}$" = c(theVal),
                                         "$\\R_0$" = c(toString(fixed_R0)))
        }
        
        df <- dplyr::bind_rows(df, holdertibble)
      }
    }
    library("ggplot2")
    library("ggpubr")
    library("latex2exp")
    ##Rename bm column
    if (vary == "bm"){
      p <- ggplot(df, mapping = aes(`$\\reducemild$`, peak, color = factor(`$\\R_0$`, levels = R0Seq))) +
        geom_line(lwd = 1.5, size = 1.2)
      if (use_directlabels){
        p <- p + directlabels::geom_dl(mapping = aes(label = `$\\R_0$`) , method = "last.bumpup")
      }
      else{
      }
      p <- p + theme_bw() +
        theme(legend.position = "none")
      p <- p +  labs(x = "relative transmission of mild cases ($\\frac{\\betamild}{\\betasevere}$)", y = "severe peak prevalence", title = "peak prevalence (severe cases)")
    }
    else{
      ##vary == gammaratio
      p <- ggplot(df, mapping = aes(`$\\frac{\\gammasevere}{\\gammamild}$`, peak, color = factor(`$\\R_0$`, levels = R0Seq))) +
        geom_smooth(lwd = 1.5, size = 1.2, se = FALSE)
      if (use_directlabels){
        p <- p + directlabels::geom_dl(mapping = aes(label = `$\\R_0$`) , method = "last.bumpup")
      }
      else{
        
      }
      p <- p + theme_bw() + theme(legend.position = "none")
      p <- p +  labs(x = "relative length of mild infections ($\\frac{\\gammasevere}{\\gammamild}$)", y = "severe peak prevalence", title = "peak prevalence (severe cases)")
    }
  }
  else{
    ##Do nothing.
  }
  p
}

final_size_formula <- function(R0){
  if (R0 <= 0)
    warning("R0 = ", R0)
  tmp <- emdbook::lambertW(-R0 * exp(-R0))
  return(1 + (1/R0) * tmp)
}
##To work in main_results_plot, order of params needs to be R0Seq,  "vary", params, use_unequal_transmission
severe_attack_rate_plot <- function(
  R0Seq = c(3, 1.5*3, 1.5*1.5*3),
  vary = "rho",
  use_directlabels = FALSE,
  params = make_params(),
  let_R0_vary = FALSE,
  use_unequal_transmission = TRUE,
  to = 1,
  by = 0.01,
  from = 0.01,
  theSeq = seq(from = from,
               to = to,
               by = by),
  gamma_to = 2,
  gamma_seq = seq(from = from,
                  to = gamma_to,
                  by = by)){
  if (vary != "gammaratio"){
    fix_R0_at <- params[[vary]]
  }
  else{
    fix_R0_at <- params[["gamma2"]]/params[["gamma1"]]
  }
  ##Match the names of the other functions (this one uses base_params instead of params) without rewriting
  base_params <- params
  if (vary == "gammaratio"){
    ##new vals
    theSeq <- gamma_seq
  }
  else{
    
  }
  df <- data.frame("ZSevere" = c(), "$\\R_0$"  = c(), "theVal" = c())
  for (fixed_R0 in R0Seq){
    for (rhoVal in theSeq){
      if (!let_R0_vary){
        if (vary == "rho"){
          base_params["rho"] <- rhoVal
        }
        else if (vary == "bm"){
          base_params["bm"] <- rhoVal
        }
        else{
          ##vary == gammaratio
          ##We've already bumped up the iterating sequence of values to consider
          base_params["gamma2"] <- as.numeric(base_params["gamma1"])*rhoVal
        }
        ##now fix beta so that we hit the right R0 value.
        correction_beta <- compute_beta(R0 = fixed_R0, params = base_params, use_unequal_transmission = use_unequal_transmission)
        base_params <- make_params(base_params = base_params,
                                   change_pars = c("Beta"),
                                   Beta = correction_beta)
        ##We don't need to compute R0 since we explicitly set it.
        Z <- final_size_formula(fixed_R0)
        severe_attack_rate <- (1-as.numeric(base_params["rho"]))*Z
      }
      else{
        if (vary == "rho"){
          ##compute the Beta so that at vary = fix_R0_at, we hit R0 = R0val
          fix_R0_params <- make_params(base_params = params,
                                       change_pars = c("rho"),
                                       rho = fix_R0_at)
          correction_beta <- compute_beta(R0 = fixed_R0,  use_unequal_transmission = use_unequal_transmission, params = fix_R0_params)
          adjust_params <- make_params(base_params = params,
                                       change_pars = c("Beta"),
                                       Beta = correction_beta)
          adjust_params["rho"] <- rhoVal
          base_params <- adjust_params
        }
        else if (vary == "bm"){
          fix_R0_params <- make_params(base_params = params,
                                       change_pars = c("bm"),
                                       bm = fix_R0_at)
          correction_beta <- compute_beta(R0 = fixed_R0,  use_unequal_transmission = use_unequal_transmission, params = fix_R0_params)
          adjust_params <- make_params(base_params = params,
                                       change_pars = c("Beta"),
                                       Beta = correction_beta)
          adjust_params["bm"] <- rhoVal
          base_params <- adjust_params
          
        }
        else{
          ##vary == gammaratio
          fix_R0_params <- make_params(base_params = params,
                                       change_pars = c("gamma1"),
                                       gamma1 = params[["gamma2"]]/fix_R0_at)
          correction_beta <- compute_beta(R0 = fixed_R0,  use_unequal_transmission = use_unequal_transmission, params = fix_R0_params)
          adjust_params <- make_params(base_params = params,
                                       change_pars = c("Beta"),
                                       Beta = correction_beta)
          adjust_params["gamma1"] <- adjust_params[["gamma2"]]/rhoVal
          base_params <- adjust_params
        }
        Z <- final_size_formula(compute_R0(params = adjust_params, use_unequal_transmission = use_unequal_transmission))
        severe_attack_rate <- (1-as.numeric(base_params["rho"]))*Z
      }
      tempdf <- data.frame("ZSevere" = c(severe_attack_rate), "$\\R_0$" = c(fixed_R0), "theVal" = c(rhoVal))
      df <- dplyr::bind_rows(df, tempdf)
    }
  }
  if (vary == "rho"){
    colnames(df) <- c("Zsevere", "$\\R_0$", "$\\mildprop$")
  }
  else if (vary == "bm"){
    ##vary == bm
    colnames(df) <- c("Zsevere", "$\\R_0$",  "$\\frac{\\betamild}{\\betasevere}$")
  }
  else{
    #vary == "gammaratio"
    colnames(df) <- c("Zsevere", "$\\R_0$", "$\\frac{\\gammasevere}{\\gammamild}$")
  }
  library("ggplot2")
  library("ggpubr")
  library("tibble")
  if (vary == "bm"){
    p <- ggplot(df, mapping = aes(`$\\frac{\\betamild}{\\betasevere}$`, Zsevere, color = factor(`$\\R_0$`, levels = R0Seq)))
    if (use_directlabels){
      p <- p +  directlabels::geom_dl(mapping = aes(label = `$\\R_0$`) , method = "last.bumpup")
    }
    else{
    }
    p <- p + geom_line(lwd = 1.5, size = 1.2, aes(group = `$\\R_0$`)) +
      theme_bw() +
      theme(legend.position = "none")
    p <- p +  labs(x = "relative transmission of mild cases ($\\frac{\\betamild}{\\betasevere}$)", y = "severe attack rate", title = "attack rate (severe cases)")
  }
  else if (vary == "rho"){
    p <- ggplot(df, mapping = aes(`$\\mildprop$`, Zsevere, color = factor(`$\\R_0$`, levels = R0Seq)))
    if (use_directlabels){
      p <- p + directlabels::geom_dl(mapping = aes(label = `$\\R_0$`) , method = "last.bumpup")
    }
    p <- p + geom_line(lwd = 1.5, size = 1.2, aes(group = `$\\R_0$`)) +
      theme_bw() +
      theme(legend.position = "none")
    p <- p +  labs(x = "probability of mild infection ($\\mildprop$)", y = "severe attack rate", title = "attack rate (severe cases)")
  }
  else{
    ##vary == gammaratio
    p <- ggplot(df, mapping = aes(`$\\frac{\\gammasevere}{\\gammamild}$`, Zsevere, color = factor(`$\\R_0$`, levels = R0Seq))) +
      geom_line(lwd = 1.5, size = 1.2, aes(group = `$\\R_0$`))
    if (use_directlabels){
      p <- p + directlabels::geom_dl(mapping = aes(label = `$\\R_0$`) , method = "last.bumpup")
    }
    else{
    }
    p <- p +  theme_bw() +
      theme(legend.position = "none")
    p <- p +  labs(x = "relative length of mild infections ($\\frac{\\gammasevere}{\\gammamild}$)", y = "severe attack rate", title = "attack rate (severe cases)")
  }
  p
}
R0_plot <- function(
  vary = "rho",
  params = make_params(),
  R0Seq = c(3, 1.5*3, 1.5*1.5*3, 4*1.5*1.5*3),
  R0labs = c("WT", "Alpha", "Delta", "Omicron"),
  use_unequal_transmission  = TRUE,
  from = 0,
  to = 1,
  by = 0.01,
  theSeq = seq(from = from,
               to = to,
               by = by),
  gamma_to = 2,
  gammaSeq = seq(from = from,
                 to = gamma_to,
                 by = by),
  legendSize = 0.1,
  do_sqrt = FALSE){
  names(R0Seq) <- R0labs
  if (vary == "gammaratio"){
    ###If we're varying the gammaratio, switch to using the sequence with higher values to iterate over.
    theSeq <- gammaSeq
  }
  else{
  }
  if (vary != "gammaratio"){
    fix_R0_at <- params[[vary]]
  }
  else{
    fix_R0_at <- params[["gamma2"]]/params[["gamma1"]]
  }
  df <- data.frame("Variant" = c(), "$\\R_0$" = c(), "theVal" = c())
  for (fixed_R0 in R0Seq){
    for (rhoVal in theSeq){
      if (vary == "rho"){
        ##compute the Beta so that at vary = fix_R0_at, we hit R0 = R0val
        fix_R0_params <- make_params(base_params = params,
                                     change_pars = c("rho"),
                                     rho = fix_R0_at)
        correction_beta <- compute_beta(R0 = fixed_R0,  use_unequal_transmission = use_unequal_transmission, params = fix_R0_params)
        adjust_params <- make_params(base_params = params,
                                     change_pars = c("Beta"),
                                     Beta = correction_beta)
        adjust_params["rho"] <- rhoVal
      }
      else if (vary == "bm"){
        fix_R0_params <- make_params(base_params = params,
                                     change_pars = c("bm"),
                                     bm = fix_R0_at)
        correction_beta <- compute_beta(R0 = fixed_R0,  use_unequal_transmission = use_unequal_transmission, params = fix_R0_params)
        adjust_params <- make_params(base_params = params,
                                     change_pars = c("Beta"),
                                     Beta = correction_beta)
        adjust_params["bm"] <- rhoVal
        
      }
      else{
        ##vary == gammaratio
        fix_R0_params <- make_params(base_params = params,
                                     change_pars = c("gamma1"),
                                     gamma1 = params[["gamma2"]]/ fix_R0_at)
        correction_beta <- compute_beta(R0 = fixed_R0,  use_unequal_transmission = use_unequal_transmission, params = fix_R0_params)
        adjust_params <- make_params(base_params = params,
                                     change_pars = c("Beta"),
                                     Beta = correction_beta)
        adjust_params["gamma1"] <- adjust_params[["gamma2"]]/rhoVal
      }
      R0 <- compute_R0(params = adjust_params, use_unequal_transmission = use_unequal_transmission)
      tempdf <- data.frame("Variant" = c(names(R0Seq)[R0Seq == fixed_R0]), "$\\R_0$" = c(R0), "theVal" = c(rhoVal))
      df <- dplyr::bind_rows(df, tempdf)
    }
  }
  library("ggplot2")
  library("tibble")
  if (vary == "rho"){
    ##use a tibble for proper column names, easier this way.
    colnames(df) <- c("Variant", "$\\R_0$", "$\\mildprop$")
    if (do_sqrt){
      df$"$\\R_0$" <- sqrt(df$"$\\R_0$")
    }
    else{
    }
    p <- ggplot2::ggplot(df, mapping = aes(x = `$\\mildprop$`, y = `$\\R_0$`, color = Variant)) + ##factor lets us order the entries in the legend
      geom_line(lwd = 1.5, size = 1, aes(group = as.factor(Variant)))
    ##Apply ggppubr publication style theme
    p <- p + labs(title = "basic reproduction number ($\\R_0$)", x = "probability of mild infection $(\\mildprop)$")
  }
  else if (vary == "bm"){
    colnames(df) <- c("Variant", "$\\R_0$", "$\\reducemild$")
    if (do_sqrt){
      df$"$\\R_0$" <- sqrt(df$"$\\R_0$")
    }
    else{
    }
    p <- ggplot2::ggplot(df, mapping = aes(x = `$\\reducemild$`, y = `$\\R_0$`, color = Variant))  +
      geom_line(lwd = 1.5, size = 1, aes(group = as.factor(Variant)))
    ##Apply ggppubr publication style theme
    p <- p + labs(title = "basic reproduction number ($\\R_0$)", x = "relative transmission of mild cases ($\\frac{\\betamild}{\\betasevere}$)")
  }
  else if (vary == "gammaratio"){
    colnames(df) <- c("Variant", "$\\R_0$", "$\\frac{\\gammasevere}{\\gammamild}$")
    if (do_sqrt){
      df$"$\\R_0$" <- sqrt(df$"$\\R_0$")
    }
    else{
    }
    p <- ggplot2::ggplot(df, mapping = aes(x = `$\\frac{\\gammasevere}{\\gammamild}$`, y = `$\\R_0$`, color = Variant)) +
      geom_line(lwd = 1.5, size = 1, aes(group = as.factor(Variant)))
    ##Apply ggppubr publication style theme
    p <- p  + labs(title = "basic reproduction number ($\\R_0$)", x = "relative length of mild infections ($\\frac{\\gammasevere}{\\gammamild}$)")
  }
  else{
  }
  ##legend time
  p <- theme_legend_manager(p = p, vary = vary, legendSize = legendSize, nothing = FALSE)
  return(p)
}

##Compute SIR exact peak prevalence as a function of R0
SIR_compute_exact_peak_prevalence <- function(R0){
  peak <- 1-(1/R0)*(1+log(R0))
  return(peak)
}

peak_R0_plot <- function(use_rel_error = FALSE,
                         ##Vary can be rho or R0, depending on what we want to show on the x-axis
                         vary = "rho"){
  rhoSeq <- seq(from = 0, to = 1, by = 0.01)
  ##Limits of R0 should be between R0m and R0s
  variolation_peaks <- sapply(rhoSeq, function(rho){
    params <- make_params(mu = 0, nu = 0)
    params["rho"] <- rho
    ##Small by argument avoids spikes in the plots.
    sim <- simulate_solve(by = 0.1, params = params, use_unequal_transmission = TRUE)
    peak <- max(sim$I1 + sim$I2)
    peak
  })
  variolation_R0s <- sapply(rhoSeq, function(rho){
    params <- make_params(mu = 0, nu = 0)
    params["rho"] <- rho
    compute_R0(params, use_unequal_transmission = TRUE)
  })
  SIR_peaks <- sapply(variolation_R0s, FUN = SIR_compute_exact_peak_prevalence)
  if (vary == "R0"){
    library("tibble")
    library("ggplot2")
    library("ggpubr")
    if (!use_rel_error){
      df <- tibble("Curve" = c(rep("SIR Standard Formula", length(variolation_R0s)), rep("Variolation Model", length(variolation_R0s))), "Value" = c(SIR_peaks, variolation_peaks), "$\\R_0$" = c(variolation_R0s, variolation_R0s))
      p <-  ggplot(data = df, mapping = aes(x = `$\\R_0$`, y = Value, color = Curve)) +
        geom_line(lwd = 1.5) + labs(y = "peak prevalence", color = "Model", title = "peak prevalence approximation: accuracy") + theme_pubr()
    }
    else{
      ##Use relative error
      df <- tibble("Relative Error" = (SIR_peaks - variolation_peaks)/variolation_peaks, "$\\R_0$" = variolation_R0s)
      p <-  ggplot(data = df, mapping = aes(x = `$\\R_0$`, y = `Relative Error`)) +
        geom_line(lwd = 1.5) + theme_pubr() + labs(y = "relative error", title = "peak prevalence approximation: relative error")
    }
  }
  else{
    library("tibble")
    library("ggplot2")
    library("ggpubr")
    if (!use_rel_error){
      df <- tibble("Curve" = c(rep("approximation (using standard SIR formula)", length(variolation_R0s)), rep("exact (numerical solution of variolation model)", length(variolation_R0s))), "Value" = c(SIR_peaks, variolation_peaks), "$\\mildprop$" = c(rhoSeq, rhoSeq))
      p <-  ggplot(data = df, mapping = aes(x = `$\\mildprop$`, y = Value, color = Curve)) +
        geom_line(lwd = 1.5) + theme_pubr() + labs(y = "peak prevalence", color = "", title = "peak prevalence approximation: accuracy")
      p <- p +  theme(legend.direction = "vertical", legend.margin = margin(-15,0,0,0))
    }
    else{
      ##Use relative error
      df <- tibble("Relative Error" = (SIR_peaks - variolation_peaks)/variolation_peaks, "$\\mildprop$" = rhoSeq)
      p <-  ggplot(data = df, mapping = aes(x = `$\\mildprop$`, y = `Relative Error`)) +
        geom_line(lwd = 1.5) + theme_pubr() + labs(y = "relative error", title = "peak prevalence approximation: relative error")
    }
  }
  p
}

prevalence_time_plot <- function(rho = as.numeric(make_params()["rho"]),
                                 gamma1 = as.numeric(make_params()["gamma1"]),
                                 gamma2 = as.numeric(make_params()["gamma2"]),
                                 total_prevalence_version = FALSE,
                                 params = make_params(rho = rho, gamma1 = gamma1, gamma2 = gamma2),
                                 length = 158,
                                 by = 1,
                                 gamma_I_sum_over_FOI_curve_colour = "#FF8300",
                                 Imild_ratio_curve_colour = "#000000",
                                 showS = TRUE,
                                 dashedHorizontalLines = TRUE,
                                 showItot = FALSE,
                                 line_ending = 170,
                                 Itot_colour = "#F90505",
                                 time = seq(from = 0, to = length, by = by) ,
                                 use_unequal_transmission = FALSE,
                                 sim = simulate_solve(params = params, length = length, by = by, use_unequal_transmission = use_unequal_transmission)){
  
  yint <- 1/compute_R0(params = params, use_unequal_transmission = use_unequal_transmission)
  library("ggpubr")
  library("ggplot2")
  if (!total_prevalence_version){
    if (!use_unequal_transmission){
      ##The formula for the coefficient we're checking uses bm, set it to one if we're using equal transmission to avoid using the default one in make_params.
      bm <- params["bm"] <- 1
    }
    else{
    }
    ##To end the curve where we want it
    factordf <- tibble::tibble("Curve" = "$\\frac{\\gammamild \\Imild + \\gammasevere \\Isevere}{\\betamild \\Imild + \\betasevere \\Isevere}$", "Value" = (params["gamma1"] * sim$I1 + params["gamma2"]*sim$I2)/(params["Beta"]*bm*sim$I1 + params["Beta"]*sim$I2), "Time" = time)
    imdf <- tibble::tibble("Curve" = "Variolated", "Value" = sim$I1, "Time" = time)
    isdf <- tibble::tibble("Curve" = "Severely Infectious", "Value" = sim$I2, "Time" = time)
    ifracdf <- tibble::tibble("Curve" = "$\\frac{\\Imild}{\\Imild + \\Isevere}$", "Value" = sim$I1/(sim$I1 + sim$I2), "Time" = time)
    tempdf <- dplyr::bind_rows(factordf, imdf, isdf, ifracdf)
    Sdf <- tibble::tibble("Curve" = "$S$", "Value" = sim$S, "Time" = time)
    tempdf <- dplyr::bind_rows(tempdf, Sdf)
    if (showItot){
      Itotdf <- tibble::tibble("Curve" = "Total Prevalence", "Value" = sim$I1 + sim$I2, "Time" = time)
      tempdf <- dplyr::bind_rows(tempdf, Itotdf)
    }
    else{
    }
    if (showItot){
      thePal <- c(gamma_I_sum_over_FOI_curve_colour, Imild_ratio_curve_colour,  all_curve_sim_palette[4], all_curve_sim_palette[3], Itot_colour, all_curve_sim_palette[1])
    }
    else{
      thePal <- c(gamma_I_sum_over_FOI_curve_colour, Imild_ratio_curve_colour,  all_curve_sim_palette[4], all_curve_sim_palette[c(3,1)])
    }
    p <- ggplot(data = tempdf, mapping = aes(x = Time, y = Value, color = Curve)) + scale_x_continuous(expand = c(0,0), limits = c(0, 180)) +
      scale_color_manual(values = thePal)
    ##Place the horizontal lines first so they're at the back
    if (!dashedHorizontalLines){
      p <- p + geom_segment( x = 0, y = params["rho"], xend = line_ending, yend = params["rho"], color = "#808080", lwd = 1.5) +
        geom_segment( x = 0, y = yint, xend = line_ending, yend = yint, color = "#808080", lwd = 1.5)
    }
    else{
      ##Dash the horizontal lines
      p <- p + geom_segment(linetype ="dashed", x = 0, y = params["rho"], xend = line_ending, yend = params["rho"], color = "#808080", lwd = 1) +
        geom_segment(linetype ="dashed", x = 0, y = yint, xend = line_ending, yend = yint, color = "#808080", lwd = 1)
    }
    p <- p + geom_line(lwd = 1.5, size = 2) +
      labs(title = "Time Series and Peak Prevalence Plot",
           x = "Time",
           y = "Proportion of population") +
      theme_pubr() +
      ##Labels
      geom_text(aes(line_ending + 3, params["rho"], label = "$\\mildprop$"), color =  "#808080")
    yint <- 1/compute_R0(params = params, use_unequal_transmission = use_unequal_transmission)
    ##Labels
    p <- p + geom_text(aes(line_ending + 3, yint, label = "$\\frac{1}{\\R_0}$"), color =  "#808080")
    p
  }
  else{
    ##Make an alternate plot
    totaldf <- tibble::tibble("Curve" = "$\\Imild + \\Isevere$", "Value" = sim$I1 + sim$I2, "Time" = time)
    if (!use_unequal_transmission){
      ##The formula for the coefficient we're checking uses bm, set it to one if we're using equal transmission to avoid using the default one in make_params.
      params["bm"] <- 1
    }
    else{
    }
    approx_df <- tibble::tibble("Curve" = "$\\frac{1}{\\R_0} (\\frac{\\betamild}{\\gammamild}\\frac{\\Imild}{N} + \\frac{\\betasevere}{\\gammasevere}\\frac{\\Isevere}{N})$", "Value" = yint*((params["Beta"]*params["bm"]/params["gamma1"])*(sim$I1) + (params["Beta"]/params["gamma2"])*(sim$I2)), "Time" = time)
    tempdf <- dplyr::bind_rows(totaldf, approx_df)
    p <- ggplot(data = tempdf, mapping = aes(x = Time, y = Value, color = Curve)) +
      geom_line(lwd = 1.5,size = 2) +
      labs(title = paste0("$\\mildprop = ", rho, "$, ", "$\\frac{1}{\\gammamild} = ", 1/gamma1, "$, ", "$\\frac{1}{\\gammasevere} = ", 1/gamma2, " $ days"),
           x = "Time",
           y = "Proportion of population") + 
      theme_pubr()
    p
  }
}

##Make a single panel of the useful plot.
peak_approximation_plot <- function(gamma_multiplier = 1.5,
                                    show_labels = FALSE,
                                    length = 350,
                                    ##Make curves smooth. The time step shouldn't be too far away from m_by or we'll end up with spiky lines
                                    by = 0.01,
                                    showS = TRUE,
                                    time = seq(from = 0, to = length, by = by),
                                    use_unequal_transmission = TRUE,
                                    ##m_by to be <= 0.005 so that we get close enough to 1 at the right endpoint to hit a relative error of 0.
                                    m_by = 0.005,
                                    ##Change the axis limits
                                    ylims = c(0.8, 1.6),
                                    ##For title and m vertical shift above the y axis limits set by ylims
                                    shift = 0.6,
                                    multi_overlay = FALSE,
                                    gammaratiovals = c(0.25, 0.5, 2, 4),
                                    justPeak = FALSE){
  library("ggplot2")
  library("ggpubr")
  library("directlabels")
  if (!multi_overlay){
    gamma1 <- as.numeric(make_params()["gamma1"]) 
    gamma2 <- gamma_multiplier*gamma1
    if (justPeak){
      ##If we're just storing the peaks, then we don't need a time column
      bigdf <- tibble::tibble("Value" = c(), "Curve" = c(), "$\\mildprop$" = c())
    }
    else{
      bigdf <- tibble::tibble("Time" = c(), "Value" = c(), "Curve" = c(), "$\\mildprop$" = c())
    }
    peaks <- c(x = c(), y = c())
    verticals <- c()
    mVal <- 0
    while (mVal <= 1){
      ##no vital dynamics.
      params <- make_params(rho = mVal, gamma1 = gamma1, gamma2 = gamma2, mu = 0, nu = 0)
      sim <- simulate_solve(params = params, length = length, by = by, use_unequal_transmission = use_unequal_transmission)
      yint <- 1/compute_R0(params = params, use_unequal_transmission = use_unequal_transmission)
      totaldf <- tibble::tibble("Curve" = "$\\Imild + \\Isevere$", "Value" = sim$I1 + sim$I2, "Time" = time)
      if (!use_unequal_transmission){
        ##The formula for the coefficient we're checking uses bm, set it to one if we're using equal transmission to avoid using the default one in make_params.
        params["bm"] <- 1
      }
      else{
      }
      approx_df <- tibble::tibble("Curve" = "$\\frac{1}{\\R_0} (\\frac{\\betamild}{\\gammamild}\\frac{\\Imild}{N} + \\frac{\\betasevere}{\\gammasevere}\\frac{\\Isevere}{N})$", "Value" = yint*((params["Beta"]*params["bm"]/params["gamma1"])*(sim$I1) + (params["Beta"]/params["gamma2"])*(sim$I2)), "Time" = time)
      peakVal <- ((approx_df$Value - totaldf$Value)/totaldf$Value)[which.max(totaldf$Value)]
      peakTime <- totaldf[[which.max(totaldf$Value), "Time"]]
      peaks <- dplyr::bind_rows(peaks, data.frame("x" = peakTime, "y" = peakVal))
      ##May not be exact, minimize the distance to find the best guess.
      Stime <- time[which.min(abs(sim$S - yint))]
      verticals <- c(verticals, Stime)
      if (!justPeak){
        red_to_blue_df <- tibble::tibble("Value" = (approx_df$Value - totaldf$Value)/totaldf$Value, "Curve" = "Relative error", "Time" = time, "$\\mildprop$" = mVal)
      }
      else{
        ##Just store the peak value
        red_to_blue_df <- tibble::tibble("Value" = peakVal, "Curve" = "$\\frac{\\frac{1}{\\R_0} (\\frac{\\betamild}{\\gammamild}\\frac{\\Imild}{N} + \\frac{\\betasevere}{\\gammasevere}\\frac{\\Isevere}{N})}{\\Imild + \\Isevere}$", "$\\mildprop$" = mVal)
      }
      bigdf <- dplyr::bind_rows(bigdf, red_to_blue_df)
      mVal <- mVal + m_by
    }
    verticals <- data.frame("x" = verticals)
    ##For common multipiliers, use more natural english for the titles
    #  if (gamma_multiplier == 1){
    #    gamma_multiplier <- "the same"
    #  }
    #  else if (gamma_multiplier == 2){
    #    gamma_multiplier <- "twice"
    #  }
    #  else if (gamma_multiplier > 2){
    #    gamma_multiplier <- paste0(xfun::numbers_to_words(gamma_multiplier), " times")
    #  }
    #  else if (gamma_multiplier == 1/2){
    #    gamma_multiplier <-   "one quarter"
    #  }
    #  else if (gamma_multiplier == 1/3){
    #    gamma_multiplier <-   "one third"
    #  }
    #  else if (gamma_multiplier == 1/4){
    #    gamma_multiplier <-   "one fourth"
    #  }
    #  else if (gamma_multiplier == 1/5){
    #    gamma_multiplier <-  "one fifth"
    #  }
    #  else{
    #  }
    theTitle <- paste0("$T_{\\rm gen,m} / T_{\\rm gen,s} = ", gamma_multiplier, "$")
    if (!justPeak){
      ##Can maybe use to remove garbage points and lines
      ##Shouldn't have to use this anymore, fixed.
      ##The simulation needs to be long enough to detect and grab the peak.
      #   peaks <- as.data.frame(peaks[peaks$x >= 1,])
      #  verticals <- data.frame("x"  = verticals[verticals$x >= 1,])
      p <- ggplot(data = bigdf, mapping = aes(x = Time, y = Value, color = `$\\mildprop$`)) +
        geom_line(lwd = 1.5, size = 1, aes(group = `$\\mildprop$`)) +
        geom_point(data = peaks, mapping = aes(x,y), color = "darkred", size = 2) +
        geom_vline(data = verticals, aes(xintercept = x), lwd = 1, color = "red") +
        labs(x = "Time",
             y = "Relative error") +
        theme_pubr() +
        theme(legend.position = "none") +
        ##Put the title as a geom_text object to save space
        geom_text(aes(40, max(ylims) + shift, label = theTitle), color =  "black") +
        ##Then we set clipping off so that we overwrite this scales commmand with the one above
        ggplot2::coord_cartesian(clip="off", ylim = ylims)
      if (show_labels){
        ##Then we place the direct labels label
        p <- p + directlabels::geom_dl(aes(label = `$\\mildprop$`), method = "last.bumpup")
        ##Then we widen the upper margin, order is (top, right, bottom, and left)
        ##The we place the label for mildprop
        p <-  p + theme(plot.margin = unit(c(2,1,1,1), "lines")) + geom_text(aes(length +3, max(ylims) + shift + 0.2, label = "$\\mildprop$"), color =  "black")
      }
      p
    }
    else{
      p <- ggplot(data = bigdf, mapping = aes(x = `$\\mildprop$`, y = Value)) +
        geom_line(lwd = 1.5,size = 1, color = "red") +
        labs(title = theTitle,
             x = "$\\mildprop$",
             y = "Relative error") + 
        theme_pubr()
      p <- p + lims(y = ylims) + scale_x_continuous(expand = c(0,0), breaks = seq(from = 0, to = 1, by = 0.1))
      p
    }
  }
  else{
    bigdf <- tibble::tibble("Value" = c(), "Curve" = c(), "$\\mildprop$" = c(), "$\\frac{\\gammasevere}{\\gammamild}$" = c())
    peaks <- c(x = c(), y = c())
    verticals <- c()
    for (gammaMult in gammaratiovals){
      ##Fix one gamma to be a multplier of the other
      gamma1 <- as.numeric(make_params()["gamma1"]) 
      gamma2 <- gammaMult*gamma1
      ##If we're just storing the peaks, then we don't need a time column
      mVal <- 0
      while (mVal <= 1){
        params <- make_params(rho = mVal, gamma1 = gamma1, gamma2 = gamma2, nu = 0, mu = 0)
        sim <- simulate_solve(params = params, length = length, by = by, use_unequal_transmission = use_unequal_transmission)
        yint <- 1/compute_R0(params = params, use_unequal_transmission = use_unequal_transmission)
        totaldf <- tibble::tibble("Curve" = "$\\Imild + \\Isevere$", "Value" = sim$I1 + sim$I2, "Time" = time)
        if (!use_unequal_transmission){
          ##The formula for the coefficient we're checking uses bm, set it to one if we're using equal transmission to avoid using the default one in make_params.
          params["bm"] <- 1
        }
        else{
        }
        ##both the simulation and the approximation are indexed by the same time vecor
        approx_df <- tibble::tibble("Curve" = "$\\frac{1}{\\R_0} (\\frac{\\betamild}{\\gammamild}\\frac{\\Imild}{N} + \\frac{\\betasevere}{\\gammasevere}\\frac{\\Isevere}{N})$", "Value" = yint*(as.numeric(params[["Beta"]]*params[["bm"]]/params[["gamma1"]])*(sim$I1) + as.numeric(params[["Beta"]]/params[["gamma2"]])*(sim$I2)), "Time" = time)
        peakVal <- ((approx_df$Value - totaldf$Value)/totaldf$Value)[which.max(totaldf$Value)]
        red_to_blue_df <- tibble::tibble("Value" = peakVal, "$\\mildprop$" = mVal, "$\\frac{\\gammasevere}{\\gammamild}$" = gammaMult)
        bigdf <- dplyr::bind_rows(bigdf, red_to_blue_df)
        mVal <- mVal + m_by
      }
    }
    theTitle <- "Accuracy of peak prevalence approximation for different $T_{\\rm gen,m} / T_{\\rm gen,s}$"
    p <- ggplot(data = bigdf, mapping = aes(x = `$\\mildprop$`, y = Value, color = `$\\frac{\\gammasevere}{\\gammamild}$`)) +
      geom_line(lwd = 1.5, size = 1, aes(group = as.factor(`$\\frac{\\gammasevere}{\\gammamild}$`))) +
      labs(title = theTitle,
           x = "$\\mildprop$",
           y = "Relative error") + 
      theme_pubr() +
      lims(y = ylims)
    p <- p + geom_dl(method = "last.bumpup", mapping = aes(label  = `$\\frac{\\gammasevere}{\\gammamild}$`))
    p <- p + theme(legend.position = "none")# +
    scale_x_continuous(expand = c(0,0), breaks = seq(from = 0, to = 1, by = 0.1))
    p
  }
}

##Make sure use_unequal_transmission agrees with the call to compute_EE if comparing the output
new_EE_formula <- function(params = make_params(), use_unequal_transmission  = TRUE){
  e_m <- as.numeric(params["mu"]/(params["mu"] + params["gamma1"]))
  e_s <- as.numeric(params["mu"]/(params["mu"] + params["gamma2"]))
  rho <- as.numeric(params["rho"])
  R0 <- compute_R0(params = params, use_unequal_transmission = use_unequal_transmission)
  S <- as.numeric(1/R0)
  if (params[["delta"]] != 0){
    e <- params[["rho"]]*e_m + (1-params[["rho"]])*e_s
    eta <- params[["mu"]]/params[["delta"]]
    I1 <- as.numeric(rho * e_m *(1 - 1/R0)*((eta + 1)/(eta + e)))
    I2 <- as.numeric((1-rho) * e_s *(1 - 1/R0)*(eta + 1)/(eta + e))
    
  }
  else{
    I1 <- as.numeric(rho * e_m *(1 - 1/R0))
    I2 <- as.numeric((1-rho) * e_s *(1 - 1/R0))
  }
  ##Faster than using the params environment
  EE <- c("S" = S, "I1" = I1, "I2" = I2, "R" = 1-S-I1-I2)
  return(EE)
}

##vary can be bm, rho, or gammaratio

##To work in main_results_plot, order of params needs to be R0Seq,  "vary", params, use_unequal_transmission
EE_plot_multiple_R0 <- function(R0seq = c(3, 3*1.5, 3*1.5*1.5),
                                vary = "rho",
                                params = make_params(),
                                use_directlabels = FALSE,
                                use_unequal_transmission = TRUE,
                                theSeq = seq(),
                                let_R0_vary = FALSE,
                                rho_to = 1,
                                rho_by = 0.01,
                                ##Use the same sequence for rho and bm
                                rhoseq = seq(from = 0.01,
                                             to = rho_to,
                                             by = rho_by),
                                gamma_to = 2,
                                gamma_by = 0.01,
                                ##Add some new values for the gammasequence
                                gammaseq = seq(from = 0.01,
                                               to = gamma_to,
                                               by = gamma_by)){
  library("dplyr")
  if (vary == "gammaratio"){
    ##new sequence
    rhoseq <- gammaseq
  }
  else{
  }
  if (vary != "gammaratio"){
    fix_R0_at <- params[[vary]]
  }
  else{
    fix_R0_at <- params[["gamma2"]]/params[["gamma1"]]
  }
  df <- data.frame("$\\R_0$" = c(), theVal = c(), "EE" = c())
  for (R0 in R0seq){
    for (theVal in rhoseq){
      if (!let_R0_vary){
        if (vary == "rho"){
          correction_params <- make_params(rho = theVal)
        }
        else if (vary == "bm"){
          correction_params <- make_params(bm = theVal)
        }
        else{
          ##vary == gammaratio
          correction_params <- make_params(gamma2 = as.numeric(params["gamma1"])*theVal)
        }
        ##Fix R0 stuff
        correction_beta <- compute_beta(R0 = R0, params = correction_params, use_unequal_transmission = use_unequal_transmission)
        params <- make_params(base_params = correction_params,
                              change_pars = c("Beta"),
                              Beta = correction_beta)
      }
      else{
        if (vary == "rho"){
          ##compute the Beta so that at vary = fix_R0_at, we hit R0 = R0val
          fix_R0_params <- make_params(base_params = params,
                                       change_pars = c("rho"),
                                       rho = fix_R0_at)
          correction_beta <- compute_beta(R0 = R0,  use_unequal_transmission = use_unequal_transmission, params = fix_R0_params)
          adjust_params <- make_params(base_params = params,
                                       change_pars = c("Beta"),
                                       Beta = correction_beta)
          adjust_params["rho"] <- theVal
          params <- adjust_params
        }
        else if (vary == "bm"){
          fix_R0_params <- make_params(base_params = params,
                                       change_pars = c("bm"),
                                       bm = fix_R0_at)
          correction_beta <- compute_beta(R0 = R0,  use_unequal_transmission = use_unequal_transmission, params = fix_R0_params)
          adjust_params <- make_params(base_params = params,
                                       change_pars = c("Beta"),
                                       Beta = correction_beta)
          adjust_params["bm"] <- theVal
          params <- adjust_params
          
        }
        else{
          ##vary == gammaratio
          fix_R0_params <- make_params(base_params = params,
                                       change_pars = c("gamma1"),
                                       gamma1 = params[["gamma2"]]/fix_R0_at)
          correction_beta <- compute_beta(R0 = R0,  use_unequal_transmission = use_unequal_transmission, params = fix_R0_params)
          adjust_params <- make_params(base_params = params,
                                       change_pars = c("Beta"),
                                       Beta = correction_beta)
          adjust_params["gamma1"] <- adjust_params[["gamma2"]]/theVal
          params <- adjust_params
        }
      }
      EE <- new_EE_formula(params = params, use_unequal_transmission = use_unequal_transmission)
      EE_severe <- as.numeric(EE["I2"])
      temp_df <- data.frame("EE" = c(EE_severe), "$\\R_0$" = c(R0), theVal = c(theVal))
      df <- bind_rows(temp_df, df)
    }
  }
  if (vary == "rho"){
    colnames(df) <- c("EE", "$\\R_0$", "$\\mildprop$")
  }
  else if (vary == "bm"){
    ##vary == bm
    colnames(df) <- c("EE", "$\\R_0$",  "$\\frac{\\betamild}{\\betasevere}$")
  }
  else{
    #vary == "gammaratio"
    colnames(df) <- c("EE", "$\\R_0$", "$\\frac{\\gammasevere}{\\gammamild}$")
  }
  library("ggplot2")
  if (vary == "bm"){
    p <- ggplot(df, mapping = aes(`$\\frac{\\betamild}{\\betasevere}$`, EE, color = factor(`$\\R_0$`, levels = R0seq)))
    if (use_directlabels){
      p <- p + directlabels::geom_dl(mapping = aes(label = `$\\R_0$`) , method = "last.bumpup")
    }
    p  <- p + geom_line(lwd = 1.5, size = 1.2, aes(group = `$\\R_0$`)) +
      theme_bw() +
      theme(legend.position = "none")
    p <- p +  labs(x = "relative transmission of mild cases ($\\frac{\\betamild}{\\betasevere}$)", y = "severe prevalence", title = "equilibrium prevalence (severe cases)")
  }
  else if (vary == "rho"){
    p <- ggplot(df, mapping = aes(`$\\mildprop$`, EE, color = factor(`$\\R_0$`, levels = R0seq)))
    if (use_directlabels){
      p <- p + directlabels::geom_dl(mapping = aes(label = `$\\R_0$`) , method = "last.bumpup")
    }
    p <- p +  geom_line(lwd = 1.5, size = 1.2, aes(group = `$\\R_0$`)) +
      theme_bw() +
      theme(legend.position = "none")
    p <- p +  labs(x = "probability of mild infection ($\\mildprop$)", y = "severe prevalence", title = "equilibrium prevalence (severe cases)")
  }
  else{
    ##vary == gammaratio
    p <- ggplot(df, mapping = aes(`$\\frac{\\gammasevere}{\\gammamild}$`, EE, color = factor(`$\\R_0$`, levels = R0seq))) +
      geom_line(lwd = 1.5, size = 1.2, aes(group = `$\\R_0$`))
    if (use_directlabels){
      p <- p + directlabels::geom_dl(mapping = aes(label = `$\\R_0$`) , method = "last.bumpup")
    }
    p <- p + theme_bw() +
      theme(legend.position = "none")
    p <- p +  labs(x = "relative length of mild infections ($\\frac{\\gammasevere}{\\gammamild}$)", y = "severe prevalence", title = "equilibrium prevalence (severe cases)")
    ##Works better because we have a different range for the data here
    ##Don't use pretty_axis_labels.
    # p <- pretty_axis_labels(p, at = seq(from = 0, to = max(df$EE), by = 0.005))
    return(p)
  }
  ##Do first
  p <- p + ggplot2::coord_cartesian(clip="off") 
  ##Don't call pretty axis labels (breaks the plot)
  # p <- pretty_axis_labels(p, to = 6-04)
  p
}


compute_exact_growth_rate <- function(params = make_params(), use_unequal_transmission = TRUE){
  if (!use_unequal_transmission){
    params["bm"] <- 1
  }
  else{
  }
  with(as.list(c(params)), {
    newBeta <- rho*bm*Beta + (1-rho)*Beta
    r <- (1/2)*(N*newBeta - (gamma1 + gamma2 + 2*mu) + sqrt((N*newBeta + gamma1 - gamma2)^2 - 4*N*rho*bm*Beta*(gamma1 - gamma2)))
    return(r)
  })
}

##start with delta
estimate_R0_from_growth_rate <- function(r_mult = 3.3, base_r =  compute_exact_growth_rate(params  = make_params(base_params = make_params(), change_pars = c("Beta"), Beta = compute_beta(R0 = 6.75, use_unequal_transmission = TRUE, params = make_params()))), params = make_params(), R0Seq = seq(from = 1, to = 100, by = 0.01), eps = 1e-3){
  r = r_mult * base_r
  for (R0 in R0Seq){
    correction_beta <- compute_beta(R0 = R0,  use_unequal_transmission = TRUE, params = make_params())
    adjust_params <- make_params(base_params = make_params(),
                                 change_pars = c("Beta"),
                                 Beta = correction_beta)
    
    estimated_r <- compute_exact_growth_rate(params = adjust_params)
    if (abs(estimated_r - r) < eps){
      return(R0)
    }
  }
}


theme_legend_manager <- function(p, vary = "rho", legendSize, nothing, directlabels_instead = TRUE){
  library("ggplot2")
  library("ggpubr")
  if (directlabels_instead == TRUE){
    nothing  <- TRUE
  }
  if (nothing){
    ##no direct labels or legend instead
    return(p + theme_bw() + theme(legend.position = "none")) ##theme pubr puts the legend back, so order matters
  }
  else{
  }
  if (directlabels_instead == FALSE){
    if (vary == "rho"){
      legendPos <- c(0.92, 0.72)
    }
    else{
      ## move  the legends to the left in the 2nd and 3rd figures
      legendPos <- c(0.08, 0.72)
    }
    ##don't add direct labels and don't remove the legend
    ##We need to call the legend.position theme command after theme_pubr, so it doesn't get obliterated
    ##With no legend, the call does nothing.
    ##With a legend, it puts the legend where we want it to be
    ##Put the legend inside the plot area by specifying coordinates.
    ##(1,1 gives top right)
    p <- p + theme_bw()  + theme(legend.position = legendPos, legend.key.size = unit(legendSize, "cm"), legend.key.width = unit(0.75, 'cm')) ##order matters, see comments asbove.
  }
  else{
  }
  p
}

##To work in main_results_plot, order of params needs to be R0Seq,  "vary", params, use_unequal_transmission
growth_rate_plot <- function(
  R0Seq = c(3, 1.5*3, 1.5*1.5*3),
  vary = "rho",
  legendSize = 0.1,
  nothing = TRUE,
  params = make_params(),
  use_unequal_transmission  = TRUE,
  let_R0_vary = FALSE,
  from = 0,
  to = 1,
  by = 0.01,
  ##Labels need to match the order of R0Seq
  R0labs = c("WT", "Alpha", "Delta", "Omicron"),
  use_doubling_time_instead = TRUE,
  theSeq = seq(from = from,
               to = to,
               by = by),
  gamma_to = 2,
  gammaSeq = seq(from = from,
                 to = gamma_to,
                 by = by)){
  names(R0Seq) <- R0labs
  if (vary != "gammaratio"){
    fix_R0_at <- params[[vary]]
  }
  else{
    fix_R0_at <- params[["gamma2"]]/params[["gamma1"]]
  }
  ##if letting R0 vary, fix the value of vary to set R0 to (by using beta)
  library("ggplot2")
  library("ggpubr")
  library("tibble")
  if (!use_doubling_time_instead){
    return("this function has been rewritten to only plot with doubling times.")
  }
  if (vary == "rho"){
    df <- tibble("$\\mildprop$" = c(), "Variant" = c(), "Doubling time" = c())
    for (R0 in R0Seq){
      for (rhoVal in theSeq){
        if (!let_R0_vary){
          params["rho"] <- rhoVal
          correction_beta <- compute_beta(R0 = R0,  use_unequal_transmission = use_unequal_transmission, params = params)
          adjust_params <- make_params(base_params = params,
                                       change_pars = c("Beta"),
                                       Beta = correction_beta)
        }
        else{
          ##compute the Beta so that at vary = fix_R0_at, we hit R0 = R0val
          fix_R0_params <- make_params(base_params = params,
                                       change_pars = c("rho"),
                                       rho = fix_R0_at)
          correction_beta <- compute_beta(R0 = R0,  use_unequal_transmission = use_unequal_transmission, params = fix_R0_params)
          adjust_params <- make_params(base_params = params,
                                       change_pars = c("Beta"),
                                       Beta = correction_beta)
          adjust_params["rho"] <- rhoVal
        }
        dbl_time <- log(2)/compute_exact_growth_rate(params = adjust_params, use_unequal_transmission = use_unequal_transmission)
        tempdf <- tibble("$\\mildprop$" = c(rhoVal), "Variant" = c(names(R0Seq)[R0Seq == R0]), "Doubling time" = c(dbl_time))
        df <- dplyr::bind_rows(df, tempdf)
      }
    }
    p <- ggplot2::ggplot(data = df, mapping = aes(x = `$\\mildprop$`, y = `Doubling time`, color = factor(Variant, levels = c("WT", "Alpha", "Delta", "Omicron")))) +
      geom_line(lwd = 1.5)
    p <- theme_legend_manager(p, legendSize = legendSize, nothing = nothing)
    p <-  p + labs(title = "initial doubling time", y = "doubling time (days)", x = "probability of mild infection ($\\mildprop$)")
  }
  else if (vary == "gammaratio"){
    theSeq <- gammaSeq
    df <- tibble("$\\frac{\\gammasevere}{\\gammamild}$" = c(), "Variant" = c(), "Doubling time" = c())
    for (R0 in R0Seq){
      for (rhoVal in theSeq){
        if (!let_R0_vary){
          params["gamma2"] <- as.numeric(params["gamma1"])*rhoVal
          correction_beta <- compute_beta(R0 = R0,  use_unequal_transmission = use_unequal_transmission, params = params)
          adjust_params <- make_params(base_params = params,
                                       change_pars = c("Beta"),
                                       Beta = correction_beta)
        }
        else{
          ##compute the Beta so that at vary = fix_R0_at, we hit R0 = R0val
          fix_R0_params <- make_params(base_params = params,
                                       change_pars = c("gamma1"),
                                       gamma1 = params[["gamma2"]]/fix_R0_at)
          correction_beta <- compute_beta(R0 = R0,  use_unequal_transmission = use_unequal_transmission, params = fix_R0_params)
          adjust_params <- make_params(base_params = params,
                                       change_pars = c("Beta"),
                                       Beta = correction_beta)
          adjust_params["gamma1"] <- as.numeric(params[["gamma2"]])/rhoVal
        }
        dbl_time <- log(2)/compute_exact_growth_rate(params = adjust_params, use_unequal_transmission = use_unequal_transmission)
        tempdf <- tibble("$\\frac{\\gammasevere}{\\gammamild}$" = c(rhoVal), "Variant" = c(names(R0Seq)[R0Seq == R0]), "Doubling time" = c(dbl_time))
        df <- dplyr::bind_rows(df, tempdf)
      }
    }
    p <- ggplot2::ggplot(data = df, mapping = aes(x = `$\\frac{\\gammasevere}{\\gammamild}$`, y = `Doubling time`, color = factor(Variant, levels = c("WT", "Alpha", "Delta", "Omicron")))) +
      geom_line(lwd = 1.5)
    p <- theme_legend_manager(p, legendSize = legendSize, nothing = nothing)
    p <-  p + labs(title = "initial doubling time", y = "doubling time (days)", color = "variant", x = "relative length of mild infections ($\\frac{\\gammasevere}{\\gammamild}$)")
  }
  else{
    ##vary == bm
    df <- tibble("$\\frac{\\betamild}{\\betasevere}$" = c(), "Variant" = c(), "Doubling time" = c())
    for (R0 in R0Seq){
      for (rhoVal in theSeq){
        if (!let_R0_vary){
          params["bm"] <- rhoVal
          correction_beta <- compute_beta(R0 = R0,  use_unequal_transmission = use_unequal_transmission, params = params)
          adjust_params <- make_params(base_params = params,
                                       change_pars = c("Beta"),
                                       Beta = correction_beta)
        }
        else{
          ##compute the Beta so that at vary = fix_R0_at, we hit R0 = R0val
          fix_R0_params <- make_params(base_params = params,
                                       change_pars = c("bm"),
                                       bm = fix_R0_at)
          correction_beta <- compute_beta(R0 = R0,  use_unequal_transmission = use_unequal_transmission, params = fix_R0_params)
          adjust_params <- make_params(base_params = params,
                                       change_pars = c("Beta"),
                                       Beta = correction_beta)
          adjust_params["bm"] <- rhoVal
        }
        dbl_time <- log(2)/compute_exact_growth_rate(params = adjust_params, use_unequal_transmission = use_unequal_transmission)
        tempdf <- tibble("$\\frac{\\betamild}{\\betasevere}$" = c(rhoVal), "Variant" = c(names(R0Seq)[R0Seq == R0]), "Doubling time" = c(dbl_time))
        df <- dplyr::bind_rows(df, tempdf)
      }
    }
    p <- ggplot2::ggplot(data = df, mapping = aes(x = `$\\frac{\\betamild}{\\betasevere}$`, y = `Doubling time`, color = factor(Variant, levels = c("WT", "Alpha", "Delta", "Omicron")))) +
      geom_line(lwd = 1.5)
    p <- theme_legend_manager(p,  legendSize = legendSize, nothing = nothing)
    p <-  p + labs(title = "initial doubling time", y = "doubling time (days)", x = "relative transmission of mild cases ($\\frac{\\betamild}{\\betasevere}$)")
  }
  return(p)
}

##Accepts all the relevant params for simulate_solve, as well as complete simulations via the sim argument, which overwrite all passed simulation argument
estimate_growth_rate <- function(params = make_params(),
                                 init = make_inits(useExpandedModel = useExpandedModel),
                                 length = 360,
                                 by = 1,
                                 use_unequal_transmission = TRUE,
                                 sim = simulate_solve(params = params, useExpandedModel = FALSE, length = length, use_unequal_transmission = use_unequal_transmission, by = by, init = init)){
  ##Fit a straight line to log prevalence and estimate the slope
  sim$total <- sim$I1 + sim$I2
  ##Just fit to the data up to the peak
  sim <- sim[1:which.max(sim$total),]
  fit <- lm(data = sim, formula = log(total) ~ time)
  estimated_r <- as.numeric(coef(fit)["time"])
  return(estimated_r)
}

set_bottom_ylim_zero <- function(p, vary = "rho", set_ymax = FALSE, ymax = NA){
  if (!set_ymax){
    ##Setting the upper limit to NA keeps the default ggplot2 selection.
    return(p + ggplot2::ylim(0, NA))
  }
  else if (vary != "rho"){
    ##Do the above but also add an upper y limit
    return(p + ggplot2::ylim(0, ymax))
  }
  else{
    ### retain the tick labels at 0,3,6,12 but have grid lines
    ### at all the integers, with the grid lines at 1,2, 4,5 etc being fainter
    ###  than the grid lines associated with labelled ticks.
    labelSeq <- c(0,3, 6, 9, 12)
    breakSeq <- labelSeq
    minor_break_seq <- c(1,2, 4,5, 7,8, 10,11)
    p <- p + ggplot2::scale_y_continuous(breaks = breakSeq, minor_breaks = minor_break_seq, limits = c(0, ymax)) ##make the minor grid lines fainter
    p <- p + theme(panel.grid.minor.y = element_line(size = 0.2))
    return(p)
  }
}

set_vertical_gridlines <- function(p){
  library("ggplot2")
  return(p + coord_cartesian(clip = "off", expand = c(0,0)) + scale_x_continuous(breaks=seq(from = 0, to = 1, by = 0.2),  minor_breaks =seq(from = 0, to = 1, by = 0.1))) ##use the breaks argument for the labels at 0.2, 0.4,... and use minor breaks to add the labels in increments of 0.1.
}

##by default, if we let R0 vary, beta will be set so that for the varying parameter, when vary = the default value (i.e rho = make_params()[["rho]], R0 will be set to the variant R0)
##To set beta based on a different value of fix_R0_at, pass that value as the default value for vary, which will be thrown out after that, to replace by the iterating sequence
##By default, for values of m, rho will be fixed at a fix_R0_at of 0.
main_results_plot <- function(vary = "rho",
                              params = make_params(),
                              ##Any pallete of three colours will work
                              thePal = three_curve_sim_palette,
                              init = make_inits(useExpandedModel = useExpandedModel),
                              length = 360,
                              gammaratio_bm_R0_yMax = 25,
                              gammaratio_bm_dbl_time_yMax = 11,
                              gammaratio_bm_peak_yMax = 0.40,
                              gammaratio_bm_attack_rate_yMax = 0.41,
                              gammaratio_bm_EE_rate_yMax = 0.065,
                              fix_vary_at = NULL,
                              R0seq = c(3, 1.5*3, 1.5*1.5*3, 28.4),
                              use_unequal_transmission = TRUE,
                              let_R0_vary = FALSE, ##if TRUE, let R0 change as we fix the parameter, and fix the beta so that at the default value of the parameter (rho = make_params()["rho], R0 = 3, 0r 1.5, or 1.5.1.5*3)
                              singleType = NULL, ##If not set to NULL, set to a single plot type (one of "dbl_time", "peak_height", "severe_attack_rate", or "EE") to display a three panel plot for each variable on the x axis for only that variable
                              set_doubling_time_ymax = vary == "rho",
                              doubling_time_ymax = 25,
                              bm_doubling_time_limits = vary == "bm",
                              all_bottom_y_lim_zero = FALSE,
                              do_sqrt = FALSE,
                              sim = simulate_solve(params = params, useExpandedModel = FALSE, length = length, use_unequal_transmission = use_unequal_transmission, by = by)){
  library("patchwork")
  if (!is.null(fix_vary_at)){
    if (vary == "rho" || vary ==  "bm"){
      params[[vary]] <- fix_vary_at
    }
    else{
      ##vary == gammaratio
      params[["gamma2"]] <- params[["gamma1"]]*fix_vary_at
    }
  }
  
  if (is.null(singleType)){
    ##Epidemic doubling time as a function of m.
    ##flip to get colouring right for the R0 plot
    thepal_flipped <- thePal[c(2,3,4,1)]
    p0 <- R0_plot(R0Seq = R0seq, vary = vary, params = params, use_unequal_transmission = use_unequal_transmission, do_sqrt = do_sqrt) +
      coord_cartesian(clip = "off", expand = c(0,0)) +
      scale_color_manual(values = thepal_flipped)
    ## ggplot2::coord_cartesian() is the ggplot equivalent to xpd = NA, so it will make sure nothing, like the direct labels, is cut off by the plotting scale.
    p1 <- growth_rate_plot(R0Seq = R0seq, vary = vary, params = params, use_unequal_transmission = use_unequal_transmission, let_R0_vary = let_R0_vary, use_doubling_time_instead = TRUE)  +
      ggplot2::coord_cartesian(clip="off", expand = c(0,0))   +
      ##As factor in the original ggplot call lets scale_color_manual work.
      scale_color_manual(values = thePal)
    p2 <- peak_R0_bm_plot(R0Seq = R0seq, vary = vary, params = params, maximize = "severe", use_unequal_transmission = use_unequal_transmission, let_R0_vary = let_R0_vary)  +
      ggplot2::coord_cartesian(clip="off", expand = c(0,0))   +
      scale_color_manual(values = thePal)
    ##Attack rate of severe infections as a function of m.
    p3 <-  severe_attack_rate_plot(R0Seq = R0seq, vary = vary, params = params, use_unequal_transmission = use_unequal_transmission, let_R0_vary = let_R0_vary)  +
      ggplot2::coord_cartesian(clip="off", expand = c(0,0))   +
      scale_color_manual(values = thePal) ##in addition to adding grid lines, theme light put the legend back. get rid of it.
    ##Equilibrium prevalence of severe infections as a function of m.
    ##The I1 column gets renamed to a verbose form, so to actually drop it, we have to pass the full name to drop
    p4 <- EE_plot_multiple_R0(R0seq = R0seq, vary = vary, params, use_unequal_transmission = use_unequal_transmission, let_R0_vary = let_R0_vary)  +
      ##Don't coord_cartesian this here (called within the function instead) so that we don'tt overwrite the pretty axis labels command
      scale_color_manual(values = thePal) +
      ggplot2::coord_cartesian(clip = "off", expand = c(0,0))
    if (all_bottom_y_lim_zero){
      ##should only be for rho
      p1 <- set_bottom_ylim_zero(p1, vary = vary, set_ymax = set_doubling_time_ymax, ymax = doubling_time_ymax)
      p2 <- set_bottom_ylim_zero(p2, vary = vary)
      p3 <- set_bottom_ylim_zero(p3, vary = vary)
      p4 <- set_bottom_ylim_zero(p4, vary = vary)
    }
    else if (bm_doubling_time_limits){
      ##the lower ylim for doubling time should be 0 in fig 6
      p1 <- set_bottom_ylim_zero(p1, vary = vary)
    }
    else{
    }
    if (vary == "rho" || vary == "bm"){
      #In the first two figures, change the vertical gridlines to be at
      #0.1, 0.2, 0.3, etc, and make the tick labels at 0, 0.2, 0.4,
      p0 <- set_vertical_gridlines(p0)
      p1 <- set_vertical_gridlines(p1)
      p2 <- set_vertical_gridlines(p2)
      p3 <- set_vertical_gridlines(p3)
      p4 <- set_vertical_gridlines(p4)
    }  
    ##In the R0(m) plot, please make the lower ylim 1 rather than the min value
    ##achieved. Then there will be a tick label at R0=1, which will be helpful.
    ##
    if (vary == "rho"){
      R0_yMax <- 25
      labels_goby <- 5
      ##do this here to override what was set with theme_bw.
    }
    else {
      ##just the same scales for the gammaratio and bm plots
      R0_yMax <- gammaratio_bm_R0_yMax <- 40
      labels_goby <- 6
      ##p1 has been done separetely in set_bottom_ylim_zero for rho.
      if (vary != "rho"){
        p1 <- p1 + scale_y_continuous(breaks = seq(from = 0, to = gammaratio_bm_dbl_time_yMax, by = 2), limits = c(0, gammaratio_bm_dbl_time_yMax))
      }
      p2 <- p2 + scale_y_continuous(breaks = seq(from = 0.05, to = gammaratio_bm_peak_yMax, by = 0.05), limits = c(0.05, gammaratio_bm_peak_yMax))
      p3 <- p3 + scale_y_continuous(breaks = seq(from = 0.31, to = gammaratio_bm_attack_rate_yMax, by = 0.02), limits = c(0.32, gammaratio_bm_attack_rate_yMax))
      p4 <- p4 + scale_y_continuous(breaks = seq(from = 0.03, to = gammaratio_bm_EE_rate_yMax, by = 0.005), limits = c(0.03, gammaratio_bm_EE_rate_yMax))
    }
    ##do p0 here
    #p0 <- p0 + scale_y_continuous(breaks = seq(from = 1, to = R0_yMax, by = 2),  minor_breaks =seq(from = 0, to = 1, by = 1))##use the breaks argument for the labels at 0.2, 0.4,... and use minor breaks to add the labels in increments of 0.1.
    # p0 <- p0 + scale_y_continuous(breaks = seq(from = 0, to = R0_yMax, by = labels_goby), labels = seq(from = 0, to = R0_yMax, by = labels_goby))##use the breaks argument for the labels at 0.2, 0.4,... and use minor breaks to add the labels in increments of 0.1.
    library("ggh4x")
    p0 <- p0 +  theme(legend.position = "none", legend.text = element_text(size = 7))
    if (vary == "rho"){
      p0 <- p0 + scale_y_sqrt(breaks = c(10,20,30), labels = c(10,20,30), minor_breaks = c(5,15,25, 35), guide = "axis_minor")
    }
    else{
      p0 <- p0 + scale_y_sqrt(breaks = c(10,20,30,40), labels = c(10,20,30, 40), minor_breaks = c(5,15,25, 35), guide = "axis_minor", limits = c(NA, 40) )
    }
    ##Match all the colours of the plots
    if (vary == "rho"){
      p3 <- p3 + ylim(0, 1)
    }
    else{
    }
    ##Revised layout for the plots
    ## p0. R0
    ## p1. Initial doubling time
    ## p2. Peak prevalence (severe cases)
    ## p3. Attack rate (severe cases)
    ## p4. Equilibrium prevalence (severe cases)
    p0 <- manual_direct_labels(vary, p0, use_sqrt = TRUE)
    (p0) /
      (p1) /
      (p2) /
      (p3) /
      (p4)
  }
  else{
    ##The order of params is the same in each function, so it doesn't matter which we chose and this will work
    if (singleType == "dbl_time"){
      fun <- growth_rate_plot
    }
    else if (singleType == "peak_height"){
      fun <- peak_R0_bm_plot
    }
    else if (singleType == "severe_attack_rate"){
      fun <- severe_attack_rate_plot
    }
    else{
      ##singletype == EE
      fun <- EE_plot_multiple_R0
    }
    ##Only order matters, don't use names as the different functions used accept different named parameters.
    ##Names match now but do this for safety's sake.
    rhoplot <- fun(R0seq,  "rho", params, use_unequal_transmission)  +
      ggplot2::coord_cartesian(clip="off")   +
      ##As factor in the original ggplot call lets scale_color_manual work.
      scale_color_manual(values = thePal)
    bmplot <- fun(R0seq, "bm",  params,  use_unequal_transmission)  +
      ggplot2::coord_cartesian(clip="off")   +
      scale_color_manual(values = thePal)
    gammaratioplot <- fun(R0seq, "gammaratio", params, use_unequal_transmission)  +
      ggplot2::coord_cartesian(clip="off")   +
      scale_color_manual(values = thePal)
    (rhoplot) /
      (bmplot) /
      (gammaratioplot)
  }
}

check_r_plot <- function(){
  exact_rs <- c()
  approx_rs <- c()
  i <- 1
  R0seq <- seq(from  = 1, to = 9, by  = 0.01)
  while (i <= length(R0seq)){
    R0 <- R0seq[i]
    correction_beta <- compute_beta(R0 = R0,  use_unequal_transmission = TRUE, params = make_params())
    adjust_params <- make_params(base_params = make_params(),
                                 change_pars = c("Beta"),
                                 Beta = correction_beta)
    
    exact_r <- compute_exact_growth_rate(params = adjust_params)
    approx_r <- mean(make_params()["gamma1"], make_params()["gamma2"])*(compute_R0(adjust_params) - 1)
    exact_rs <- c(exact_rs, exact_r)
    approx_rs <- c(approx_rs, approx_r)
    i <- i + 1
  }
  df <- data.frame("exact"  = exact_rs, "approx" = approx_rs, "R0" = R0seq)
  library("tidyr")
  df_long <- pivot_longer(df, cols = !R0)
  library("ggplot2")
  library("ggpubr")
  
  p <- ggplot(df_long, mapping = aes(R0, value, color = name)) + geom_point() + labs(y = "r", title = "validity of formula for r")
  p <- p + theme_pubr()
  p
}

white_out_string_ggplot2 <- function(x, y, variant, p, col="white", textcolour = "black", k = 6, use_sqrt = FALSE){
  library("ggplot2")
  ##we break the ploting functions if we call strwidth inside them
  ##Presaved strwidths and fugde factors
  if (variant == "Omicron"){
    sw = 0.1*k
  }
  else if (variant == "Delta"){
    sw = 0.06383266*k
  }
  else if (variant == "Alpha"){
    sw = 0.07*k
  }
  else{
    ##variant == WT
    sw = 0.04*k
  }
  ##only the vertical scaling is affected by the sqrt y scale.
  if (use_sqrt == FALSE){
    sh = 0.015 ##Only determined by line spacing, 0.031 is a big big
    
  }
  else{
    sh = 0.005
  }
  p <- p + geom_tile(fill = col,
                     color = col,
                     mapping = aes(x = x,
                                   y = y,
                                   width = sw/7.5,
                                   height = sh*90))
  p <- p + geom_text(aes(x,y), color = textcolour, label = variant, size = 3)
  p
}

manual_direct_labels <- function(vary, p0, use_sqrt = FALSE){
  if (use_sqrt == FALSE){
    if (vary == "rho"){
      ##axis scale is the same
      p0 <- white_out_string_ggplot2(0.3, 11.5, variant = "Omicron", textcolour = "darkred", p0, use_sqrt = use_sqrt)
      p0 <- white_out_string_ggplot2(0.35, 3.8, variant = "Alpha", textcolour  = "darkorange", p0, use_sqrt = use_sqrt)
      p0 <- white_out_string_ggplot2(0.2, 5, variant = "Delta", textcolour = "firebrick3", p0, use_sqrt = use_sqrt)
      p0 <- white_out_string_ggplot2(0.1, 2, variant = "WT", textcolour = "black", p0, use_sqrt = use_sqrt)
    }
    else if (vary == "bm"){
      ##axis scale is the same
      p0 <- white_out_string_ggplot2(0.4, 13.2, variant = "Omicron",textcolour = "darkred", p0, use_sqrt = use_sqrt)
      p0 <- white_out_string_ggplot2(0.6, 7, variant = "Delta", textcolour  = "firebrick3", p0, use_sqrt = use_sqrt)
      p0 <- white_out_string_ggplot2(0.7, 5, variant = "Alpha", textcolour = "darkorange", p0, use_sqrt = use_sqrt)
      p0 <- white_out_string_ggplot2(0.8, 3.8, variant = "WT", textcolour = "black", p0, use_sqrt = use_sqrt)
    }
    else{
      ##vary == "gammaratio
      p0 <- white_out_string_ggplot2(0.5, 12.8, variant = "Omicron", textcolour = "darkred", p0, k = 12, use_sqrt = use_sqrt)
      p0 <- white_out_string_ggplot2(1.2, 5.3, variant = "Alpha", textcolour = "darkorange" , p0, k = 12, use_sqrt = use_sqrt)
      p0 <- white_out_string_ggplot2(0.5, 6, variant = "Delta", textcolour = "firebrick3", p0, k = 12, use_sqrt = use_sqrt)
      p0 <- white_out_string_ggplot2(1.5, 4, variant = "WT", textcolour = "black", p0, k = 12, use_sqrt = use_sqrt)
    }
  }
  else{
    if (vary == "rho"){
      ##axis scale is the same
      p0 <- white_out_string_ggplot2(sqrt(0.3), 18, variant = "Omicron", textcolour = "darkred", p0, use_sqrt = use_sqrt)
      p0 <- white_out_string_ggplot2(sqrt(0.35), 3.05, variant = "Alpha", textcolour  = "darkorange", p0, use_sqrt = use_sqrt)
      p0 <- white_out_string_ggplot2(sqrt(0.2), 5, variant = "Delta", textcolour = "firebrick3", p0, use_sqrt = use_sqrt)
      p0 <- white_out_string_ggplot2(sqrt(0.1), 2.3, variant = "WT", textcolour = "black", p0, use_sqrt = use_sqrt)
    }
    else if (vary == "bm"){
      ##axis scale is the same
      p0 <- white_out_string_ggplot2(sqrt(0.4), 30, variant = "Omicron",textcolour = "darkred", p0, use_sqrt = use_sqrt)
      p0 <- white_out_string_ggplot2(sqrt(0.6), 8.7, variant = "Delta", textcolour  = "firebrick3", p0, use_sqrt = use_sqrt)
      p0 <- white_out_string_ggplot2(sqrt(0.7), 6, variant = "Alpha", textcolour = "darkorange", p0, use_sqrt = use_sqrt)
      p0 <- white_out_string_ggplot2(sqrt(0.8), 4, variant = "WT", textcolour = "black", p0, use_sqrt = use_sqrt)
    }
    else{
      ##vary == "gammaratio
      p0 <- white_out_string_ggplot2(sqrt(0.5), 28, variant = "Omicron", textcolour = "darkred", p0, k = 12, use_sqrt = use_sqrt)
      p0 <- white_out_string_ggplot2(sqrt(1.2), 5, variant = "Alpha", textcolour = "darkorange" , p0, k = 12, use_sqrt = use_sqrt)
      p0 <- white_out_string_ggplot2(sqrt(0.5), 7, variant = "Delta", textcolour = "firebrick3", p0, k = 12, use_sqrt = use_sqrt)
      p0 <- white_out_string_ggplot2(sqrt(1.5), 3.5, variant = "WT", textcolour = "black", p0, k = 12, use_sqrt = use_sqrt)
    }
  }
  p0
}
