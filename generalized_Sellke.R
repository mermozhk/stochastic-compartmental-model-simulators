rm(list = ls())
generalized_Sellke <- function(edges,initial,parameters,final_time, z,susceptibles, functions,plot=F){
  
  
  draw <- function(s,j){
    set.seed(s)
    out <- runif(j)
    set.seed(NULL)
    return(out)
  }
  
  edges_reshaped <- matrix(data = edges,ncol=2,byrow = T)
  suscep_trans_compartments <- unique(edges_reshaped[,1])
  
  initialize_variables <- function(which_compartment,initial,z,states,parameters){
    output <- rep(NA,sum(initial$initial_value))
    index <- which(initial$compartment==which_compartment)
    if(initial$initial_value[index] >0){
      U <- draw(s = z[index],j=initial$initial_value[index])
      q_dtrib_name <- initial$distributions[index]
      output[1:initial$initial_value[index]] <- sort(do.call(what = q_dtrib_name,args = list(u=U,parameters=parameters)),na.last = T)
    }
    output
  }
  compute_pressure_derivatives <- function(what_susceptible,functions,states,parameters){
    index <- which(edges_reshaped[,1]==what_susceptible)
    pres <- sapply(X=1:length(index),FUN = function(i,functions,states,parameters) do.call(what = functions[i],args = list(states=states,parameters=parameters)),functions=functions,states=states,parameters=parameters)
    sum(pres)
  }
  
  next_exit_times_function <- function(which_compartment,instant,initial_variables,pressure_levels,states,parameters){
   
    colum_names <- colnames(initial_variables)
    index_2 <- which(colum_names==which_compartment)
    #print(index_2)
    cat("variable=",initial_variables[1,index_2],"\n")
    if(which_compartment%in% susceptibles){
      index <- which(susceptibles==which_compartment)
      denom <- compute_pressure_derivatives(what_susceptible = which_compartment,functions=functions,states=states,parameters=parameters)
      resistance <- initial_variables[1,index_2]
      if(is.na(resistance)){
        return(Inf)
      }else{
        return(instant + (resistance-pressure_levels[[length(pressure_levels)]][index])/denom)
      }
      
      
    }else{
      if(is.na(initial_variables[1,index_2])){
        return(Inf)
      }else{
        return(initial_variables[1,index_2])
      }
      
    }
    
  }
  
  
  
  pressure_levels <- list(rep(0,length(susceptibles)))
  
  compartments <- initial$compartment
  
  
  t <- 0
  current_state <- initial$initial_value
  Tps <- list(t)
  output <- list(current_state)
  
  i <- 1   
  
  initial_variables <- as.data.frame(sapply(X = suscep_trans_compartments,FUN = initialize_variables,initial=initial,z=z,states=current_state,parameters=parameters))
  
  while(t<final_time){
    cat("iteration",i,"\n")
    i <- i+1
    
    #head(initial_variables)
    print(current_state)
    
    next_exit_times <- matrix(sapply(X = suscep_trans_compartments,FUN = next_exit_times_function,instant=t,initial_variables=initial_variables,pressure_levels=pressure_levels,states=current_state,parameters=parameters),nrow=1)
    
    colnames(next_exit_times) <- suscep_trans_compartments
    next_exit_times <- as.data.frame(next_exit_times)
    next_exit_event_index <- which.min(next_exit_times[1,])
    next_event_time <- next_exit_times[1,next_exit_event_index]
    
    #print("exit times")
    #print(next_exit_times)
    
    if(identical(next_event_time,Inf)){
      y <- as.data.frame(matrix(unlist(output),ncol=length(compartments),byrow = T))
      colnames(y) <- compartments
      return(list(Tps=unlist(Tps),output=y))
    }else{
      
    compartment_index <- names(next_exit_event_index)
    transition_indices <-which(edges_reshaped[,1]==compartment_index)
    possible_transitions <- edges_reshaped[transition_indices,]
    possible_transition_rate_functions <- functions[transition_indices]
    
    rates <- sapply(X = 1:length(transition_indices),FUN = function(i,func) do.call(what = func[i],args = list(states=current_state,parameters=parameters)),func=possible_transition_rate_functions)
    cat("rates=",rates,"\n")
    proba <- rates/sum(rates)
    proba_sorted <- sort(proba)
    U_decision <- draw(s = z[length(z)],j=i)[i]
    next_event_index <- which(proba==proba_sorted[min(which((U_decision<=cumsum(proba_sorted))==TRUE))])
    
    if(is.vector(possible_transitions)){
      involved_compartments <- possible_transitions
    }else{
      involved_compartments <- possible_transitions[next_event_index,]
    }
    
    #print(involved_compartments)
    
    involved_compartment_indices <- sapply(X = 1:length(involved_compartments),FUN = function(i,initial,involved_compartments) which(initial$compartment==involved_compartments[i]),initial=initial,involved_compartments=involved_compartments)
    #print(involved_compartment_indices)
   
    cat("next_time=",next_event_time,"\n")
   
    A <- sapply(X = susceptibles,FUN = compute_pressure_derivatives,functions=functions,states=current_state,parameters=parameters)
    cat("A=",A,"\n")
    pressure_levels[[i]] <- pressure_levels[[i-1]] + (next_event_time-t)*A
      
    current_state[involved_compartment_indices] <- current_state[involved_compartment_indices]+c(-1,+1)
    
    initial_variables[1,which(suscep_trans_compartments==involved_compartments[1])] <- NA
    destination_compartment <- involved_compartments[2]
    if(destination_compartment%in% suscep_trans_compartments){
      destination_index <- which(initial$compartment==destination_compartment)
      u <- draw(s=z[destination_index],j=i)[i]
      if(destination_compartment%in%susceptibles){
       
        new_resistance <- do.call(what = initial$distributions[destination_index],args = list(u=u,parameters=parameters))
        initial_variables[sum(current_state),which(suscep_trans_compartments==destination_compartment)] <- pressure_levels[[i]][which(susceptibles==destination_compartment)]+new_resistance
      }else{
       
        new_duration <- do.call(what = initial$distributions[destination_index],args = list(u=u,parameters=parameters))
        initial_variables[sum(current_state),which(suscep_trans_compartments==destination_compartment)] <- next_event_time+new_duration
      }
    }
    initial_variables <- apply(X = initial_variables,MARGIN = 2,FUN = function(x) sort(x,na.last = T))
    t <- next_event_time
    Tps[[i]] <- t
    output[[i]] <- current_state
    
      }
  
  
  }
  y <- as.data.frame(matrix(unlist(output),ncol=length(compartments),byrow = T))
  colnames(y) <- compartments
  return(list(Tps=unlist(Tps),output=y))
}




parameters <- c(.1,.2,1/30,1/30)

f1 <- function(u,parameters){
  qexp(u,rate = 1)
}

f2 <- function(u,parameters){
  qexp(u,rate = parameters[3])
}

f3 <- function(u,parameters){
  qexp(u,rate = parameters[4])
}

f_na <- function(u,parameters){
  NA
}

pressure_function1 <- function(states,parameters){
  parameters[2]*(sum(states[2:3])*parameters[1])/sum(states)
}
pressure_function2 <- function(states,parameters){
  (1-parameters[2])*(sum(states[2:3])*parameters[1])/sum(states)
}

rate_function1 <- function(states,parameters){
  parameters[3]
}

rate_function2 <- function(states,parameters){
  parameters[4]
}

init <- as.data.frame(c("S","A","I","R"))
colnames(init) <- "compartment"
init$initial_value <- c(10000,0,5,0)
edges <- c("S","I","S","A","I","R","A","R")
init$distributions <- c("f1","f2","f3","f_na")
z <- sample(1:1000,5)
functions <- c("pressure_function1","pressure_function2","rate_function1","rate_function2")
final_time <- 100

y <- generalized_Sellke(edges=edges,initial=init,parameters=parameters,final_time=final_time, z=z,susceptibles="S", functions=functions,plot=F)

par(mfrow=c(2,2))
plot(y$Tps,y$output$S,type = "s")
plot(y$Tps,y$output$A,type = "s")
plot(y$Tps,y$output$I,type = "s")
plot(y$Tps,y$output$R,type = "s")