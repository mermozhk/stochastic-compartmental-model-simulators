rm(list=ls())
gillespie_per_transitions<- function(edges,initial,parameters, rate_functions,final_time, z,plot=F){
  
  draw <- function(s,j){
    set.seed(s)
    out <- runif(j)
    set.seed(NULL)
    return(out)
  }
  
  compartments <- initial$compartment
  
  t <- 0
  current_state <- initial$initial_value
  Tps <- list(t)
  output <- list(current_state)
  
  events <- 1:length(rate_functions)
  i <- 1
  states <- current_state
  rates <- sapply(X=1:length(rate_functions),FUN = function(i,rate_functions,states,parameters) do.call(what=rate_functions[i],args=list(states=states,parameters=parameters)),states=states,parameters=parameters,rate_functions=rate_functions)
  #print(rates)
  while(t<final_time && !identical(rates,rep(0,length(rates)))){
    #print(i)
    if(i==1){
      U <- sapply(z,FUN = draw,j=i)
    }else{
      U <- sapply(z,FUN = draw,j=i)[i,]
    }
    
    i <- i+1
    # cat("i=",i,"\n")
    states <- current_state
    #print(rates/sum(rates))
    #next_event <- sample(x=events,size=1,prob = rates/sum(rates))
    events_time <- -log(U)/rates
    
    next_event <- which.min(events_time)
    next_event_time <- min(events_time)
    involved_compartments <- edges[c(2*next_event-1,2*next_event)]
    #print(involved_compartments)
    compartment_indices <- sapply(X = 1:length(involved_compartments),FUN = function(i,compartments,involved_compartments) which(compartments==involved_compartments[i]),compartments=compartments,involved_compartments=involved_compartments)
    #other_compartments <- sapply(X = 1:length(involved_compartments),FUN = function(i,compartments,involved_compartments) which(compartments!=involved_compartments[i]),compartments=compartments,involved_compartments=involved_compartments)
    #print(compartment_indices)
    current_state[compartment_indices] <- current_state[compartment_indices]+c(-1,+1)
    #current_state[other_compartments] <- 0
    t <- t + next_event_time
    cat("Time=",t,"\n")
    Tps[[i]] <- t
    output[[i]] <- current_state
    #print(current_state)
    rates <- sapply(X=1:length(rate_functions),FUN = function(i,rate_functions,states,parameters) do.call(rate_functions[i][[1]],args=list(states,parameters)),states=current_state,parameters=parameters,rate_functions=rate_functions)
    
  }
  
  if(plot==T){
    require(igraph)
    g <- graph(edges = edges,directed = T)
    plot(g)
  }
  y <- as.data.frame(matrix(unlist(output),ncol=length(compartments),byrow = T))
  colnames(y) <- compartments
  list(Tps=unlist(Tps),output=y)
}

init <- as.data.frame(c("S","I","R"))
colnames(init) <- "compartment"
init$initial_value <- c(10000,1,0)
parameters <- c(.1,1/50)

infection_A <- function(states,parameters){
  states[1]*sum(states[2])*parameters[1]*parameters[2]/sum(states)
}
infection_I <- function(states,parameters){
  states[1]*sum(states[2])*parameters[1]/sum(states)
}

removal_A <-function(states,parameters){
  states[2]*parameters[3]
}
removal_I <-function(states,parameters){
  states[2]*parameters[2]
}
reinfection <-function(states,parameters){
  states[4]*parameters[5]
}

edges <- c("S","I","I","R")
z=sample(x = 1:1000,length(edges)/2)
Tf <- 400
#y <- model(edges = c("S","A","S","I","A","R","I","R"),initial = init,parameters = parameters, rate_functions = c("infection_A","infection_I","removal_A","removal_I"),final_time = 200,plot = T)
y <- gillespie_per_transitions(edges = edges,initial = init,parameters = parameters, rate_functions = c("infection_I","removal_I"),final_time = Tf,z=z,plot = F)
par(mfrow=c(2,2))
plot(y$Tps,y$output$S,type = "l",ylab = "S",xlab = "Time")
plot(y$Tps,y$output$I,type = "l",ylab = "I",xlab = "Time")
plot(y$Tps,y$output$R,type = "l",ylab = "R",xlab = "Time")

