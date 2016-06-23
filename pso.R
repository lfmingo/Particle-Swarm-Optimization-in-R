fitness <- function(pop, calculate_fitness) {
  pop[,dimensions+1] <- calculate_fitness(pop[,1:dimensions]);
  return(pop)
}

pso <- function(verbose, individuals, dimensions, iterations, c1, c2, inertia, calculate_fitness) {
  if (verbose)
    cat ("Startint PSO\n");
  
  #########################
  # Init population (actual + fitness + local_best + fitness)
  population <- array (runif(2*(dimensions+1)*individuals), c(individuals, 2*(dimensions+1)));
  
  # Calculate fitness
  population <- fitness (population, calculate_fitness);
  # Init local_best
  population[,(dimensions+2):(2*(dimensions+1))] <- population[,1:(dimensions+1)];
  # Init global_best
  global_best <- population[max.col(t(population[,dimensions+1])),1:(dimensions+1)];
  
  # Init velocity
  velocity <- array (rnorm(individuals*dimensions), c(individuals, dimensions));
  inertia_coeff <- inertia/iterations;
  
  # Init historic data
  fitness_historic <- array(0, iterations);
  global_best_historic <- array(0, c(iterations,dimensions + 1));
  
  for (iteration in 1:iterations) {
    uniform_coeffs <- runif(individuals*2);
    # Update velocity
    velocity <- inertia * velocity +
      c1 * uniform_coeffs[1:individuals] * (population[,(dimensions+2):(2*(dimensions+1)-1)] - population[,1:dimensions]) +
      c2 * uniform_coeffs[(individuals+1):(2*individuals)] * (global_best[1:dimensions] - population[,1:dimensions]);
    # Update inertia
    inertia <- inertia - inertia_coeff;
    
    # Update population with fitness
    population[,1:dimensions] <- population[,1:dimensions] + velocity;
    population <- fitness (population, calculate_fitness);
    
    # Update local_best
    for (individual in 1:individuals) {
      if (population[individual, dimensions+1] > population[individual, 2*(dimensions+1)]) {
        population[individual,(dimensions+2):(2*(dimensions+1))] <- population[individual,1:(dimensions+1)];
      } 
    }
    
    # Update global_best
    if (max (population[,dimensions+1] > global_best[dimensions+1])) {
      global_best <- population[max.col(t(population[,dimensions+1])),1:(dimensions+1)];
    }
    
    if (verbose)
      cat("Iteration ", iteration, " Fitness ", global_best[dimensions+1], "\n");
    
    # Keep historic fitness data
    fitness_historic[iteration] <- global_best[dimensions+1];
    
    # Keep historic global_best
    global_best_historic[iteration,] <- global_best;
  }
  
  #########################
  if (verbose)
    cat ("Finish PSO\n");
  
  return(list(best=global_best, fitness_evol=fitness_historic, global_best=global_best_historic));
}






