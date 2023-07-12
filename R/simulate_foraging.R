body_mass <- 1 # kg
food_consumption_rate <- 0.3 # 30% of body mass
time_limit <- 720 # 12 hours in minutes
memory_capacity <- 100 # maximum number of food locations remembered
weight_loss_rate <- 0.1 # 0.1% of body mass per minute

landscape_size <- 100
food_concentration <- sample(1:10, size = landscape_size^2, replace = TRUE)/100
landscape <- matrix(food_concentration, nrow = landscape_size, ncol = landscape_size, byrow = TRUE)

simulate_foraging <- function(start_location, body_mass, food_consumption_rate, time_limit, memory_capacity, landscape, weight_loss_rate) {

  # initialize variables
  location <- start_location
  time_elapsed <- 0
  food_consumed <- 0
  memory <- matrix(nrow = memory_capacity, ncol = 2)
  memory_count <- 0
  weight_loss <- body_mass * weight_loss_rate / 100

  # foraging loop
  while (food_consumed < body_mass * food_consumption_rate & time_elapsed < time_limit) {

    # choose a random direction
    direction <- sample(1:4, 1)
    if (direction == 1) {
      location[1] <- location[1] - 1
    } else if (direction == 2) {
      location[1] <- location[1] + 1
    } else if (direction == 3) {
      location[2] <- location[2] - 1
    } else if (direction == 4) {
      location[2] <- location[2] + 1
    }

    # check if the animal has found food
    if (landscape[location[1], location[2]] > 0) {
      food_found <- landscape[location[1], location[2]] * body_mass
      food_consumed <- food_consumed + food_found
      landscape[location[1], location[2]] <- 0
      # check if the animal can remember this food location
      if (memory_count < memory_capacity) {
        memory_count <- memory_count + 1
        memory[memory_count, ] <- location
      } else {
        # overwrite the oldest memory location
        memory[1:(memory_capacity-1), ] <- memory[2:memory_capacity, ]
        memory[memory_capacity, ] <- location
      }
    }

    # update time elapsed and weight loss
    time
