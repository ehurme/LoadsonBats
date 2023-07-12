simulate_foraging_trajectory <- function(body_mass, time_interval, max_food_consumption_rate, time_threshold) {
  # Define initial variables
  start_location <- c(0, 0) # Starting location
  current_location <- start_location # Current location of the animal
  previous_location <- start_location # Previous location of the animal
  trajectory <- data.frame(x = current_location[1], y = current_location[2], time = 0, food_consumption = 0, weight_loss = 0) # Data frame to store the trajectory
  food_consumed <- 0 # Amount of food consumed
  food_required <- 0.3 * body_mass # Total amount of food required
  current_weight <- body_mass # Current weight of the animal
  food_patches <- data.frame(x = runif(100, -50, 50), y = runif(100, -50, 50), food = runif(100, 0.01 * food_required, 0.1 * food_required)) # Randomly generated food patches

  # Simulate foraging trajectory
  for (t in seq(0, 60 * time_threshold, time_interval)) { # Loop for the time threshold or until the animal reaches its food requirements
    # Check if the animal has reached its food requirements
    if (food_consumed >= food_required) {
      break
    }

    # Calculate weight loss
    weight_loss <- 0.1 * current_weight * (time_interval / 60)
    current_weight <- current_weight - weight_loss

    # Check if the animal is able to eat at the current location
    current_food_patch <- food_patches[sqrt((food_patches$x - current_location[1])^2 + (food_patches$y - current_location[2])^2) < 2, ]
    if (nrow(current_food_patch) > 0) { # If there is food at the current location
      food_available <- sum(current_food_patch$food) # Total amount of food available
      food_consumed_at_location <- min(food_available, max_food_consumption_rate * current_weight) # Amount of food consumed at the location
      food_consumed <- food_consumed + food_consumed_at_location
      food_patches <- food_patches[-which(sqrt((food_patches$x - current_location[1])^2 + (food_patches$y - current_location[2])^2) < 2), ] # Remove food from the patch
      current_weight <- current_weight + (food_consumed_at_location / current_weight) * 100 # Update current weight
      current_location <- current_food_patch[sample(nrow(current_food_patch), 1), 1:2] # Move to the new location
    } else { # If there is no food at the current location
      # Move to a random location within a certain distance
      new_location <- c(runif(1, -50, 50), runif(1, -50, 50))
      while (sqrt((new
