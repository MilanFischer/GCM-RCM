# Rescaling function
rescale <- function(x, Range){
  new_min = Range[1]
  new_max = Range[2]
  old_min = 0
  old_max = 1
  new_min + (new_max - new_min) * ((x - old_min) / (old_max - old_min))
}

# Function to write linear regression the equation
LM_equation <- function(slope, intercept){
  ifelse(intercept >= 0, 
         sprintf("y = %.2fx + %.2f", slope, intercept), 
         sprintf("y = %.2fx - %.2f", slope, -intercept))
  # sprintf("y = %.2fx \u2013 %.2f", slope, -intercept)) # Here \u2013 is the Unicode character code for en dash
}

XY_range <- function(x, y, round=0.05, offset=0.04){
  X_range <- c(min(x, na.rm=TRUE), max(x, na.rm=TRUE))
  Y_range <- c(min(y, na.rm=TRUE), max(y, na.rm=TRUE))
  
  X_range[1] <- floor(X_range[1]/round)*round
  Y_range[1] <- floor(Y_range[1]/round)*round
  X_range[2] <- ceiling(X_range[2]/round)*round
  Y_range[2] <- ceiling(Y_range[2]/round)*round
  
  x_offs <- (X_range[2] - X_range[1])*offset
  y_offs <- (Y_range[2] - Y_range[1])*offset

  assign(x = "X_range" , value = X_range + c(-x_offs, x_offs), envir = parent.frame())
  assign(x = "Y_range" , value = Y_range + c(-y_offs, y_offs), envir = parent.frame())
}

# Function to convert a single-row data frame to a named numeric vector
df2nam_vec <- function(df) {
  # Check if the input is a data frame
  if (!is.data.frame(df)) {
    stop("Input must be a data frame")
  }
  
  # Check if the data frame has exactly one row
  if (nrow(df) != 1) {
    stop("Data frame must have exactly one row")
  }
  
  # Step 1: Extract the values from the data frame
  values <- as.numeric(df[1, ])
  
  # Step 2: Extract the names of the columns
  names <- colnames(df)
  
  # Step 3: Combine these into a named numeric vector
  named_vector <- setNames(values, names)
  
  return(named_vector)
}