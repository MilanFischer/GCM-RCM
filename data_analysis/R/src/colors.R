darken_col <- function(COL, darkening_factor=1){
  
  # Convert the hexadecimal color to RGB
  original_color_rgb <- col2rgb(COL_RCMs)
  
  # Darken the color by reducing the RGB values
  # You can adjust the factor to control how much darker you want the color to be
  
  darker_color_rgb <- original_color_rgb * darkening_factor
  
  # Ensure values are within the valid range [0, 255]
  darker_color_rgb <- pmax(0, pmin(255, darker_color_rgb))
  
  # Convert back to hexadecimal format
  darker_color_hex <- rgb(darker_color_rgb[1], darker_color_rgb[2], darker_color_rgb[3], maxColorValue = 255)
  
  return(darker_color_hex)
  
}


lighten_col <- function(COL, lightening_factor=1){
  # Convert the hexadecimal color to RGB
  original_color_rgb <- col2rgb(COL)
  
  # Lighten the color by increasing the RGB values
  # Adjust the factor to control how much lighter you want the color to be
  lighter_color_rgb <- original_color_rgb + (255 - original_color_rgb) * lightening_factor
  
  # Ensure values are within the valid range [0, 255]
  lighter_color_rgb <- pmax(0, pmin(255, lighter_color_rgb))
  
  # Convert back to hexadecimal format
  lighter_color_hex <- rgb(lighter_color_rgb[1], lighter_color_rgb[2], lighter_color_rgb[3], maxColorValue = 255)
  
  return(lighter_color_hex)
  
}