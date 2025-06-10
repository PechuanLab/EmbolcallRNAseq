# Function: plot_pca_loadings
# Description:
# This function plots the PCA loadings ranked in descending order and identifies inflection points where the second derivative of a fitted spline changes sign. It also allows for the addition of vertical dashed lines and labels for specified variables.

# Parameters:
# - loadings_df: DataFrame containing the PCA loadings. Must include columns Rank and Variable.
# - plot_variables: (Optional) A vector of variable names to be highlighted on the plot with vertical dashed lines and labels.
# - spar: (Optional) Smoothing parameter for the spline fitting. Default is 0.999.

# Returns:
# - inflection_points: A vector of indices where the second derivative of the spline changes sign, indicating inflection points.

plot_pca_loadings <- function(loadings_df, plot_variables = NULL, spar = 0.999) {
  # Sort the loadings in descending order
  loadings_df <- loadings_df %>% arrange(desc(Rank))
  
  # Create a rank array
  loadings_df$RankIndex <- 1:nrow(loadings_df)
  
  # Fit a spline to the data
  spline_fit <- smooth.spline(loadings_df$RankIndex, loadings_df$Rank, spar = spar)
  
  # Calculate the first and second derivatives of the spline
  spline_base <- predict(spline_fit, loadings_df$RankIndex)$y
  first_derivative <- predict(spline_fit, loadings_df$RankIndex, deriv = 1)$y
  second_derivative <- predict(spline_fit, loadings_df$RankIndex, deriv = 2)$y
  
  # Find the inflection points where the second derivative changes sign
  inflection_points <- which(diff(sign(second_derivative)) != 0) + 1
  
  # Create the plot
  p <- ggplot(loadings_df, aes(x = RankIndex, y = Rank)) +
    geom_point(color="#FED16A") +
    geom_line() +
    geom_point(data = loadings_df[inflection_points, ], aes(x = RankIndex, y = Rank), color = "#16610E", size = 3) +
    geom_line(aes(y = spline_base), color = "#F97A00", linetype = "dashed") +
    ggtitle("PCA Loadings Ranked") +
    xlab("Rank") +
    ylab("Loading Value") +
    theme_pubr()
  
  # Add vertical dashed lines and labels for specified variables
  if (!is.null(plot_variables)) {
    for (var in plot_variables) {
      var_index <- which(loadings_df$Variable == var)
      p <- p + geom_vline(xintercept = var_index, color = "#FF3F33", linetype = "dashed") +
        annotate("text", x = var_index, y = max(loadings_df$Rank), label = var, 
                 color = "black", angle = 90, vjust = -0.5)
    }
  }
  
  # Print the plot
  print(p)
  
  # Return inflection points
  return(inflection_points)
}
