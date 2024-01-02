plot_centrality <- function(obj,
                            # accepts overview or density
                            plot_type = "tiefighter",
                            # accepts: instrength, outstrength, strength, density_beta, density_pcor
                            density_type = "outstrength",
                            # only accepts one for now
                            cis = 0.95){
  
  # TODO could add centrality difference
  
  
  if (!is.numeric(cis) || any(cis <= 0) || any(cis >= 1)) {
    stop("cis must be a numeric vector with values between 0 and 1 (exclusive)")
  }
  
  
  #--- Preparation
  # Combine all centrality measures
  create_centrality_df <- function(centrality, suffix) {
    df <- as.data.frame(centrality)
    colnames(df) <- paste0(colnames(df), "_", suffix)
    return(df)
  }
  
  instrength <- create_centrality_df(obj$instrength, "instrength")
  outstrength <- create_centrality_df(obj$outstrength, "outstrength")
  strength <- create_centrality_df(obj$strength, "strength")
 
  
  df_centrality <- cbind(instrength, 
                         outstrength, 
                         strength)
  
  #--- Overview
  if (plot_type == "tiefighter") {
     
    #--- Plot
    overview_plot <- df_centrality %>% 
      tidyr::pivot_longer(cols = everything(),
                          names_to = "measure",
                          values_to = "value") |>
      dplyr::group_by(.data$measure) |>
      dplyr::summarize(mean_value = mean(.data$value, 
                                      na.rm = TRUE),
                    lb = quantile(.data$value, 
                                  probs = (1 - cis) / 2,
                                  na.rm = TRUE),
                    ub = quantile(.data$value,
                                  probs = 1 -((1 - cis) / 2)),
                                  na.rm = TRUE) |>
      dplyr::ungroup() |>
      tidyr::separate(col = measure,
                      into = c("variable", "centrality"),
                      sep = "_") |>
      ggplot2::ggplot(aes(x = .data$mean_value, y = .data$variable)) +  
      ggplot2::geom_point() + 
      ggplot2::geom_errorbarh(aes(
        xmin = lb , 
        xmax = ub)) +
      ggdist::theme_ggdist() + 
      ggplot2::facet_grid(. ~ .data$centrality) + 
      labs(x = "Centrality",
           y = "Variable")

      return(overview_plot)
  
  }
  #--- Density
  if(plot_type == "density"){
    

    #--- Plot
    # browser()
    density_plot <- df_centrality %>% 
      tidyr::pivot_longer(cols = everything(),
                          names_to = "measure",
                          values_to = "value") |>
      dplyr::group_by(.data$measure) |>
      dplyr::mutate(mean_value = mean(.data$value, 
                                         na.rm = TRUE)) |>
      dplyr::ungroup() |>
      tidyr::separate(col = measure,
                      into = c("variable", "centrality"),
                      sep = "_") |>
      ggplot2::ggplot(aes(
                          x = .data$value, 
                          y = .data$variable)) +  
      ggdist::stat_slab(aes(
                        fill = after_stat(.data$level)), 
                        # fixed width for now
                        .width = c(0.8, 0.9, 0.95, 1)) +
      ggdist::stat_pointinterval(aes(x = .data$mean_value), 
                                 size = 1) +
      scale_alpha(guide = "none") +
      ggdist::theme_ggdist() + 
      ggplot2::facet_grid(. ~ .data$centrality) + 
      scale_fill_brewer() + 
      labs(x = "Centrality",
           y = "Variable", 
           fill = "Credible Interval")
    
    return(density_plot)
    
    
    
    
  }
}  
