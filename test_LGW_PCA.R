#ANPC PCA quality control visualisation

# FUNC_data = data containing individual lipid data - MUST CONTAIN sampleID column 
# FUNC_colour_by = how to colour the plot (e.g. sample class, is_ltr or cohort)
# FUNC_plot label = what to label the scores plot with (e.g. sampleID)
# FUNC_scaling = UV or Pareto

lgw_pca <- function(FUNC_data, 
                    FUNC_metabolite_list, 
                    FUNC_colour_by, 
                    FUNC_plot_label, 
                    FUNC_scaling,
                    FUNC_title,
                    FUNC_project_colours,
                    FUNC_option_invert_y,
                    FUNC_option_invert_x,
                    FUNC_option_plot_qc){
  
  pca_output <- list()
  
  title_text <- FUNC_title
  
  qc_idx <- which(FUNC_data[["sample_type"]] == "qc")
  
  if(FUNC_option_plot_qc == FALSE){
    FUNC_data <- FUNC_data %>% 
      filter(sample_type != "qc")
  }
  
  #create data matrix for PCA
  pca_x <- FUNC_data %>%  select(all_of(FUNC_metabolite_list)) %>% as.matrix() 
  pca_x[pca_x == 0] <- NA #remove all 0 values
  pca_x[is.infinite(pca_x)] <- NA #remove all infinite values
  min_value <- min(pca_x, na.rm = TRUE) # find the lowest value in the matrix
  pca_x[is.na(pca_x)] <- min_value/100 # replace NA, Inf, and 0 values
  
  pca_x <- log(pca_x + 1) #log values
  
  #create PCA model
  pca_output$pca_model <- pca(pca_x, pc=3, scale = paste(FUNC_scaling), center = TRUE)
  
  # extract score values
  PC1 <- as.numeric(as.matrix(pca_output$pca_model@t[,1]))
  PC2 <- as.numeric(as.matrix(pca_output$pca_model@t[,2]))
  
  # extract loadings values for plotting
  plotly_loadings_data <- pca_output$pca_model@p %>% 
    as_tibble(rownames = "variable") %>% 
    dplyr::rename(PC1 = V1, PC2 = V2)
  
  #create plot data
  pca_colour <- FUNC_data %>% select(all_of(FUNC_colour_by)) %>% pull()
  pca_plot_label <- FUNC_data %>% select(all_of(FUNC_plot_label)) %>% pull()
  
  plot_Val <- tibble(PC1, PC2, pca_colour, pca_plot_label)
  
  #plot settings
  x_axis_settings_scores <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    title = paste("PC1 (", round(pca_output$pca_model@Parameters$R2[1]*100,1), " %)", sep = "")
  )
  
  y_axis_settings_scores <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    title = paste("PC2 (", round(pca_output$pca_model@Parameters$R2[2]*100,1), " %)", sep = "")
  )
  
  pca_output$plot_scores <- plot_ly(type = "scatter", 
                                    mode = "markers", 
                                    data = plot_Val, 
                                    x = ~PC1, 
                                    y = ~PC2, 
                                    text = ~pca_plot_label, 
                                    color = ~pca_colour, 
                                    colors = FUNC_project_colours,
                                    legendgroup = ~pca_colour,
                                    showlegend = NULL,
                                    marker = list(size = 10, 
                                                  opacity = 1,
                                                  line = list(
                                                    color = '#000000',
                                                    width = 1)
                                    )) %>% 
    layout(
      title = paste("PCA - ", title_text, sep = ""),
      xaxis = x_axis_settings_scores,
      yaxis = y_axis_settings_scores,
      margin = list(l = 65, r = 50, b=65, t=85),
      title = paste0(FUNC_title, ": PCA Scores", "\n", nrow(plot_Val), " samples; ", nrow(plotly_loadings_data), " features")
    )
  
  if(FUNC_option_invert_x){
    pca_output$plot_scores <- pca_output$plot_scores %>%
      layout(xaxis = list(autorange = "reversed"))
  }
  
  if(FUNC_option_invert_y){
    pca_output$plot_scores <- pca_output$plot_scores %>%
      layout(yaxis = list(autorange = "reversed"))
  }
  
  # create loadings plot
  x_axis_settings_loading <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    title = paste("")
  )
  
  y_axis_settings_loading <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    title = paste("")
  )
  
  pca_output$plot_loadings <- plot_ly(type = "scatter", 
                                      mode = "markers", 
                                      data = plotly_loadings_data, 
                                      x = ~PC1, 
                                      y = ~PC2, 
                                      text = ~variable, 
                                      marker = list(size = 10, color = '#808080', opacity = 0.5,
                                                    line = list(color = '#000000', width = 1)
                                      )) %>% 
    layout(
      xaxis = x_axis_settings_loading,
      yaxis = y_axis_settings_loading,
      showlegend = TRUE, 
      margin = list(l = 65, r = 50, b=65, t=85),
      title = paste0(FUNC_title, ": PCA Loadings", "\n", nrow(plot_Val), " samples; ", nrow(plotly_loadings_data), " features")
    )
  
  pca_output$plot_combined <- subplot(pca_output$plot_scores, 
                                      pca_output$plot_loadings, 
                                      nrows = 1,
                                      margin = 0.05,
                                      titleX = TRUE,
                                      titleY = TRUE
  ) %>% layout(showlegend = TRUE, 
               margin = list(l = 65, r = 50, b=65, t=85),
               title = paste0(FUNC_title, "\n", nrow(plot_Val), " samples; ", nrow(plotly_loadings_data), " features")
  )
  
  pca_output
}