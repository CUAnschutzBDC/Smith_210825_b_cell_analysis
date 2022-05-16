groupContinuousPlots <- function(group, plot_df, col_by, color = NULL,
                                 limits = NULL, axis_names = c("dim1", "dim2"),
                                 save_plot = NULL, show_legend = TRUE,
                                 size = 0.25, ggrastr = FALSE,
                                 raster_scale = 1, raster_res = 300) {
  plot_name_comb <- paste(group, collapse = "_")
  plot1 <- plot_df[plot_df$all == "all_samples", ]
  plot2 <- plot_df[plot_df$all != "all_samples", ]
  
  base_plot <- ggplot2::ggplot(data = plot2, ggplot2::aes(dim1, dim2)) +
    ggplot2::xlab(axis_names[1]) +
    ggplot2::ylab(axis_names[2])
  
  if(ggrastr){
    base_plot <- base_plot + ggrastr::geom_point_rast(data = plot1, 
                                                      ggplot2::aes(dim1, dim2), 
                                                      color = "#DCDCDC",
                                                      size = size,
                                                      show.legend = FALSE,
                                                      scale = raster_scale,
                                                      raster.dpi = raster_res)
    base_plot <- base_plot + ggrastr::geom_point_rast(data = plot2,
                                                      ggplot2::aes(dim1, dim2,
                                                                    color = colour_metric),
                                                      size = size,
                                                      show.legend = show_legend,
                                                      scale = raster_scale,
                                                      raster.dpi = raster_res)    
  } else {
    base_plot <- base_plot + ggplot2::geom_point(data = plot1, 
                                                 ggplot2::aes(dim1, dim2), 
                                                 color = "#DCDCDC",
                                                 size = size,
                                                 show.legend = FALSE)
    base_plot <- base_plot + ggplot2::geom_point(data = plot2,
                                                 ggplot2::aes(dim1, dim2,
                                                               color = colour_metric),
                                                 size = size,
                                                 show.legend = show_legend)
  }
  
  base_plot <- base_plot + #ggplot2::theme_classic() +
    ggplot2::ggtitle(paste0(plot_name_comb, " ", col_by)) +
    ggplot2::xlab(axis_names[1]) +
    ggplot2::ylab(axis_names[2])
  
  if (is.null(color)) {
    base_plot <- base_plot + ggplot2::scale_color_viridis_c(option = "magma") +
      labs(color = col_by)
  } else {
    low <- color[1]
    high <- color[2]
    if(is.null(limits)){
      base_plot <- base_plot + 
        ggplot2::scale_color_gradient(low = low, high = high, name = col_by)
    } else {
      base_plot <- base_plot + ggplot2::scale_color_gradient(low = low,
                                                             high = high, 
                                                             name = col_by,
                                                             limits = limits)
    }
  }
  
  return(base_plot)
}


groupDiscretePlots <- function(group, plot_df, col_by, axis_names = c("dim1", "dim2"),
                               color = NULL, show_legend = TRUE,
                               size = 0.25, ggrastr = FALSE,
                               raster_scale = 1, raster_res = 300) {
  plot1 <- plot_df[plot_df$all == "all_samples", ]
  plot2 <- plot_df[plot_df$all != "all_samples", ]
  
  base_plot <- ggplot2::ggplot(data = plot2, ggplot2::aes(dim1,
                                                          dim2)) +
    ggplot2::xlab(axis_names[1]) +
    ggplot2::ylab(axis_names[2])
  
  if(ggrastr){
    base_plot <- base_plot + ggrastr::geom_point_rast(data = plot1, 
                                                      ggplot2::aes(dim1, dim2), 
                                                      color = "#DCDCDC",
                                                      size = size,
                                                      show.legend = FALSE,
                                                      scale = raster_scale,
                                                      raster.dpi = raster_res)
    base_plot <- base_plot + ggrastr::geom_point_rast(data = plot2,
                                                      ggplot2::aes(dim1, dim2,
                                                                    color = colour_metric),
                                                      size = size,
                                                      show.legend = show_legend,
                                                      scale = raster_scale,
                                                      raster.dpi = raster_res)  
  } else {
    base_plot <- base_plot + ggplot2::geom_point(data = plot1, 
                                                 ggplot2::aes(dim1, dim2), 
                                                 color = "#DCDCDC",
                                                 size = size,
                                                 show.legend = FALSE)
    base_plot <- base_plot + ggplot2::geom_point(data = plot2,
                                                 ggplot2::aes_(dim1, dim2,
                                                               color = colour_metric),
                                                 size = size,
                                                 show.legend = show_legend)
  }
  
  base_plot <- base_plot + #ggplot2::theme_classic() + 
    ggplot2::ggtitle(paste(group, collapse = "_")) +
    ggplot2::xlab(axis_names[1]) +
    ggplot2::ylab(axis_names[2]) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3)))
  if (is.null(color)) {
    if(!is.null(levels(plot2$color_metric))){
      nColors <- length(levels(plot2$colour_metric))
    } else {
      nColors <- length(unique(plot2$colour_metric))
      
    }
    base_plot <- base_plot + ggplot2::scale_color_manual(
      values = grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(9, "Set1"))(nColors), name = col_by)
  } else {
    base_plot <- base_plot +
      ggplot2::scale_color_manual(values = color, name = col_by)
  }
  
  return(base_plot)
}
