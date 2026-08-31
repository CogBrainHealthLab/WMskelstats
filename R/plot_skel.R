#' Interactive 3D and 2D Skeleton Data Visualizer
#'
#' Launches an interactive Shiny application to explore 3D spatial skeleton data,
#' overlay values, and template backgrounds across 3D interactive plots and 2D cross-sectional slices.
#' Includes dynamic thresholding to filter overlay values below a user-selected threshold.
#'
#' @param coords A matrix or data frame with at least 3 numeric columns representing 3D spatial coordinates (X, Y, and Z).
#'   The number of rows must match the length of \code{data}.
#' @param data A numeric vector of overlay values mapped to each 3D coordinate point. May contain \code{NA} values.
#' @param template An optional numeric vector of background or template values matching the length of \code{data}.
#'   Used to display background structure where \code{data} contains \code{NA} values. Defaults to \code{rep(NA, length(data))}.
#'
#' @details
#' The interactive Shiny application includes:
#' \itemize{
#'   \item \strong{Value Thresholding:} Dynamic slider that recodes overlay values below the threshold to \code{NA}.
#'   \item \strong{Auto-Scaling Density Sparklines:} Located directly beneath coordinate sliders, dynamically updating
#'     and scaling vertical height to fit remaining valid data distributions even under high sparsity.
#'   \item \strong{3D Interactive Plot:} Rendered via Plotly, allowing full 3D rotation, zooming, and point selection.
#'     Clicking a 3D point dynamically updates the 2D slice coordinate sliders.
#'   \item \strong{2D Cross-Sectional Slice Views:} Displays orthogonal Sagittal (Y-Z), Coronal (X-Z), and Axial (X-Y) views 
#'     for the currently selected slice plane.
#'   \item \strong{Export Capabilities:} Includes option to save publication-ready high-resolution PNG renders of the 2D slices.
#' }
#'
#' @return A \code{\link[shiny]{shinyApp}} object representing the interactive application.
#'
#' @importFrom shiny fluidPage titlePanel sidebarLayout sidebarPanel sliderInput plotOutput hr h4
#'   selectInput checkboxInput conditionalPanel br downloadButton mainPanel fluidRow column uiOutput
#'   observeEvent req updateSliderInput renderUI renderPlot reactive downloadHandler shinyApp
#' @importFrom plotly plotlyOutput renderPlotly plot_ly add_trace event_register layout event_data
#' @importFrom grDevices colorRampPalette adjustcolor png dev.off rainbow
#' @importFrom graphics par plot polygon abline grid points layout image axis box mtext
#' @importFrom stats density median bw.nrd0
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' n=5000
#' coords=cbind(X = runif(n, 1, 100), Y = runif(n, 1, 100), Z = runif(n, 1, 100))
#' data_vec=runif(n, 10, 50)
#' data_vec[c(1:1500, 3000:4000)]=NA
#' tmpl_vec=runif(n, 5, 25)
#' 
#' plot_interactive_skeleton(coords, data_vec, tmpl_vec)
#' }
##################################################################################################################
##################################################################################################################

plot_interactive_skeleton=function(coords, data, template = rep(NA, length(data))) {
  stopifnot(is.matrix(coords) || is.data.frame(coords))
  stopifnot(ncol(coords) >= 3)
  stopifnot(length(data) == nrow(coords))
  
  if (is.null(template)) {template=rep(NA, length(data))}
  stopifnot(length(template) == length(data))
  
  df=data.frame(X = coords[, 1],Y = coords[, 2],Z = coords[, 3],val = data,tmpl = template)
  
  # Compute unified range across template and overlay data
  all_vals=c(df$val, df$tmpl)
  all_vals=all_vals[!is.na(all_vals)]
  
  if (length(all_vals) == 0) {val_range=c(0, 1)} 
  else {
    val_range=range(all_vals)
    if (val_range[1] == val_range[2]) {val_range=c(val_range[1] - 0.5, val_range[2] + 0.5)}
  }
  
  # Determine bounds for data threshold slider
  val_vector=df$val[!is.na(df$val)]
  thresh_min=if (length(val_vector) > 0) min(val_vector) else 0
  thresh_max=if (length(val_vector) > 0) max(val_vector) else 1
  thresh_step=if (thresh_max > thresh_min) (thresh_max - thresh_min) / 100 else 0.1
  
  x_range=range(df$X, na.rm = TRUE)
  y_range=range(df$Y, na.rm = TRUE)
  z_range=range(df$Z, na.rm = TRUE)
  
  ui=fluidPage(
    titlePanel("Interactive Skeleton Visualizer"),
    
    sidebarLayout(
      sidebarPanel(width = 3,h4("Slice Coordinates"),
        
        sliderInput("x_slider", "X Coordinate (Sagittal):", min = floor(x_range[1]), max = ceiling(x_range[2]), value = round(median(df$X)), step = 1),
        plotOutput("x_density", height = "25px"),
        sliderInput("y_slider", "Y Coordinate (Coronal):", min = floor(y_range[1]), max = ceiling(y_range[2]), value = round(median(df$Y)), step = 1),
        plotOutput("y_density", height = "25px"),
        sliderInput("z_slider", "Z Coordinate (Axial):",min = floor(z_range[1]), max = ceiling(z_range[2]), value = round(median(df$Z)), step = 1),
        plotOutput("z_density", height = "25px"),
        
        hr(),
        h4("Data Controls"),
        sliderInput("thresh_val", "Value Threshold (Min):",min = round(thresh_min, 2), max = round(thresh_max, 2),value = round(thresh_min, 2), step = round(thresh_step, 2)),
        
        hr(),
        h4("Appearance Controls"),
        checkboxInput("black_bg", "Black Background", value = FALSE),
        selectInput("colorscale", "Color Palette:",choices = c("viridis", "plasma", "inferno", "magma", "cividis", "jet", "rainbow", "hot"),selected = "viridis"),
        sliderInput("pt_size", "Point Size:", min = 1, max = 10, value = 3, step = 0.5),
        sliderInput("na_opacity", "Template Opacity:", min = 0, max = 1, value = 0.5, step = 0.1),
        
        hr(),
        checkboxInput("show_2d", "Show 2D Slice Views", value = TRUE),
        conditionalPanel(condition = "input.show_2d == true",
        downloadButton("download_2d_png", "Save 2D Slices as PNG", class = "btn-primary", style = "width: 100%;")
        )
      ),
      
      mainPanel(width = 9,
        fluidRow(column(12, uiOutput("plot_3d_container"))),
        conditionalPanel(condition = "input.show_2d == true",
          br(),
          fluidRow(
            column(4, plotOutput("plot_sagittal", height = "450px")),
            column(4, plotOutput("plot_coronal", height = "450px")),
            column(4, plotOutput("plot_axial", height = "450px"))
          )
        )
      )
    )
  )
  
  server=function(input, output, session) {
    
    # Dynamic subsetting based on value thresholding
    df_filtered=reactive({
      req(input$thresh_val)
      df_mod=df
      if (!is.na(input$thresh_val)) {df_mod$val[!is.na(df_mod$val) & df_mod$val < input$thresh_val]=NA}
      df_mod
    })
    
    df_valid=reactive({ d=df_filtered(); d[!is.na(d$val), ] })
    df_na   =reactive({ d=df_filtered(); d[is.na(d$val) & !is.na(d$tmpl), ] })
    df_grey =reactive({ d=df_filtered(); d[is.na(d$val) & is.na(d$tmpl), ] })
    
    # Reactive local densities updated when threshold changes
    dens_x=reactive({ compute_slice_valid_ratio(df$X, df_valid()$X, x_range) })
    dens_y=reactive({ compute_slice_valid_ratio(df$Y, df_valid()$Y, y_range) })
    dens_z=reactive({ compute_slice_valid_ratio(df$Z, df_valid()$Z, z_range) })
    
    # 3D Click Observer
    observeEvent(suppressWarnings(event_data("plotly_click", source = "skeleton_3d")), {
      click_data=suppressWarnings(event_data("plotly_click", source = "skeleton_3d"))
      req(click_data)
      
      updateSliderInput(session, "x_slider", value = round(click_data$x))
      updateSliderInput(session, "y_slider", value = round(click_data$y))
      updateSliderInput(session, "z_slider", value = round(click_data$z))
    })
    
    # Dynamic 3D Plot Container Height
    output$plot_3d_container=renderUI({
      h=if (isTRUE(input$show_2d)) "500px" else "800px"
      plotlyOutput("plot_3d", height = h)
    })
    
    # 1D Density Plots
    output$x_density=renderPlot({ render_density_slice(dens_x(), input$x_slider, x_range) }, bg = "transparent")
    output$y_density=renderPlot({ render_density_slice(dens_y(), input$y_slider, y_range) }, bg = "transparent")
    output$z_density=renderPlot({ render_density_slice(dens_z(), input$z_slider, z_range) }, bg = "transparent")
    
    # Reactive slice extractions
    sagittal_valid=reactive({ req(input$x_slider); v=df_valid(); v[abs(v$X - input$x_slider) <= 0.5, ] })
    sagittal_na   =reactive({ req(input$x_slider); n=df_na();    n[abs(n$X - input$x_slider) <= 0.5, ] })
    sagittal_grey =reactive({ req(input$x_slider); g=df_grey();  g[abs(g$X - input$x_slider) <= 0.5, ] })
    
    coronal_valid =reactive({ req(input$y_slider); v=df_valid(); v[abs(v$Y - input$y_slider) <= 0.5, ] })
    coronal_na    =reactive({ req(input$y_slider); n=df_na();    n[abs(n$Y - input$y_slider) <= 0.5, ] })
    coronal_grey  =reactive({ req(input$y_slider); g=df_grey();  g[abs(g$Y - input$y_slider) <= 0.5, ] })
    
    axial_valid   =reactive({ req(input$z_slider); v=df_valid(); v[abs(v$Z - input$z_slider) <= 0.5, ] })
    axial_na      =reactive({ req(input$z_slider); n=df_na();    n[abs(n$Z - input$z_slider) <= 0.5, ] })
    axial_grey    =reactive({ req(input$z_slider); g=df_grey();  g[abs(g$Z - input$z_slider) <= 0.5, ] })
    
    # Modular 2D Slice Drawing Delegates
    draw_sagittal=function() {
      render_slice_2d(sagittal_valid(), sagittal_na(), sagittal_grey(), "Y", "Z",
                      y_range, z_range, "Y (Posterior - Anterior)", "Z (Inferior - Superior)",
                      paste("Sagittal View (X =", input$x_slider, ")"),
                      input$colorscale, val_range, input$pt_size, input$na_opacity,black_bg = input$black_bg)
    }
    
    draw_coronal=function() {
      render_slice_2d(coronal_valid(), coronal_na(), coronal_grey(), "X", "Z",
                      x_range, z_range, "X (Left - Right)", "Z (Inferior - Superior)",
                      paste("Coronal View (Y =", input$y_slider, ")"),
                      input$colorscale, val_range, input$pt_size, input$na_opacity,black_bg = input$black_bg)
    }
    
    draw_axial=function() {
      render_slice_2d(axial_valid(), axial_na(), axial_grey(), "X", "Y",
                      x_range, y_range, "X (Left - Right)", "Y (Posterior - Anterior)",
                      paste("Axial View (Z =", input$z_slider, ")"),
                      input$colorscale, val_range, input$pt_size, input$na_opacity,black_bg = input$black_bg)
    }
    
    # 3D Plotly View
    output$plot_3d=renderPlotly({
      input$show_2d
      
      d_valid=df_valid()
      d_na   =df_na()
      d_grey =df_grey()
      
      bg_3d      =if (isTRUE(input$black_bg)) "#000000" else "#FFFFFF"
      fg_3d      =if (isTRUE(input$black_bg)) "#FFFFFF" else "#000000"
      grid_3d    =if (isTRUE(input$black_bg)) "#333333" else "#E5E5E5"
      zeroline_3d=if (isTRUE(input$black_bg)) "#555555" else "#CCCCCC"
      grey_3d_col=if (isTRUE(input$black_bg)) "#555555" else "lightgrey"
      
      p=plot_ly(source = "skeleton_3d")
      
      has_valid=nrow(d_valid) > 0
      has_na   =nrow(d_na) > 0
      has_grey =nrow(d_grey) > 0
      
      show_cb_valid=has_valid
      show_cb_na   =has_na && !has_valid
      
      cb_style=list(title = list(text = "Value", font = list(color = fg_3d)),tickfont = list(color = fg_3d))
      
      if (has_grey) 
        {
        p=p %>% add_trace(x = d_grey$X, y = d_grey$Y, z = d_grey$Z,type = "scatter3d", mode = "markers",showscale = FALSE,showlegend = FALSE,
                          marker = list(
                            size = max(1, input$pt_size * 0.6),
                            color = grey_3d_col,
                            opacity = input$na_opacity)
                          )
        }
      
      if (has_na) 
        {
        p=p %>% add_trace(x = d_na$X, y = d_na$Y, z = d_na$Z,type = "scatter3d", mode = "markers",showscale = show_cb_na,showlegend = FALSE,
                          marker = list(
                            size = max(1, input$pt_size * 0.6), 
                            color = d_na$tmpl,
                            colorscale = make_plotly_colorscale(input$colorscale),
                            cmin = val_range[1],
                            cmax = val_range[2],
                            opacity = input$na_opacity,
                            colorbar = if (show_cb_na) cb_style else NULL)
                          )
        }
      
      if (has_valid) 
        {
        p=p %>% add_trace(x = d_valid$X, y = d_valid$Y, z = d_valid$Z,type = "scatter3d", mode = "markers",showscale = show_cb_valid,name = "Data Overlay",
                          marker = list(
                            size = input$pt_size,
                            opacity = 0.85,
                            color = d_valid$val,
                            colorscale = make_plotly_colorscale(input$colorscale),
                            cmin = val_range[1],
                            cmax = val_range[2],
                            colorbar = if (show_cb_valid) cb_style else NULL)
                        )
        }
      
      make_axis=function(name) {
        list(
          title = list(text = name, font = list(color = fg_3d)),
          tickfont = list(color = fg_3d),
          gridcolor = grid_3d,
          zerolinecolor = zeroline_3d,
          backgroundcolor = bg_3d,
          showbackground = TRUE
        )
      }
      
      axx = list(ticketmode = 'array',ticktext =  c("Left","Right"),tickvals =x_range)
      axy = list(ticketmode = 'array',ticktext = c("Posterior","Anterior"),tickvals = y_range)
      axz = list(ticketmode = 'array',ticktext = c("Inferior","Superior"),tickvals = z_range)
      
      p %>% 
        event_register("plotly_click") %>%
        layout(title = list(text = "3D Skeleton Overlay", font = list(color = fg_3d)),paper_bgcolor = bg_3d,plot_bgcolor = bg_3d,autosize = TRUE,
               scene = list(aspectmode = "data",xaxis = axx,yaxis = axy,zaxis = axz))
    })
    
    # 2D Render outputs
    output$plot_sagittal=renderPlot({ req(input$show_2d); draw_sagittal() })
    output$plot_coronal =renderPlot({ req(input$show_2d); draw_coronal() })
    output$plot_axial   =renderPlot({ req(input$show_2d); draw_axial() })
    
    # PNG Download handler
    output$download_2d_png=downloadHandler(filename = function() {
        paste0("skeleton_2d_slices_X", input$x_slider, "_Y", input$y_slider, "_Z", input$z_slider, ".png")
      },
      content = function(file) {
        bg_col=if (isTRUE(input$black_bg)) "black" else "white"
        fg_col=if (isTRUE(input$black_bg)) "white" else "black"
        
        png(file, width = 2400, height = 950, res = 200)
        
        layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE), heights = c(4.2, 1))
        
        draw_sagittal()
        draw_coronal()
        draw_axial()
        
        cols=get_palette_colors(input$colorscale, 256)
        x_vals=seq(val_range[1], val_range[2], length.out = 256)
        
        par(mar = c(3.5, 12, 1.5, 12))
        image(x = x_vals,y = c(0, 1),z = matrix(rep(x_vals, 2), nrow = 256, ncol = 2),col = cols,axes = FALSE,xlab = "", ylab = "")
        axis(1)
        box()
        mtext("Value", side = 3, line = 0.3, cex = 0.9, font = 2)
        
        dev.off()
      }
    )
  }
  
  shinyApp(ui = ui, server = server)
}
##################################################################################################################
##################################################################################################################
## helper functions
# Global palette generator
get_palette_colors=function(palette_name, n = 256) {
  pal_func=switch(
    palette_name,
    "viridis" = grDevices::colorRampPalette(c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725")),
    "plasma"  = grDevices::colorRampPalette(c("#000004", "#6A00A8", "#B12A90", "#E16462", "#FCA636", "#F0F921")),
    "inferno" = grDevices::colorRampPalette(c("#000004", "#420A68", "#932667", "#DD513A", "#FCA50A", "#FCFFA4")),
    "magma"   = grDevices::colorRampPalette(c("#000004", "#3B0F70", "#8C2981", "#DE4968", "#FE9F6D", "#FCFDBF")),
    "cividis" = grDevices::colorRampPalette(c("#002051", "#2C456B", "#576B71", "#89926B", "#C3BC67", "#FBEA55")),
    "jet"     = grDevices::colorRampPalette(c("blue", "cyan", "green", "yellow", "red")),
    "rainbow" = grDevices::colorRampPalette(rainbow(7)),
    "hot"     = grDevices::colorRampPalette(c("black", "red", "yellow", "white")),
    grDevices::colorRampPalette(c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"))
  )
  pal_func(n)
}

# Convert palette vector to Plotly colorscale format
make_plotly_colorscale=function(palette_name) {
  cols=get_palette_colors(palette_name, n = 256)
  n=length(cols)
  lapply(seq_len(n), function(i) {
    list((i - 1) / (n - 1), cols[i])
  })
}

# Base R color indexing helper for 2D slice plots
get_base_colors=function(vals, palette_name, limits) {
  if (length(vals) == 0) return(character(0))
  
  range_span=limits[2] - limits[1]
  norm_vals=if (range_span == 0) rep(0.5, length(vals)) else (vals - limits[1]) / range_span
  norm_vals=pmin(pmax(norm_vals, 0), 1)
  
  cols=get_palette_colors(palette_name, n = 256)
  col_idx=round(norm_vals * 255) + 1
  cols[col_idx]
}

# Compute localized fraction of non-NA points along an axis
compute_slice_valid_ratio=function(x_all, x_valid, x_rng) {
  if (length(x_all) == 0) return(NULL)
  
  if (x_rng[1] == x_rng[2]) {x_rng=c(x_rng[1] - 0.5, x_rng[2] + 0.5)}
  
  from_val=x_rng[1]
  to_val  =x_rng[2]
  
  if (length(x_valid) == length(x_all)) {
    x_grid=seq(from_val, to_val, length.out = 512)
    res=list(x = x_grid, y = rep(1, 512))
    class(res)="density"
    return(res)
  }
  
  if (length(x_valid) == 0) {
    x_grid=seq(from_val, to_val, length.out = 512)
    res=list(x = x_grid, y = rep(0, 512))
    class(res)="density"
    return(res)
  }
  
  dens_all=density(x_all, from = from_val, to = to_val, n = 512, cut = 0)
  bw_val=min(dens_all$bw, bw.nrd0(x_valid))
  dens_valid=density(x_valid, bw = bw_val, from = from_val, to = to_val, n = 512, cut = 0)
  
  cnt_all  =dens_all$y * length(x_all)
  cnt_valid=dens_valid$y * length(x_valid)
  
  ratio=ifelse(cnt_all > 1e-6, cnt_valid / cnt_all, 0)
  dens_all$y=pmin(pmax(ratio, 0), 1)
  
  return(dens_all)
}

# Helper to render 1D density sparklines with auto-scaled vertical height
render_density_slice=function(dens_obj, slider_val, coord_range) {
  req(slider_val, dens_obj)
  par(mar = c(0, 0, 0, 0), bg = NA)
  
  max_y=max(dens_obj$y, na.rm = TRUE)
  max_y=if (is.finite(max_y) && max_y > 0) max_y else 1
  
  plot(dens_obj$x, dens_obj$y, type = "n", xlim = coord_range, ylim = c(0, max_y),xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
  
  px=c(dens_obj$x[1], dens_obj$x, dens_obj$x[length(dens_obj$x)])
  py=c(0, dens_obj$y, 0)
  polygon(px, py, col = adjustcolor("#4682B4", alpha.f = 0.3), border = "#4682B4")
  
  abline(v = slider_val, col = "#D9534F", lty = 2, lwd = 2)
}

# Helper to render 2D slice scatter plots with optional dark background
render_slice_2d=function(valid_pts, tmpl_pts, grey_pts, x_var, y_var,x_range, y_range, xlab, ylab, main_title,colorscale, val_range, pt_size, na_opacity,black_bg = FALSE) 
  {
  cex_val=pt_size * 0.35
  
  bg_col     =if (black_bg) "black" else "white"
  fg_col     =if (black_bg) "white" else "black"
  grid_col   =if (black_bg) "#333333" else "#E5E5E5"
  grey_pt_col=if (black_bg) "#555555" else "lightgrey"
  
  par(mar = c(4, 4, 2.5, 1), bg = bg_col, col.axis = fg_col, col.lab = fg_col, col.main = fg_col, fg = fg_col)
  plot(NA, xlim = x_range, ylim = y_range, asp = 1, xaxs = "i", yaxs = "i",xlab = xlab, ylab = ylab, main = main_title, col.axis = fg_col, col.lab = fg_col)
  grid(col = grid_col, lty = "solid")
  
  if (nrow(grey_pts) > 0) {
    points(grey_pts[[x_var]], grey_pts[[y_var]],col = adjustcolor(grey_pt_col, alpha.f = na_opacity),pch = 16, cex = cex_val * 0.8)
  }
  if (nrow(tmpl_pts) > 0) {
    tmpl_cols=get_base_colors(tmpl_pts$tmpl, colorscale, val_range)
    points(tmpl_pts[[x_var]], tmpl_pts[[y_var]],col = adjustcolor(tmpl_cols, alpha.f = na_opacity),pch = 16, cex = cex_val * 0.8)
  }
  if (nrow(valid_pts) > 0) {
    cols=get_base_colors(valid_pts$val, colorscale, val_range)
    points(valid_pts[[x_var]], valid_pts[[y_var]],col = adjustcolor(cols, alpha.f = 0.9),pch = 16, cex = cex_val)
  }
}

# Helper to render 1D density sparklines with auto-scaled vertical height
render_density_slice=function(dens_obj, slider_val, coord_range) {
  req(slider_val, dens_obj)
  par(mar = c(0, 0, 0, 0), bg = NA)
  
  # Dynamically calculate maximum density peak for auto-scaling vertical plot bounds
  max_y=max(dens_obj$y, na.rm = TRUE)
  max_y=if (is.finite(max_y) && max_y > 0) max_y else 1
  
  plot(dens_obj$x, dens_obj$y, type = "n", xlim = coord_range, ylim = c(0, max_y),xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
  
  px=c(dens_obj$x[1], dens_obj$x, dens_obj$x[length(dens_obj$x)])
  py=c(0, dens_obj$y, 0)
  polygon(px, py, col = adjustcolor("#4682B4", alpha.f = 0.3), border = "#4682B4")
  abline(v = slider_val, col = "#D9534F", lty = 2, lwd = 2)
}
