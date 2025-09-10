# Resize selected elements of an existing ggplot object
# - repel_size: size (mm) for ggrepel text/labels
# - text_size:  size (mm) for annotate("text") / annotate("label")
# - axis_title_size: base size (pt) for axis titles
# - axis_text_size:  base size (pt) for tick labels
# - axis_title_size_x / _y and axis_text_size_x / _y override per-axis (optional)
# - point_size: fixed size (mm) for point-like geoms
# - line_width: fixed width (mm) for line-like geoms (handles both size & linewidth fields)
# - tick_length: unit for tick length, e.g. unit(3, "mm")  (optional)
resize_plot_elements <- function(
    p,
    repel_size = NULL,
    text_size  = NULL,
    axis_title_size = NULL,
    axis_text_size  = NULL,
    axis_title_size_x = NULL,
    axis_title_size_y = NULL,
    axis_text_size_x  = NULL,
    axis_text_size_y  = NULL,
    point_size = NULL,
    line_width = NULL,
    tick_length = NULL
) {
  # ----- 1) LAYERS: text & symbols -----
  for (i in seq_along(p$layers)) {
    lyr <- p$layers[[i]]
    g   <- class(lyr$geom)[1]
    
    set_size <- function(lyr, s) {
      lyr$aes_params$size  <- s
      lyr$geom_params$size <- s
      if (!is.null(lyr$mapping) && "size" %in% names(lyr$mapping)) {
        m <- lyr$mapping; m$size <- NULL; lyr$mapping <- m
      }
      lyr
    }
    set_linewidth <- function(lyr, w) {
      # ggplot2 >= 3.4 uses linewidth; many geoms still accept size too.
      lyr$aes_params$linewidth  <- w
      lyr$geom_params$linewidth <- w
      lyr$aes_params$size  <- w
      lyr$geom_params$size <- w
      if (!is.null(lyr$mapping) && "linewidth" %in% names(lyr$mapping)) {
        m <- lyr$mapping; m$linewidth <- NULL; lyr$mapping <- m
      }
      if (!is.null(lyr$mapping) && "size" %in% names(lyr$mapping)) {
        m <- lyr$mapping; m$size <- NULL; lyr$mapping <- m
      }
      lyr
    }
    
    # ggrepel layers
    if (!is.null(repel_size) && grepl("Geom(Text|Label)Repel", g)) {
      lyr <- set_size(lyr, repel_size)
    }
    
    # plain text/label layers (incl. annotate("text"/"label"))
    if (!is.null(text_size) && g %in% c("GeomText", "GeomLabel")) {
      lyr <- set_size(lyr, text_size)
    }
    
    # point-like geoms
    if (!is.null(point_size) && g %in% c("GeomPoint", "GeomJitter", "GeomDotplot")) {
      lyr <- set_size(lyr, point_size)
    }
    
    # line-like geoms (common ones)
    if (!is.null(line_width) && g %in% c("GeomLine", "GeomPath", "GeomAbline",
                                         "GeomSegment", "GeomHline", "GeomVline",
                                         "GeomSmooth", "GeomErrorbar", "GeomErrorbarh",
                                         "GeomLinerange", "GeomCrossbar", "GeomCurve",
                                         "GeomSpoke", "GeomContour", "GeomDensity")) {
      lyr <- set_linewidth(lyr, line_width)
    }
    
    p$layers[[i]] <- lyr
  }
  
  # ----- 2) THEME: axes & ticks -----
  # Build theme updates selectively (don’t clobber user theme)
  thm <- theme()
  
  # Axis titles
  if (!is.null(axis_title_size)) {
    thm <- thm + theme(axis.title = element_text(size = axis_title_size))
  }
  if (!is.null(axis_title_size_x)) {
    thm <- thm + theme(axis.title.x = element_text(size = axis_title_size_x))
  }
  if (!is.null(axis_title_size_y)) {
    thm <- thm + theme(axis.title.y = element_text(size = axis_title_size_y))
  }
  
  # Tick labels
  if (!is.null(axis_text_size)) {
    thm <- thm + theme(axis.text = element_text(size = axis_text_size))
  }
  if (!is.null(axis_text_size_x)) {
    thm <- thm + theme(axis.text.x = element_text(size = axis_text_size_x))
  }
  if (!is.null(axis_text_size_y)) {
    thm <- thm + theme(axis.text.y = element_text(size = axis_text_size_y))
  }
  
  # Tick length (optional)
  if (!is.null(tick_length)) {
    thm <- thm + theme(axis.ticks.length = tick_length)
  }
  
  p + thm
}
