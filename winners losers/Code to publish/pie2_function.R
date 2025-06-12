pie2 = function (x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE, 
          init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
          col = NULL, border = NULL, lty = NULL, main = NULL, line_length = 1.05, text_center = 1.1, textcol = 1, doubleplot_text=FALSE, ...) 
{
  if (!is.numeric(x) || any(is.na(x) | x < 0)) 
    stop("'x' values must be positive.")
  if (is.null(labels)) 
    labels <- as.character(seq_along(x))
  else labels <- as.graphicsAnnot(labels)
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L]) 
    xlim <- (pin[1L]/pin[2L]) * xlim
  else ylim <- (pin[2L]/pin[1L]) * ylim
  dev.hold()
  on.exit(dev.flush())
  plot.window(xlim, ylim, "", asp = 1)
  if (is.null(col)) 
    col <- if (is.null(density)) 
      c("white", "lightblue", "mistyrose", 
        "lightcyan", "lavender", "cornsilk")
  else par("fg")
  if (!is.null(col)) 
    col <- rep_len(col, nx)
  if (!is.null(border)) 
    border <- rep_len(border, nx)
  if (!is.null(lty)) 
    lty <- rep_len(lty, nx)
  angle <- rep(angle, nx)
  if (!is.null(density)) 
    density <- rep_len(density, nx)
  twopi <- if (clockwise) 
    -2 * pi
  else 2 * pi
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
    polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i])
    P <- t2xy(mean(x[i + 0:1]))
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)) {
      #lines(c(1, line_length) * P$x, c(1, line_length) * P$y)
      if(doubleplot_text) {
        text(text_center * P$x, text_center * P$y, labels[i], xpd = TRUE, 
             adj = ifelse(P$x < 0, 1, 0), cex = 1.1, col = 1, ...)
      }
      text(text_center * P$x, text_center * P$y, labels[i], xpd = TRUE, 
           #adj = ifelse(P$x < 0, 1, 0),
           col = textcol, ...)
    }
  }
  title(main = main, ...)
  invisible(NULL)
}

pie(c(1,1,1), lab = letters[1:3])
pie2(c(1,1,1), lab = letters[1:3], line_length = 1, text_center = 0.5, textcol = "white")
pie2(c(1,1,1), lab = "+", line_length = 1, text_center = 0.5, textcol = "white", doubleplot_text = TRUE)
