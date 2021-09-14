library(parallel)

intersect_bed <- function(a, b){
  library(GenomicRanges)
  my_hit <- findOverlaps(a, b)
  class(my_hit)

  da = as.data.frame(a[queryHits(my_hit)])
  db = as.data.frame(b[subjectHits(my_hit)])

  #colnames(db) = paste0('b_', colnames(db))

  inter  <- cbind(da, db)
  return(inter)
}

get_divergent_color_set <- function()
{
  colors37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d",
               "#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22",
               "#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d",
               "#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")
  return(colors37)
}

get_marker_genes <- function(context)
{
  marker_genes_list = list()
  marker_genes_list[['leuk']] = c('MME', 'CD19', 'MS4A1', 'CD22', 'FCER2', 'CD40', 'IGHM', 'PAX5', 'EBF1'  #B Cell dev
                                  ,'HOXA7', 'HOXA9','MEIS1', 'RUNX2', 'MEF2C', 'FLT3' # MLLr-program
                                  , 'CD3D', 'CD3G', 'CD27', 'CD28'    # T cell
                                  ,'CD8A', 'CD8B',  'CD4', 'CCR7' , 'CCL5', 'TBX21',   'CXCR6'  #T cell
                                  ,'FCGR3A', 'NCAM1', 'NCR3' , 'NCR1', 'GNLY',   'NKG7' , 'GZMA', 'GZMB' #NK Cell
                                  ,'CD68', 'CSF1R' #Monocytes
                                  , 'KIT' , 'CD34' )  #Progenitors


  marker_genes_list[['hgg']] = c('SPARC', 'SPARCL1', 'GFAP', 'APOE' #Astro like
                                 ,'PDGFRA', 'ETV1','MBP' # OC-like
                                 , 'TOP2A', 'CDK1' #Cell cycle
                                 , 'MSR1', 'CD163' ,'TLR1', 'TLR2', 'CD80'     # Glioma-infiltration microglia
  )

  return(as.character(marker_genes_list[[context]]))
}

#pie with labels inside and without ticks
nice.pie <- function (x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE,
                   init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45,
                   col = NULL, border = NULL, lty = NULL, main = NULL, text_col ='white', ...)
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
      c("white", "lightblue", "mistyrose", "lightcyan",
        "lavender", "cornsilk")
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
      #lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)
      text(0.7 * P$x, 0.7 * P$y, labels[i], xpd = TRUE, col = text_col, font = 1,
           adj = ifelse(P$x < 0, 1, 0), ...)
    }
  }
  title(main = main, ...)
  invisible(NULL)
}



# An mc-version of the sapply function.
mcsapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
  FUN <- match.fun(FUN)
  answer <- parallel::mclapply(X = X, FUN = FUN, ...)
  if (USE.NAMES && is.character(X) && is.null(names(answer)))
    names(answer) <- X
  if (!isFALSE(simplify) && length(answer))
    simplify2array(answer, higher = (simplify == "array"))
  else answer
}


