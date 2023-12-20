##############
# new plot
##############

##code written by Adam Clark

make_boxplot = function(toplot_data = toplot, # data to plot
                        trt.labels_data = trt.labels, # treatment labels
                        trait.labels_data=trait.labels, # trait labels
                        p_alpha = 0.025, # alpha value for significance tests
                        groupbytrait = FALSE, # should individual panels include multiple traits?
                        legend_line_length = 1.55, # length of lines in legend
                        legend_cex = 0.9, # text size in legend
                        lower_margin = 5.5, # margin for plotting category names
                        sigadj = 0, # shift the asterisk up or down,
                        legend_textwidth = 0.2, # text spacing for legend
                        legend_yadj = 0.35, # margin spacing for legend
                        legend_xadj = 0, # horizontal adjustment for the legend 
                        flipaxes = TRUE, # plot with quantitative values on the x-axis
                        traitorder = NULL, # optional order in which to plot traits
                        group_colors = NULL,  # optional vector of colors for backgrounds of individual panels
                        xtit = "Effect Size, DCi diff. vs. Standardized ln(Trait)", # axis title
                        xlm = NULL,  # optional x-axis limits
                        collst = NULL, # optional vector of colors for the bars
                        autopar = TRUE, # should par settings be handled automatically?
                        axisside = 1+flipaxes, # which side should the axis labels be plotted on?
                        n_table = NULL, ncex = 0.6, yadj_text = 0) { # optional table of n values to be plotted; size for plotted text; vertical adjustment for text
  lwd_lines = 2.5
  lwd_diff_fact = 0.5
  
  # get significance based on estimate, se, and desired alpha
  tmp = pnorm(0, abs(toplot_data$Estimate), toplot_data$SE) <= p_alpha
  sig = tmp; sig[tmp] = "*"; sig[!tmp] = ""
  
  trt_type2_levels = as.character(sort(unique(toplot_data$trt_type2[!is.na(toplot_data$trt_type2)])))
  trt_type2_levels = trt_type2_levels[match(names(trt.labels_data), trt_type2_levels)]
  
  trait_levels = sort(unique(toplot_data$trait))
  trait_levels = trait_levels[match(names(trait.labels_data), trait_levels)]
  if(!is.null(traitorder) & !groupbytrait) {
    trait_levels = trait_levels[traitorder]
    trait.labels_data = trait.labels_data[traitorder]
  }
  
  if(groupbytrait) {
    v1 = trt_type2_levels
    v2 = trait_levels
  } else {
    v1 = trait_levels
    v2 = trt_type2_levels
  }
  if(is.null(collst)) {
    collst = RColorBrewer::brewer.pal(n = length(v2), name = "Dark2")
  }
  
  xps = seq(1, length(v1))
  dx = seq(0.35, -0.35, length=length(v2))
  whisker_length = min(diff(dx))/5
  
  yrng = c(min(toplot_data$Estimate-2*toplot_data$SE),
           max(toplot_data$Estimate+2*toplot_data$SE))
  
  if(autopar) {
    if(flipaxes) {
      par(mar=c(4,lower_margin,3.5,1))
    } else {
      par(mar=c(lower_margin,4,3,1))
    }
  }
  
  if(flipaxes) {
    if(!is.null(xlm)) {
      plot(yrng, c(min(xps)-0.5, max(xps)+0.5), xlab = "", ylab = "", type = "n", axes = FALSE, yaxs = "i", xlim = xlm)
    } else {
      plot(yrng, c(min(xps)-0.5, max(xps)+0.5), xlab = "", ylab = "", type = "n", axes = FALSE, yaxs = "i")
    }
    axis(1)
    if(groupbytrait) {
      axis(axisside, at = xps, labels = trt.labels_data, las = 2)
    } else {
      axis(axisside, at = xps, labels = trait.labels_data, las = 2)
    }
  } else {
    if(!is.null(xlm)) {
      plot(c(min(xps)-0.5, max(xps)+0.5), yrng, xlab = "", ylab = "", type = "n", axes = FALSE, xaxs = "i", xlim = xlm)
    } else {
      plot(c(min(xps)-0.5, max(xps)+0.5), yrng, xlab = "", ylab = "", type = "n", axes = FALSE, xaxs = "i")
    }
    axis(2)
    if(groupbytrait) {
      axis(axisside, at = xps, labels = trt.labels_data, las = 2)
    } else {
      axis(axisside, at = xps, labels = trait.labels_data, las = 2)
    }
  }
  box()
  
  if(!is.null(group_colors)) {
    for(i in 1:length(v1)) {
      if(flipaxes) {
        polygon(rep(c(yrng[1]-diff(range(yrng)), yrng[2]+diff(range(yrng))), each=2),
                xps[i]+c(-0.5, 0.5, 0.5, -0.5),
                col = group_colors[i], border = NA)
      }
    }
  }
  
  if(flipaxes) {
    abline(v=0, lty=2)
  } else {
    abline(h=0, lty=2)
  }
  
  for(i in 1:length(v1)) {
    for(j in 1:length(v2)) {
      if(groupbytrait) {
        ps = which(toplot_data$trait == trait_levels[j] & 
                     toplot_data$trt_type2 == trt_type2_levels[i])
      } else {
        ps = which(toplot_data$trait == trait_levels[i] & 
                     toplot_data$trt_type2 == trt_type2_levels[j])
      }
      
      xv = xps[i] + dx[j]
      yv = c(toplot_data$Estimate[ps]-toplot_data$SE[ps],
             toplot_data$Estimate[ps],
             toplot_data$Estimate[ps]+toplot_data$SE[ps],
             qnorm(p_alpha, toplot_data$Estimate[ps], toplot_data$SE[ps]),
             qnorm(1-p_alpha, toplot_data$Estimate[ps], toplot_data$SE[ps]))
      
      if(flipaxes) {
        segments(yv[1], xv, yv[3], xv,
                 col = collst[j], lwd = lwd_lines, lend = 2)
        segments(yv[4], xv, yv[5], xv,
                 col = collst[j], lwd = lwd_lines*lwd_diff_fact, lend = 2)
        segments(yv[c(4,5)], xv+whisker_length, yv[c(4,5)], xv-whisker_length,
                 col = collst[j], lwd = lwd_lines*lwd_diff_fact, lend = 2)
        points(yv[2], xv, col = collst[j], pch = 16, cex = 0.7)
      } else {
        segments(xv, yv[1], xv, yv[3],
                 col = collst[j], lwd = lwd_lines, lend = 2)
        segments(xv, yv[4], xv, yv[5],
                 col = collst[j], lwd = lwd_lines*lwd_diff_fact, lend = 2)
        segments(xv-whisker_length, yv[c(4,5)], xv+whisker_length, yv[c(4,5)],
                 col = collst[j], lwd = lwd_lines*lwd_diff_fact, lend = 2)
        points(xv, yv[2], col = collst[j], pch = 16, cex = 0.7)
      }
      ## add significance marker
      if(sum(ps)>0) {
        if(flipaxes) {
          if(all(yv[1:3]<0)) {
            text(yv[4], xv+sigadj, sig[ps], pos = 2,
                 cex = 0.8, offset = 0.1)
          } else {
            text(yv[5], xv+sigadj, sig[ps], pos = 4,
                 cex = 0.8, offset = 0.1)
          }
        } else {
          text(xv+sigadj, yv[3], sig[ps], pos = 3,
               cex = 0.8, offset = 0)
        }
      }
      
      #add n
      if(flipaxes) {
        if(!is.null(n_table)) {
          row_ps = which(row.names(n_table)==v2[j])
          col_ps = which(colnames(n_table)==v1[i])
          tmp_xrng = range(c(yrng, xlm))
          xps_text = tmp_xrng[2]-diff(tmp_xrng)*0.03
          yps_text = xv+yadj_text
          
          text(xps_text, yps_text, paste("(", n_table[row_ps, col_ps], ")", sep =""), cex=ncex)
        }
      }
      
    }
  }
  if(flipaxes) {
    abline(h = seq(1.5, length(v1)-0.5, by = 1), col = "darkgrey")
  } else {
    abline(v = seq(1.5, length(v1)-0.5, by = 1), col = "darkgrey")
  }
  
  if(flipaxes) {
    if(groupbytrait) {
      legend(mean(yrng)+legend_xadj, max(xps)+diff(range(xps))*legend_yadj, legend = trait.labels_data, lty = 1, lwd = lwd_lines*1.5,
             col = collst, ncol = ceiling(length(trait.labels_data)/2), xpd = NA, xjust = 0.5,
             seg.len = legend_line_length, cex = legend_cex, x.intersp = 0.5, text.width = legend_textwidth, bty = "n")
    } else {
      legend(mean(yrng)+legend_xadj, max(xps)+diff(range(xps))*legend_yadj, legend = trt.labels_data, lty = 1, lwd = lwd_lines*1.5,
             col = collst, ncol = ceiling(length(trt.labels_data)/2), xpd = NA, xjust = 0.5,
             seg.len = legend_line_length, cex = legend_cex, x.intersp = 0.5, text.width = legend_textwidth, bty = "n")
    }
    mtext(xtit, 1, line = 2.4)
  } else {
    if(groupbytrait) {
      legend(mean(xps)+legend_xadj, yrng[2]+diff(range(yrng))*legend_yadj, legend = trait.labels_data, lty = 1, lwd = lwd_lines*1.5,
             col = collst, ncol = ceiling(length(v2)/2), xpd = NA, xjust = 0.5,
             seg.len = legend_line_length, cex = legend_cex, x.intersp = 0.5, text.width = legend_textwidth, bty = "n")
    } else {
      legend(mean(xps)+legend_xadj, yrng[2]+diff(range(yrng))*legend_yadj, legend = trt.labels_data, lty = 1, lwd = lwd_lines*1.5,
             col = collst, ncol = ceiling(length(v2)/2), xpd = NA, xjust = 0.5,
             seg.len = legend_line_length, cex = legend_cex, x.intersp = 0.5, text.width = legend_textwidth, bty = "n")
    }
    mtext(xtit, 2, line = 2.4)
  }
}


