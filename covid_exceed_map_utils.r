# Script Title: covid_exceed_map_utils.r
# Author: Massimo Cavallaro
# Description: Defines a number of graphical function and utils to plot color
# maps.
# Dependencies: sf, viridisLite

map = sf::st_geometry(read_sf('Data/LTLA.shp'))

float2int<-function(x, n=10){
  ifelse(x<1, as.integer(x * n) + 1, n)
}

wilson.ci<-function(prop, n){
  if (n>0){
    return(suppressWarnings(c(prop.test(prop * n, n)$conf.int)))    
  }else{
    return(c(0, 1)) 
  }
}

Wilson.ci<-Vectorize(wilson.ci)

plot_exceed_map<-function(ws, baseline.matrix, observation.matrix, column, label, map,
                          region_names=LTLA_names, n_colors=10, n.matrix=NULL){  
  par(mfrow = c(1,3))
  color_names = c("#ff7f0e","#2ca02c","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
  
  # tabulate color palette:
  if (is.null(n.matrix)){
    Palette = colorRampPalette(c("#1f77b488","#d6272888"), alpha = TRUE)(n_colors)
    colors = Palette[float2int(ws[,column], n_colors)]
  }else{
    Palette = colorRampPalette(c("#1f77b4","#d62728"), alpha = F)(n_colors)
    main_colors = Palette[float2int(ws[,column], n_colors)]
    
    # Alphas = c("F2","D9","BF","A6","8C","73","59","40","26","0D")
    Alphas = c("F2","D9","BF","A6","8C","73","59","40","26","26")
    opacities = Alphas[float2int(diff(Wilson.ci(ws[,column], n.matrix[,column])))]
    
    colors = sapply(1:nrow(ws), function(x){paste0(main_colors[x], opacities[x])})
  }
  
  if(is.null(n.matrix)){
    #take the index of the highest warning scores (mean over last few days),
    ws_mean = rowMeans(ws[,(column-3):column])
    idx = order(ws_mean, decreasing = T)[1:3]
  }else{
    #take the index of the highest warning scores, but only if statistically significant
    ss = n.matrix[,column] > 10
    ws_mean = rowMeans(ws[,(column-3):column])
    ws_min = min(ws_mean)
    ws_mean = ifelse(ss, ws_mean, ws_min)
    idx = order(ws_mean, decreasing = T)
    nn = min(3, length(idx))
    if(nn>0){
      idx = idx[1:nn]      
    }else{
      idx = NULL
    }
  }
  day = as.integer(colnames(ws)[i])
  x = (day-10):day
  
  # plot warning scores
  plot(range(x), c(0,1), col='white', xlim=, ylim=c(0,1), ylab = expression(w%*%100), xlab = 'Day')
  title(main=sprintf('Top %d warning scores for %s (3-day average)', length(idx), label), cex.main = 0.8)
  if(!is.null(idx)){
    if (!is.null(n.matrix)){
      for(I in 1:length(idx)){
        ci = Wilson.ci(ws[idx[I], as.character(x)], n.matrix[idx[I], as.character(x)])
        polygon(c(x, rev(x)), c(ci[1,] ,rev(ci[2,])), col = paste0(color_names[I], '50'), border = NA)      
      }
    }
    for (I in 1:length(idx)){
      lines(x, ws[idx[I],as.character(x)], col=color_names[I])
      points(x, ws[idx[I],as.character(x)], col=color_names[I])
    }
    legend('bottomleft', region_names[idx,], lty=c(1,1), pch=c(1,1), col = color_names[1:length(idx)])
  }
  # plot baseline and counts
  # print(range(c(baseline.matrix[idx, as.character(x)], observation.matrix[idx, as.character(x)])))
  plot(range(x),
       range(c(baseline.matrix[idx, as.character(x)], observation.matrix[idx, as.character(x)])),
       col='white', ylab='No. of cases', xlab = 'Day')
  if(!is.null(idx)){
    for (I in 1:length(idx)){
      # print(idx[I])
      lines(x, baseline.matrix[idx[I], as.character(x)], col=color_names[I], lty=1)
      lower = qpois(0.05, baseline.matrix[idx[I], as.character(x)])
      upper = qpois(0.95, baseline.matrix[idx[I], as.character(x)])
      polygon(c(x, rev(x)), c(lower ,rev(upper)), col = paste0(color_names[I], '50'), border = NA)
      lines(x, observation.matrix[idx[I],as.character(x)], col=color_names[I], lty=2)
    }      
  }
  legend('topleft', c('Baseline', 'Observations'), lty=c(1,2), cex=0.4)
  
  # plot map
  title = sprintf("Day: %d", day)
  plot(map, col = colors, border = 'white', lwd = 0.2, main = title)
  if(!is.null(idx)){
    plot(map[idx], border = color_names[1:length(idx)], col=NA, add=T, lwd=0.5)    
  }
  par(fig=c(0.70, 0.76, 0.4, 0.6), new=T)
  par(mar=c(0.2, 2.2, 0.2, 0.2))
  plot(c(0, 1), c(0,1), ylab = 'w', yaxt='n',
       xaxt='n', col='white', frame.plot=F)
  axis(2, at = c(0,0.5,1), lwd=0.5)
  axis(1, at = c(0,0.5,1), lwd=0.5, cex=0.3)
  mtext(side=2, 'w', line = 2)
  if(!is.null(n.matrix)){
    mtext(side=1, 'CI width', line = 2, cex = 0.6)
    Cb = sapply(1:10, function(i){sapply(1:10, function(x){paste0(Palette[x], Alphas[i])}) })
  }else{
    Cb = matrix(rep(Palette, 10), nrow = length(Palette), ncol = 10)
  }
  rasterImage(Cb[seq(length(Palette), 1, -1),], 0, 0, 1, 1)
  rect(0, 0, 1, 1, lwd=0.5)
}



plot_all<-function(Positives, Negatives, vline){
  plot(colSums(Positives), type = 'l', xlab = '', ylab = '# cases', col='blue', xaxt='n')
  
  xtick = c(1, 367, 732)
  axis(side=1, at=xtick, labels = FALSE)
  text(x=c(183, 518, 775),  par("usr")[3],
       labels = c('2020', '2021', '2022'),
       pos = 1, xpd = TRUE, cex = 1.1)
  
  lines(colSums(Negatives),  col='orange')
  mtext(side=3, "Whole-England observations", cex=0.7, adj=0)
  abline(v=vline, col='grey')
  legend('topleft', c('Positives','Negatives'),lty=1,col=c('blue','orange'),  box.col = 'grey')
  VOC = c(353, 497, 694)
  text(x=c(353, 497, 684)+50, 37000,
       labels = c("Alpha","Delta","Omicron"),
       pos = 1, xpd = TRUE, cex = 1.1)

  abline(v=VOC, lwd=1, lty='dashed', col=tab.gray)
  
}


int2daydate<-function(x){
  day_names = as.POSIXct('2022-2-21', '%Y-%m-%d', tz='GMT') - as.difftime(783 - x,units='days')
  return(as.character(day_names))
}


int2daydate.lubridate<-function(x){
  day_names = as.POSIXct('2022-2-21', '%Y-%m-%d', tz='GMT') - as.difftime(783 - x,units='days')
  Day = as.character(lubridate::day(day_names))
  Month = as.character(lubridate::month(day_names,label = T))
  Year = paste0("'",substring(as.character(lubridate::year(day_names)), 3))
  ret = format( paste(Day, Month, Year, sep=' '), width=10, justify='r')
  return(ret)
}


plot_exceed_map2<-function(ws.p, ws.f,
                           baseline.matrix.p, observation.matrix.p,
                           baseline.matrix.f, observation.matrix.f, column, map,
                          region_names=LTLA_names, n_colors=10, n.matrix=NULL){  
  par(mfrow = c(2,3))
  color_names = c("#ff7f0e","#2ca02c","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
  labels = c('Positives', 'Negatives')
  
  ws.list = list(ws.p, ws.f)
  baseline.matrix.list = list(baseline.matrix.p, baseline.matrix.f)
  observation.matrix.list = list(observation.matrix.p, observation.matrix.f)
  for (iii in 1:2){
    ws = ws.list[[iii]]
    baseline.matrix = baseline.matrix.list[[iii]]
    observation.matrix = observation.matrix.list[[iii]]
    label=labels[iii]
  # tabulate color palette:
    
    if (is.null(n.matrix)){
      Palette = colorRampPalette(c("#1f77b488","#d6272888"), alpha = TRUE)(n_colors)
      colors = Palette[float2int(ws[,column], n_colors)]
    }else{
      Palette = colorRampPalette(c("#1f77b4","#d62728"), alpha = F)(n_colors)
      main_colors = Palette[float2int(ws[,column], n_colors)]
      
      Alphas = c("F2","D9","BF","A6","8C","73","59","40","26","0D")
      opacities = Alphas[float2int(diff(Wilson.ci(ws[,column], n.matrix[,column])))]
      
      colors = sapply(1:nrow(ws), function(x){paste0(main_colors[x], opacities[x])})
    }
    
    if(is.null(n.matrix)){
      #take the index of the highest warning scores (mean over last few days),
      ws_mean = rowMeans(ws[,(column-3):column])
      idx = order(ws_mean, decreasing = T)[1:3]
    }else{
      #take the index of the highest warning scores, but only if statistically significant
      ss = n.matrix[,column] > 10
      ws_mean = rowMeans(ws[,(column-3):column])
      ws_min = min(ws_mean)
      ws_mean = ifelse(ss, ws_mean, ws_min)
      idx = order(ws_mean, decreasing = T)
      nn = min(3, length(idx))
      if(nn>0){
        idx = idx[1:nn]      
      }else{
        idx = NULL
      }
    }
    day = as.integer(colnames(ws)[column])
    x = (day-10):day
    
    # plot warning scores
    plot(range(x), c(0,1), type='n',  xlim=, ylim=c(0,1), ylab = 'w', xlab = 'Day', cex.lab=1.2)
    title(main=sprintf('Top %d warning scores for %s (3-day average)', length(idx), label), cex.main = 0.8)
    if(!is.null(idx)){
      if (!is.null(n.matrix)){
        for(I in 1:length(idx)){
          a = ws[idx[I], as.character(x)]
          b = n.matrix[idx[I], as.character(x)]
          ci = Wilson.ci(a,b)
          polygon(c(x, rev(x)), c(ci[1,] ,rev(ci[2,])), col = paste0(color_names[I], '50'), border = NA)      
        }
      }
      for (I in 1:length(idx)){
        lines(x, ws[idx[I],as.character(x)], col=color_names[I])
        points(x, ws[idx[I],as.character(x)], col=color_names[I])
      }
      legend('bottomleft', region_names[idx,], lty=c(1,1), pch=c(1,1), col = color_names[1:length(idx)])
    }
    # plot baseline and counts
    # print(range(c(baseline.matrix[idx, as.character(x)], observation.matrix[idx, as.character(x)])))
    Range = range(c(baseline.matrix[idx, as.character(x)], observation.matrix[idx, as.character(x)]))
    if (Range[2] <= 1){
      Range[2] = 1
    }
    plot(range(x),
         Range,
         type='n', ylab='No. of cases', xlab = 'day')
    if(!is.null(idx)){
      for (I in 1:length(idx)){
        # print(idx[I])
        lines(x, baseline.matrix[idx[I], as.character(x)], col=color_names[I], lty=1)
        lower = qpois(0.05, baseline.matrix[idx[I], as.character(x)])
        upper = qpois(0.95, baseline.matrix[idx[I], as.character(x)])
        polygon(c(x, rev(x)), c(lower ,rev(upper)), col = paste0(color_names[I], '50'), border = NA)
        lines(x, observation.matrix[idx[I],as.character(x)], col=color_names[I], lty=2)
      }      
    }
    legend('topleft', c('Baseline', 'Observations'), lty=c(1,2) )
    
    # plot map
    title = sprintf("%s.", as.character(int2daydate.lubridate(day)))
    plot(map, col = colors, border = 'white', lwd = 0.1, main = title)
    if(!is.null(idx)){
      plot(map[idx], border = color_names[1:length(idx)], col=NA, add=T, lwd=0.5)    
    }
  }
  par(fig=c(0.70, 0.76, 0.8, 0.9), new=T)
  par(mar=c(0.2, 2.2, 0.2, 0.2))
  plot(c(0, 1), c(0,1), ylab = 'w', yaxt='n',
         xaxt='n', type='n', frame.plot=F)
  axis(2, at = c(0,0.5,1), lwd=0.5)
  axis(1, at = c(0,0.5,1), lwd=0.5, cex=0.3)
  mtext(side=2, expression(w%*%100), line = 2)
  if(!is.null(n.matrix)){
    mtext(side=1, 'CI width', line = 2, cex = 0.6)
    Cb = sapply(1:10, function(i){sapply(1:10, function(x){paste0(Palette[x], Alphas[i])}) })
  }else{
    Cb = matrix(rep(Palette, 10), nrow = length(Palette), ncol = 10)
  }
  rasterImage(Cb[seq(length(Palette), 1, -1),], 0, 0, 1, 1)
  rect(0, 0, 1, 1, lwd=0.5)

  par(fig=c(0.65, 0.82, 0.22, 0.4), new=T)
  par(mar=c(0.2, 2.2, 0.2, 0.2))
  plot_all(observation.matrix.p,  observation.matrix.f, as.integer(colnames(ws.p)[column]))
}


plot_exceed_map3<-function(ws.p, ws.f, averages.p, averages.f,
                           baseline.matrix.p, observation.matrix.p,
                           baseline.matrix.f, observation.matrix.f, column, map,
                           region_names=LTLA_names, n_colors=10,
                           n.matrix.p=NULL,
                           n.matrix.f=NULL,
                           method = 'rancover'){
  
  TMP = 999
  
  par(mfrow = c(2, 3), mar = c(10,8,4,4)/2)
  color_names = c("#ff7f0e","#2ca02c","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
  labels = c('S-gene positives', 'S-gene negatives')

  Y1 = c(0.75, 0.25) - 0.02
  Y2 = c(0.93, 0.43)
  x1 = 0.05
  x2 = 0.2
  
  ws.list = list(ws.p, ws.f)
  baseline.matrix.list = list(baseline.matrix.p, baseline.matrix.f)
  observation.matrix.list = list(observation.matrix.p, observation.matrix.f)
  averages.list = list(averages.p, averages.f)
  n.matrix.list = list(n.matrix.p, n.matrix.f)
  #
  for (iii in 1:2){
    ws = ws.list[[iii]]
    baseline.matrix = baseline.matrix.list[[iii]]
    observation.matrix = observation.matrix.list[[iii]]
    averages = averages.list[[iii]]
    n.matrix = n.matrix.list[[iii]]
    label=labels[iii]
    y1=Y1[iii] #0.25
    y2=Y2[iii] #0.43
    
    # tabulate color palette:
    
    if (is.null(n.matrix)){
      Palette = colorRampPalette(c("#1f77b4","#d62728"), alpha = F)(n_colors)
      colors = Palette[float2int(ws[,column], n_colors)]
      #
      Palette1 = viridisLite::viridis(n_colors)
      colors1 = Palette1[float2int(averages[,column], n_colors)]

      
    }else{
      Palette = colorRampPalette(c("#1f77b4","#d62728"), alpha = F)(n_colors)
      main_colors = Palette[float2int(ws[,column], n_colors)]
      
#      Alphas = c("F2","D9","BF","A6","8C","73","59","40","26","0D")
      Alphas = c("F2","D9","BF","A6","8C","73","59","40","40","40")
      opacities = Alphas[float2int(diff(Wilson.ci(ws[,column], n.matrix[,column])))]
      colors = sapply(1:nrow(ws), function(x){paste0(main_colors[x], opacities[x])})
      
        Palette1 = substring(viridisLite::viridis(n_colors), 1, 7)
        main_colors1 = Palette1[float2int(averages[,column], n_colors)]
        opacities1 = Alphas[float2int(diff(Wilson.ci(ws[,column], n.matrix[,column])))]
        colors1 = sapply(1:nrow(averages), function(x){paste0(main_colors1[x], opacities[x])})
    }
    if(is.null(n.matrix)){
      #take the index of the highest warning scores (mean over last few days),
      ws_mean = rowMeans(ws[,(column-3):column])
      idx = order(ws_mean, decreasing = T)[1:3]
    }else{
      #take the index of the highest warning scores, but only if statistically significant
      ss = n.matrix[,column] > 10
      ws_mean = rowMeans(ws[,(column-3):column])
      ws_min = min(ws_mean)
      ws_mean = ifelse(ss, ws_mean, ws_min)
      idx = order(ws_mean, decreasing = T)
      nn = min(3, length(idx))
      if(nn>0){
        idx = idx[1:nn]      
      }else{
        idx = NULL
      }
    }
    day = as.integer(colnames(ws)[column])
    x = (day-10):day
    y.label = ifelse(method=='rancover', 'w', 'Exceedance')
    # plot warning scores
    rx = range(x)
    plot(rx, c(0,1), type='n', xlab=NA, ylim=c(0,1), ylab = y.label,  xaxt='n')
    xtick<-seq(rx[1],rx[2],2)
    axis(side=1, at=xtick, labels = FALSE)
    text(x=xtick,  par("usr")[3],
         labels = int2daydate.lubridate(xtick), pos = 1, xpd = TRUE, cex = 0.9)
    
    
    title_label = ifelse(method == 'rancover', 'warning scores', 'exceedances')
    
    title(main=sprintf('Top %d %s for %s (last 3-day average)', length(idx), title_label, label), cex.main = 0.8)
    if(!is.null(idx)){
      if (!is.null(n.matrix)){
        for(I in 1:length(idx)){
          a = ws[idx[I], as.character(x)]
          b = n.matrix[idx[I], as.character(x)]
          ci = Wilson.ci(a,b)
          polygon(c(x, rev(x)), c(ci[1,] ,rev(ci[2,])), col = paste0(color_names[I], '50'), border = NA)      
        }
      }
      for (I in 1:length(idx)){
        lines(x, ws[idx[I],as.character(x)], col=color_names[I])
        points(x, ws[idx[I],as.character(x)], col=color_names[I])
      }
      legend('bottomleft', region_names[idx,], lty=c(1,1), pch=c(1,1),
             box.col = "grey",
             col = color_names[1:length(idx)])
    }
    rect(min(x)-0.1, 0.3, max(x) - 4, 1, col = "white", lwd=0)
    
    # plot map
    title = sprintf("%s - %s", label, int2daydate.lubridate(day))
    plot(map, col = colors, border = 'white', lwd = 0.1)
    if(!is.null(idx)){
      plot(map[idx], border = color_names[1:length(idx)], col=NA, add=T, lwd=0.5)
    }
    if(TMP > 900){
      mtext(side=3, title)
    }else{
      mtext(side=1, title)      
    }
    plot(map, col = colors1, border = 'white', lwd = 0.1)
    if(!is.null(idx)){
      plot(map[idx], border = color_names[1:length(idx)], col=NA, add=T, lwd=0.5)
    }
    if(TMP > 900){
      mtext(side=3, title)
      TMP = 0
    }else{
      mtext(side=1, title)      
    }
    
  }
  
  # plot baseline and counts in insets
  for(iii in 1:2){
    
    baseline.matrix = baseline.matrix.list[[iii]]
    observation.matrix = observation.matrix.list[[iii]]
    label=labels[iii]
    
    # print(range(c(baseline.matrix[idx, as.character(x)], observation.matrix[idx, as.character(x)])))
    Range = range(c(baseline.matrix[idx, as.character(x)], observation.matrix[idx, as.character(x)]))
    if (Range[2] <= 1){
      Range[2] = 1
    }
    y1=Y1[iii] #0.25
    y2=Y2[iii] #0.43
    par(fig=c(x1, x2, y1, y2),  new=T, mar=c(1.25, 2.2, 0.2, 0.2) )
    rx = range(x)
    plot(rx,
         Range,
         col='white', ylab='No. of cases', xaxt='n')
    # mtext(side=1, "Day", adj=0, cex=0.6)
    xtick<-seq(rx[1],rx[2], 4)
    axis(side=1, at=xtick, labels = FALSE)
    text(x=xtick,  par("usr")[3],
         labels = int2daydate.lubridate(xtick), pos = 1, xpd = TRUE, cex = 0.9)
    mtext(side=3, "No. of cases", padj=0, cex=0.5)
    
    if(!is.null(idx)){
      for (I in 1:length(idx)){
        # print(idx[I])
        lines(x, baseline.matrix[idx[I], as.character(x)], col=color_names[I], lty=1)
        lower = qpois(0.05, baseline.matrix[idx[I], as.character(x)])
        upper = qpois(0.95, baseline.matrix[idx[I], as.character(x)])
        polygon(c(x, rev(x)), c(lower ,rev(upper)), col = paste0(color_names[I], '50'), border = NA)
        lines(x, observation.matrix[idx[I],as.character(x)], col=color_names[I], lty=2)
      }      
    }
    legend('topleft', c('Baseline', 'Observations'), lty=c(1,2), cex=0.8, box.col = 'grey')
  }
  y.label = ifelse(method == 'rancover', 'w', 'Exceedance')
  
  # plot color Palette
  par(fig=c(0.40, 0.46, 0.8, 0.9), new=T)
  par(mar=c(0.2, 2.2, 0.2, 0.2))
  plot(c(0, 1), c(0,1), ylab = y.label, yaxt='n',
       xaxt='n', col='white', frame.plot=F)
  axis(2, at = c(0,0.5,1), lwd=0.5)
  
  
  if(!is.null(n.matrix)){
    mtext(side=1, 'CI width', line = 2, cex = 0.6)
    Cb = sapply(1:10, function(i){sapply(1:10, function(x){paste0(Palette[x], Alphas[i])}) })
    axis(1, at = c(0,0.5,1), lwd=0.5, cex=0.3)
    mtext(side=2, y.label, line = 2)
    
  }else{
    Cb = matrix(rep(Palette, 10), nrow = length(Palette), ncol = 10)
    mtext(side=2, y.label, line = 2, cex = 0.6)
    
  }
  rasterImage(Cb[seq(length(Palette), 1, -1),], 0, 0, 1, 1)
  rect(0, 0, 1, 1, lwd=0.5)
  
  
  # plot color Palette1
  par(fig=c(0.40 + 0.333, 0.46 + 0.333, 0.8, 0.9), new=T)
  par(mar=c(0.2, 2.2, 0.2, 0.2))
  plot(c(0, 1), c(0,1), yaxt='n',
       xaxt='n', col='white', frame.plot=F)
  axis(2, at = c(0,0.5,1), lwd=0.5)
  mtext(side=2, 'Avrg. ratio', line = 2, cex = 0.6)
  if(!is.null(n.matrix)){
    mtext(side=1, 'CI width', line = 2, cex = 0.6)
    axis(1, at = c(0,0.5,1), lwd=0.5, cex=0.3)
    Cb = sapply(1:10, function(i){sapply(1:10, function(x){paste0(Palette1[x], Alphas[i])}) })
  }else{
    Cb = matrix(rep(Palette1, 10), nrow = length(Palette), ncol = 10)
  }
  rasterImage(Cb[seq(length(Palette1), 1, -1),], 0, 0, 1, 1)
  rect(0, 0, 1, 1, lwd=0.5)
  
  par(fig=c(0.55, 0.8, 0.23 + 0.1, 0.4 + 0.2), new=T)
  par(mar=c(1.2, 2.2, 0.2, 0.2))
  plot_all(observation.matrix.p,  observation.matrix.f, as.integer(colnames(ws.p)[column]))
}

