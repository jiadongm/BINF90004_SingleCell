### These functions are useful for visualising dimension reduction results
# Author: Dr Jiadong Mao, Melbourne Integrative Genomics & School of Mathematics and Statistics, The Univeristy of Melbourne
# Github: https://github.com/jiadongm


## Feature mean plots
marginMeanPlot <- function(plot_dat, log = T, byvar = F){
  
  plot_dat <- as.matrix(plot_dat)
  if(byvar){
    marginMeans <- apply(plot_dat, 2, var)
  } else {
    marginMeans <- colMeans(plot_dat)
  }
  marginMeans <- sort(marginMeans, decreasing = T)
  
  if(log){
    plot(log10(marginMeans + 1))
  } else {
    plot(marginMeans)
  }
  
  return(marginMeans)
}


## Marginal density plots
marginPlot <- function(plot_dat, featIdx, 
                       plotNrows = NULL, plotNcols = NULL,
                       xlim = NULL,
                       marginMeans = NULL,
                       groupKey = NULL,
                       pointSize = 3, pointAlpha = NULL,
                       manualCol = NULL, manualAlpha = NULL,
                       omitZeros = F,
                       graph.out = T){
  
  outPlots <- vector("list", length(featIdx))
  for(i in 1:length(featIdx)){
    
    var2plot <- featIdx[i]
    
    xx <- plot_dat[,var2plot]
    if(omitZeros) xx <- xx[xx != 0]
    # Plug-in bw selection
    bw <- try(hpi(xx, binned = F), silent = T)
    
    if(inherits(bw, "try-error")){
      
      if(omitZeros){
        # If bw failed (due to eg sparsity), plot histogram
        out <- plot_dat %>%
          filter(!! sym(var2plot) != 0) %>%
          ggplot(aes(x = !! sym(var2plot))) + 
          geom_histogram()
      } else {
        out <- plot_dat %>%
          ggplot(aes(x = !! sym(var2plot))) + 
          geom_histogram()
      }
      
    } else {
      
      if(omitZeros){
        out <- plot_dat %>%
          filter(!! sym(var2plot) != 0) %>%
          ggplot(aes(x = !! sym(var2plot))) + 
          geom_density(bw = 'SJ')
      } else {
        out <- plot_dat %>%
          ggplot(aes(x = !! sym(var2plot))) + 
          geom_density(bw = 'SJ')
      }
      
      
    }
    
    if(is.null(groupKey)){ # if groupkey has been specified
      
      out <-
        out +
        geom_jitter(aes(y=0), 
                    height = diff(layer_scales(out)$y$range$range)/20, 
                    size = pointSize,
                    shape = 16)
      
    } else {
      
      out <-
        out +
        geom_jitter(aes(y=0, colour = !! sym(groupKey), alpha = !! sym(groupKey)), 
                    height = diff(layer_scales(out)$y$range$range)/20, 
                    size = pointSize,
                    shape = 16)
      
    }
    
    if(!is.null(marginMeans)){
      
      out <- out +
        labs(x = paste0(var2plot, "(", signif(marginMeans[var2plot], 4), ")"))
      
    }
    
    if(!is.null(manualCol)){
      out <- out + 
        scale_color_manual(values = manualCol)
    }
    
    if(is.null(manualAlpha)){
      
      Ngroups <- length(unique(plot_dat[,groupKey]))
      
      if(is.null(pointAlpha)) pointAlpha <- 1
      
      out <- out +
        scale_alpha_manual(values = rep(pointAlpha, Ngroups))
      
    } else {
      out <- out +
        scale_alpha_manual(values = manualAlpha)
    }
    
    if(!is.null(xlim)){
      out <- out +
        scale_x_continuous(limits = xlim)
    }
    
    outPlots[[i]] <- out
    
  }
  
  if(graph.out){
    print(ggarrange(plotlist = outPlots, 
                    nrow = plotNrows, 
                    ncol = plotNcols,
                    common.legend = T))
  }
  
  
  return(outPlots)
  
}




#' Score matrix plot
#'
#' Density and pairwise scatter plots for visualising score matrices.
#'
#' @param plot_dat Data frame like objects. Score matrix to plot.
#' @param max_ncomp Default NULL. Number of first components to plot. If specified, will override comp_idx.
#' @param comp_idx Default NULL. Indices of components to plot.
#' @param groupKey Column name of plot_dat storing group information.
#' @param showGraph Logic. Whether to print matrix plot or not.
#'
#' @examples
#' 
#'
#' @export
matrixPlot <- function(plot_dat, max_ncomp = NULL, comp_idx = NULL, groupKey = NULL, 
                       pointAlpha = NULL,
                       pointSize = 1, 
                       manualCol = NULL, manualAlpha = NULL){
  
  if(is.null(max_ncomp) & is.null(comp_idx)){
    stop("Need to specify either max_ncomp or comp_idx.")
  }
  
  if(!is.null(max_ncomp)) comp_idx <- 1:max_ncomp
  
  
  if(!all(paste0("comp", comp_idx) %in% colnames(plot_dat))){
    
    missingComps <-
      paste0("comp", comp_idx)[!(paste0("comp", comp_idx) %in% colnames(plot_dat))]
    stop(paste0("These components are missing from plot_dat: ", missingComps))
    
  }
  
  
  if(length(comp_idx) >= 3){
    
    # Density plots on diagonal
    out_diag <- vector("list", length(comp_idx))
    
    for(comp_i in 1:length(comp_idx)){
      
      var2plot <- paste0("comp", comp_idx[comp_i])
      
      p <- 
        plot_dat %>%
        ggplot(aes(x = !! sym(var2plot))) +
        geom_density(bw = "SJ") +
        scale_y_continuous(name = NULL)
      
      
      
      if(is.null(groupKey)){
        
        p <- p + 
          geom_jitter(aes(y = 0), 
                      height = diff(layer_scales(p)$y$range$range)/20) 
        out_diag[[comp_i]] <- p
        
      } else {
        
        p <- p + 
          geom_jitter(aes(y = 0, colour = !! sym(groupKey)), 
                      height = diff(layer_scales(p)$y$range$range)/20) 
        
        if(!is.null(manualCol)){
          p <- p + 
            scale_color_manual(values = manualCol)
        }
        
        out_diag[[comp_i]] <- p +
          theme(legend.position = "none") 
        
        
        
      }
      
      
      
    }
    
    # Get the legend 
    if(!is.null(groupKey)) p_legend <- cowplot::get_legend(p)
    
    
    # Scatter plots on non-diagonal
    combs <- combn(comp_idx, 2) %>% t()
    out_nondiag <- vector("list", nrow(combs))
    for(comb in 1:nrow(combs)){
      
      var1 <- paste0("comp", combs[comb, 1])
      var2 <- paste0("comp", combs[comb, 2])
      
      if(is.null(groupKey)){
        
        p <-
          plot_dat %>%
          ggplot(aes(x = !! sym(var1), y = !! sym(var2))) +
          geom_point(size = pointSize) +
          theme(legend.position = "none", 
                axis.ticks = element_blank()) +
          scale_x_continuous(name = NULL, labels = NULL) +
          scale_y_continuous(name = NULL, labels = NULL)
        
        
      } else {
        
        p <-
          plot_dat %>%
          ggplot(aes(x = !! sym(var1), y = !! sym(var2))) +
          geom_point(aes(colour = !! sym(groupKey)), size = pointSize) +
          theme(legend.position = "none", 
                axis.ticks = element_blank()) +
          scale_x_continuous(name = NULL, labels = NULL) +
          scale_y_continuous(name = NULL, labels = NULL)
        
        if(!is.null(manualCol)){
          p <- p + 
            scale_color_manual(values = manualCol)
        }
      }
      
      out_nondiag[[comb]] <- p
      
    }
    
    
    # Arrange plots
    out <- c(out_nondiag, out_diag)
    if(!is.null(groupKey)){
      out[[length(out) + 1]] <- p_legend
    }
    layoutM <- matrix(NA, length(comp_idx), length(comp_idx))
    layoutM[lower.tri(layoutM, diag = F)] <- 1:nrow(combs)
    diag(layoutM) <- (nrow(combs)+1):(nrow(combs)+length(comp_idx))
    if(!is.null(groupKey)) layoutM[ceiling(length(comp_idx)/2)-1, ceiling(length(comp_idx)/2)+1] <- length(out)
    p <- grid.arrange(grobs = out, layout_matrix = layoutM)
    
    print(p)
    
  } else if (length(comp_idx) == 1){
    
    var2plot <- paste0("comp", comp_idx)
    
    ## Only density plot is given
    out <- 
      plot_dat %>%
      ggplot(aes(x = !! sym(var2plot))) +
      geom_density(bw = "SJ") 
    if(is.null(groupKey)){
      out <-
        out +
        geom_jitter(aes(y=0), 
                    height = diff(layer_scales(out)$y$range$range)/20, 
                    size = pointSize,
                    shape = 16)
    } else {
      out <-
        out +
        geom_jitter(aes(y=0, colour = !! sym(groupKey), alpha = !! sym(groupKey)), 
                    height = diff(layer_scales(out)$y$range$range)/20, 
                    size = pointSize,
                    shape = 16)
    }
    
    if(!is.null(manualCol)){
      out <- out + 
        scale_color_manual(values = manualCol)
    }
    
    if(is.null(manualAlpha)){
      
      Ngroups <- length(unique(plot_dat[,groupKey]))
      
      if(is.null(pointAlpha)) pointAlpha <- 1
      
      out <- out +
        scale_alpha_manual(values = rep(pointAlpha, Ngroups))
      
    } else {
      out <- out +
        scale_alpha_manual(values = manualAlpha)
    }
    
    
  } else {
    
    ## Plot 2 components
    
    var1 <- paste0("comp", comp_idx[1])
    var2 <- paste0("comp", comp_idx[2])
    
    if(is.null(groupKey)){
      
      p <-
        plot_dat %>%
        ggplot(aes(x = !! sym(var1), y = !! sym(var2))) +
        geom_point(size = pointSize) 
      
      
    } else {
      
      p <-
        plot_dat %>%
        ggplot(aes(x = !! sym(var1), y = !! sym(var2))) +
        geom_point(aes(colour = !! sym(groupKey)), size = pointSize) 
      
      if(!is.null(manualCol)){
        p <- p + 
          scale_color_manual(values = manualCol)
      }
    }
    
    print(p)
    out <- p
    
  }
  
  return(out)
  
}


## Sankey plot
plotSankey <- function(class1, class2, fontsize = 12, calcErrs = F){
  confTab <- table(class1, class2)
  # Plot Sankey Network Diagram
  colLinks <- c('source', 'target', 'value')
  sankeyLinks <- matrix(ncol = length(colLinks), nrow = 0)
  colnames(sankeyLinks) <- colLinks
  for(i in 1:nrow(confTab)){
    for(j in 1:ncol(confTab)){
      if(confTab[i, j]>0){
        sankeyLinks <- rbind(sankeyLinks,
                             # - 1 is for zero index
                             c(i - 1, nrow(confTab) + j - 1, confTab[i, j]))
      }
    }
  }
  sankeyNodes <- data.frame(name = c(rownames(confTab), colnames(confTab)))
  p <- networkD3::sankeyNetwork(Links = as.data.frame(sankeyLinks),
                                Nodes = sankeyNodes, NodeID = 'name',
                                Source = 'source', Target = 'target',
                                Value = 'value',
                                fontSize = fontsize)
  print(p)
  if(calcErrs) return(classErr(class1, class2)$err)
}

# Sankey diagram for 3 classification results
plotSankey3 <- function(class1, class2, class3, add = T, fontsize = 12, calcErrs = F){
  if(calcErrs){
    er_re12 <- classErr(class1, class2)
    er_re23 <- classErr(class3, class2)
  }
  if(add){
    class1 <- paste0(class1, "-")
    class2 <- paste0(class2, "_")
  }
  
  
  # Define a flow incidence matrix, from rows to cols
  rowNams <- colNams <- unique(c(class1, class2, class3))
  
  incidMat <- matrix(0, nrow = length(rowNams), ncol = length(colNams),
                     dimnames = list(rowNams, colNams))
  
  
  confTab12 <- table(class1, class2)
  confTab23 <- table(class2, class3)
  
  confTab <- confTab12
  for(i in 1:nrow(confTab)){
    fromClass <- rownames(confTab)[i]
    for(j in 1:ncol(confTab)){
      toClass <- colnames(confTab)[j]
      incidMat[fromClass, toClass] <- confTab[i, j]
    }
  }
  confTab <- confTab23
  for(i in 1:nrow(confTab)){
    fromClass <- rownames(confTab)[i]
    for(j in 1:ncol(confTab)){
      toClass <- colnames(confTab)[j]
      incidMat[fromClass, toClass] <- confTab[i, j]
    }
  }
  
  # Following code from https://r-graph-gallery.com/321-introduction-to-interactive-sankey-diagram-2.html
  sankeyLinks <- incidMat %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="source") %>% 
    tidyr::gather(key="target", value="value", -1) %>%
    dplyr::filter(value != 0)
  
  sankeyNodes <- data.frame(
    name=c(as.character(sankeyLinks$source), as.character(sankeyLinks$target)) %>% 
      unique())
  
  sankeyLinks$IDsource <- match(sankeyLinks$source, sankeyNodes$name)-1 
  sankeyLinks$IDtarget <- match(sankeyLinks$target, sankeyNodes$name)-1
  
  p <- 
    networkD3::sankeyNetwork(Links = sankeyLinks, Nodes = sankeyNodes,
                  Source = "IDsource", Target = "IDtarget",
                  Value = "value", NodeID = "name", 
                  sinksRight = FALSE, fontSize = fontsize)
  print(p)
  
  if(calcErrs) return(list(er12 = er_re12$err, er23 = er_re23$err))
}