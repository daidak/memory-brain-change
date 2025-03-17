
brain_plots <- function(df, fill, significance_column = NULL, scale_name = "", na_color = "ivory", clrscale = "viridis", color.group = NULL) {
  library(ggseg)
  library(ggsegDesterieux)
  df <- df %>% 
    rename(label = features) %>%
    ungroup()
  
  # Filter based on significance if a significance_column is provided
  # Significant_colum; p.value, fdr, p.value.age, fdr.age etc.
  if (!is.null(significance_column) && !significance_column %in% names(df)) {
    stop("wrong variable name in significance_column (e.g. fdr, p.value, fdr.age, p.value.age.")
  }
  
  if (!fill %in% names(df)) {
    stop("wrong variable name for fill. (e.g. estimate, log, log_fdr..")
  }
  
  # Filter based on significance if a significance_column is provided
  if (!is.null(significance_column)) {
    df <- df %>% filter(.[[significance_column]] < 0.05)
  }
  
  # infinite values can appear due to permutation testing 
  df[[fill]][is.infinite(df[[fill]])] <- 4
  df[[fill]][is.infinite(df[[fill]])] <- -4
  
  
  # Calculate min and max fill values
  min_fill <- min(df[[fill]], na.rm = TRUE)
  max_fill <- max(df[[fill]], na.rm = TRUE)
  
  if (clrscale == "viridis") {
  colors <- viridis::viridis(100)
  } else if (clrscale == "fs") {
    my_palette <- colorRampPalette(c("lightblue", "blue", "white", "red", "yellow"))
    colors <- my_palette(100)
    max_fill <- max(c(abs(min_fill), abs(max_fill)))
    min_fill <-  -max_fill
  }
  
  # Plot 1
  g1 <- ggplot(df) +
    geom_brain(atlas = desterieux,
               position = position_brain(hemi ~ side),
               aes(fill = !!sym(fill)), color = "black", size = 0.8) + 
    scale_fill_gradientn(colors = colors, 
                         limits = c(min_fill, max_fill),
                         name = scale_name, 
                         na.value = na_color)+
    theme_brain2(text.colour = "black") +
    labs(title = "cortical thickness") +
    theme(legend.position = 'none',
          title = element_text(hjust = 0.5, color = "black"), 
          text=element_text(family="Arial")) 
  
  
  # Plot 2
  g2 <- ggplot(df) +
    geom_brain(atlas = aseg,
               side = "coronal",
               aes(fill = !!sym(fill)), color = "black", size = 0.8) + 
    scale_fill_gradientn(colors = colors, 
                         limits = c(min_fill, max_fill),
                         name = scale_name) +
    theme_brain2(text.colour = "black")+
    labs(title = "subcortical volume") +
    theme(legend.position = 'bottom',
          legend.title = element_text(size = 14),    # Increase legend title size
          legend.text = element_text(size = 12),     # Increase legend text size
          legend.key.size = unit(1.5, "lines"), 
          title = element_text(hjust = 0.5, color = "black"),
          text=element_text(family="Arial")) 
  
  
  # Return the plots
  return(list(g1 = g1, g2 = g2))
}


wrapper_brain_plots <- function(df, fill, na_color = "ivory", mode = "main") {
  library(ggseg)
  library(ggsegDesterieux)
  
  if(mode == "main") {
  df <- df %>% 
    rename(label = features) %>%
    ungroup()  %>% 
    mutate(color.groupF = if_else(gamm.w.pval.bootstrap > 0.05, 0,
                                  if_else(gamm.w.pval.bootstrap.fdr > .05, 1, 2)) %>% as.factor(),
           color.name = if_else(color.groupF == 2, "orange", "ivory"), 
           size.I = if_else(color.groupF == 2, 1.1, .5))
  } else if (mode == "int") {
    df <- df %>% 
      rename(label = features) %>%
      ungroup()  %>% 
      mutate(color.groupF = if_else(gamm.int.w.pval.bootstrap > 0.05, 0,
                                    if_else(gamm.int.w.pval.bootstrap.fdr > .05, 1, 2)) %>% as.factor(),
             color.name = if_else(color.groupF == 2, "orange", "ivory"), 
             size.I = if_else(color.groupF == 2, 1.1, .5))
  }
  
  # Calculate min and max fill values
  min_fill <- min(df[[fill]], na.rm = TRUE)
  max_fill <- max(df[[fill]], na.rm = TRUE)
  
  colors <- viridis::viridis(100)
  
  
  # Plot 1
  g1 <- ggplot(df, color = "black") +
    geom_brain(atlas = desterieux,
               position = position_brain(hemi ~ side),
               aes(fill = !!sym(fill), 
                   color = I(color.name), 
                   size = I(size.I), 
                   alpha = color.groupF)) + 
    scale_fill_gradientn(colors = colors, 
                         limits = c(min_fill, max_fill),
                         name = "", 
                         na.value = na_color)+
    theme_brain2(text.colour = "black") +
    labs(title = "cortical thickness") +
    theme(legend.position = 'none',
          title = element_text(hjust = 0.5, color = "black"),
          text=element_text(family="Arial")) + 
    #scale_color_manual(values = c("ivory","ivory", "orange")) + 
    scale_alpha_manual( values = c(1,1,1))
  
  
  # Plot 2
  g2 <- ggplot(df, color = "black") +
    geom_brain(atlas = aseg,
               side = "coronal",
               aes(fill = !!sym(fill), 
                   color = I(color.name), 
                   size = I(size.I), 
                   alpha = color.groupF))+ 
    scale_fill_gradientn(colors = colors, 
                         limits = c(min_fill, max_fill),
                         name = "") +
    theme_brain2(text.colour = "black")+
    labs(title = "subcortical volume") +
    theme(legend.position = 'bottom',
          legend.title = element_text(size = 14),    # Increase legend title size
          legend.text = element_text(size = 12),     # Increase legend text size
          legend.key.size = unit(1.5, "lines"), 
          title = element_text(hjust = 0.5, color = "black"),
          text=element_text(family="Arial")) + 
    #scale_color_manual(values = c("ivory","ivory", "orange"), guide = "none") + 
    scale_alpha_manual( values = c(1,1,1), guide = "none")
  
  
  # Return the plots
  return(list(g1 = g1, g2 = g2))
}


inferiority.test = function(lm, con, rhs, alternative) {
  x = glht(model = lm,
           linfct = matrix(con, nrow = 1),
           rhs = rhs,
           alternative = alternative)
  return(summary(x)$test$pvalues[[1]])
}


plot_consensus_matrix = function(mod, clusdir) {
  library(ComplexHeatmap)
  ccl <- list()
  x <- c(
    "#FDB4C7",  # Pinkish
    "#B5EAD7",  # Greenish
    "#A7D8F0",  # Blueish
    "#FFECB3",  # Yellowish
    "#FDCFE8",  # Light pink
    "#B5D0EA",  # Soft blue-green
    "#C2F2E5",  # Aqua
    "#FFDAC1",  # Peach
    "#FFE5A3",  # Light yellow-orange
    "#D3E8F5"   # Sky blue
  )
  names(x) <- as.character(seq(1,10,by=1))
  for (i in seq(2,10)){
    # get cc matrix and labels
    ccmatrix <- mod$realdataresults[[i]]$consensus_matrix
    annon <- mod$realdataresults[[i]]$ordered_annotation
    # do heatmap
    n <- 10
    seq <- rev(seq(0,255,by=255/(n)))
    palRGB <- cbind(seq,seq,255)
    mypal <- rgb(palRGB,maxColorValue=255)
    #my_palette <- colorRampPalette(c("white", "lightblue", "blue"))
    #mypal <- my_palette(100)
    ha = HeatmapAnnotation(
      df= data.frame(Cluster=as.character(annon[,1])), col = list(Cluster=x), show_legend = F, show_annotation_name = F)
    ra = rowAnnotation(
      df= data.frame(Cluster=as.character(annon[,1])), col = list(Cluster=x), show_legend = F, show_annotation_name = F)
    ccl[[i]] <- Heatmap(ccmatrix, 
                        name = "Consensus_index", 
                        top_annotation = ha,
                        left_annotation = ra,
                        col=mypal, 
                        show_row_dend = FALSE,
                        show_column_dend = FALSE, 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE,
                        show_column_names = FALSE, 
                        show_heatmap_legend = TRUE,
                        border_gp = gpar(col = "black", lwd = 1.5),
                        border = T, 
                        rect_gp = gpar(col = "grey90", lwd = 1),
                        heatmap_legend_param = list(
                          legend_direction = "horizontal"))
  }
  save(ccl, file = file.path(clusdir, "consensus.heatmaps.rda"))
  return(ccl)
}






wrapper_plot_trajectories_main= function(mega, outdir.mega.main){
  df.plot.deriv = 
    mega %>%  
    filter(gamm.w.pval.bootstrap.fdr < .05)  %>% 
    select(c(features, datamain)) %>% 
    unnest() %>% 
    dplyr::select(features, 
                  .derivative, 
                  .se, 
                  .lower_ci, 
                  .upper_ci, 
                  delta_brain) %>% 
    mutate(type = "derivative") %>% 
    rename(.estimate = .derivative)
  
  df.plot = mega %>%  
    filter(gamm.w.pval.bootstrap.fdr < .05)  %>% 
    select(c(features, datamainsmooth)) %>% 
    unnest() %>% 
    mutate(.lower_ci = .estimate - 1.96*.se, 
           .upper_ci = .estimate + 1.96*.se) %>% 
    select(features, 
           .estimate, 
           .se,
           .lower_ci, 
           .upper_ci, 
           delta_brain) %>% 
    mutate(type = "smooth")
  
  df.plot = rbind(df.plot, df.plot.deriv) 
  
  df.plot = 
    df.plot %>% 
    filter(features %in% 
             c("Left-Caudate", 
               "Left-Hippocampus" , 
               "Left-Inf-Lat-Vent", 
               "lh_G_oc-temp_med-Parahip")) %>% 
    mutate(type = 
             factor(type, levels = c("smooth", "derivative"), 
                    labels = c(expression(Delta ~ "Memory"), 
                               expression(Delta ~ "Memory/yr"))))
  nlabels = c("left\ncaudate",
              "left\nhippocampus", 
              "left inferior\nlateral ventricle", 
              "left \nparahippocampal gyrus") 
  
  df.plot = df.plot %>% filter(between(delta_brain, -2.5, 2.5))
  gs = 
    ggplot(df.plot, 
           aes(delta_brain,.estimate, group = features, fill = features)) + 
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey50") + 
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = .6)+ 
    geom_line(size = 1.5) +
    facet_grid(rows = vars(type), cols = vars(features), scales = "free_y", switch = "y", labeller = label_parsed) + 
    theme_classic() +
    scale_fill_manual(values =   c("#FDB4C7", "#B5EAD7", "#A7D8F0", "#FFECB3"), 
                      labels = nlabels) +
    theme(text=element_text(family="Arial", size = 16),
          strip.placement = "outside", 
          legend.position = "bottom", 
          legend.title = element_blank(), 
          strip.background = element_blank(), 
          strip.text.x = element_blank(),
          axis.title.y = element_blank(), 
          axis.title.x = element_text(face = "bold", size = 16),
          strip.text = element_text(size = 16, face = "bold"),
          axis.text = element_text(size = 14, face = "bold"), 
          legend.text = element_text(size = 14, face = "bold"),
          panel.spacing.y = unit(1, "lines")) +
    labs(x = expression(Delta~"Brain")) 
  ggsave(file.path(outdir.mega.main, "main.trajectories.selected.png"), 
         plot = gs, 
         width = 9,
         height = 6)  
  return(gs)
}


wrapper_plot_trajectories_int= function(mega){
  
  df.plot = 
    mega %>%  
    filter(gamm.int.w.pval.bootstrap.fdr < .05)  %>% 
    select(c(features, datainteraction)) %>% 
    unnest() %>% 
    mutate(lower_ci = estimate - 1.96*se, 
           upper_ci = estimate + 1.96*se) 
  
  
  df.plot = 
    df.plot %>% 
    filter(features %in% 
             c("Left-Hippocampus" , 
               "lh_G_insular_short",
               "Right-Inf-Lat-Vent"
               )) %>% 
    filter(xage %in% c(40,50,60,70,80)) 
  
  nlabels = c("left\nhippocampus",
              "left short\ninsular gyri", 
              "right inferior\nlateral ventricle")
  df.plot = df.plot %>% filter(between(delta_brain, -2.5, 2.5))
  gs = 
    ggplot(df.plot, aes(delta_brain,estimate, group = features, fill = features)) + 
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey50") + 
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = .6)+ 
    geom_line(size = 1.5) +
    geom_line(mapping = aes(y = derivative*2), color = "lightblue", linewidth = 1.5, alpha = .5) +
    facet_grid(cols = vars(xage), rows = vars(features),  scales = "free_y", switch = "y", labeller = label_parsed) + 
    theme_classic() +
    scale_fill_manual(values =   c("#FDB4C7", "#B5EAD7","#FFECB3"), 
                      labels = nlabels) +
    scale_y_continuous(
      name = expression(Delta~"Memory"), 
      sec.axis = sec_axis(~ . / 2, expression(Delta~"Memory/year"))
    ) +
    theme(strip.placement = "outside", 
          legend.position = "bottom", 
          legend.title = element_blank(), 
          strip.background = element_blank(), 
          strip.text.y = element_blank(),
          #axis.title.y = element_blank(), 
          axis.title.y = element_text(face = "bold", size = 16),,
          axis.title.x = element_text(face = "bold", size = 16),
          strip.text = element_text(size = 16, face = "bold"),
          axis.text = element_text(size = 14, face = "bold"), 
          legend.text = element_text(size = 14, face = "bold"), 
          text=element_text(family="Arial", face = "bold"), 
          panel.spacing.y = unit(1, "lines")) +
    labs(x = expression(Delta~"Brain"))
  
  
  
  ggsave(file.path(outdir.mega.int, "int.trajectories.selected.png"), 
         plot = gs, 
         width = 9,
         height = 6.2)  
  return(gs)
}

wrapper_density_tile_plot = function(wd, outdir) {
  outdir.mega.dim = file.path(outdir, "dimensionality")
  load(file.path(wd,"mega","consensus_cluster" ,"cluster_mod.rda"))
  
  # CORR PLOT
  
  corM = cor(m3c$M)
  distM = corM[upper.tri(corM, diag = F)]
  clim = range(corM)
  dlim = range(distM)
  annon <- m3c$mod$realdataresults[[8]]$ordered_annotation
  idx = match(gsub("\\.","-",rownames(annon)), rownames(corM))
  corM = corM[idx,idx]
  
  df.changenames = retrieve_orig2plot_names()
  old.names = gsub("\\.","-",rownames(annon))
  new.names = df.changenames$plot.names[which(df.changenames$orig.names %in% old.names )]
  
  
  colnames(corM) = new.names
  rownames(corM) = new.names
  
  distM = reshape::melt(corM) %>% 
    mutate(value = if_else(value == 1, NaN, value))
  
  my_palette <- colorRampPalette(c("lightblue", "blue", "#dfefff", "white","lightpink", "red", "yellow"))
  colors <- my_palette(100)
  max_fill <- max(c(abs(dlim), abs(dlim)))
  min_fill <-  -max_fill
  
  gs= 
    ggplot(distM, aes(x = X1, y = X2, fill = value)) + 
    geom_tile(color = "white",
              lwd = 1,
              linetype = 1) +
    scale_fill_gradientn(name = "Corr.",
                         colors = colors, 
                         limits = c(min_fill, max_fill),
                         na.value = "ivory") +
    coord_fixed() +
    #theme_classic() +
    theme( panel.border = element_rect(colour = "black", fill=NA, size=2), 
           axis.title = element_blank(), 
           axis.text.y = element_text(size = 14, face = "bold"),
           axis.text.x = element_blank(), 
           axis.ticks.x = element_blank(),
           legend.text = element_text(size = 14, face = "bold"), 
           text=element_text(family="Arial"),
    ) 
  ggsave(file.path(outdir.mega.dim, "correlation_matrix.png"), 
         plot = gs)  
  
  # DENSITY PLOT
  distM= corM[upper.tri(corM, diag = F)]
  df.density = data.frame("X" = "slopeCorr", 
                          distM = distM)
  
  gs1 = ggplot(df.density, aes(x = distM)) + 
    geom_density(aes(fill = "#A7D8F0"), size = 1.5, color = "black",alpha = 0.5) + 
    scale_fill_identity() + 
    theme_void() +  
    scale_x_continuous(position = "top", limits = dlim) +
    theme(text = element_text(family = "Arial"), 
          aspect.ratio = .25,
          legend.position = "none",
          #axis.line.x.top = element_line(color = "black"), 
          axis.text.x.top =  element_text(size = 14, face = "bold"))+
    scale_y_reverse()
  
  ggsave(file.path(outdir.mega.dim, "densityplot.png"), 
         plot = gs1)  
  
  
  
  gs.mix = gs /gs1
  ggsave(file.path(outdir.mega.dim, "combinedplot.png"), 
         plot = gs.mix)  
  
  out = list(gs = gs, 
             gs1 = gs1, 
             gs.mix = gs.mix)
  return(out)
}


wrapper_consensus_cluster_plot = function(wd, outdir) {
  outdir.mega.dim = file.path(outdir, "dimensionality")
  load(file.path(wd,"mega","consensus_cluster", "cluster_mod.rda"))
  ccl = plot_consensus_matrix(m3c$mod, file.path(wd,"mega","consensus_cluster"))
  
  png(file.path(outdir.mega.dim, "consensus_cluster.jpeg"),4000,4000, res = 300) 
  draw(ccl[[8]], merge_legend = TRUE, heatmap_legend_side = "bottom")
  dev.off()
  
  
  
  ###
  annon <- m3c$mod$realdataresults[[8]]$ordered_annotation
  annon$label = gsub("\\.","-",rownames(annon))
  annon$ccf = as.factor(annon$consensuscluster)
  annon = annon %>% 
    mutate(label = if_else(label == "Left-Inf-Lat-Vent", "Left-Lateral-Ventricle", label))
  
  
  
  mycols = c(
    "1" = "#FDB4C7",  # Pinkish
    "2" = "#B5EAD7",  # Greenish
    "3" = "#A7D8F0",  # Blueish
    "4" = "#FFECB3",  # Yellowish
    "5" = "#FDCFE8",  # Light pink
    "6" = "#B5D0EA",  # Soft blue-green
    "7" = "#C2F2E5",  # Aqua
    "8" = "#FFDAC1"  # Peach
  )
  c1 <- ggplot(annon) + 
    geom_brain(atlas = desterieux, 
               position = position_brain(hemi ~ side),
               aes(fill = ccf), colour = "black", size = 0.6) +
    scale_fill_manual(name = "clusters", values =  mycols,   na.value = "ivory")+
    theme_brain2(text.colour = "black") + 
    theme(title = element_blank(),
          text=element_text(family="Arial"), 
          legend.position = "none")
  
  c2 <- ggplot(annon) + 
    geom_brain(atlas = aseg, 
               side = "coronal", 
               aes(fill =ccf), colour = "black", size = 0.7) +
    scale_fill_manual(name = "clusters", values = mycols,   na.value = "ivory")+
    theme_brain2(text.colour = "black") + 
    theme(title = element_blank(),
          text=element_text(family="Arial"), 
          legend.position = "none")
  
  
  
  legend_only <- ggplot(annon, aes(x = label, y = ccf,group = ccf, fill = ccf)) +
    geom_point(size = 4, shape = 21) +
    scale_fill_manual(name = "Cluster", values = mycols) +
    theme_void() +  # Remove all plot elements
    theme(text=element_text(family="Arial", size = 16), 
          legend.position = "bottom") +
    labs(color = "Factor Levels")
  
  
  
  ggsave(file.path(outdir.mega.dim, "cluster.cth.png"), 
         plot = c1)  
  ggsave(file.path(outdir.mega.dim, "cluster.cth.jpeg"), 
         plot = c2)  
  ggsave(file.path(outdir.mega.dim, "legend.jpeg"), 
         plot = legend_only)  
  
  out = list(
    gs1 = ccl,
    gs2 = c1, 
    gs3 = c2, 
    gs4 = legend_only
  )
  return(out)
}




retrieve_orig2plot_names = function() {
orig2plot_names = 
tribble(~orig.names, ~plot.names, 
        "lh_G_oc-temp_med-Parahip", "l parahippocampal gyrus",
        "lh_G_temporal_inf", "l inferior temporal gyrus",
        "lh_S_pericallosal","l pericallosal sulcus",
        "rh_G_Ins_lg_and_S_cent_ins","r long insular gyrus",
        "rh_G_pariet_inf-Supramar","r supramarginal gyrus",
        "rh_G_rectus", "r gyrus rectus",
        "rh_G_subcallosal","r subcallosal gyrus",
        "rh_G_temp_sup-Plan_polar","r planum polare",
        "rh_G_temporal_inf","r inferior temporal gyrus",
        "rh_S_collat_transv_ant","r ant. transverse col. sulcus",
        "rh_S_temporal_transverse", "r transverse temporal sulcus",
        "Left-Inf-Lat-Vent","l inf. lat. vent.",
        "Left-Thalamus-Proper", "l thalamus",
        "Left-Caudate","l caudate", 
        "Left-Hippocampus", "l hippocampus",
        "Left-Amygdala","l amygdala",
        "Right-Putamen","r putamen",
        "Right-Hippocampus","r hippocampus",
        "Right-Amygdala", "r amygdala")
return(orig2plot_names)
}
      
wrapper_plot_pca = function(wd, outdir) {
  outdir.mega.dim = file.path(outdir, "dimensionality")
  load(file.path(wd,"mega","PCA" ,"pca_result.rda"))
  
  df.importance =  summary(pca_result)$importance  %>% t() %>% data.frame()
  df.importance$PC = gsub("PC", "", rownames(df.importance)) %>% 
    as.numeric() %>% as.factor()
  df.loading = data.frame(features = rownames(pca_result$rotation), 
                          loading = -pca_result$rotation[,1])
  
  
  
  df.changenames = retrieve_orig2plot_names()
  df.loading$plot.names = df.changenames$plot.names[which(df.changenames$orig.names %in% rownames(pca_result$rotation) )]
  df.loading = 
    df.loading %>% arrange(desc(loading)) %>% 
    mutate(id = 1:nrow(.) %>% as.factor())
  
  
  gs1 = 
    ggplot(df.importance, aes(x = PC, y = Proportion.of.Variance)) +
    geom_bar(stat = "identity", width = 0.8, color = "black", fill = "#96C7C1", size = 0.8) +  # Bars with outlines
    #geom_text(aes(label = PC), vjust = -0.5, size = 5, fontface = "bold") +  # Add labels above bars
    labs(
      x = "Principal component",
      y = "Proportion of Variance"
    ) +
    theme_classic() +  # Minimal theme with custom base size
    theme(
      text=element_text(family="Arial", face = "bold"),
      plot.title = element_blank(),  # Centered title, 
      axis.line.x = element_blank(),
      axis.text.x = element_text(face = "bold", size = 14,vjust = 5.5),
      axis.text.y = element_text(face = "bold", size = 14),
      axis.title.x = element_text(face = "bold", size = 16, vjust = 4.5),
      axis.title.y = element_text(face = "bold", size = 16),
      axis.ticks.x = element_blank(),
      legend.position = "none" )
  
  
  
  gs2 = 
    ggplot(df.loading, aes(x = id, y = loading)) +
    geom_bar(stat = "identity", width = 0.8, color = "black", fill = "#96C7C1", size = 0.8) +  # Bars with outlines
    #geom_text(aes(label = PC), vjust = -0.5, size = 5, fontface = "bold") +  # Add labels above bars
    labs(
      x = "Region",
      y = "Loadings"
    ) +
    #scale_x_discrete(labels = df.loading$plot.names,expand = expansion(mult = c(0.05, 0))) + 
    theme_classic() +  # Minimal theme with custom base size
    theme(
      text=element_text(family="Arial", face = "bold"),
      plot.title = element_blank(),  # Centered title, 
      axis.line.x = element_blank(),
      axis.text.x = element_text(vjust = 5.5),
      axis.ticks.x = element_blank(), 
      axis.text = element_text(face = "bold", size = 14),
      axis.title.x = element_text(face = "bold", size = 16, vjust = 4.5),
      axis.title.y = element_text(face = "bold", size = 16),
      #axis.ticks.x = element_blank(),
      legend.position = "none" )
  
  ggsave(file.path(outdir.mega.dim, "pca.variance.png"), 
         plot = gs1)  
  ggsave(file.path(outdir.mega.dim, "pca.loading.jpeg"), 
         plot = gs2)  
  
  gs = list(
    gs1 = gs1,
    gs2 = gs2
  )
  return(gs)
}



wrapper_hemi_plot_main = function(wd, outdir) {
  
  library(see)
  library(gghalves)
  
  outdir.mega.main = file.path(outdir, "main")
  load(file.path(wd,"mega","hemispheric_diff", "hemi.rda"))
  
  # CONVERT TO LONG
  data_long = 
    hemi.diff$dat2 %>% 
    select(feature2, gamm.w.beta.m5.0.weighted_rh, gamm.w.beta.m5.0.weighted_lh) %>% 
    pivot_longer(-feature2, 
                 names_to = "hemi", 
                 values_to = "values") %>% 
    mutate(hemi = if_else(hemi == "gamm.w.beta.m5.0.weighted_rh", "right", "left") %>% 
             as.factor)
  
  
  data_long_nudged <- data_long %>%
    mutate(adjusted_x = ifelse(hemi == "left", 1.2, 
                               1.8)) %>% 
    rowwise() %>% 
    mutate(adjusted_x = adjusted_x + rnorm(1, sd = .04)) 
  
  
  
  gs = ggplot(data_long, aes(x = hemi, y = values, fill = hemi)) +
    geom_half_violin(
      side = c("l", "r"), 
      nudge = .12, 
      trim = TRUE, 
      width = .3,
      size = 1
    )  +
    geom_half_boxplot(
      data = data_long %>% filter(hemi == "left"), 
      width = .01,
      nudge = -0.09,
      side = "r",
      size = 1,
      outlier.alpha = 0
    ) +
    geom_half_boxplot(
      data = data_long %>% filter(hemi == "right"), 
      width = .01,
      nudge = -0.09,
      side = "l",
      size = 1,
      outlier.alpha = 0
    ) +
    # Add points from the nudged data
    geom_point(data = data_long_nudged, 
               aes(x = adjusted_x, y = values, fill = hemi), color = "black", shape = 21, alpha = .5) +
    # Add lines based on the original data (no nudge)
    geom_line(data = data_long_nudged,
              aes(x = adjusted_x,group = feature2), 
              color = "grey70", 
              size = 0.7, 
              linetype = "dashed", 
              lineend = "butt") +
    scale_color_manual(values = c("#FDB4C7", "#B5EAD7")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 1) +
    scale_fill_manual(values = c("#FDB4C7", "#B5EAD7")) +
    scale_x_discrete(name = "hemisphere") + #
    scale_y_continuous(name =expression(beta * "eta")) +
    theme_minimal() +
    theme(
      #aspect.ratio = 1.25,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 1),
      axis.text.x = element_text(size = 18, family = "Arial", color = "black"),
      axis.text.y = element_text(size = 18, family = "Arial", color = "black"),
      axis.title= element_text(size = 20, family = "Arial", color = "black"),
      #axis.title.y = element_blank(),
      plot.title = element_text(size = 24, face = "bold", family = "Arial", hjust = 0.5),
      legend.position = "none"  # Hide legend in individual plot
    ) 
  ggsave(file.path(outdir.mega.main, "hemi_paired.png"), 
         plot = gs)  
  
  
  df.diff = 
    hemi.diff$dat2 %>% 
    mutate(diff = gamm.w.beta.m5.0.weighted_lh - gamm.w.beta.m5.0.weighted_rh,
           g = "g", 
           adjusted_x = 1.2) %>% 
    rowwise() %>% 
    mutate(adjusted_x = adjusted_x + rnorm(1, sd = .04)) 
  
  
  gs1 = ggplot(df.diff, aes(x = g, y = diff, fill = g)) +
    geom_half_violin(
      side = c("l"), 
      nudge = .12, 
      trim = TRUE, 
      width = .3,
      size = 1
    )  +
    geom_half_boxplot(
      width = .01,
      nudge = -0.09,
      side = "r",
      size = 1,
      outlier.alpha = 0
    ) +
    geom_hline(
      yintercept = 0, 
      linetype = "dashed", 
      color = "grey50"
    ) +
    # Add points from the nudged data
    geom_point(mapping =  
                 aes(x = adjusted_x, y = diff, fill = g), color = "black",shape = 21, alpha = .5) + 
    scale_fill_manual(values = "#A7D8F0")+
    theme_minimal() +
    theme(
      #aspect.ratio = 1.25,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 1),
      #axis.text.x =  element_blank(),
      axis.text = element_text(size = 18, family = "Arial", color = "black"),
      axis.title = element_text(size = 18, family = "Arial", color = "black"),
      plot.title = element_text(size = 24, face = "bold", family = "Arial", hjust = 0.5),
      legend.position = "none",  # Hide legend in individual plot
      plot.margin = margin(5, 5, 5, 5)) + 
    scale_x_discrete(name = "difference", expand = c(.2, .2),labels = c(expression(Delta))) + 
    scale_y_continuous(position = "right", expression(beta * "eta (left - right)"))
  
  ggsave(file.path(outdir.mega.main, "hemi_diff.png"), 
         plot = gs1)  
  gsall = gs + gs1 + plot_layout(widths = c(3, 1.3))
  ggsave(file.path(outdir.mega.main, "hemi_combined.png"), 
         plot = gsall, width = 10, height = 7)  
  
  out = list()
  out$gs = gs
  out$gs1 = gs1
  out$gsall = gsall
  return(out)
}



wrapper_hemi_plot_int = function(wd, outdir) {
  outdir.mega.main = file.path(outdir, "main")
  load(file.path(wd,"mega","hemispheric_diff", "hemi.rda"))
  hemi.diff$x = hemi.diff$x %>% 
    mutate(lower_CI = map_dbl(ttest.age50, ~.$conf.int[[1]]), 
           upper_CI = map_dbl(ttest.age50, ~.$conf.int[[2]]))
  
  df.plot = hemi.diff$x %>% filter(xage %in% c(40,50,60,70,80)) %>% 
    mutate(age = as.factor(xage))
  
  (p.adjust(df.plot$t.test.age50.p.value, method = "fdr"))
  
  gs.errorbar = 
    ggplot(df.plot, aes(age, t.test.age50.mean_difference, color = age, group = age, fill = age)) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") + 
    geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = .3, size = 1.5) +
    geom_point(shape = 21, size = 5, fill = "ivory", stroke = 1.5) + 
    scale_fill_manual(values = c("#FFECB3", "#FDCFE8", "#B5D0EA","#C2F2E5", "#FFDAC1")) + 
    scale_color_manual(values = c("#FFECB3", "#FDCFE8", "#B5D0EA","#C2F2E5", "#FFDAC1")) +
    scale_x_discrete(name = "Age") + 
    scale_y_continuous(name = expression(beta * "eta (left - right)")) + 
    theme_minimal() +
    theme(
      #aspect.ratio = 1.25,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 1),
      #axis.text.x =  element_blank(),
      axis.text = element_text(size = 18, family = "Arial", color = "black"),
      axis.title = element_text(size = 18, family = "Arial", color = "black"),
      plot.title = element_text(size = 24, face = "bold", family = "Arial", hjust = 0.5),
      legend.position = "none") +
    annotate("text", x = 2, y = .032, label = "†", size = 10, color = "grey50") +
    annotate("text", x = 1, y = .034, label = "*", size = 12, color = "grey50")
  
  ggsave(file.path(outdir.mega.main, "hemi_int.png"), 
         plot = gs.errorbar)  
  return(gs.errorbar)
}



wrapper_memory_plot_main = function(wd, outdir) {
  library(see)
  library(gghalves)
  outdir.mega.mem = file.path(outdir, "memory_apoe")
  
  load(file.path(wd,"mega/apoe/output/temp_models_mem",
                 paste0("models.Rda")))
  
  
  res = predict(models$mod.lme.w.apoe.memory, re.form = ~0) + residuals(models$mod.lme.w.apoe.memory)
  df = models$mod.lme.w.apoe.memory@frame %>% 
    mutate(apoe_statusF = factor(apoe_status, levels = c(0,1), labels = c("non-carrier", "carrier")))
  df$res = res
  
  gs = 
    ggplot(df, aes(x = apoe_statusF, y = res, fill = apoe_statusF)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")+
    geom_half_violin(
      side = c("l", "r"), 
      nudge = -.05, 
      trim = TRUE, 
      width = .45,
      size = 1
    ) +
    geom_boxplot(width = 0.25, position = position_nudge(x = -0.25),, outlier.shape = NA, 
                 data = df %>% filter(apoe_statusF == "carrier")) +
    geom_boxplot(width = 0.25, position = position_nudge(x = 0.25),, outlier.shape = NA, 
                 data = df %>% filter(apoe_statusF == "non-carrier")) +
    theme_minimal() +
    scale_color_manual(values = c("#B5EAD7", "#FDB4C7")) +
    scale_fill_manual(values = c( "#B5EAD7", "#FDB4C7")) +
    scale_x_discrete(name = expression("APOE " * epsilon * "4")) + #
    scale_y_continuous(name =expression(Delta * "memory")) +
    theme_minimal() +
    theme(
      #aspect.ratio = 1.25,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 1),
      axis.text.x = element_text(size = 18, family = "Arial", color = "black"),
      axis.text.y = element_text(size = 18, family = "Arial", color = "black"),
      axis.title= element_text(size = 20, family = "Arial", color = "black"),
      #axis.title.y = element_blank(),
      plot.title = element_text(size = 24, face = "bold", family = "Arial", hjust = 0.5),
      legend.position = "none"  # Hide legend in individual plot
    ) 
  
  ggsave(file.path(outdir.mega.mem, "maineffect.png"), 
         plot = gs)  
  
  return(gs)
}


wrapper_memory_plot_age = function(wd, outdir) {
  outdir.mega.mem = file.path(outdir, "memory_apoe")
  load(file.path(wd,"mega/apoe/output/temp_models_mem",
                 paste0("models.Rda")))
  
  df.plot = predictions[[2]] %>% 
    mutate(apoe_statusF = factor(apoe_statusO, levels = c(0,1), labels = c("non-carrier", "carrier"))) %>% 
    filter(between(age, 35,85)) # rounding 2.5%  and 97.5%
  
  gs = 
    ggplot(df.plot, aes(age,estimate, group = apoe_statusF, fill = apoe_statusF)) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") + 
    geom_ribbon(aes(ymin = estimate - 1.96*se, ymax = estimate + 1.96*se), alpha = .6)+ 
    geom_line(size = 1.5) +
    theme_classic() +
    scale_fill_manual( values = c( "#B5EAD7", "#FDB4C7")) +
    theme(#legend.position = "right", 
          #legend.box.margin = margin(30, 0, 0, 0), 
          legend.position= c(0.6, .9), 
          axis.line = element_line(color = "black", size = 1),
          legend.title = element_blank(), 
          axis.text.x = element_text(size = 18, family = "Arial", color = "black"),
          axis.text.y = element_text(size = 18, family = "Arial", color = "black"),
          axis.title= element_text(size = 20, family = "Arial", color = "black"),
          strip.text = element_text(size = 16, face = "bold"),
          axis.text = element_text(size = 16), 
          legend.text = element_text(size = 18, family = "Arial"), 
          text=element_text(family="Arial")) + #
    scale_y_continuous(name =expression(Delta * "memory")) +
    labs(y = expression(Delta~"Brain"))
  ggsave(file.path(outdir.mega.mem, "ageinteraction.png"), 
         plot = gs)  
  
  return(gs)
}

wrapper_plot_apoe_brain_main = function(wd, outdir) {
  library(see)
  library(gghalves)
  outdir.mega.brain= file.path(outdir, "brain_apoe")
  load( file.path(wd,"mega","apoe", "df.mega.apoe.brain.rda"))
  
  
  f = c("Left-Hippocampus", "Right-Hippocampus", "Right-Amygdala")
  grot = list()
  for (i in f) {
    load(file.path(wd,"mega/apoe/output/temp_models_mem",
                   paste0("models.",i,".Rda")))
    df = models$mod.lme.w.apoe.brain@frame %>% 
      mutate(apoe_statusF = factor(apoe_status, levels = c(0,1), labels = c("non-carrier", "carrier")))
    df$res = predict(models$mod.lme.w.apoe.brain, re.form = ~0) + residuals(models$mod.lme.w.apoe.brain)
    df$features = i
    grot[[i]] = df %>% 
      select(apoe_statusF, 
             res, 
             features)
    
  }
  
  #quantile(models$mod.gam.w.apoe.brain$gam$model$age, probs = c(.025,.975))
  df.plot = data.table::rbindlist(grot)
  df.plot = 
    df.plot %>%  
    mutate(features = factor(features, 
                             levels = c("Left-Hippocampus",  "Right-Hippocampus", "Right-Amygdala"), 
                             labels = c("Left Hippocampus",  "Right Hippocampus", "Right Amygdala")))
  
  gs = 
    ggplot(df.plot, aes(x = apoe_statusF, y = res, fill = apoe_statusF)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", size = 1.5)+
    geom_half_violin(
      side = c("l", "r"), 
      nudge = -.10, 
      trim = TRUE, 
      width = .45,
      size = 1
    ) +
    geom_boxplot(width = 0.25, position = position_nudge(x = -0.3),, outlier.shape = NA, 
                 data = df.plot %>% filter(apoe_statusF == "carrier")) +
    geom_boxplot(width = 0.25, position = position_nudge(x = 0.3),, outlier.shape = NA, 
                 data = df.plot %>% filter(apoe_statusF == "non-carrier")) +
    theme_minimal() +
    scale_color_manual(values = c("#B5EAD7", "#FDB4C7")) +
    scale_fill_manual(values = c("#B5EAD7", "#FDB4C7")) +
    scale_x_discrete(name = expression("APOE " * epsilon * "4"), expand = c(0.3,0)) + #
    scale_y_continuous(name =expression(Delta * "brain")) +
    theme_minimal() +
    theme(
      #aspect.ratio = 1.25,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 1),
      axis.text.x = element_text(size = 18, family = "Arial", color = "black"),
      axis.text.y = element_text(size = 18, family = "Arial", color = "black"),
      axis.title= element_text(size = 20, family = "Arial", color = "black"),
      strip.text = element_text(size = 20, family = "Arial", color = "black"), 
      plot.title = element_text(size = 24, face = "bold", family = "Arial", hjust = 0.5),
      legend.position = "none") +
    facet_wrap(vars(features), ncol = 3)
  
  ggsave(file.path(outdir.mega.brain, "maineffect.png"), 
         plot = gs, 
         width = 20, 
         heigh = 7)  
  return(gs)
}


wrapper_plot_apoe_brain_main_across= function(wd, outdir) {
  outdir.mega.brain= file.path(outdir, "brain_apoe")
  load( file.path(wd,"mega","apoe", "df.mega.apoe.brain.rda"))  
  #mega.apoe.brain$lme.apoe.brain_beta
  df.plot1 = 
    mega.apoe.brain %>% 
    mutate(g = "g", 
           adjusted_x = 1.2) %>% 
    rowwise() %>% 
    mutate(adjusted_x = adjusted_x + rnorm(1, sd = .04)) 
  
  
  gs1 = ggplot(df.plot1, aes(x = g, y = lme.apoe.brain_beta, fill = g)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", size = 1.5)+
    geom_half_violin(
      side = c("l"), 
      nudge = .12, 
      trim = TRUE, 
      width = .3,
      size = 1
    )  +
    geom_half_boxplot(
      width = .01,
      nudge = -0.09,
      side = "r",
      size = 1,
      outlier.alpha = 0
    ) +
  geom_point(mapping =  
                 aes(x = adjusted_x, y = lme.apoe.brain_beta, fill = g), color = "black",shape = 21, alpha = .5) + 
    scale_fill_manual(values = "#A7D8F0")+
    theme_minimal() +
    theme(
      #aspect.ratio = 1.25,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 1),
      #axis.text.x =  element_blank(),
      axis.text = element_text(size = 18, family = "Arial", color = "black"),
      axis.title = element_text(size = 18, family = "Arial", color = "black"),
      plot.title = element_text(size = 20, family = "Arial", hjust = 0.5),
      legend.position = "none", 
      plot.subtitle = element_text(hjust = 0.5, size = 20, family = "Arial", vjust = -6)) + 
    scale_x_discrete(name = expression("APOE " * epsilon * "4"),expand = c(0.2,0.15), labels = expression(beta)) + #labels = c(expression(epsilon * "4 effect"))) + 
    scale_y_continuous(position = "right", expression(Delta * "brain")) +
    #ggtitle("Across regions") +
    labs(subtitle = "Across regions")
  
  ggsave(file.path(outdir.mega.brain, "allregions.png"), 
         plot = gs1,
         width = 6, height = 7)  
  
  return(gs1)
}


#######PLot
wrapper_plot_apoe_brain_age = function(wd, outdir) {
  f = c("Left-Hippocampus", "Right-Hippocampus")
  grot = list()
  grot2 = list()
  for (i in f) {
    load(file.path(wd,"mega/apoe/output/temp_models_mem",
                   paste0("models.",i,".Rda")))
    
    #df.predict.tensory.interaction by apoe_group
    grid <- expand.grid(age = seq(25, 85, by = .1), 
                        apoe_statusO = c(0,1))
    pred.apoe <- predict(models$mod.gam.w.apoe.brain$gam, newdata = grid, se.fit = TRUE)
    
    sm1 = 
      data.frame(
        age = grid$age, 
        apoe_statusO = grid$apoe_statusO, 
        estimate = pred.apoe$fit, 
        se = pred.apoe$se.fit)
    sm1 = 
      sm1 %>% 
      group_by(apoe_statusO) %>% 
      nest() %>% 
      mutate(data = map(data, ~ .x %>%
                          mutate(derivative = c(NA, diff(.$estimate) / diff(.$age))))) %>% 
      unnest() 
    
    
    grot[[i]] = sm1
    grot2[[i]] = smooth_estimates(models$mod.gam.w.apoe.brain) |>
      add_confint() %>% filter(!is.na(apoe_statusO)) %>% 
      mutate(.estimate = .estimate + summary(models$mod.gam.w.apoe.brain$gam)$p.coeff[[2]])
  }
  df.age = data.table::rbindlist(grot, idcol = "features") %>% 
    mutate(apoe_statusF = factor(apoe_statusO, levels = c(0,1), labels = c("non-carrier", "carrier"))) %>% 
    filter(between(age, 38.5,84))  %>% 
    mutate(features = factor(features, 
                             levels = c("Left-Hippocampus", 
                                        "Right-Hippocampus"),
                             labels = c("Left Hippocampus", 
                                        "Right Hippocampus")))
  # rounding the first/3rd quarticle quartile
  df.term= data.table::rbindlist(grot2, idcol = "features") %>% 
    mutate(apoe_statusF = factor(apoe_statusO, levels = c(0,1), labels = c("non-carrier", "carrier"))) %>% 
    filter(between(age, 38.5,84))  %>% 
    mutate(features = factor(features, 
                             levels = c("Left-Hippocampus", 
                                        "Right-Hippocampus"),
                             labels = c("Left Hippocampus", 
                                        "Right Hippocampus"))) # rounding the first/3rd quarticle quartile)
  
  
  gs = 
    ggplot(df.age, aes(age,estimate, group = apoe_statusF, fill = apoe_statusF)) + 
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey50", size = 1.5) + 
    geom_ribbon(aes(ymin = estimate - se*1.96, ymax = estimate + se*1.96), alpha = .6)+ 
    geom_line(size = 1.5) +
    facet_grid(cols = vars(features), scales = "free_y", switch = "y") + 
    theme_classic() +
    scale_fill_manual(values =   c("#B5EAD7", "#FDB4C7")) +
    theme(strip.placement = "outside", 
          legend.position = "right", 
          legend.title = element_blank(), 
          axis.title = element_text(face = "bold", size = 18, family = "Arial"),
          strip.background = element_blank(), 
          strip.text.x =element_text(size = 18, face = "bold", family = "Arial"),
          strip.text = element_text(size = 18, face = "bold", family = "Arial"),
          axis.text = element_text(size = 18, face = "bold", family = "Arial"),
          legend.text = element_text(size = 18, face = "bold", family = "Arial"), 
          #axis.line.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(), 
          text=element_text(family="Arial")) +
    labs(y = expression(Delta~"Brain"))
  
  
  gs1 = 
    ggplot(df.term, aes(age,.estimate), color = "black", fill = "lightblue") + 
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey50", size = 1.5) + 
    geom_ribbon(aes(ymin = .estimate - .se*1.96, ymax = .estimate + .se*1.96), alpha = .6, fill = "lightblue")+ 
    geom_line(size = 1.5) +
    facet_grid(cols = vars(features), scales = "free_y", switch = "y") + 
    theme_classic() +
    theme(strip.placement = "none", 
          legend.position = "bottom", 
          legend.title = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_blank(), 
          axis.title = element_text(face = "bold", size = 18, family = "Arial"),
          #axis.title.x = element_blank(),
          #strip.text = element_text(size = 18, face = "bold", family = "Arial"),
          axis.text = element_text(size = 18, face = "bold", family = "Arial"),
          legend.text = element_text(size = 18, face = "bold", family = "Arial"), 
          text=element_text(family="Arial")) +
    labs(y = expression(Delta~"Brain/year"))
  
  
  gsall =  gs / gs1 + 
    plot_layout(heights = c(2, 1)) 
  
  ggsave(file.path(outdir.mega.brain, "ageintraction.png"), 
         plot = gsall,
         width = 10, height = 7)  
  
  return(gsall)
}


wrapper_plot_apoe_main = function(wd, outdir) {
  outdir.mega.apoe= file.path(outdir, "apoe")
  
  f = c("Left-Hippocampus", "Right-Hippocampus")
  grot = list()
  grot2 = list()
  for (i in f) {
    load(file.path(wd,"mega/apoe/output/temp_models",
                   paste0("models_",i,".Rda")))
    
    
    grid <- expand.grid(apoe_statusO = c(0,1), 
                        delta_brain = seq(-4.5,4.5, length.out = 100))
    predict.apoeO = predict(models$mod.gam_main.w.apoeO$gam, grid, se.fit = T)
    df.1 = data.frame(apoe_statusO = grid$apoe_statusO, 
                      delta_brain = grid$delta_brain, 
                      estimate = predict.apoeO$fit, 
                      se = predict.apoeO$fit)
    
    df.1 = df.1  %>% 
      mutate(apoe_statusF = factor(apoe_statusO, levels = c(0,1), labels = c("non-carrier", "carrier")), 
             features = i) 
    
    grot[[i]] = df.1
    grot2[[i]] = predictions[[5]] %>% filter(!is.na(apoe_statusO))
  }
  
  df.1 = data.table::rbindlist(grot) %>% 
    mutate(features = factor(features, 
                             levels = c("Left-Hippocampus", 
                                        "Right-Hippocampus"),
                             labels = c("Left Hippocampus", 
                                        "Right Hippocampus")))
  # rounding the first/3rd quarticle quartile
  df.term= data.table::rbindlist(grot2) %>% 
    mutate(features = factor(features, 
                             levels = c("Left-Hippocampus", 
                                        "Right-Hippocampus"),
                             labels = c("Left Hippocampus", 
                                        "Right Hippocampus"))) # rounding the first/3rd quarticle quartile)
  
  df.1 = df.1 %>% filter(between(delta_brain, -2.5, 2.5))
  gs = 
    ggplot(df.1, aes(delta_brain,estimate, group = apoe_statusF, fill = apoe_statusF)) + 
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey50", size = 1.5) + 
    geom_ribbon(aes(ymin = estimate - se*1.96, ymax = estimate + se*1.96), alpha = .6)+ 
    geom_line(aes(color = apoe_statusF), size = 1.5) +
    facet_grid(rows =  vars(features), scales = "free_y", switch = "y") + 
    theme_classic() +
    scale_color_manual(values =   c("#7CAE00", "#F8766D")) +
    scale_fill_manual(values =   c( "#B5EAD7", "#FDB4C7")) +
    theme(strip.placement = "outside", 
          legend.position = "bottom", 
          legend.title = element_blank(), 
          axis.title = element_text(face = "bold", size = 18, family = "Arial"),
          #axis.title = element_blank(),
          strip.background = element_blank(),
          strip.text.x =element_text(size = 18, face = "bold", family = "Arial"),
          #strip.text = element_text(size = 18, face = "bold", family = "Arial"),
          strip.text = element_blank(),
          axis.text = element_text(size = 18, face = "bold", family = "Arial"),
          legend.text = element_text(size = 18, face = "bold", family = "Arial"), 
          text=element_text(family="Arial")) + 
    labs(y = expression(Delta~"Memory"), x = expression(Delta~"Brain"))
  
  
  ggsave(file.path(outdir.mega.apoe, "main_hipp.png"), 
         plot = gs, 
         height = 9, 
         width = 7)  
  return(gs)
}


wrapper_plot_apoe_int = function(wd, outdir) {
  outdir.mega.apoe= file.path(outdir, "apoe")
  f = c("Left-Hippocampus", "Right-Hippocampus")
  grot = list()
  grot2 = list()
  for (i in f) {
    load(file.path(wd,"mega/apoe/output/temp_models",
                   paste0("models_",i,".Rda")))
    
    
    grid2 <- expand.grid(apoe_statusO = c(0,1), 
                         delta_brain = seq(-4.5,4.5, length.out = 100), 
                         xage = seq(40, 80, by = 10))
    predict.apoeO.int = predict(models$mod.gam_int.w.apoeO$gam, grid2, se.fit = T)
    df.2 = data.frame(apoe_statusO = grid2$apoe_statusO, 
                      delta_brain = grid2$delta_brain, 
                      age = grid2$xage,
                      estimate = predict.apoeO.int$fit, 
                      se = predict.apoeO.int$fit)
    
    df.2 = df.2  %>% 
      mutate(apoe_statusF = factor(apoe_statusO, levels = c(0,1), labels = c("non-carrier", "carrier")), 
             features = i) 
    
    grot[[i]] = df.2
    grot2[[i]] = predictions[[2]] %>% filter(!is.na(apoe_statusO))
  }
  
  df.2 = data.table::rbindlist(grot) %>% 
    mutate(features = factor(features, 
                             levels = c("Left-Hippocampus", 
                                        "Right-Hippocampus"),
                             labels = c("Left Hippocampus", 
                                        "Right Hippocampus"))) %>% 
    mutate(f_label = if_else(features == "Left Hippocampus", "left\nhippocampus","right\nhippocampus"))
  # rounding the first/3rd quarticle quartile
  df.term= data.table::rbindlist(grot2) %>% 
    mutate(features = factor(features, 
                             levels = c("Left-Hippocampus", 
                                        "Right-Hippocampus"),
                             labels = c("Left Hippocampus", 
                                        "Right Hippocampus"))) # rounding the first/3rd quarticle quartile)
  
  df.2 = df.2 %>% filter(between(delta_brain, -2.5, 2.5))
  gs1 = 
    ggplot(df.2, aes(delta_brain,estimate, group = apoe_statusF, fill = apoe_statusF)) + 
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey50", size = 1.5) + 
    geom_ribbon(aes(ymin = estimate - se*1.96, ymax = estimate + se*1.96), alpha = .6)+ 
    geom_line(aes(color = apoe_statusF), size = 1.5) +
    #facet_grid(rows =  vars(features), scales = "free_y") + 
    theme_classic() +
    scale_color_manual(values =   c("#7CAE00", "#F8766D")) +
    scale_fill_manual(values =   c("#B5EAD7", "#FDB4C7")) +
    theme(strip.placement = "outside", 
          legend.position = "bottom", 
          legend.title = element_blank(), 
          axis.title = element_text(face = "bold", size = 18, family = "Arial"),
          strip.background = element_blank(), 
          strip.text.x =element_text(size = 18, face = "bold", family = "Arial"),
          strip.text = element_text(size = 18, face = "bold", family = "Arial"),
          strip.text.y = element_text(face = "bold", size = 18, angle = 90, family = "Arial"), 
          axis.title.y = element_blank(), 
          axis.text = element_text(size = 18, face = "bold", family = "Arial"),
          legend.text = element_text(size = 18, face = "bold", family = "Arial"), 
          text=element_text(family="Arial")) + 
    facet_grid(cols = vars(age), rows = vars(features), scales = "free_y") + 
    labs(y = expression(Delta~"Memory"), x = expression(Delta~"Brain"))  
  
  
  ggsave(file.path(outdir.mega.apoe, "int_hipp.png"), 
         plot = gs1, 
         height = 9, 
         width = 18)  
  return(gs1) 
}

kable_2_table = function(tb, bold = "last", digit = 2) {
  kb = 
    tb  %>% 
    kable(digits = digit) %>% 
    kable_styling(bootstrap_options = "striped", full_width = F)
  
  if (bold == "last") {
    kb = 
      kb%>% 
      row_spec(nrow(tb), bold = TRUE)
  }
  return(kb)
}

retrieve_pboot_cc = function(wd, type ="basic") {
  
  wd = here("data_memory_long/all")
  
  if (!type == "all") {
    if (type == "pca") {
      load(file.path(wd,"mega","consensus_cluster", "pca_consensus_cluster_gam.Rda"))
      cl = paste0("cl", c(1:8) %>% as.character())
    } else if (type == "cl1") {
      load(file.path(wd,"mega","consensus_cluster", "cl1_consensus_cluster_gam.Rda"))
      cl = paste0("cl", c(2:8) %>% as.character())
    } else if (type == "basic") {
      load(file.path(wd,"mega","consensus_cluster", "basic_consensus_cluster_gam.Rda"))
      cl = paste0("cl", c(1:8) %>% as.character())
    }
    
    trues = ls(pattern = "mod.")
    
    pvalues = lapply(trues, function(x) {summary(get(x)$gam)$s.table[1,4]}) %>% 
      simplify2array()
    
    df = 
      data.frame(cl = cl, 
                 pval = pvalues)
    
    df = 
      bootstraps %>%
      filter(!is.na(.))%>% 
      group_by(cl) %>% 
      nest() %>% 
      left_join(df,.)
    
    df = 
      df %>% 
      mutate(pboot = 
               map2_dbl(pval, data, ~sum(.x > .y) / length(.y[[1]]))) %>% 
      select(-data)
  }
  
  return(df)
}


wrapper_simulation_nonlinear = function(outdir.mega.simulation) {
  load(file.path(outdir.mega.simulation, "modelling.rda"))
  ## Select analysis
  df = 
    mini.grid %>% 
    filter(between(x1,-.9, -.7), # -.8 
           sd1 == 1.2,
           skew1 == -5, 
           sd2 == .8)
  
  
  
  df.model = rbind(df$smt.model[[1]][[2]], df$smt.aging[[1]][[2]])
  out.dataset.x = 
    df$out.dataset[[1]] %>% 
    mutate(g = "g") %>% 
    pivot_longer(-g, 
                 names_to = "G", 
                 values_to = "X") %>% 
    filter(G %in% c("X", "X1", "X2"))
  
  out.dataset.y = 
    df$out.dataset[[1]] %>% 
    mutate(g = "g") %>% 
    pivot_longer(-g, 
                 names_to = "G", 
                 values_to = "mean") %>% 
    filter(G %in% c("Y", "Ynoise"))
  
  
  
  gs = ggplot(df.model, aes(X, mean, group = G, fill = G, color = G)) + 
    geom_hline(yintercept = 0, color = "grey50", linetype = "dashed")+ 
    geom_vline(xintercept = 0, color = "grey50", linetype = "dashed")+ 
    geom_line(linewidth = 2) +
    geom_ribbon(aes(ymin = li, ymax = ui), alpha = .6, linewidth = 0)+ 
    theme_classic() +
    scale_fill_manual(values =   c("#FDB4C7", "#B5EAD7", "#A7D8F0","#FFECB3","#FFDAC1", "#A7D8F0", "#FFDAC1")) +
    scale_color_manual(values =   c("#FDB4C7", "#B5EAD7", "#A7D8F0","#FFECB3","#FFDAC1", "#A7D8F0", "#FFDAC1")) +
    theme(legend.position = 'none', 
          legend.title = element_blank(), 
          axis.title = element_text(face = "bold", size = 18),
          axis.text = element_text(size = 18), 
          #legend.text = element_text(size = 14, face = "bold"), 
          text=element_text(family="Arial"), face = "bold") +
    labs(x = expression(Delta~"Brain"),
         y = expression(Delta~"Memory"))
  
  gs1 =   
    gs + 
    geom_xsidedensity(data = out.dataset.x, aes(X, group = G, color = G), alpha = .6) +
    geom_ysidedensity(data = out.dataset.y, aes(y =mean, group = G, color = G), alpha = .6) +
    ylim(c(-1.1,.8)) + 
    xlim(c(c(-4.5,1.5))) + 
    scale_ysidex_continuous(minor_breaks = NULL, breaks = NULL) + 
    scale_xsidey_continuous(minor_breaks = NULL, breaks = NULL)
  
  ggsave(file.path(outdir.mega.simulation, "main.simulation.png"), 
         plot = gs1, 
         width = 6, 
         height = 8)   
  
  gs2 = gs1 + theme(legend.position = "bottom", 
                    legend.text = element_text(size = 16), 
                    legend.title = element_text(size = 16))
  ggsave(file.path(outdir.mega.simulation, "main.simulation.legend.png"), 
         plot = gs2, 
         width = 6, 
         height = 8)   
  
  out = list(gs1, 
             gs2, 
             df) 
  return(out)
}

wrapper_simulation_meanchange = function(outdir.mega.simulation){
  
  load(file.path(outdir.mega.simulation, "modelling.rda"))
  ## Select analysis
  df = 
    mini.grid %>% 
    filter(sd1 == 1.2,
           skew1 == -5, 
           sd2 == .8)
  
  ###
  df.model = rbind(df$smt.model[[1]][[2]] %>% mutate(M ="low"),
                   df$smt.model[[2]][[2]] %>% mutate(M ="mid"),
                   df$smt.model[[3]][[2]] %>% mutate(M ="high"))
  
  out.dataset.x = 
    rbind(
      df$out.dataset[[1]] %>% 
        mutate(g = "g") %>% 
        pivot_longer(-g, 
                     names_to = "G", 
                     values_to = "X") %>% 
        filter(G %in% c("X1")) %>% 
        mutate(M = "low"), 
      df$out.dataset[[2]] %>% 
        mutate(g = "g") %>% 
        pivot_longer(-g, 
                     names_to = "G", 
                     values_to = "X") %>% 
        filter(G %in% c("X1")) %>% 
        mutate(M = "mid"), 
      df$out.dataset[[3]] %>% 
        mutate(g = "g") %>% 
        pivot_longer(-g, 
                     names_to = "G", 
                     values_to = "X") %>% 
        filter(G %in% c("X1")) %>% 
        mutate(M = "high")
    )
  
  
  out.dataset.y = 
    rbind(
      df$out.dataset[[1]] %>% 
        mutate(g = "g") %>% 
        pivot_longer(-g, 
                     names_to = "G", 
                     values_to = "X") %>% 
        filter(G %in% c("Y")) %>% 
        mutate(M = "low"), 
      df$out.dataset[[2]] %>% 
        mutate(g = "g") %>% 
        pivot_longer(-g, 
                     names_to = "G", 
                     values_to = "X") %>% 
        filter(G %in% c("Y")) %>% 
        mutate(M = "mid"), 
      df$out.dataset[[3]] %>% 
        mutate(g = "g") %>% 
        pivot_longer(-g, 
                     names_to = "G", 
                     values_to = "X") %>% 
        filter(G %in% c("Y")) %>% 
        mutate(M = "high")
    ) %>% 
    rename(mean = X)
  
  
  df.model = 
    df.model %>% 
    mutate(M = factor(M, levels = c("low", "mid", "high")))
  
  out.dataset.y = 
    out.dataset.y %>% 
    mutate(M = factor(M, levels = c("low", "mid", "high")))
  
  out.dataset.x = 
    out.dataset.x %>% 
    mutate(M = factor(M, levels = c("low", "mid", "high")))
  
  gs = ggplot(df.model, aes(X, mean, group = M, fill = M, color = M)) + 
    geom_hline(yintercept = 0, color = "grey50", linetype = "dashed")+ 
    geom_vline(xintercept = 0, color = "grey50", linetype = "dashed")+ 
    geom_line(linewidth = 2) +
    geom_ribbon(aes(ymin = li, ymax = ui), alpha = .6, linewidth = 0)+ 
    theme_classic() +
    scale_fill_manual(values =   c("#FDB4C7", "#B5EAD7", "#A7D8F0")) +
    scale_color_manual(values =   c("#FDB4C7", "#B5EAD7", "#A7D8F0")) +
    theme(legend.position = 'bottom', 
          legend.title = element_blank(), 
          axis.title = element_text(face = "bold", size = 18),
          axis.text = element_text(size = 14),
          legend.text = element_text(face = "bold", size = 18),
          text=element_text(family="Arial",face = "bold")) +
    labs(x = expression(Delta~"Brain"),
         y = expression(Delta~"Memory"))
  
  
  gs1 =   
    gs + 
    geom_xsidedensity(data = out.dataset.x, aes(X, group = M, color = M), alpha = .6) +
    geom_ysidedensity(data = out.dataset.y, aes(y =mean, group = M, color = M), alpha = .6) +
    ylim(c(-1.5,.8)) + 
    xlim(c(c(-5,2))) + 
    scale_ysidex_continuous(minor_breaks = NULL, breaks = NULL) + 
    scale_xsidey_continuous(minor_breaks = NULL, breaks = NULL)
  
  ggsave(file.path(outdir.mega.simulation, "meanchange.simulation.png"), 
         plot = gs1, 
         width = 6, 
         height = 8)     
  
  out = 
    list(gs1, 
         df)
  return(out)
}


wrapper_simulation_dispersion =function(outdir.mega.simulation ) {
  load(file.path(outdir.mega.simulation, "modelling.rda")) 
  ## Select analysis
  df = 
    mini.grid %>% 
    filter(between(x1,-.9, -.7),
           skew1 == -5, 
           sd2 == .8)
  
  
  ###
  df.model = rbind(df$smt.model[[1]][[2]] %>% mutate(SD ="low"),
                   df$smt.model[[2]][[2]] %>% mutate(SD ="mid"),
                   df$smt.model[[3]][[2]] %>% mutate(SD ="high"))
  
  out.dataset.x = 
    rbind(
      df$out.dataset[[1]] %>% 
        mutate(g = "g") %>% 
        pivot_longer(-g, 
                     names_to = "G", 
                     values_to = "X") %>% 
        filter(G %in% c("X1")) %>% 
        mutate(SD = "low"), 
      df$out.dataset[[2]] %>% 
        mutate(g = "g") %>% 
        pivot_longer(-g, 
                     names_to = "G", 
                     values_to = "X") %>% 
        filter(G %in% c("X1")) %>% 
        mutate(SD = "mid"), 
      df$out.dataset[[3]] %>% 
        mutate(g = "g") %>% 
        pivot_longer(-g, 
                     names_to = "G", 
                     values_to = "X") %>% 
        filter(G %in% c("X1")) %>% 
        mutate(SD = "high")
    )
  
  
  out.dataset.y = 
    rbind(
      df$out.dataset[[1]] %>% 
        mutate(g = "g") %>% 
        pivot_longer(-g, 
                     names_to = "G", 
                     values_to = "X") %>% 
        filter(G %in% c("Y")) %>% 
        mutate(SD = "low"), 
      df$out.dataset[[2]] %>% 
        mutate(g = "g") %>% 
        pivot_longer(-g, 
                     names_to = "G", 
                     values_to = "X") %>% 
        filter(G %in% c("Y")) %>% 
        mutate(SD = "mid"), 
      df$out.dataset[[3]] %>% 
        mutate(g = "g") %>% 
        pivot_longer(-g, 
                     names_to = "G", 
                     values_to = "X") %>% 
        filter(G %in% c("Y")) %>% 
        mutate(SD = "high")
    ) %>% 
    rename(mean = X)
  
  df.model = 
    df.model %>% 
    mutate(SD = factor(SD, levels = c("low", "mid", "high")))
  
  out.dataset.y = 
    out.dataset.y %>% 
    mutate(SD = factor(SD, levels = c("low", "mid", "high")))
  
  out.dataset.x = 
    out.dataset.x %>% 
    mutate(SD = factor(SD, levels = c("low", "mid", "high")))
  
  gs = ggplot(df.model, aes(X, mean, group = SD, fill = SD, color = SD)) + 
    geom_hline(yintercept = 0, color = "grey50", linetype = "dashed")+ 
    geom_vline(xintercept = 0, color = "grey50", linetype = "dashed")+ 
    geom_line(linewidth = 2) +
    geom_ribbon(aes(ymin = li, ymax = ui), alpha = .6, linewidth = 0)+ 
    theme_classic() +
    scale_fill_manual(values =   c("#FDB4C7", "#B5EAD7", "#A7D8F0")) +
    scale_color_manual(values =   c("#FDB4C7", "#B5EAD7", "#A7D8F0")) +
    theme(legend.position = 'bottom', 
          legend.title = element_blank(), 
          axis.title = element_text(face = "bold", size = 18),
          axis.text = element_text(size = 14),
          legend.text = element_text(face = "bold", size = 18),
          text=element_text(family="Arial",face = "bold")) +
    labs(x = expression(Delta~"Brain"),
         y = expression(Delta~"Memory"))
  
  
  
  gs1 =   
    gs + 
    geom_xsidedensity(data = out.dataset.x, aes(X, group = SD, color = SD), alpha = .6) +
    geom_ysidedensity(data = out.dataset.y, aes(y =mean, group = SD, color = SD), alpha = .6) +
    ylim(c(-1.5,.8)) + 
    xlim(c(c(-6,2))) + 
    scale_ysidex_continuous(minor_breaks = NULL, breaks = NULL) + 
    scale_xsidey_continuous(minor_breaks = NULL, breaks = NULL)
  
  ggsave(file.path(outdir.mega.simulation, "variance.simulation.png"), 
         plot = gs1, 
         width = 6, 
         height = 8)
  
  out = 
    list(gs1, 
         df)
  return(out)
}


wrapper_simulation_skew = function(outdir.mega.simulation) {
  
  load(file.path(outdir.mega.simulation, "modelling.rda"))
  load(file.path(outdir.mega.simulation, "skewgrid.rda"))
  
  ## Select analysis
  df = 
    skew.grid 
  
  ###
  df.model = rbind(df$smt.model[[1]][[2]] %>% mutate(skew ="low"),
                   df$smt.model[[2]][[2]] %>% mutate(skew ="mid"),
                   df$smt.model[[3]][[2]] %>% mutate(skew ="high"))
  
  
  out.dataset.x = 
    rbind(
      df$out.dataset[[1]] %>% 
        mutate(g = "g") %>% 
        pivot_longer(-g, 
                     names_to = "G", 
                     values_to = "X") %>% 
        filter(G %in% c("X1")) %>% 
        mutate(skew = "low"), 
      df$out.dataset[[2]] %>% 
        mutate(g = "g") %>% 
        pivot_longer(-g, 
                     names_to = "G", 
                     values_to = "X") %>% 
        filter(G %in% c("X1")) %>% 
        mutate(skew = "mid"), 
      df$out.dataset[[3]] %>% 
        mutate(g = "g") %>% 
        pivot_longer(-g, 
                     names_to = "G", 
                     values_to = "X") %>% 
        filter(G %in% c("X1")) %>% 
        mutate(skew = "high")
    )
  
  
  out.dataset.y = 
    rbind(
      df$out.dataset[[1]] %>% 
        mutate(g = "g") %>% 
        pivot_longer(-g, 
                     names_to = "G", 
                     values_to = "X") %>% 
        filter(G %in% c("Y")) %>% 
        mutate(skew = "low"), 
      df$out.dataset[[2]] %>% 
        mutate(g = "g") %>% 
        pivot_longer(-g, 
                     names_to = "G", 
                     values_to = "X") %>% 
        filter(G %in% c("Y")) %>% 
        mutate(skew = "mid"), 
      df$out.dataset[[3]] %>% 
        mutate(g = "g") %>% 
        pivot_longer(-g, 
                     names_to = "G", 
                     values_to = "X") %>% 
        filter(G %in% c("Y")) %>% 
        mutate(skew = "high")
    ) %>% 
    rename(mean = X)
  
  
  
  df.model = 
    df.model %>% 
    mutate(skew = factor(skew, levels = c("low", "mid", "high")))
  
  out.dataset.y = 
    out.dataset.y %>% 
    mutate(skew = factor(skew, levels = c("low", "mid", "high")))
  
  out.dataset.x = 
    out.dataset.x %>% 
    mutate(skew = factor(skew, levels = c("low", "mid", "high")))
  
  gs = ggplot(df.model, aes(X, mean, group = skew, fill = skew, color = skew)) + 
    geom_hline(yintercept = 0, color = "grey50", linetype = "dashed")+ 
    geom_vline(xintercept = 0, color = "grey50", linetype = "dashed")+ 
    geom_line(linewidth = 2) +
    geom_ribbon(aes(ymin = li, ymax = ui), alpha = .6, linewidth = 0)+ 
    theme_classic() +
    scale_fill_manual(values =   c("#FDB4C7", "#B5EAD7", "#A7D8F0")) +
    scale_color_manual(values =   c("#FDB4C7", "#B5EAD7", "#A7D8F0")) +
    theme(legend.position = 'bottom', 
          legend.title = element_blank(), 
          axis.title = element_text(face = "bold", size = 18),
          axis.text = element_text(size = 14),
          legend.text = element_text(face = "bold", size = 18),
          text=element_text(family="Arial",face = "bold")) +
    labs(x = expression(Delta~"Brain"),
         y = expression(Delta~"Memory"))
  
  
  
  gs1 =   
    gs + 
    geom_xsidedensity(data = out.dataset.x, aes(X, group = skew, color = skew), alpha = .6) +
    geom_ysidedensity(data = out.dataset.y, aes(y =mean, group = skew, color = skew), alpha = .6) +
    ylim(c(-1.5,.8)) + 
    xlim(c(c(-6,2))) + 
    scale_ysidex_continuous(minor_breaks = NULL, breaks = NULL) + 
    scale_xsidey_continuous(minor_breaks = NULL, breaks = NULL)
  
  ggsave(file.path(outdir.mega.simulation, "skew.simulation.png"), 
         plot = gs1, 
         width = 6, 
         height = 8)
  
  out = 
    list(gs1, 
         df)
  return(out)
}