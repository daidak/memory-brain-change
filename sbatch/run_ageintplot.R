args = commandArgs(TRUE)

feature=as.character(args[1])
wd=as.character(args[2])

# feature = "rh_Lat_Fis-post"

run_ageintplot = function(feature, wd) {
  library(here)
  library(tidyverse)
  library(mgcv)
  library(metagam)
  library(itsadug)
  library(lmerTest)
  library(simex)
  wd = here("data_memory_long/all")
  load(file.path(wd, "merged_reliability.rda"))
  options(bitmapType = "cairo")
  
  
  
  ## META
  df.metagam = 
    df.merged %>% 
    filter(features == feature) %>% 
    mutate(
      data = map(data, ~ filter(.x, abs(memory_slope) < 5 & abs(delta_brain) < 5, ))
    )
  
  df.metagam = 
    df.metagam %>% 
    mutate(mod3 = map(data, ~ gam(memory_slope ~  te(xage, delta_brain), data = .x, gamma = 1)),
           mod3_strip = map(mod3, ~strip_rawdata(.x)))
  
  
  models.mod3 = df.metagam$mod3_strip
  metafit.mod3 <- metagam(models.mod3, 
                          terms = "te(xage,delta_brain)", 
                          grid_size = 25,
                          nsim = 10000, 
                          ci_alpha = .05) 
  
  grid <- list(xage = seq(25, 85, by = 5), 
               delta_brain = seq(-3,3, length.out = 20))
  metafit.mod3.grid <- metagam(models.mod3, 
                               terms = "te(xage,delta_brain)", 
                               grid = grid,
                               nsim = 10000, 
                               ci_alpha = .05)
  
  metamod = 
    metafit.mod3.grid$meta_models$`te(xage,delta_brain)`$predictions %>%  
    group_by(xage) %>% 
    nest() %>% 
    mutate(mod = map(data, ~ lm(estimate ~ delta_brain, data = .x)),
           tidy = map(mod, ~ broom::tidy(.x))) %>% 
    unnest(tidy) %>% 
    filter(term == "delta_brain")
  
  
  
  ## MEGA
  rel_parameters = "/ess/p274/cluster/projects/p039_image_brain/data_normative_long/df_mri/all/reliability/reliability_parameters.rda"
  load(rel_parameters)
  
  
  df.megagam =
    df.merged %>% 
    filter(features == feature) %>% 
    unnest(data) 
  
  
  df.megagam = inner_join(df.megagam, df.features)
  df.megagam = 
    df.megagam %>% 
    mutate(bs = seD^2,
           ws = pct.err.mean^2*((time^2*n*(n+1))/ (12*(n-1)))^-1,
           icc.sub = bs /(ws + bs),
           weight.sub = 1/(1-icc.sub), 
           weight.sub2 = icc.sub^2) %>% 
    mutate(weight.sub = if_else(weight.sub > 10*min(weight.sub), 10*min(weight.sub), weight.sub))
  
  
  mega.mod = 
    gamm(memory_slope ~ te(xage, delta_brain, k = c(5, 5)), 
         data = df.megagam, random=list(dataset=~1), weights =weight.sub2)
  
  

  grid2 <- expand.grid(xage = seq(25, 85, by = 5), 
                      delta_brain = seq(-3,3, length.out = 20))
  
  pred <- predict(mega.mod$gam, newdata = grid2, se.fit = TRUE)
  
  mod2 <- data.frame(
    xage = grid2$xage, 
    delta_brain = grid2$delta_brain, 
    estimate = pred$fit, 
    se = pred$se.fit
  )
  
  megamod = 
    mod2 %>%  
    group_by(xage) %>% 
    nest() %>% 
    mutate(mod = map(data, ~ lm(estimate ~ delta_brain, data = .x)),
           tidy = map(mod, ~ broom::tidy(.x))) %>% 
    unnest(tidy) %>% 
    filter(term == "delta_brain")

  
  
  ##merge and plot
  ##merge
  # mod.mega <- megamod %>% select(-c(data,mod))
  # mod.mega$source <- "Mega"
  # 
  # mod.meta <- metamod %>% select(-c(data,mod))
  # mod.meta$source <- "Meta"
  # 
  # merged <- rbind(mod.mega, mod.meta)
  # 
  # merged <- merged%>%
  #   mutate(
  #     lower_ci = estimate - 1.96 * std.error,
  #     upper_ci = estimate + 1.96 * std.error,
  #     significant = ifelse(p.value < 0.05, "Significant", "Not Significant")
  #   )
  # 
  
  # mergedplot <- ggplot(merged, aes(x = xage, y = estimate)) +
  #   geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = source), alpha = 0.2) +
  #   geom_line(aes(color=source)) +
  #   geom_point(aes(color=significant)) +
  #   scale_color_manual(values = c("Significant" = "green", "Not Significant" = "red")) +
  #   scale_fill_manual(values = c("Mega" = "red", "Meta" = "blue")) +
  #   labs(
  #     title = "Effect of delta_brain on Memory across Age Groups",
  #     x = "Age",
  #     y = "Effect Estimate",
  #     color = "Significance",
  #     fill = "Source"
  #   ) +
  #   theme_minimal() +
  #   theme(legend.position = "bottom")
  

  
  


  
  
  #interaction plot
  #for mega - cant save this plot. So save mega.mod so it can be plotted
  #mega.int <- fvisgam(mod$gam, view=c("xage", "delta_brain"))
  
  
  #for meta
  meta.int <- plot(metafit.mod3, ci = "pointwise", only_meta = TRUE) + 
    theme_bw() + 
    theme(panel.grid = element_blank())
  
  #save plots
  save(mega.mod,
       meta.int,
       file = file.path(wd,"megameta",
                        "temp_models_age",
                        paste0(feature, ".Rda")))
  
  
}

run_ageintplot(feature, wd)