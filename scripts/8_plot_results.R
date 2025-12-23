# make coefficients plots and partial effect plots

# required library
library(dplyr)
library(ggplot2)
library(ggpubr) # ggarrange

#-------------------------------------------------
# Functions
#-------------------------------------------------
if (TRUE) {
  # ----------- Functions for coefficient plots -----------
  get_coef <- function(filepath) {
    models <- readRDS(filepath)
    coef <- NULL
    for (i in seq_along(models)) {
      if (length(models[[i]]$warnings) == 0) {
        coef <- rbind(coef, models[[i]]$model$coef)
      }
    }
    return(coef)
  }

  df_for_coef_plot <- function(coef, group, nest, ci) {
    x <- (1 - 0.01*ci)/2
    dt <- data.frame(mean=apply(coef, 2, mean), 
                     lower=apply(coef, 2, quantile, x), 
                     upper=apply(coef, 2, quantile, 1-x))
    # arrange the columns as the order to show up in the figure
    if (model_version == "full") {
      dt <- dt[c("pred", "pred:Ground", "pred:Cooperative", "pred:Clutch",
                 "PC1", "PC1:Egg", "PC2", "PC2:Egg", "I(PC2^2)", "I(PC2^2):Egg",
                 "Ground", "Cooperative", "Clutch", "Egg", "Migration"), ]
    } else if (model_version == "pred") {
      dt <- dt[c("pred", "pred:Ground", "pred:Cooperative", "pred:Clutch",
                 "Ground", "Cooperative", "Clutch", "Migration"), ]
    } else if (model_version == "clim") {
      dt <- dt[c("PC1", "PC1:Egg", "PC2", "PC2:Egg", "I(PC2^2)", "I(PC2^2):Egg",
                 "Egg", "Migration"), ]
    }
    
    dt$y <- factor(rownames(dt), levels = rownames(dt))
    dt$sig <- ifelse(dt$lower * dt$upper < 0, 0, 1) # whether significant
    dt$group <- factor(group, levels=c("Passerines", "Non-passerines", "All"))
    dt$nest <- nest
    return(dt)
  }
  
  plot_coef <- function(dt.plot, xlim, ylab) {
    dodge <- position_dodge(width = .5)
    p <- ggplot(dt.plot) +
      geom_vline(aes(xintercept = 0), color = "black", linetype = 1, linewidth = .2) +
      geom_linerange(aes(xmin = lower, xmax = upper, y = y, color = group),
                     position = dodge, linewidth = 2.5, lineend = "round") +
      scale_color_manual(values = c("#50a9b3", "#b37350", "gray35")) +
      geom_point(aes(x = mean, y = y, color = group, fill = factor(sig)),
                 position = dodge, shape = 21, size = 1.8, stroke = 0) +
      scale_fill_manual(values = c("0" = "black", "1" = "white")) +
      scale_y_discrete(limits = rev) +
      coord_cartesian(xlim = xlim) +
      theme_classic() +
      theme(panel.grid = element_blank(),
            text = element_text(size = 12),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title = element_blank(),
            legend.position = "none",
            panel.background = element_rect(fill = 'transparent'),
            plot.background = element_rect(fill = 'transparent', color = NA))
    if (ylab == F) p <- p + theme(axis.text.y = element_blank())
    return(p)
  }
  
  # ----------- Functions for partial effect plots -----------
  process_dataset <- function(data) {
    out <- data %>%
      filter(Y+N > 0) %>%
      select(Migration, Egg, Ground, Cooperative, Clutch, pred, PC1, PC2) %>%
      mutate_at(c("Migration","Egg","Ground","Cooperative","Clutch","pred","PC1","PC2"), ~(scale(.) %>% as.vector))
    return(out)
  }
  
  subset_data <- function(data, nest) {
    if (nest == "enclosed") {
      rs <- data %>% mutate(Y = Enclosed>0, N = Open>0 & Enclosed==0) %>% process_dataset
    } else if (nest == "dome") {
      rs <- data %>% mutate(Y = Dome>0, N = Open>0 & Enclosed==0) %>% process_dataset
    } else if (nest == "cavity") {
      rs <- data %>% mutate(Y = Cavity>0, N = Open>0 & Enclosed==0) %>% process_dataset
    }
    return(rs)
  }
  
  pred_partial <- function(coef, vm, pred = vm$pred, PC1 = vm$PC1, PC2 = vm$PC2, 
                           ground = vm$Ground, coop = vm$Cooperative, 
                           migrate = vm$Migration, mass = vm$Egg, clutch = vm$Clutch) {
    sigmoid <- function(x) 1/(1 + exp(-x))
    # the dataframe columes must be in the same order as the coefficient data
    if (model_version == "pred") {
      df <- data.frame(1, pred, ground, coop, migrate, clutch,
                       pred*ground, pred*coop, pred*clutch) %>% as.matrix
    } else if (model_version == "clim") {
      df <- data.frame(1, PC1, PC2, PC2^2, migrate, mass,
                       PC1*mass, PC2*mass, (PC2^2)*mass) %>% as.matrix
    } else {
      df <- data.frame(1, pred, PC1, PC2, PC2^2,
                       ground, coop, migrate, mass, clutch,
                       pred*ground, pred*coop, pred*clutch,
                       PC1*mass, PC2*mass, (PC2^2)*mass) %>% as.matrix
    }
    rs <- sigmoid(df %*% t(coef))
    mean <- apply(rs, 1, mean)
    upper <- apply(rs, 1, quantile, 0.975)
    lower <- apply(rs, 1, quantile, 0.025)
    return(cbind(mean, upper, lower))
  }
  
  plot_curve_single <- function(df, x, y1, y2) {
    mycol <- c("a" = "dodgerblue", "b" = "darkorange")
    dfplot <- data.frame(x = df[,x], y1 = df[,y1][,1], y2 = df[,y2][,1], 
                         y1max = df[,y1][,2], y1min = df[,y1][,3], 
                         y2max = df[,y2][,2], y2min = df[,y2][,3])
    ggplot(dfplot, mapping = aes(x = x)) + 
      geom_ribbon(mapping = aes(ymax = y1max, ymin = y1min, fill = "a"), alpha = 0.15, show.legend = F) +
      geom_ribbon(mapping = aes(ymax = y2max, ymin = y2min, fill = "b"), alpha = 0.15, show.legend = F) +
      geom_line(mapping = aes(y = y1, color = "a"), linewidth = 1) +
      geom_line(mapping = aes(y = y2, color = "b"), linewidth = 1) +
      scale_color_manual(labels = c(y1, y2), values = mycol) +
      scale_fill_manual(values = mycol) +
      coord_cartesian(ylim = c(0,1), expand = F) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            legend.position = "none",
            legend.background = element_blank(),
            legend.title = element_blank(),
            text = element_text(size = 12))
  }
  
  plot_curve_all <- function(nest, ds, coefs) {
    allplot <- list()
    for (i in 1:3) {
      coef <- coefs[[i]]
      data <- subset_data(ds[[i]], nest)
      if (model_version == "full") {
        data <- data %>% rename(pred = npp)
      }
      vm <- data.frame(t(apply(data, 2, mean))) # mean values of the variables
      
      df_pred <- df_pc1 <- df_pc2 <- data %>% select(pred, PC1, PC2)
      if (model_version != "clim") {
        df_pred$`Off-ground` <-      pred_partial(coef = coef, vm = vm, pred = df_pred$pred, ground = min(data$Ground))
        df_pred$`Ground` <-          pred_partial(coef = coef, vm = vm, pred = df_pred$pred, ground = max(data$Ground))
        df_pred$`Non-cooperative` <- pred_partial(coef = coef, vm = vm, pred = df_pred$pred, coop = min(data$Cooperative))
        df_pred$`Cooperative` <-     pred_partial(coef = coef, vm = vm, pred = df_pred$pred, coop = max(data$Cooperative))
        df_pred$`Small clutch` <-    pred_partial(coef = coef, vm = vm, pred = df_pred$pred, clutch = quantile(data$Clutch, 0.1))
        df_pred$`Large clutch` <-    pred_partial(coef = coef, vm = vm, pred = df_pred$pred, clutch = quantile(data$Clutch, 0.9))
        
        p1 <- plot_curve_single(df_pred, 'pred', 'Ground', 'Off-ground')
        p2 <- plot_curve_single(df_pred, 'pred', 'Non-cooperative', 'Cooperative')
        p3 <- plot_curve_single(df_pred, 'pred', 'Large clutch', 'Small clutch')
      }
      if (model_version != "pred") {
        df_pc1$`Small mass` <-       pred_partial(coef = coef, vm = vm, PC1 = df_pc1$PC1, mass = quantile(data$Egg, 0.1))
        df_pc1$`Large mass` <-       pred_partial(coef = coef, vm = vm, PC1 = df_pc1$PC1, mass = quantile(data$Egg, 0.9))
        df_pc2$`Small mass` <-       pred_partial(coef = coef, vm = vm, PC2 = df_pc2$PC2, mass = quantile(data$Egg, 0.1))
        df_pc2$`Large mass` <-       pred_partial(coef = coef, vm = vm, PC2 = df_pc2$PC2, mass = quantile(data$Egg, 0.9))
        
        p4 <- plot_curve_single(df_pc1, 'PC1', 'Small mass', 'Large mass')
        p5 <- plot_curve_single(df_pc2, 'PC2', 'Small mass', 'Large mass')
      }
      
      if (model_version == "full") allplot <- c(allplot, list(p1, p2, p3, p4, p5))
      if (model_version == "pred") allplot <- c(allplot, list(p1, p2, p3))
      if (model_version == "clim") allplot <- c(allplot, list(p4, p5))
    }
    if (model_version == "full") plot <- ggarrange(plotlist = allplot, ncol = 5, nrow = 3)
    if (model_version == "pred") plot <- ggarrange(plotlist = allplot, ncol = 3, nrow = 3)
    if (model_version == "clim") plot <- ggarrange(plotlist = allplot, ncol = 2, nrow = 3)
    
    return(plot)
  }

  # ----------- Function for preparing data for the plots -----------
  process_all_data <- function(data_version, model_version, multitree) {
    # output folder and filename
    path_output <- paste0("output/model", data_version, "_", model_version, "/")

    # Import bird trait data 
    path_trait <- paste0("data/processed/trait", data_version, "/")
    dt_all <- read.csv(paste0(path_trait ,"sp_traits_pca.csv"), row.names=1) # all apseices
    dt_np <- dt_all %>% filter(Order != "Passeriformes") # non-passerines
    dt_psr <- dt_all %>% filter(Order == "Passeriformes") # passerines
    ds <<- list(dt_all, dt_np, dt_psr)
    
    if (multitree) {
      # Import data (supplementary analysis with 1,000 trees)
      coef_Ae <- get_coef(paste0(path_output, "multitrees_model_all_enclosed.rds"))
      coef_Ad <- get_coef(paste0(path_output, "multitrees_model_all_dome.rds"))
      coef_Ac <- get_coef(paste0(path_output, "multitrees_model_all_cavity.rds"))
      coef_Ne <- get_coef(paste0(path_output, "multitrees_model_np_enclosed.rds"))
      coef_Nd <- get_coef(paste0(path_output, "multitrees_model_np_dome.rds"))
      coef_Nc <- get_coef(paste0(path_output, "multitrees_model_np_cavity.rds"))
      coef_Pe <- get_coef(paste0(path_output, "multitrees_model_psr_enclosed.rds"))
      coef_Pd <- get_coef(paste0(path_output, "multitrees_model_psr_dome.rds"))
      coef_Pc <- get_coef(paste0(path_output, "multitrees_model_psr_cavity.rds"))
    } else {
      ncol <- ifelse(model_version %in% c("full", "full_npp"), 17, 10)
      # Import data (main analysis with consensus tree)
      coef_Ae <- readRDS(paste0(path_output, "model_all_enclosed.rds"))$bootstrap[,-ncol]
      coef_Ad <- readRDS(paste0(path_output, "model_all_dome.rds"))$bootstrap[,-ncol]
      coef_Ac <- readRDS(paste0(path_output, "model_all_cavity.rds"))$bootstrap[,-ncol]
      coef_Ne <- readRDS(paste0(path_output, "model_np_enclosed.rds"))$bootstrap[,-ncol]
      coef_Nd <- readRDS(paste0(path_output, "model_np_dome.rds"))$bootstrap[,-ncol]
      coef_Nc <- readRDS(paste0(path_output, "model_np_cavity.rds"))$bootstrap[,-ncol]
      coef_Pe <- readRDS(paste0(path_output, "model_psr_enclosed.rds"))$bootstrap[,-ncol]
      coef_Pd <- readRDS(paste0(path_output, "model_psr_dome.rds"))$bootstrap[,-ncol]
      coef_Pc <- readRDS(paste0(path_output, "model_psr_cavity.rds"))$bootstrap[,-ncol]
    }
    
    # Set confidence interval for coefficient plots
    ci <- 95
    # ci <- 70
    
    # Get coefficient plot datasets
    dt_enclosed_all <<- rbind(df_for_coef_plot(coef_Ae, "All", "Enclosed", ci),
                              df_for_coef_plot(coef_Ne, "Non-passerines", "Enclosed", ci),
                              df_for_coef_plot(coef_Pe, "Passerines", "Enclosed", ci))
    dt_dome_all <<-     rbind(df_for_coef_plot(coef_Ad, "All", "Dome", ci),
                              df_for_coef_plot(coef_Nd, "Non-passerines", "Dome", ci),
                              df_for_coef_plot(coef_Pd, "Passerines", "Dome", ci))
    dt_cavity_all <<-   rbind(df_for_coef_plot(coef_Ac, "All", "Cavity", ci),
                              df_for_coef_plot(coef_Nc, "Non-passerines", "Cavity", ci),
                              df_for_coef_plot(coef_Pc, "Passerines", "Cavity", ci))

    # Correction of round-end line range for plotting purpose (output into 8x6.6in pdf image)
    adj <- 0.012
    dt_enclosed_all$lower <<- dt_enclosed_all$lower + adj
    dt_enclosed_all$upper <<- dt_enclosed_all$upper - adj
    dt_dome_all$lower <<- dt_dome_all$lower + adj
    dt_dome_all$upper <<- dt_dome_all$upper - adj
    dt_cavity_all$lower <<- dt_cavity_all$lower + adj
    dt_cavity_all$upper <<- dt_cavity_all$upper - adj
    
    # Prepare partial effect datasets
    coefs_enc <<- list(coef_Ae, coef_Ne, coef_Pe)
    coefs_dom <<- list(coef_Ad, coef_Nd, coef_Pd)
    coefs_cav <<- list(coef_Ac, coef_Nc, coef_Pc)
  }
}

#-------------------------------------------------
# MAIN
#-------------------------------------------------
# model version
# Main model (global, consensus tree)
data_version = ""; model_version = "pred"; multitree = FALSE
data_version = ""; model_version = "clim"; multitree = FALSE
# Global, multiple trees
data_version = ""; model_version = "pred"; multitree = TRUE
data_version = ""; model_version = "clim"; multitree = TRUE
# Northern hemisphere only, consensus tree
data_version = "_north"; model_version = "pred"; multitree = FALSE
data_version = "_north"; model_version = "clim"; multitree = FALSE
# # Northern hemisphere only, multiple trees
data_version = "_north"; model_version = "pred"; multitree = TRUE
data_version = "_north"; model_version = "clim"; multitree = TRUE

# plot and save
if (TRUE) {
  process_all_data(data_version, model_version, multitree)
  suffix <- ifelse(isTRUE(multitree), "_multi.pdf", ".pdf")

  # Forest plots
  xlim <- if (data_version == "_north") c(-.5, .5) else c(-.3, .3)
  p1 <- plot_coef(dt_enclosed_all, xlim, F)
  p2 <- plot_coef(dt_dome_all, xlim, F)
  p3 <- plot_coef(dt_cavity_all, xlim, F)
  ggarrange(p1, p2, p3, ncol = 3, nrow = 1)
  ggsave(filename = paste0("output/img/forest_", model_version, data_version, suffix), width = 8, height = 4)
  # ggsave(filename = paste0("output/img/forest_70_", model_version, data_version, suffix), width = 8, height = 4)

  # Partial effect plots
  pep_width <- ifelse(model_version == "pred", 6, ifelse(model_version == "clim", 4, 10))
  plot_curve_all("enclosed", ds, coefs_enc)
  ggsave(filename = paste0("output/img/curve_enclosed_", model_version, data_version, suffix), width = pep_width, height = 6)
  plot_curve_all("dome", ds, coefs_dom)
  ggsave(filename = paste0("output/img/curve_dome_", model_version, data_version, suffix), width = pep_width, height = 6)
  plot_curve_all("cavity", ds, coefs_cav)
  ggsave(filename = paste0("output/img/curve_cavity_", model_version, data_version, suffix), width = pep_width, height = 6)
}
 
