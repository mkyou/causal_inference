library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

N_VALS = c(200, 500, 2000)
T_VALS = c(64, 128, 256)
OUT_DIR = "output"
PLOT_DIR = file.path(OUT_DIR, "consolidated_plots")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# STANDARD COLORS (STRICT)
N_SIM = 1000 # Hardcode simulation count for SE calc
COLOR_MAP_STD = c(
  "Wavelet (Causal/Hybrid)" = "royalblue",
  "Wavelet (Offline)" = "forestgreen",
  "Spline" = "firebrick"
)

# COLOR MAP for raw mapping (just to be safe/consistent)
COLOR_MAP = c(
  "Hybrid" = "royalblue",
  "hybrid" = "royalblue",
  "Causal Wavelet" = "royalblue",
  "rlifting_offline" = "royalblue",
  "Online Wavelet" = "royalblue",
  
  "Ecker (Spline)" = "firebrick",
  "Spline" = "firebrick",
  "spline" = "firebrick",
  
  "wavethresh" = "forestgreen",
  "Offline Wavelet" = "forestgreen"
)

# Standardize Method Names for Plotting
standardize_methods = function(method_vec) {
  sapply(method_vec, function(m) {
    # RoyalBlue Group
    if (m %in% c("hybrid", "Causal Wavelet", "rlifting_offline", "Online Wavelet", "Hybrid")) 
      return("Wavelet (Causal/Hybrid)")
    
    # Firebrick Group
    if (m %in% c("ecker_spline", "Spline", "spline", "Ecker (Spline)")) 
      return("Spline")
    
    # ForestGreen Group
    if (m %in% c("wavethresh", "Offline Wavelet", "Offline")) 
      return("Wavelet (Offline)")
      
    return(m)
  })
}



get_sim1_pointwise = function() {
  results = list()
  
  for (n in N_VALS) {
    for (t_val in T_VALS) {
      path = file.path(
        OUT_DIR, 
        sprintf("sim1_N%d_T%d", n, t_val),
        sprintf("sim1_time_series_results_N%d_T%d.csv", n, t_val)
        )
      if (file.exists(path)) {
        dt = fread(path)
        # Aggregate over time to get global scalar metrics per Scenario/Method
        # Using RMSE average over t as proxy for global performance
        dt_agg = dt[, .(
          NRMSE = mean(NRMSE_t, na.rm=TRUE), 
          # Calculate SE from the SD of the simulations (N_SIM=1000)
          # NRMSE_sd in source is SD across Sims at time t, averaged over t here
          NRMSE_SE = mean(NRMSE_sd, na.rm=TRUE) / sqrt(1000), 
          MAE = mean(MAE_t, na.rm=TRUE)
        ), by = .(Scenario, Method)]
        
        dt_agg[, `:=`(N = n, T_size = t_val)]

        results[[length(results) + 1]] = dt_agg
      }
    }
  }
  return(rbindlist(results))
}

get_sim1_inference = function() {
  results = list()
  
  for (n in N_VALS) {
    for (t_val in T_VALS) {
      path = file.path(
        OUT_DIR, 
        sprintf("sim1_N%d_T%d", n, t_val), "inference",
        sprintf("inference_full_N%d_T%d.csv", n, t_val)
        )
      if (file.exists(path)) {
        dt = fread(path)
        # RAW DATA RETURNED NOW (contains SimID)
        dt[, `:=`(N = n, T_size = t_val)]
        results[[length(results) + 1]] = dt
      }
    }
  }
  return(rbindlist(results))
}

get_sim3_leakage = function() {
  results = list()
  
  for (n in N_VALS) {
    for (t_val in T_VALS) {
      path = file.path(
        OUT_DIR, 
        sprintf("sim3_N%d_T%d", n, t_val), "artifacts",
        sprintf("sim3_raw_results_N%d_T%d.rds", n, t_val)
        )
      
      if (file.exists(path)) {
        raw_list = readRDS(path)
        all_metrics = rbindlist(lapply(raw_list, `[[`, "metrics"))
        
        # RAW DATA RETURNED NOW
        # We need to ensure we have SimID or equivalent if valid, 
        # but here we just need the distribution of Leakage_T1/T2 per Method
        
        all_metrics[, `:=`(N = n, T_size = t_val)]
        results[[length(results) + 1]] = all_metrics
      }
    }
  }
  return(rbindlist(results))
}

create_composite_plot = function(
    data, y_var, y_label, x_var, 
    facet_var, title, method_col = "Method",
    sd_var = NULL # New: Variable containing SD/SE for error bars
    ) {
  
  # Standardize Method Names
  data$MethodStd = standardize_methods(data[[method_col]])
  
  # Scenarios present
  scenarios = unique(data$Scenario)
  if (is.null(scenarios)) scenarios = "Mixed Effects" 
  
  plot_list = list()
  
  # DETERMINE PLOT TYPE
  # 1. Raw Data -> Boxplot (if sim data observed)
  # 2. Aggregated with SD -> Point + ErrorBar
  # 3. Aggregated without SD -> Line + Point
  
  check_grp = data[1, .N, by = c("MethodStd", "Scenario", x_var, facet_var)]
  is_raw_data = any(check_grp$N > 1) 
  
  for (scen in scenarios) {
    dt_sub = data
    if ("Scenario" %in% names(data)) {
      dt_sub = data[Scenario == scen]
    }
    
    if (nrow(dt_sub) == 0) next
    
    # Base Plot (Clean, no global grouping that messes up Boxplots)
    p = ggplot(
      dt_sub, 
      aes_string(
        x = paste0("factor(", x_var, ")"), y = y_var
        )
      ) +
      labs(
        title = scen, 
        x = ifelse(x_var == "T_size", "T (Resolution)", "N (Sample Size)"),
        y = y_label,
        color = "Method",
        fill = "Method"
      ) +
      theme_bw() +
      theme(legend.position = "bottom") + 
      facet_wrap(as.formula(paste("~", facet_var)), labeller = label_both, scales = "fixed")

    
    if (!is.null(sd_var)) {
      # ERROR BAR Plot (Mean +/- 2*SE/SD) -> Takes precedence if SE/SD is provided
      
      p = p + 
        # Error Bars (Black)
        geom_errorbar(
          aes_string(
            ymin = paste0(y_var, " - 2 * ", sd_var),
            ymax = paste0(y_var, " + 2 * ", sd_var),
            group = "MethodStd" 
          ),
          width = 0.2,
          linewidth = 0.5,
          position = position_dodge(0.5), 
          color = "black" 
        ) +
        # Points (Colored by Method)
        geom_point(
          aes_string(color = "MethodStd", group = "MethodStd"),
          size = 3, 
          position = position_dodge(0.5)
        ) +
        scale_color_manual(values = COLOR_MAP_STD)
        
    } else if (is_raw_data) {
      # BOXPLOT (Black Outline, Colored Fill)
      # Only if no sd_var and data is detected as raw
      p = p + geom_boxplot(
        aes_string(fill = "MethodStd"), 
        color = "black",        
        outlier.shape = NA, 
        alpha = 0.7,
        lwd = 0.3
      ) +
      scale_fill_manual(values = COLOR_MAP_STD)
      
    } else {
      # PLAIN LINE PLOT (Fallback)
      p = p + geom_line(
        aes_string(group = "MethodStd", color = "MethodStd"),
        linewidth = 1
      ) +
      geom_point(
        aes_string(group = "MethodStd", color = "MethodStd"),
        size = 3
      ) +
      scale_color_manual(values = COLOR_MAP_STD)
    }


      
    plot_list[[scen]] = p
  }
  
  # Patchwork aggregation
  if (length(plot_list) > 0) {
    final_plot = wrap_plots(plot_list, ncol = 2) + 
      plot_annotation(
        title = NULL,
        theme = theme(plot.title = element_blank())
      ) +
      plot_layout(guides = "collect") & theme(legend.position = "bottom")
    
    return(final_plot)
  } else {
    return(NULL)
  }
}

# --- 3. Execution ---

cat("Reading Sim 1 Pointwise (Aggregated)...\n")
dt_sim1_pt = get_sim1_pointwise()

cat("Reading Sim 1 Inference (Raw)...\n")
dt_sim1_inf = get_sim1_inference()

cat("Reading Sim 3 Leakage (Raw)...\n")
dt_sim3_leak = get_sim3_leakage()

# --- PLOTS START ---

# 1. Sim 1 Pointwise NRMSE (Aggregated -> Mean +/- 2SD)
if (nrow(dt_sim1_pt) > 0) {
  # Reverted to SD per user request to visualize variability
  
  p1 = create_composite_plot(
    dt_sim1_pt[N %in% c(200, 500, 2000)], # Full Grid
    y_var = "NRMSE", y_label = "NRMSE (Mean +/- 2 SE)", 
    x_var = "T_size", facet_var = "N",

    sd_var = "NRMSE_SE",
    title = "Sim 1 Pointwise Error: Fix N, Vary T"
  )
  ggsave(
    file.path(
      PLOT_DIR, "sim1_nrmse_fixN_varyT.png"
      ), p1, width = 14, height = 8
    )
}




# 2. Sim 1 Inference (Coverage & Width)
if (nrow(dt_sim1_inf) > 0) {
  # A) Width - Raw Data -> Box Plot
  p3 = create_composite_plot(
    dt_sim1_inf[N %in% c(200, 500, 2000)],
    y_var = "Width", y_label = "Band Width",
    x_var = "T_size", facet_var = "N",
    title = "Sim 1 Inference Width: Fix N, Vary T"
  )
  ggsave(
    file.path(
      PLOT_DIR, "sim1_width_fixN_varyT.png"), p3, width = 14, height = 8
    )
  
  # B) Coverage - Coverage Rate (Mean +/- 2*SE)
  # We aggregate manually to get Mean and SE (Standard Error of Proportion)
  dt_inf_agg = dt_sim1_inf[, .(
    Coverage_Mean = mean(Cov_Simul, na.rm=TRUE),
    N_Sims = .N
  ), by = .(Scenario, Method, N, T_size)]
  
  # Calculate Standard Error: sqrt(p(1-p)/n)
  dt_inf_agg[, Coverage_SE := sqrt(Coverage_Mean * (1 - Coverage_Mean) / N_Sims)]
  
  p4 = create_composite_plot(
    dt_inf_agg[N %in% c(200, 500, 2000)],
    y_var = "Coverage_Mean", y_label = "Coverage Rate (Mean +/- 2 SE)",
    x_var = "T_size", facet_var = "N",

    sd_var = "Coverage_SE",
    title = "Sim 1 Inference Coverage: Fix N, Vary T"
  )
  
  # Add 0.95 reference line
  p4 = p4 & geom_hline(
    yintercept = 0.95, linetype = "dashed", color = "gray50"
    )
  ggsave(
    file.path(
      PLOT_DIR, "sim1_coverage_fixN_varyT.png"), p4, width = 14, height = 8
    )
}


# 3. Sim 3 Leakage (Raw Data -> Box Plot)
if (nrow(dt_sim3_leak) > 0) {
  
  # Reshape to Long Format to have T1 and T2 as Scenarios
  dt_sim3_long = melt(
    dt_sim3_leak,
    id.vars = c("Method", "N", "T_size"),
    measure.vars = c("Leakage_T1", "Leakage_T2"),
    variable.name = "ScenarioRaw",
    value.name = "LeakageValue"
  )
  
  # Rename for nicer plot titles
  dt_sim3_long[, Scenario := fcase(
    ScenarioRaw == "Leakage_T1", "Leakage T1 (Smooth/Sigmoid)",
    ScenarioRaw == "Leakage_T2", "Leakage T2 (Abrupt/Step)"
  )]
  
  p5 = create_composite_plot(
    dt_sim3_long[N %in% c(200, 500, 2000)],
    y_var = "LeakageValue", y_label = "Leakage (Pre-event Bias)",
    x_var = "T_size", facet_var = "N",
    title = "Sim 3 Leakage: Fix N, Vary T"
  )
  ggsave(
    file.path(
      PLOT_DIR, "sim3_leakage_fixN_varyT.png"
      ), p5, width = 14, height = 7
    )
}

cat("Done. Plots saved to", PLOT_DIR, "\n")
