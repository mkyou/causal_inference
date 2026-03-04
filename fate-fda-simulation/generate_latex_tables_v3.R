library(data.table)
s1 = fread('tmp_s1.csv')
s2 = fread('tmp_s2.csv')

# Helper to format numbers with comma as decimal
fmt = function(x, d=3) {
  if (is.na(x)) return("-")
  format(round(x, d), nsmall=d, decimal.mark=',', big.mark='.')
}

# 1. SIMPLIFIED TABLE (Main Text)
s1_simple = s1[Scenario %in% c("Sigmoid", "Doppler", "Escada") & N_fac %in% c(200, 500, 2000) & T_fac %in% c(128, 512)]

sink('table1_simple_rows_v3.txt')
for (scen in c("Sigmoid", "Doppler", "Escada")) {
  sub_scen = s1_simple[Scenario == scen]
  if(nrow(sub_scen) == 0) next
  cat(sprintf('\\midrule\n\\multicolumn{8}{c}{\\textbf{Cenário: %s}} \\\\\n\\midrule\n', scen))
  for (n in c(200, 500, 2000)) {
    for (t in c(128, 512)) {
      sub = sub_scen[N_fac == n & T_fac == t]
      methods_order = c("Ondaletas (Online)", "Ondaletas (Offline)", "Splines de Ecker")
      for (m in methods_order) {
        row = sub[Method == m]
        if(nrow(row) == 0) next
        nrmse = fmt(row$NRMSE_Mean, 3)
        rmse = fmt(row$RMSE_Mean, 3)
        bias = fmt(row$Bias_Mean, 3)
        cov = paste0(fmt(row$Cov_Simul_Mean*100, 1), '\\%')
        m_label = switch(m, 
          "Ondaletas (Online)" = "Ond. (On)",
          "Ondaletas (Offline)" = "Ond. (Off)",
          "Splines de Ecker" = "Ecker"
        )
        cat(sprintf('%s & %s & %d & %d & %s & %s & %s & %s \\\\\n', "", m_label, n, t, rmse, nrmse, bias, cov))
      }
      cat('\\cmidrule{2-8}\n')
    }
  }
}
sink()

# 2. FULL TABLE (Appendix)
sink('table1_full_rows_v3.txt')
for (scen in sort(unique(s1$Scenario))) {
  cat(sprintf('\\midrule\n\\multicolumn{8}{c}{\\textbf{Cenário: %s}} \\\\\n\\midrule\n', scen))
  for (n in sort(unique(s1$N_fac))) {
    for (t in sort(unique(s1$T_fac))) {
      sub = s1[Scenario == scen & N_fac == n & T_fac == t]
      methods_order = c("Ondaletas (Online)", "Ondaletas (Offline)", "Splines de Ecker")
      for (m in methods_order) {
        row = sub[Method == m]
        if(nrow(row) == 0) next
        nrmse = fmt(row$NRMSE_Mean, 3)
        rmse = fmt(row$RMSE_Mean, 3)
        bias = fmt(row$Bias_Mean, 3)
        cov = paste0(fmt(row$Cov_Simul_Mean*100, 1), '\\%')
        m_label = switch(m, 
          "Ondaletas (Online)" = "Ond. (On)",
          "Ondaletas (Offline)" = "Ond. (Off)",
          "Splines de Ecker" = "Ecker"
        )
        cat(sprintf('%s & %s & %d & %d & %s & %s & %s & %s \\\\\n', "", m_label, n, t, rmse, nrmse, bias, cov))
      }
      cat('\\cmidrule{2-8}\n')
    }
  }
}
sink()

# 3. SIMULATION 2 TABLE
sink('table2_rows_v3.txt')
for (scen in sort(unique(s2$Scenario))) {
  cat(sprintf('\\midrule\n\\multicolumn{7}{c}{\\textbf{Cenário: %s}} \\\\\n\\midrule\n', scen))
  for (t in sort(unique(s2$T_fac))) {
    for (n in sort(unique(s2$N_fac))) {
      row = s2[Scenario == scen & T_fac == t & N_fac == n]
      if(nrow(row) == 0) next
      nrmse = fmt(row$NRMSE_Mean, 3)
      rmse = fmt(row$RMSE_Mean, 3)
      bias = fmt(row$Bias_Mean, 3)
      lat = fmt(row$Latency_Mean, 1)
      cat(sprintf('%s & %d & %d & %s & %s & %s & %s \\\\\n', "", t, n, rmse, nrmse, bias, lat))
    }
    cat('\\cmidrule{2-7}\n')
  }
}
sink()
