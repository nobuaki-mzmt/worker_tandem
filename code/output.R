# output.R
# load format data created by procesing.R
# create plots and perform stats

# config ----
{
  source("code/source.R")
  
  treat_labels <- c(
    FM  = "Male - Female",
    MW  = "Male - Worker (same)",
    MW2 = "Male - Worker (other)",
    SM  = "Male - Soldier",
    MM  = "Male - Male"
  )
  
  treat_colors <- c(
    FM  = "#D55E00",  # vermillion
    MM  = "#0072B2",  # blue
    MW  = "#009E73",  # green (same colony worker)
    SM  = "#66A61E",  # olive green (same colony soldier)
    MW2 = "#CC79A7"   # purple (different colony worker)
  )

  load("data_fmt/df_pair.rda")
  load("data_fmt/df_dish.rda")
  
  df_attack <- read.csv("data_raw/master_well.csv")
  df_attack <- df_attack |> filter(!is.na(tandem_before_attack ), treat == "MW2") |> 
    mutate(pair_id = paste(Video, well, sep = "_")) |> 
    dplyr::select(pair_id, treat, tandem_before_attack, severe_attack_sec, fetal_attack_sec, tandem_after_attack)
    
  # remove sec after the fetal attack
  df_pair <- df_pair |> 
    left_join(df_attack |> dplyr::select(pair_id, fetal_attack_sec), by = "pair_id") |> 
    mutate(fetal_attack_sec = if_else(is.na(fetal_attack_sec), 1801, fetal_attack_sec)) |> 
    filter(time_sec < fetal_attack_sec) |> dplyr::select(-fetal_attack_sec)
}

# plot relative positions ----
{
  df_rel <- df_pair |>
    mutate(rx0 = partner_dis * cos(dir_1to0),
           ry0 = partner_dis * sin(dir_1to0),
           rx1 = partner_dis * cos(dir_0to1),
           ry1 = partner_dis * sin(dir_0to1))
  
  video_list <- unique(df_rel$video )
  
  for(i_v in video_list){
    df_temp <- df_rel |> filter(video == i_v)
    p1 <- ggplot(df_temp, aes(x = rx0, y = ry0)) +
      geom_bin_2d(binwidth = 0.1) +
      facet_wrap(~well_id) +
      scale_fill_viridis() +
      coord_cartesian(xlim = c(-2,2), ylim = c(-2,2)) +
      theme(legend.position = "none", aspect.ratio = 1) +
      labs(x = "", y = "", title = i_v)
    ggsave(p1, file = sprintf("output/relative_pos_each/%s.pdf", i_v), width = 4, height = 3)
  }
  
  ggplot(df_rel, aes(x = rx0, y = ry0)) +
    geom_bin_2d(binwidth = 0.1) +
    facet_wrap(~treat) +
    scale_fill_viridis() +
    coord_cartesian(xlim = c(-2,2), ylim = c(-2,2)) +
    theme(legend.position = "none", aspect.ratio = 1) +
    labs(x = "", y = "")
  ggsave(file = sprintf("output/relative_pos.pdf"), width = 4, height = 3)
}

# tandem prop time development ----
{
  # data prep + plot
  {
    plot_tandem_prop_time <- function(df, time_bin = 60){ 
      df_tandem_time <- df |>
        mutate(time_min = floor(time_sec / time_bin)) |>
        group_by(time_min, treat, pair_id, colony) |>
        summarize(
          interact_prop = mean(interact),
          tandem_prop = mean(tandem),
          non_tandem_interact_prop = mean(!tandem & interact),
          .groups = "drop"
        )
      p <- ggplot(df_tandem_time,
             aes(x = time_min, y = tandem_prop, col = treat)) +
        stat_summary(fun = mean, geom = "line", linewidth = 1) +
        stat_summary(fun.data = mean_se, geom = "ribbon",
                     aes(fill = treat), alpha = 0.2, color = NA) +
        scale_color_manual(values = treat_colors, labels = treat_labels)+
        scale_fill_manual(values = treat_colors, labels = treat_labels)  +
        scale_y_continuous(limits = c(0, 1), breaks = c(0,0.5,1), labels = c(0,0.5,1)) +
        theme_classic(base_size = 10) +
        labs(x = "Time (minutes)", y = "Propotion of tandem")+
        theme(aspect.ratio = 1,
              legend.position = c(0.7,0.9),
              legend.title = element_blank(),
              legend.background = element_blank())
      list(df_tandem_time, p)
    }
    
    df_out <- plot_tandem_prop_time(df_pair, time_bin = 60)
    df_tandem_prop_well = df_out[[1]]
    ggsave(df_out[[2]], filename = "output/tandem_prop_well.pdf", 
            device = cairo_pdf, family = "Arial",
            width = 3, height = 3)
    
    df_out <- plot_tandem_prop_time(df_dish, time_bin = 60)
    df_tandem_prop_dish = df_out[[1]]
    ggsave(df_out[[2]], filename = "output/tandem_prop_well.pdf", 
           device = cairo_pdf, family = "Arial",
           width = 3, height = 3)
    
  }
  
  # stat
  {
    # well
    {
      df_tandem_prop_well |> filter(tandem_prop > 0, tandem_prop < 1) |> pull(tandem_prop) |> range()
      df_tandem_prop_well <- df_tandem_prop_well |> 
        mutate(logit_tandem_prop = car::logit(tandem_prop, adjust = 0.003, percents = FALSE))
      r <- lmer(logit_tandem_prop ~ time_min * treat + (1 | colony/pair_id), data = df_tandem_prop_well)
      df_res <- tibble(tidy(Anova(r)), formula = deparse(formula(r) ))
      res <- summary(r)
      pairwise_res <- emtrends(r, pairwise ~ treat, var = "time_min")
      
      wb <- createWorkbook()
      addWorksheet(wb, "anova")
      addWorksheet(wb, "slopes")
      addWorksheet(wb, "contrasts")
      
      writeData(wb, "anova",     df_res)
      writeData(wb, "slopes",    as.data.frame(pairwise_res$emtrends))
      writeData(wb, "contrasts", as.data.frame(pairwise_res$contrasts))
      
      saveWorkbook(wb, "output/tandem_prop_well_lmer.xlsx", overwrite = TRUE)
    }
    # dish
    {
      df_tandem_prop_dish |> filter(tandem_prop > 0, tandem_prop < 1) |> pull(tandem_prop) |> range()
      df_tandem_prop_dish <- df_tandem_prop_dish |> 
        mutate(logit_tandem_prop = car::logit(tandem_prop, adjust = 0.003, percents = FALSE))
      r <- lmer(logit_tandem_prop ~ time_min * treat + (1 | colony/pair_id), data = df_tandem_prop_dish)
      df_res <- tibble(tidy(Anova(r)), formula = deparse(formula(r) ))
      res <- summary(r)
      pairwise_res <- emtrends(r, pairwise ~ treat, var = "time_min")
      slope_tests <- test(pairwise_res$emtrends)
      
      wb <- createWorkbook()
      addWorksheet(wb, "anova")
      addWorksheet(wb, "slopes")
      addWorksheet(wb, "slopes_P")
      addWorksheet(wb, "contrasts")
      
      writeData(wb, "anova",     df_res)
      writeData(wb, "slopes",    as.data.frame(pairwise_res$emtrends))
      writeData(wb, "slopes_P",    as.data.frame(slope_tests ))
      writeData(wb, "contrasts", as.data.frame(pairwise_res$contrasts))
      
      saveWorkbook(wb, "output/tandem_prop_well_lmer.xlsx", overwrite = TRUE)
    }
  }
  
}

# tandem duration ----
{
  # create tandem event data
  tandem_df <- function(df,
                        gap_max = 1, # number of frames allowed as gap 
                        min_len = 1 ){ # minimum frames to count as tandem
    df_tandem <- df|>
      group_by(pair_id) |>
      arrange(time_sec, .by_group = TRUE) |>
      mutate(
        run_id = rleid(tandem),
        run_len = ave(tandem, run_id, FUN = length)
      ) |>
      mutate(
        tandem = if_else(tandem & run_len >= min_len, TRUE, FALSE)
      ) |>
      mutate(
        cluster_start = tandem & lag(!tandem, default = TRUE),
        cluster_idx = cumsum(cluster_start),
        cluster_idx = if_else(tandem, cluster_idx, NA_integer_),
        pair_event = if_else(
          !is.na(cluster_idx),
          sprintf("%s_%02d", pair_id, cluster_idx),
          NA_character_
        )
      ) |>
      dplyr::select(-cluster_start, -run_id, -run_len) |>
      ungroup()
    
    df_pair_dur <- df_tandem |> filter(time_sec < 1800.1) |>
      filter(!is.na(pair_event)) |>
      group_by(pair_event, treat, pair_id) |>
      summarise(
        start_time = first(time_sec),
        end_time = last(time_sec),
        tandem_event = first(cluster_idx),
        duration = (n() * 0.2), 
        cens = !(dplyr::near(last(time_sec), max(df$time_sec)) | dplyr::near(first(time_sec), 0)),
        leader_moved_dis = sum(post_step_0, na.rm = T),
        follower_moved_dis = sum(post_step_1, na.rm = T),
        .groups = "drop"
      )
    df_pair_dur
  }
  
  # dish analysis
  df_tandem <- tandem_df(df_dish, gap_max = 1, min_len = 1)
  
  df_first_tandem <- df_tandem |> group_by(pair_id) |>
    slice_min(start_time, n = 1, with_ties = TRUE) |>
    ungroup()
  
  ggplot(df_first_tandem, aes(x = treat, y = start_time)) +
    geom_boxplot() +
    geom_jitter()
    
  ggplot(df_first_tandem, aes(x = treat, y = duration)) +
    geom_boxplot() +
    geom_jitter()
  
  ggsurvplot(
    survfit(Surv(duration, cens) ~ treat, data = df_first_tandem),
    data = df_first_tandem,
    censor = FALSE,
    conf.int.style = "ribbon",
    conf.int.fill = TRUE,
    ggtheme = theme_classic(),
    color = "treat"
  )$plot + 
    scale_x_continuous(trans = "pseudo_log", breaks = c(0, 1, 10, 100, 1000))+
    scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7"), labels = treat_labels)  +
    scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7"), labels = treat_labels) +
    scale_y_continuous(breaks = c(0,1))+
    labs(x = "Duration (sec)", y = "Tandem probability") +
    theme(aspect.ratio = 3/4,
          legend.title = element_blank(),
          legend.position = c(0.8,0.8))
  
  
  df_non_first_tandem <- df_tandem |>
    group_by(pair_id) |>
    filter(start_time != min(start_time)) |>
    ungroup()
  
  ggsurvplot(
    survfit(Surv(duration, cens) ~ treat, data = df_non_first_tandem),
    data = df_longest_tandem,
    censor = TRUE,
    conf.int.style = "ribbon",
    conf.int.fill = TRUE,
    ggtheme = theme_classic(),
    color = "treat"
  )$plot + 
    scale_x_continuous(trans = "pseudo_log", breaks = c(0, 1, 10, 100, 1000))+
    scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7"), labels = treat_labels)  +
    scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7"), labels = treat_labels) +
    scale_y_continuous(breaks = c(0,1))+
    labs(x = "Duration (sec)", y = "Tandem probability") +
    theme(aspect.ratio = 3/4,
          legend.title = element_blank(),
          legend.position = c(0.8,0.8))
  
  df_combined <- bind_rows(
    "non_first" = df_non_first_tandem,
    "first"   = df_first_tandem,
    .id = "tandem_type"
  )
  
  ggsurvplot(
    survfit(Surv(duration, cens) ~ tandem_type + treat, data = df_combined),
    data = df_combined,
    censor = T,
    conf.int.style = "ribbon",
    conf.int.fill = TRUE,
    ggtheme = theme_classic()
  )$plot + 
    scale_x_continuous(trans = "pseudo_log", breaks = c(0, 1, 10, 100, 1000))+
    scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7"), labels = treat_labels)  +
    scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7"), labels = treat_labels) +
    scale_y_continuous(breaks = c(0,1))+
    labs(x = "Duration (sec)", y = "Tandem probability") +
    theme(aspect.ratio = 3/4,
          legend.title = element_blank(),
          legend.position = c(0.8,0.8))
  
  df_combined$treat <- as.factor(df_combined$treat)
  fit_cox <- coxme(Surv(duration, cens) ~ tandem_type + (1|pair_id), 
                   data = df_combined |> filter(treat == "FM"))
  Anova(fit_cox)
  
  fit_cox <- coxme(Surv(duration, cens) ~ tandem_type + (1|pair_id), 
                   data = df_combined |> filter(treat == "MM"))
  Anova(fit_cox)
  
  
  
  df_longest_tandem <- df_tandem |> group_by(pair_id) |>
    slice_max(duration, n = 1, with_ties = TRUE) |>
    ungroup()
  
  ggplot(df_longest_tandem, aes(x = treat, y = start_time)) +
    geom_boxplot() +
    geom_jitter()
  
  ggplot(df_longest_tandem, aes(x = treat, y = tandem_event)) +
    geom_boxplot() +
    geom_jitter()
  
  
  ggplot(df_longest_tandem, aes(x = duration, fill = treat)) +
    geom_density(alpha = 0.25)+
    scale_x_log10(
      breaks = c(1, 10, 100, 1000),
      labels = scales::label_number()
    ) 
  
  ggplot(df_longest_tandem, aes(x = duration, color = treat)) +
    stat_ecdf(linewidth = 0.8) +
    scale_color_viridis_d(end = 0.8) +
    labs(x = "Duration (s)", y = "Cumulative proportion",
         color = "Treatment") +
    theme_classic() +
    theme(aspect.ratio = 1, legend.position = c(0.8, 0.3))+scale_x_log10(
      breaks = c(1, 10, 100, 1000),
      labels = scales::label_number()
    ) 
  
  
  ggplot(df_longest_tandem, aes(x = duration, y = treat, fill = treat, color = treat)) +
    stat_halfeye( width = 0.4, position = position_nudge(y = 0.11),
                  .width = 0, point_colour = NA, alpha = 0.7, scale = 0.5) +
    stat_dots( side = "bottom", dotsize = 0.75, binwidth = 0.08,
      alpha = 0.6, position = position_nudge(y = -0.1) ) +
    geom_boxplot( width = 0.06, outlier.shape = NA,
                  fill = "white", alpha = 0.6, color = "grey30", linewidth = 0.4 ) +
    scale_x_log10( breaks = c(1, 10, 100, 1000), labels = scales::label_number()) +
    coord_cartesian(ylim = c(1,  2.2)) +
    scale_fill_viridis_d(end = 0.8) +
    scale_color_viridis_d(end = 0.8) +
    labs(x = "Longest tandem duration (s)", y = NULL) +
    theme_classic(base_size = 10) +
    theme(aspect.ratio = 1.25, legend.position = "none")
  
  # OVL from KDE
  d1 <- density(log10(df_longest_tandem$duration[df_longest_tandem$treat == "FM"]))
  d2 <- density(log10(df_longest_tandem$duration[df_longest_tandem$treat == "MM"]))
  # Interpolate to common grid
  x_grid <- seq(min(d1$x, d2$x), max(d1$x, d2$x), length.out = 1000)
  f1 <- approx(d1$x, d1$y, x_grid)$y
  f2 <- approx(d2$x, d2$y, x_grid)$y
  ovl <- sum(pmin(f1, f2, na.rm = TRUE)) / sum(pmax(f1, f2, na.rm = TRUE))
  
  
  
  library(effsize)   # cliff.delta()
  
  cd <- cliff.delta(duration ~ treat, data = df_longest_tandem)
  
  ggsurvplot(
    survfit(Surv(duration, cens) ~ treat, data = df_longest_tandem),
    data = df_longest_tandem,
    censor = T,
    conf.int.style = "ribbon",
    conf.int.fill = TRUE,
    ggtheme = theme_classic()
  )$plot + 
    scale_x_continuous(trans = "pseudo_log", breaks = c(0, 1, 10, 100, 1000))+
    scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7"), labels = treat_labels)  +
    scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7"), labels = treat_labels) +
    scale_y_continuous(breaks = c(0,1))+
    labs(x = "Duration (sec)", y = "Tandem probability") +
    theme(aspect.ratio = 3/4,
          legend.title = element_blank(),
          legend.position = c(0.8,0.8))
  
  
  
  df_combined <- bind_rows(
    "longest" = df_longest_tandem,
    "first"   = df_first_tandem,
    .id = "tandem_type"
  )
  
  ggsurvplot(
    survfit(Surv(duration, cens) ~ treat + tandem_type, data = df_combined),
    data = df_combined,
    censor = T,
    conf.int.style = "ribbon",
    conf.int.fill = TRUE,
    ggtheme = theme_classic()
  )$plot + 
    scale_x_continuous(trans = "pseudo_log", breaks = c(0, 1, 10, 100, 1000))+
    scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7"), labels = treat_labels)  +
    scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7"), labels = treat_labels) +
    scale_y_continuous(breaks = c(0,1))+
    labs(x = "Duration (sec)", y = "Tandem probability") +
    theme(aspect.ratio = 3/4,
          legend.title = element_blank(),
          legend.position = c(0.8,0.8))
  
  
  posthoc <- summary(glht(fit_cox, linfct = mcp(treat = "Tukey")))
  
  
  # tandem duration
  {
    ggsurvplot(
      survfit(Surv(duration, cens) ~ treat, data = df_tandem),
      data = df_plot,
      censor = FALSE,
      conf.int.style = "ribbon",
      conf.int.fill = TRUE,
      ggtheme = theme_classic(),
      color = "treat"
    )$plot + 
      scale_x_continuous(trans = "pseudo_log", breaks = c(0, 1, 10, 100, 1000))+
      scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7"), labels = treat_labels)  +
      scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7"), labels = treat_labels) +
      scale_y_continuous(breaks = c(0,1))+
      labs(x = "Duration (sec)", y = "Tandem probability") +
      theme(aspect.ratio = 3/4,
            legend.title = element_blank(),
            legend.position = c(0.8,0.8))
    
    ggsave( filename = "output/tandem_duration.pdf", 
            device = cairo_pdf, family = "Arial",
            width = 4, height = 3)
  }
  
  # cumulative hazard plot
  {
    ggsurvplot(
      survfit(Surv(duration, cens) ~ treat, data = df_plot),
      data = df_plot,
      fun = "cumhaz",
      censor = TRUE,
      ggtheme = theme_classic(),
      color = "treat"
    )$plot+ 
      scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7"), labels = treat_labels)  +
      scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7"), labels = treat_labels) +
      scale_y_continuous(breaks = c(0, 5))+
      coord_cartesian(ylim = c(0, 6.5)) +
      labs(x = "Duration (sec)", y = "Cumulative hazard") +
      theme(aspect.ratio = 3/4,
            legend.title = element_blank(),
            legend.position = c(0.8,0.4))
    
    ggsave( filename = "output/tandem_duration_cum_hazard.pdf", 
            device = cairo_pdf, family = "Arial",
            width = 4, height = 3)
  }
  
  # stat
  {
    df_plot$treat <- as.factor(df_plot$treat)
    fit_cox <- coxme(Surv(duration, cens) ~ treat + (1|video/well_id), data = df_plot)
    res <- Anova(fit_cox)
    posthoc <- summary(glht(fit_cox, linfct = mcp(treat = "Tukey")))
    
    df_res <- tibble(tibble(res), formula = deparse(formula(fit_cox) ))
    write.csv(df_res, "output/tandem_duration_coxme.csv")
    df_res <- tidy(posthoc)
    write.csv(df_res, "output/tandem_duration_coxme_posthoc_Tukey.csv")
  }
}

# interaction transition ---- 
{
  
  df <- df_dish
  
  
  df_choice <- df %>%
    filter(!is.na(next_state)) |>
    filter(state == "interact", next_state %in% c("tandem", "sep"), state != next_state) |>
    mutate(time_min = floor(time_sec / 300)) |>
    group_by(time_min, treat, pair_id, colony) |>
    summarize(
      sep_prop = mean(next_state == "sep"),
      tandem_prop = mean(next_state == "tandem"),
      .groups = "drop"
    )
  
  ggplot(df_choice, aes(x = as.factor(time_min), y = sep_prop, fill = treat, col = treat)) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5) +
    scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7"))  +
    labs(y = "P(tandem | interact)", x = NULL) +
    scale_x_discrete(labels = treat_labels) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none", aspect.ratio = 1)
  
  
  %>%
    count(treat, pair_id, next_state) %>%
    mutate(prob = n / sum(n),.by = c(treat, pair_id))
  
  df_choice_tandem <- df_choice %>%
    filter(next_state == "tandem")
  
  # plot
  {
    ggplot(df_choice_tandem, aes(treat, prob, fill = treat, col = treat)) +
      geom_jitter(width = 0.1, alpha = 0.6, size = 1) +
      stat_summary(fun = median, geom = "crossbar", width = 0.5) +
      scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7"))  +
      labs(y = "P(tandem | interact)", x = NULL) +
      scale_x_discrete(labels = treat_labels) +
      coord_cartesian(ylim = c(0, 1)) +
      scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
      theme_classic(base_size = 10) +
      theme(legend.position = "none", aspect.ratio = 1)
    
    ggsave( filename = "output/Interaction_transition.pdf", 
            device = cairo_pdf, family = "Arial",
            width = 3.5, height = 3.5)
  }
  
  # stat
  df_stat <- df_state %>%
    filter(state == "interact", next_state %in% c("tandem", "sep"), state != next_state) %>%
    mutate(pair_id = paste(video, well_id, sep = "_"),
           tandem_transition = next_state == "tandem")
  
  r <- glmer(tandem_transition ~ treat + (1 | pair_id), data = df_stat,
             family = binomial)
  r_sum <- summary(r)
  res <- Anova(r)
  posthoc <- summary(glht(r, linfct = mcp(treat = "Tukey")))
  
  df_res <- tibble(tibble(res), formula = deparse(formula(r) ))
  write.csv(df_res, "output/Interaction_transition_glmer.csv")
  df_res <- tidy(posthoc)
  write.csv(df_res, "output/Interaction_transition_glmer_posthoc_Tukey.csv")
  
  
  ##
  df_choice <- df_dish %>%
    filter(!is.na(next_state)) |>
    filter(state == "interact", next_state %in% c("tandem", "sep"), state != next_state) %>%
    count(treat, video, next_state) %>%
    mutate(prob = n / sum(n),.by = c(treat, video))
  
  df_choice_tandem <- df_choice %>%
    filter(next_state == "tandem")
  
  # plot
  {
    ggplot(df_choice_tandem, aes(treat, prob, fill = treat, col = treat)) +
      geom_jitter(width = 0.1, alpha = 0.6, size = 1) +
      stat_summary(fun = median, geom = "crossbar", width = 0.5) +
      scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7"))  +
      labs(y = "P(tandem | interact)", x = NULL) +
      scale_x_discrete(labels = treat_labels) +
      coord_cartesian(ylim = c(0, 1)) +
      scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
      theme_classic(base_size = 10) +
      theme(legend.position = "none", aspect.ratio = 1)
    
    ggsave( filename = "output/Interaction_transition.pdf", 
            device = cairo_pdf, family = "Arial",
            width = 3.5, height = 3.5)
    }
  
  # stat
  df_stat <- df_state %>%
    filter(state == "interact", next_state %in% c("tandem", "sep"), state != next_state) %>%
    mutate(pair_id = video,
           tandem_transition = next_state == "tandem")
  
  r <- glmer(tandem_transition ~ treat + (1 | pair_id), data = df_stat,
             family = binomial)
  r_sum <- summary(r)
  res <- Anova(r)
  posthoc <- summary(glht(r, linfct = mcp(treat = "Tukey")))
  
}  



# speed ----
{
  df_dis_dist <- df_pair |> pivot_longer(cols = starts_with("post_step")) |> 
    dplyr::select(name, value, treat)
  
  ggplot(df_dis_dist |> filter(value < 1)) +
    geom_density_ridges(aes(x = value, y = treat, fill= name), stat = "binline", 
                        alpha = 0.5, binwidth = 0.04, scale = 0.85) +
    labs(x = "Step length (BL)", y = "") +
    scale_y_discrete(expand = c(0, 0.1),labels = treat_labels) +
    scale_x_continuous(breaks = c(0,0.5,1), labels =  c(0,0.5,1)) +
    scale_fill_manual(values = c(post_step_0 = "#1B7837", post_step_1 = "#D8B58A")) +
    coord_cartesian(xlim = c(0,1)) +
    theme_classic(base_size = 10) +
    theme(aspect.ratio = 3,
          legend.position = "none")
  
  ggsave("output/step_distribution.pdf", device = cairo_pdf, family = "Arial",
         width = 3, height = 5.5)
  
  df_dis_dist <- df_dish |> pivot_longer(cols = starts_with("post_step")) |> 
    dplyr::select(name, value, treat)
  
  
  ggplot(df_dis_dist |> filter(value < 1)) +
    geom_density_ridges(aes(x = value, y = treat, fill= name), stat = "binline", 
                        alpha = 0.5, binwidth = 0.04, scale = 0.85) +
    labs(x = "Step length (BL)", y = "") +
    scale_y_discrete(expand = c(0, 0.1),labels = treat_labels) +
    scale_x_continuous(breaks = c(0,0.5,1), labels =  c(0,0.5,1)) +
    scale_fill_manual(values = c(post_step_0 = "#1B7837", post_step_1 = "#D8B58A")) +
    coord_cartesian(xlim = c(0,1)) +
    theme_classic(base_size = 10) +
    theme(aspect.ratio = 3,
          legend.position = "none")
  
  
}

# acc ----
{
  df_dis_acc <- df_pair |> dplyr::select(follow_dis, acc_0, acc_1, treat) |>
    mutate(follow_dis_bin = round(follow_dis, 1)) |>
    pivot_longer(cols = starts_with("acc"))
  
  ggplot(df_dis_acc |> filter(follow_dis_bin < 1.5), 
         aes(x = follow_dis_bin, y = value)) + 
    stat_summary(geom = "ribbon", fun.data = mean_se, alpha = 0.2, aes(fill = name) ) +
    stat_summary(geom = "line", fun = mean, linewidth = 1, aes(col = name)) +
    scale_color_manual(values = c(acc_0 = "#1B7837", acc_1 = "#D8B58A")) +
    scale_fill_manual(values = c(acc_0 = "#1B7837", acc_1 = "#D8B58A")) +
    geom_hline(yintercept = 0)+
    scale_x_continuous(breaks = c(0,0.5,1), labels =  c(0,0.5,1)) +
    scale_y_continuous(breaks = c(-0.03, 0, 0.03), 
                       labels = c(-0.03, 0, 0.03)) +
    coord_cartesian(xlim = c(0, 1.1), ylim = c(-0.035, 0.035))+
    theme_bw(base_size = 11) +
    facet_wrap(~ treat, labeller = labeller(treat = treat_labels)) +
    labs(x = "Distance (BL)", y = "Accerelation (BL/sec2)") +
    theme(strip.placement = "outside",
          strip.background = element_blank(),
          legend.position = "none",
          legend.title = element_blank(),
          aspect.ratio = 3/4
          )
  
  ggsave("output/accerelation.pdf", device = cairo_pdf, family = "Arial",
         width = 5, height = 4)
}
