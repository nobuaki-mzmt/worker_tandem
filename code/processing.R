# ------------------------------------------------------------------------------
# packages
# ------------------------------------------------------------------------------
{
  library(arrow)
  library(stringr)
  library(data.table)
  library(dplyr)
  library(tidyr)
  
  library(ggplot2)
  library(viridis)
  library(ggridges)
  library(patchwork)
  
  scale_factor <- 112/1440 #mm/pix
}
# ------------------------------------------------------------------------------
library(survminer)    
library(survival)

library(scales)

# data prep ----
{
  df_FM <-     arrow::read_feather("data_fmt/trajectory/FM_df.feather")
  df_alates <- arrow::read_feather("data_fmt/trajectory/alates_df.feather")
  df_worker <- arrow::read_feather("data_fmt/trajectory/worker_df.feather")
  df_soldier <- arrow::read_feather("data_fmt/trajectory/soldier_df.feather")
  
  pattern <- "^.*\\d+-\\d+"
  df_FM <- df_FM %>% mutate(video = str_extract(video, pattern))
  df_alates <- df_alates %>% mutate(video = str_extract(video, pattern))
  df_worker <- df_worker %>% mutate(video = str_extract(video, pattern))
  df_soldier <- df_soldier %>% mutate(video = str_extract(video, pattern))
  
  # skip_list
  df_alates <- df_alates %>% filter(!(video == "Ret_ama_FW_F_1-6" & ind_id %in% c(4, 5)))
  df_alates <- df_alates %>% filter(!(video == "Ret_ama_FW_H_7-12" & ind_id %in% c(4)))
  df_alates <- df_alates %>% filter(!(video == "Ret_ama_FW_I_1-6" & ind_id %in% c(1,4)))
  df_alates <- df_alates %>% filter(!(video == "Ret_ama_FW_I_7-12" & ind_id %in% c(1))) 
  df_alates <- df_alates %>% filter(!(video == "Ret_ama_MW_G_7-12" & ind_id %in% c(4,5)))
  df_alates <- df_alates %>% filter(!(video == "Ret_ama_SM_I_7-12"& ind_id %in% c(0, 3)))
  
  df_worker <- df_worker %>% filter(!(video == "Ret_ama_FW_F_1-6" & ind_id %in% c(4, 5)))
  df_worker <- df_worker %>% filter(!(video == "Ret_ama_FW_G_1-6" & ind_id %in% c(1,2,4,5)))
  df_worker <- df_worker %>% filter(!(video == "Ret_ama_FW_H_7-12" & ind_id %in% c(4)))
  df_worker <- df_worker %>% filter(!(video == "Ret_ama_FW_I_1-6" & ind_id %in% c(1,4)))
  df_worker <- df_worker %>% filter(!(video == "Ret_ama_FW_I_7-12" & ind_id %in% c(1))) 
  df_worker <- df_worker %>% filter(!(video == "Ret_ama_MW_G_7-12" & ind_id %in% c(4,5)))
  df_soldier <- df_soldier %>% filter(!(video == "Ret_ama_SM_I_7-12"& ind_id %in% c(0, 3)))
}

# plot all trajectories ----
{
  plot_traj <- function(df, ...){
    ggplot(df, aes(x = x_body, y = y_body, col = as.factor(ind_id) ))+
      scale_color_viridis(discrete = T, option = "D")+
      geom_path(alpha = 1)+
      coord_cartesian(xlim = c(0, 2100), ylim=c(0, 1400)) +
      scale_y_reverse() +
      facet_wrap(~video)+
      theme_classic()+
      theme(aspect.ratio = 2/3, legend.position = "none")+
      labs(...)
  }
  
  save_comparison_plot <- function(df1, df2, video_name) {
    p1 <- plot_traj(df1, title = video_name)
    p2 <- plot_traj(df2)
    ggsave(
      filename = file.path("output/trajectory", paste0(video_name, ".png")),
      plot = p1 + p2, width = 6, height = 4
    )
  }
  
  # FM
  video_list <- unique(df_FM$video)
  for(i in 1:length(video_list)){
    save_comparison_plot(
      df1 = df_FM %>% filter(video == video_list[i] & ind_id %% 2 == 0),
      df2 = df_FM %>% filter(video == video_list[i] & ind_id %% 2 == 1),
      video_name = video_list[i])
  }
  
  # alate-worker
  video_list_a <- unique(df_alates$video)
  video_list_w <- unique(df_worker$video)
  video_list <- video_list_a[video_list_a %in% video_list_w]
  for(i in 1:length(video_list)){
    save_comparison_plot(
      df1 = df_alates %>% filter(video == video_list[i]),
      df2 = df_worker %>% filter(video == video_list[i]),
      video_name = video_list[i])
  }
  
  # alate-soldier
  video_list_s <- unique(df_soldier$video)
  video_list <- video_list_a[video_list_a %in% video_list_s]
  for(i in 1:length(video_list)){
    save_comparison_plot(
      df1 = df_alates %>% filter(video == video_list[i]),
      df2 = df_soldier %>% filter(video == video_list[i]),
      video_name = video_list[i])
  }
}

# data pairing ----
{
  df_FM <- df_FM |> mutate(
    colony = str_split_i(video, "_", 4),
    treat  = str_split_i(video, "_", 3),
    well_id = ind_id %/% 2,
    role = if_else( ind_id %% 2 == 0, "female", "male" ),
    role_id = if_else( ind_id %% 2 == 0, 0, 1)
  )
  
  df_alates <- df_alates |> group_by(video) |>
    mutate(
    colony = str_split_i(video, "_", 4),
    treat  = str_split_i(video, "_", 3),
    well_id = dense_rank(ind_id),
    role = if_else( treat == "FW", "female", "male"),
    role_id = 1
  )
  
  df_worker <- df_worker |> group_by(video) |>
    mutate(
      colony = str_split_i(video, "_", 4),
      treat  = str_split_i(video, "_", 3),
      well_id = dense_rank(ind_id),
      role = "worker", role_id = 0
    )
  
  df_soldier <- df_soldier |> group_by(video) |>
    mutate(
      colony = str_split_i(video, "_", 4),
      treat  = str_split_i(video, "_", 3),
      well_id = dense_rank(ind_id),
      role = "soldier", role_id = 0
    )
  
  df_all <- bind_rows(df_FM, df_alates, df_worker, df_soldier)
  df_all <- df_all |> mutate(ind_name = paste(video, well_id, role, sep = "_"))

  euclid_dis <- function(x0, y0, x1, y1){
    sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1))
  }
  
  df_body <- df_all |> mutate(body_length = euclid_dis(x_head, y_head, x_body, y_body) +
                     euclid_dis(x_tip, y_tip, x_body, y_body)) |>
    group_by(video, well_id, role, ind_name) |>
    summarise(body_length = mean(body_length, na.rm = T), .groups = "drop")
  
  
}
save(df_all, df_body, file = "data_fmt/df_all.rda")
load("data_fmt/df_all.rda")

df_pair <- df_all |> pivot_wider(
  id_cols = c(video, well_id, time_sec, colony, treat), 
  names_from = role_id,
  values_from = starts_with("x_") | starts_with("y_")
)
df_body <- df_body |> group_by(video, well_id) |> 
  summarise(body_length = mean(body_length))

df_pair <- df_pair |> left_join(df_body, by = c("video", "well_id")) |>
  mutate(across(starts_with("x_"), ~ .x / body_length),
         across(starts_with("y_"), ~ .x / body_length),
         partner_dis = euclid_dis(x_body_0, y_body_0, x_body_1, y_body_1),
         follow_dis = euclid_dis(x_tip_0, y_tip_0, x_head_1, y_head_1),
         lead_dis   = euclid_dis(x_tip_1, y_tip_1, x_head_0, y_head_0),
         ang_0 = atan2(y_head_0 - y_tip_0, x_head_0 - x_tip_0),
         ang_1 = atan2(y_head_1 - y_tip_1, x_head_1 - x_tip_1),
         ang_1to0 = atan2(y_body_0 - y_body_1, x_body_0 - x_body_1),
         ang_0to1 = ang_1to0 + pi,
         dir_1to0 = ang_1to0 - ang_1,
         dir_0to1 = ang_0to1 - ang_0
  ) 


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
  ggsave(p1, file = sprintf("output/relative_pos/%s.pdf", i_v), width = 4, height = 3)
}

# ----

df_tandem_time <- df_pair |> dplyr::select(video, well_id, time_sec, colony, treat,
                                           partner_dis, follow_dis, lead_dis) |>
  mutate(tandem = follow_dis < 1  & lead_dis > 1,
         time_bin = round(time_sec)) |>
  group_by(time_bin, treat) |>
  summarize(tandem_prop = mean(tandem)) |> ungroup()

ggplot(df_tandem_time, aes(x = time_bin, y = tandem_prop, 
                           group = treat, col = treat)) +
  geom_line(alpha=0.5)


gap_max <- 10  # number of frames allowed as gap (e.g., 1 sec)

min_len <- 10  # minimum frames to count as tandem (e.g., 1 sec)

df_tandem <- df_pair|> dplyr::select(video, well_id, time_sec, colony, treat,
                                     partner_dis, follow_dis, lead_dis) |>
  group_by(video, well_id) |>
  arrange(time_sec, .by_group = TRUE) |>
  mutate(
    #tandem_raw = partner_dis < 2,
    tandem_raw = follow_dis < 1 & lead_dis > 1,
    
    # identify runs
    run_id = rleid(tandem_raw),
    
    # run length
    run_len = ave(tandem_raw, run_id, FUN = length)
  ) |>
  
  # keep only sufficiently long tandem runs
  mutate(
    tandem = if_else(tandem_raw & run_len >= min_len, TRUE, FALSE)
  ) |>
  
  # define clusters
  mutate(
    cluster_start = tandem & lag(!tandem, default = TRUE),
    cluster_idx = cumsum(cluster_start),
    cluster_idx = if_else(tandem, cluster_idx, NA_integer_),
    pair_event = if_else(
      !is.na(cluster_idx),
      sprintf("%s_%d_%02d", video, well_id, cluster_idx),
      NA_character_
    )
  ) |>
  select(-cluster_start, -cluster_idx, -run_id, -run_len, -tandem_raw) |>
  ungroup()

df_pair_dur <- df_tandem |> filter(time_sec < 1800.1) |>
  filter(!is.na(pair_event)) |>
  group_by(pair_event, treat, well_id, video) |>
  summarise(
    start_time = first(time_sec),
    end_time = last(time_sec),
    duration = (n() * 0.2), 
    cens = !(dplyr::near(last(time_sec), 1800) | dplyr::near(first(time_sec), 0)),
    .groups = "drop"
  )

df_plot <- df_pair_dur |> filter(duration > 0) |>
  filter(treat != "FW")

ggsurvplot(
  survfit(Surv(duration, cens) ~ treat, data = df_plot),
  data = df_plot,
  censor = FALSE,
  conf.int.style = "ribbon",
  conf.int.fill = TRUE,
  ggtheme = theme_classic(),
  
)$plot + 
  scale_x_continuous(trans = "pseudo_log", breaks = c(0, 1, 10, 100, 1000))

#cumulative hazard plot
ggsurvplot(
  survfit(Surv(duration, cens) ~ treat, data = df_plot),
  data = df_plot,
  fun = "cumhaz",
  censor = TRUE,
  ggtheme = theme_classic()
)$plot

df_plot |> ggplot(aes(x = start_time, y = 1*(duration > 5))) +
  stat_smooth(method = glm, method.args = list(family = "binomial")) + 
  facet_wrap(~treat)

# hazard is proportional? 
library(coxme)
fit_cox <- coxme(Surv(duration, cens) ~ treat + (1|video/well_id), data = df_plot)
zph <- cox.zph(fit_cox)

plot(zph, resid = TRUE, se = TRUE, col = "steelblue")
abline(h = 0, lty = 2, col = "red")


ggsurvplot(
  survfit(Surv(duration, cens) ~ treat, data = df_plot),
  data = df_plot,
  fun = "cumhaz",
  censor = TRUE,
  ggtheme = theme_classic()
)$plot 


## speed ----

df_dis <- df_pair |> filter(treat != "FW") |>
  group_by(video, well_id) |>  mutate(
  pre_step_0  = euclid_dis(x_body_0, y_body_0, lag(x_body_0),  lag(y_body_0)),
  pre_step_1  = euclid_dis(x_body_1, y_body_1, lag(x_body_1),  lag(y_body_1)),
  post_step_0 = euclid_dis(x_body_0, y_body_0, lead(x_body_0), lead(y_body_0)),
  post_step_1 = euclid_dis(x_body_1, y_body_1, lead(x_body_1), lead(y_body_1)),
  pre_step_head_0 = euclid_dis(x_head_0, y_head_0, lag(x_head_0), lag(y_head_0)),
  pre_step_head_1 = euclid_dis(x_head_1, y_head_1, lag(x_head_1), lag(y_head_1))
)

df_dis |>
  ggplot(aes(x = partner_dis)) + 
  stat_smooth(aes(y = pre_step_0), col = "red") +
  stat_smooth(aes(y = pre_step_1), col = "blue") +
  facet_wrap(~treat)


library(data.table)

df_pair_dur <- tibble(df_pair_dur)
df_check <- df_pair_dur |> filter(cens) |> filter(duration < 10)

setDT(df_dis)
setDT(df_pair_dur)

df_dis_split <- split(df_dis, list(df_dis$video, df_dis$well_id), drop = TRUE)

list_df <- vector("list", nrow(df_check))

for(i in seq_len(nrow(df_check))){
  key <- paste(df_check$video[i], df_check$well_id[i], sep = ".")
  df_sub <- df_dis_split[[key]]
  if (is.null(df_sub)) next
  st <- df_check$start_time[i]
  list_df[[i]] <- df_sub[
    df_sub$time_sec > st - 5 & df_sub$time_sec < st + 5, 
  ][, `:=`(
    start_time = st,
    pair_event = df_check$pair_event[i],
    rel_time = time_sec - st
  )]
}

df_out <- rbindlist(list_df)

df_out <- tibble(df_out)

df_out |> ggplot(aes(x = rel_time)) + 
  stat_smooth(aes(y = pre_step_0), col = "red") +
  stat_smooth(aes(y = pre_step_1), col = "blue") +
  facet_wrap(~treat)
