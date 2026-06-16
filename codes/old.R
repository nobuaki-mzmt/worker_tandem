## analysis performed on onset of separation events.
## after spending hours, I decided not to go this direction.
## this is a code record:


df_cleaned <- df_pair |> 
  group_by(pair_id) |> 
  dplyr::select(pair_id, time_sec, treat, tandem, post_step_0, post_step_1, acc_0, acc_1, follow_dis, partner_dis, dir_0to1) |>
  arrange(time_sec, .by_group = TRUE) |>
  mutate(run_id = rleid(tandem)) |> 
  group_by(pair_id, run_id) |> 
  mutate(
    run_duration = max(time_sec) - min(time_sec),
    tandem = if_else(tandem & run_duration < 1.5, FALSE, tandem)
  ) |> 
  group_by(pair_id) |> 
  mutate(
    run_id = rleid(tandem),
    is_offset = !tandem & lag(tandem, default = FALSE)
  ) |> 
  ungroup()


df_indexed <- df_cleaned |>
  group_by(pair_id) |>
  mutate(row_idx = row_number()) |>
  ungroup()

# Get the row index of each offset event
offset_idx <- df_indexed |>
  filter(is_offset) |>
  group_by(pair_id) |>
  mutate(offset_event_id = row_number()) |>
  ungroup() |>
  dplyr::select(pair_id, offset_row_idx = row_idx, offset_event_id)

# For each offset event, expand to +/-20 rows, then join back
window_size <- 50

df_offsets <- offset_idx |>
  reframe(
    row_idx       = (offset_row_idx - window_size):(offset_row_idx + window_size),
    relative_step = -window_size:window_size,
    offset_event_id = paste(pair_id, offset_event_id, sep ="_"),   # carry through
    .by = c(pair_id, offset_event_id)
  ) |>
  left_join(df_indexed, by = c("pair_id", "row_idx")) |>
  filter(!is.na(time_sec))


df_offsets <- df_offsets |> filter( ((relative_step < 0) & tandem) | ((relative_step >= 0) & !tandem) ) 

df_plot_offset <- df_offsets |>
  #pivot_longer(cols = c(post_step_0, post_step_1), names_to = "step_type", values_to = "step_length") |>
  pivot_longer(cols = c(acc_0, acc_1), names_to = "step_type", values_to = "step_length") |>
  group_by(relative_step, step_type, treat) |>
  summarise(mean_step = mean(step_length, na.rm = TRUE), 
            mean_partner_dis = mean(follow_dis, na.rm = T), 
            .groups = 'drop')

ggplot(df_plot_offset, aes(x = relative_step, y = mean_step, color = step_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_line(size = 1) +
  facet_wrap(~treat) +
  labs(
    title = "Step Length Dynamics Heading Into Tandem Run",
    x = "Relative Steps (0 = Tandem Starts)",
    y = "Mean Step Length"
  ) +
  theme_minimal()+
  geom_hline(yintercept = 0)


df_offsets |> 
  summarise(sep_duration = max(relative_step) / 5, .by = c(offset_event_id, treat)) |>
  ggplot(aes(x = sep_duration, y = treat)) +
  geom_density_ridges(stat = "binline", binwidth = 0.2) +
  xlim(c(0,5))


df_offsets |> dplyr::select(tandem, partner_dis, offset_event_id, run_duration, is_offset,
                            "relative_step") |> View()




# longer than 2 seconds
df_offset_long_separation <- df_offsets |>
  group_by(offset_event_id) |>
  filter(max(relative_step) > 40) |>
  ungroup()

df_plot_offset <- df_offset_long_separation |>
  #pivot_longer(cols = c(post_step_0, post_step_1), names_to = "step_type", values_to = "step_length") |>
  pivot_longer(cols = c(acc_0, acc_1), names_to = "step_type", values_to = "step_length") |>
  group_by(relative_step, step_type, treat) |>
  summarise(mean_step = mean(step_length, na.rm = TRUE), 
            mean_partner_dis = mean(follow_dis, na.rm = T), 
            .groups = 'drop')

ggplot(df_plot_offset, aes(x = relative_step, y = mean_step, color = step_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(-20,10))+
  facet_wrap(~treat) +
  labs(
    title = "Step Length Dynamics Heading Into Tandem Run",
    x = "Relative Steps (0 = Tandem Starts)",
    y = "Mean Step Length"
  ) +
  theme_minimal()+
  geom_hline(yintercept = 0)







df_analysis <- df_pair |> 
  group_by(pair_id) |> 
  dplyr::select(pair_id, time_sec, treat, tandem, post_step_0, post_step_1, acc_0, acc_1, follow_dis, partner_dis, dir_0to1) |>
  arrange(time_sec, .by_group = TRUE) |>
  mutate(
    run_id = rleid(tandem)
  ) |> 
  group_by(pair_id, run_id) |> 
  mutate(
    run_duration = max(time_sec) - min(time_sec),
    tandem = if_else(tandem & run_duration < 1.5, FALSE, tandem)
  ) |> 
  group_by(pair_id) |> 
  mutate(
    run_id = rleid(tandem), # Recalculate IDs since some TRUE states became FALSE
    is_onset = tandem & lag(!tandem, default = FALSE),
    is_offset = !tandem & lag(tandem, default = FALSE)
  ) |> 
  ungroup()

extract_window <- function(df, flag_col, window_size = 10) {
  event_indices <- which(df[[flag_col]])
  if(length(event_indices) == 0) return(tibble())
  
  lapply(seq_along(event_indices), function(i) {
    idx <- event_indices[i]
    
    start_idx <- idx - window_size
    end_idx <- idx + window_size
    
    if (start_idx < 1 || end_idx > nrow(df)) {
      return(tibble())
    }
    
    df[start_idx:end_idx, ] |> 
      mutate(
        event_instance = i,
        relative_step = (start_idx:end_idx) - idx
      )
  }) |> bind_rows()
}

df_onsets <- df_analysis |> 
  group_split() |> 
  lapply(extract_window, flag_col = "is_onset", window_size = 5) |> 
  bind_rows()

df_onsets <- df_onsets |> filter( ((relative_step < 0) & !tandem) | ((relative_step >= 0) & tandem) )   

df_plot_onset <- df_onsets |>
  pivot_longer(cols = c(post_step_0, post_step_1), names_to = "step_type", values_to = "step_length") |>
  #pivot_longer(cols = c(acc_0, acc_1), names_to = "step_type", values_to = "step_length") |>
  group_by(relative_step, step_type, treat) |>
  summarise(mean_step = mean(step_length, na.rm = TRUE), .groups = 'drop')

ggplot(df_plot_onset, aes(x = relative_step, y = mean_step, color = step_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_line(size = 1) +
  facet_wrap(~treat) +
  labs(
    title = "Step Length Dynamics Heading Into Tandem Run",
    x = "Relative Steps (0 = Tandem Starts)",
    y = "Mean Step Length"
  ) +
  theme_minimal()

df_offsets <- df_analysis |> 
  group_split() |> 
  lapply(extract_window, flag_col = "is_offset", window_size = 20) |> 
  bind_rows()


df_offsets <- df_offsets |> filter( ((relative_step < 0) & tandem) | ((relative_step >= 0) & !tandem) ) |>
  mutate(event_id = paste0(pair_id, event_instance, sep = ""))

df_dis_dist <- df_offsets |> pivot_longer(cols = starts_with("acc_")) |> 
  dplyr::select(name, value, treat, tandem)

ggplot(df_dis_dist) +
  geom_density_ridges(aes(x = value, y = treat, fill= name), stat = "binline", 
                      alpha = 0.5, binwidth = 0.04, scale = 0.85) +
  labs(x = "Step length (BL)", y = "") +
  scale_y_discrete(expand = c(0, 0.1),labels = treat_labels) +
  scale_x_continuous(breaks = c(0,0.5,1), labels =  c(0,0.5,1)) +
  scale_fill_manual(values = c(post_step_0 = "#1B7837", post_step_1 = "#D8B58A")) +
  coord_cartesian(xlim = c(0,1.5)) +
  theme_classic(base_size = 10) +
  theme(aspect.ratio = 3,
        legend.position = "none")


df_plot_offset <- df_offsets |>
  #pivot_longer(cols = c(post_step_0, post_step_1), names_to = "step_type", values_to = "step_length") |>
  pivot_longer(cols = c(acc_0, acc_1), names_to = "step_type", values_to = "step_length") |>
  group_by(relative_step, step_type, treat) |>
  summarise(mean_step = mean(step_length, na.rm = TRUE), 
            mean_partner_dis = mean(follow_dis, na.rm = T), 
            .groups = 'drop')

ggplot(df_plot_offset, aes(x = relative_step, y = mean_step, color = step_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_line(size = 1) +
  #ylim(c(0,0.4))+
  facet_wrap(~treat) +
  labs(
    title = "Step Length Dynamics Heading Into Tandem Run",
    x = "Relative Steps (0 = Tandem Starts)",
    y = "Mean Step Length"
  ) +
  theme_minimal()+
  geom_hline(yintercept = 0)

ggplot(df_plot_offset, aes(x = relative_step, y = mean_partner_dis, color = step_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_line(size = 1) +
  #ylim(c(0,0.4))+
  facet_wrap(~treat) +
  labs(
    title = "Step Length Dynamics Heading Into Tandem Run",
    x = "Relative Steps (0 = Tandem Starts)",
    y = "Mean Step Length"
  ) +
  theme_minimal()+
  geom_hline(yintercept = 0)




df_offsets <- df_analysis |> 
  group_split() |> 
  lapply(extract_window, flag_col = "is_offset", window_size = 20) |> 
  bind_rows() |>
  mutate(event_id = paste(pair_id, event_instance, sep = "_"))

df_offsets <- df_offsets |> filter( ((relative_step < 0) & tandem) | ((relative_step >= 0) & !tandem) ) 

df_offsets |> filter(event_id == "Ret_ama_MW2_G_1-6_4_2259") |> View()
df_analysis |> filter(pair_id == "Ret_ama_MW2_G_1-6_4") |> View()

df_offsets |> filter(event_id == "Ret_ama_MW2_G_1-6_4_2259") |> dplyr::select(tandem, is_offset, relative_step)

df_offsets |> 
  group_by(event_id, treat) |>
  summarise(sep_duration = max(relative_step)) |> 
  ungroup() |> filter(sep_duration < 0)


df_offsets |> 
  group_by(event_id, treat) |>
  summarise(sep_duration = max(relative_step)) |>
  ungroup() |> filter(sep_duration < 10) |>
  ggplot(aes(x = sep_duration, y = treat)) +
  geom_density_ridges(stat = "binline", binwidth = 0.2) 

# longer than 2 seconds
df_offset_long_separation <- df_offsets |>
  group_by(event_id) |>
  filter(max(relative_step) > 8) |>
  ungroup()

#df_offsets |> dplyr::select(tandem, partner_dis, event_id, run_duration, is_offset, is_onset, "relative_step") |> View()


df_plot_offset <- df_offset_long_separation |>
  #pivot_longer(cols = c(post_step_0, post_step_1), names_to = "step_type", values_to = "step_length") |>
  pivot_longer(cols = c(acc_0, acc_1), names_to = "step_type", values_to = "step_length") |>
  group_by(relative_step, step_type, treat) |>
  summarise(mean_step = mean(step_length, na.rm = TRUE), 
            mean_partner_dis = mean(follow_dis, na.rm = T), 
            .groups = 'drop')

ggplot(df_plot_offset, aes(x = relative_step, y = mean_step, color = step_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_line(size = 1) +
  #ylim(c(0,0.4))+
  facet_wrap(~treat) +
  labs(
    title = "Step Length Dynamics Heading Into Tandem Run",
    x = "Relative Steps (0 = Tandem Starts)",
    y = "Mean Step Length"
  ) +
  theme_minimal()+
  geom_hline(yintercept = 0)



df_sep_move <- df_offsets |> filter(relative_step > 0) |> mutate(event_id = paste0(pair_id, event_instance, sep = ""))

df_sep_move_rel <- df_sep_move |>
  mutate(rx = partner_dis * cos(dir_0to1), ry = partner_dis * sin(dir_0to1))

ggplot(df_sep_move_rel,
       aes(x = rx, y = ry, group = event_id, col = relative_step)) +
  geom_path(alpha = 0.25) +
  facet_wrap(~treat)


df_pair |> pull(partner_dis) |> max()

df_pair |> filter(partner_dis > 8) |> pull(pair_id)

