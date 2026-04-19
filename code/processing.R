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

# data prep
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

# plot all trajectories
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

# data pairing
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
  
  save(df_all, df_body, file = "data_fmt/df_all.rda")
}

df_pair <- df_all |> pivot_wider(
  id_cols = c(video, well_id, time_sec, colony, treat), 
  names_from = role_id,
  values_from = starts_with("x_") | starts_with("y_")
)

df_all |> filter(video == "Ret_ama_SM_F_1-4", time_sec == 0.2) |> View()




# ------------------------------------------------------------------------------
# compute relative positioning of termites and termitophiles with kinetic parameters
# ------------------------------------------------------------------------------
Plot = F
{
  df_termite <- arrow::read_feather("data_fmt/lon_tra/termite_df.feather")
  meta_df <- data.frame(fread("individual_list.csv"))[,1:6]
  
  f.namesplace <- list.files("data_fmt/lon_tra", 
                             pattern=".feather", full.names=T)
  f.names <- basename(list.files("data_fmt/lon_tra", 
                                 pattern=".feather", full.names=F))
  
  list_df_plot_sum <- NULL
  list_df_interact <- NULL
  list_i <- 1
  
  for(i_s in 1:(length(f.names)-1)){
    
    species <- str_replace(f.names[i_s], "_df.feather", "")
    data_name <- f.namesplace[i_s]
    df <- arrow::read_feather(data_name)
    
    print(species)
    
    #df$video <- str_replace(df$video, "_[^_]+$", "")
    df$video <- sub("_([^_]+)$", "", df$video)
    
    video_names <- unique(df$video)
    
    for(i_video in video_names){
      v_name <- i_video
      print(i_video)
      meta_df_temp <- subset(meta_df, video == i_video & type == species)
      dftemp <- df %>% filter(video == i_video)
      df_termitetemp <- df_termite %>% filter(video == i_video)
      
      # exceptions TODO
      if(i_video == "Lon_lon_NM23074_sf1-w1_01-06-1"){
        v_name <- "Lon_lon_NM23074_sf1-w1_01-06"
        meta_df_temp <- subset(meta_df, video == v_name)
        dftemp1 <- subset(df, video == i_video)
        dftemp2 <- subset(df, video == "Lon_lon_NM23074_sf1-w1_01-06-2")
        dftemp <- rbind(dftemp1, dftemp2)
        
        df_termitetemp1 <- subset(df_termite, video == i_video)
        df_termitetemp2 <- subset(df_termite, 
                                  video == "Lon_lon_NM23074_sf1-w1_01-06-2")
        df_termitetemp <- rbind(df_termitetemp1, df_termitetemp2)
      }
      if(i_video == "Lon_lon_NM23074_sf1-w1_01-06-2"){ next; }
      if(i_video == "Lon_lon_NM23086_sf1-w1_01-03-1"){
        v_name <- "Lon_lon_NM23086_sf1-w1_01-03"
        meta_df_temp <- subset(meta_df, video == v_name)
        dftemp1 <- subset(df, video == i_video)
        dftemp2 <- subset(df, video == "Lon_lon_NM23086_sf1-w1_01-03-2")
        dftemp <- rbind(dftemp1, dftemp2)
        dftemp$fill[dftemp$fill == 2] <- 3
        dftemp$fill[dftemp$fill == 1] <- 2
        
        df_termitetemp1 <- subset(df_termite, video == i_video)
        df_termitetemp2 <- subset(df_termite, 
                                  video == "Lon_lon_NM23086_sf1-w1_01-03-2")
        df_termitetemp <- rbind(df_termitetemp1, df_termitetemp2)
        
        
      }
      if(i_video == "Lon_lon_NM23086_sf1-w1_01-03-2"){ next; }
      
      if(length(unique(dftemp$fill)) != length(meta_df_temp$well)){
        print("something wrong")
        browser()
      }
      sym_fill_list <- unique(dftemp$fill)
      
      for(i in 1:length(sym_fill_list)){
        sym_fill  <- sym_fill_list[i]
        real_fill <- meta_df_temp[i, "well"]
        if( meta_df_temp[i,]$analyze == 0){ next; }
        
        df_m <- dftemp %>% filter(fill == sym_fill) %>%
          mutate(index = row_number()) %>%
          filter(index%%6 == 0)
        
        df_t <- df_termitetemp %>% filter(fill == real_fill) %>%
          mutate(index = row_number()) %>%
          filter(index%%6 == 0)
        
        # allocate other termite (next) in the same well
        next_well <- meta_df %>% filter(video == v_name) 
        x <- next_well$well == (real_fill)
        next_well <- next_well$well[x[c(length(x), 1:(length(x)-1))]]
        if(i_video == "Lon_lon_NM23087_termitophile1-w1_01"){
          ## this video only record one cell, so take it from other vieo
          df_ran <- df_termitetemp <- df_termite %>% 
            filter(video == video_names[1]) %>% mutate(index = row_number())
        } else {
          df_ran <- df_termitetemp %>% filter(fill == next_well) %>%
            mutate(index = row_number()) %>%
            filter(index%%6 == 0)
        }
        
        if(dim(df_ran)[1] > dim(df_m)[1]){
          df_ran <- df_ran[1:dim(df_m)[1],]
        }
        
        df_t <- df_t %>% mutate(across(matches("^(x|y)_"), ~ .x * scale_factor))
        df_ran <- df_ran %>% mutate(across(matches("^(x|y)_"), ~ .x * scale_factor))
        df_m <- df_m %>% mutate(across(c(x, y), ~ .x * scale_factor))
        
        dir_vec_x = df_t$x_head - df_t$x_tip
        dir_vec_y = df_t$y_head - df_t$y_tip
        tbl <- mean(sqrt( dir_vec_x^2 + dir_vec_y^2))
        
        calc_rel_geometry <- function(mx, my, tx, ty, dir_x, dir_y) {
          rel_x <- mx - tx
          rel_y <- my - ty
          list(
            ang = atan2(rel_y, rel_x) - atan2(dir_y, dir_x) + pi/2,
            dis = sqrt(rel_x^2 + rel_y^2)
          )
        }
        
        # relative to body
        body_geo <- calc_rel_geometry(df_m$x, df_m$y, df_t$x_body, df_t$y_body,
                                      dir_vec_x, dir_vec_y)
        
        # relative to tail
        tail_geo <- calc_rel_geometry(df_m$x, df_m$y, df_t$x_tip, df_t$y_tip,
                                      dir_vec_x, dir_vec_y)
        
        # relative to head
        head_geo <- calc_rel_geometry(df_m$x, df_m$y, df_t$x_tip, df_t$y_tip,
                                      dir_vec_x, dir_vec_y)
        
        # random data
        df_ran$x_body <- df_ran$x_body + mean(df_m$x) - mean(df_ran$x_body)
        df_ran$y_body <- df_ran$y_body + mean(df_m$y) - mean(df_ran$y_body)
        random_geo <- calc_rel_geometry(df_m$x, df_m$y,
                                        df_ran$x_body, df_ran$y_body,
                                        dir_vec_x, dir_vec_y)
        
        
        m_moved_dis <- c(NA, sqrt( diff(df_m$x)^2 + diff(df_m$y)^2))
        t_moved_dis <- c(NA, sqrt( diff(df_t$x_body)^2 + diff(df_t$y_body)^2))
        
        m_acc <- c(diff(m_moved_dis), NA)
        t_acc <- c(diff(t_moved_dis), NA)
        
        if(Plot){
          df_plot_temp <- data.frame(
            frame = df_m$index,
            t_moved_dis,
            m_moved_dis
          )
          
          # check if termites and termitophile in the same well by ploting 
          g1 <- ggplot()+
            geom_path(data = subset(df_m, index < 18000), 
                      aes(x=x, y=y), col=1, alpha = 0.5) +
            geom_path(data = subset(df_t, index < 18000), 
                      aes(x=x_body, y=y_body), col=2, alpha = 0.5) +
            ggtitle(paste(i_video, real_fill)) +
            theme(aspect.ratio = 1, plot.title = element_text(size = 7))
          g2 <- ggplot()+
            geom_path(data = subset(df_m, index > 18000 & index < 36000), 
                      aes(x=x, y=y), col=1, alpha = 0.5) +
            geom_path(data = subset(df_t, index > 18000 & index < 36000), 
                      aes(x=x_body, y=y_body), col=2, alpha = 0.5) +
            ggtitle(paste(i_video, real_fill)) +
            theme(aspect.ratio = 1, plot.title = element_text(size = 7))
          g3 <- ggplot()+
            geom_path(data = subset(df_m, index > 36000), 
                      aes(x=x, y=y), col=1, alpha = 0.5) +
            geom_path(data = subset(df_t, index > 36000), 
                      aes(x=x_body, y=y_body), col=2, alpha = 0.5) +
            ggtitle(paste(i_video, real_fill)) +
            theme(aspect.ratio = 1, plot.title = element_text(size = 7))
          g4 <- ggplot(data = df_plot_temp) +
            geom_density(aes(x = m_moved_dis), col = 1) +
            geom_density(aes(x = t_moved_dis), col = 2) +
            xlim(c(0,5))
          g1+g2+g3+g4
          ggsave(paste0("output/trajectory/", species, "_", i_video, "_", real_fill, ".png"),
                 width = 6, height = 4)
          
          
        }
        
        scale_cols <- c(
          "x","y",
          "dis","tail_dis","head_dis",
          "m_moved_dis","t_moved_dis",
          "m_acc","t_acc"
        )
        
        df_plot = data.frame(
          species,
          x = body_geo$dis * cos(body_geo$ang)/tbl,
          y = body_geo$dis * sin(body_geo$ang)/tbl,
          dis = body_geo$dis/tbl,
          tail_dis = tail_geo$dis/tbl,
          head_dis = head_geo$dis/tbl,
          ran_dis = random_geo$dis/tbl,
          m_moved_dis = m_moved_dis/tbl,
          t_moved_dis = t_moved_dis/tbl,
          m_acc = m_acc/tbl,
          t_acc = t_acc/tbl,
          video = i_video,
          fill  = real_fill
        )
        list_df_plot_sum[[list_i]] <- df_plot
        
        # interaction
        interact <- body_geo$dis < tbl * 1.5
        interact.end <- which(interact)[c(diff(which(interact))>1,T)]
        interact.sta <- which(interact)[c(T, diff(which(interact))>1)]
        interact_duration <- (interact.end - interact.sta)
        
        df_interact_temp <- data.frame(
          species,
          video = i_video,
          fill = real_fill,
          interact_duration,
          cens = (interact.end != length(interact)) & (interact.sta != 0)
        )
        list_df_interact[[list_i]] <- df_interact_temp
        
        
        list_i <- list_i + 1
      }
    }
  }
  
  df_plot_sum <- rbindlist(list_df_plot_sum)
  df_interact <- rbindlist(list_df_interact)
  
  
  df_plot_sum <- df_plot_sum %>%
    mutate(across(all_of(scale_cols), round, 4))
  
  save(df_plot_sum, file = "data_fmt/df_rel_pos.rda")
  save(df_interact, file = "data_fmt/df_interact.rda")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# preprocessing tandem data in R. amamianu
# ------------------------------------------------------------------------------
Plot= T
{
  df_tandem <- arrow::read_feather("data_fmt/tandem_df.feather")
  
  df_tandem$sex <- df_tandem$fill%%2
  df_tandem$fill <- floor(df_tandem$fill/2)
  
  df_tandem[,4:9] <- df_tandem[,4:9] * scale_factor
  colnames(df_tandem)
  
  df_tandem <- subset(df_tandem, index%%6 == 0) # 5FPS
  
  video_list <- unique(df_tandem$video)
  df_res <- NULL
  for(i_v in 1:length(video_list)){
    print(video_list[i_v])
    for(i_f in 0:5){
      df_temp <- subset(df_tandem, fill == i_f & video == video_list[i_v])
      df_temp_f <- subset(df_temp, sex == 0)
      df_temp_m <- subset(df_temp, sex == 1)
      
      fbl <- sqrt( (df_temp_f$x_head - df_temp_f$x_body)^2 +
                     (df_temp_f$y_head - df_temp_f$y_body)^2 ) +
        sqrt( (df_temp_f$x_body - df_temp_f$x_tip)^2 +
                (df_temp_f$y_body - df_temp_f$y_tip)^2 ) 
      
      mbl <- sqrt( (df_temp_m$x_head - df_temp_m$x_body)^2 +
                     (df_temp_m$y_head - df_temp_m$y_body)^2 ) +
        sqrt( (df_temp_m$x_body - df_temp_m$x_tip)^2 +
                (df_temp_m$y_body - df_temp_m$y_tip)^2 ) 
      
      tbl <- (mean(fbl) + mean(mbl))/2
      
      dis <- sqrt( (df_temp_f$x_body - df_temp_m$x_body)^2 +
                     (df_temp_f$y_body - df_temp_m$y_body)^2 )
      tail_dis <- sqrt( (df_temp_f$x_tip - df_temp_m$x_body)^2 +
                          (df_temp_f$y_tip - df_temp_m$y_body)^2 )
      head_dis <- sqrt( (df_temp_f$x_head - df_temp_m$x_body)^2 +
                          (df_temp_f$y_head - df_temp_m$y_body)^2 )
      
      f_step <- c(NA, sqrt((diff(df_temp_f$x_body))^2 + (diff(df_temp_f$y_body))^2))
      m_step <- c(NA, sqrt((diff(df_temp_m$x_body))^2 + (diff(df_temp_m$y_body))^2))
      
      f_acc  <- c(diff(f_step), NA)
      m_acc  <- c(diff(m_step), NA)
      
      dir_vec_x = df_temp_f$x_head - df_temp_f$x_tip
      dir_vec_y = df_temp_f$y_head - df_temp_f$y_tip
      
      rel_vec_x = df_temp_m$x_body - df_temp_f$x_body
      rel_vec_y = df_temp_m$y_body - df_temp_f$y_body
      
      sta_ang <- atan2(rel_vec_y, rel_vec_x) - atan2(dir_vec_y, dir_vec_x) + pi/2
      dis <- sqrt(rel_vec_x^2 + rel_vec_y^2)
      
      df_temp <- data.frame(
        video = i_v,
        fill = i_f,
        x = dis * cos(sta_ang)/tbl,
        y = dis * sin(sta_ang)/tbl,
        dis = dis/tbl,
        tail_dis = tail_dis/tbl,
        head_dis = head_dis/tbl,
        f_step = f_step/tbl,
        m_step = m_step/tbl,
        f_acc = f_acc/tbl,
        m_acc = m_acc/tbl
      )
      
      df_res <- rbind(df_res, df_temp)
      
      if(Plot){
        df_plot_temp <- data.frame(
          frame = df_temp_f$index, f_step, m_step
        )
        
        # check if termites and termitophile in the same well by ploting 
        g1 <- ggplot()+
          geom_path(data = subset(df_temp_f, index < 18000), 
                    aes(x=x_body, y=y_body), col=1, alpha = 0.5) +
          geom_path(data = subset(df_temp_m, index < 18000), 
                    aes(x=x_body, y=y_body), col=2, alpha = 0.5) +
          ggtitle(paste(video_list[i_v], i_f)) +
          theme(aspect.ratio = 1, plot.title = element_text(size = 7))
        g2 <- ggplot()+
          geom_path(data = subset(df_temp_f, index > 18000 & index < 36000), 
                    aes(x=x_body, y=y_body), col=1, alpha = 0.5) +
          geom_path(data = subset(df_temp_m, index > 18000 & index < 36000), 
                    aes(x=x_body, y=y_body), col=2, alpha = 0.5) +
          ggtitle(paste(video_list[i_v], i_f)) +
          theme(aspect.ratio = 1, plot.title = element_text(size = 7))
        g3 <- ggplot()+
          geom_path(data = subset(df_temp_f, index > 36000), 
                    aes(x=x_body, y=y_body), col=1, alpha = 0.5) +
          geom_path(data = subset(df_temp_m, index > 36000), 
                    aes(x=x_body, y=y_body), col=2, alpha = 0.5) +
          ggtitle(paste(video_list[i_v], i_f)) +
          theme(aspect.ratio = 1, plot.title = element_text(size = 7))
        g4 <- ggplot(data = df_plot_temp) +
          geom_density(aes(x = f_step), col = 1) +
          geom_density(aes(x = m_step), col = 2) +
          xlim(c(0,5))
        g1+g2+g3+g4
        ggsave(paste0("output/trajectory/", "tandem_", 
                      video_list[i_v], "_", i_f, ".png"),
               width= 6, height = 4)
      }
    }
  }
  save(df_res, file = "data_fmt/tandem.rda")
}
# ------------------------------------------------------------------------------






