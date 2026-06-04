df <- read.csv("data_raw/master_well.csv")
df_stat <- df |> filter(!is.na(tandem_before_attack )) |> dplyr::select(tandem_before_attack, severe_attack_sec, fetal_attack_sec, tandem_after_attack)
dim(df |> filter(treat == "MW2"))

r <- glm(
  tandem_before_attack ~ log(severe_attack_sec + 1),
  family = binomial,
  data = df_stat
)
Anova(r)
