install.packages('tidysynth')

require(tidysynth)
data("smoking")
smoking %>% dplyr::glimpse()

smoking_out = smoking %>%
  synthetic_control(outcome = cigsale, 
                    unit = state, 
                    time = year, 
                    i_unit = "California", 
                    i_time = 1988, 
                    generate_placebos=T) %>%
  generate_predictor(time_window = 1980:1988,
                     ln_income = mean(lnincome, na.rm = T),
                     ret_price = mean(retprice, na.rm = T),
                     youth = mean(age15to24, na.rm = T)) %>%
  generate_predictor(time_window = 1984:1988,
                     beer_sales = mean(beer, na.rm = T)) %>%
  generate_predictor(time_window = 1975,
                     cigsale_1975 = cigsale) %>%
  generate_predictor(time_window = 1980,
                     cigsale_1980 = cigsale) %>%
  generate_predictor(time_window = 1988,
                     cigsale_1988 = cigsale) %>%
  generate_weights(optimization_window = 1970:1988, 
                   margin_ipop = .02,sigf_ipop = 7,bound_ipop = 6) %>%
  generate_control()

smoking_out %>% plot_trends()

smoking_out %>% plot_differences()

smoking_out %>% plot_weights()

smoking_out %>% grab_balance_table()

smoking_out %>% plot_placebos()

smoking_out %>% plot_placebos(prune = FALSE)

smoking_out %>% plot_mspe_ratio()

smoking_out %>% grab_significance()

smoking_out

smoking_out %>% grab_synthetic_control(placebo = T)

smoking_out %>% 
  tidyr::unnest(cols = c(.outcome)) 
