
df = as.data.frame(y)

df = df %>% 
  rename(`Treament Series` = c(1),
         `Donor 1` = c(2),
         `Donor 2` = c(3),
         `Donor 3` = c(4),
         `Donor 4` = c(5)) %>% 
  mutate(Time = 1:30) %>% 
  gather(Series, Y, -(Time))

color_codes <- viridis::viridis_pal(option = "mako")(10)

# Create a plot to visualize the colors
plot(1, type = "n", xlim = c(0, 1), ylim = c(0, length(color_codes) + 1),
     xlab = "", ylab = "", axes = "none")
for (i in seq_along(color_codes)) {
  rect(0.2, i, 0.8, i + 0.8, col = color_codes[i], border = NA)
  text(0.9, i + 0.4, color_codes[i], adj = c(0, 0.5))
}

color_codes <- viridis::viridis_pal(option = "mako")(10)
custom_colors <- c(color_codes[5], color_codes[6], color_codes[8], 
                   color_codes[9], "black")

p1 = ggplot(df) +
  aes(x = Time, y = Y, colour = Series, linetype = Series) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = custom_colors) +
  scale_linetype_manual(values = c("solid", "solid", "solid", "solid", "dotted"))+
  geom_vline(xintercept = 20, linetype = "dashed", color = "red", linewidth = .8)+
  theme_minimal()+
  theme(
    axis.title.y = element_text(size = 12L),
    axis.title.x = element_text(size = 12L),
    plot.caption = element_text(size = 10L))



png('~/Diss/Topics/Synthetic Control/Latex/Paper/images/p1.png', width = 7, height = 3.5, units = 'in', res = 1000)
p1
dev.off()
