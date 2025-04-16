library(tidyverse)

data <- read.table("/Users/aarondesouza/Desktop/ArabTTAGGG_Three_1M.txt", 
                   header = TRUE,
                   sep = "\t",
                   stringsAsFactors = FALSE) 

data <- data %>%
  mutate(group = as.factor(group), observation = as.factor(observation), as.integer(id), as.integer(value))

data$group <- factor(data$group)

empty_bar <- 3
twolevels <- nlevels(data$observation)
to_add <- data.frame(matrix(NA, empty_bar * nlevels(data$group) * twolevels, ncol(data)))
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each = empty_bar * twolevels)
 data <- rbind(data, to_add)
 data <- data %>% arrange(group)
 data$id <- rep(seq(1, nrow(data)/twolevels) , each=twolevels)

base_data <- data %>%
  group_by(group) %>%
  summarize(
    start = min(id),
    end = max(id) - empty_bar,
    .groups = 'drop'
  ) %>%
  rowwise() %>%
  mutate(title = mean(c(start, end)))

grid_data <- base_data
grid_data$end <- grid_data$end[c(nrow(grid_data), 1:nrow(grid_data) - 1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1, ]

custom_colors <- c("forward_hits" = "red","reverse_hits" = "blue")

max_tot <- max(data$value, na.rm = TRUE)
upper_limit <- (max_tot + 100)
data$value <- replace(data$value, is.na(data$value), 0)

p <- ggplot(data) +
  geom_bar(aes(x = id, y = value, fill = observation), stat = "identity", alpha = 1) +
  scale_fill_manual(values = custom_colors) +
  geom_segment(data=grid_data, aes(x = end, y = -100, xend = start, yend = -100), colour = "grey", alpha=1, linewidth=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = -50, xend = start, yend = -50), colour = "black", alpha=0.81, linewidth=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, linewidth=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "black", alpha=1, linewidth=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=0.81, linewidth=0.3 , inherit.aes = FALSE ) +
  ggplot2::annotate("text", x = rep(max(data$id),5), y = c(-100, -50, 0, 50, 100), label = c("-100", "-50", "0", "50", "100") , color="black", size=4 , angle=0, fontface="bold", hjust=1) +
  
  ylim(-500, upper_limit) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,5), "cm")
  ) +
  coord_polar() +
  geom_segment(data=base_data, aes(x = start, y = -0, xend = end, yend = -0), colour = "black", alpha=0.4, size=0.2 , inherit.aes = FALSE )  +
  geom_text(data = base_data, aes(x = title, y = +25, label = group), 
            hjust = 0,  # Apply left alignment to all labels
            colour = "red", alpha = 0.5, size = 5, fontface = "bold", inherit.aes = FALSE) 

ggsave(p, filename = "/Users/aarondesouza/Desktop/output.pdf", width = 30, height = 30, units = "cm") 
