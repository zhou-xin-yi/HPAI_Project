# packages
library(dplyr)
library(tidyr)
library(phylotools)

library(ggplot2)
library(usmap)
library(viridis)
library(gridExtra)

### bird flu detections in wild birds
birds <- read.csv("./data/hpai-wild-birds.csv")
summary(birds)
birds[birds$State == "michigan", c("State")] = "Michigan"
colnames(birds)[1] = "state"

birds$Collection.Date <- as.Date(birds$Collection.Date, "%m/%d/%Y")

birds_mean <- birds %>%
  filter(Collection.Date > '2023-07-31') %>% # from August 2023
  group_by(state) %>%
  summarize(count = n()) 


pA <- plot_usmap(data = birds_mean, values = "count", labels = TRUE) + 
  scale_fill_continuous(
    low = "white", high = "red", name = "wild birds", label = scales::comma
  ) + 
  theme(legend.position = "right") +
  labs(tag = "A")

### cattle outbreak
cattle <- read.csv("./data/cattle.csv")
colnames(cattle)[2] = "state"
cattle$Confirmed <- as.Date(cattle$Confirmed, "%d-%b-%y")
colnames(cattle)[1] = "Collection.Date"

cattle_mean <- cattle %>%
  group_by(state) %>%
  summarize(count = n()) 


pB <- plot_usmap(data = cattle_mean, values = "count", labels = TRUE) + 
  scale_fill_continuous(
    low = "white", high = "red", name = "cattle", label = scales::comma
  ) + 
  theme(legend.position = "right") +
  labs(tag = "B")

### correlation
birds_cattle <- inner_join(birds_mean, cattle_mean, by = "state")
pC <- ggplot(birds_cattle, aes(x=count.x, y=count.y)) + 
  geom_point(mapping = aes(colour = state)) +
  geom_smooth(method = lm, formula=y~0+x) +
  coord_equal() +
  xlab("Reported no. wild bird cases") +
  ylab("Confirmed no. cattle cases") +
  theme_bw() + 
  labs(tag = "C")

p2 <- grid.arrange(pA, pB, pC, nrow = 2)

ggsave("../figures/figure2.png", p2, width = 16, height = 9, dpi = 400)

cor.test(birds_cattle$count.x, birds_cattle$count.y)

# Pearson's product-moment correlation
# 
# data:  birds_cattle$count.x and birds_cattle$count.y
# t = 4.9288, df = 12, p-value = 0.0003486
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.5081417 0.9404725
# sample estimates:
#       cor 
# 0.8181422  

### time series
birds_time <- birds %>%
  group_by(Collection.Date) %>%
  summarise(n_birds = n())

cattle_time <- cattle %>%
  group_by(Collection.Date) %>%
  summarise(n_cattle = n())

n_time <- full_join(birds_time, cattle_time, by = c("Collection.Date"))

p1 <- ggplot(n_time, aes(x = Collection.Date)) +
  geom_line(aes(y = n_birds, color = "Wild birds")) + 
  geom_line(aes(y = n_cattle, color = "Cattle")) + 
  geom_vline(xintercept = as.numeric(as.Date("2023-08-01")), 
             linetype = "dashed", 
             color = "blue") + # start of 2023 high prevalence in wild birds
  geom_vline(xintercept = as.numeric(as.Date("2024-03-25")), 
             linetype = "dotdash", 
             color = "blue") + # first confirmed cattle case
  ylab("Daily no. reported cases") +
  xlab("") +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 month") +
  theme_bw()

ggsave("../figures/figure1.png", p1, width = 16, height = 9, dpi = 400)


