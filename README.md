# ozone-and-drought-project-

This is a joint project with Dr. Jesse Berman at UMN and other collaborate. In this project I am responsible for data analysis. 

############ Figure 1 code ##################
#############################################

m(list = ls())
require("rgdal")
require("sf")
#require("sfheaders")
require("ggplot2")

data.2 <- read.csv("Table 2 Ozone by GEOID.csv")
names(data.2)
data.3 <- data.2[ , c("share_ModerateDrought", "lat", "long")]
data.3$type = "moderate"
data.4 <- data.2[ ,c("share_SevereDrought", "lat", "long")]
data.4$type = "severe"
names(data.3) <- names(data.4)
data.all <- rbind(data.3, data.4)
names(data.all)[1] <- c("share")

data.all$id <- 1:NROW(data.all)
data.all.sf <- sf_point( obj = data.all, x = "long", y = "lat")
data.all.sf$id <- 1:NROW(data.all.sf)
data.all.sf <- merge(data.all.sf, data.all, by = "id")

map.sf <- st_read("noaa.4326.shp")
st_crs(data.all.sf) <- st_crs(map.sf)
#plot
plot <- ggplot() + 
  geom_sf(data = map.sf) + 
  # geom_sf(data = dat.sf) + 
  geom_sf(data = data.all.sf, aes(size = share, color = type,
                               geometry = geometry), pch = 1, fill = NA) +
  scale_color_manual(values = c("#000000", "blue")) +  # <---------------- add the colors you prefer
  scale_size(range = c(0, 5), guide = "legend") + 
  labs(size = "% of drought", color = "Type")
plot  

#export
ggsave("plot.circles.png", plot, width = 3150, height = 2000, units = "px")

![0](https://user-images.githubusercontent.com/27037723/218239960-5b800007-a34e-49ae-93b3-ebd13bb4b7cb.jpeg)



######################################################
######################################################
Figure 2 

##############################
# ---- Libreries ----
library(tidyverse)
library(dplyr)
library(ggplot2)
library(plm)
library(fixest)
library(broom)
library(forestplot)

# ---- Data ----
setwd("C:/Users/Chandra/Downloads/")
df <- readRDS("ozone_daily_drought_dat.rds")

# ---- Analysis ----
density_plot <- df %>% 
  ggplot(aes(x=Ozone.Mean,
             fill=USDM.categorical)) +
  geom_density() +
  scale_fill_manual(
    values = c(rgb(1, 255, 27, maxColorValue = 255, alpha = 80),
               rgb(25, 85, 255, maxColorValue = 255, alpha = 80),
               rgb(255, 25, 100, maxColorValue = 255, alpha = 80)),
    labels = c("No Drought", "Moderate Drought", "Severe Drought"),
    name = NULL
  ) +
  theme_classic() + 
  xlab("Ozone mean") + ylab("Density") +
  labs(title="Ozone mean distribution across levels of droughtiness") +
  scale_y_continuous(expand=expansion(mult = c(0,0.1)))

density_plot

df$wday <- lubridate::wday(df$date)

df <- df %>% arrange(date) %>% 
  group_by(date) %>% 
  mutate(
    time = cur_group_id()
  )

# Baseline
model0 <- df %>% 
  lm(data=.,
    formula=Ozone.Mean ~ USDM.categorical
  )

# No fixed effects
model1 <- df[complete.cases(df),] %>% 
  lm(data=.,
     formula=Ozone.Mean ~ USDM.categorical + splines::ns(tmean, 4) + precip + factor(wday) +
       mean_elev + factor(holiday) + I(long^2) + I(lat^2) + long + lat
  )

# Fixed effects
fe_model <- plm(
  Ozone.Mean ~ USDM.categorical + splines::ns(tmean, 4) + precip + factor(wday) +
    mean_elev + factor(holiday) + I(long^2) + I(lat^2) + long + lat,
  data=df,
  index = c("GEOID", "date"),
  model = "within",
#  vcov. = vcovHC, type = "HC1"
)

df$noaa_region <- factor(df$noaa_region)

fe_model_interaction <- plm(
  Ozone.Mean ~ 1 + USDM.categorical + noaa_region + USDM.categorical:noaa_region + splines::ns(tmean,4) + 
    precip + factor(wday) + mean_elev + factor(holiday) + I(long^2) + I(lat^2) + long + lat,
  data=df,
  index = c("GEOID", "date"),
  model = "within",
  vcov. = vcovHC, type = "HC1"
)

#summary(fe_model_interaction)
lmtest::coeftest(fe_model_interaction)

coefs <- fe_model_interaction$coefficients

regions <- as.character(unique(df$noaa_region))
drought_vals <- c("NoDrought", "ModerateDrought")

interaction_values_moderate <- data.frame()
interaction_values_severe <- data.frame()

vcov <- as.data.frame(fe_model_interaction$vcov)

for (x in regions[regions!="northeast"]) {
  name_val_Moderatedrought <- paste0("USDM.categoricalModerateDrought:noaa_region", x)
  name_val_Severedrought <- paste0("USDM.categoricalSevereDrought:noaa_region", x)
  
  pos_ModerateDrought.region <- which(names(coefs)==name_val_Moderatedrought)
  pos_SevereDrought.region <- which(names(coefs)==name_val_Severedrought)
 
  val_ModerateDrought.region <- coefs[pos_ModerateDrought.region]
  val_SevereDrought.region <- coefs[pos_SevereDrought.region]
  
  val_ModerateDrought <- coefs[1]
  val_SevereDrought <- coefs[2]
  
  se_ModerateDrought.region <- vcov[name_val_Moderatedrought,name_val_Moderatedrought]
  cov_ModerateDrought.region <- vcov["USDM.categoricalModerateDrought", name_val_Moderatedrought]
  
  se_ModerateDrought <- vcov["USDM.categoricalModerateDrought","USDM.categoricalModerateDrought"]
  
  se_SevereDrought.region <- vcov[name_val_Severedrought,name_val_Severedrought]
  cov_SevereDrought.region <- vcov["USDM.categoricalSevereDrought", name_val_Severedrought]
  
  se_SevereDrought <- vcov["USDM.categoricalSevereDrought","USDM.categoricalSevereDrought"]
  
  se_interaction_Moderate <- sqrt(se_ModerateDrought + se_ModerateDrought.region +
                                    2*cov_ModerateDrought.region)
  
  se_interaction_Severe <- sqrt(se_SevereDrought + se_SevereDrought.region +
                                    2*cov_SevereDrought.region)
  
  df_values_moderate <- data.frame(
    "region" = x,
    "coef" = val_ModerateDrought.region + val_ModerateDrought,
    "se" = se_interaction_Moderate
  )
  
  df_values_severe <- data.frame(
    "region" = x,
    "coef" = val_SevereDrought.region + val_SevereDrought,
    "se" = se_interaction_Severe
  )
  
  interaction_values_moderate <- bind_rows(interaction_values_moderate,
                                           df_values_moderate)
  
  interaction_values_severe <- bind_rows(interaction_values_severe,
                                           df_values_severe)
}

df_values_moderate <- data.frame(
  "region" = "northeast",
  "coef" = coefs[1],
  "se" = vcov[1, 1]
)

df_values_severe <- data.frame(
  "region" = "northeast",
  "coef" = coefs[2],
  "se" = vcov[2, 2]
)

interaction_values_moderate <- bind_rows(interaction_values_moderate,
                                         df_values_moderate)

interaction_values_severe <- bind_rows(interaction_values_severe,
                                       df_values_severe)

row.names(interaction_values_moderate) <- 1:9
row.names(interaction_values_severe) <- 1:9

interaction_values_moderate$t <- interaction_values_moderate$coef/interaction_values_moderate$se
interaction_values_severe$t <- interaction_values_severe$coef/interaction_values_severe$se

interaction_values_moderate$sig <- interaction_values_moderate$t>1.96
interaction_values_severe$sig <- interaction_values_severe$t>1.96

tab_moderate <- interaction_values_moderate %>% 
  mutate(
    p = pt(t, 1970303, lower.tail = FALSE)
  ) %>% 
  select(-sig, -t)

tab_severe <- interaction_values_severe %>% 
  mutate(
    p = pt(t, 1970303, lower.tail = FALSE)
  ) %>% 
  select(-sig, -t)

rm(list = setdiff(ls(), c("density_plot", "interaction_values_moderate",
                          "interaction_values_severe", "model0", "model1",
                          "fe_model", "fe_model", "tab_moderate", "tab_severe")))

tab2 <- tidy(lmtest::coeftest(fe_model))

tab2$term <- c("Moderate drought", "Severe Drought", "spline 1", "spline 2",
               "spline 3", "spline 4", "Precipitation", "Tuesday", "Wednesday",
              "Thursday", "Friday", "Saturday", "Sunday",
              "Holiday")

interaction_values_moderate$lower <- (interaction_values_moderate$coef - 1.96*interaction_values_moderate$se)
interaction_values_moderate$upper <- (interaction_values_moderate$coef + 1.96*interaction_values_moderate$se)

interaction_values_severe$lower <- (interaction_values_severe$coef - 1.96*interaction_values_severe$se)
interaction_values_severe$upper <- (interaction_values_severe$coef + 1.96*interaction_values_severe$se)

tab_regions <- interaction_values_moderate %>% 
  mutate(Type = "Moderate") %>% 
  bind_rows(
    interaction_values_severe %>% mutate(Type="Severe")
  ) %>% 
  mutate(p = round(pnorm(t, se, lower.tail = FALSE), 3)) %>% 
  select(region, Type, coef, p,lower, upper) %>% 
  mutate(
    val = paste(round(coef, 1), "(", round(lower,1), round(upper, 1), ")")
  ) %>% 
  select(region, Type, val, p) %>% 
  arrange(region, Type)

tab_regions[seq(2, 18, 2), "region"] <- NA

tab_regions <- matrix(c(tab_regions[,1], tab_regions[,2], tab_regions[,3], tab_regions[,4]), ncol=4)
tab_regions <- rbind(c("Region", "Type", "Estimations", "P-Value"), tab_regions)

fig_regions <- interaction_values_moderate %>% 
  bind_rows(
    interaction_values_severe
  ) %>% 
  select(region, coef, lower, upper) %>% 
  rename(
    mean = coef
  ) %>% 
  arrange(region) %>% 
  select(-region)

fig_regions <- rbind(c(NA, NA, NA), fig_regions)

tab_regions[2:19,4] <- tab_regions[2:19,4]

forestplot(tab_regions, fig_regions,
           align="l", 
           is.summary=c(TRUE, rep(FALSE,18)), 
           boxsize=0.01, 
           vertices=T,
           col=fpColors(box="black", line="black"))

save.image("output.RData")



######################################## 
## MY VER
tab_regions <- interaction_values_moderate %>% 
  mutate(Type = "Moderate") %>% 
  bind_rows(
    interaction_values_severe %>% mutate(Type="Severe")
  ) %>% 
  mutate(p = round(pnorm(t, se, lower.tail = FALSE), 3)) %>% 
  select(region, Type, coef, p,lower, upper) %>% 
  mutate(
    val = paste(round(coef, 1), "(", round(lower,1), round(upper, 1), ")")
  ) %>% 
  select(region, Type, val, p) %>% 
  arrange(region, Type)

tab_regions[seq(2, 18, 2), "region"] <- NA

#tab_regions <- matrix(c(tab_regions[,1], tab_regions[,2], tab_regions[,3], tab_regions[,4]), ncol=4)
tab_regions <- rbind(c("Region", "Type", "Estimations", "P-Value"), tab_regions)

tab_regions<-data.frame(tab_regions)
tab_regions$uniq<-LETTERS[seq( from = 19, to = 1 )]
tab_regions[is.na(tab_regions)]<-""
#colnames(tab_regions)<-tab_regions[1,]
#tab_regions<-tab_regions[-1,]
tab_regions$Coll<-rep("white",19)
tab_regions$region<-factor(tab_regions$region)
tab_regions$Type<-factor(tab_regions$Type)




fdf<-ggplot(tab_regions,aes(y=uniq))+
  geom_hline(aes(yintercept=uniq),col=uniq,size=5)+
  geom_text(aes(x=0,label=region),hjust=0)+
  geom_text(aes(x=1,label=Type))+
  geom_text(aes(x=2,label=val),hjust=1)+
  geom_text(aes(x=2.5,label=p),hjust=1)+
  theme_void()

##################################################
fig_regions <- interaction_values_moderate %>% 
  bind_rows(
    interaction_values_severe
  ) %>% 
  select(region, coef, lower, upper) %>% 
  rename(
    mean = coef
  ) %>% 
  arrange(region) %>% 
  select(-region)


fig_regions
fig_regions <- rbind(c(NA,NA,NA), fig_regions)
fig_regions$uniq<-LETTERS[seq( from = 19, to = 1 )]


fp1<-ggplot(fig_regions,aes(x=mean,xmin=lower,xmax=upper,y=uniq))+
  geom_hline(aes(yintercept=factor(-1)),col="grey",size=5)+
  geom_pointrange(shape=22,fill="black",size=.2)+
  geom_vline(xintercept = 0,linetype=3)+
  xlab("Estimations")+
  ylab("")+
  theme_classic()+
  scale_colour_identity()+
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  scale_y_discrete(limits = rev(fig_regions$uniq)) +
  scale_x_continuous(limits = c(-1, 5), 
                     breaks = c(-1, 0, 1 ,2,3, 4,5), 
                     labels = c("-1", "0", "1", "2","3", "4","5"), expand = c(0,0))

fp1

library(gridExtra)
grid.arrange(arrangeGrob(fdf, ncol = 1), 
             fp1,# Second row with 2 plots in 2 different columns
             ncol = 2)
############
##########






