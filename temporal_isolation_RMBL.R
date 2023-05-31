# Analysis of phenological data from the RMBL Phenology Project by Sébastien Rivest (2021),
# for the manuscript "Warmer springs increase potential for temporal reproductive isolation among habitat patches in subalpine flowering plants".
# Three datasets are used: flowering phenology and abundance data from 1974 to 2020 ("pheno.dat"),
# temperature data ("weather.dat"),
# and plots topography and microclimate data ("plots_data").

# Load libraries
Packages <- c("tidyr", "tidyverse", "ggplot2", "ggthemes", "lattice",
              "magrittr", "ggpubr", "extrafont", "gridExtra", "brms",
              "rstan", "bayesplot", "tidybayes", "glue", "ggridges",
              "pracma", "stringr", "forcats", "kulife", "mgcv", "geosphere")
lapply(Packages, library, character.only = TRUE)

# Load phenology dataset
pheno.dat <- read.csv("C:/Users/seb69/OneDrive/Documents/Doc/Assortative_mating/Datasets/pheno_dat.csv")
pheno.dat$species <- as.factor(pheno.dat$species)
pheno.dat <- pheno.dat[!is.na(pheno.dat$doy),]

# Load weather dataset
weather.dat <- read.csv("C:/Users/seb69/OneDrive/Documents/Doc/Assortative_mating/Datasets/Crested_Butte_weather_data.CSV") 
weather.dat$date <- substring(weather.dat$date, 1, 4)
weather.dat <- weather.dat %>% 
  rename(year = date)

# Load plots topography and microclimate dataset
plots_data <- read.csv("C:/Users/seb69/OneDrive/Documents/Doc/Assortative_mating/Datasets/plots_data.csv")

# Measure synchrony among plots.
# Create empty data frame in which the data on flowering synchrony among plots will go.
assort <- data.frame()
for (i in pheno.dat$species %>% droplevels %>% levels){ # Loop for species (approx running time: 15 mins)
  for (year.f in pheno.dat$year[pheno.dat$species == i] 
       %>% as.factor %>% droplevels %>% levels %>% as.numeric){ # Loop for year
    
    for (focal in subset(pheno.dat,species == i & year == year.f) %>%
         droplevels %>% .$plot %>% as.factor %>% levels){ # Loop for plot1
      
      for (donor in subset(pheno.dat,species == i & year == year.f) %>% 
           droplevels %>% .$plot %>% as.factor %>% levels){ # Loop for plot2
        # Measure synchrony if the donor and focal plots are not the same and...
        # if they do not constitute a pair of plots for which synchrony has already been calculated.
        if (donor != focal & 
            any(assort$focal[assort$species == i
                             & assort$year.f == year.f
                             & assort$donor == focal] == donor) == FALSE) {
          # Flower abundance data for the focal plot
          focal.dat <- pheno.dat[pheno.dat$species == i 
                                 & pheno.dat$year == year.f 
                                 & pheno.dat$plot == focal,] %>%
            subset(.,substring(.$date,1,4) ==
                     year.f) %>% unique
          # Flower abundance data for the donor plot
          donor.dat <- pheno.dat[pheno.dat$species == i 
                                 & pheno.dat$year == year.f
                                 & pheno.dat$plot == donor,] %>% 
            subset(.,substring(.$date,1,4) ==
                     year.f) %>% unique
          # Range of days for which the flower curves need to be calculated
          day <- c(min(focal.dat$doy, donor.dat$doy):
                     max(focal.dat$doy, donor.dat$doy))
          # Create flowering curves
          # Estimated values of floral abundance lower than 0 are reported to 0.
          focal.curve <- pchip(focal.dat$doy, focal.dat$floralcount, day)
          0 -> focal.curve[is.na(focal.curve)| focal.curve < 0]
          donor.curve <- pchip(donor.dat$doy, donor.dat$floralcount, day)
          0 -> donor.curve[is.na(donor.curve)| donor.curve < 0]
          
          # Determine the flowering peak of the focal and donor plots.
          flower.peak.focal <- mean(which(focal.curve == max(focal.curve))) + min(day)
          flower.peak.donor <- mean(which(donor.curve == max(donor.curve))) + min(day)
          # Measure overlap.
          overlap <- auc(day, pmin(donor.curve/auc(day,donor.curve),
                                   focal.curve/auc(day,focal.curve)))
          # Measure distance between flowering peaks
          peak.dist <- max(flower.peak.focal,flower.peak.donor) - 
            min(flower.peak.focal,flower.peak.donor)
          # Measure distance between first flowerings dates.
          first.dist <- max(min(day[focal.curve > 0.01]), min(day[donor.curve > 0.01])) -
            min(min(day[focal.curve > 0.01]), min(day[donor.curve > 0.01]))
          # Measure spatial distance between the focal and donor plots.
          spatial.dist <- distm(c(plots_data$longitude[plots_data$plots == focal],
                                  plots_data$latitude[plots_data$plots == focal]),
                                c(plots_data$longitude[plots_data$plots == donor],
                                  plots_data$latitude[plots_data$plots == donor]),
                                fun = distHaversine)   
          
          # Input information for the new data frame
          species <- i
          av.temp <- weather.dat$temp_max[weather.dat$year == year.f & 
                                            weather.dat$month %in% c(4,5,6)] %>%
            na.omit %>% mean
          assort.f <- data.frame(species, focal, donor, year.f, av.temp, 
                                 overlap, peak.dist, first.dist, spatial.dist)
          assort <- rbind(assort, assort.f)
        }
      }
    }
  }
}
# For each species, remove pairs of plots with less than 15 years of observation.
assort$pair <- paste(assort$donor, assort$focal, sep = "-")
rem.dat <- numeric(0)
for (i in 1:nrow(assort)){
  if(assort$year.f[assort$pair == assort$pair[i] &
                   assort$species == assort$species[i]] %>% 
     as.factor %>% levels %>% length < 15) {
    rem.dat <- rbind(rem.dat, i)
  }
} 
assort <- assort[-rem.dat,]

# Defining the habitat type of each plot.
assort$habtype <- paste(substring(assort$focal, 1,2), 
                             substring(assort$donor, 1,2), sep = "-")

# Scaling fixed effect variables.
scaled_data <- assort %>%
  mutate(av.temp_scaled = scale(av.temp),
         year_scaled = scale(year.f),
         spatial.dist_scaled = scale(spatial.dist))

# Model of spring temperature on overlap (approx running time: 10 hours)
model.over <- brm(overlap ~ year_scaled + av.temp_scaled * spatial.dist_scaled + 
                    (-1 + year_scaled + av.temp_scaled|species) +
                    (1 |species/habtype/pair),
                  family = zero_one_inflated_beta(), data = scaled_data, 
                  warmup = 2500, iter = 10000, chains = 4, thin = 10, 
                  control = list(adapt_delta = 0.95),
                  cores = 4, seed = 1234)
print(summary(model.over), digits = 4)
# Verifying if the chains have converged.
plot(model.over)
# Looking for the presence of autocorrelation in the chains.
mcmc_acf_bar(model.over, regex_pars = c("sd"))

# Plots
# Standardized effect sizes of spring temperature on overlap for each species
draws<-spread_draws(model.over, r_species[species,term], b_av.temp_scaled) %>% 
  filter(term == "av.temp_scaled") %>% 
  mutate(b_av.temp_scaled = r_species + b_av.temp_scaled)

pooled.draws <- spread_draws(model.over, b_av.temp_scaled) %>% 
  mutate(species = "Overall effect size")

pheno.data <- bind_rows(draws, pooled.draws) %>% 
  ungroup() %>%
  mutate(species = str_replace_all(species, "[.]", " ")) %>% 
  mutate(species = reorder(species, b_av.temp_scaled))

pheno.data.summary <- group_by(pheno.data, species) %>% 
  mean_qi(b_av.temp_scaled)

ggplot(aes(b_av.temp_scaled, relevel(species, "Overall effect size", after = Inf),
           fill = factor(stat(quantile))), data = pheno.data) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantiles = c(0.025, 0.975), scale = 5
  ) +
  scale_fill_manual(
    name = "Probability", values = c("#ffffff", "#7D96C496", "#ffffff"),
    labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")
  ) +
  geom_vline(xintercept = 0, color = "black", linetype="dashed", si5ze = 1, alpha = 1) +
  coord_cartesian(xlim = c(-0.45, 0.15)) +
  labs(x = "Standardized change in overlap",
       y = element_blank()) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(face = ifelse((
          c(1,rep(0,scaled_data$species %>% as.factor %>% levels %>% length)) == 0),
          "italic","bold")))

# Overall relationship between spring temperature and overlap
ce.overt<-conditional_effects(model.over, method = "fitted", effects = "av.temp_scaled")
# Scale back spring temperature.
ce.overt$av.temp_scaled$effect1__ <- ce.overt$av.temp_scaled$effect1__ *
  sd(assort$av.temp) + mean(assort$av.temp)
# plot
plot(ce.overt, plot = FALSE, col = "black")[[1]] + 
  geom_line(colour = "#000000", size = 1.75) +
  coord_cartesian(ylim = c(0, 0.7)) +
  scale_x_continuous(breaks = round(seq(min(11.5), max(19), by = 1.5),1)) +
  labs(x = "Spring temperature (°C)", y = "Overlap") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# Plot of the relationship for all species
mycolors <- colorRampPalette(c("#2e5ba3","#91a9cf"))(assort$species %>% as.factor %>% levels %>% length)
mycolors2 <- colorRampPalette("#ffffff")(assort$species %>% as.factor %>% levels %>% length)
ce.over <- conditional_effects(model.over,
                               method = "fitted",
                               effects = "av.temp_scaled:species",
                               re_formula = NULL)
ce.over$`av.temp_scaled:species`$effect1__ <- ce.over$`av.temp_scaled:species`$effect1__ *
  sd(assort$av.temp) + mean(assort$av.temp)
plot(ce.over, plot = FALSE, col = "black")[[1]] +
  scale_fill_manual(values = mycolors2) +
  scale_color_manual(values = mycolors) +
  geom_line(size = 1, alpha = 1) +
  coord_cartesian(ylim = c(0, 0.75)) +
  scale_x_continuous(breaks = round(seq(min(11.5), max(19), by = 1.5),1)) +
  labs(x = "Spring temperature (°C)", y = "Overlap") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.position = "none")  

# Model of spring temperature on distance between flowering peaks (approx running time: 5 hours)
model.pd <- brm(peak.dist ~ year_scaled + av.temp_scaled * spatial.dist_scaled + 
                  (-1 + year_scaled + av.temp_scaled|species) +
                  (1 |species/habtype/pair),
                family = hurdle_gamma(), data = scaled_data, 
                warmup = 2500, iter = 10000, chains = 4, thin = 10, 
                control = list(adapt_delta = 0.95, max_treedepth = 12),
                cores = 4, seed = 1234)
print(summary(model.pd), digits = 4)
plot(model.pd)
mcmc_acf_bar(model.pd, regex_pars = c("sd"))

# Model of spring temperature on first flowering dates (approx running time: 5 hours)
model.first <- brm(first.dist ~ year_scaled + av.temp_scaled * spatial.dist_scaled + 
                     (-1 + year_scaled + av.temp_scaled|species) +
                     (1 |species/habtype/pair),
                   family = hurdle_gamma(), data = scaled_data, 
                   warmup = 2500, iter = 10000, chains = 4, thin = 10, 
                   control = list(adapt_delta = 0.95, max_treedepth = 12),
                   cores = 4, seed = 1234)
print(summary(model.first), digits = 4)
plot(model.first)
mcmc_acf_bar(model.first, regex_pars = c("sd"))

# Measure flowering duration for each species-plot-year combination.
duration  <- data.frame()
for (i in assort$species %>% as.factor %>% droplevels %>% levels){ # approx running time: 2 minutes
  for (year.p in assort$year.f[assort$species == i] %>%
       as.factor %>% droplevels %>% levels %>% as.numeric){
    for (micro in c(subset(assort,species == i & year.f == year.p) %>% 
                    droplevels  %>% .$focal %>% as.factor %>% levels,
                    subset(assort,species == i & year.f == year.p) %>% 
                    droplevels  %>% .$donor %>% as.factor %>% levels) %>% unique){
      
      focal.dat <- pheno.dat[pheno.dat$species == i 
                             & pheno.dat$year == year.p 
                             & pheno.dat$plot == micro,]
      focal.dat <- subset(focal.dat,substring(focal.dat$date,1,4) == year.p)
      focal.dat <- unique(focal.dat)
      
      day <- c(min(focal.dat$doy):
                 max(focal.dat$doy))
      
      # Calculate the flowering curve of the focal plot.
      val <- data.frame(doy = seq(min(day), max(day), by = 1))
      
      focal.curve <- pchip(focal.dat$doy, focal.dat$floralcount, day)
      0 -> focal.curve[is.na(focal.curve)| focal.curve <= 0]
      
      flowering.length <- max(day[focal.curve > 0]) - 
        min(day[focal.curve > 0])
      
      flower.peak.micro <- mean(which(focal.curve == max(focal.curve))) + min(day)
      
      species <- i 
      av.temp <- weather.dat$temp_max[weather.dat$year == year.p & 
                                        weather.dat$month %in% c(4,5,6)] %>%
        na.omit %>% mean
      
      duration.f <- data.frame(species, year.p, micro, av.temp, flowering.length, flower.peak.micro)
      duration <- rbind(duration, duration.f)
    }
  }
}
scaled_duration <- duration %>%
  mutate(av.temp_scaled = scale(av.temp),
         year_scaled = scale(year.p))
scaled_duration$habtype <- paste(substring(scaled_duration$micro, 1,2))

# Model of spring temperature on flowering duration (approx running time: 1 hour)
model.d <- brm(flowering.length ~ year_scaled + av.temp_scaled +
                 (-1 + year_scaled + av.temp_scaled |species) +
                 (1 |species/habtype/micro),
               family = Gamma(link = log), data = scaled_duration, 
               warmup = 2500, iter = 10000, chains = 4, thin = 10, 
               control = list(adapt_delta = 0.95),
               cores = 4, seed = 1234)
print(summary(model.d), digits = 3)
plot(model.d)
mcmc_acf_bar(model.d, regex_pars = c("sd"))

# Plot of the effect size of spring temperature on days between first flowerings,
# days between flowering peaks and flowering duration
effect.sizes <- data.frame(response=c("Days between first flowerings",
                                      "Days between flowering peaks",
                                      "Flowering duration"), index=1:3,
                           effect=c(fixef(model.first, summary = TRUE)[3,1],
                                    fixef(model.pd, summary = TRUE)[3,1],
                                    fixef(model.d, summary = TRUE)[3,1]),
                           lower=c(fixef(model.first, summary = TRUE)[3,3],
                                   fixef(model.pd, summary = TRUE)[3,3],
                                   fixef(model.d, summary = TRUE)[3,3]),
                           upper=c(fixef(model.first, summary = TRUE)[3,4]
                                   , fixef(model.pd, summary = TRUE)[3,4],
                                   fixef(model.d, summary = TRUE)[3,4]))
ggplot(data=effect.sizes, aes(y=index, x=effect, xmin=lower, xmax=upper)) +
  geom_pointrange(position = position_dodge(width = 1.1)) +
  geom_errorbarh(height = 0.1, size = 0.8) +
  geom_point(size = 4, col =c("#94a6c5","#4982d9","#2d4f8c")) +
  geom_point(shape = 1, size = 4, col = "black") +
  scale_y_continuous(name = "", breaks=1:nrow(effect.sizes), labels=effect.sizes$response) +
  labs(x='Effect Size', y = 'Study') +
  geom_vline(xintercept = 0, color='black', linetype='dashed', alpha=1) +
  coord_cartesian(xlim = c(-0.08, 0.08)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 13))

# Plot of the effect size of spring temperature on flowering overlap
effect.sizes <- data.frame(response=c("Flowering overlap"), index=1,
                           effect=c(fixef(model.over, summary = TRUE)[3,1]),
                           lower=c(fixef(model.over, summary = TRUE)[3,3]),
                           upper=c(fixef(model.over, summary = TRUE)[3,4]))
ggplot(data=effect.sizes, aes(y=index, x=effect, xmin=lower, xmax=upper)) +
  geom_pointrange(position = position_dodge(width = 1.1)) +
  geom_errorbarh(height = 0.05, size = 0.8) +
  geom_point(size = 4, col = c("#d08889")) +
  geom_point(shape = 1, size = 4, col = "black") +
  scale_y_continuous(name = "", breaks=1:nrow(effect.sizes), labels=effect.sizes$response) +
  labs(x='Effect Size', y = 'Study') +
  geom_vline(xintercept = 0, color='black', linetype='dashed', alpha=1) +
  coord_cartesian(xlim = c(-0.15, 0.15)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 13))

# Test the effect of plot and species charecteristics on overlap
# Measure differences in microclimatic variables between plots
assort$snowmelt.diff <- "NA"
assort$slope.diff <- "NA"
assort$aspect.diff <- "NA"
for (focal in assort$focal %>% as.factor %>% levels){ # Loop for plot1
  for (donor in assort$donor %>% as.factor %>% levels){ # Loop for plot2
    if (donor != focal) {

assort$snowmelt.diff[assort$focal == focal & assort$donor == donor] <- abs(plots_data$snowmelt_date[plots_data$plot == focal] - 
                    plots_data$snowmelt_date[plots_data$plot == donor])
assort$slope.diff[assort$focal == focal & assort$donor == donor] <- abs(plots_data$slope[plots_data$plot == focal] - 
                    plots_data$slope[plots_data$plot == donor])
assort$aspect.diff[assort$focal == focal & assort$donor == donor] <- abs(plots_data$aspect[plots_data$plot == focal] - 
                     plots_data$aspect[plots_data$plot == donor])
    }
  }
}

# Measure species mean phenology
assort$phenology <- "NA"
assort$pollination_mode <- "NA"
for (i in duration$species %>% as.factor %>% levels){
  assort$phenology[assort$species == i] <- flowering_info$flower.peak.micro[flowering_info$species == i]
  assort$pollination_mode[assort$species == i] <- flowering_info$pollination_mode[flowering_info$species == i]
} 

scaled_data <- assort %>%
  mutate(av.temp_scaled = scale(av.temp),
         year_scaled = scale(year.f),
         spatial.dist_scaled = scale(spatial.dist),
         phenology_scaled = scale(as.numeric(phenology)),
         snowmelt.diff_scaled = scale(as.numeric(snowmelt.diff)),
         slope.diff_scaled = scale(as.numeric(slope.diff)),
         aspect.diff_scaled = scale(as.numeric(aspect.diff)))

# Model of overlap in function of temperature, plot and species variables, and their interaction
model.over.micro <- brm(overlap ~ year_scaled + av.temp_scaled * spatial.dist_scaled + 
                          av.temp_scaled * snowmelt.diff_scaled + av.temp_scaled * slope.diff_scaled +
                          av.temp_scaled * aspect.diff_scaled + av.temp_scaled * phenology_scaled +
                          (-1 + year_scaled + av.temp_scaled|species) +
                          (1 |species/habtype/pair),
                          family = zero_one_inflated_beta(), data = scaled_data, 
                          warmup = 2500, iter = 10000, chains = 4, thin = 10, 
                          control = list(adapt_delta = 0.95),
                          cores = 4, seed = 1234)
print(summary(model.over.micro), digits = 4)
plot(model.over.micro)
mcmc_acf_bar(model.over.micro, regex_pars = c("sd"))

# plots of the effect of spatial distance and microclimatic variables on overlap and change in overlap with temperature
# Spatial distance
# Determine range of values of spatial distance
condition <- data.frame(spatial.dist_scaled = c((50 - mean(assort$spatial.dist))/sd(assort$spatial.dist),
                                                (200 - mean(assort$spatial.dist))/sd(assort$spatial.dist),
                                                (800 - mean(assort$spatial.dist))/sd(assort$spatial.dist)))
ce.spatial<-conditional_effects(model.over6,
                              method = "fitted",
                              effects = "av.temp_scaled:spatial.dist_scaled",
                              int_conditions = condition)
ce.spatial$`av.temp_scaled:spatial.dist_scaled`$av.temp_scaled <- ce.spatial$`av.temp_scaled:spatial.dist_scaled`$av.temp_scaled *
  sd(assort$av.temp) + mean(assort$av.temp)
df <- as.data.frame(ce.spatial$`av.temp_scaled:spatial.dist_scaled`)
# Plot
ggplot(df,aes(x=av.temp_scaled,y=estimate__, group=as.factor(spatial.dist_scaled)))+
  geom_ribbon(aes(ymin=lower__, ymax=upper__, fill = as.factor(spatial.dist_scaled)), alpha=0.15)+
  geom_line(size=1, position=position_dodge(0.05), aes(color=as.factor(spatial.dist_scaled), linetype=as.factor(spatial.dist_scaled))) +
  scale_color_manual(values = c("#030303","#616161","#858585")) +
  scale_fill_manual(values = c("#030303","#616161","#858585")) +
  coord_cartesian(ylim = c(0.30, 0.62)) +
  labs(x = "Spring temperature (°C)", y = "Overlap") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# Slope
condition <- data.frame(slope.diff_scaled = c(mean(scaled_data$slope.diff_scaled)-sd(scaled_data$slope.diff_scaled),
                                              mean(scaled_data$slope.diff_scaled),
                                              mean(scaled_data$slope.diff_scaled)+sd(scaled_data$slope.diff_scaled)))
ce.slope<-conditional_effects(model.over6,
                              method = "fitted",
                              effects = "av.temp_scaled:slope.diff_scaled",
                              int_conditions = condition)
ce.slope$`av.temp_scaled:slope.diff_scaled`$av.temp_scaled <- ce.slope$`av.temp_scaled:slope.diff_scaled`$av.temp_scaled *
  sd(assort$av.temp) + mean(assort$av.temp)
df <- as.data.frame(ce.slope$`av.temp_scaled:slope.diff_scaled`)

ggplot(df,aes(x=av.temp_scaled,y=estimate__, group=as.factor(slope.diff_scaled)))+
  geom_ribbon(aes(ymin=lower__, ymax=upper__, fill = as.factor(slope.diff_scaled)), alpha=0.15)+
  geom_line(size=1, position=position_dodge(0.05), aes(color=as.factor(slope.diff_scaled), linetype=as.factor(slope.diff_scaled))) +
  scale_color_manual(values = c("#284b41","#4b8879","#5da994")) +
  scale_fill_manual(values = c("#284b41","#4b8879","#5da994")) +
  coord_cartesian(ylim = c(0.30, 0.62)) +
  labs(x = "Spring temperature (°C)", y = "Overlap") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# Snowmelt date
condition <- data.frame(snowmelt.diff_scaled = c(mean(scaled_data$snowmelt.diff_scaled)-sd(scaled_data$snowmelt.diff_scaled),
                                              mean(scaled_data$snowmelt.diff_scaled),
                                              mean(scaled_data$snowmelt.diff_scaled)+sd(scaled_data$snowmelt.diff_scaled)))
ce.snowmelt<-conditional_effects(model.over6,
                              method = "fitted",
                              effects = "av.temp_scaled:snowmelt.diff_scaled",
                              int_conditions = condition)
ce.snowmelt$`av.temp_scaled:snowmelt.diff_scaled`$av.temp_scaled <- ce.snowmelt$`av.temp_scaled:snowmelt.diff_scaled`$av.temp_scaled *
  sd(assort$av.temp) + mean(assort$av.temp)
df <- as.data.frame(ce.snowmelt$`av.temp_scaled:snowmelt.diff_scaled`)

ggplot(df,aes(x=av.temp_scaled,y=estimate__, group=as.factor(snowmelt.diff_scaled)))+
  geom_ribbon(aes(ymin=lower__, ymax=upper__, fill = as.factor(snowmelt.diff_scaled)), alpha=0.15)+
  geom_line(size=1, position=position_dodge(0.05), aes(color=as.factor(snowmelt.diff_scaled), linetype=as.factor(snowmelt.diff_scaled))) +
  scale_color_manual(values = c("#3f5664","#5f839b","#82b3d6")) +
  scale_fill_manual(values = c("#3f5664","#5f839b","#82b3d6")) +
  coord_cartesian(ylim = c(0.30, 0.62)) +
  labs(x = "Spring temperature (°C)", y = "Overlap") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# Aspect
condition <- data.frame(aspect.diff_scaled = c(mean(scaled_data$aspect.diff_scaled)-sd(scaled_data$aspect.diff_scaled),
                                               mean(scaled_data$aspect.diff_scaled),
                                               mean(scaled_data$aspect.diff_scaled)+sd(scaled_data$aspect.diff_scaled)))
ce.aspect<-conditional_effects(model.over6,
                              method = "fitted",
                              effects = "av.temp_scaled:aspect.diff_scaled",
                              int_conditions = condition)
ce.aspect$`av.temp_scaled:aspect.diff_scaled`$av.temp_scaled <- ce.aspect$`av.temp_scaled:aspect.diff_scaled`$av.temp_scaled *
  sd(assort$av.temp) + mean(assort$av.temp)
df <- as.data.frame(ce.aspect$`av.temp_scaled:aspect.diff_scaled`)

ggplot(df,aes(x=av.temp_scaled,y=estimate__, group=as.factor(aspect.diff_scaled)))+
  geom_ribbon(aes(ymin=lower__, ymax=upper__, fill = as.factor(aspect.diff_scaled)), alpha=0.15)+
  geom_line(size=1, position=position_dodge(0.05), aes(color=as.factor(aspect.diff_scaled), linetype=as.factor(aspect.diff_scaled))) +
  scale_color_manual(values = c("#362f14","#7a6a2b","#bca245")) +
  scale_fill_manual(values = c("#362f14","#7a6a2b","#bca245")) +
  coord_cartesian(ylim = c(0.30, 0.6)) +
  labs(x = "Spring temperature (°C)", y = "Overlap") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

