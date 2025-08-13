
library(PaleoSpec)
library(RColorBrewer)
library(tidyverse)
library(cmdstanr)
library(scales)

# Data --------------------------------------------------------------------

#General function for finding spectra from potentially unevenly spaced data
find_spec <- function(x, d18O, smooth = NULL){
  
  reg_x <- seq(x[1], x[length(x)], length.out = length(x))
  
  reg_d18O <- approx(x, d18O, reg_x)$y
  
  dx <- reg_x[2] - reg_x[1]
  ts <- ts(reg_d18O, start = reg_x[1], deltat = dx)
  sp <- spectrum(ts, plot = FALSE, spans = smooth, taper = 0.1, pad = 0, fast = T, detrend = T)
  
  freq <- sp$freq
  spec <- sp$spec
  
  return(list(freq = freq, spec = spec))
}

#Load Dome C data as a tibble
DomeC <- read_tsv("./Vasileios_etal-2021.tab")

#Transform the data
DomeC <- DomeC[-c(1:37),] %>%
  separate(names(DomeC), into = c("depth_top", "depth_bot", "depth", "age", "d18O"), sep = "\t") %>%
  transmute(depth_top = as.numeric(depth_top),
            depth_bot = as.numeric(depth_bot),
            depth = as.numeric(depth),
            age = as.numeric(age),
            d18O = as.numeric(d18O)
  )

#Linearly interpolate gaps in data
DomeC$d18O <- approx(DomeC$depth, DomeC$d18O, DomeC$depth)$y

#Define time-series for analysis
MIS1 <- DomeC[1:3150, ]

dt <- diff(DomeC$age)

bin_size <- c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)

t_max <- vector()
binned_sec <- list()
spec <- list()
spectra <- list()

for(i in 1:length(bin_size)){
  
  #Finds the age cut-off 't_max' where the temporal resolution is greater than the current bin size (and therefore cannot be binned without potentially creating gaps)
  if(bin_size[i] < max(dt)){
    t_max[i] <- min(which(dt >= bin_size[i])) - 1
  }else{
    t_max[i] <- length(DomeC$d18O)
  }
  
  #Bins the data to current bin size. Only includes data which has temporal resolution > bin size
  binned_sec[[i]] <- AvgToBin(x = DomeC$age[1:t_max[i]], y = DomeC$d18O[1:t_max[i]], breaks = seq(DomeC$age[1], DomeC$age[t_max[i]], by = bin_size[i]))
  
  #Computes the power spectrum of the binned section
  spec[[i]] <- find_spec(binned_sec[[i]]$centers, binned_sec[[i]]$avg)
  
  #Converts into format compatible with MeanSpectrum
  spectra[[i]] <- list(freq = spec[[i]]$freq, spec = spec[[i]]$spec, dof = rep(1, length(spec[[i]]$freq)))
  
}

#Computes the composite mean spectrum of all the binned spectra
mean_P0 <- MeanSpectrum(specList = spectra)


pdf("./fig_A1.pdf",
    width = 10,
    height = 3)

par(mar = c(5,5,1,7))

plot(DomeC$age, DomeC$d18O, 'l', xlab = 'Age (kya)', ylab = expression(paste(delta^{18}, "O")))

colfun <- brewer.pal(n = length(bin_size), 'Accent')

for(i in 1:length(bin_size)){
  
  abline(v = max(binned_sec[[i]]$centers), col = colfun[i], lwd = 2)
  
}

legend(850, -43.5, legend = (as.character(bin_size)), title = 'bin size (kyr)', col = colfun, lwd = 2, bty = 'n', xpd = TRUE)

dev.off()

pdf("./fig_A2.pdf",
    width = 6,
    height = 5)

plot(spec[[1]]$freq, spec[[1]]$spec, 'l', log = 'xy', col = colfun[1], lwd = 1, xlim = c(1e-3, 100), ylim = c(1e-6, 1e3),
     xlab = 'Frequency (1/ky)', ylab = 'PSD')

for(i in 2:length(bin_size)){
  
  lines(spec[[i]]$freq, spec[[i]]$spec, col = colfun[i], lwd = 1)
  
}

lines(mean_P0$freq, mean_P0$spec, 'l', lwd = 2)

legend(7e-4, 1, legend = as.character(bin_size), title = 'bin size (ky)', col = colfun, lwd = 2, bty = 'n')
legend(7e-4, 5e-6, 'mean', col = 'black', bty = 'n', lwd = 2)

dev.off()


#Model universal P0 fit
P0_mod <- cmdstan_model("./Universal_P0_fit.stan")

stan_data <- list(N = length(mean_P0$freq[which(mean_P0$freq >= 0.01 & mean_P0$freq <= 10)]), freq = mean_P0$freq[which(mean_P0$freq >= 0.01 & mean_P0$freq <= 10)], spec = mean_P0$spec[which(mean_P0$freq >= 0.01 & mean_P0$freq <= 10)],
                  alpha_mu = 1, alpha_sigma = 10, beta_mu = 1, beta_sigma = 1, epsilon_mu = 0.01, epsilon_sigma = 1, phi_scale = 1)

mean_P0_params <- P0_mod$sample(data = stan_data, chains = 4, parallel_chains = 4, refresh = 500)

alpha_est <- mean(mean_P0_params$draws('alpha'))
alpha_sd <- sd(mean_P0_params$draws('alpha'))
beta_est <- mean(mean_P0_params$draws('beta'))
beta_sd <- sd(mean_P0_params$draws('beta'))
epsilon_est <- mean(mean_P0_params$draws('epsilon'))
epsilon_sd <- sd(mean_P0_params$draws('epsilon'))

mean_P0_fit <- alpha_est*mean_P0$freq[which(mean_P0$freq >= 0.01 & mean_P0$freq <= 10)]^-beta_est + epsilon_est
mean_P0_fit_min <- mean_P0_params$summary()$q5[-(1:5)]
mean_P0_fit_max <- mean_P0_params$summary()$q95[-(1:5)]

mean_P0_df <- tibble(freq = mean_P0$freq,
                     spec = mean_P0$spec)

mean_P0_fit_df <- tibble(freq = mean_P0$freq[which(mean_P0$freq >= 0.01 & mean_P0$freq <= 10)],
                         fit = mean_P0_fit,
                         fit_min = mean_P0_fit_min,
                         fit_max = mean_P0_fit_max)

pdf(file = "./fig_2.pdf",
    width = 8,
    height = 5)

ggplot(data = mean_P0_fit_df) +
  geom_rect(aes(xmin = freq[1], xmax = freq[length(freq)], ymin = 0, ymax = Inf, fill = 'Fitted Region'), alpha = 1) +
  geom_line(data = mean_P0_df, aes(freq, spec, colour = 'spectrum')) +
  geom_line(aes(freq, fit, colour = 'fit')) +
  geom_ribbon(aes(x = freq, ymin = fit_min, ymax = fit_max), colour = NA, fill = 'red', alpha = 0.5) +
  scale_x_continuous(trans = 'log', breaks = 10^(-2:3)) +
  scale_y_continuous(trans = 'log', breaks = 10^(-6:2)) +
  scale_colour_manual(name = expression(mean~P[0]), values = c('black', 'red'), breaks = c('spectrum', 'fit')) +
  scale_fill_manual(name = '', values = 'grey92') +
  xlab("Frequency (1/kyr)") +
  ylab("PSD") +
  theme_bw()

dev.off()


# Load LDC age-depth model and diffusion length data ----------------------

LDC <- read_tsv("./LDC-VHF_updated.txt") #Updated version
LDC <- LDC[-(1:10), ] %>% 
  separate(names(LDC), into = c("depth",	"age",	"sigma_age",	"age_density",	"accumulation", "thinning",	"vertical_velocity"), sep = "\t")
LDC <- as.tibble(sapply(LDC, as.numeric))
LDC <- LDC[which(LDC$age < 1.5e6), ]

vas <- read_tsv("./vas_sigmas.txt")


# Finding required noise for a given frequency ---------------------------------------------------------------

#LDC layer thickness
LDC_bdot <- diff(LDC$depth)/diff(LDC$age) #In metres per year
LDC_bdot <- LDC_bdot*1000 #In metres per ky

#Shortening age-depth model timeseries by 1 to match length of LDC_bdot
LDC <- list(depth = LDC$depth[-1] - diff(LDC$depth),
            age = LDC$age[-1] - diff(LDC$age))

dz <- 0.11 #In metres
dt <- dz/LDC_bdot #In kiloyears

#Interpolate diffusion lengths to age-depth model
s <- approx(vas$Depth, vas$sigma_total, LDC$depth)$y

#Choose range of frequencies for figure (desired temporal frequencies)
ft_des <- 10^(seq(-3, 4, 0.001)) #per ky

#Define P0 using estimated parameters from Dome C
P0 <- alpha_est*ft_des^-beta_est + epsilon_est

A <- 0.25 #Fraction of initial power after Wiener deconvolution

#f_max_H is f_max_diff from Eq. 13. Diffusion length is converted into time domain
f_max_H <- sqrt(- (log(A)) / (2*pi^2 * (s/LDC_bdot)^2 ) )

noise_req <- list()
frac_ind <- vector()

for(i in 1:length(LDC_bdot)){
  
  H <- exp(-0.5 * (2*pi * (ft_des) * s[i]/LDC_bdot[i])^2)
  
  N_req <- P0 * H^2 * (1/A - 1) #Eq. 17
  
  noise_req[[i]] <- sqrt(N_req/dt[i]) #Convert from variance to SD
  
  frac_ind[i] <- max(which(H > A)) #Finds indices for f_max_diff (vertical dashed lines in Fig. 5)
  
  print(i)
}

select_ages <- c(5, 10, 20, 50, 100, 200, 500, 1000, 1500) #Selected ages in ky

select_inds <- vector()

for(i in 1:length(select_ages)){
  
  select_inds[i] <- which.min(abs(LDC$age/1000 - select_ages[i])) #Find indices for selected ages
  
}

inds <- select_inds

# BE-OIC f_max figures (Figs. 5a + 5b) -------------------------------------------------

#Co-ordinates of dots in Fig. 5a
H_x <- ft_des[frac_ind]
H_y1 <- noise_req[[inds[1]]][frac_ind[inds[1]]]
H_y2 <- noise_req[[inds[2]]][frac_ind[inds[2]]]
H_y3 <- noise_req[[inds[3]]][frac_ind[inds[3]]]
H_y4 <- noise_req[[inds[4]]][frac_ind[inds[4]]]
H_y5 <- noise_req[[inds[5]]][frac_ind[inds[5]]]
H_y6 <- noise_req[[inds[6]]][frac_ind[inds[6]]]
H_y7 <- noise_req[[inds[7]]][frac_ind[inds[7]]]
H_y8 <- noise_req[[inds[8]]][frac_ind[inds[8]]]
H_y9 <- noise_req[[inds[9]]][frac_ind[inds[9]]]

#Reformatting data for ggplot
fgain_df_Age_5 <- tibble(ft_des = ft_des, ft_des_rel = ft_des/f_max_H[inds[1]], noise_req = noise_req[[inds[1]]], H_x = H_x[inds[1]], H_y = H_y1, Age = '5')
fgain_df_Age_10 <- tibble(ft_des = ft_des, ft_des_rel = ft_des/f_max_H[inds[2]], noise_req = noise_req[[inds[2]]], H_x = H_x[inds[2]], H_y = H_y2, Age = '10')
fgain_df_Age_20 <- tibble(ft_des = ft_des, ft_des_rel = ft_des/f_max_H[inds[3]], noise_req = noise_req[[inds[3]]], H_x = H_x[inds[3]], H_y = H_y3, Age = '20')
fgain_df_Age_50 <- tibble(ft_des = ft_des, ft_des_rel = ft_des/f_max_H[inds[4]], noise_req = noise_req[[inds[4]]], H_x = H_x[inds[4]], H_y = H_y4, Age = '50')
fgain_df_Age_100 <- tibble(ft_des = ft_des, ft_des_rel = ft_des/f_max_H[inds[5]], noise_req = noise_req[[inds[5]]], H_x = H_x[inds[5]], H_y = H_y5, Age = '100')
fgain_df_Age_200 <- tibble(ft_des = ft_des, ft_des_rel = ft_des/f_max_H[inds[6]], noise_req = noise_req[[inds[6]]], H_x = H_x[inds[6]], H_y = H_y6, Age = '200')
fgain_df_Age_500 <- tibble(ft_des = ft_des, ft_des_rel = ft_des/f_max_H[inds[7]], noise_req = noise_req[[inds[7]]], H_x = H_x[inds[7]], H_y = H_y7, Age = '500')
fgain_df_Age_1000 <- tibble(ft_des = ft_des, ft_des_rel = ft_des/f_max_H[inds[8]], noise_req = noise_req[[inds[8]]], H_x = H_x[inds[8]], H_y = H_y8, Age = '1000')
fgain_df_Age_1500 <- tibble(ft_des = ft_des, ft_des_rel = ft_des/f_max_H[inds[9]], noise_req = noise_req[[inds[9]]], H_x = H_x[inds[9]], H_y = H_y9, Age = '1500')

fgain_df <- bind_rows(fgain_df_Age_5,
                      fgain_df_Age_10,
                      fgain_df_Age_20,
                      fgain_df_Age_50,
                      fgain_df_Age_100,
                      fgain_df_Age_200,
                      fgain_df_Age_500,
                      fgain_df_Age_1000,
                      fgain_df_Age_1500)
fgain_df$Age <- factor(fgain_df$Age, levels = c('5', '10', '20', '50', '100', '200', '500', '1000', '1500'))

colfun1 <- viridisLite::viridis(n = length(inds))

fig_5a <- ggplot(data = fgain_df) +
  geom_line(aes(x = ft_des, y = noise_req, colour = Age)) +
  geom_hline(aes(yintercept = 0.1, size = '0.1'), colour = '#fb8500', linetype = 'dashed') +
  geom_hline(aes(yintercept = 0.01, size = '0.01'), colour = '#219ebc', linetype = 'dashed') +
  geom_point(aes(x = H_x, y = H_y, colour = Age), shape = 19) +
  geom_vline(aes(xintercept = H_x, colour = Age), linetype = 'dotted') +
  coord_cartesian(xlim=c(1e-2, 1e3), ylim=c(1e-3, 10)) +
  scale_x_continuous(trans = 'log', breaks = 10^(-2:3), labels = scales::number_format(accuracy = c(0.01, 0.1, 1, 1, 1, 1), big.mark = ','), 
                     sec.axis = sec_axis(~.^(-1)*1000, name = 'Timescale (years)', breaks = 10^(5:0), labels = scales::number_format(accuracy = c(1, 1, 1, 1, 1, 1), big.mark = ','))) +
  scale_y_continuous(trans = 'log', breaks = 10^(-3:1), labels = scales::number_format(accuracy = c(0.001, 0.01, 0.1, 1, 1), big.mark = ',')) +
  scale_colour_manual(name = 'Age (kya)', values = colfun1) +
  scale_size_manual(name = 'Noise level (permil)', values = c(0.5, 0.5), breaks = c('0.1', '0.01')) +
  labs(title = '(a)') +
  xlab(expression(paste('f'['max_decon'], ' (ky'^{-1}, ')'))) +
  ylab(expression(paste('Measurement noise level (permil)'))) +
  theme_bw() +
  theme(legend.position = "none",
        plot.tag = element_text())

fig_5b <- ggplot(data = fgain_df) +
  geom_line(aes(x = ft_des_rel, y = noise_req, colour = Age, linetype = 'solid')) +
  geom_hline(aes(yintercept = 0.1, size = '0.1'), colour = '#fb8500', linetype = 'dashed') +
  geom_hline(aes(yintercept = 0.01, size = '0.01'), colour = '#219ebc', linetype = 'dashed') +
  geom_vline(aes(xintercept = 1), linetype = 'dotted') +
  geom_point(aes(x = rep(1, length(fgain_df$ft_des)), y = H_y, colour = Age), shape = 19) +
  geom_line(aes(x = -2, y = 1, colour = Age, linetype = 'dotted')) +
  coord_cartesian(xlim=c(1, 7), ylim=c(1e-3, 10)) + #FOR A = 0.9
  scale_x_continuous(breaks = seq(1, 7, by = 1), labels = scales::number_format(accuracy = 1, big.mark = ',')) + #FOR A = 0.9
  scale_y_continuous(trans = 'log', breaks = 10^(-3:1), labels = scales::number_format(accuracy = c(0.001, 0.01, 0.1, 1, 1), big.mark = ',')) +
  scale_colour_manual(name = 'Age (kya)', values = colfun1) +
  scale_size_manual(name = 'Noise level (permil)', values = c(0.5, 0.5), breaks = c('0.1', '0.01')) +
  scale_linetype_manual(name = '', values = c('dotted' = 'dotted', 'solid' = 'solid'), labels = expression(f[max_diff], f[max_decon])) +
  guides(colour = guide_legend(order = 2),
         shape = guide_legend(order = 1),
         size = guide_legend(order = 0)) +
  labs(title = '(b)') +
  xlab(expression(Frequency~gain~(relative~to~f[max_diff]))) +
  ylab('Measurement noise level (permil)') +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.tag = element_text())

pdf(file = "./fig_5.pdf",
    width = 10,
    height = 5)

grid::grid.draw(cbind(ggplotGrob(fig_5a), ggplotGrob(fig_5b)))

dev.off()



# Diffusion length profile (Fig. B1) --------------------------------------

#Labels for x axis
diff_axis_ages <- c(0, 5, 20, 100, 300, 1500)
diff_axis_inds <- vector()

for(i in 1:length(diff_axis_ages)){
  
  diff_axis_inds[i] <- which.min(abs(LDC$age/1000 - diff_axis_ages[i]))
  
}

#Chosen indices for firn and deep ice cases
diff_sim_inds <- c(12, 858)


pdf(file = "./fig_B1.pdf",
    width = 6,
    height = 5)

plot(LDC$depth, s, 'l', xlab = 'Depth (m)', ylab = 'Diffusion length (m)', lwd = 2)
points(LDC$depth[diff_sim_inds[1]], s[diff_sim_inds[1]], pch = 19, col = 'blue')
lines(c(-100, LDC$depth[diff_sim_inds[1]]), rep(s[diff_sim_inds[1]], 2), col = 'blue', lwd = 2, lty = 2)
abline(v = LDC$depth[diff_sim_inds[1]], col = 'blue', lwd = 2, lty = 2)
points(LDC$depth[diff_sim_inds[2]], s[diff_sim_inds[2]], pch = 19, col = 'red')
lines(c(-100, LDC$depth[diff_sim_inds[2]]), rep(s[diff_sim_inds[2]], 2), col = 'red', lwd = 2, lty = 2)
abline(v = LDC$depth[diff_sim_inds[2]], col = 'red', lwd = 2, lty = 2)
axis(side = 3, at = LDC$depth[diff_axis_inds], labels = diff_axis_ages)
mtext('Age (kya)', side = 3, padj = -4)
legend(50, 0.2, c('diffusion length profile', 'shallow (white noise) value', 'deep (power law) value'), col = c('black', 'blue', 'red'), lwd = 2, lty = c(1, 2, 2), bty = 'n')

dev.off()

