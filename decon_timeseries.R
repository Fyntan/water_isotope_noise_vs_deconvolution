
library(PaleoSpec)

#Function which applies Wiener deconvolution
decon <- function(z, rec, sigma, alpha, beta, noise){
  
  dz <- z[2] - z[1]
  N <- length(z)
  
  f <- seq(1/(N*dz), 1/(2*dz), 1/(N*dz))
  f <- c(0, f, -rev(f)[-1])
  
  H <- exp(-0.5 * (2*pi * f * sigma)^2)
  
  P0 <- alpha*abs(f)^-beta
  SNR <- P0 / (noise^2*dz)
  
  G <- ifelse(SNR == 0|abs(H)^2 == 0, 0, 1/H * 1 / (1 + 1 / (abs(H)^2 * SNR)) )
  
  rec_decon <- Re(fft(G * fft(rec), inverse = TRUE)/N)
  
  return(rec_decon)
}


# Load LDC data --------------------------------------------------------

LDC <- read_tsv("./LDC-VHF_updated.txt") #Updated version
LDC <- LDC[-(1:10), ] %>% 
  separate(names(LDC), into = c("depth",	"age",	"sigma_age",	"age_density",	"accumulation", "thinning",	"vertical_velocity"), sep = "\t")
LDC <- as.tibble(sapply(LDC, as.numeric))
LDC <- LDC[which(LDC$age < 1.5e6), ]

#LDC layer thickness
LDC_bdot <- diff(LDC$depth)/diff(LDC$age) #In metres per year
LDC_bdot <- LDC_bdot*1000 #In metres per ky

#Diffusion lengths
vas <- read_tsv("./vas_sigmas.txt")
s <- approx(vas$Depth, vas$sigma_total, LDC$depth)$y



# Firn case (Fig. 6) ------------------------------------------------------

dz <- 0.05
N <- 1000
z <- seq(dz, N*dz, dz)

shallow_ind <- 12

alpha_true <- epsilon_est*LDC_bdot[shallow_ind]
beta_true <- 0
sigma_true <- s[shallow_ind]
noise_true <- 0.1

seed <- 5738 #This seed was used for Fig. 6

set.seed(seed)
rec_ori <- SimProxySeries(a = alpha_true, b = beta_true, nt = N, f.scl = 1/dz, t.smpl = z)

set.seed(seed)
rec_diff_n1 <- SimProxySeries(a = alpha_true, b = beta_true, nt = N, f.scl = 1/dz, t.smpl = z, smth.arch = list(type = 'diffusion', tau = sigma_true), var.noise = noise_true^2)

set.seed(seed)
rec_diff_n2 <- SimProxySeries(a = alpha_true, b = beta_true, nt = N, f.scl = 1/dz, t.smpl = z, smth.arch = list(type = 'diffusion', tau = sigma_true), var.noise = (noise_true/10)^2)

rec_decon_n1 <- decon(z = z, rec = rec_diff_n1, sigma = sigma_true, alpha = alpha_true, beta = beta_true, noise = noise_true)

rec_decon_n2 <- decon(z = z, rec = rec_diff_n2, sigma = sigma_true, alpha = alpha_true, beta = beta_true, noise = (noise_true/10))

rec_anom_n1 <- abs(rec_ori - rec_decon_n1)
rec_anom_n2 <- abs(rec_ori - rec_decon_n2)


#Plotting

zlims <- c(LDC$depth[shallow_ind], LDC$depth[shallow_ind] + 5)

recs_df <- tibble(z = z + LDC$depth[shallow_ind] - 5, rec_ori = rec_ori, rec_diff_n1 = rec_diff_n1, rec_diff_n2 = rec_diff_n2, rec_decon_n1 = rec_decon_n1, rec_decon_n2 = rec_decon_n2, rec_anom_n1 = rec_anom_n1, rec_anom_n2 = rec_anom_n2)

fig_6a <- ggplot(recs_df) +
  geom_line(aes(z, rec_ori, col = 'original')) +
  geom_line(aes(z, rec_diff_n1, col = 'noise = 0.1')) +
  geom_line(aes(z, rec_diff_n2, col = 'noise = 0.01')) +
  coord_cartesian(xlim = zlims, ylim = range(rec_ori[which(recs_df$z > zlims[1] & recs_df$z < zlims[2])])) +
  scale_colour_manual(values = c('grey', '#fb8500', '#219ebc'), breaks = c('original', 'noise = 0.1', 'noise = 0.01')) +
  labs(title = '(a)') +
  xlab('Depth (m)') +
  ylab(expression(paste(delta^{18}, "O anomolies"))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        plot.tag = element_text())

fig_6b <- ggplot(recs_df) +
  geom_line(aes(z, rec_ori, col = 'original')) +
  geom_line(aes(z, rec_decon_n1, col = 'noise = 0.1')) +
  geom_line(aes(z, rec_decon_n2, col = 'noise = 0.01')) +
  coord_cartesian(xlim = zlims, ylim = range(rec_ori[which(recs_df$z > zlims[1] & recs_df$z < zlims[2])])) +
  scale_colour_manual(name = '', values = c('grey', '#fb8500', '#219ebc'), breaks = c('original', 'noise = 0.1', 'noise = 0.01')) +
  labs(title = '(b)') +
  xlab('Depth (m)') +
  ylab(expression(paste(delta^{18}, "O anomolies"))) +
  theme_bw() +
  theme(axis.title.x=element_blank())

fig_6c <- ggplot(recs_df) +
  geom_line(aes(z, rec_anom_n1, col = 'noise = 0.1')) +
  geom_line(aes(z, rec_anom_n2, col = 'noise = 0.01')) +
  coord_cartesian(xlim = zlims) +
  scale_colour_manual(values = c('#fb8500', '#219ebc'), breaks = c('noise = 0.1', 'noise = 0.01')) +
  labs(title = '(c)') +
  xlab('Depth (m)') +
  ylab(expression(paste(delta^{18}, "O residuals"))) +
  theme_bw() +
  theme(legend.position = "none",
        plot.tag = element_text())

pdf(file = "./fig_6.pdf",
    width = 8,
    height = 5)

grid::grid.draw(rbind(ggplotGrob(fig_6a), ggplotGrob(fig_6b), ggplotGrob(fig_6c)))

dev.off()



# Deep ice case (Fig. 7) --------------------------------------------------

dz <- 0.05
N <- 250
z <- seq(dz, N*dz, dz)

deep_ind <- 858

alpha_true <- alpha_est*LDC_bdot[deep_ind]^(1-beta_est)
beta_true <- beta_est
sigma_true <- s[deep_ind]
noise_true <- 0.1

seed <- 5241 #This seed was used for Fig. 6

set.seed(seed)
rec_ori <- SimProxySeries(a = alpha_true, b = beta_true, nt = N, f.scl = 1/dz, t.smpl = z)

set.seed(seed)
rec_diff_n1 <- SimProxySeries(a = alpha_true, b = beta_true, nt = N, f.scl = 1/dz, t.smpl = z, smth.arch = list(type = 'diffusion', tau = sigma_true), var.noise = noise_true^2)

set.seed(seed)
rec_diff_n2 <- SimProxySeries(a = alpha_true, b = beta_true, nt = N, f.scl = 1/dz, t.smpl = z, smth.arch = list(type = 'diffusion', tau = sigma_true), var.noise = (noise_true/10)^2)

rec_decon_n1 <- decon(z = z, rec = rec_diff_n1, sigma = sigma_true, alpha = alpha_true, beta = beta_true, noise = noise_true)

rec_decon_n2 <- decon(z = z, rec = rec_diff_n2, sigma = sigma_true, alpha = alpha_true, beta = beta_true, noise = (noise_true/10))

rec_anom_n1 <- abs(rec_ori - rec_decon_n1)
rec_anom_n2 <- abs(rec_ori - rec_decon_n2)


#Plotting

zlims <- c(LDC$depth[deep_ind], LDC$depth[deep_ind] + 5)

recs_df <- tibble(z = z + LDC$depth[deep_ind] - 5, rec_ori = rec_ori, rec_diff_n1 = rec_diff_n1, rec_diff_n2 = rec_diff_n2, rec_decon_n1 = rec_decon_n1, rec_decon_n2 = rec_decon_n2, rec_anom_n1 = rec_anom_n1, rec_anom_n2 = rec_anom_n2)

fig_7a <- ggplot(recs_df) +
  geom_line(aes(z, rec_ori, col = 'original')) +
  geom_line(aes(z, rec_diff_n1, col = 'noise = 0.1')) +
  geom_line(aes(z, rec_diff_n2, col = 'noise = 0.01')) +
  coord_cartesian(xlim = zlims, ylim = range(rec_ori[which(recs_df$z > zlims[1] & recs_df$z < zlims[2])])) +
  scale_colour_manual(values = c('grey', '#fb8500', '#219ebc'), breaks = c('original', 'noise = 0.1', 'noise = 0.01')) +
  labs(title = '(a)') +
  xlab('Depth (m)') +
  ylab(expression(paste(delta^{18}, "O anomolies"))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        plot.tag = element_text())

fig_7b <- ggplot(recs_df) +
  geom_line(aes(z, rec_ori, col = 'original')) +
  geom_line(aes(z, rec_decon_n1, col = 'noise = 0.1')) +
  geom_line(aes(z, rec_decon_n2, col = 'noise = 0.01')) +
  coord_cartesian(xlim = zlims, ylim = range(rec_ori[which(recs_df$z > zlims[1] & recs_df$z < zlims[2])])) +
  scale_colour_manual(name = '', values = c('grey', '#fb8500', '#219ebc'), breaks = c('original', 'noise = 0.1', 'noise = 0.01')) +
  labs(title = '(b)') +
  xlab('Depth (m)') +
  ylab(expression(paste(delta^{18}, "O anomolies"))) +
  theme_bw() +
  theme(axis.title.x=element_blank())

fig_7c <- ggplot(recs_df) +
  geom_line(aes(z, rec_anom_n1, col = 'noise = 0.1')) +
  geom_line(aes(z, rec_anom_n2, col = 'noise = 0.01')) +
  coord_cartesian(xlim = zlims) +
  scale_colour_manual(values = c('#fb8500', '#219ebc'), breaks = c('noise = 0.1', 'noise = 0.01')) +
  labs(title = '(c)') +
  xlab('Depth (m)') +
  ylab(expression(paste(delta^{18}, "O residuals"))) +
  theme_bw() +
  theme(legend.position = "none",
        plot.tag = element_text())

pdf(file = "./fig_7.pdf",
    width = 8,
    height = 5)

grid::grid.draw(rbind(ggplotGrob(fig_7a), ggplotGrob(fig_7b), ggplotGrob(fig_7c)))

dev.off()
