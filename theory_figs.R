#Quick simulated ts and power spectrum

library(PaleoSpec)
library(tidyverse)

alpha <- 1
beta <- 1
sigma <- 3
noise <- 0.2

N <- 1000
dz <- 0.1

diff_theory <- SimProxySeries(a = alpha, b = beta, nt = N, smth.arch = list(type = 'diffusion', tau = sigma), val = 1)
noise_1 <- 0.07
noise_2 <- 0.0007

pdf(file = "./fig_1.pdf",
    width = 6,
    height = 4)

par(mar = c(1, 1, 1, 16))

plot(diff_theory$fax/dz, diff_theory$psd*dz, 'l', log = 'xy', lwd = 2, ylim = c(0.0000001, 10), xlim = c(0.2, 2),
     xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
abline(h = noise_1^2, col = '#fb8500', lwd = 2, lty = 2)
abline(h = noise_2^2, col = '#219ebc', lwd = 2, lty = 2)
abline(v = diff_theory$fax[which.min(abs(diff_theory$psd*dz - noise_1^2))]/dz, col = '#fb8500', lwd = 2)
abline(v = diff_theory$fax[which.min(abs(diff_theory$psd*dz - noise_2^2))]/dz, col = '#219ebc', lwd = 2)

legend(2.2, 0.03, c('diffused record', 'initial noise', 'initial frequency limit', 'reduced noise', 'improved frequency limit'), col = c('black', '#fb8500', '#fb8500', '#219ebc', '#219ebc'), xpd = TRUE, lty = c(1, 2, 1, 2, 1), lwd = 2, bty = 'n')

dev.off()




alpha <- 1
beta <- 1
sigma <- 3
noise <- 0.07^2

freq <- seq(0.001, 0.5, 0.001)
P_0 <- alpha*freq^-beta

H <- exp(-0.5*(2*pi*freq*sigma)^2)
G <- 1/H * (1/(1 + 1/(H^2*P_0/noise)))

frac <- 0.25

def_SNR <- max(which(H^2*P_0/noise > 1))
def_GH <- max(which(G*H > frac))
def_H <- max(which(H > frac))

pdf("./fig_3.pdf",
    width = 10,
    height = 4.3)

par(mfrow = c(1,2))

plot(freq, P_0, 'l', log = 'xy', col = 'grey', lwd = 2, ylim = c(1e-7, 1e4),
     xlab = 'Frequency (arbitrary units)', ylab = 'PSD (arbitrary units)')
lines(freq, P_0*H^2, col = 'green4', lwd = 2)
lines(freq, P_0*H^2*G^2, col = 'blue', lwd = 2)
abline(h = noise, col = 'red3', lwd = 2)
points(freq[def_SNR], noise, pch = 19, col = 'green4')
legend('bottomleft', c(expression(P[0]), expression(noise~P[n]), expression(diffused~P[0]~(P[0]*H^2)), expression(deconvolved~P[0]~(P[0]*H^2*G^2))), col = c('grey', 'red3', 'green4', 'blue'), lty = 1, lwd = 2, bty = 'n')             
legend('topleft', 'SNR = 1', col = 'green4', pch = 19, bty = 'n')
title('(a)', adj = 0)

plot(freq, H, 'l', log = 'xy', col = 'orange', lwd = 2, ylim = c(1e-5, 1e2),
     xlab = 'Frequency (arbitrary units)', ylab = 'Relative amplitude (unitless)')
abline(h = 0.25, lwd = 2, col = 'grey30', lty = 2)
abline(v = freq[def_GH], col = 'skyblue3', lwd = 2, lty = 2)
abline(v = freq[def_H], col = 'orange', lwd = 2, lty = 2)
lines(freq, G, col = 'purple2', lwd = 2)
lines(freq, G*H, col = 'skyblue3', lwd = 2)
points(freq[def_GH], frac, pch = 19, col = 'skyblue3')
points(freq[def_H], frac, pch = 19, col = 'orange')
legend('bottomleft', c(expression(A[diff]~(H)), 'G', expression(A[decon]~(HG)), expression(f[max_diff]), expression(f[max_decon]), expression(A[min]~'='~0.25)), col = c('orange', 'purple2', 'skyblue3', 'orange', 'skyblue3', 'grey30'), lty = c(1, 1, 1, 1, 2, 2, 2), lwd = 2, bty = 'n')
title('(b)', adj = 0)

dev.off()

