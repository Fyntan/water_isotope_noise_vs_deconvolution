
f_gain <- function(power_law, A_frac = 0.25, a = alpha, b = beta, s = sigma, n1 = noise_sd){
  
  f <- seq(0.001, 600, 0.001)
  
  N1 <- n1^2*dz
  
  if(power_law == TRUE){
    P_0 <- a*f^-b
  }else{
    P_0 <- a
  }
  
  H <- exp(-0.5*(2*pi*f*s)^2)
  G1 <- 1/H * (1/(1 + 1/(H^2*P_0/N1)))
  
  f_H <- f[max(which(H > 0.25))]
  
  P <- P_0*H^2
  
  f1 <- f[which.min(abs(G1*H - A_frac))]
  
  eta <- seq(1, 100, by = 1)
  
  N2 <- N1/eta
  
  f2 <- vector(length = length(eta))
  
  for(i in 1:length(eta)){
    G2 <- 1/H * (1/(1 + 1/(H^2*P_0/N2[i])))
    f2[i] <- f[which.min(abs(G2*H - A_frac))]
  }
  
  fgain <- f2/f1
  
  list(f = f, P = P, eta = eta, fgain = fgain, N1 = N1, N2 = N2, f1 = f1, f2 = f2, H = H, f_H = f_H)
}

dz <- 0.11
beta <- 1.4 #Not necessary, as this is white-noise case
sigma <- 0.2
noise_sd <- 0.1

#Three different alphas for comparison
alpha_low <- noise_sd^2*dz*5
alpha_mid <- noise_sd^2*dz*10
alpha_high <- noise_sd^2*dz*20

#Relative amplitude
A_frac <- 0.25

alpha_low_results <- f_gain(power_law = FALSE, A_frac = A_frac, a = alpha_low)
alpha_mid_results <- f_gain(power_law = FALSE, A_frac = A_frac, a = alpha_mid)
alpha_high_results <- f_gain(power_law = FALSE, A_frac = A_frac, a = alpha_high)

alpha_low_theory <- sqrt(1 + log(alpha_low_results$eta)/log((alpha_low*(1/A_frac - 1))/alpha_low_results$N1))
alpha_mid_theory <- sqrt(1 + log(alpha_mid_results$eta)/log((alpha_mid*(1/A_frac - 1))/alpha_mid_results$N1))
alpha_high_theory <- sqrt(1 + log(alpha_high_results$eta)/log((alpha_high*(1/A_frac - 1))/alpha_high_results$N1))




pdf(file = "./fig_4.pdf",
   width = 5,
   height = 4)

par(mfrow = c(1, 1), mar = c(5, 5, 1, 1))

fgain_colfun <- c('mediumorchid1',
                  'turquoise1',
                  'springgreen')

plot(sqrt(alpha_low_results$eta), alpha_low_results$fgain, 'l', lwd = 4, col = 'black',
     xlab = expression(paste('Noise SD reduction ', sqrt(eta))), ylab = expression(paste('Frequency gain ', f[gain])),
     cex.lab = 1.2, cex.axis = 1.2, ylim = range(alpha_low_results$fgain))
lines(sqrt(alpha_low_results$eta), alpha_low_theory, col = fgain_colfun[1], lwd = 4, lty = 2)
lines(sqrt(alpha_mid_results$eta), alpha_mid_results$fgain, lwd = 4, col = 'black')
lines(sqrt(alpha_mid_results$eta), alpha_mid_theory, col = fgain_colfun[2], lwd = 4, lty = 2)
lines(sqrt(alpha_high_results$eta), alpha_high_results$fgain, lwd = 4, col = 'black')
lines(sqrt(alpha_high_results$eta), alpha_high_theory, col = fgain_colfun[3], lwd = 4, lty = 2)
legend('topleft', c('5', '10', '20'), title = expression(paste(SNR[0])), col = fgain_colfun, lwd = 4, lty = c(1, 1, 1), bty = 'n')

dev.off()

