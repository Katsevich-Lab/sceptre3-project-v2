gamma_left <- Vectorize(function(n, m, t, L) {
  x <- (n * t / (m / L + n) - (n * t)/(m + n))/(sqrt(t * m * n)/(m + n))
  pnorm(x)
})

n_hyp <- 10000
t <- MASS::rnegbin(n = n_hyp, mu = 800, theta = 1)
m <- rpois(n_hyp, 40)
n <- rpois(n_hyp, 500)
L <- 1/2
p <- gamma_left(n, m, t, L)
p[is.na(p)] <- 1
df <- data.frame(hyp_no = seq_len(n_hyp), p = p)
df <- df |> dplyr::arrange(p)

L <- 1/3
p <- gamma_left(n, m, t, L)
p[is.na(p)] <- 1
df_2 <- data.frame(hyp_no = seq_len(n_hyp), p = p)
df_2 <- df_2 |> dplyr::arrange(p)
