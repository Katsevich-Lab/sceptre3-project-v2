library(ggplot2)
n <- 10
m <- 50
mean_exp <- 0.5
L <- 1/4
t <- mean_exp * (n + m)

# the range of Y is 0 to tot_umi_count; either none of the UMIs or all of the UMIs are contained within the treatment cells
mu <- (t * n)/(m + n)
sigma <- sqrt((t * m * n)/((m + n)^2))
get_extreme_y <- function(n, m, t, L = L) {
  n * t /(n + m / L)
}

# plot the normal approximation to the null distribution, alongside the minimum and maximum values of y (n_umis_trt)
extreme_decrease <- get_extreme_y(n, m, t, L = L)
extreme_increase <- get_extreme_y(n, m, t, L = 1/L)
my_min <- min(qnorm(p = 1e-4, mean = mu, sd = sigma), extreme_decrease)
my_max <- max(qnorm(p = 1 - 1e-4, mean = mu, sd = sigma), extreme_increase)

x_grid <- seq(my_min, my_max, length.out = 1000)
dens <- dnorm(x = x_grid, mean = mu, sd = sigma)
df <- data.frame(x_grid = x_grid, dens = dens)
no_change <- get_extreme_y(n, m, t, L = 1)
ggplot(data = df, mapping = aes(x = x_grid, y = dens)) +
  geom_line() +
  geom_vline(xintercept = extreme_decrease, col = "darkblue") +
  geom_vline(xintercept = extreme_increase, col = "darkblue") +
  geom_vline(xintercept = no_change, col = "darkred") +
  theme_bw() +
  xlab("Y")

# compute the smallest possible left-tailed p-value
gamma_left <- function(n, m, t, L) {
  x <- (n * t / (m / L + n) - (n * t)/(m + n))/(sqrt(t * m * n)/(m + n))
  pnorm(x)
}
gamma_right <- function(n, m, t, L) {
  x <- (n * t / (m / L + n) - (n * t)/(m + n))/(sqrt(t * m * n)/(m + n))
  pnorm(x, lower.tail = FALSE)
}
gamma_both <- Vectorize(function(n, m, t, L) {
  min(1, 2 * min(gamma_left(n, m, t, L), gamma_right(n, m, t, 1/L)))
})


# If we are doing a Bonf. adjustment, this is probably how we would select the cutoff; order the hypotheses according to their candidate p-values, then select the K most promising hypotheses according to some rule.
# suppose m, n fixed across hypotheses (as in bulk rna seq), but t varies
n_hyp <- 500
t <- MASS::rnegbin(n = n_hyp, mu = 800, theta = 1)
m <- 40
n <- 950
L <- 2

# compute pilot alternative p-values under a given fold change
# p <- gamma_right(n = n, m = m, t = t, L = L)
# p <- gamma_left(n = n, m = m, t = t, L = 1/L)
p <- gamma_both(n = n, m = m, t = t, L = 1/L)
p[is.na(p)] <- 1
p_sort <- sort(p, decreasing = FALSE)
for (i in seq_along(p_sort)) {
  if (p_sort[i] >= alpha/i) break
}
thresh <- alpha / i
ggplot(data = data.frame(x_grid = seq_along(p_sort),
                              p = p_sort,
                              retain = seq_along(p_sort) <= i),
       mapping = aes(x = x_grid, y = p_sort, col = retain)) +
  geom_point(size = 0.1) +
  theme_bw() +
  scale_y_continuous(trans = sceptre:::revlog_trans(10)) +
  ylab("Pilot p-value (under alternative)") +
  xlab("Rank") +
  geom_function(fun = function(x) alpha/x, color = "black")

# If we are doing a BH adjustment, then it's more complicated: we do not know the strength of the signal.
# optimistic: fold change is large, 
pi <- 0.05
alpha <- 0.1
pi * alpha /(1 - alpha)
