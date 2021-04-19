library(here)
library(tidyverse)
library(wesanderson)
library(patchwork)
library(animation)
set.seed(34729)
pal <- wes_palette("Darjeeling1")

n_obs <- 1000
delta_grid <- seq(-2, 2, 0.5)
delta_grid <- c(0, delta_grid)
w <- rbinom(n_obs, 1, 0.5)
a <- rnorm(n_obs, mean = 0.5 * w, sd = 0.5)

qbar_and_shift <- lapply(delta_grid, function(delta) {
  a_shifted <- a + delta
  qbar_aplusdelta <- plogis(-3 - 3 * a_shifted + w)
  return(as_tibble(cbind(a_shifted, qbar_aplusdelta)))
})

saveGIF(
  {
    for (iter in seq_along(qbar_and_shift)) {
      df_input <- qbar_and_shift[[iter]]
      a_mean <- mean(df_input$a_shifted)
      p_a_shifted <- df_input %>%
        ggplot(aes(x = a_shifted)) +
        geom_histogram(
          colour = "white", fill = pal[[5]], alpha = 0.6,
          binwidth = 0.2
        ) +
        geom_vline(
          xintercept = a_mean, size = 0.5, linetype = "dashed",
          colour = "black"
        ) +
        xlim(-4, 4) +
        xlab("") +
        ylab("") +
        ggtitle("Exposure profile") +
        theme_bw() +
        theme(
          text = element_text(size = 22),
          axis.text.x = element_text(colour = "black", size = 22),
          axis.text.y = element_text(colour = "black", size = 22)
        )

      qbar_mean <- mean(df_input$qbar_aplusdelta)
      p_qbar_shifted <- df_input %>%
        ggplot(aes(y = qbar_aplusdelta)) +
        geom_hline(
          yintercept = qbar_mean, size = 2, linetype = "dashed",
          colour = "black"
        ) +
        xlab("") +
        ylab("") +
        ylim(0, 1) +
        ggtitle("Counterfactual mean risk") +
        theme_bw() +
        theme(
          text = element_text(size = 22),
          axis.text.x = element_text(colour = "black", size = 22),
          axis.text.y = element_text(colour = "black", size = 22)
        )

      p_out <- p_a_shifted + p_qbar_shifted
      print(p_out)
    }
  },
  movie.name = here("img", "gif", "shift_animation.gif"),
  interval = 1,
  ani.width = 1400,
  ani.height = 1000
)
