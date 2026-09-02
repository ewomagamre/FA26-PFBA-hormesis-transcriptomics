
suppressPackageStartupMessages({
  library(ggplot2)
})

dir.create("results/phenotype", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/main", recursive = TRUE, showWarnings = FALSE)

dat <- data.frame(
  Treatment = rep(c("Control", "PFBA 0.01", "PFBA 1"), each = 12),
  Time_hr = rep(
    c(rep(2,4), rep(4,4), rep(6,4)),
    3
  ),
  Replicate = rep(1:4, 9),
  Weight_change = c(
    20,36,40,21,
    4,-3,-1,33,
    16,19,-6,-6,

    1,27,35,6,
    13,5,18,3,
    10,3,26,26,

    21,24,61,22,
    31,48,33,20,
    40,37,12,27
  )
)

write.csv(
  dat,
  "results/phenotype/FA26_weight_change.csv",
  row.names = FALSE
)

groups <- split(
  dat,
  interaction(dat$Time_hr, dat$Treatment, drop = TRUE)
)

summary_list <- lapply(groups, function(x) {
  data.frame(
    Time_hr = x$Time_hr[1],
    Treatment = x$Treatment[1],
    n = nrow(x),
    mean = mean(x$Weight_change),
    sd = sd(x$Weight_change),
    sem = sd(x$Weight_change) / sqrt(nrow(x))
  )
})

sumdat <- do.call(rbind, summary_list)
rownames(sumdat) <- NULL

sumdat$Treatment <- factor(
  sumdat$Treatment,
  levels = c("Control", "PFBA 0.01", "PFBA 1")
)

sumdat <- sumdat[
  order(sumdat$Time_hr, sumdat$Treatment),
]

write.csv(
  sumdat,
  "results/phenotype/FA26_weight_change_summary.csv",
  row.names = FALSE
)

cols <- c(
  "Control"   = "#E69F00",
  "PFBA 0.01" = "#1675B9",
  "PFBA 1"    = "#E6007E"
)

pd <- position_dodge(width = 0.78)

p <- ggplot(
  sumdat,
  aes(
    x = factor(Time_hr),
    y = mean,
    fill = Treatment
  )
) +

  geom_hline(
    yintercept = 0,
    linewidth = 0.45,
    linetype = "dashed",
    colour = "grey55"
  ) +

  geom_col(
    position = pd,
    width = 0.68,
    colour = "black",
    linewidth = 0.35
  ) +

  geom_errorbar(
    aes(
      ymin = mean - sem,
      ymax = mean + sem
    ),
    position = pd,
    width = 0.16,
    linewidth = 0.65
  ) +

  scale_fill_manual(
    values = cols,
    labels = c(
      "Control",
      "PFBA 0.01 \u00b5g/g",
      "PFBA 1 \u00b5g/g"
    )
  ) +

  labs(
    x = "Time after exposure (h)",
    y = "Weight change",
    fill = NULL
  ) +

  theme_classic(base_size = 15) +

  theme(
    axis.title = element_text(
      face = "bold",
      size = 16
    ),
    axis.text = element_text(
      size = 13,
      colour = "black"
    ),
    axis.line = element_line(
      linewidth = 0.7,
      colour = "black"
    ),
    axis.ticks = element_line(
      linewidth = 0.6
    ),
    legend.position = "top",
    legend.text = element_text(
      size = 13
    ),
    plot.margin = margin(
      12, 15, 12, 12
    )
  )

ggsave(
  "figures/main/Figure1_WeightPhenotype.pdf",
  plot = p,
  width = 7.2,
  height = 5.4,
  units = "in"
)

ggsave(
  "figures/main/Figure1_WeightPhenotype.tiff",
  plot = p,
  width = 7.2,
  height = 5.4,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

cat("\nFigure 1 updated successfully.\n")
cat("figures/main/Figure1_WeightPhenotype.pdf\n")
cat("figures/main/Figure1_WeightPhenotype.tiff\n")

