# =============================================================================
# _common.R  —  Shared setup for the T-helper / exercise meta-analysis
# Sourced by every chapter of the Quarto book AND by the single-file qmd.
# Holds: package loading, global ggplot theme, the leave-one-out outlier
# functions, a project-root resolver, a data loader, and styled-plot helpers
# so the analysis chapters stay short and consistent.
# =============================================================================

## ---- Project root (works from the book/ dir or from the project root) ------
proj_root <- if (file.exists("effect_size")) "." else ".."
es_path   <- function(...) file.path(proj_root, ...)

## ---- Packages --------------------------------------------------------------
.required <- c(
  "metafor", "tidyr", "dplyr", "data.table", "stringr", "readr",
  "ggplot2", "reshape2", "cowplot", "gridExtra",
  "forestplot", "metaviz", "meta",
  "readxl", "patchwork", "scales", "kableExtra"
)
.missing <- .required[!(.required %in% rownames(installed.packages()))]
if (length(.missing)) {
  install.packages(.missing, repos = "https://cloud.r-project.org")
}

invisible(suppressWarnings(suppressMessages(
  lapply(.required, require, character.only = TRUE)
)))

## ---- Global ggplot2 theme (matches the document styling) -------------------
ggplot2::theme_set(
  ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      text             = ggplot2::element_text(color = "#1E1E1E"),
      plot.title       = ggplot2::element_text(face = "bold", size = 13, margin = ggplot2::margin(b = 8)),
      plot.subtitle    = ggplot2::element_text(size = 11, color = "grey40"),
      axis.title       = ggplot2::element_text(size = 11),
      axis.text        = ggplot2::element_text(size = 10, color = "#1E1E1E"),
      legend.text      = ggplot2::element_text(size = 10),
      legend.title     = ggplot2::element_text(size = 10, face = "bold"),
      legend.position  = "bottom",
      panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#f0f4f8", color = NA),
      strip.text       = ggplot2::element_text(face = "bold", size = 10)
    )
)

## ---- Data loader -----------------------------------------------------------
# Reads a pre-computed effect-size CSV (";"-separated, period decimals) and
# maps the pipeline column names to the conventions used in the analyses.
# Only columns that exist in a given file are coerced (the +2h files, for
# example, carry no cytokine columns).
load_es <- function(rel_path) {
  d <- read.csv(es_path(rel_path), header = TRUE, sep = ";", stringsAsFactors = FALSE)
  d <- as.data.frame(d)
  num <- function(x) suppressWarnings(as.numeric(x))
  d$Short_reference <- as.character(d$Study)
  d$g    <- num(d$g)
  d$V_g  <- num(d$V_g)
  d$VE_g <- num(d$V_g)        # variance, kept under both names used in the script
  d$SE_g <- num(d$SE_g)
  for (col in c("age", "bmi", "duration", "dose", "met",
                "IFN", "IL_10", "IL_2", "IL_4", "IL_6")) {
    if (col %in% names(d)) d[[col]] <- num(d[[col]])
  }
  for (col in c("RoB", "intensity", "sex", "staining")) {
    if (col %in% names(d)) d[[col]] <- as.character(d[[col]])
  }
  # Capitalised aliases used in prose / forest labels
  d$Sex <- if ("sex" %in% names(d)) as.character(d$sex) else NA_character_
  d$Age <- if ("age" %in% names(d)) num(d$age) else NA_real_
  d$BMI <- if ("bmi" %in% names(d)) num(d$bmi) else NA_real_
  d
}

# Convenience: drop the study flagged as unrealistic in the main analyses.
drop_kostrzewa19 <- function(d) {
  dplyr::filter(d, !Study == "Kostrzewa-Nowak et al. (2019)")
}

## ---- Leave-one-out outlier detection (Hedges & Olkin 1985; Viechtbauer &
##      Cheung 2010). A study is an outlier if |std. residual| >= 3. ----------
metaoutliers <- function(y, s2, model) {
  if (length(y) != length(s2) | any(s2 < 0)) stop("error in the input data.")
  w <- 1 / s2
  y.p <- sum(y * w) / sum(w)
  n <- length(y)

  if (missing(model)) {
    hetmeasure <- metahet.base(y, s2)
    Ir2 <- hetmeasure$Ir2
    if (Ir2 < 0.3) {
      model <- "FE"
      cat("This function uses fixed-effect meta-analysis because Ir2 < 30%.\n")
    } else {
      model <- "RE"
      cat("This function uses random-effects meta-analysis because Ir2 >= 30%.\n")
    }
  }

  if (!is.element(model, c("FE", "RE"))) stop("wrong input for the argument model.")

  y.p.i <- res <- std.res <- numeric(n)
  if (model == "FE") {
    for (i in 1:n) {
      w.temp <- w[-i]
      y.temp <- y[-i]
      y.p.i[i] <- sum(y.temp * w.temp) / sum(w.temp)
      res[i] <- y[i] - y.p.i[i]
      var.res.i <- 1 / sum(w.temp) + s2[i]
      std.res[i] <- res[i] / sqrt(var.res.i)
    }
  } else {
    for (i in 1:n) {
      s2.temp <- s2[-i]
      y.temp <- y[-i]
      tau2.temp <- metahet.base(y.temp, s2.temp)$tau2.DL
      w.temp <- 1 / (s2.temp + tau2.temp)
      y.p.i[i] <- sum(y.temp * w.temp) / sum(w.temp)
      res[i] <- y[i] - y.p.i[i]
      var.res.i <- 1 / sum(w.temp) + s2[i] + tau2.temp
      std.res[i] <- res[i] / sqrt(var.res.i)
    }
  }

  outliers <- which(abs(std.res) >= 3)
  if (length(outliers) == 0) outliers <- "All the standardized residuals are smaller than 3"

  out <- NULL
  out$model <- model
  out$std.res <- std.res
  out$outliers <- outliers
  class(out) <- "metaoutliers"
  return(out)
}

metahet.base <- function(y, s2) {
  if (length(y) != length(s2) | any(s2 < 0)) stop("error in the input data.")
  n <- length(y)
  w <- 1 / s2
  mu.bar <- sum(w * y) / sum(w)
  out <- NULL
  out$weighted.mean <- mu.bar

  # conventional methods
  Q <- sum(w * (y - mu.bar)^2)
  H <- sqrt(Q / (n - 1))
  I2 <- (Q - n + 1) / Q
  tau2.DL <- (Q - n + 1) / (sum(w) - sum(w^2) / sum(w))
  tau2.DL <- max(c(0, tau2.DL))
  out$Q <- Q; out$H <- H; out$I2 <- I2; out$tau2.DL <- tau2.DL

  # absolute deviation based on weighted mean
  Qr <- sum(sqrt(w) * abs(y - mu.bar))
  Hr <- sqrt((3.14159 * Qr^2) / (2 * n * (n - 1)))
  Ir2 <- (Qr^2 - 2 * n * (n - 1) / 3.14159) / (Qr^2)
  tau2.r <- tau2.r.solver(w, Qr)
  out$Qr <- Qr; out$Hr <- Hr; out$Ir2 <- Ir2; out$tau2.r <- tau2.r

  # absolute deviation based on weighted median
  expit <- function(x) ifelse(x >= 0, 1 / (1 + exp(-x / 0.0001)), exp(x / 0.0001) / (1 + exp(x / 0.0001)))
  psi <- function(x) sum(w * (expit(x - y) - 0.5))
  mu.med <- uniroot(psi, c(min(y) - 0.001, max(y) + 0.001))$root
  out$weighted.median <- mu.med
  Qm <- sum(sqrt(w) * abs(y - mu.med))
  Hm <- sqrt(3.14159 / 2) * Qm / n
  Im2 <- (Qm^2 - 2 * n^2 / 3.14159) / Qm^2
  tau2.m <- tau2.m.solver(w, Qm)
  out$Qm <- Qm; out$Hm <- Hm; out$Im2 <- Im2; out$tau2.m <- tau2.m
  return(out)
}

tau2.r.solver <- function(w, Qr) {
  f <- function(tau2) {
    sum(sqrt(1 - w / sum(w) + tau2 * (w - 2 * w^2 / sum(w) + w * sum(w^2) / (sum(w))^2))) - Qr * sqrt(3.14159 / 2)
  }
  f <- Vectorize(f)
  tau.upp <- Qr * sqrt(3.14159 / 2) / sum(sqrt(w - 2 * w^2 / sum(w) + w * sum(w^2) / (sum(w))^2))
  tau2.upp <- tau.upp^2
  if (f(0) * f(tau2.upp) > 0) 0 else uniroot(f, interval = c(0, tau2.upp))$root
}

tau2.m.solver <- function(w, Qm) {
  f <- function(tau2) sum(sqrt(1 + w * tau2)) - Qm * sqrt(3.14159 / 2)
  f <- Vectorize(f)
  n <- length(w)
  tau2.upp <- sum(1 / w) * (Qm^2 / n * 2 / 3.14159 - 1)
  tau2.upp <- max(c(tau2.upp, 0.01))
  if (f(0) * f(tau2.upp) > 0) 0 else uniroot(f, interval = c(0, tau2.upp))$root
}

## ---- PDF output: hide code (code-folding is HTML-only) & narrow printing ---
if (requireNamespace("knitr", quietly = TRUE)) {
  if (isTRUE(knitr::is_latex_output())) {
    options(width = 60)
    knitr::opts_chunk$set(echo = FALSE)
  }
}

## ---- Styled-plot helpers (consistent look across all chapters) -------------

# Contour-enhanced funnel plot with the significance-region labels.
# The labels are right-aligned to an individual x-value at the right edge:
#   - label_x : the x-value the labels are pinned to (default = right plot edge).
#   - label_y : top y-value where the stack starts (default = near the top edge).
#   - label_cex : text size.
# Pass label_x = <number> to place them at a specific Hedges' g value.
funnel_ce <- function(m, ref = m$TE.random, xlim = NULL, ylim = NULL,
                      label_x = NULL, label_y = NULL, label_cex = 0.9) {
  meta::funnel(
    m, random = TRUE, xlim = xlim, ylim = ylim,
    xlab = "Hedges' g", ylab = "Standard error", level = 0.95,
    contour = c(0.9, 0.95, 0.99),
    shade = c("white", "gray", "darkgray"),
    col = "black", bg = "darkgray", pch = 16, ref = ref, lwd = 1, cex = 1.2
  )

  # Plotting region in user coordinates: c(x1, x2, y1, y2)
  # (x1,y1) = bottom-left, (x2,y2) = top-right — works for the reversed SE axis too.
  usr <- par("usr")
  if (is.null(label_x)) label_x <- usr[2] - 0.02 * (usr[2] - usr[1])  # just inside the right edge
  if (is.null(label_y)) label_y <- usr[4] + 0.04 * (usr[3] - usr[4])  # just below the top edge

  legend(
    x = label_x, y = label_y,
    xjust = 1, yjust = 1,                       # anchor the box's top-RIGHT corner at (label_x, label_y)
    legend = c("0.1 > p > 0.05", "0.05 > p > 0.01", "p < 0.01"),
    fill   = c("white", "gray", "darkgray"),
    border = "black", bty = "n", cex = label_cex,
    x.intersp = 0.8, y.intersp = 1
  )
}

# Main random-effects forest plot (RevMan5 layout, document palette).
forest_main <- function(m, xlim = c(-2, 5),
                        left = "Lower at post", right = "Higher at post") {
  meta::forest(
    m, sortvar = -TE, layout = "RevMan5",
    common = FALSE, random = TRUE, prediction = FALSE,
    label.left = left, label.right = right, xlab = "Hedges' g",
    fontsize = 8, xlim = xlim, method.tau = "HE", weight.study = "same",
    digits.mean = 2, digits.sd = 2, digits.tau2 = 2, digits.I2 = 1,
    col.study = "black", col.square = "#1f4e79", col.square.lines = "#1f1f1f",
    col.diamond = "#7fb3d5", col.diamond.random = "#7fb3d5",
    col.diamond.lines = "#1f1f1f", col.label = "black", col.lines = "grey40"
  )
}

# Subgroup forest plot (RevMan5 layout, document palette).
forest_sub <- function(m, xlim = c(-2, 5),
                       left = "Lower at post", right = "Higher at post") {
  meta::forest(
    m, sortvar = -TE, layout = "RevMan5",
    common = FALSE, random = TRUE, prediction = FALSE, subgroup = TRUE,
    label.left = left, label.right = right, xlab = "Hedges' g",
    xlim = xlim, fontsize = 8, spacing = 1.2, colgap = "4mm",
    colgap.forest.left = "10mm", colgap.forest.right = "3mm", plotwidth = "9cm",
    method.tau = "HE", weight.study = "same",
    rightcols = c("effect", "ci"), print.tau2 = TRUE, digits.tau2 = 2, digits.I2 = 1,
    ff.random = "bold", print.subgroup.labels = TRUE, subgroup.name = NULL,
    test.subgroup.random = TRUE, test.effect.subgroup.random = TRUE,
    col.study = "black", col.square = "#1f4e79", col.square.lines = "#1f1f1f",
    col.diamond = "#7fb3d5", col.diamond.random = "#7fb3d5",
    col.diamond.lines = "#1f1f1f", col.label = "black", col.lines = "grey40",
    col.subgroup = "#1f4e79"
  )
}

# Meta-regression bubble plot (metafor::regplot, document palette).
bubble_reg <- function(mr, mod, xlab) {
  metafor::regplot(
    mr, mod = mod, xlab = xlab, ylab = "Hedges' g", ci = TRUE,
    shade = grDevices::adjustcolor("#7fb3d5", alpha.f = 0.20), pi = FALSE,
    pch = 21, col = "#1f1f1f", bg = grDevices::adjustcolor("#7fb3d5", alpha.f = 0.45),
    lcol = c("#1f4e79", "#4f81bd", NA, "grey75"),
    lwd = c(3, 1.5, NA, 1.2), lty = c(1, 2, NA, 2), refline = 0, grid = "grey92"
  )
}

# ---- Nicely formatted table (HTML + PDF compatible) ------------------------
# Drop-in replacement for printing a data.frame. Uses knitr::kable(), which
# auto-targets HTML or LaTeX depending on the output format; in HTML it adds
# light kableExtra styling, in PDF it stays clean booktabs (needs \usepackage
# {booktabs}, already loaded in the book header).
#
#   x         : a data.frame / tibble (works with the pipe).
#   digits    : decimal places for numeric columns (default 3).
#   caption   : optional table caption.
#   col_names : optional nicer column headers (defaults to the data's names).
#   align     : column alignment, e.g. "l", "c", or a per-column string.
#
# Usage:
#   fullData |> select(Short_reference, g) |>
#     mutate(std_res = outlier$std.res) |>
#     nice_table(caption = "Standardised residuals")
nice_table <- function(x, digits = 3, caption = NULL, col_names = NULL,
                       align = "l") {
  x <- as.data.frame(x)
  if (is.null(col_names)) col_names <- colnames(x)

  tab <- knitr::kable(
    x,
    digits    = digits,
    caption   = caption,
    col.names = col_names,
    align     = align,
    booktabs  = TRUE,
    linesep   = ""        # no extra space every 5 rows in LaTeX
  )

  # Extra styling only in HTML (keeps PDF as plain booktabs).
  if (knitr::is_html_output() && requireNamespace("kableExtra", quietly = TRUE)) {
    tab <- kableExtra::kable_styling(
      tab,
      bootstrap_options = c("striped", "hover", "condensed", "responsive"),
      full_width = FALSE,
      position   = "left",
      font_size  = 14
    )
  }
  tab
}

# ---- Standardised-residual plot (outlier diagnostics) ----------------------
# ggplot replacement for `plot(stdres, type = "o", ...)`.
# Draws the connected residual series with RED dashed threshold lines at
# +/- `threshold` (default 3); points with |z| >= threshold are highlighted
# in red. Optional grey dotted "inspection" lines at +/- `inspect` (1.96).
# Returns a ggplot -> works in HTML and PDF, and can go into patchwork.
#
#   stdres    : numeric vector of standardised residuals (outlier$std.res).
#   threshold : outlier cut-off (default 3).
#   inspect   : inspection cut-off, NULL to hide (default 1.96).
#   labels    : optional study labels; flagged points (|z| >= threshold) get
#               a text label.
#
# Usage:
#   stdres <- outlier$std.res
#   plot_stdres(stdres)
#   plot_stdres(stdres, labels = fullData$Short_reference)
plot_stdres <- function(stdres, threshold = 3, inspect = 1.96, labels = NULL) {
  df <- data.frame(id = seq_along(stdres), z = as.numeric(stdres))
  df$flag <- abs(df$z) >= threshold
  if (!is.null(labels)) df$label <- ifelse(df$flag, as.character(labels), NA_character_)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = id, y = z)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey55", linewidth = 0.3)

  if (!is.null(inspect)) {
    p <- p + ggplot2::geom_hline(
      yintercept = c(-inspect, inspect),
      color = "grey70", linetype = "dotted", linewidth = 0.4
    )
  }

  p <- p +
    ggplot2::geom_hline(
      yintercept = c(-threshold, threshold),
      color = "red", linetype = "dashed", linewidth = 0.7
    ) +
    ggplot2::geom_line(color = "grey60", linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(color = flag), size = 2.8) +
    ggplot2::scale_color_manual(
      values = c(`FALSE` = "#1f4e79", `TRUE` = "red"), guide = "none"
    ) +
    ggplot2::scale_x_continuous(breaks = df$id) +
    ggplot2::expand_limits(y = c(-threshold - 0.3, threshold + 0.3)) +
    ggplot2::labs(x = "Outcome ID", y = "Standardized residual")

  if (!is.null(labels)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = label),
      color = "red", size = 3, vjust = -0.8, na.rm = TRUE
    )
  }
  p
}
