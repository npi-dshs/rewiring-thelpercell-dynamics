# =============================================================================
# _chapter_helpers.R  -  r-switch, model ledger, result tables, interactive plots
# Sourced AFTER _common.R by every analysis chapter (Th1 / Th2 / Ratio).
# Each chapter defines a `cfg` list; the r-dependent sections are generated from
# it with the same model code as meta_Th1Th2_neu.R and knitted once per r file.
# =============================================================================

suppressMessages(library(plotly))

R_VALUES <- c("030", "040", "050", "060", "070")
R_LABEL  <- function(r) sprintf("r = .%s", substr(r, 2, 3))
R_MAIN   <- "050"

## ---- r-variant blocks -----------------------------------------------------
# One environment per r value; every r-dependent section is knitted once per
# environment from the SAME chunk text, so code and R output stay identical to
# the original chapter. Only the visible block is switched in the browser.
r_envs <- setNames(lapply(R_VALUES, function(r) {
  e <- new.env(parent = globalenv()); e$r_file <- r; e
}), R_VALUES)

r_block <- function(template) {
  if (length(template) == 1 && file.exists(template))
    template <- paste(readLines(template, warn = FALSE), collapse = "\n")
  out <- vapply(R_VALUES, function(r) {
    txt <- gsub("{r}", r, template, fixed = TRUE)
    md  <- knitr::knit_child(text = txt, envir = r_envs[[r]], quiet = TRUE)
    sprintf('\n\n::: {.r-variant data-r="%s"%s}\n%s\n:::\n\n',
            r, if (r == R_MAIN) " .active" else "", md)
  }, character(1))
  cat(out, sep = "\n")
}

# Copy the main-r environment into the global env (moderators: r = .50 only).
use_main_r <- function() {
  e <- r_envs[[R_MAIN]]
  for (n in ls(e)) assign(n, get(n, envir = e), envir = globalenv())
  invisible(NULL)
}

## ---- Model badge ----------------------------------------------------------
badge <- function(model, label, k, r = NULL) {
  txt <- paste0("Based on: ", model, " - ", label, " - k = ", k,
                if (!is.null(r)) paste0(" - ", R_LABEL(r)) else "")
  cls <- if (is.null(r)) "model-badge" else "model-badge fixed-r"
  sprintf('<span class="%s" data-k="%s">%s</span>', cls, k, txt)
}

## ---- Plotly theme bits ----------------------------------------------------
.pl_font   <- list(family = "Arial, Helvetica, Roboto, sans-serif", size = 13, color = "#1E1E1E")
.pl_config <- function(p, name) {
  plotly::config(
    p, displaylogo = FALSE,
    modeBarButtons = list(list("zoom2d", "pan2d", "resetScale2d", "toImage")),
    toImageButtonOptions = list(format = "svg", filename = name, width = 1000, height = 800)
  )
}
.axis <- function(title, ...) list(
  title = list(text = title, font = .pl_font), showline = TRUE, linecolor = "#000",
  mirror = TRUE, ticks = "outside", zeroline = FALSE, showgrid = FALSE, ...
)

## ---- Interactive contour-enhanced funnel (style of funnel_ce) ------------
funnel_ce_plotly <- function(m, ref = m$TE.random, name = "funnel") {
  d <- data.frame(Study = m$studlab, g = m$TE, SE = m$seTE)
  ymax <- max(d$SE) * 1.0
  xr   <- range(c(d$g, ref + c(-1, 1) * qnorm(.995) * ymax))
  xr   <- xr + c(-1, 1) * 0.02 * diff(xr)
  z <- qnorm(c(.95, .975, .995))
  band <- function(z_in, z_out, s) {
    list(x = c(ref, ref + s * z_out * ymax, ref + s * z_in * ymax), y = c(0, ymax, ymax))
  }
  p <- plotly::plot_ly(height = 620)
  # outside 99 % (light), 95-99 (medium), 90-95 (dark), inside (white)
  p <- plotly::add_polygons(p, x = c(xr[1], xr[2], xr[2], xr[1]), y = c(0, 0, ymax, ymax),
                            fillcolor = "#E0E0E0", line = list(width = 0),
                            hoverinfo = "skip", name = "p < 0.01")
  shades <- list(c(z[2], z[3], "#B2B2B2", "0.05 > p > 0.01"),
                 c(z[1], z[2], "#878787", "0.1 > p > 0.05"))
  for (s in shades) for (side in c(-1, 1)) {
    b <- band(as.numeric(s[1]), as.numeric(s[2]), side)
    p <- plotly::add_polygons(p, x = b$x, y = b$y, fillcolor = s[3],
                              line = list(width = 0), hoverinfo = "skip",
                              name = s[4], legendgroup = s[4], showlegend = side == 1)
  }
  p <- plotly::add_polygons(p, x = c(ref, ref - z[1] * ymax, ref + z[1] * ymax),
                            y = c(0, ymax, ymax), fillcolor = "#FFFFFF",
                            line = list(width = 0), hoverinfo = "skip", showlegend = FALSE)
  # 95 % pseudo-CI (dotted) and pooled estimate
  p <- plotly::add_lines(p, x = c(ref - z[2] * ymax, ref, ref + z[2] * ymax), y = c(ymax, 0, ymax),
                         line = list(color = "#000", dash = "dot", width = 1),
                         hoverinfo = "skip", showlegend = FALSE)
  p <- plotly::add_lines(p, x = c(ref, ref), y = c(0, ymax),
                         line = list(color = "#000", dash = "dot", width = 1),
                         hoverinfo = "skip", showlegend = FALSE)
  p <- plotly::add_markers(
    p, data = d, x = ~g, y = ~SE, marker = list(color = "#000", size = 9),
    text = ~sprintf("<b>%s</b><br>g = %.2f<br>SE = %.3f", Study, g, SE),
    hoverinfo = "text", showlegend = FALSE
  )
  p <- plotly::layout(
    p, font = .pl_font, plot_bgcolor = "#fff", paper_bgcolor = "#fff",
    xaxis = .axis("Hedges' g", range = xr),
    yaxis = .axis("Standard error", autorange = "reversed"),
    legend = list(x = 0.99, y = 0.99, xanchor = "right", bgcolor = "#fff",
                  bordercolor = "#000", borderwidth = 0.6,
                  traceorder = "reversed"),
    hoverlabel = list(bgcolor = "#fff", font = .pl_font),
    margin = list(l = 70, r = 20, t = 20, b = 60)
  )
  .pl_config(p, name)
}

## ---- Interactive standardised residuals (style of plot_stdres) -----------
plot_stdres_plotly <- function(stdres, labels, threshold = 3, inspect = 1.96,
                               name = "stdres") {
  d <- data.frame(id = seq_along(stdres), z = as.numeric(stdres), Study = labels)
  d$col <- ifelse(abs(d$z) >= threshold, "red", "#1f4e79")
  lim <- max(abs(d$z), threshold + 0.3) * 1.05
  hl <- function(y, col, dash, w) list(type = "line", xref = "paper", x0 = 0, x1 = 1,
                                       y0 = y, y1 = y, line = list(color = col, dash = dash, width = w))
  p <- plotly::plot_ly(d, x = ~id, y = ~z, height = 460)
  p <- plotly::add_lines(p, line = list(color = "grey", width = 1.2), hoverinfo = "skip")
  p <- plotly::add_markers(
    p, marker = list(color = ~col, size = 11),
    text = ~sprintf("<b>%s</b><br>Outcome ID %d<br>z = %.2f", Study, id, z),
    hoverinfo = "text"
  )
  p <- plotly::layout(
    p, font = .pl_font, showlegend = FALSE, plot_bgcolor = "#fff", paper_bgcolor = "#fff",
    shapes = list(hl(0, "#8c8c8c", "solid", 0.6),
                  hl(inspect, "#b3b3b3", "dot", 1), hl(-inspect, "#b3b3b3", "dot", 1),
                  hl(threshold, "red", "dash", 1.6), hl(-threshold, "red", "dash", 1.6)),
    xaxis = list(title = "Outcome ID", tickmode = "array", tickvals = d$id,
                 showline = TRUE, linecolor = "#000", ticks = "outside", zeroline = FALSE,
                 gridcolor = "#ebebeb"),
    yaxis = list(title = "Standardized residual", range = c(-lim, lim),
                 showline = TRUE, linecolor = "#000", ticks = "outside", zeroline = FALSE,
                 gridcolor = "#ebebeb"),
    hoverlabel = list(bgcolor = "#fff", font = .pl_font),
    margin = list(l = 70, r = 20, t = 10, b = 60)
  )
  .pl_config(p, name)
}

## ---- Interactive forest plot (test) ---------------------------------------
forest_plotly <- function(m, xlim = c(-2, 5), left = "Lower at post",
                          right = "Higher at post", name = "forest") {
  w <- m$w.random / sum(m$w.random)
  d <- data.frame(Study = m$studlab, g = m$TE, lo = m$lower, up = m$upper, w = w)
  d <- d[order(-d$g), ]
  n <- nrow(d); d$y <- rev(seq_len(n)) + 1
  est <- m$TE.random; elo <- m$lower.random; eup <- m$upper.random
  lab <- sprintf("%.2f [%.2f; %.2f]", d$g, d$lo, d$up)

  p <- plotly::plot_ly(height = 120 + 34 * (n + 3))
  for (i in seq_len(n)) {
    p <- plotly::add_segments(p, x = max(d$lo[i], xlim[1]), xend = min(d$up[i], xlim[2]),
                              y = d$y[i], yend = d$y[i],
                              line = list(color = "#1f1f1f", width = 1.2),
                              hoverinfo = "skip", showlegend = FALSE)
  }
  p <- plotly::add_markers(
    p, data = d, x = ~g, y = ~y, showlegend = FALSE,
    marker = list(symbol = "square", color = "#1f4e79", size = 7 + 12 * d$w / max(d$w),
                  line = list(color = "#1f1f1f", width = 1)),
    text = sprintf("<b>%s</b><br>g = %s<br>Weight = %.1f%%", d$Study, lab, 100 * d$w),
    hoverinfo = "text"
  )
  p <- plotly::add_polygons(
    p, x = c(elo, est, eup, est), y = c(0, 0.3, 0, -0.3),
    fillcolor = "#7fb3d5", line = list(color = "#1f1f1f", width = 1), showlegend = FALSE,
    text = sprintf("<b>Random effects</b><br>g = %.2f [%.2f; %.2f]<br>I² = %.1f%%",
                   est, elo, eup, 100 * m$I2),
    hoverinfo = "text", hoveron = "fills"
  )
  ann <- c(
    lapply(seq_len(n), function(i) list(x = 1.0, xref = "paper", xanchor = "left", y = d$y[i],
      text = sprintf("%s   %4.1f%%", lab[i], 100 * d$w[i]), showarrow = FALSE,
      font = list(family = "Arial, Helvetica, Roboto, sans-serif", size = 12))),
    list(list(x = 1.0, xref = "paper", xanchor = "left", y = 0,
              text = sprintf("<b>%.2f [%.2f; %.2f]   100.0%%</b>", est, elo, eup),
              showarrow = FALSE, font = list(size = 12)),
         list(x = 1.0, xref = "paper", xanchor = "left", y = n + 2.3,
              text = "<b>g [95% CI]          Weight</b>", showarrow = FALSE, font = list(size = 12)),
         list(x = -0.15, y = 0, yref = "paper", yanchor = "top", yshift = -26, xanchor = "right", text = left,
              showarrow = FALSE, font = list(size = 11)),
         list(x = 0.15, y = 0, yref = "paper", yanchor = "top", yshift = -26, xanchor = "left", text = right,
              showarrow = FALSE, font = list(size = 11))),
    lapply(seq_len(n), function(i) .lab(d$y[i], d$Study[i])),
    list(.lab(n + 2.3, "Study", bold = TRUE), .lab(0, "Total (95% CI)", bold = TRUE),
         .lab(-1.2, sprintf("Heterogeneity: Tau² = %.2f; Chi² = %.2f, df = %d (P = %.4f); I² = %.1f%%",
                            m$tau2, m$Q, m$df.Q, m$pval.Q, 100 * m$I2), size = 10.5))
  )
  p <- plotly::layout(
    p, font = .pl_font, plot_bgcolor = "#fff", paper_bgcolor = "#fff",
    shapes = list(
      list(type = "line", x0 = 0, x1 = 0, y0 = -0.6, y1 = n + 1.6, line = list(color = "grey", width = 1)),
      list(type = "line", x0 = est, x1 = est, y0 = -0.6, y1 = n + 1.6,
           line = list(color = "#000", width = 1, dash = "dot"))
    ),
    annotations = ann,
    xaxis = list(title = list(text = "Hedges' g", standoff = 28), range = xlim, showline = TRUE, linecolor = "#000",
                 ticks = "outside", zeroline = FALSE, showgrid = FALSE, dtick = 1),
    yaxis = list(title = "", showticklabels = FALSE, range = c(-2.8, n + 2.8),
                 showgrid = FALSE, zeroline = FALSE),
    hoverlabel = list(bgcolor = "#fff", font = .pl_font),
    margin = list(l = .LEFT_COL, r = 175, t = 10, b = 80)
  )
  .pl_config(p, name)
}

## ---- Left-aligned label column for plotly forests -------------------------
.LEFT_COL <- 245
.lab <- function(y, text, bold = FALSE, color = "#1E1E1E", size = 12) list(
  x = 0, xref = "paper", xanchor = "left", xshift = -.LEFT_COL + 4, y = y, showarrow = FALSE,
  text = if (bold) paste0("<b>", text, "</b>") else text, font = list(size = size, color = color))

## ---- Format helpers -------------------------------------------------------
fmt_ci <- function(est, lo, up) sprintf("%.2f [%.2f; %.2f]", est, lo, up)
fmt_p  <- function(p) ifelse(p < .001, "< .001", sub("^0", "", sprintf("%.3f", p)))

## ---- Format-aware wrappers (HTML -> plotly, PDF -> original static) -------
.funnel_ce_static  <- funnel_ce
.plot_stdres_static <- plot_stdres
funnel_ce <- function(m, ..., name = "funnel") {
  if (knitr::is_html_output()) funnel_ce_plotly(m, name = name) else .funnel_ce_static(m, ...)
}
plot_stdres <- function(stdres, ..., labels = NULL, name = "stdres") {
  if (knitr::is_html_output() && !is.null(labels)) {
    plot_stdres_plotly(stdres, labels = labels, name = name)
  } else .plot_stdres_static(stdres, ..., labels = labels)
}

## ---- Subgroup forest: tight auto-sized static version ---------------------
.forest_sub_static <- forest_sub
forest_sub <- function(m, xlim = c(-2, 5),
                       left = "Lower at post", right = "Higher at post") {
  if (is.null(knitr::opts_current$get("label"))) return(.forest_sub_static(m, xlim, left, right))
  rel <- knitr::fig_path(".png", number = 1)
  fn  <- file.path(knitr::opts_knit$get("output.dir"), rel)
  dir.create(dirname(fn), recursive = TRUE, showWarnings = FALSE)
  meta::forest(
    m, sortvar = -TE, layout = "RevMan5",
    common = FALSE, random = TRUE, prediction = FALSE, subgroup = TRUE,
    label.left = left, label.right = right, xlab = "Hedges' g",
    xlim = xlim, fontsize = 9, fontfamily = "Helvetica",
    method.tau = "HE", weight.study = "same",
    print.tau2 = TRUE, digits.tau2 = 2, digits.I2 = 1,
    ff.random = "bold", print.subgroup.labels = TRUE, subgroup.name = NULL,
    test.subgroup.random = TRUE, test.effect.subgroup.random = TRUE,
    col.study = "black", col.square = "#1f4e79", col.square.lines = "#1f1f1f",
    col.diamond = "#7fb3d5", col.diamond.random = "#7fb3d5",
    col.diamond.lines = "#1f1f1f", col.label = "black", col.lines = "grey40",
    col.subgroup = "#1f4e79", filename = fn
  )
  knitr::include_graphics(rel, rel_path = FALSE, error = FALSE)
}

## ---- Interactive subgroup forest ------------------------------------------
forest_sub_plotly <- function(m, xlim = c(-2, 5), left = "Lower at post",
                              right = "Higher at post") {
  w <- m$w.random / sum(m$w.random)
  d <- data.frame(Study = m$studlab, g = m$TE, lo = m$lower, up = m$upper,
                  w = w, grp = as.character(m$subgroup))
  lv <- m$subgroup.levels
  rows <- list(); y <- 0
  # build bottom-up: overall at y = 0, then groups above
  add <- function(...) rows[[length(rows) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)
  y <- 1.6
  for (gname in rev(lv)) {
    s <- d[d$grp == gname, ]; s <- s[order(s$g), ]
    add(type = "diamond", label = "Total (95% CI)", g = m$TE.random.w[[gname]],
        lo = m$lower.random.w[[gname]], up = m$upper.random.w[[gname]],
        w = sum(s$w), y = y, grp = gname)
    y <- y + 1
    for (i in seq_len(nrow(s))) {
      add(type = "study", label = s$Study[i], g = s$g[i], lo = s$lo[i], up = s$up[i],
          w = s$w[i], y = y, grp = gname); y <- y + 1
    }
    add(type = "header", label = gname, g = NA, lo = NA, up = NA, w = NA, y = y, grp = gname)
    y <- y + 1.6
  }
  r <- do.call(rbind, rows)
  st <- r[r$type == "study", ]; dm <- r[r$type == "diamond", ]
  est <- m$TE.random; top <- max(r$y) + 1
  fmt <- function(g, lo, up) sprintf("%.2f [%.2f; %.2f]", g, lo, up)

  p <- plotly::plot_ly(height = 140 + 30 * top)
  for (i in seq_len(nrow(st))) {
    p <- plotly::add_segments(p, x = max(st$lo[i], xlim[1]), xend = min(st$up[i], xlim[2]),
                              y = st$y[i], yend = st$y[i], line = list(color = "#1f1f1f", width = 1.2),
                              hoverinfo = "skip", showlegend = FALSE)
  }
  p <- plotly::add_markers(
    p, x = st$g, y = st$y, showlegend = FALSE,
    marker = list(symbol = "square", color = "#1f4e79", size = 7 + 12 * st$w / max(st$w),
                  line = list(color = "#1f1f1f", width = 1)),
    text = sprintf("<b>%s</b><br>%s<br>g = %s<br>Weight = %.1f%%", st$label, st$grp,
                   fmt(st$g, st$lo, st$up), 100 * st$w),
    hoverinfo = "text"
  )
  dia <- function(p, g, lo, up, yy, txt, h = 0.32) plotly::add_polygons(
    p, x = c(lo, g, up, g), y = c(yy, yy + h, yy, yy - h),
    fillcolor = "#7fb3d5", line = list(color = "#1f1f1f", width = 1), showlegend = FALSE,
    text = txt, hoverinfo = "text", hoveron = "fills")
  # subgroup estimates: CI line + diamond marker (arrow when CI exceeds the axis)
  for (i in seq_len(nrow(dm))) {
    lo <- max(dm$lo[i], xlim[1]); up <- min(dm$up[i], xlim[2])
    p <- plotly::add_segments(p, x = lo, xend = up, y = dm$y[i], yend = dm$y[i],
                              line = list(color = "#1f4e79", width = 2.2),
                              hoverinfo = "skip", showlegend = FALSE)
    if (dm$lo[i] < xlim[1]) p <- plotly::add_markers(p, x = xlim[1], y = dm$y[i], hoverinfo = "skip",
      showlegend = FALSE, marker = list(symbol = "triangle-left", size = 9, color = "#1f4e79"))
    if (dm$up[i] > xlim[2]) p <- plotly::add_markers(p, x = xlim[2], y = dm$y[i], hoverinfo = "skip",
      showlegend = FALSE, marker = list(symbol = "triangle-right", size = 9, color = "#1f4e79"))
  }
  p <- plotly::add_markers(
    p, x = dm$g, y = dm$y, showlegend = FALSE,
    marker = list(symbol = "diamond-wide", size = 18, color = "#7fb3d5",
                  line = list(color = "#1f1f1f", width = 1)),
    text = vapply(seq_len(nrow(dm)), function(i) { gname <- dm$grp[i]; sprintf(
      "<b>%s</b><br>g = %s<br>k = %d<br>Tau² = %.2f; I² = %.1f%%<br>p = %.4f",
      gname, fmt(dm$g[i], dm$lo[i], dm$up[i]), m$k.w[[gname]], m$tau2.w[[gname]],
      100 * m$I2.w[[gname]], m$pval.random.w[[gname]]) }, character(1)),
    hoverinfo = "text"
  )
  p <- dia(p, est, m$lower.random, m$upper.random, 0, sprintf(
    "<b>Overall</b><br>g = %s<br>Tau² = %.2f; I² = %.1f%%<br>Subgroup differences: Chi² = %.2f, df = %d, p = %.4f",
    fmt(est, m$lower.random, m$upper.random), m$tau2, 100 * m$I2,
    m$Q.b.random, m$df.Q.b.random, m$pval.Q.b.random))

  bold <- r$type != "study"
  ticktext <- ifelse(r$type == "header", sprintf("<b><span style='color:#1f4e79'>%s</span></b>", r$label),
               ifelse(r$type == "diamond", sprintf("<b><span style='color:#1f4e79'>%s</span></b>", r$label), r$label))
  ann_right <- lapply(which(r$type != "header"), function(i) list(
    x = 1.0, xref = "paper", xanchor = "left", y = r$y[i], showarrow = FALSE,
    text = sprintf(if (bold[i]) "<b>%s   %5.1f%%</b>" else "%s   %5.1f%%",
                   fmt(r$g[i], r$lo[i], r$up[i]), 100 * r$w[i]),
    font = list(size = 12, color = if (bold[i]) "#1f4e79" else "#1E1E1E")))
  ann <- c(ann_right, list(
    list(x = 1.0, xref = "paper", xanchor = "left", y = 0, showarrow = FALSE,
         text = sprintf("<b>%s   100.0%%</b>", fmt(est, m$lower.random, m$upper.random)),
         font = list(size = 12)),
    list(x = 1.0, xref = "paper", xanchor = "left", y = top, showarrow = FALSE,
         text = "<b>g [95% CI]          Weight</b>", font = list(size = 12)),
    list(x = -0.15, y = 0, yref = "paper", yanchor = "top", yshift = -26, xanchor = "right", text = left,
         showarrow = FALSE, font = list(size = 11)),
    list(x = 0.15, y = 0, yref = "paper", yanchor = "top", yshift = -26, xanchor = "left", text = right,
         showarrow = FALSE, font = list(size = 11)),
    .lab(top, "Study or Subgroup", bold = TRUE),
    .lab(0, "Total (95% CI)", bold = TRUE),
    .lab(-1.0, sprintf("Heterogeneity: Tau² = %.2f; Chi² = %.2f, df = %d (P = %.4f); I² = %.1f%%",
                       m$tau2, m$Q, m$df.Q, m$pval.Q, 100 * m$I2), size = 10.5),
    .lab(-1.8, sprintf("Test for subgroup differences: Chi² = %.2f, df = %d (P = %.4f)",
                       m$Q.b.random, m$df.Q.b.random, m$pval.Q.b.random), size = 10.5)
  ),
  lapply(seq_len(nrow(r)), function(i) {
    if (r$type[i] == "header") .lab(r$y[i], r$label[i], bold = TRUE, color = "#1f4e79")
    else if (r$type[i] == "diamond") .lab(r$y[i], r$label[i], bold = TRUE, color = "#1f4e79")
    else .lab(r$y[i], paste0("&nbsp;&nbsp;", r$label[i]))
  }))
  p <- plotly::layout(
    p, font = .pl_font, plot_bgcolor = "#fff", paper_bgcolor = "#fff",
    shapes = list(
      list(type = "line", x0 = 0, x1 = 0, y0 = -0.6, y1 = top - 0.5, line = list(color = "grey", width = 1)),
      list(type = "line", x0 = est, x1 = est, y0 = -0.6, y1 = top - 0.5,
           line = list(color = "#000", width = 1, dash = "dot"))
    ),
    annotations = ann,
    xaxis = list(title = list(text = "Hedges' g", standoff = 28), range = xlim, showline = TRUE, linecolor = "#000",
                 ticks = "outside", zeroline = FALSE, showgrid = FALSE, dtick = 1),
    yaxis = list(title = "", showticklabels = FALSE, range = c(-2.4, top + 0.8),
                 showgrid = FALSE, zeroline = FALSE),
    hoverlabel = list(bgcolor = "#fff", font = .pl_font),
    margin = list(l = .LEFT_COL, r = 175, t = 10, b = 80)
  )
  .pl_config(p, "forest-subgroup")
}

## ---- Clean tables ---------------------------------------------------------
.col_labels <- c(Short_reference = "Study", g = "g", std_res = "Std. residual",
                 dfbetas = "DFBETAS", dfbetas.intrcpt = "DFBETAS")
nice_table <- function(x, digits = 3, caption = NULL, col_names = NULL, align = NULL) {
  x <- do.call(data.frame, c(as.list(as.data.frame(x)), check.names = FALSE))
  if (is.null(col_names)) {
    col_names <- colnames(x)
    hit <- col_names %in% names(.col_labels)
    col_names[hit] <- .col_labels[col_names[hit]]
  }
  if (is.null(align)) align <- ifelse(vapply(x, is.numeric, logical(1)), "r", "l")
  tab <- knitr::kable(x, digits = digits, caption = caption, col.names = col_names,
                      align = align, booktabs = TRUE, linesep = "")
  if (knitr::is_html_output()) {
    tab <- kableExtra::kable_styling(tab, bootstrap_options = "condensed",
                                     full_width = FALSE, position = "left")
  }
  tab
}

## ---- Result tables + collapsible raw R output -----------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a
.f2 <- function(x, d = 2) ifelse(is.na(x), "-", formatC(x, format = "f", digits = d))
.ci <- function(lo, up, d = 2) sprintf("[%s; %s]", .f2(lo, d), .f2(up, d))
.pp <- function(p) ifelse(is.na(p), "-", ifelse(p < .001, "&lt; .001", sub("^0", "", sprintf("%.3f", p))))
# APA: p italic; significant p values bold
.pcell <- function(p) ifelse(!is.na(p) & p < .05, paste0("<b>", .pp(p), "</b>"), .pp(p))
.ptxt  <- function(p) {
  v <- .pp(p); s <- paste0("<i>p</i> ", if (grepl("&lt;", v)) v else paste("=", v))
  ifelse(!is.na(p) & p < .05, paste0("<b>", s, "</b>"), s)
}
.SIG <- "Bold: <i>p</i> &lt; .05."
.PH  <- "<i>p</i>"
.note <- function(tab, txt) {
  txt <- txt[nzchar(txt)]
  if (!length(txt)) return(tab)
  kableExtra::footnote(tab, general_title = "Note.", title_format = "italic",
                       general = paste(txt, collapse = " "), footnote_as_chunk = TRUE, escape = FALSE)
}
.html_table <- function(df, caption, align, cls = "table") {
  tab <- knitr::kable(df, format = "html", caption = caption, align = align, escape = FALSE,
                      row.names = FALSE, table.attr = sprintf('class="%s"', cls))
  kableExtra::kable_styling(tab, bootstrap_options = "condensed", full_width = FALSE, position = "left")
}

.kv_table <- function(rows, groups, caption, note = NULL, sig = FALSE) {
  df <- do.call(rbind, lapply(rows, function(r) data.frame(Parameter = r[1], Estimate = r[2],
                                                           `95% CI` = r[3], Test = r[4],
                                                           check.names = FALSE)))
  keep <- c(TRUE, vapply(df[-1], function(v) any(nzchar(v)), logical(1)))
  df <- df[, keep, drop = FALSE]
  tab <- .html_table(df, caption, c("l", "r", "l", "l")[keep], "table kv-table")
  start <- 1
  for (g in names(groups)) {
    tab <- kableExtra::pack_rows(tab, g, start, start + groups[[g]] - 1,
                                 label_row_css = "border-bottom: 0; color: #0053A0; font-weight: 600;")
    start <- start + groups[[g]]
  }
  .note(tab, c(note %||% "", if (isTRUE(sig)) .SIG else ""))
}

.tbl_meta <- function(m, caption) {
  pooled <- list(
    c("<i>k</i>", m$k, "", ""),
    c("Hedges' <i>g</i> (random effects)", .f2(m$TE.random), .ci(m$lower.random, m$upper.random),
      sprintf("<i>t</i>(%d) = %s, %s", m$k - 1, .f2(m$statistic.random), .ptxt(m$pval.random))),
    c("Prediction interval", "", .ci(m$lower.predict, m$upper.predict), "")
  )
  het <- list(
    c("&tau;&sup2;", .f2(m$tau2), .ci(m$lower.tau2, m$upper.tau2), ""),
    c("&tau;", .f2(m$tau), .ci(m$lower.tau, m$upper.tau), ""),
    c("<i>I</i>&sup2;", paste0(.f2(100 * m$I2, 1), "%"),
      sprintf("[%s%%; %s%%]", .f2(100 * m$lower.I2, 1), .f2(100 * m$upper.I2, 1)), ""),
    c("<i>H</i>", .f2(m$H), .ci(m$lower.H, m$upper.H), ""),
    c("<i>Q</i>", .f2(m$Q), "", sprintf("<i>df</i> = %d, %s", m$df.Q, .ptxt(m$pval.Q)))
  )
  note <- sprintf("Random effects, inverse variance; &tau;&sup2; estimator: %s; CI: %s.",
                  m$method.tau, m$method.random.ci)
  .kv_table(c(pooled, het), list(`Pooled effect` = 3, Heterogeneity = 5), caption, note,
            sig = any(c(m$pval.random, m$pval.Q) < .05, na.rm = TRUE))
}

.tbl_meta_sub <- function(m, caption) {
  lv <- m$subgroup.levels
  df <- data.frame(
    Subgroup = lv, k = as.integer(m$k.w[lv]),
    g = .f2(m$TE.random.w[lv]), `95% CI` = .ci(m$lower.random.w[lv], m$upper.random.w[lv]),
    p = .pcell(m$pval.random.w[lv]), `τ²` = .f2(m$tau2.w[lv]),
    `I²` = paste0(.f2(100 * m$I2.w[lv], 1), "%"), check.names = FALSE
  )
  names(df) <- c("Subgroup", "<i>k</i>", "<i>g</i>", "95% CI", .PH, "&tau;&sup2;", "<i>I</i>&sup2;")
  tab <- .html_table(df, caption, c("l", "r", "r", "l", "r", "r", "r"))
  qb <- sprintf("Test for subgroup differences (random effects): <i>Q</i> = %s, <i>df</i> = %d, %s.",
                .f2(m$Q.b.random), m$df.Q.b.random, .ptxt(m$pval.Q.b.random))
  .note(tab, c(qb, if (any(c(m$pval.random.w[lv], m$pval.Q.b.random) < .05, na.rm = TRUE)) .SIG else ""))
}

.tbl_rma <- function(x, caption) {
  if (isTRUE(x$int.only)) {
    rows <- list(
      c("<i>k</i>", x$k, "", ""),
      c("Estimate", .f2(x$b[1]), .ci(x$ci.lb, x$ci.ub),
        sprintf("%s = %s, %s", if (x$test == "t") sprintf("<i>t</i>(%d)", x$dfs) else "<i>z</i>",
                .f2(x$zval), .ptxt(x$pval))),
      c("&tau;&sup2;", .f2(x$tau2), "", ""), c("<i>I</i>&sup2;", paste0(.f2(x$I2, 1), "%"), "", ""),
      c("<i>H</i>&sup2;", .f2(x$H2), "", ""),
      c("<i>Q</i>", .f2(x$QE), "", sprintf("<i>df</i> = %d, %s", x$k - x$p, .ptxt(x$QEp)))
    )
    if (!is.null(x$k0)) rows <- c(list(c("Imputed studies (<i>k</i>0)", x$k0, "", paste("side:", x$side))), rows)
    grp <- if (!is.null(x$k0)) list(`Trim and fill` = 1, `Adjusted model` = 6) else list(`Pooled effect` = 2, Heterogeneity = 4)
    return(.kv_table(rows, grp, caption, sig = any(c(x$pval, x$QEp) < .05, na.rm = TRUE)))
  }
  df <- data.frame(Term = rownames(x$b), b = .f2(x$b[, 1], 3), SE = .f2(x$se, 3),
                   t = .f2(x$zval), p = .pcell(x$pval), ci = .ci(x$ci.lb, x$ci.ub))
  names(df) <- c("Term", "<i>b</i>", "<i>SE</i>", if (x$test == "t") "<i>t</i>" else "<i>z</i>", .PH, "95% CI")
  tab <- .html_table(df, caption, c("l", "r", "r", "r", "r", "l"))
  .note(tab, c(
    sprintf("<i>k</i> = %d; residual &tau;&sup2; = %s; residual <i>I</i>&sup2; = %s%%; <i>R</i>&sup2; = %s%%.",
            x$k, .f2(x$tau2), .f2(x$I2, 1), ifelse(is.null(x$R2) || is.na(x$R2), "-", .f2(x$R2, 1))),
    sprintf("Test of moderators: %s = %s, %s; residual heterogeneity: <i>Q</i>(%d) = %s, %s.",
            if (x$test == "t") sprintf("<i>F</i>(%d, %d)", x$QMdf[1], x$QMdf[2]) else sprintf("<i>QM</i>(%d)", x$QMdf[1]),
            .f2(x$QM), .ptxt(x$QMp), x$k - x$p, .f2(x$QE), .ptxt(x$QEp)),
    if (any(c(x$pval, x$QMp, x$QEp) < .05, na.rm = TRUE)) .SIG else ""))
}

.tbl_bias <- function(x, caption, msg) {
  if (is.null(x$statistic)) {
    kmin <- suppressWarnings(as.integer(sub(".*k\\.min=(\\d+).*", "\\1", msg %||% "")))
    rows <- list(c("<i>k</i>", x$k %||% "", "", ""), c("<i>k</i>.min", if (is.na(kmin)) "-" else kmin, "", ""),
                 c("Status", "Not computed", "", ""))
    return(.kv_table(rows, list(`Egger's regression test` = 3), caption, note = msg))
  }
  rows <- list(
    c("<i>k</i>", x$k, "", ""),
    c("Bias (intercept)", .f2(x$estimate[1]), "", sprintf("<i>SE</i> = %s", .f2(x$estimate[2]))),
    c("Test", "", "", sprintf("<i>t</i>(%d) = %s, %s", x$df, .f2(x$statistic), .ptxt(x$pval)))
  )
  .kv_table(rows, list(`Egger's regression test` = 3), caption, sig = isTRUE(x$pval < .05))
}

report <- function(x, caption = NULL, ...) {
  msg <- NULL
  raw <- withCallingHandlers(
    capture.output(print(x, ...)),
    warning = function(w) { msg <<- conditionMessage(w); invokeRestart("muffleWarning") }
  )
  if (!knitr::is_html_output()) { cat(raw, sep = "\n"); return(invisible(x)) }
  if (inherits(x, "metabias") && is.null(msg) && is.null(x$statistic)) msg <- attr(x, "warning")
  tab <- if (inherits(x, "metabias")) .tbl_bias(x, caption %||% "Egger's regression test", msg)
    else if (inherits(x, "meta") && !is.null(x$subgroup)) .tbl_meta_sub(x, caption %||% "Subgroup analysis")
    else if (inherits(x, "meta")) .tbl_meta(x, caption %||% "Random-effects meta-analysis")
    else if (inherits(x, "rma.uni.trimfill")) .tbl_rma(x, caption %||% "Trim and fill")
    else if (inherits(x, "rma")) .tbl_rma(x, caption %||% if (isTRUE(x$int.only)) "Random-effects model" else "Meta-regression")
  esc <- function(s) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", s)))
  if (!length(raw) || all(!nzchar(raw))) raw <- msg %||% ""
  knitr::asis_output(paste0(
    "\n\n```{=html}\n", as.character(tab),
    '\n<details class="raw-output"><summary>R output</summary><pre><code>',
    esc(paste(raw, collapse = "\n")), "</code></pre></details>\n```\n\n"))
}


## ---- Sensitivity overview plot (same type size as the text) --------------
sens_plotly <- function(d) {
  d$r <- factor(d$r, levels = rev(d$r))
  p <- plotly::plot_ly(d, height = 230)
  p <- plotly::add_markers(
    p, x = ~g, y = ~r, error_x = list(type = "data", symmetric = FALSE, array = d$up - d$g,
                                      arrayminus = d$g - d$lo, color = "#1f1f1f", thickness = 1, width = 4),
    marker = list(symbol = "square", size = 9, color = ifelse(d$main, "#1f4e79", "#7fb3d5"),
                  line = list(color = "#1f1f1f", width = 1)),
    text = sprintf("<b>%s</b><br>g = %.2f [%.2f; %.2f]", d$r, d$g, d$lo, d$up), hoverinfo = "text"
  )
  p <- plotly::layout(
    p, font = list(family = "Arial, Helvetica, Roboto, sans-serif", size = 12, color = "#1E1E1E"),
    plot_bgcolor = "rgba(0,0,0,0)", paper_bgcolor = "rgba(0,0,0,0)", showlegend = FALSE,
    shapes = list(list(type = "line", x0 = 0, x1 = 0, yref = "paper", y0 = 0, y1 = 1,
                       line = list(color = "grey", width = 1))),
    xaxis = list(title = "Hedges' g [95% CI]", zeroline = FALSE, showline = TRUE, linecolor = "#000",
                 ticks = "outside", gridcolor = "#ebebeb", rangemode = "tozero"),
    yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, ticks = ""),
    margin = list(l = 60, r = 10, t = 5, b = 45)
  )
  plotly::config(p, displayModeBar = FALSE)
}

## =============================================================================
## Chapter configuration + section generators
## =============================================================================
# cfg fields
#   id        chunk-label prefix, e.g. "th1-imm"
#   outcome   "Th1" | "Th2" | "Ratio";  win  "pre_post" | "pre_post+2h" | "pre_post17h"
#   M         model object name used by the original code (m1, m3, m5, m2, ...)
#   drop_kost TRUE -> M1 - Primary (Kostrzewa-Nowak et al. (2019) removed)
#   tabs      order of the primary tabs, from "M0-meta", "M0-metafor", "M1-meta", "M1-metafor"
#   m1_note   optional callout text shown in the M1 tab (verbatim from the R script)
#   refit     NULL or list(study = "...", M2 = "m2_1")  -> M2 - Outlier-adjusted
#   forest    list(w, h, args)  args = extra arguments of forest_main() as text
#   asym      "egger_tf" | "egger" | "none";  tf_data  data used for trim-and-fill
#   lancaster TRUE -> Lancaster et al. (2005) crossover sensitivity section

ES_DIR   <- "book/data/effect_size"
LANC_DIR <- "book/data/effect_size_lancaster"
es_file  <- function(cfg, r, dir = ES_DIR)
  sprintf("%s/%s/%s/%s_%s_r_%s.csv", dir, cfg$outcome, cfg$win, cfg$outcome, cfg$win, r)

.fit <- function(d) meta::metagen(g, seTE = SE_g, data = d, studlab = paste(Study),
                                  random = TRUE, method.tau = "HE", hakn = TRUE,
                                  prediction = TRUE, sm = "SMD", fixed = FALSE)

model_ids <- function(cfg) c("M0", if (isTRUE(cfg$drop_kost)) "M1", if (!is.null(cfg$refit)) "M2")
MODEL_LABEL <- c(M0 = "Full", M1 = "Primary", M2 = "Outlier-adjusted")
reported_id <- function(cfg) tail(model_ids(cfg), 1)
primary_id  <- function(cfg) if (isTRUE(cfg$drop_kost)) "M1" else "M0"
mlab <- function(id) paste(id, "-", MODEL_LABEL[[id]])

# data behind each model
model_data <- function(cfg, d0, id) {
  d <- d0
  if (id %in% c("M1", "M2") && isTRUE(cfg$drop_kost)) d <- drop_kostrzewa19(d)
  if (id == "M2") d <- dplyr::filter(d, !Study == cfg$refit$study)
  d
}

chapter_chain <- function(cfg) {
  lapply(setNames(R_VALUES, R_VALUES), function(r) {
    d0 <- load_es(es_file(cfg, r))
    lapply(setNames(model_ids(cfg), model_ids(cfg)), function(id) .fit(model_data(cfg, d0, id)))
  })
}

.row_of <- function(m) data.frame(
  k = m$k, g = .f2(m$TE.random), ci = .ci(m$lower.random, m$upper.random),
  p = .pcell(m$pval.random), tau2 = .f2(m$tau2), I2 = paste0(.f2(100 * m$I2, 1), "%"),
  pi = .ci(m$lower.predict, m$upper.predict)
)

## ---- Summary panel --------------------------------------------------------
sens_overview_table <- function(cfg, chain) {
  rep <- reported_id(cfg)
  rows <- vapply(R_VALUES, function(r) {
    x <- .row_of(chain[[r]][[rep]])
    sprintf('<tr data-r="%s"%s><td>%s%s</td><td>%d</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>',
            r, if (r == R_MAIN) ' class="is-main"' else "", R_LABEL(r), if (r == R_MAIN) " (main)" else "",
            x$k, x$g, x$ci, x$p, x$tau2, x$I2, x$pi)
  }, character(1))
  sig <- any(vapply(chain, function(ch) ch[[rep]]$pval.random < .05, logical(1)))
  knitr::asis_output(paste0(
    '\n```{=html}\n<div class="table-scroll"><table class="table table-sm sens-table">',
    '<caption>Sensitivity analysis - ', mlab(rep), '</caption><thead><tr>',
    '<th>Effect size file</th><th><i>k</i></th><th><i>g</i></th><th>95% CI</th><th><i>p</i></th>',
    '<th>&tau;&sup2;</th><th><i>I</i>&sup2;</th><th>95% PI</th></tr></thead><tbody>',
    paste(rows, collapse = ""), '</tbody>',
    if (sig) paste0('<tfoot><tr><td colspan="8"><i>Note.</i> ', .SIG, '</td></tr></tfoot>') else "",
    '</table></div>\n```\n'))
}

sens_overview_plot <- function(cfg, chain) {
  rep <- reported_id(cfg)
  sens_plotly(do.call(rbind, lapply(R_VALUES, function(r) {
    m <- chain[[r]][[rep]]
    data.frame(r = R_LABEL(r), g = m$TE.random, lo = m$lower.random, up = m$upper.random,
               main = r == R_MAIN)
  })))
}

ledger_table <- function(cfg, chain) {
  ids <- model_ids(cfg); ch <- chain[[R_MAIN]]
  excl <- c(M0 = "-", M1 = "Kostrzewa-Nowak et al. (2019)",
            M2 = paste0(if (isTRUE(cfg$drop_kost)) "Kostrzewa-Nowak et al. (2019), " else "", cfg$refit$study %||% ""))
  df <- do.call(rbind, lapply(ids, function(id) {
    x <- .row_of(ch[[id]])
    data.frame(Model = mlab(id), Excluded = excl[[id]], k = x$k, g = x$g, ci = x$ci, p = x$p,
               tau2 = x$tau2, I2 = x$I2, Status = if (id == reported_id(cfg)) "Reported" else "")
  }))
  names(df) <- c("Model", "Excluded", "<i>k</i>", "<i>g</i>", "95% CI", .PH, "&tau;&sup2;", "<i>I</i>&sup2;", "Status")
  tab <- .html_table(df, paste0("Model ledger - ", R_LABEL(R_MAIN)),
                     c("l", "l", "r", "r", "l", "r", "r", "r", "l"))
  tab <- .note(tab, if (any(vapply(ch, function(m) m$pval.random < .05, logical(1)))) .SIG else "")
  knitr::asis_output(paste0('\n```{=html}\n<div class="ledger-wrap">', as.character(tab), '</div>\n```\n'))
}

r_switch <- function(cfg) {
  pre <- sprintf("%s_%s_r_", cfg$outcome, cfg$win)
  opts <- paste(vapply(R_VALUES, function(r) sprintf('<option value="%s"%s>r = .%s%s</option>',
    r, if (r == R_MAIN) " selected" else "", substr(r, 2, 3), if (r == R_MAIN) " (main)" else ""),
    character(1)), collapse = "")
  knitr::asis_output(sprintf(paste0(
    '\n```{=html}\n<div class="r-switch" id="r-switch" data-prefix="%s">',
    '<label for="r-select">Effect size file</label>',
    '<select id="r-select" class="form-select form-select-sm">%s</select>',
    '<code class="r-file">%s050.csv</code>',
    '<span class="r-fixed-note">Moderators - r = .50</span></div>\n```\n'), pre, opts, pre))
}

## ---- Section templates ----------------------------------------------------
.chunk <- function(label, code, opts = character(0))
  paste0("```{r ", label, "}\n", paste0(opts, collapse = ""), code, "\n```\n")
.fig <- function(w = NULL, h = NULL)
  paste0(if (!is.null(w)) sprintf("#| fig-width: %s\n", w), if (!is.null(h)) sprintf("#| fig-height: %s\n", h),
         "#| out-width: 100%\n")
.badge_inline <- function(id, k_expr, r = NULL)
  sprintf('`r badge("%s", "%s", %s%s)`\n\n', id, MODEL_LABEL[[id]], k_expr,
          if (is.null(r)) "" else sprintf(', r = "%s"', r))

sec_data <- function(cfg)
  .chunk(paste0(cfg$id, "-data-{r}"),
         sprintf('fullData <- load_es("%s")', es_file(cfg, "{r}")))

sec_primary <- function(cfg) {
  M <- cfg$M
  fit_full <- sprintf(paste0(
    '%s <- meta::metagen(TE = g, seTE = SE_g, data = fullData, studlab = paste(Study),\n',
    '                    random = TRUE, method.tau = "HE", fixed = FALSE,\n',
    '                    hakn = TRUE, prediction = TRUE, sm = "SMD")\n'), M)
  fit_adapt <- sprintf(paste0(
    'fullData <- drop_kostrzewa19(fullData)\n',
    '%s <- meta::metagen(g, seTE = SE_g, data = fullData, studlab = paste(Study),\n',
    '                    random = TRUE, method.tau = "HE", hakn = TRUE,\n',
    '                    prediction = TRUE, sm = "SMD", fixed = FALSE)\n'), M)
  metafor <- function(id) sprintf(paste0(
    'report(metafor::rma(yi = g, sei = SE_g, data = fullData, slab = paste(Study),\n',
    '                    method = "HE", test = "knha"), caption = "%s - <i>metafor</i>")'), mlab(id))
  tab <- function(key, i) {
    id <- sub("-.*", "", key); pkg <- sub(".*-", "", key)
    head <- sprintf("### %s - {%s}\n\n", mlab(id), pkg)
    lab <- sprintf("%s-main-%d-{r}", cfg$id, i)
    if (pkg == "metafor")
      return(paste0(head, .badge_inline(id, "nrow(fullData)"), .chunk(lab, metafor(id))))
    if (id == "M0")
      return(paste0(head, .badge_inline(id, "nrow(fullData)"),
                    .chunk(lab, paste0(fit_full, sprintf('report(%s, caption = "%s - <i>meta</i>")', M, mlab(id))))))
    paste0(head, .badge_inline(id, "nrow(drop_kostrzewa19(fullData))"),
           if (!is.null(cfg$m1_note)) sprintf("::: {.callout-caution}\n%s\n:::\n\n", cfg$m1_note) else "",
           .chunk(lab, paste0(fit_adapt, sprintf('report(%s, caption = "%s - <i>meta</i>")', M, mlab(id)))))
  }
  paste0("::: panel-tabset\n", paste(mapply(tab, cfg$tabs, seq_along(cfg$tabs)), collapse = "\n"), ":::\n")
}

sec_outlier <- function(cfg) {
  M <- cfg$M; id <- primary_id(cfg)
  paste0(.badge_inline(id, paste0(M, "$k")),
    .chunk(paste0(cfg$id, "-outliers-{r}"), sprintf(
      'outlier <- metaoutliers(y = %s$TE, s2 = %s$seTE^2, model = "RE")\nstdres <- outlier$std.res', M, M)),
    "\n::: panel-tabset\n### Plot\n",
    .chunk(paste0(cfg$id, "-outliers-plot-{r}"), "plot_stdres(stdres, labels = fullData$Short_reference)", .fig(h = 6)),
    "\n### Table\n",
    .chunk(paste0(cfg$id, "-resid-{r}"), paste0(
      'fullData |> select(Short_reference, g) |> mutate(std_res = outlier$std.res) |>\n',
      '  nice_table(caption = "Standardised residuals")')),
    ":::\n")
}

sec_funnel <- function(cfg) {
  M <- cfg$M; id <- primary_id(cfg)
  paste0(.badge_inline(id, paste0(M, "$k")),
    .chunk(paste0(cfg$id, "-funnel-{r}"), sprintf("funnel_ce(%s)", M), .fig(10, 8)), "\n",
    .chunk(paste0(cfg$id, "-dfbetas-{r}"), paste0(
      'overall_metafor <- rma.mv(g, VE_g, tdist = TRUE, data = fullData)\n',
      'df_betas <- dfbetas(overall_metafor)')),
    "\n::: panel-tabset\n### DFBETAS - Plot\n",
    .chunk(paste0(cfg$id, "-dfbetas-plot-{r}"),
           "plot_dfbetas(df_betas[[1]], labels = fullData$Short_reference)", .fig(h = 5)),
    "\n### DFBETAS - Table\n",
    .chunk(paste0(cfg$id, "-dfbetas-tab-{r}"), paste0(
      'fullData |> select(Short_reference, g) |> mutate(dfbetas = df_betas) |>\n',
      '  nice_table(caption = "DFBETAS influence diagnostics")')),
    ":::\n")
}

sec_refit <- function(cfg) {
  M <- cfg$M; M2 <- cfg$refit$M2
  paste0(.badge_inline("M2", paste0(M, "$k - 1")),
    .chunk(paste0(cfg$id, "-refit-{r}"), sprintf(paste0(
      '# Withdraw the outlier/influential case and re-fit\n',
      'fullData_bis <- fullData |> dplyr::filter(!Study == "%s")\n',
      '%s <- meta::metagen(g, seTE = SE_g, data = fullData_bis, studlab = paste(Study),\n',
      '                      random = TRUE, method.tau = "HE", hakn = TRUE,\n',
      '                      prediction = TRUE, sm = "SMD", fixed = FALSE)\n',
      'report(%s, caption = "M2 - Outlier-adjusted - <i>meta</i>")'), cfg$refit$study, M2, M2)), "\n",
    .chunk(paste0(cfg$id, "-funnel2-{r}"), sprintf("funnel_ce(%s)", M2), .fig(10, 8)))
}

rep_obj <- function(cfg) if (!is.null(cfg$refit)) cfg$refit$M2 else cfg$M

sec_forest <- function(cfg) {
  R <- rep_obj(cfg); a <- cfg$forest$args
  args <- if (nzchar(a %||% "")) paste0(", ", a) else ""
  paste0(.badge_inline(reported_id(cfg), paste0(R, "$k")),
    "::: panel-tabset\n### Interactive\n::: forest-scroll\n",
    .chunk(paste0(cfg$id, "-forest-plotly-{r}"), sprintf("forest_plotly(%s%s)", R, args), .fig(h = 6)),
    ":::\n\n### RevMan5\n",
    .chunk(paste0(cfg$id, "-forest-{r}"), sprintf("forest_main(%s%s)", R, args),
           .fig(cfg$forest$w %||% 8, cfg$forest$h %||% 7)),
    ":::\n")
}

sec_asym <- function(cfg) {
  R <- rep_obj(cfg)
  egger <- .chunk(paste0(cfg$id, "-egger-{r}"),
                  sprintf('report(metabias(%s, method = "linreg"), caption = "Egger\'s regression test")', R))
  if (cfg$asym == "egger") return(paste0(.badge_inline(reported_id(cfg), paste0(R, "$k")), egger))
  paste0(.badge_inline(reported_id(cfg), paste0(R, "$k")),
    "::: {.panel-tabset}\n\n### Egger´s regression test\n", egger,
    "\n### Trim and Fill\n",
    .chunk(paste0(cfg$id, "-trim_fill-{r}"), sprintf(paste0(
      '# Trim-and-fill\n',
      'res <- metafor::rma(yi = g, vi = V_g, data = %s, method = "HE")\n',
      'taf <- metafor::trimfill(res)\n',
      'report(taf, caption = "Trim and fill")'), cfg$tf_data %||% "fullData")),
    "\n::: {.callout-note}\n# No missing studies\n:::\n\n[r = .50]{.callout-r}\n\n:::\n")
}

## ---- Lancaster et al. (2005) crossover sensitivity ------------------------
sec_lancaster <- function(cfg)
  .chunk(paste0(cfg$id, "-sens-lancaster-{r}"), 'sens_lancaster(cfg, "{r}")')

sens_lancaster <- function(cfg, r) {
  id <- reported_id(cfg)
  dm <- model_data(cfg, load_es(es_file(cfg, r)), id)
  dl <- model_data(cfg, load_es(es_file(cfg, r, LANC_DIR)), id)
  fits <- list(
    `Main analysis` = .fit(dm),
    `Lancaster et al. (2005) - <i>n</i> / 3` = .fit(dl),
    `Pre/post designs only` = .fit(dplyr::filter(dm, !grepl("^Lancaster", Study)))
  )
  df <- do.call(rbind, lapply(names(fits), function(n) cbind(Analysis = n, .row_of(fits[[n]])[, 1:6])))
  names(df) <- c("Analysis", "<i>k</i>", "<i>g</i>", "95% CI", .PH, "&tau;&sup2;", "<i>I</i>&sup2;")
  tab <- .html_table(df, paste0("Sensitivity analyses - ", mlab(id)), c("l", "r", "r", "l", "r", "r", "r"))
  tab <- .note(tab, if (any(vapply(fits, function(m) m$pval.random < .05, logical(1)))) .SIG else "")
  knitr::asis_output(paste0("\n```{=html}\n", as.character(tab), "\n```\n"))
}

## ---- Shortcuts used in the chapter files ----------------------------------
sec <- function(fun) r_block(fun(cfg))

## ---- Interactive DFBETAS (style of plot_stdres) ---------------------------
plot_dfbetas_plotly <- function(dfb, labels, threshold = 1, name = "dfbetas") {
  d <- data.frame(id = seq_along(dfb), v = as.numeric(dfb), Study = labels)
  d$col <- ifelse(abs(d$v) > threshold, "red", "#1f4e79")
  lim <- max(abs(d$v), threshold + 0.2) * 1.08
  hl <- function(y, col, dash, w) list(type = "line", xref = "paper", x0 = 0, x1 = 1,
                                       y0 = y, y1 = y, line = list(color = col, dash = dash, width = w))
  p <- plotly::plot_ly(d, x = ~id, y = ~v, height = 420)
  p <- plotly::add_lines(p, line = list(color = "grey", width = 1.2), hoverinfo = "skip")
  p <- plotly::add_markers(
    p, marker = list(color = ~col, size = 11),
    text = ~sprintf("<b>%s</b><br>Outcome ID %d<br>DFBETAS = %.2f", Study, id, v), hoverinfo = "text"
  )
  p <- plotly::layout(
    p, font = .pl_font, showlegend = FALSE, plot_bgcolor = "#fff", paper_bgcolor = "#fff",
    shapes = list(hl(0, "#8c8c8c", "solid", 0.6), hl(threshold, "red", "dash", 1.6), hl(-threshold, "red", "dash", 1.6)),
    xaxis = list(title = "Outcome ID", tickmode = "array", tickvals = d$id, showline = TRUE,
                 linecolor = "#000", ticks = "outside", zeroline = FALSE, gridcolor = "#ebebeb"),
    yaxis = list(title = "DFBETAS", range = c(-lim, lim), showline = TRUE, linecolor = "#000",
                 ticks = "outside", zeroline = FALSE, gridcolor = "#ebebeb"),
    hoverlabel = list(bgcolor = "#fff", font = .pl_font),
    margin = list(l = 70, r = 20, t = 10, b = 60)
  )
  .pl_config(p, name)
}

plot_dfbetas <- function(dfb, labels, threshold = 1) {
  if (knitr::is_html_output()) return(plot_dfbetas_plotly(dfb, labels, threshold))
  d <- data.frame(id = seq_along(dfb), v = as.numeric(dfb))
  d$flag <- abs(d$v) > threshold
  ggplot2::ggplot(d, ggplot2::aes(id, v)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey55", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = c(-threshold, threshold), color = "red", linetype = "dashed", linewidth = 0.7) +
    ggplot2::geom_line(color = "grey60", linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(color = flag), size = 2.8) +
    ggplot2::scale_color_manual(values = c(`FALSE` = "#1f4e79", `TRUE` = "red"), guide = "none") +
    ggplot2::scale_x_continuous(breaks = d$id) +
    ggplot2::labs(x = "Outcome ID", y = "DFBETAS")
}

## ---- Page resources (shown under the TOC by r-switch.html) ----------------
REPO_URL <- "https://github.com/npi-dshs/rewiring-thelpercell-dynamics"
DL_DIR   <- "downloads"   # relative to book/, copied to _book/downloads (resources)

# links: list of c(label, href, icon[, pattern]); href "#" = placeholder,
# pattern = href with {r} that the r-switch keeps in sync with the selected file
page_resources <- function(links) {
  a <- vapply(links, function(l) {
    ext <- grepl("^https?://", l[2])
    sprintf('<a href="%s"%s data-icon="%s"%s>%s</a>', l[2],
            if (l[2] == "#") ' class="is-ph"' else if (ext) ' target="_blank" rel="noopener"'
            else sprintf(' download="%s"', basename(l[2])),
            l[3], if (length(l) > 3) sprintf(' data-pattern="%s"', l[4]) else "", l[1])
  }, character(1))
  knitr::asis_output(paste0('\n```{=html}\n<div id="page-resources" hidden>', paste(a, collapse = ""), '</div>\n```\n'))
}

# data ZIP (all r, incl. Lancaster sensitivity files) + results XLSX for one chapter
chapter_downloads <- function(cfg, chain) {
  dir.create(DL_DIR, showWarnings = FALSE)
  stem <- sprintf("%s_%s", cfg$outcome, cfg$win)
  zipf <- file.path(DL_DIR, paste0(stem, "_all_r.zip"))
  csv  <- vapply(R_VALUES, function(r) es_path(es_file(cfg, r)), character(1))
  lanc <- vapply(R_VALUES, function(r) es_path(es_file(cfg, r, LANC_DIR)), character(1))
  lanc <- lanc[file.exists(lanc)]
  tmp <- file.path(tempdir(), stem); unlink(tmp, recursive = TRUE)
  dir.create(file.path(tmp, "effect_size"), recursive = TRUE)
  file.copy(csv, file.path(tmp, "effect_size"))
  if (length(lanc)) { dir.create(file.path(tmp, "effect_size_lancaster")); file.copy(lanc, file.path(tmp, "effect_size_lancaster")) }
  unlink(zipf)
  zip::zip(file.path(normalizePath(DL_DIR), basename(zipf)), files = list.files(tmp, recursive = TRUE), root = tmp)

  num <- function(m) data.frame(k = m$k, g = m$TE.random, ci_lower = m$lower.random, ci_upper = m$upper.random,
                                p = m$pval.random, tau2 = m$tau2, I2 = m$I2, pi_lower = m$lower.predict,
                                pi_upper = m$upper.predict)
  ledger <- do.call(rbind, lapply(R_VALUES, function(r) do.call(rbind, lapply(names(chain[[r]]), function(id)
    cbind(r = as.numeric(r) / 100, model = mlab(id), reported = id == reported_id(cfg), num(chain[[r]][[id]]))))))
  sheets <- list(models_by_r = ledger)
  if (isTRUE(cfg$lancaster)) {
    id <- reported_id(cfg)
    sheets$sensitivity <- do.call(rbind, lapply(R_VALUES, function(r) {
      dm <- model_data(cfg, load_es(es_file(cfg, r)), id)
      dl <- model_data(cfg, load_es(es_file(cfg, r, LANC_DIR)), id)
      rbind(cbind(r = as.numeric(r) / 100, analysis = "Main analysis", num(.fit(dm))),
            cbind(r = as.numeric(r) / 100, analysis = "Lancaster et al. (2005) - n / 3", num(.fit(dl))),
            cbind(r = as.numeric(r) / 100, analysis = "Pre/post designs only",
                  num(.fit(dplyr::filter(dm, !grepl("^Lancaster", Study))))))
    }))
  }
  xlsx <- file.path(DL_DIR, paste0(stem, "_results.xlsx"))
  writexl::write_xlsx(sheets, xlsx)
  c(zip = zipf, xlsx = xlsx)
}

chapter_resources <- function(cfg, chain) {
  f <- chapter_downloads(cfg, chain)
  csv <- function(r) paste0("../", es_file(cfg, r) |> sub(pattern = "^book/", replacement = ""))
  page_resources(list(
    c("GitHub repository", REPO_URL, "github"),
    c("Data - r = .50", csv(R_MAIN), "filetype-csv", csv("{r}")),
    c("Data - all r (ZIP)", paste0("../", f[["zip"]]), "file-earmark-zip"),
    c("R script", "../R/meta_Th1Th2.R", "file-earmark-code"),
    c("Figures (ZIP)", sprintf("../%s/%s_figures.zip", DL_DIR, sub("\\..*$", "", basename(knitr::current_input()))), "images"),
    c("Results (XLSX)", paste0("../", f[["xlsx"]]), "file-earmark-spreadsheet")
  ))
}
