# =========================
# ROBINS-I plots (A/B/C) + patchwork
# =========================

## Packages laden
library(tidyverse)
library(patchwork)
library(scales)
library(readxl)

# -------------------------
## Data laden
# -------------------------
data_RoB <- read_xlsx(
  path = "04_risk_of_bias/results.xlsx",
  sheet = "Overall",
  range = "A1:H10"
)

data_ConFac <- read_xlsx(
  path = "04_risk_of_bias/results.xlsx",
  sheet = "ConFac",
  range = "A1:G10"
)

data_ConFac_notes <- tribble(
  ~Item,   ~Note,
  "D1.1",  "No exercise 24h-48h pre-intervention?",
  "D1.2",  "Standardized food intake among participants (e.g. all fasted or standardized breakfast)?",
  "D1.3",  "No caffeine, acute medication, or alcohol 24h before exercise?",
  "D1.4",  "Identical warm-up between participants?",
  "D1.5",  "good way to assess the controlled intensity (HR control, VO2, VT)",
  "D1.6",  "Homogeneity of the group (sex, age, BMI, training status)"
)

# =========================================================================
## EINSTELLUNGEN (editierbare Vektoren am Anfang)
# =========================================================================

## 1. STUDIENREIHENFOLGE (alphabetisch absteigend: Z→A, A ganz unten)
study_order <- sort(unique(data_RoB$Study), decreasing = TRUE)  # Z→A

## 2. PLOT C - Domain-Reihenfolge (editierbar, "Overall" ganz unten)
domain_order_C <- c(
  "Bias due to confounding",
  "Bias due to selection of participants",
  "Bias due to deviations from intended interventions",
  "Bias due to missing data",
  "Bias in measurement of outcomes",
  "Bias in selection of the reported result",
  "Overall risk of bias"  # <-- ganz unten
)

## 3. BALKEN-REIHENFOLGE in Plot C (Low → Moderate → Serious, links→rechts)
risk_levels_bar <- c("Low", "Moderate", "Serious")
risk_levels_all <- c("Low", "Moderate", "Serious", "Critical")

## 4. Farben (robvis-Standard)
robvis_cols <- c(
  Low      = "#00C000",
  Moderate = "#E6E600",
  Serious  = "#E31A1C",
  Critical = "#7F0000"
)

## 5. Critical anzeigen ja/nein?
show_critical <- FALSE
active_levels <- if (show_critical) risk_levels_all else risk_levels_bar
active_cols   <- robvis_cols[active_levels]

## 6. Schriftgrößen
base_font   <- 12
legend_font <- 10

# =========================================================================
## DOMAIN-LABELS
# =========================================================================
domain_labels <- c(
  Overall = "Overall risk of bias",
  D1      = "Bias due to confounding",
  D3      = "Bias due to selection of participants",
  D4      = "Bias due to deviations from intended interventions",
  D5      = "Bias due to missing data",
  D6      = "Bias in measurement of outcomes",
  D7      = "Bias in selection of the reported result"
)

# =========================================================================
## THEMES
# =========================================================================

common_theme <- theme(
  text         = element_text(size = base_font, color = "black"),
  axis.text    = element_text(size = base_font),
  legend.text  = element_text(size = legend_font),
  legend.title = element_text(size = legend_font),
  plot.title   = element_text(face = "bold", size = base_font + 2, hjust = 0),
  plot.tag     = element_text(face = "bold", size = base_font + 4)
)

top_axis_theme <- theme(
  panel.grid        = element_blank(),
  axis.text         = element_text(color = "black"),
  axis.text.x       = element_text(margin = margin(b = 6)),
  axis.ticks        = element_blank(),
  axis.line         = element_blank(),
  axis.line.x.top   = element_line(color = "black", linewidth = 0.8),
  axis.line.y       = element_blank(),
  plot.margin       = margin(10, 10, 10, 10),
  legend.position   = "bottom",
  legend.direction  = "horizontal",
  legend.background = element_blank(),
  legend.key        = element_blank()
)

# =========================================================================
## PLOT A: RoB Matrix (Punkte je Studie × Domain)
## Studienreihenfolge: Z→A (Z oben, A unten)
# =========================================================================

df_matrix <- data_RoB |>
  pivot_longer(cols = c(D1, D3, D4, D5, D6, D7, Overall),
               names_to = "Domain", values_to = "Risk") |>
  mutate(
    Risk   = factor(trimws(Risk), levels = risk_levels_all),
    Domain = factor(Domain, levels = c("D1", "D3", "D4", "D5", "D6", "D7", "Overall")),
    Study  = factor(Study, levels = study_order)  # Z→A
  )

p_A <- ggplot(df_matrix, aes(x = Domain, y = Study)) +
  geom_point(aes(fill = Risk), shape = 21, size = 9,
             color = "black", stroke = 1, show.legend = TRUE) +
  scale_fill_manual(
    values = robvis_cols,
    breaks = risk_levels_all,
    drop   = !show_critical
  ) +
  scale_x_discrete(position = "top") +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_minimal() +
  top_axis_theme

# =========================================================================
## PLOT B: Confounding Factors Matrix (Punkte je Studie × Item)
## Studienreihenfolge: synchron mit Plot A
# =========================================================================

df_conf_long <- data_ConFac |>
  pivot_longer(cols = -Study, names_to = "Item", values_to = "Answer") |>
  mutate(
    Item   = factor(Item, levels = data_ConFac_notes$Item),
    Study  = factor(Study, levels = study_order),  # Z→A (synchron mit A)
    Answer = factor(Answer, levels = c("Unsure", "Sure"))
  )

p_B <- ggplot(df_conf_long, aes(x = Item, y = Study)) +
  geom_point(aes(fill = Answer), shape = 21, size = 9,
             color = "black", stroke = 1, show.legend = TRUE) +
  scale_x_discrete(position = "top") +
  scale_fill_manual(
    values = c(Unsure = "orange1", Sure = "#00C000"),
    breaks = c("Sure", "Unsure"),
    drop   = TRUE
  ) +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_minimal() +
  top_axis_theme

# =========================================================================
## PLOT C: Risk of Bias Balkendiagramm (proportional)
## Domain-Reihenfolge: Overall visuell ganz unten / als letztes
## Balken-Reihenfolge: Low → Moderate → Critical
# =========================================================================

risk_levels_C <- c("Low", "Moderate", "Serious")
active_cols_C <- robvis_cols[risk_levels_C]

df_long_C <- data_RoB |>
  pivot_longer(
    cols = c(D1, D3, D4, D5, D6, D7, Overall),
    names_to = "Domain",
    values_to = "Risk"
  ) |>
  mutate(
    Risk = factor(trimws(Risk), levels = risk_levels_C),

    Domain = recode(Domain, !!!domain_labels),

    # wichtig: wegen coord_flip() muss hier rev() rein,
    # damit Overall im fertigen Plot unten steht
    Domain = factor(Domain, levels = rev(domain_order_C))
  ) |>
  filter(!is.na(Risk))

p_C <- ggplot(df_long_C, aes(x = Domain, fill = Risk)) +
  geom_bar(
    position = position_fill(reverse = TRUE),
    width = 0.55,
    color = "black",
    linewidth = 0.5
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, .25, .50, .75, 1),
    labels = percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = active_cols_C,
    breaks = risk_levels_C,
    limits = risk_levels_C,
    drop   = FALSE
  ) +
  guides(
    fill = guide_legend(
      reverse = FALSE,
      override.aes = list(colour = "black")
    )
  ) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_classic() +
  theme(
    panel.grid      = element_blank(),
    axis.line.y     = element_blank(),
    axis.line.x     = element_line(color = "black", linewidth = 0.8),
    axis.ticks.y    = element_blank(),
    axis.ticks.x    = element_line(color = "black", linewidth = 0.8),
    legend.position = "bottom",
    legend.key      = element_rect(fill = "white", colour = "black"),
    plot.margin     = margin(10, 30, 10, 10)
  )

# =========================================================================
## FINAL LAYOUT: A+B oben, C unten
# =========================================================================

final_plot <- ((p_A | p_B) + plot_layout(widths = c(7, 6)) & common_theme) /
  (p_C & common_theme) +
  plot_layout(heights = c(2, 1.2)) +
  plot_annotation(tag_levels = list(c("A", "B", "C")))

# =========================================================================
## EXPORT
# =========================================================================
#ggsave("04_risk_of_bias/ROBINS-I V4.1.png", final_plot, width = 14, height = 10, dpi = 600, bg = "white")
