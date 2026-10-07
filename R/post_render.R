# Post-render: zip the static figures of every chapter into _book/downloads/<chapter>_figures.zip
out <- Sys.getenv("QUARTO_PROJECT_OUTPUT_DIR", "_book")
dl  <- file.path(out, "downloads")
dir.create(dl, showWarnings = FALSE, recursive = TRUE)
for (d in Sys.glob(file.path(out, "qmd", "*_files", "figure-html"))) {
  chapter <- sub("_files$", "", basename(dirname(d)))
  files <- list.files(d, pattern = "\\.(png|svg|pdf)$")
  if (!length(files)) next
  zipf <- normalizePath(file.path(dl, paste0(chapter, "_figures.zip")), mustWork = FALSE)
  unlink(zipf)
  zip::zip(zipf, files = files, root = d)
}
