# generate_manuscript_pdf.R
# Compile the manuscript LaTeX file to PDF.
# Requires: A LaTeX distribution (TeX Live, MiKTeX, or TinyTeX)
#           with pdflatex and bibtex.
#
# Usage: Rscript data-raw/generate_manuscript_pdf.R
#
# The .tex source lives at manuscript/manuscript.tex
# The PDF will be placed at manuscript/manuscript.pdf

tex_file <- "manuscript/manuscript.tex"
tex_dir  <- dirname(tex_file)

if (!file.exists(tex_file)) {
  stop("Manuscript .tex not found at: ", normalizePath(tex_file, mustWork = FALSE))
}

# Check for pdflatex
if (Sys.which("pdflatex") == "") {
  stop("pdflatex not found. Install a LaTeX distribution (e.g., TinyTeX: tinytex::install_tinytex())")
}

old_wd <- getwd()
on.exit(setwd(old_wd))
setwd(tex_dir)

message("Compiling ", basename(tex_file), " ...")

# pdflatex -> bibtex -> pdflatex -> pdflatex
for (i in 1:2) {
  system2("pdflatex", c("-interaction=nonstopmode", basename(tex_file)))
}
system2("bibtex", basename(tools::file_path_sans_ext(basename(tex_file))))
for (i in 1:2) {
  system2("pdflatex", c("-interaction=nonstopmode", basename(tex_file)))
}

# Clean up auxiliary files
aux_exts <- c("aux", "bbl", "blg", "log", "out", "toc")
for (ext in aux_exts) {
  f <- paste0(tools::file_path_sans_ext(basename(tex_file)), ".", ext)
  if (file.exists(f)) file.remove(f)
}

setwd(old_wd)

pdf_file <- file.path(tex_dir, paste0(tools::file_path_sans_ext(basename(tex_file)), ".pdf"))
if (file.exists(pdf_file)) {
  message("PDF generated: ", normalizePath(pdf_file))
} else {
  stop("PDF generation failed; check manuscript/manuscript.log for errors")
}
