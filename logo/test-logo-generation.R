args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) {
  sub("^--file=", "", file_arg[[1]])
} else {
  "logo/test-logo-generation.R"
}
repo_root <- normalizePath(
  file.path(dirname(script_path), ".."),
  mustWork = TRUE
)

generator <- new.env(parent = globalenv())
sys.source(
  file.path(repo_root, "logo", "generate-logo.R"),
  envir = generator
)

required_functions <- c(
  "logo_palette",
  "build_emblem_svg",
  "build_logo_svg",
  "build_hex_svg",
  "output_manifest",
  "generate_logo_files"
)
stopifnot(all(vapply(
  required_functions,
  exists,
  logical(1),
  envir = generator,
  inherits = FALSE
)))

emblem <- generator$build_emblem_svg()
full_logo <- generator$build_logo_svg()
hex <- generator$build_hex_svg()

stopifnot(
  grepl("<svg", emblem, fixed = TRUE),
  length(gregexpr("<rect", emblem, fixed = TRUE)[[1]]) >= 25L,
  grepl("id=\"roc-glow\"", emblem, fixed = TRUE),
  grepl("id=\"roc-gradient\"", emblem, fixed = TRUE),
  grepl(
    "ROC analysis for biomarker matrices",
    full_logo,
    fixed = TRUE
  ),
  grepl("aucmat", hex, fixed = TRUE),
  !grepl(
    "AUC analysis for biomarker matrices",
    hex,
    fixed = TRUE
  )
)

build_root <- tempfile("aucmat-logo-")
dir.create(build_root)
on.exit(
  unlink(build_root, recursive = TRUE, force = TRUE),
  add = TRUE
)
generator$generate_logo_files(build_root)

manifest <- generator$output_manifest(build_root)
stopifnot(all(file.exists(manifest$path)))
stopifnot(all(file.info(manifest$path)$size > 100L))

png_rows <- manifest[manifest$format == "png", ]
for (i in seq_len(nrow(png_rows))) {
  info <- magick::image_info(
    magick::image_read(png_rows$path[[i]])
  )
  stopifnot(
    info$width[[1]] == png_rows$width[[i]],
    info$height[[1]] == png_rows$height[[i]]
  )
}

svg_rows <- manifest[manifest$format == "svg", ]
stopifnot(all(vapply(
  svg_rows$path,
  function(path) {
    identical(
      readBin(path, "raw", n = 4L),
      charToRaw("<svg")
    )
  },
  logical(1)
)))

ico_path <- file.path(
  build_root,
  "pkgdown",
  "favicon",
  "favicon.ico"
)
stopifnot(identical(
  readBin(ico_path, "raw", n = 4L),
  as.raw(c(0, 0, 1, 0))
))

message("All logo generation checks passed.")
