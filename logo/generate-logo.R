# Neon ROC Flow logo generator -----------------------------------------------
#
# Run from any directory:
#   Rscript path/to/aucmat/logo/generate-logo.R
#
# Development dependencies:
#   install.packages(c("rsvg", "magick"))

logo_palette <- function() {
  c(
    ink = "#102A56",
    navy = "#173D72",
    indigo = "#4B45C6",
    violet = "#8559E8",
    blue = "#238AD2",
    cyan = "#2CD5DB",
    green = "#25B88A",
    yellow = "#F6C84A",
    orange = "#F58A35",
    coral = "#EF5D68",
    ice = "#F5FAFF",
    mist = "#DFEEFA"
  )
}

xml_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

svg_document <- function(width, height, body, defs = "") {
  paste0(
    '<svg xmlns="http://www.w3.org/2000/svg" ',
    'width="', width, '" height="', height, '" ',
    'viewBox="0 0 ', width, " ", height, '" ',
    'role="img">',
    defs,
    body,
    "</svg>"
  )
}

svg_defs <- function() {
  paste0(
    "<defs>",
    '<radialGradient id="surface" cx="28%" cy="18%" r="92%">',
    '<stop offset="0" stop-color="#FFFFFF"/>',
    '<stop offset=".62" stop-color="#F5FAFF"/>',
    '<stop offset="1" stop-color="#E4F2FF"/>',
    "</radialGradient>",
    '<linearGradient id="roc-gradient" x1="0" y1="1" x2="1" y2="0">',
    '<stop offset="0" stop-color="#173D72"/>',
    '<stop offset=".55" stop-color="#4B45C6"/>',
    '<stop offset="1" stop-color="#2CD5DB"/>',
    "</linearGradient>",
    '<linearGradient id="word-gradient" x1="0" x2="1">',
    '<stop offset="0" stop-color="#238AD2"/>',
    '<stop offset=".48" stop-color="#2CD5DB"/>',
    '<stop offset="1" stop-color="#8559E8"/>',
    "</linearGradient>",
    '<linearGradient id="hex-edge" x1="0" y1="0" x2="1" y2="1">',
    '<stop offset="0" stop-color="#173D72"/>',
    '<stop offset=".52" stop-color="#4B45C6"/>',
    '<stop offset="1" stop-color="#2CD5DB"/>',
    "</linearGradient>",
    '<linearGradient id="hex-surface" x1="0" y1="0" x2="1" y2="1">',
    '<stop offset="0" stop-color="#FFFFFF"/>',
    '<stop offset=".68" stop-color="#F3F8FF"/>',
    '<stop offset="1" stop-color="#E7F4FF"/>',
    "</linearGradient>",
    '<filter id="panel-shadow" x="-20%" y="-20%" width="140%" height="155%">',
    '<feDropShadow dx="0" dy="18" stdDeviation="20" ',
    'flood-color="#102A56" flood-opacity=".18"/>',
    "</filter>",
    '<filter id="roc-glow" x="-30%" y="-30%" width="160%" height="160%">',
    '<feGaussianBlur stdDeviation="15" result="blur"/>',
    '<feFlood flood-color="#2CD5DB" flood-opacity=".42"/>',
    '<feComposite in2="blur" operator="in"/>',
    "<feMerge>",
    "<feMergeNode/>",
    '<feMergeNode in="SourceGraphic"/>',
    "</feMerge>",
    "</filter>",
    '<filter id="hex-shadow" x="-20%" y="-20%" width="140%" height="150%">',
    '<feDropShadow dx="0" dy="22" stdDeviation="18" ',
    'flood-color="#102A56" flood-opacity=".22"/>',
    "</filter>",
    "</defs>"
  )
}

emblem_parts <- function() {
  p <- logo_palette()
  cell_colors <- unname(p[c(
    "violet", "blue", "cyan", "green",
    "yellow", "orange", "coral"
  )])

  cells <- character(25)
  cell_index <- 1L
  for (row in 0:4) {
    for (column in 0:4) {
      x <- 180 + column * 132
      y <- 180 + row * 132
      color_index <- (row + column) %% length(cell_colors) + 1L
      opacity <- c(.92, .98, .86)[(row * 5 + column) %% 3 + 1L]
      cells[[cell_index]] <- paste0(
        '<rect x="', x, '" y="', y,
        '" width="112" height="112" rx="22" ',
        'fill="', cell_colors[[color_index]],
        '" fill-opacity="', opacity,
        '" stroke="#FFFFFF" stroke-opacity=".72" stroke-width="7"/>'
      )
      cell_index <- cell_index + 1L
    }
  }

  path <- paste0(
    "M 168 806 ",
    "C 205 640 262 468 372 340 ",
    "C 520 190 690 160 836 146"
  )
  area <- paste0(
    path,
    " L 836 836 L 168 836 Z"
  )

  knots_xy <- matrix(
    c(
      208, 650,
      294, 444,
      480, 260,
      704, 176
    ),
    ncol = 2,
    byrow = TRUE
  )
  knots <- apply(knots_xy, 1, function(xy) {
    paste0(
      '<circle cx="', xy[[1]], '" cy="', xy[[2]],
      '" r="23" fill="#FFFFFF" stroke="#173D72" stroke-width="9"/>'
    )
  })

  list(
    cells = paste0(cells, collapse = ""),
    path = path,
    area = area,
    knots = paste0(knots, collapse = "")
  )
}

emblem_body <- function() {
  parts <- emblem_parts()
  paste0(
    '<rect width="1000" height="1000" rx="170" fill="url(#surface)"/>',
    '<circle cx="790" cy="174" r="120" fill="#2CD5DB" opacity=".11"/>',
    '<circle cx="156" cy="796" r="150" fill="#8559E8" opacity=".10"/>',
    '<path d="M 735 90 C 835 112 900 188 920 284" ',
    'fill="none" stroke="#FFFFFF" stroke-opacity=".72" stroke-width="18" ',
    'stroke-linecap="round"/>',
    '<g filter="url(#panel-shadow)">',
    '<rect x="128" y="128" width="744" height="744" rx="82" ',
    'fill="#FFFFFF" fill-opacity=".80" stroke="#FFFFFF" stroke-width="6"/>',
    parts$cells,
    '<path d="', parts$area, '" fill="#173D72" opacity=".13"/>',
    '<path d="', parts$path, '" fill="none" stroke="#FFFFFF" ',
    'stroke-width="62" stroke-linecap="round" stroke-linejoin="round"/>',
    '<path d="', parts$path, '" fill="none" stroke="#2CD5DB" ',
    'stroke-opacity=".52" stroke-width="48" stroke-linecap="round" ',
    'stroke-linejoin="round" filter="url(#roc-glow)"/>',
    '<path d="', parts$path, '" fill="none" stroke="url(#roc-gradient)" ',
    'stroke-width="32" stroke-linecap="round" stroke-linejoin="round"/>',
    parts$knots,
    "</g>"
  )
}

build_emblem_svg <- function() {
  svg_document(
    width = 1000,
    height = 1000,
    body = paste0(
      "<title>aucmat matrix and ROC emblem</title>",
      emblem_body()
    ),
    defs = svg_defs()
  )
}

wordmark_body <- function(y = 1024, descriptor_y = 1113) {
  paste0(
    '<text x="500" y="', y, '" text-anchor="middle" ',
    'font-family="Arial, sans-serif" font-size="154" ',
    'font-weight="800" letter-spacing="-9">',
    '<tspan fill="#102A56">auc</tspan>',
    '<tspan fill="url(#word-gradient)">mat</tspan>',
    "</text>",
    '<rect x="392" y="', y + 36,
    '" width="216" height="10" rx="5" fill="url(#word-gradient)"/>',
    '<text x="500" y="', descriptor_y, '" text-anchor="middle" ',
    'font-family="Arial, sans-serif" font-size="34" font-weight="600" ',
    'letter-spacing="2.2" fill="#58718F">',
    xml_escape("ROC analysis for biomarker matrices"),
    "</text>"
  )
}

build_logo_svg <- function() {
  body <- paste0(
    "<title>aucmat logo</title>",
    '<rect width="1000" height="1273" fill="#FFFFFF"/>',
    '<g transform="translate(90 20) scale(.82)">',
    emblem_body(),
    "</g>",
    wordmark_body()
  )
  svg_document(
    width = 1000,
    height = 1273,
    body = body,
    defs = svg_defs()
  )
}

build_hex_svg <- function() {
  points <- "500,26 954,288 954,812 500,1074 46,812 46,288"
  body <- paste0(
    "<title>aucmat hex sticker</title>",
    '<rect width="1000" height="1155" fill="none"/>',
    '<polygon points="', points, '" fill="url(#hex-surface)" ',
    'stroke="#FFFFFF" stroke-width="58" stroke-linejoin="round" ',
    'filter="url(#hex-shadow)"/>',
    '<polygon points="', points, '" fill="url(#hex-surface)" ',
    'stroke="url(#hex-edge)" stroke-width="32" stroke-linejoin="round"/>',
    '<g transform="translate(140 90) scale(.72)">',
    emblem_body(),
    "</g>",
    '<text x="500" y="956" text-anchor="middle" ',
    'font-family="Arial, sans-serif" font-size="126" ',
    'font-weight="800" letter-spacing="-7">',
    '<tspan fill="#102A56">auc</tspan>',
    '<tspan fill="url(#word-gradient)">mat</tspan>',
    "</text>",
    '<rect x="412" y="988" width="176" height="9" rx="4.5" ',
    'fill="url(#word-gradient)"/>'
  )
  svg_document(
    width = 1000,
    height = 1155,
    body = body,
    defs = svg_defs()
  )
}

output_manifest <- function(root_dir = ".") {
  data.frame(
    path = file.path(
      root_dir,
      c(
        "logo/aucmat_emblem.svg",
        "logo/aucmat_logo.svg",
        "logo/aucmat_hex.svg",
        "logo/aucmat_emblem.png",
        "logo/aucmat_logo.png",
        "logo/aucmat_hex.png",
        "man/figures/logo.png",
        "docs/logo.png",
        "pkgdown/favicon/favicon.svg",
        "pkgdown/favicon/favicon-96x96.png",
        "pkgdown/favicon/apple-touch-icon.png",
        "pkgdown/favicon/web-app-manifest-192x192.png",
        "pkgdown/favicon/web-app-manifest-512x512.png",
        "pkgdown/favicon/favicon.ico"
      )
    ),
    format = c(
      rep("svg", 3),
      rep("png", 5),
      "svg",
      rep("png", 4),
      "ico"
    ),
    width = c(
      1000, 1000, 1000,
      2400, 2200, 691, 1200, 1200,
      1000,
      96, 180, 192, 512,
      NA
    ),
    height = c(
      1000, 1273, 1155,
      2400, 2800, 800, 1200, 1200,
      1000,
      96, 180, 192, 512,
      NA
    ),
    stringsAsFactors = FALSE
  )
}

write_utf8 <- function(text, path) {
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )
  connection <- file(path, open = "wb")
  on.exit(close(connection), add = TRUE)
  writeBin(charToRaw(enc2utf8(text)), connection)
}

render_png <- function(svg_path, png_path, width, height) {
  dir.create(
    dirname(png_path),
    recursive = TRUE,
    showWarnings = FALSE
  )
  rsvg::rsvg_png(
    svg = svg_path,
    file = png_path,
    width = width,
    height = height
  )
}

generate_favicon_ico <- function(source_path, output_path) {
  source <- magick::image_read(source_path)
  sizes <- c(16L, 32L, 48L)
  frames <- lapply(sizes, function(size) {
    magick::image_resize(
      source,
      paste0(size, "x", size, "!")
    )
  })
  icon <- do.call(c, frames)
  magick::image_write(
    icon,
    path = output_path,
    format = "ico"
  )
}

generate_logo_files <- function(root_dir = ".") {
  packages <- c("rsvg", "magick")
  available <- vapply(
    packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
  if (any(!available)) {
    stop(
      "Install required logo packages: ",
      paste(packages[!available], collapse = ", "),
      call. = FALSE
    )
  }

  root_dir <- normalizePath(root_dir, mustWork = TRUE)
  emblem_path <- file.path(
    root_dir,
    "logo",
    "aucmat_emblem.svg"
  )
  logo_path <- file.path(
    root_dir,
    "logo",
    "aucmat_logo.svg"
  )
  hex_path <- file.path(
    root_dir,
    "logo",
    "aucmat_hex.svg"
  )
  favicon_path <- file.path(
    root_dir,
    "pkgdown",
    "favicon",
    "favicon.svg"
  )

  write_utf8(build_emblem_svg(), emblem_path)
  write_utf8(build_logo_svg(), logo_path)
  write_utf8(build_hex_svg(), hex_path)
  write_utf8(build_emblem_svg(), favicon_path)

  jobs <- list(
    list(
      emblem_path,
      file.path(root_dir, "logo", "aucmat_emblem.png"),
      2400L, 2400L
    ),
    list(
      logo_path,
      file.path(root_dir, "logo", "aucmat_logo.png"),
      2200L, 2800L
    ),
    list(
      hex_path,
      file.path(root_dir, "logo", "aucmat_hex.png"),
      691L, 800L
    ),
    list(
      emblem_path,
      file.path(root_dir, "man", "figures", "logo.png"),
      1200L, 1200L
    ),
    list(
      emblem_path,
      file.path(root_dir, "docs", "logo.png"),
      1200L, 1200L
    ),
    list(
      emblem_path,
      file.path(
        root_dir,
        "pkgdown",
        "favicon",
        "favicon-96x96.png"
      ),
      96L, 96L
    ),
    list(
      emblem_path,
      file.path(
        root_dir,
        "pkgdown",
        "favicon",
        "apple-touch-icon.png"
      ),
      180L, 180L
    ),
    list(
      emblem_path,
      file.path(
        root_dir,
        "pkgdown",
        "favicon",
        "web-app-manifest-192x192.png"
      ),
      192L, 192L
    ),
    list(
      emblem_path,
      file.path(
        root_dir,
        "pkgdown",
        "favicon",
        "web-app-manifest-512x512.png"
      ),
      512L, 512L
    )
  )

  for (job in jobs) {
    render_png(
      svg_path = job[[1]],
      png_path = job[[2]],
      width = job[[3]],
      height = job[[4]]
    )
  }

  generate_favicon_ico(
    source_path = file.path(
      root_dir,
      "pkgdown",
      "favicon",
      "web-app-manifest-512x512.png"
    ),
    output_path = file.path(
      root_dir,
      "pkgdown",
      "favicon",
      "favicon.ico"
    )
  )

  invisible(output_manifest(root_dir))
}

script_repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (!length(file_arg)) {
    return(normalizePath(".", mustWork = TRUE))
  }
  script_path <- sub("^--file=", "", file_arg[[1]])
  normalizePath(
    file.path(dirname(script_path), ".."),
    mustWork = TRUE
  )
}

if (sys.nframe() == 0L) {
  generate_logo_files(script_repo_root())
  message("All Neon ROC Flow logo assets generated.")
}
