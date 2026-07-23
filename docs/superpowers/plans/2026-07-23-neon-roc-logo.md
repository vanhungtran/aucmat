# Neon ROC Flow Logo Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the crowded logo with a vibrant, modern matrix-plus-ROC identity and generate every vector, PNG, and favicon asset reproducibly from R.

**Architecture:** `logo/generate-logo.R` is a self-contained build tool. It assembles SVG from deterministic geometry, writes three vector masters, rasterizes exact output sizes with `rsvg`, and creates the multi-resolution ICO with `magick`. A standalone base-R verification script checks the SVG semantics, output manifest, file signatures, and dimensions.

**Tech Stack:** R 4.x, base R, SVG 1.1, `rsvg`, `magick`

---

## File structure

- Modify: `logo/generate-logo.R` — palette, SVG builders, raster export, ICO export, and command-line entry point.
- Create: `logo/test-logo-generation.R` — isolated generator contract and output validation.
- Regenerate: `logo/aucmat_emblem.svg`, `logo/aucmat_emblem.png`, `logo/aucmat_logo.svg`, `logo/aucmat_logo.png`, `logo/aucmat_hex.svg`, `logo/aucmat_hex.png`.
- Regenerate: `man/figures/logo.png`, `docs/logo.png`.
- Regenerate: `pkgdown/favicon/favicon.svg`, `pkgdown/favicon/favicon.ico`, `pkgdown/favicon/favicon-96x96.png`, `pkgdown/favicon/apple-touch-icon.png`, `pkgdown/favicon/web-app-manifest-192x192.png`, `pkgdown/favicon/web-app-manifest-512x512.png`.

### Task 1: Lock the generator contract with a failing verification script

**Files:**
- Create: `logo/test-logo-generation.R`
- Test: `logo/test-logo-generation.R`

- [ ] **Step 1: Create the verification script**

```r
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "logo/test-logo-generation.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

generator <- new.env(parent = baseenv())
sys.source(file.path(repo_root, "logo", "generate-logo.R"), envir = generator)

required_functions <- c(
  "logo_palette", "build_emblem_svg", "build_logo_svg", "build_hex_svg",
  "output_manifest", "generate_logo_files"
)
stopifnot(all(vapply(required_functions, exists, logical(1), envir = generator, inherits = FALSE)))

emblem <- generator$build_emblem_svg()
full_logo <- generator$build_logo_svg()
hex <- generator$build_hex_svg()

stopifnot(
  grepl("<svg", emblem, fixed = TRUE),
  length(gregexpr("<rect", emblem, fixed = TRUE)[[1]]) >= 25L,
  grepl("id=\"roc-glow\"", emblem, fixed = TRUE),
  grepl("id=\"roc-gradient\"", emblem, fixed = TRUE),
  grepl("ROC analysis for biomarker matrices", full_logo, fixed = TRUE),
  grepl("aucmat", hex, fixed = TRUE),
  !grepl("AUC analysis for biomarker matrices", hex, fixed = TRUE)
)

build_root <- tempfile("aucmat-logo-")
dir.create(build_root)
generator$generate_logo_files(build_root)

manifest <- generator$output_manifest(build_root)
stopifnot(all(file.exists(manifest$path)))
stopifnot(all(file.info(manifest$path)$size > 100L))

png_rows <- manifest[manifest$format == "png", ]
for (i in seq_len(nrow(png_rows))) {
  info <- magick::image_info(magick::image_read(png_rows$path[[i]]))
  stopifnot(
    info$width[[1]] == png_rows$width[[i]],
    info$height[[1]] == png_rows$height[[i]]
  )
}

svg_rows <- manifest[manifest$format == "svg", ]
stopifnot(all(vapply(svg_rows$path, function(path) {
  identical(readBin(path, "raw", n = 4L), charToRaw("<svg"))
}, logical(1))))

ico_path <- file.path(build_root, "pkgdown", "favicon", "favicon.ico")
stopifnot(identical(readBin(ico_path, "raw", n = 4L), as.raw(c(0, 0, 1, 0))))

message("All logo generation checks passed.")
```

- [ ] **Step 2: Run the verification script and confirm it fails**

Run:

```powershell
& 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe' logo/test-logo-generation.R
```

Expected: failure because `build_emblem_svg` and the other generator functions do not exist in the current script.

- [ ] **Step 3: Commit the failing contract**

```powershell
git add -- logo/test-logo-generation.R
git commit -m "test: define reproducible logo generation contract"
```

### Task 2: Implement deterministic SVG builders

**Files:**
- Modify: `logo/generate-logo.R`
- Test: `logo/test-logo-generation.R`

- [ ] **Step 1: Replace the package-dependent ggplot generator with pure SVG helpers**

Define the public generator interface exactly as follows:

```r
logo_palette <- function() {
  c(
    ink = "#102A56", navy = "#173D72", indigo = "#4B45C6",
    violet = "#8559E8", blue = "#238AD2", cyan = "#2CD5DB",
    green = "#25B88A", yellow = "#F6C84A", orange = "#F58A35",
    coral = "#EF5D68", ice = "#F5FAFF", mist = "#DFEEFA"
  )
}

svg_document <- function(width, height, body, defs = "") {
  paste0(
    '<svg xmlns="http://www.w3.org/2000/svg" width="', width,
    '" height="', height, '" viewBox="0 0 ', width, " ", height,
    '" role="img">', defs, body, "</svg>"
  )
}

xml_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}
```

Create `emblem_parts()` to return:

- a 5 x 5 set of rounded 118 x 118 cell rectangles at x/y positions `180 + index * 132`;
- cell colors selected by `(row + column - 2) %% 7 + 1` from violet, blue, cyan, green, yellow, orange, coral;
- a translucent navy AUC fill bounded by `M 168 806 C 205 640 262 468 372 340 C 520 190 690 160 836 146 L 836 836 L 168 836 Z`;
- the ROC stroke path `M 168 806 C 205 640 262 468 372 340 C 520 190 690 160 836 146`;
- four circular ROC knots at `(208,650)`, `(294,444)`, `(480,260)`, and `(704,176)`.

Create `svg_defs()` containing:

```r
paste0(
  "<defs>",
  '<radialGradient id="surface" cx="30%" cy="20%" r="85%">',
  '<stop offset="0" stop-color="#FFFFFF"/><stop offset="1" stop-color="#E8F4FF"/>',
  "</radialGradient>",
  '<linearGradient id="roc-gradient" x1="0" y1="1" x2="1" y2="0">',
  '<stop offset="0" stop-color="#173D72"/><stop offset=".55" stop-color="#4B45C6"/>',
  '<stop offset="1" stop-color="#2CD5DB"/></linearGradient>',
  '<linearGradient id="word-gradient" x1="0" x2="1">',
  '<stop offset="0" stop-color="#238AD2"/><stop offset=".5" stop-color="#2CD5DB"/>',
  '<stop offset="1" stop-color="#8559E8"/></linearGradient>',
  '<linearGradient id="hex-edge" x1="0" y1="0" x2="1" y2="1">',
  '<stop offset="0" stop-color="#173D72"/><stop offset=".55" stop-color="#4B45C6"/>',
  '<stop offset="1" stop-color="#2CD5DB"/></linearGradient>',
  '<filter id="panel-shadow" x="-20%" y="-20%" width="140%" height="150%">',
  '<feDropShadow dx="0" dy="18" stdDeviation="20" flood-color="#102A56" flood-opacity=".18"/>',
  "</filter>",
  '<filter id="roc-glow" x="-30%" y="-30%" width="160%" height="160%">',
  '<feGaussianBlur stdDeviation="15" result="blur"/>',
  '<feFlood flood-color="#2CD5DB" flood-opacity=".42"/>',
  '<feComposite in2="blur" operator="in"/><feMerge><feMergeNode/><feMergeNode in="SourceGraphic"/></feMerge>',
  "</filter>",
  "</defs>"
)
```

Implement:

```r
build_emblem_svg <- function() {
  p <- logo_palette()
  parts <- emblem_parts()
  body <- paste0(
    '<rect width="1000" height="1000" rx="170" fill="url(#surface)"/>',
    '<circle cx="780" cy="180" r="115" fill="#2CD5DB" opacity=".10"/>',
    '<circle cx="155" cy="780" r="145" fill="#8559E8" opacity=".09"/>',
    '<g filter="url(#panel-shadow)">',
    '<rect x="128" y="128" width="744" height="744" rx="82" fill="#FFFFFF" fill-opacity=".82"/>',
    parts$cells,
    '<path d="', parts$area, '" fill="#173D72" opacity=".13"/>',
    '<path d="', parts$path, '" fill="none" stroke="#FFFFFF" stroke-width="56" stroke-linecap="round"/>',
    '<path d="', parts$path, '" fill="none" stroke="#2CD5DB" stroke-opacity=".50" stroke-width="42" ',
    'stroke-linecap="round" filter="url(#roc-glow)"/>',
    '<path d="', parts$path, '" fill="none" stroke="url(#roc-gradient)" stroke-width="30" ',
    'stroke-linecap="round" stroke-linejoin="round"/>',
    parts$knots,
    "</g>"
  )
  svg_document(1000, 1000, body, svg_defs())
}
```

`build_logo_svg()` must embed the emblem in a `1000 x 1280` document at
`x=90, y=20, width=820, height=820`, then place `auc` and `mat` on one baseline
at y=1015 and the descriptor at y=1105. `build_hex_svg()` must embed the emblem
at `x=126, y=120, width=748, height=748` inside a `1000 x 1155` hexagon with
points `500,26 954,288 954,812 500,1074 46,812 46,288`, a 32-pixel gradient
border, and centered `aucmat` text at y=955.

- [ ] **Step 2: Parse the script to catch R syntax errors**

Run:

```powershell
& 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe' -e "parse(file='logo/generate-logo.R'); cat('parse ok\n')"
```

Expected: `parse ok`.

- [ ] **Step 3: Run the contract and confirm the builder assertions advance to output generation**

Run:

```powershell
& 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe' logo/test-logo-generation.R
```

Expected: failure at `generate_logo_files` until Task 3 implements the exporter.

- [ ] **Step 4: Commit the vector builders**

```powershell
git add -- logo/generate-logo.R
git commit -m "feat: build Neon ROC Flow vector artwork"
```

### Task 3: Implement complete asset export

**Files:**
- Modify: `logo/generate-logo.R`
- Test: `logo/test-logo-generation.R`
- Regenerate: all paths listed in the file structure

- [ ] **Step 1: Add the exact output manifest**

```r
output_manifest <- function(root_dir = ".") {
  data.frame(
    path = file.path(root_dir, c(
      "logo/aucmat_emblem.svg", "logo/aucmat_logo.svg", "logo/aucmat_hex.svg",
      "logo/aucmat_emblem.png", "logo/aucmat_logo.png", "logo/aucmat_hex.png",
      "man/figures/logo.png", "docs/logo.png",
      "pkgdown/favicon/favicon.svg",
      "pkgdown/favicon/favicon-96x96.png",
      "pkgdown/favicon/apple-touch-icon.png",
      "pkgdown/favicon/web-app-manifest-192x192.png",
      "pkgdown/favicon/web-app-manifest-512x512.png",
      "pkgdown/favicon/favicon.ico"
    )),
    format = c(rep("svg", 3), rep("png", 5), "svg", rep("png", 4), "ico"),
    width = c(1000, 1000, 1000, 2400, 2200, 691, 1200, 1200, 1000, 96, 180, 192, 512, NA),
    height = c(1000, 1280, 1155, 2400, 2800, 800, 1200, 1200, 1000, 96, 180, 192, 512, NA),
    stringsAsFactors = FALSE
  )
}
```

- [ ] **Step 2: Add SVG writing and PNG rendering**

```r
write_utf8 <- function(text, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw(enc2utf8(text)), con)
}

render_png <- function(svg_path, png_path, width, height) {
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
  rsvg::rsvg_png(svg_path, png_path, width = width, height = height)
}
```

- [ ] **Step 3: Add the orchestrator and CLI entry point**

```r
generate_logo_files <- function(root_dir = ".") {
  missing <- c("rsvg", "magick")[!vapply(c("rsvg", "magick"), requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Install required logo packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  root_dir <- normalizePath(root_dir, mustWork = TRUE)
  emblem_svg <- build_emblem_svg()
  logo_svg <- build_logo_svg()
  hex_svg <- build_hex_svg()

  emblem_path <- file.path(root_dir, "logo", "aucmat_emblem.svg")
  logo_path <- file.path(root_dir, "logo", "aucmat_logo.svg")
  hex_path <- file.path(root_dir, "logo", "aucmat_hex.svg")
  favicon_path <- file.path(root_dir, "pkgdown", "favicon", "favicon.svg")
  write_utf8(emblem_svg, emblem_path)
  write_utf8(logo_svg, logo_path)
  write_utf8(hex_svg, hex_path)
  write_utf8(emblem_svg, favicon_path)

  jobs <- list(
    c(emblem_path, file.path(root_dir, "logo", "aucmat_emblem.png"), 2400, 2400),
    c(logo_path, file.path(root_dir, "logo", "aucmat_logo.png"), 2200, 2800),
    c(hex_path, file.path(root_dir, "logo", "aucmat_hex.png"), 691, 800),
    c(emblem_path, file.path(root_dir, "man", "figures", "logo.png"), 1200, 1200),
    c(emblem_path, file.path(root_dir, "docs", "logo.png"), 1200, 1200),
    c(emblem_path, file.path(root_dir, "pkgdown", "favicon", "favicon-96x96.png"), 96, 96),
    c(emblem_path, file.path(root_dir, "pkgdown", "favicon", "apple-touch-icon.png"), 180, 180),
    c(emblem_path, file.path(root_dir, "pkgdown", "favicon", "web-app-manifest-192x192.png"), 192, 192),
    c(emblem_path, file.path(root_dir, "pkgdown", "favicon", "web-app-manifest-512x512.png"), 512, 512)
  )
  for (job in jobs) render_png(job[[1]], job[[2]], as.integer(job[[3]]), as.integer(job[[4]]))

  icon <- magick::image_read(c(
    file.path(root_dir, "pkgdown", "favicon", "favicon-96x96.png"),
    file.path(root_dir, "pkgdown", "favicon", "web-app-manifest-192x192.png"),
    file.path(root_dir, "pkgdown", "favicon", "web-app-manifest-512x512.png")
  ))
  magick::image_write(icon, file.path(root_dir, "pkgdown", "favicon", "favicon.ico"), format = "ico")
  invisible(output_manifest(root_dir))
}

if (sys.nframe() == 0L) {
  generate_logo_files(normalizePath(".", mustWork = TRUE))
  message("All Neon ROC Flow logo assets generated.")
}
```

- [ ] **Step 4: Run the full isolated verification**

Run:

```powershell
& 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe' logo/test-logo-generation.R
```

Expected: `All logo generation checks passed.`

- [ ] **Step 5: Regenerate repository assets**

Run:

```powershell
& 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe' logo/generate-logo.R
```

Expected: `All Neon ROC Flow logo assets generated.`

- [ ] **Step 6: Commit generator and all generated outputs**

```powershell
git add -- logo/generate-logo.R logo/test-logo-generation.R logo/aucmat_emblem.svg logo/aucmat_emblem.png logo/aucmat_logo.svg logo/aucmat_logo.png logo/aucmat_hex.svg logo/aucmat_hex.png man/figures/logo.png pkgdown/favicon/favicon.svg pkgdown/favicon/favicon.ico pkgdown/favicon/favicon-96x96.png pkgdown/favicon/apple-touch-icon.png pkgdown/favicon/web-app-manifest-192x192.png pkgdown/favicon/web-app-manifest-512x512.png
git add -f -- docs/logo.png
git commit -m "feat: redesign logo with Neon ROC Flow identity"
```

### Task 4: Visual and repository verification

**Files:**
- Inspect: `logo/aucmat_emblem.png`
- Inspect: `logo/aucmat_logo.png`
- Inspect: `logo/aucmat_hex.png`
- Inspect: `pkgdown/favicon/favicon-96x96.png`

- [ ] **Step 1: Inspect the four representative raster outputs**

Open the four files at native resolution. Confirm:

- all 25 matrix cells are distinct;
- the ROC curve reads from lower-left to upper-right;
- glow does not obscure the white separator or curve knots;
- the wordmark does not overlap the emblem;
- the hex has even optical padding;
- `aucmat` is legible at hex size;
- the 96-pixel favicon still reads as matrix plus ROC.

- [ ] **Step 2: Run the package test suite**

Run:

```powershell
& 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe' -e "testthat::test_local(reporter='summary')"
```

Expected: all existing package tests pass.

- [ ] **Step 3: Confirm the intended git scope**

Run:

```powershell
git status --short
git diff --stat HEAD~1
```

Expected: the logo generator, its verification script, and generated logo/favicon assets are included; the pre-existing `Rplots.pdf` remains untracked and is not committed.
