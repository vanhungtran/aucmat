# aucmat ggplot2 theme

A clean, professional theme for all aucmat plots. Uses Fira Sans if
available via showtext, otherwise falls back to system defaults.

## Usage

``` r
theme_aucmat(base_size = 12, base_family = "")
```

## Arguments

- base_size:

  Base font size (default 12).

- base_family:

  Font family. If `"fira"`, checks for Fira Sans registration via
  showtext.

## Value

A ggplot2 theme object.
