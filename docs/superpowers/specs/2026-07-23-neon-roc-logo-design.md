# Neon ROC Flow logo design

## Objective

Redesign the `aucmat` logo so it feels vibrant, modern, and immediately
recognizable while preserving the existing matrix-plus-ROC concept. The artwork
must be reproducible from code and remain legible as a README logo, package hex
sticker, and small favicon.

## Visual design

The primary emblem is a simplified 5 x 5 rounded biomarker matrix. Its cells
progress from violet through blue, cyan, green, yellow, orange, and coral,
creating energy without the visual noise of the existing 6 x 6 grid. A smooth,
upward ROC curve crosses the matrix. The curve uses a deep-navy core, a narrow
white separator, and a restrained cyan glow so it remains crisp rather than
looking blurred.

The matrix sits on a subtle blue-white radial background with a soft shadow.
Grid cells use modest corner radii and slight tonal variation to add depth.
Decorative effects are limited to the glow, gradient, and shadow; no sparkles,
3D extrusion, or extra chart annotations are included.

The standalone emblem contains no text. The full lockup places the `aucmat`
wordmark below the emblem, with `auc` in navy and `mat` in a cyan-to-violet
gradient. The descriptor reads `ROC analysis for biomarker matrices`. Text
never overlaps the icon.

## Hex sticker

The hex sticker uses a dark navy-to-indigo border and a pale cool background.
The emblem is centered and scaled to fill the upper-middle region while leaving
clear optical padding around all six edges. The package name sits in a dedicated
lower band. The long curved tagline is removed because it becomes illegible at
README and thumbnail sizes.

## Reproducible implementation

`logo/generate-logo.R` will be the single source of truth. It will construct SVG
markup from fixed numeric geometry, colors, gradients, masks, and filters, then
write the vector masters and render PNG outputs with `rsvg`.

The script will generate:

- `logo/aucmat_emblem.svg` and `logo/aucmat_emblem.png`
- `logo/aucmat_logo.svg` and `logo/aucmat_logo.png`
- `logo/aucmat_hex.svg` and `logo/aucmat_hex.png`
- `man/figures/logo.png`
- `docs/logo.png`
- PNG favicon sizes used by the pkgdown site

All paths are relative to the repository root. The script will check for the
required renderer, create missing output directories, use an SVG-safe generic
sans-serif font stack, and fail with a clear installation message when `rsvg`
is unavailable.

## Validation

The generator must run from a clean R session without interactive input. Each
SVG must parse successfully, every expected output must exist and be non-empty,
and PNG dimensions must match their declared targets. Visual validation will
inspect the emblem, full lockup, hex sticker, and the smallest favicon to confirm
that the ROC curve, matrix, and wordmark remain distinct.

## Scope

This change is limited to logo-generation code and generated logo/favicon
assets. It does not alter package functionality, plotting themes, or the public
R API.
