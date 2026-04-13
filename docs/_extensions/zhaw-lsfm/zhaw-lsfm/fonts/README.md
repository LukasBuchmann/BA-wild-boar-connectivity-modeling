# Fonts Directory

This directory can contain local font files (.ttf or .otf) that will be used by Typst in addition to system fonts.

## Font Configuration

The template does not set a default font — Typst will use its built-in default. ZHAW recommends Arial 11pt, but does not require it. To use Arial or another font, set `mainfont` in your `_quarto.yml`:

```yaml
format:
  zhaw-lsfm-typst:
    mainfont: "Arial"
```

If the font is not installed on your system, place the font files (.ttf or .otf) in this directory and Typst will find them automatically.
