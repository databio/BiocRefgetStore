#!/usr/bin/env Rscript
# Render the package vignettes to a standalone HTML site for GitHub Pages.
#
# This is intentionally NOT pkgdown: it just knits each vignette to HTML with
# rmarkdown and writes a small landing page. The benchmark vignette executes
# its real harness (run_benchmarks) at knit time, so the deployed page shows
# live timing tables + the ggplot, not pre-baked numbers. That requires the
# package + gtars + the Bioconductor Suggests (Biostrings/BSgenome/rtracklayer/
# Rsamtools) + bench + ggplot2 to be installed in the CI job before this runs.

out_dir <- "site"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Build the package itself first so system.file() resolves the bundled FASTA
# and the benchmark harness is loadable.
message("Rendering vignettes to ", normalizePath(out_dir))

# (vignette file, output html, human title, one-line description)
vignettes <- list(
  list(src = "vignettes/benchmark.Rmd",
       html = "benchmark.html",
       title = "Benchmark: BiocRefgetStore vs BSgenome/Biostrings",
       desc = "Live, executed benchmark comparing RefgetStore-backed sequence extraction against the classic Bioconductor stack. Real timing tables and a comparison plot, produced at build time."),
  list(src = "vignettes/getting-started.Rmd",
       html = "getting-started.html",
       title = "Getting Started",
       desc = "Connect to a refget store and access sequences with a BSgenome-compatible interface."),
  list(src = "vignettes/reference.Rmd",
       html = "reference.html",
       title = "Reference",
       desc = "Complete reference for every exported constructor, method, and class.")
)

rendered <- character(0)
for (v in vignettes) {
  message("== rendering ", v$src)
  ok <- tryCatch({
    rmarkdown::render(
      input = v$src,
      output_format = "html_document",
      output_file = v$html,
      output_dir = out_dir,
      output_options = list(
        toc = TRUE, toc_float = TRUE, self_contained = TRUE,
        theme = "cosmo"
      ),
      quiet = FALSE,
      envir = new.env()
    )
    TRUE
  }, error = function(e) {
    message("!! FAILED to render ", v$src, ": ", conditionMessage(e))
    FALSE
  })
  if (ok) rendered <- c(rendered, v$html)
}

# The benchmark vignette is the headline; fail the build if it didn't render,
# so a broken benchmark page never silently deploys.
if (!"benchmark.html" %in% rendered) {
  stop("benchmark.html did not render; aborting so a broken page is not deployed.")
}

# Landing page. Benchmark first and called out, since that is the point of the site.
links <- vapply(vignettes, function(v) {
  status <- if (v$html %in% rendered) {
    sprintf('<a href="%s">%s</a>', v$html, v$title)
  } else {
    sprintf('<span style="color:#999">%s (failed to render)</span>', v$title)
  }
  sprintf("<li><strong>%s</strong><br/><span class=\"desc\">%s</span></li>",
          status, v$desc)
}, character(1))

index <- sprintf('<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>BiocRefgetStore</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
         max-width: 760px; margin: 3rem auto; padding: 0 1rem; line-height: 1.55; color: #1a1a1a; }
  h1 { margin-bottom: 0.2rem; }
  .sub { color: #555; margin-top: 0; }
  ul { list-style: none; padding: 0; }
  li { margin: 1.1rem 0; padding: 1rem 1.2rem; border: 1px solid #e2e2e2; border-radius: 8px; }
  li:first-of-type { border-color: #2b6cb0; background: #f5f9ff; }
  a { color: #2b6cb0; text-decoration: none; font-size: 1.1rem; }
  a:hover { text-decoration: underline; }
  .desc { color: #555; font-size: 0.95rem; }
  footer { margin-top: 2.5rem; color: #888; font-size: 0.85rem; }
</style>
</head>
<body>
<h1>BiocRefgetStore</h1>
<p class="sub">A BSgenome-compatible interface to sequences backed by a refget store (Rust-backed <code>gtars</code>).</p>
<p>The <strong>benchmark vignette</strong> below runs the package&#39;s real benchmark
harness at build time, so its timing tables and comparison plot reflect an actual
execution on the CI runner &mdash; not pre-baked numbers.</p>
<ul>
%s
</ul>
<footer>
Source: <a href="https://github.com/databio/BiocRefgetStore">github.com/databio/BiocRefgetStore</a>.
Pages rendered from the package vignettes by a GitHub Actions workflow.
</footer>
</body>
</html>
', paste(links, collapse = "\n"))

writeLines(index, file.path(out_dir, "index.html"))
message("Wrote ", file.path(out_dir, "index.html"))
message("Rendered: ", paste(rendered, collapse = ", "))
