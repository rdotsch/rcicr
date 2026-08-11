# Regenerates man/figures/logo.png. Run from the package root:
#
#   Rscript tools/make-logo.R
#
# The logo is a real classification image, computed by this package from a
# simulated 2IFC experiment -- the same construction the walkthrough vignette
# uses, so the picture on the tin is the thing in the tin.
#
# The observer is simulated with low internal noise so the CI recovers its
# template cleanly. A CI from noisy responses is a grey smudge, which is true to
# life but illegible at favicon size.
#
# The base face is synthetic, never a photograph (CONTRIBUTING.md, "Base images in
# tests and vignettes are always synthetic").

suppressMessages(devtools::load_all(".", quiet = TRUE))

size <- 256
trials <- 400

# ncores = 1 is required, not a preference. generateStimuli2IFC()'s workers
# library(rcicr) themselves, which load_all() does not satisfy: on a fresh
# checkout they fail with "there is no package called 'rcicr'", and where an
# older rcicr happens to be installed they quietly build the logo from *that*
# instead of the working tree. Running in-process is also the only way this
# script can be trusted to reflect the code beside it. See AGENTS.md, "Common
# commands".
ncores <- 1

# --- synthetic base face ------------------------------------------------------

ax <- seq(-1, 1, length.out = size)
x <- matrix(ax, size, size, byrow = TRUE)
y <- -matrix(ax, size, size)

face <- exp(-(x^2 / 0.45 + y^2 / 0.75))
face <- face - 0.55 * exp(-(((x + 0.3)^2 + (y - 0.25)^2) / 0.012))
face <- face - 0.55 * exp(-(((x - 0.3)^2 + (y - 0.25)^2) / 0.012))
face <- face - 0.35 * exp(-(x^2 / 0.10 + (y + 0.42)^2 / 0.006))
face <- (face - min(face)) / (max(face) - min(face))

work <- tempfile("logo")
dir.create(work)
base_face <- file.path(work, "face.png")
png::writePNG(face, base_face)

# --- run the pipeline for real ------------------------------------------------

generateStimuli2IFC(
  base_face_files = list(face = base_face),
  n_trials = trials,
  img_size = size,
  stimulus_path = work,
  seed = 1,
  nscales = 5,
  ncores = ncores,
  save_as_png = FALSE
)

rdata <- list.files(work, pattern = "\\.Rdata$", full.names = TRUE)[1]

e <- new.env()
load(rdata, envir = e)
params <- e$stimuli_params[["face"]]

set.seed(99)
template <- generateNoiseImage(rnorm(max(e$p$patchIdx)), e$p)

stack <- vapply(
  seq_len(nrow(params)),
  function(i) generateNoiseImage(params[i, ], e$p),
  matrix(0, size, size)
)
evidence <- apply(stack, 3, function(z) base::sum(z * template))
evidence <- evidence / sd(evidence)

set.seed(7)
responses <- ifelse(evidence + rnorm(length(evidence), 0, 0.4) > 0, 1, -1)

ci <- generateCI(
  stimuli = seq_len(trials),
  responses = responses,
  baseimage = "face",
  rdata = rdata,
  save_as_png = FALSE
)

# combined is what save_as_png writes: the scaled CI overlaid on the base image,
# which is the output researchers actually look at.
img <- ci$combined

# Stretch to the 1st-99th percentile. Straight min/max normalisation leaves the
# logo a flat mid-grey -- combined averages the CI against the base image, so its
# values bunch in the middle and the structure is invisible below about 150px.
q <- stats::quantile(img, c(0.01, 0.99))
img <- pmin(pmax((img - q[1]) / (q[2] - q[1]), 0), 1)

# --- draw at 4x with the package name -----------------------------------------

out_size <- 1024

tmp <- tempfile(fileext = ".png")
grDevices::png(tmp, width = out_size, height = out_size, bg = "transparent", type = "cairo")
graphics::par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
graphics::rasterImage(img, 0, 0, 1, 1, interpolate = TRUE)
# A pointy-top hexagon is only about 47% as wide at this height as at its
# middle, so the name has to stay well under half the image width or it clips
# through the edge. Drawn twice: a dark offset copy first, so white letters stay
# legible over the light patches of the CI.
for (d in list(c(0.006, -0.006), c(0, 0))) {
  graphics::text(0.5 + d[1], 0.15 + d[2], "rcicr",
    col = if (all(d == 0)) "#FFFFFF" else "#2A0033",
    cex = 7, font = 2, family = "sans"
  )
}
grDevices::dev.off()

img <- png::readPNG(tmp)[, , 1:3]
size <- out_size

# --- hexagon mask -------------------------------------------------------------

# Same idea as applyMask(): decide per pixel whether it is inside the shape and
# zero the alpha outside. Pointy-top, the orientation R hex stickers use.
px <- seq(-1, 1, length.out = size)
gx <- matrix(px, size, size, byrow = TRUE)
gy <- matrix(px, size, size)

hex <- function(r) {
  q <- sqrt(3) / 2
  abs(gx) <= r * q & (q * abs(gy) + 0.5 * abs(gx)) <= r * q
}

inside <- hex(1.0)
border <- inside & !hex(0.93)

rgba <- array(0, dim = c(size, size, 4))
rgba[, , 1:3] <- img
rgba[, , 4] <- 1

edge <- c(0.27, 0.00, 0.33) # viridis dark end, matching the z-map palette
for (k in 1:3) {
  channel <- rgba[, , k]
  channel[border] <- edge[k]
  rgba[, , k] <- channel
}
alpha <- rgba[, , 4]
alpha[!inside] <- 0
rgba[, , 4] <- alpha

dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)
png::writePNG(rgba, "man/figures/logo.png")

cat("wrote man/figures/logo.png (", size, "x", size, ")\n")

# --- favicons -----------------------------------------------------------------

# Regenerated here rather than by hand: they are derived from the logo, so a new
# logo with the old favicons beside it is a silent mismatch.
pkgdown::build_favicons(overwrite = TRUE)

# build_favicons() writes root-absolute icon paths ("/web-app-manifest-...png"),
# which assume the site is served from the domain root. This one is a project
# page under /rcicr/, so those resolve to https://rdotsch.github.io/... and 404
# for anyone installing the site as an app. Relative paths work under any
# prefix. The generator also leaves the names empty, which shows as a blank
# label on the installed icon.
manifest_path <- "pkgdown/favicon/site.webmanifest"
manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
manifest$name <- "rcicr"
manifest$short_name <- "rcicr"
manifest$icons <- lapply(manifest$icons, function(icon) {
  icon$src <- sub("^/", "", icon$src)
  icon
})
jsonlite::write_json(manifest, manifest_path,
  auto_unbox = TRUE, pretty = TRUE
)

cat("wrote pkgdown/favicon/ (icon paths made relative)\n")
