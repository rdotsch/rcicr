# Run from the repository root with this checkout installed: Rscript tools/audit-runtime.R REPORT.md
# Diagnostic observations are separate from regression tests of intended behavior.
library(rcicr)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
report_path <- args[[1]]
report <- c('# Audit runtime measurements', '',
  paste('Commit:', Sys.getenv('GITHUB_SHA', unset = 'local checkout')), '',
  paste('Package:', as.character(packageVersion('rcicr'))), '',
  'These are synthetic reproductions, not estimates of affected-study prevalence.', '')
failures <- character()

observe <- function(expr) {
  warnings <- character()
  error <- NULL
  output <- capture.output(value <- tryCatch(
    withCallingHandlers(expr, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart('muffleWarning')
    }), error = function(e) {
      error <<- conditionMessage(e)
      NULL
    }
  ))
  list(value = value, warnings = warnings, error = error, output = output)
}

checked <- function(expr) {
  result <- observe(expr)
  if (!is.null(result$error)) stop(result$error)
  result$value
}

metricText <- function(x) {
  if (!length(x)) return('(none)')
  if (is.numeric(x)) x <- format(x, digits = 17)
  gsub('[\r\n\t]+', ' ', paste(x, collapse = ', '))
}

runCase <- function(issue, fn) {
  cat('\nAUDIT_START\t', issue, '\n', sep = '')
  result <- observe(fn())
  if (!is.null(result$error)) {
    failures <<- c(failures, issue)
    values <- list(status = 'DIAGNOSTIC ERROR', error = result$error)
  } else {
    values <- result$value
  }
  values$uncaught_warnings <- result$warnings
  lines <- vapply(names(values), function(key) {
    value <- metricText(values[[key]])
    cat('AUDIT\t', issue, '\t', key, '\t', value, '\n', sep = '')
    paste0('- **', key, ':** ', value)
  }, character(1))
  report <<- c(report, paste('##', issue), '', lines, '')
  writeLines(report, report_path)
}

readFixture <- function(path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  as.list(env)
}

makeFixture <- function(dir, bases = 1L, n_trials = 24L, nscales = 2L,
                        shared = TRUE, uniform = FALSE, ncores = 1L,
                        return_frame = FALSE, save_png = FALSE) {
  dir.create(dir, recursive = TRUE)
  image_paths <- setNames(file.path(dir, paste0('base', seq_len(bases), '.png')),
    paste0('base', seq_len(bases)))
  for (i in seq_len(bases)) {
    img <- if (uniform) matrix(0.5, 32, 32) else {
      matrix(((seq_len(1024) * (2 * i + 1)) %% 1024) / 1023, 32, 32)
    }
    png::writePNG(img, image_paths[[i]])
  }
  generated <- checked(generateStimuli2IFC(
    as.list(image_paths), n_trials = n_trials, img_size = 32,
    stimulus_path = file.path(dir, 'stimuli'), seed = 31,
    use_same_parameters = shared, maximize_baseimage_contrast = !uniform,
    nscales = nscales, ncores = ncores, return_as_dataframe = return_frame,
    save_as_png = save_png, save_rdata = TRUE
  ))
  path <- list.files(file.path(dir, 'stimuli'), pattern = '\\.Rdata$', full.names = TRUE)
  stopifnot(length(path) == 1L)
  list(path = path, saved = readFixture(path), frame = generated,
    images = image_paths, dir = dir)
}

noiseFor <- function(fixture, base) {
  params <- fixture$saved$stimuli_params[[base]]
  stopifnot(is.matrix(params), nrow(params) == fixture$saved$n_trials)
  vapply(seq_len(nrow(params)), function(i) {
    as.vector(generateNoiseImage(params[i, ], fixture$saved$p))
  }, numeric(fixture$saved$img_size^2))
}

normsFor <- function(noise, draws) {
  # A column is one complete response draw; the noise selector is supplied independently.
  vapply(seq_len(ncol(draws)), function(i) {
    ci <- noise %*% draws[, i] / ncol(noise)
    sqrt(sum(ci^2))
  }, numeric(1))
}

maxDiff <- function(a, b) max(abs(a - b))
near <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-12, check.attributes = FALSE))

main <- function() {
  root <- tempfile('rcicr-audit-')
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  runCase('#299', function() {
    f <- makeFixture(file.path(root, 'independent'), bases = 2L, shared = FALSE)
    first <- noiseFor(f, 'base1')
    second <- noiseFor(f, 'base2')
    stopifnot(maxDiff(first, second) > 1e-8)
    iter <- 512L
    set.seed(97)
    draws <- matrix(2 * (runif(f$saved$n_trials * iter) > 0.5) - 1,
      nrow = f$saved$n_trials)
    oracle1 <- normsFor(first, draws)
    oracle2 <- normsFor(second, draws)
    stopifnot(stats::mad(oracle1) > 0, stats::mad(oracle2) > 0)
    actual <- checked(generateReferenceDistribution2IFC(f$path, iter = iter,
      ncores = 1, response_seed = 97, save_rdata = TRUE))
    parallel <- checked(generateReferenceDistribution2IFC(f$path, iter = iter,
      ncores = 2, response_seed = 97, save_rdata = FALSE))
    stopifnot(identical(actual, parallel), near(actual, oracle1))

    # Seeded synthetic choices favor one nonzero pixel of the second base's noise.
    pixel <- which.max(rowSums(second^2))
    responses <- ifelse(second[pixel, ] >= 0, 1, -1)
    target <- checked(generateCI(seq_len(f$saved$n_trials), responses,
      'base2', f$path, save_as_png = FALSE, scaling = 'none', n_cores = 1))
    direct_ci <- as.vector(second %*% responses / ncol(second))
    stopifnot(near(as.vector(target$ci), direct_ci))
    observed_z <- checked(computeInfoVal2IFC(target, f$path))
    target_norm <- sqrt(sum(target$ci^2))
    correct_z <- (target_norm - stats::median(oracle2)) / stats::mad(oracle2)
    first_z <- (target_norm - stats::median(oracle1)) / stats::mad(oracle1)
    stopifnot(near(observed_z, first_z))

    shared <- makeFixture(file.path(root, 'shared'), bases = 2L)
    shared_noise <- noiseFor(shared, 'base2')
    shared_ref <- checked(generateReferenceDistribution2IFC(shared$path,
      iter = iter, ncores = 1, response_seed = 97, save_rdata = FALSE))
    stopifnot(identical(noiseFor(shared, 'base1'), shared_noise),
      near(shared_ref, normsFor(shared_noise, draws)))

    rebuilt <- checked(generateStimuli2IFC(f$saved$base_face_files,
      n_trials = f$saved$n_trials, img_size = 32, seed = f$saved$seed,
      nscales = f$saved$nscales, ncores = 1, return_as_dataframe = TRUE,
      save_as_png = FALSE, save_rdata = FALSE))
    default_draws <- matrix(2 * (runif(f$saved$n_trials * iter) > 0.5) - 1,
      nrow = f$saved$n_trials)
    default_oracle <- normsFor(first, default_draws)
    default_ref <- checked(generateReferenceDistribution2IFC(f$path, iter = iter,
      ncores = 1, save_rdata = FALSE))
    default_parallel <- checked(generateReferenceDistribution2IFC(f$path, iter = iter,
      ncores = 2, save_rdata = FALSE))
    stopifnot(near(as.matrix(rebuilt), first), near(default_ref, default_oracle),
      identical(default_ref, default_parallel))
    list(status = if (near(actual, oracle1) && !near(actual, oracle2)) 'CONFIRMED' else 'NOT REPRODUCED',
      stimulus_seed = 31, response_seed = 97, trials = f$saved$n_trials, iterations = iter,
      first_vs_second_noise_max = maxDiff(first, second),
      actual_vs_first_reference_max = maxDiff(actual, oracle1),
      actual_vs_second_reference_max = maxDiff(actual, oracle2),
      actual_reference_median = stats::median(actual), actual_reference_mad = stats::mad(actual),
      second_reference_median = stats::median(oracle2), second_reference_mad = stats::mad(oracle2),
      second_base_reported_infoval = observed_z, second_base_oracle_infoval = correct_z,
      infoval_difference = observed_z - correct_z,
      serial_parallel_identical = identical(actual, parallel),
      shared_parameter_control = near(shared_ref, normsFor(shared_noise, draws)),
      default_stream_control = near(default_ref, default_oracle),
      default_serial_parallel_identical = identical(default_ref, default_parallel))
  })

  runCase('#301', function() {
    moved <- makeFixture(file.path(root, 'moved'))
    baseline <- checked(generateReferenceDistribution2IFC(moved$path,
      iter = 8, ncores = 1, save_rdata = FALSE))
    stopifnot(all(is.finite(baseline)), all(is.finite(noiseFor(moved, 'base1'))))
    unlink(moved$images)
    stopifnot(!any(file.exists(moved$images)), file.exists(moved$path))
    missing <- observe(generateReferenceDistribution2IFC(moved$path,
      iter = 8, ncores = 1, save_rdata = FALSE))
    grey <- makeFixture(file.path(root, 'grey'), uniform = TRUE)
    stopifnot(all(is.finite(noiseFor(grey, 'base1'))), all(file.exists(grey$images)))
    uniform <- observe(generateReferenceDistribution2IFC(grey$path,
      iter = 8, ncores = 1, save_rdata = FALSE))
    confirmed <- !is.null(missing$error) && grepl('could not be read', missing$error) &&
      !is.null(uniform$error) && grepl('has no contrast', uniform$error)
    list(status = if (confirmed) 'CONFIRMED' else 'NOT REPRODUCED',
      original_path_control_finite = all(is.finite(baseline)),
      saved_noise_finite_in_both_cases = TRUE,
      missing_source_error = missing$error, uniform_source_error = uniform$error)
  })

  runCase('#300', function() {
    f <- makeFixture(file.path(root, 'indices'), n_trials = 6L, nscales = 1L)
    stopifnot(ncol(f$saved$stimuli_params$base1) == 12L)
    ciFor <- function(ids, responses) generateCI(ids, responses, 'base1', f$path,
      save_as_png = FALSE, scaling = 'none', n_cores = 1)
    integer <- checked(ciFor(c(1, 2), c(1, -1)))
    fractional <- observe(ciFor(c(1.2, 2.8), c(1, -1)))
    zero <- observe(ciFor(c(0, 1, 2), c(1, -1, 1)))
    aligned <- checked(ciFor(c(1, 2), c(-1, 1)))
    params <- f$saved$stimuli_params$base1[c(0, 1, 2), ]
    stopifnot(nrow(params) == 2L, length(params) %% 3L == 0L)
    recycled <- generateCINoise(params, c(1, -1, 1), f$saved$p)
    fraction_equal <- is.null(fractional$error) && near(fractional$value$ci, integer$ci)
    zero_misaligned <- is.null(zero$error) && near(zero$value$ci, recycled) &&
      !near(zero$value$ci, aligned$ci)
    list(status = if (fraction_equal && zero_misaligned) 'CONFIRMED' else 'NOT REPRODUCED',
      fractional_equals_truncated_ids = fraction_equal, fractional_warnings = fractional$warnings,
      fractional_error = fractional$error, zero_warnings = zero$warnings, zero_error = zero$error,
      zero_selected_rows = nrow(params), responses_supplied = 3L,
      zero_matches_recycled_multiplication = is.null(zero$error) && near(zero$value$ci, recycled),
      zero_vs_aligned_ci_max = if (is.null(zero$error)) maxDiff(zero$value$ci, aligned$ci) else NA_real_)
  })

  runCase('#294', function() {
    f <- makeFixture(file.path(root, 'participants'), n_trials = 20L)
    responses <- rep(c(1, 1, -1, 1, -1), 4)
    params <- f$saved$stimuli_params$base1
    desired_ids <- c(rep('a', 4), rep('b', 16))
    short_ids <- c('a', 'b')
    recycled_ids <- rep(short_ids, 10)
    perPerson <- function(ids) lapply(c('a', 'b'), function(id) {
      rows <- ids == id
      generateCINoise(params[rows, ], responses[rows], f$saved$p)
    })
    intended_individual <- perPerson(desired_ids)
    recycled_individual <- perPerson(recycled_ids)
    intended_mean <- (intended_individual[[1]] + intended_individual[[2]]) / 2
    stopifnot(maxDiff(intended_mean,
      (recycled_individual[[1]] + recycled_individual[[2]]) / 2) > 1e-8)
    observations <- lapply(c(1L, 2L), function(cores) {
      out <- file.path(f$dir, paste0('cis-', cores))
      short <- observe(generateCI(seq_len(20), responses, 'base1', f$path,
        participants = short_ids, save_individual_cis = TRUE, targetpath = out,
        individual_scaling = 'constant', individual_scaling_constant = 1,
        save_as_png = FALSE, scaling = 'none', n_cores = cores))
      explicit <- checked(generateCI(seq_len(20), responses, 'base1', f$path,
        participants = recycled_ids, save_as_png = FALSE, scaling = 'none', n_cores = cores))
      intended <- checked(generateCI(seq_len(20), responses, 'base1', f$path,
        participants = desired_ids, save_as_png = FALSE, scaling = 'none', n_cores = cores))
      stopifnot(near(intended$ci, intended_mean))
      individual_match <- FALSE
      if (is.null(short$error)) {
        actual_png <- png::readPNG(file.path(out, 'individual_cis', 'ci_a.png'))
        expected_png_path <- file.path(out, 'expected-a.png')
        png::writePNG(((recycled_individual[[1]] + 1) / 2 + f$saved$base_faces$base1) / 2,
          expected_png_path)
        individual_match <- identical(actual_png, png::readPNG(expected_png_path))
      }
      list(short = short, explicit = explicit, intended = intended,
        individual_matches_recycled = individual_match)
    })
    serial <- observations[[1]]
    parallel <- observations[[2]]
    confirmed <- all(vapply(observations, function(x) {
      is.null(x$short$error) && near(x$short$value$ci, x$explicit$ci) &&
        !near(x$short$value$ci, x$intended$ci) && x$individual_matches_recycled
    }, logical(1)))
    list(status = if (confirmed) 'CONFIRMED' else 'NOT REPRODUCED',
      supplied_participants = 2L, trials = 20L, intended_group_sizes = c(4L, 16L),
      serial_error = serial$short$error, parallel_error = parallel$short$error,
      serial_warnings = serial$short$warnings, parallel_warnings = parallel$short$warnings,
      group_vs_intended_max = if (is.null(serial$short$error)) {
        maxDiff(serial$short$value$ci, serial$intended$ci)
      } else NA_real_,
      individual_png_matches_recycled = serial$individual_matches_recycled,
      serial_parallel_identical = is.null(serial$short$error) && is.null(parallel$short$error) &&
        identical(serial$short$value$ci, parallel$short$value$ci))
  })

  runCase('#302', function() {
    values <- list()
    confirmed <- logical()
    for (shared in c(TRUE, FALSE)) {
      for (cores in c(1L, 2L)) {
        key <- paste0(if (shared) 'shared' else 'independent', '_cores', cores)
        f <- makeFixture(file.path(root, key), bases = 2L, n_trials = 2L,
          shared = shared, ncores = cores, return_frame = TRUE, save_png = TRUE)
        files <- list.files(file.path(f$dir, 'stimuli'), pattern = '\\.png$')
        stopifnot(identical(dim(f$frame), c(1024L, 2L)), near(as.matrix(f$frame), noiseFor(f, 'base1')))
        expected <- unlist(lapply(c('base1', 'base2'), function(base) {
          paste0('rcic_', base, '_31_', rep(sprintf('%05d', seq_len(2)), each = 2),
            rep(c('_ori.png', '_inv.png'), 2))
        }))
        missing <- setdiff(expected, files)
        confirmed <- c(confirmed, length(files) == 4L && length(missing) == 4L &&
          all(grepl('_base2_', missing)))
        values[[paste0(key, '_png_count')]] <- length(files)
        values[[paste0(key, '_missing')]] <- missing
        values[[paste0(key, '_first_base_return_max_difference')]] <-
          maxDiff(as.matrix(f$frame), noiseFor(f, 'base1'))
      }
    }
    c(list(status = if (all(confirmed)) 'CONFIRMED' else 'NOT REPRODUCED',
      expected_png_count_per_call = 8L), values)
  })

  runCase('#303', function() {
    f <- makeFixture(file.path(root, 'zero'), n_trials = 2L)
    zero <- checked(generateCI(c(1, 1, 2, 2), c(1, -1, 1, -1), 'base1', f$path,
      save_as_png = FALSE, n_cores = 1))
    stopifnot(all(zero$ci == 0), all(is.finite(zero$ci)))
    mask <- matrix(1, 32, 32)
    mask[seq_len(16), ] <- 0
    masked <- checked(generateCI(c(1, 1, 2, 2), c(1, -1, 1, -1), 'base1', f$path,
      save_as_png = FALSE, n_cores = 1, mask = mask))
    stopifnot(all(is.na(masked$ci[mask == 0])), all(masked$ci[mask == 1] == 0))
    all_zero <- checked(autoscale(list(zero = zero, masked = masked), save_as_pngs = FALSE))
    signal <- checked(generateCI(c(1, 2), c(1, -1), 'base1', f$path,
      save_as_png = FALSE, n_cores = 1))
    stopifnot(max(abs(signal$ci)) > 0)
    mixed <- checked(autoscale(list(zero = zero, signal = signal), save_as_pngs = FALSE))
    stopifnot(all(mixed$zero$scaled == 0.5), all(is.finite(mixed$signal$scaled)),
      identical(mixed$zero$combined, zero$combined))
    matched <- checked(generateCI(c(1, 1, 2, 2), c(1, -1, 1, -1), 'base1', f$path,
      save_as_png = FALSE, n_cores = 1, scaling = 'matched'))
    confirmed <- all(is.nan(zero$scaled)) && all(is.nan(zero$combined)) &&
      all(is.nan(all_zero$zero$scaled)) && all(is.nan(masked$scaled[mask == 1]))
    list(status = if (confirmed) 'CONFIRMED' else 'NOT REPRODUCED',
      raw_zero_pixels = sum(zero$ci == 0), scaled_nan_pixels = sum(is.nan(zero$scaled)),
      combined_nan_pixels = sum(is.nan(zero$combined)),
      mask_intended_na_pixels = sum(is.na(masked$ci)),
      masked_unmasked_nan_pixels = sum(is.nan(masked$scaled[mask == 1])),
      autoscale_all_zero_nan_pixels = sum(is.nan(all_zero$zero$scaled)),
      autoscale_masked_unmasked_nan_pixels = sum(is.nan(all_zero$masked$scaled[mask == 1])),
      mixed_zero_scaled_value = unique(as.vector(mixed$zero$scaled)),
      autoscale_combined_unchanged = identical(mixed$zero$combined, zero$combined),
      matched_nan_pixels = sum(is.nan(matched$scaled)))
  })
}

main()
report <- c(report, '## Session', '', '```', capture.output(sessionInfo()), '```')
writeLines(report, report_path)
cat('\nAUDIT_REPORT_BEGIN\n', paste(report, collapse = '\n'), '\nAUDIT_REPORT_END\n', sep = '')
if (length(failures)) stop('Diagnostic errors: ', paste(failures, collapse = ', '))
