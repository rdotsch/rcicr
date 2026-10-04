# Plan: less memory per parallel worker (#405)

## What changes

Every `foreach` worker holds its own copy of the noise basis (`patches` 120 MB and `patchIdx` 120 MB at 512 px) plus everything else its loop body names, and each render adds about 250 MB of temporaries. Three changes, none of which moves a number:

1. **`patchIdx` is integer in memory.** A copy of the basis used for rendering gets `storage.mode(patchIdx) <- "integer"`, where the basis is built (`generateStimuli2IFC()`) and where it is loaded (`loadStimulusParams()`, `referenceNoise()`). The basis saved to the `.Rdata` file stays double, so the file contract is unchanged. `params[patchIdx]` selects the same elements with either index type; integer indexing skips the per-element conversion R does for a double index.
2. **The participant loop sends each participant only its own rows.** `computeParticipantCIs()` iterates over per-participant parameter and response subsets, so they travel with each task instead of the whole `params` matrix going to every worker.
3. **Never more workers than tasks.** `startBackend()` takes the number of iterations and starts `min(ncores, iterations)` workers. Four participants on eight cores start four.
`generateCI()`'s `n_cores` is documented as the z-map's core count; it also runs the participant loop. Its `@param`, and the others' `ncores`, say what a worker costs.

## Not changed

- No numeric output. Workers draw no random numbers (`startBackend()`'s comment), so where a CI is computed cannot change it; every change above is checked with `identical()`.
- No `.Rdata` field or type, no argument, no default `ncores`.
- **Not a serial threshold.** Parallel loses with few participants (table below) but wins with many; where they cross depends on the machine, so a fixed threshold would be wrong elsewhere.
- **Not fork workers.** `doSNOW` cannot drive them: `snow`'s `sendData()` has no method for `parallel`'s fork nodes, so a fork cluster fails on its first `clusterCall()`. Through `doParallel` they run, but every free variable of the loop body is still serialized to each worker, so the basis is copied per worker anyway. Only with the basis kept out of the export, and read from a shared environment instead, did memory fall: 4 workers rendering 40 trials, 733 MB per worker and 3.5 GB in total with socket workers, 589 MB and 2.6 GB with fork workers (`fork2.R`). That needs a new dependency, two forms of every loop body, and gives up the `.options.snow` progress callback `doSNOW` was chosen for (#178).
- **Not chunked rendering.** Rendering in blocks of pixels is bit-identical and faster, but did not lower peak memory: R's collector lets the blocks' temporaries accumulate. It belongs in an issue of its own, since it changes the code every render runs.

## Measured

Synthetic stimulus file: 512 px, 5 scales, 1400 trials (`make.R` below). Container with 4 cores and 16 GB, R 4.3.3. Elapsed seconds and peak RSS summed over all R processes (`measure.sh`), mean of two runs. "Prototype" is changes 1 and 2, socket workers.

| individual CIs saved | participants | cores | `main` s | `main` MB | prototype s | prototype MB |
|---|---|---|---|---|---|---|
| no | 4 | 1 | 5.3 | 953 | 4.5 | 919 |
| no | 4 | 4 | 11.6 | 3573 | 8.9 | 3343 |
| no | 20 | 1 | 11.1 | 965 | 9.2 | 1073 |
| no | 20 | 4 | 14.1 | 4614 | 10.8 | 3861 |
| no | 100 | 1 | 42.3 | 1335 | 32.1 | 1516 |
| no | 100 | 4 | 22.8 | 4266 | 18.9 | 4176 |
| yes | 100 | 1 | 47.2 | 1332 | 37.0 | 1502 |
| yes | 100 | 4 | 23.6 | 4622 | 19.0 | 4102 |

The prototype's `generateCI()` result is `identical()` to `main`'s for 4, 20 and 100 participants at 1 and 4 cores, individual CIs saved.

Integer `patchIdx`, 20 renders of the file's own trials: `identical()` output; 7.25 s and 7.22 s with double, 5.27 s and 5.34 s with integer.

Fork against socket workers through `foreach` (`fork2.R`), 4 workers rendering 40 trials, proportional set size (PSS), which splits shared pages between the processes sharing them:

```
psock      : workers 733/733/733/733 MB PSS, total 3465 MB, 40 renders 8.5s, checksum -71.85336191
fork       : workers 883/883/883/883 MB PSS, total 3759 MB, 40 renders 8.0s, checksum -71.85336191
fork-shared: workers 589/589/589/589 MB PSS, total 2583 MB, 40 renders 3.6s, checksum -71.85336191
```

`psock` registers `doSNOW`; `fork` and `fork-shared` register `doParallel`, since `doSNOW` with a fork cluster stops with `no applicable method for 'sendData' applied to an object of class "c('forknode', 'SOCK0node')"`.

The serial prototype peaks higher than `main` at 20 and 100 participants, because the per-participant subsets are a second copy of `params`. The implementation builds them only when workers run.

The other three loops (`generateStimuli2IFC()`, the t-test z-map, `referenceNoise()`) render one image per trial, hundreds per call, so parallel pays off there; they get changes 1 and 3 and are measured once each, `main` against the branch, in the analysis file.

## Tests

- `identical()` results serial against parallel, for `generateCI(participants = )`, `generateStimuli2IFC()` parameters, the t-test z-map and an `"images"` reference.
- The saved `.Rdata` file's `p$patchIdx` is still double.
- `startBackend()` starts `min(ncores, iterations)` workers.
- The golden master unchanged, and the release gate against `main` reports 0 deviations (`tools/compare-release-output.R --ref="$(git rev-parse origin/main)"`).

## Documentation

`analyses/worker-memory.Rmd`, knitted, with the scripts below and the numbers above re-measured on the branch. `NEWS.md` "Performance and dependencies" gets the speed and memory change; `DECISIONS.md` records the rejected fork workers and serial threshold. `@param n_cores`/`ncores` corrected.

## Most likely to fail

The worker cap and the per-participant iteration in `computeParticipantCIs()` interact with its `on.exit()` cluster teardown and the serial progress bar; both paths need the `identical()` test, and the cap needs a run where `ncores` exceeds the participants.

## The scripts

`make.R` writes the stimulus file:

```r
library(rcicr)
out <- commandArgs(TRUE)[1]
n <- 512
b <- file.path(out, "base.png"); dir.create(out, showWarnings = FALSE)
x <- matrix(seq(-1, 1, length.out = n), n, n, byrow = TRUE); y <- -matrix(seq(-1, 1, length.out = n), n, n)
face <- exp(-(x^2 / 0.45 + y^2 / 0.75)); png::writePNG((face - min(face)) / diff(range(face)), b)
invisible(capture.output(generateStimuli2IFC(list(face = b), n_trials = 1400, img_size = n, stimulus_path = out, seed = 1, ncores = 4, save_as_png = FALSE)))
```

`run.R` times one `generateCI()` call; arguments are the file, participants, cores, whether to save individual CIs, and optionally a path to save the result for `identical()`:

```r
suppressMessages(library(rcicr))
a <- commandArgs(TRUE); rdata <- a[1]; npids <- as.integer(a[2]); ncores <- as.integer(a[3]); save <- a[4] == "1"
trials_per <- 1400 %/% npids
stimuli <- rep(seq_len(trials_per), npids)
set.seed(9); responses <- sample(c(-1, 1), length(stimuli), TRUE)
participants <- rep(sprintf("p%03d", seq_len(npids)), each = trials_per)
out <- tempfile(); dir.create(out)
t <- system.time(invisible(capture.output(ci <- generateCI(stimuli, responses, "face", rdata, participants = participants,
  save_individual_cis = save, save_as_png = FALSE, targetpath = out, n_cores = ncores))))
if (!is.na(a[5])) saveRDS(ci, a[5])
cat(sprintf("%.2f %s\n", t[["elapsed"]], digest <- format(sum(ci$ci), digits = 17)))
```

`measure.sh` adds the peak memory:

```sh
#!/bin/bash
# usage: measure.sh rdata npids ncores save -> "elapsed checksum peak_rss_MB"
# Peak RSS summed over every R process, workers included, sampled every 0.1 s.
out=$(mktemp)
Rscript "$(dirname "$0")/run.R" "$@" > "$out" 2>/dev/null &
pid=$!
peak=0
while kill -0 $pid 2>/dev/null; do
  s=$(ps -C R -C Rscript -o rss= 2>/dev/null | awk '{t+=$1} END {print int(t/1024)}')
  [ "${s:-0}" -gt "$peak" ] && peak=$s
  sleep 0.1
done
echo "$(grep -v '^starting worker' "$out") $peak"
rm -f "$out"
```

The table is `measure.sh <file> <participants> <cores> <saved>` over participants 4, 20, 100, cores 1, 2, 4 and saved 0, 1, twice each, with `R_LIBS` naming a library holding `main` or the prototype.

`fork2.R` compares the cluster types through `foreach` (`Rscript fork2.R <dir> psock|fork|fork-shared`):

```r
suppressMessages({library(rcicr); library(foreach)}); a <- commandArgs(TRUE); mode <- a[2]
e <- new.env(); load(list.files(a[1], "Rdata$", full.names = TRUE), envir = e)
pss <- function(pids) sum(vapply(pids, function(pid) { l <- readLines(sprintf("/proc/%d/smaps_rollup", pid)); as.numeric(sub("^Pss:\\s+(\\d+).*", "\\1", grep("^Pss:", l, value = TRUE))) / 1024 }, 1))
shared <- new.env()
run <- function() {
  p <- e$p; storage.mode(p$patchIdx) <- "integer"; par <- e$stimuli_params$face[1:40, ]
  if (mode == "fork-shared") assign("p", p, envir = shared)
  cl <- if (mode == "psock") parallel::makeCluster(4) else parallel::makeForkCluster(4)
  if (mode == "psock") doSNOW::registerDoSNOW(cl) else doParallel::registerDoParallel(cl); wp <- unlist(parallel::clusterCall(cl, Sys.getpid))
  peak <- 0
  t <- system.time(if (mode == "fork-shared") {
    r <- foreach(i = 1:40, .combine = c, .packages = "rcicr", .noexport = "p", .export = character(0)) %dopar% { x <- generateNoiseImage(par[i, ], get("p", envir = shared)); peak <<- 0; sum(x) }
  } else {
    r <- foreach(i = 1:40, .combine = c, .packages = "rcicr") %dopar% sum(generateNoiseImage(par[i, ], p))
  })[["elapsed"]]
  w <- vapply(wp, pss, 1); tot <- pss(c(Sys.getpid(), wp))
  parallel::stopCluster(cl)
  cat(sprintf("%-11s: workers %s MB PSS, total %.0f MB, 40 renders %.1fs, checksum %.10g\n", mode, paste(round(w), collapse = "/"), tot, t, sum(r)))
}
run()
```

`idx.R` times rendering with either index type (`Rscript idx.R <dir>`):

```r
suppressMessages(library(rcicr)); e <- new.env(); load(list.files(commandArgs(TRUE)[1], "Rdata$", full.names = TRUE), envir = e)
p <- e$p; pi <- p; storage.mode(pi$patchIdx) <- "integer"
cat("patchIdx MB, double vs integer:", format(object.size(p$patchIdx), units = "MB"), "/", format(object.size(pi$patchIdx), units = "MB"), "\n")
par <- e$stimuli_params$face
same <- all(vapply(1:20, function(i) identical(generateNoiseImage(par[i, ], p), generateNoiseImage(par[i, ], pi)), TRUE))
cat("20 trials identical:", same, "\n")
for (r in 1:2) { td <- system.time(for (i in 1:20) generateNoiseImage(par[i, ], p))[["elapsed"]]; ti <- system.time(for (i in 1:20) generateNoiseImage(par[i, ], pi))[["elapsed"]]; cat(sprintf("20 renders: double %.2fs, integer %.2fs\n", td, ti)) }
```

```
patchIdx MB, double vs integer: 120 Mb / 60 Mb
20 trials identical: TRUE
20 renders: double 7.25s, integer 5.27s
20 renders: double 7.22s, integer 5.34s
```

The prototype is `main` with `computeParticipantCIs()` iterating `foreach(obs = seq_len(npids), obs_params = pid_params, obs_responses = pid_responses, .noexport = c('pid_params', 'pid_responses'))` over `lapply(seq_len(npids), function(obs) params[pids == obs, ])` and the matching responses, calling `generateCINoise(obs_params, obs_responses, p)` with `p$patchIdx` made integer first.
