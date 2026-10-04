What a parallel worker costs
================

- [Two installed versions](#two-installed-versions)
- [The measurement](#the-measurement)
- [The stimulus files](#the-stimulus-files)
- [One render, double and integer
  indices](#one-render-double-and-integer-indices)
- [The participant loop](#the-participant-loop)
- [The other loops](#the-other-loops)
- [What a worker costs](#what-a-worker-costs)
- [What this shows, and what it does
  not](#what-this-shows-and-what-it-does-not)

Every parallel loop in rcicr starts worker processes, and each holds its
own copy of the noise basis plus the working memory of one render. This
document measures what that costs at 512 pixels, on `main` and on the
change for [\#405](https://github.com/rdotsch/rcicr/issues/405):

1.  renders use the basis with its patch indices as integer, while the
    `.Rdata` file keeps the basis as built;
2.  in parallel, each participant’s rows travel with that participant’s
    task;
3.  no more workers start than the loop has tasks.

Each configuration runs once, in a fresh R process started in a session
of its own. Memory is the peak resident set size summed over that
session’s processes, workers included, sampled every 0.1 s; it counts
pages shared between processes once per process. The machine and R
version are printed below; absolute times and sizes are this machine’s.

## Two installed versions

``` r
repo <- normalizePath("..")
main_commit <- system2("git", c("-C", shQuote(repo), "rev-parse", "origin/main"), stdout = TRUE)
branch_commit <- system2("git", c("-C", shQuote(repo), "rev-parse", "HEAD"), stdout = TRUE)
install <- function(commit) {
  source_dir <- tempfile("src")
  lib <- tempfile("lib")
  dir.create(source_dir)
  dir.create(lib)
  system2("sh", c("-c", shQuote(sprintf("git -C %s archive %s | tar -x -C %s", shQuote(repo), commit,
                                        shQuote(source_dir)))))
  status <- system2(file.path(R.home("bin"), "R"),
                    c("CMD", "INSTALL", "--no-test-load", "-l", shQuote(lib), shQuote(source_dir)),
                    stdout = FALSE, stderr = FALSE)
  stopifnot(status == 0)
  lib
}
libs <- c(main = install(main_commit), branch = install(branch_commit))
c(main = main_commit, branch = branch_commit)
                                      main                                     branch 
"d4c8adcdfaba5d0e52e7f4ba3e8913a0da6b876d" "7c8878a2453ea48fe591b051ec089893c2c34df7" 
c(R = R.version.string, cores = parallel::detectCores(),
  memory_GB = round(as.numeric(sub("\\D+(\\d+).*", "\\1",
                                   grep("MemTotal", readLines("/proc/meminfo"), value = TRUE))) / 2^20))
                             R                          cores                      memory_GB 
"R version 4.3.3 (2024-02-29)"                            "4"                           "16" 
```

## The measurement

`job.R` runs one call and saves its result, so every result can be
compared with `identical()`. `measure.sh` samples memory while it runs.

``` r
job <- tempfile(fileext = ".R")
writeLines(c(
  'suppressMessages(library(rcicr))',
  'a <- commandArgs(TRUE); what <- a[1]; rdata <- a[2]; cores <- as.integer(a[3]); result <- a[4]',
  'n <- as.integer(a[5])',
  'quiet <- function(e) { utils::capture.output(v <- suppressMessages(e)); v }',
  'set.seed(9)',
  'out <- tempfile(); dir.create(out)',
  'elapsed <- system.time(value <- switch(what,',
  '  participants = {',
  '    per <- 1400 %/% n',
  '    quiet(generateCI(rep(seq_len(per), n), sample(c(-1, 1), per * n, TRUE), "face", rdata,',
  '                     participants = rep(seq_len(n), each = per), save_as_png = FALSE,',
  '                     n_cores = cores))',
  '  },',
  '  stimuli = {',
  '    quiet(generateStimuli2IFC(list(face = rdata), n_trials = n, img_size = 512, stimulus_path = out,',
  '                              seed = 1, ncores = cores, save_as_png = FALSE, save_rdata = FALSE,',
  '                              return_as_dataframe = TRUE))',
  '  },',
  '  zmap = {',
  '    quiet(generateCI(seq_len(n), sample(c(-1, 1), n, TRUE), "face", rdata, save_as_png = FALSE,',
  '                     zmap = TRUE, zmapmethod = "t.test", zmapdecoration = FALSE,',
  '                     zmaptargetpath = out, n_cores = cores))$zmap',
  '  },',
  '  reference = {',
  '    quiet(generateReferenceDistribution2IFC(rdata, iter = 100, ncores = cores, response_seed = 1,',
  '                                            save_rdata = FALSE, reference_method = "images"))',
  '  }))[["elapsed"]]',
  'saveRDS(value, result)',
  'cat(elapsed, "\\n")'
), job)

sampler <- tempfile(fileext = ".sh")
writeLines(c(
  '#!/bin/bash',
  'out=$(mktemp)',
  'setsid Rscript "$@" > "$out" 2>/dev/null &',
  'pid=$!',
  'peak=0',
  'while kill -0 $pid 2>/dev/null; do',
  '  s=$(ps -o rss= --sid $pid 2>/dev/null | awk \'{t+=$1} END {print int(t/1024)}\')',
  '  [ "${s:-0}" -gt "$peak" ] && peak=$s',
  '  sleep 0.1',
  'done',
  'echo "$(grep -v "^starting worker" "$out" | tail -1) $peak"',
  'rm -f "$out"'
), sampler)

measure <- function(version, what, rdata, cores, n) {
  result <- tempfile(fileext = ".rds")
  line <- system2("bash", c(sampler, job, what, shQuote(rdata), cores, result, n),
                  env = paste0("R_LIBS=", libs[[version]]), stdout = TRUE)
  fields <- strsplit(trimws(tail(line, 1)), " +")[[1]]
  list(seconds = as.numeric(fields[1]), peak_MB = as.numeric(fields[2]), result = readRDS(result))
}
```

## The stimulus files

Two synthetic files at 512 pixels and the default 5 scales: 1400 trials
for the participant loop, 200 for the others.

``` r
n <- 512
x <- matrix(seq(-1, 1, length.out = n), n, n, byrow = TRUE)
y <- -matrix(seq(-1, 1, length.out = n), n, n)
face <- exp(-(x^2 / 0.45 + y^2 / 0.75))
base_png <- tempfile(fileext = ".png")
png::writePNG((face - min(face)) / diff(range(face)), base_png)

.libPaths(c(libs[["branch"]], .libPaths()))
Sys.setenv(R_LIBS = libs[["branch"]])
library(rcicr)
stimulus_file <- function(n_trials) {
  dir <- tempfile("stimuli")
  invisible(capture.output(generateStimuli2IFC(list(face = base_png), n_trials = n_trials, img_size = n,
                                               stimulus_path = dir, seed = 1, ncores = 4,
                                               save_as_png = FALSE)))
  list.files(dir, "\\.Rdata$", full.names = TRUE)
}
rdata_1400 <- stimulus_file(1400)
rdata_200 <- stimulus_file(200)
```

## One render, double and integer indices

``` r
stored <- new.env()
load(rdata_200, envir = stored)
double_basis <- stored$p
integer_basis <- rcicr:::renderingBasis(double_basis)
params <- stored$stimuli_params$face
c(double_MB = as.numeric(object.size(double_basis$patchIdx)) / 2^20,
  integer_MB = as.numeric(object.size(integer_basis$patchIdx)) / 2^20)
 double_MB integer_MB 
 120.00021   60.00021 
identical(lapply(1:20, function(i) generateNoiseImage(params[i, ], double_basis)),
          lapply(1:20, function(i) generateNoiseImage(params[i, ], integer_basis)))
[1] TRUE
renders <- sapply(1:2, function(run) {
  c(double = system.time(for (i in 1:20) generateNoiseImage(params[i, ], double_basis))[["elapsed"]],
    integer = system.time(for (i in 1:20) generateNoiseImage(params[i, ], integer_basis))[["elapsed"]])
})
renders
         [,1]  [,2]
double  7.127 6.846
integer 5.031 4.956
rm(stored, double_basis, integer_basis, params)
```

## The participant loop

`generateCI(participants = )` with 1400 trials split evenly between the
participants.

``` r
grid <- expand.grid(cores = c(1, 2, 4), participants = c(4, 20, 100))
participant_rows <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  runs <- lapply(c("main", "branch"), measure, what = "participants", rdata = rdata_1400,
                 cores = grid$cores[i], n = grid$participants[i])
  data.frame(participants = grid$participants[i], cores = grid$cores[i],
             main_s = runs[[1]]$seconds, branch_s = runs[[2]]$seconds,
             main_MB = runs[[1]]$peak_MB, branch_MB = runs[[2]]$peak_MB,
             identical = identical(runs[[1]]$result, runs[[2]]$result))
}))
knitr::kable(participant_rows)
```

| participants | cores | main_s | branch_s | main_MB | branch_MB | identical |
|-------------:|------:|-------:|---------:|--------:|----------:|:----------|
|            4 |     1 |  4.185 |    3.994 |     953 |       801 | TRUE      |
|            4 |     2 |  6.124 |    5.937 |    2138 |      1968 | TRUE      |
|            4 |     4 |  6.983 |    6.482 |    3569 |      3217 | TRUE      |
|           20 |     1 | 10.006 |    8.131 |     977 |       815 | TRUE      |
|           20 |     2 |  9.435 |    8.144 |    2544 |      2272 | TRUE      |
|           20 |     4 | 10.260 |    8.265 |    4601 |      3743 | TRUE      |
|          100 |     1 | 39.929 |   30.126 |    1335 |      1296 | TRUE      |
|          100 |     2 | 25.369 |   20.481 |    2733 |      2272 | TRUE      |
|          100 |     4 | 19.117 |   15.332 |    4270 |      4033 | TRUE      |

## The other loops

200 trials each: generating the stimuli, the t-test z-map, and an
`"images"` reference distribution of 100 draws.

``` r
loop_rows <- do.call(rbind, lapply(c("stimuli", "zmap", "reference"), function(what) {
  rdata <- if (what == "stimuli") base_png else rdata_200
  do.call(rbind, lapply(c(1, 4), function(cores) {
    runs <- lapply(c("main", "branch"), measure, what = what, rdata = rdata, cores = cores, n = 200)
    data.frame(loop = what, cores = cores,
               main_s = runs[[1]]$seconds, branch_s = runs[[2]]$seconds,
               main_MB = runs[[1]]$peak_MB, branch_MB = runs[[2]]$peak_MB,
               identical = identical(runs[[1]]$result, runs[[2]]$result))
  }))
}))
knitr::kable(loop_rows)
```

| loop      | cores | main_s | branch_s | main_MB | branch_MB | identical |
|:----------|------:|-------:|---------:|--------:|----------:|:----------|
| stimuli   |     1 | 73.137 |   51.973 |    1841 |      1659 | TRUE      |
| stimuli   |     4 | 31.596 |   23.900 |    5184 |      4718 | TRUE      |
| zmap      |     1 | 87.786 |   67.970 |    1727 |      1723 | TRUE      |
| zmap      |     4 | 44.205 |   38.319 |    4573 |      4169 | TRUE      |
| reference |     1 | 85.477 |   65.235 |    1847 |      1831 | TRUE      |
| reference |     4 | 41.235 |   35.337 |    4757 |      4255 | TRUE      |

## What a worker costs

The extra peak memory of a run on 4 cores over the same run on 1 core,
divided by the 3 extra processes. Each run holds the parent too, so this
is the price of a worker, not of a process.

``` r
per_worker <- function(rows, by) {
  serial <- rows[rows$cores == 1, ]
  parallel_rows <- rows[rows$cores == 4, ]
  data.frame(run = parallel_rows[[by]], main_MB = round((parallel_rows$main_MB - serial$main_MB) / 3),
             branch_MB = round((parallel_rows$branch_MB - serial$branch_MB) / 3), row.names = NULL)
}
knitr::kable(rbind(per_worker(participant_rows, "participants"), per_worker(loop_rows, "loop")))
```

| run       | main_MB | branch_MB |
|:----------|--------:|----------:|
| 4         |     872 |       805 |
| 20        |    1208 |       976 |
| 100       |     978 |       912 |
| stimuli   |    1114 |      1020 |
| zmap      |     949 |       815 |
| reference |     970 |       808 |

## What this shows, and what it does not

These are one machine’s numbers from one run each; another machine’s
will differ, and the ratios between `main` and the branch are what
carries over. The `identical` columns are `TRUE` for every
configuration, so the change moves no number in any of them.
