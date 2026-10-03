# Plan: reject an unrecognised noise_type (#376)

## The defect

`generateNoisePattern()` treats anything but exactly `"gabor"` as sinusoid noise and records the value as given, so `noise_type = "Gabor"` writes sinusoid stimuli into a file that says Gabor.

## The change

`generateNoisePattern()` and `generateStimuli2IFC()` stop unless `noise_type` is exactly `"sinusoid"` or `"gabor"`. `generateStimuli2IFC()` checks up front, before the base images are read and anything is reserved or written.

`match.arg()` is rejected: it partially matches, so `"gab"`, sinusoid noise today, would silently become Gabor. Exact matching only turns wrong calls into errors.

## Compatibility

Every value accepted afterwards produces the same noise as before. Existing `.Rdata` files are not re-validated: nothing downstream reads `noise_type`, the saved basis `p` is used as stored (checked: `grep noise_type R/`).

## Tests

`"Gabor"`, `"gab"`, `NA`, a length-2 vector: error, and for `generateStimuli2IFC()` nothing written. `"sinusoid"`/`"gabor"`: patches identical to before.
