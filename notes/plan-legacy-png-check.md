# Plan: compare archived stimulus PNGs across the old white-pixel wrap (#419)

## Verified behavior

The current checker includes a pixel whenever both decoded PNG values lie strictly between 0 and 1, then compares `ori - inv` to `noise / 0.6`. Before 1.6.0, the stimulus writer passed values above 1 to `png::writePNG()`; the 1.6.0 NEWS entry and its measured example record that such values wrapped into the interior rather than clipping to white. Algebra gives a counterexample: with base 1 and noise +0.36, the two pre-write values are 1.05 and 0.45. The old encoded original is near 0.05, so the checker sees approximately -0.40 instead of +0.60. Its interior filter does not exclude that pixel.

The exact decoded residual depends on PNG bit depth and integer conversion. A regression test must measure it through `png::writePNG()` and `png::readPNG()`, not assert a hand-derived byte.

## Change

1. Add a test that generates a real current stimulus file with a white base and low `nscales`, rewrites its PNGs with the pre-1.6.0 expression (without `clampUnit()`), and checks that overflow pixels decode to interior values and fail the present comparison. Check that the correct file scores completely and a wrong-seed candidate remains low. Keep the current clamped archive case.
2. Compare a decoded difference with the expected difference under both the direct interpretation and the legacy writer's integer wrap. Derive the period from the PNG encoding actually used in the test (and document any unavoidable ambiguity); do not change stimulus generation or a saved `.Rdata` field.
3. Update the checker help and NEWS with the scope and limits of this recovery check. Ensure `compared` and `share` still have an honest meaning.
4. Run targeted tests, R CMD check and the release comparison. Remove this plan file before marking the PR ready.

## Main risk

Treating any difference modulo one as an agreement may give a wrong candidate extra matches. The test must measure the wrong-seed control and the exact PNG conversion; the implementation should accept only a legacy alias supported by the writer, not an arbitrary approximate integer offset. If the test disproves the proposed treatment, revise the algorithm and this PR rather than weakening the assertion.
