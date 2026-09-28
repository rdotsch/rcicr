# R setup for ChatGPT Work containers

The existing Claude container script fails in this Ubuntu 24.04 Work container: `/proc/self/uid_map` and `gid_map` each expose only ID 0, and apt cannot switch to `_apt`. A per-command `APT::Sandbox::User=root` setting allows signed package-index downloads from the configured Ubuntu snapshot. Package installation is still being tested.

Add a separate `tools/setup-work-container-r.sh` for this environment and document its invocation in CONTRIBUTING.md. Keep the Claude script unchanged. Reuse signed Ubuntu packages, keep TLS and package signature verification, and avoid persistent apt configuration changes. Scope any apt user override to the single-ID container. Install rcicr for parallel workers and make reruns safe. Prefer the real yesno package if available; otherwise clearly identify the existing fail-loud stub limitation.

The main uncertainty is whether dpkg maintainer scripts work inside this user namespace. Validate the resulting script with an initial installation and a rerun, load rcicr in both the parent and parallel workers, then run the test suite. Record environmental limitations rather than claiming a complete CRAN build environment. Package R code, numerical output, CI workflows, and release gates are outside this PR.
