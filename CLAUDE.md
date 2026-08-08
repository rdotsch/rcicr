<!--
This file exists only so that Claude Code loads this repository's conventions.
Claude Code reads CLAUDE.md and does not read AGENTS.md:
https://code.claude.com/docs/en/memory#agents-md

AGENTS.md remains the single source of truth. Put conventions there, not here.

An @-import is used rather than `ln -s AGENTS.md CLAUDE.md` because symlinks
require Administrator privileges or Developer Mode on Windows, and this package
has Windows contributors and a Windows CI runner.

Block-level HTML comments are stripped before the file is loaded, so this note
costs no context. To confirm the file is loading, run /context in a session and
look for CLAUDE.md under "Memory files".
-->

@AGENTS.md
