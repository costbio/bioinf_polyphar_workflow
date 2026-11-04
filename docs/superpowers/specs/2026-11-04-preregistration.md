# Pre-registration — moved

The pre-registration now lives at **`polyphar/data/preregistration.md`**, as package data.

It is read at run time by `polyphar stats`, which records its sha256 in every manifest. Kept
under `docs/` it was not shipped in the Docker image (`docker/Dockerfile` copies `polyphar/`
and `tests/`, not `docs/`), so `preregistration_sha256()` returned `""` and every manifest
written by a containerised run — i.e. every production run — carried an empty provenance
field. The document's only value is fixing the analysis before the numbers exist, which an
empty field destroys silently.

Editing it still requires an `## Amendments` entry and a matching `stats.PREREGISTRATION_SHA256`.
