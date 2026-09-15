# Avian FASTA wrapper: file handling review

Reviewed: 2026-09-15.

Wrapper: `/home/rasmuskopperud.riis/Coding/flu-wrappers/avianseq_fasta_wrapper.sh`.

## Scope and status

This is a limited review of the shell interface, local file safeguards, and
upload requests. It is
not an end-to-end pipeline validation or an analysis launch guide. The biological
analysis stages were not tested or modified. No analysis, reference download,
or N-drive transfer was performed.

The wrapper selects a remote pipeline revision. Changes in this local pipeline
checkout therefore do not establish what an operational run will execute.

## File-handling fixes

| Problem | Updated behavior |
| --- | --- |
| Recursive upload selected everything under `$HOME/out_fluseq`, including unrelated runs and previous results. | Select only the current run directory. |
| Report upload used both an initial remote directory and an additional relative change into that directory. | Use the initial remote directory once. |
| Missing report CSVs could go unnoticed before moving results and requesting uploads. | Stop before any result move or upload when no report CSV matches exist. |
| A timestamped archive destination could already exist. | Stop rather than move previous results into that existing destination. |
| Local upload paths were unquoted. | Quote the local paths in upload commands. |
| Help required a working Conda installation. | Handle help and argument validation before loading Conda. |
| Long options were advertised but not parsed. | Accept short options, attached short values, `--option VALUE`, and `--option=VALUE`; reject unknown options and missing values. |
| Existing temporary input was deleted before downloading. | Stop and preserve the existing directory, file, or symlink. |
| Two copies of this wrapper could update shared files concurrently. | Hold an exclusive lock for the invocation; refuse a second copy. |
| Validation mode could be mistaken for a dry run. | State in help and the execution log that report CSVs are uploaded. |

## Operator notes

- Display help without loading Conda or contacting external services:

  ```bash
  /home/rasmuskopperud.riis/Coding/flu-wrappers/avianseq_fasta_wrapper.sh --help
  ```

- `flock` is now required. The wrapper uses
  `/mnt/tempdata/fasta_fluseq/.avianseq_fasta_wrapper.lock`; the service account
  must be able to open that file. Do not delete the lock file to bypass an active
  invocation.
- If `/mnt/tempdata/fasta_fluseq/<run>` already exists, execution stops before
  downloading or updating helper scripts. Review and archive that input before
  reusing the run name. A dangling symlink is also treated as an existing path.
- New local output is expected at `$HOME/<run>/`. After the report check succeeds,
  it is moved to `$HOME/out_fluseq/<run>/`.
- Existing archived results are retained as
  `$HOME/out_fluseq/<run>.previous.<timestamp>/`.
- Report CSVs are expected in the run's `reporthuman/` directory.
- Validation mode still uploads report CSVs to its configured validation
  destination; it skips the full result-directory upload. It is not a dry run.
- If publication fails, retain the local output and console log. A failed
  transfer may leave partial remote results; this wrapper does not verify remote
  file checksums or provide a separate upload-retry command.
- Keep credentials out of logs and shared documentation.

## Other observations requiring follow-up

- Helper scripts are updated during execution, and the default pipeline branch
  is moving. The hard-coded release label does not establish the code revision
  actually executed.
- Reference downloads still overwrite a shared directory. The lock protects
  against another copy of this wrapper only; other wrappers and manual writers
  do not participate automatically. Per-run reference isolation and incomplete
  downloads remain unresolved.
- Publication still moves output into the local archive before uploading.
  Separate upload recovery and remote verification remain unresolved.

## Verification

`bash -n` passed. Four isolated checks exercised the publication block with real
temporary local files and a simulated `smbclient`: normal publication,
validation-mode publication, missing reports, and archive-name collision.
The checks confirmed preservation of unrelated local results and previous
results, plus upload requests restricted to the intended run and report CSVs.

These checks establish local shell behavior only. Actual SMB transfers, remote
permissions, pipeline execution, and scientific correctness remain unverified.

Additional isolated checks passed for help, equivalent option forms, malformed
arguments, invalid run paths, preservation of existing temporary input and
dangling symlinks, and lock contention/release. The parser and lock checks used
extracted shell sections with temporary fixtures; they did not execute the
analysis or download stages.
