# exonBlocks AGENTS.md

This file contains instructions for AI agents operating in the `exonBlocks` repository.
This is an R package with C extensions (`src/*.c`) for analyzing single-cell RNA-seq data (BAM files).

## 1. Environment & Build

**System Requirements:**
- **htslib** is REQUIRED. The C code depends on `<htslib/hts.h>`.
- Ensure `pkg-config --libs htslib` works or `HTSLIB_DIR` is set if compilation fails.

**Build Commands:**
- **Install Package:**
  ```bash
  R CMD INSTALL .
  # OR
  R -e 'devtools::install()'
  ```
- **Rebuild C Code:**
  The `src/` directory contains `Makevars`. If you modify `.c` files, you must reinstall/rebuild.

**Verification:**
- **Check Package Integrity:**
  ```bash
  R CMD check .
  ```

## 2. Testing

The project uses the `testthat` framework (Edition 3).

**Run All Tests:**
```bash
R -e 'devtools::test()'
# OR
Rscript tests/testthat.R
```

**Run Single Test File:**
```bash
# Example: running test-scan.R
R -e 'testthat::test_file("tests/testthat/test-scan.R")'
```

**Run Single Test Case:**
```bash
# Filter by test description regex
R -e 'testthat::test_file("tests/testthat/test-scan.R", filter="scan_bam works")'
```

**Creating Tests:**
- Place new tests in `tests/testthat/test-<name>.R`.
- Use `test_that("description", { ... })`.
- If testing C code, ensure the package is loaded/installed so the DLL is available.

## 3. Code Style & Conventions

### R Code
- **Style:** Mixture of `tidyverse` (dplyr/tidyr) and `data.table`.
  - Use `data.table` for high-performance operations on large datasets.
  - Use `dplyr`/`tidyr` for readability on smaller summaries.
- **Documentation:** Use **roxygen2** comments (`#'`) above functions.
  - Run `devtools::document()` to update `NAMESPACE` and `.Rd` files.
  - Do NOT edit `NAMESPACE` manually.
- **Imports:**
  - Add packages to `DESCRIPTION` (Imports/Suggests).
  - Use `@importFrom package function` in roxygen comments.

### C Code (`src/`)
- **Integration:** Accessed via `.Call()` in R.
- **Registration:** Functions must be registered in `src/init.c`.
- **Memory Management:** Use R's memory management (PROTECT/UNPROTECT) carefully.
- **Headers:** `#include <htslib/hts.h>` for SAM/BAM operations.

### Error Handling
- **R:** Use `stop()`, `warning()`.
- **C:** Use `error()` (from R API) to raise errors back to R safely. Do not use `exit()`.

## 4. Common Tasks

**Adding a new dependency:**
1. Edit `DESCRIPTION`.
2. Run `devtools::document()` if you added `@importFrom`.

**Modifying C Code:**
1. Edit `src/scan_core.c` (or other files).
2. If adding a new function, register it in `src/init.c`.
3. Reinstall package: `R CMD INSTALL .`.
4. Run tests.

**Debugging Build Failures:**
- If `htslib/hts.h` not found: Check `src/Makevars`. It usually uses `pkg-config`.
