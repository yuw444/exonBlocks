# exonBlocks AGENTS.md

This file contains instructions for AI agents operating in the `exonBlocks` repository.
This is an R package with C extensions (`src/*.c`) for analyzing single-cell RNA-seq data (BAM files).

## 0. HPC Environment Setup (IMPORTANT)

R is only available through the HPC module system. **R is NOT in the default PATH.**

### Loading R

```bash
# Check available versions
module avail R

# Load R 4.5.0 (required version)
module load R/4.5.0

# R is now available at:
which R
# /hpc/apps/R/4.5.0/bin/R
```

### User Library (No Root Required)

Install packages to your user library:

```bash
# Create user library if needed
mkdir -p ~/R/x86_64-pc-linux-gnu-library/4.5.0

# Install package
R CMD INSTALL --library=~/R/x86_64-pc-linux-gnu-library/4.5.0 .

# Or set library path in R
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.5.0")
```

### Key Paths

| Component | Path |
|-----------|------|
| R Runtime | `/hpc/apps/R/4.5.0/` |
| R Headers | `/hpc/apps/R/4.5.0/lib64/R/include` |
| htslib | `/hpc/apps/htslib/1.22.1/include` |
| User Library | `~/R/x86_64-pc-linux-gnu-library/4.5.0/` |

### VS Code Configuration

The `.vscode/` directory contains pre-configured settings:

- `.vscode/settings.json` - R LSP and C/C++ settings with correct include paths
- `.vscode/c_cpp_properties.json` - IntelliSense configuration
- `.vscode/tasks.json` - Build tasks (Ctrl+Shift+P → "R: Install Package")
- `.vscode/r.attach.sh` - R loader script for VS Code R extension

### Compile/Test Commands

```bash
# Load R module FIRST (required)
module load R/4.5.0

# Install package
R CMD INSTALL --library=~/R/x86_64-pc-linux-gnu-library/4.5.0 .

# Run tests
R -e 'devtools::test()'

# Run specific tests
R -e 'devtools::test(filter = "build_matrix")'

# Generate documentation
R -e 'devtools::document()'

# Check package
R CMD check .

# Verify C syntax
gcc -fsyntax-only -c src/scan_core.c \
    -I/hpc/apps/R/4.5.0/lib64/R/include \
    -I/hpc/apps/htslib/1.22.1/include
```

### Troubleshooting

| Error | Solution |
|-------|----------|
| "R command not found" | Run `module load R/4.5.0` first |
| "No permission to install" | Use `--library=~/R/x86_64-pc-linux-gnu-library/4.5.0` |
| "R.h not found" | Check include path has `/hpc/apps/R/4.5.0/lib64/R/include` |
| "htslib/hts.h not found" | Check include path has `/hpc/apps/htslib/1.22.1/include` |

See `COMPILE_ENV.md` for detailed environment documentation.

## 1. Environment & Build

**System Requirements:**
- **htslib** is REQUIRED. The C code depends on `<htslib/hts.h>`.
- Ensure `pkg-config --libs htslib` works or `HTSLIB_DIR` is set if compilation fails.

**Build Commands:**
- **Install Package:**
  ```bash
  module load R/4.5.0
  R CMD INSTALL --library=~/R/x86_64-pc-linux-gnu-library/4.5.0 .
  # OR
  R -e 'devtools::install()'
  ```
- **Rebuild C Code:**
  The `src/` directory contains `Makevars`. If you modify `.c` files, you must reinstall/rebuild.

**Verification:**
- **Check Package Integrity:**
  ```bash
  module load R/4.5.0
  R CMD check .
  ```

## 2. Testing

The project uses the `testthat` framework (Edition 3).

**Run All Tests:**
```bash
module load R/4.5.0
R -e 'devtools::test()'
# OR
Rscript tests/testthat.R
```

**Run Single Test File:**
```bash
# Example: running test-scan.R
module load R/4.5.0
R -e 'testthat::test_file("tests/testthat/test-scan.R")'
```

**Run Single Test Case:**
```bash
# Filter by test description regex
module load R/4.5.0
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
- If `htslib/hts.h` not found: Check include paths in `.vscode/settings.json` and `.vscode/c_cpp_properties.json`.
- On HPC: Ensure paths `/hpc/apps/R/4.5.0/lib64/R/include` and `/hpc/apps/htslib/1.22.1/include` are set.
