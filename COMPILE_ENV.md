# Compile & Test Environment Setup

This document describes the environment setup for developing exonBlocks on the HPC cluster.

## Quick Start

```bash
# Load R module
module load R/4.5.0

# Install package
R CMD INSTALL --library=~/R/x86_64-pc-linux-gnu-library/4.5.0 .

# Run tests
R -e 'devtools::test()'
```

## Environment Details

### R via Environment Modules

R is only available through the HPC module system:

```bash
# Check available R versions
module avail R

# Load R 4.5.0
module load R/4.5.0

# R is now available at:
which R
# /hpc/apps/R/4.5.0/bin/R
```

**Important**: R is NOT in the default PATH. All R commands must be run after `module load R/4.5.0`.

### R Attach Script (VS Code)

The `.vscode/r.attach.sh` script handles R loading for VS Code:

```bash
#!/bin/bash

module load R/4.5.0
R "$@"
```

Configure VS Code to use this script in `settings.json`:
```json
{
    "r.rpath.linux": "${workspaceFolder}/.vscode/r.attach.sh"
}
```

### User Library Path

Install packages to your user library (no root required):

```bash
# Create if not exists
mkdir -p ~/R/x86_64-pc-linux-gnu-library/4.5.0

# Install with
R CMD INSTALL --library=~/R/x86_64-pc-linux-gnu-library/4.5.0 .

# Or set in R
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.5.0")
```

## VS Code Configuration

### settings.json

Located at `.vscode/settings.json`:

```json
{
    // R Terminal Configuration
    "terminal.integrated.env.linux": {
        "PATH": "${env:PATH}:/home/yu89975/.local/bin:/opt/slurm/bin",
        "LD_LIBRARY_PATH": "/opt/slurm/lib64/:${env:LD_LIBRARY_PATH}"
    },
    "r.rpath.linux": "${workspaceFolder}/.vscode/r.attach.sh",
    "r.rterm.linux": "${workspaceFolder}/.vscode/radian.attach.sh",
    "r.plot.useHttpgd": true,
    "r.libPaths": [
        "/home/yu89975/R/x86_64-pc-linux-gnu-library/4.5.0",
        "/scratch/g/chlin/Yu/CD137/src/Impute/lib/"
    ],

    // C/C++ LSP Configuration
    "clangd.arguments": [
        "--header-insertion=never",
        "-I/hpc/apps/R/4.5.0/lib64/R/include",
        "-I/usr/share/R/include"
    ],
    "C_Cpp.default.compilerPath": "/usr/bin/gcc",
    "C_Cpp.default.includePath": [
        "/usr/include",
        "/usr/local/include",
        "/usr/include/c++/11",
        "/usr/lib/gcc/x86_64-linux-gnu/11/include",
        "/hpc/apps/R/4.5.0/lib64/R/include",
        "/hpc/apps/htslib/1.22.1/include",
        "${workspaceFolder}/**"
    ]
}
```

### c_cpp_properties.json

Located at `.vscode/c_cpp_properties.json`:

```json
{
    "configurations": [
        {
            "name": "R/C - exonBlocks",
            "includePath": [
                "${workspaceFolder}/**",
                "/hpc/apps/R/4.5.0/lib64/R/include",
                "/usr/share/R/include",
                "/hpc/apps/htslib/1.22.1/include"
            ],
            "defines": [
                "__STDC__",
                "__STDC_VERSION__=201710L"
            ],
            "compilerPath": "/usr/bin/gcc",
            "cStandard": "c11",
            "intelliSenseMode": "linux-clang"
        }
    ],
    "version": 4
}
```

### tasks.json

Located at `.vscode/tasks.json`:

```json
{
    "version": "2.0.0",
    "tasks": [
        {
            "label": "R: Install Package",
            "type": "shell",
            "command": "R CMD INSTALL .",
            "problemMatcher": ["$gcc"],
            "group": "build"
        },
        {
            "label": "R: Run Tests",
            "type": "shell",
            "command": "R -e 'devtools::test()'",
            "problemMatcher": ["$gcc"],
            "group": "test"
        },
        {
            "label": "R: Generate Docs",
            "type": "shell",
            "command": "R -e 'devtools::document()'",
            "problemMatcher": ["$gcc"]
        },
        {
            "label": "R: Check Package",
            "type": "shell",
            "command": "R CMD check .",
            "problemMatcher": ["$gcc"]
        }
    ]
}
```

**Keyboard Shortcut**: Press `Ctrl+Shift+P` → type task name → Enter

## Dependencies

### System Dependencies

| Dependency | Path | Purpose |
|------------|------|---------|
| R | `/hpc/apps/R/4.5.0/` | R runtime and headers |
| htslib | `/hpc/apps/htslib/1.22.1/` | BAM file handling |
| gcc | `/usr/bin/gcc` | C compiler |

### R Packages

| Package | Version | Purpose |
|---------|---------|---------|
| data.table | >= 1.17.0 | Fast data manipulation |
| dplyr | >= 1.1.4 | Data manipulation |
| tidyr | >= 1.3.1 | Data tidying |
| magrittr | >= 2.0.3 | Pipe operator |
| Matrix | >= 1.6.0 | Sparse matrices |
| GenomicRanges | >= 1.55.0 | Genomic interval handling |
| testthat | >= 3.0.0 | Unit testing |

Install required packages:
```r
install.packages(c("data.table", "dplyr", "tidyr", "magrittr", "Matrix", "GenomicRanges", "testthat"))
```

## Compilation

### C Code Compilation

The package uses R CMD for compilation, which automatically finds R headers:

```bash
# Compile C code (dry run)
R CMD SHLIB --dry-run src/scan_core.c

# Compile with explicit paths
gcc -fsyntax-only -c src/scan_core.c \
    -I/hpc/apps/R/4.5.0/lib64/R/include \
    -I/hpc/apps/htslib/1.22.1/include
```

### Package Installation

```bash
# Full installation (compiles C code)
R CMD INSTALL --library=~/R/x86_64-pc-linux-gnu-library/4.5.0 .

# Verbose output
R CMD INSTALL -v --library=~/R/x86_64-pc-linux-gnu-library/4.5.0 .
```

## Testing

### Run All Tests

```bash
module load R/4.5.0
R -e 'devtools::test()'
```

### Run Specific Tests

```bash
# Filter by file name
R -e 'devtools::test(filter = "build_matrix")'

# Filter by description
R -e 'devtools::test(filter = "cluster")'
```

### Test Results

Expected output format:
```
══ Results ════════════════════════════════════════════════
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 15 ]
```

## Troubleshooting

### "R command not found"

```bash
# Solution: Load R module first
module load R/4.5.0
```

### "No permission to install" during R CMD INSTALL

```bash
# Solution: Use user library
mkdir -p ~/R/x86_64-pc-linux-gnu-library/4.5.0
R CMD INSTALL --library=~/R/x86_64-pc-linux-gnu-library/4.5.0 .
```

### C LSP errors ("R.h not found")

Ensure VS Code settings include:
```json
{
    "clangd.arguments": [
        "-I/hpc/apps/R/4.5.0/lib64/R/include"
    ]
}
```

### htslib errors ("htslib/hts.h not found")

Ensure include path contains htslib:
```json
{
    "C_Cpp.default.includePath": [
        "/hpc/apps/htslib/1.22.1/include"
    ]
}
```

## HPC Notes

- Slurm is available at `/opt/slurm/bin/`
- Jobs can be submitted via `sbatch` for large computations
- User library path: `~/R/x86_64-pc-linux-gnu-library/4.5.0/`
- Shared data path: `/scratch/g/chlin/Yu/CD137/src/Impute/lib/`

## File Locations

```
/scratch/g/chlin/Yu/exonBlocks/
├── .vscode/
│   ├── settings.json           # VS Code settings
│   ├── c_cpp_properties.json   # C/C++ IntelliSense
│   ├── tasks.json              # Build tasks
│   └── r.attach.sh            # R loader script
├── src/
│   ├── scan_core.c            # C backend
│   └── init.c                 # R/C interface
├── R/
│   ├── extract_exon_reads_hts.R
│   ├── cb_umi_exons.R
│   └── build_matrix.R          # New: cell×cluster matrix
├── tests/
│   └── testthat/
│       └── test-build_matrix.R
└── DESCRIPTION                 # Package metadata
```
