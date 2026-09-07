#!/usr/bin/env bash
# One-time setup on macOS.
#
# The trap this script exists for: PhysiCell needs OpenMP, and Apple's `g++` is
# clang wearing a g++ name badge. Clang on macOS ships without OpenMP, so the
# build fails on -fopenmp with an error that does not mention OpenMP clearly.
# The fix is Homebrew GCC plus PHYSICELL_CPP, which PhysiCell's Makefile honours.
#
# Usage:  bash scripts/setup_macos.sh
# Then:   source .envrc   (written by this script)   and   ./bin/glioma check

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PHYSICELL_VERSION="1.14.2"
PHYSICELL_DIR="${REPO_ROOT}/PhysiCell"

say() { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
die() { printf '\n\033[31mERROR: %s\033[0m\n' "$*" >&2; exit 1; }

[[ "$(uname -s)" == "Darwin" ]] || die "This script is for macOS. On Linux the system g++ already has OpenMP: just run make."

# --- 1. a compiler that actually has OpenMP ---------------------------------
say "Locating a GCC with OpenMP support"
if ! command -v brew >/dev/null 2>&1; then
    die "Homebrew not found. Install it from https://brew.sh, then rerun this script."
fi

if ! ls "$(brew --prefix)/bin/g++-"* >/dev/null 2>&1; then
    say "Installing GCC via Homebrew (this takes a few minutes)"
    brew install gcc
fi

# Pick the highest g++-N Homebrew installed, rather than guessing a version.
GXX="$(ls "$(brew --prefix)/bin/g++-"* 2>/dev/null | sort -V | tail -1)"
[[ -n "${GXX}" ]] || die "Homebrew GCC still not found after install."
echo "    using ${GXX}"

"${GXX}" -fopenmp -x c++ -E /dev/null >/dev/null 2>&1 \
    || die "${GXX} rejected -fopenmp. Try: brew reinstall gcc"

# --- 2. PhysiCell, pinned -----------------------------------------------------
if [[ -d "${PHYSICELL_DIR}" ]]; then
    say "PhysiCell already present at ${PHYSICELL_DIR}"
else
    say "Cloning PhysiCell ${PHYSICELL_VERSION}"
    git clone --depth 1 --branch "${PHYSICELL_VERSION}" \
        https://github.com/MathCancer/PhysiCell.git "${PHYSICELL_DIR}"
fi

say "Building the template project"
(
    cd "${PHYSICELL_DIR}"
    export PHYSICELL_CPP="${GXX}"
    make template
    make -j"$(sysctl -n hw.perflevel0.logicalcpu 2>/dev/null || sysctl -n hw.ncpu)"
)
[[ -x "${PHYSICELL_DIR}/project" ]] || die "Build finished but ${PHYSICELL_DIR}/project is missing."

# --- 3. Python dependencies ---------------------------------------------------
say "Installing Python dependencies"
python3 -m pip install --user -r "${REPO_ROOT}/requirements.txt"

# --- 4. environment for later shells -----------------------------------------
cat > "${REPO_ROOT}/.envrc" <<EOF
# Written by scripts/setup_macos.sh — source this before running ./bin/glioma
export PHYSICELL_CPP="${GXX}"
export PHYSICELL_ROOT="${PHYSICELL_DIR}"
EOF

chmod +x "${REPO_ROOT}/bin/glioma"

say "Verifying the harness"
export PHYSICELL_ROOT="${PHYSICELL_DIR}"
"${REPO_ROOT}/bin/glioma" check

cat <<EOF

Setup complete.

  source .envrc          # in each new shell
  ./bin/glioma check     # should reprint the table above

Before running the real replicate set, open params/model.yaml and set
run.omp_num_threads once. PhysiCell seeds one RNG stream per thread, so changing
it partway through a set breaks seed reproducibility across that set.

  time ./bin/glioma run --replicates 1 --arms mono_TP53wt   # measure first
EOF
