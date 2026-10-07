#!/bin/bash
# Install all dependencies on Ubuntu 24.04.
# Usage: bash install_linux.sh
set -e
REPO="$(cd "$(dirname "$0")" && pwd)"
export DEBIAN_FRONTEND=noninteractive

# 1. system packages (R + HSA dependencies, GNU Octave for ChromSDE/GEM, Java for LorDG, tcsh/perl for Chromosome3D/CNS)
if [ "$(dpkg --print-architecture)" = "arm64" ]; then
  # CNS 1.3 is x86-64 only; on arm64 it needs the amd64 libraries and an x86-64 emulator
  sudo sed -i "s/^Types: deb$/Types: deb\nArchitectures: arm64/" /etc/apt/sources.list.d/ubuntu.sources
  printf 'Types: deb\nURIs: http://archive.ubuntu.com/ubuntu/\nSuites: noble noble-updates noble-security\nComponents: main universe\nArchitectures: amd64\n' \
    | sudo tee /etc/apt/sources.list.d/amd64.sources >/dev/null
  sudo dpkg --add-architecture amd64
  sudo apt-get update -qq
  sudo apt-get install -y -qq libc6:amd64 libgfortran5:amd64
fi
sudo apt-get update -qq
sudo apt-get install -y -qq r-base-core r-base-dev r-cran-mass r-cran-matrixstats octave octave-dev octave-statistics \
     tcsh openjdk-17-jre-headless perl build-essential gawk curl git
sudo update-alternatives --set java "$(update-alternatives --list java | grep 17 | head -1)"

# 2. Python environment from environment.yml
if [ ! -x ~/miniforge3/bin/conda ]; then
  curl -fsSL -o /tmp/mf.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh"
  bash /tmp/mf.sh -b -p ~/miniforge3 && rm /tmp/mf.sh
fi
[ -d ~/miniforge3/envs/gahic3d ] || ~/miniforge3/bin/conda env create -q -f "$REPO/environment.yml"

# 3. HSA (R package of the authors' scripts)
sudo R CMD INSTALL "$REPO/tools/ext/HSA"

# 4. ChromSDE mex files for Octave (MATLAB headers it expects are provided as one-line shims)
sudo mkdir -p /opt/mexshim
printf '#include "mex.h"\n' | sudo tee /opt/mexshim/matrix.h >/dev/null
printf '#include <stddef.h>\n' | sudo tee /opt/mexshim/blas.h /opt/mexshim/lapack.h >/dev/null
export CFLAGS="-O2 -fPIC -Wno-implicit-function-declaration -Wno-int-conversion -fpermissive"
for d in helperfunctions SDPT3-4.0/Solver/Mexfun PPA_Semismooth_new/solver/mexfun; do
  ( cd "$REPO/tools/ext/ChromSDE/program/$d"
    for c in *.c; do b=${c%.c}; case $b in *win|*old) continue;; esac
      sed -i "/typedef int mwIndex/d; /typedef int mwSize/d" "$c"
      mkoctfile --mex -I/opt/mexshim "$c" -o "$b.mex" -llapack -lblas >/dev/null 2>&1 || echo "mex skipped: $d/$c"
    done; rm -f ./*.o )
done

# 5. methods not redistributed here (LorDG, GEM, Chromosome3D) + instructions for CNS 1.3
bash "$REPO/tools/fetch_external_tools.sh"

# 6. check
cd "$REPO" && source env.sh && python3 tests/test_core.py
echo "Done. Run: bash run_all.sh <ncores>"
