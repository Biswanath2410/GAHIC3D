#!/bin/bash
# Download LorDG, GEM and Chromosome3D (not included in this repository).
# Usage: bash tools/fetch_external_tools.sh
set -e
TOOLS="$(cd "$(dirname "$0")" && pwd)"

fetch() {   # repo commit target
  if [ -d "$3" ]; then echo "exists: $3"; return; fi
  git clone -q https://github.com/$1.git "$3"
  git -C "$3" checkout -q $2
  rm -rf "$3/.git"
  echo "fetched $1 @ ${2:0:7} -> $3"
}

# LorDG (Trieu & Cheng 2017): uses bin/3DDistanceBaseLorentz.jar + lib/commons-math3-3.5.jar
fetch BDM-Lab/LorDG f40f1f7b88ed9458d50d89bd6ac720089a009ab8 "$TOOLS/LorDG"
# GEM (Zhu et al. 2018), MATLAB/Octave code
fetch GuangxiangZhu/GEM f464ef1c24e10d93bf27e1f57011abcce841424e "$TOOLS/ext/GEM"
# Chromosome3D (Adhikari et al. 2016), Perl driver for CNS 1.3
fetch multicom-toolbox/Chromosome3D c6e30d6d27c41eeaf312e5c56eec26ac2933c7ca "$TOOLS/ext/Chromosome3D"

# Chromosome3D: set the CNS path, create the output folder, run '&>' with bash
python3 - "$TOOLS/ext/Chromosome3D/chromosome3D.pl" "$TOOLS/ext/cns_solve_1.3" <<'EOF'
import re, sys
f, cns = sys.argv[1], sys.argv[2]
s = open(f).read()
s = re.sub(r'my \$cns_suite\s*=\s*"[^"]*";', 'my $cns_suite      = "%s";' % cns, s, count=1)
s = s.replace('mkdir $dir_out or confess "ERROR! Could not create output directory $dir_out!" if not -d $dir_out;',
              'mkdir $dir_out ,0777 or die $!; #or confess "ERROR! Could not create output directory $dir_out!" if not -d $dir_out;')
s = s.replace('system("$command &> $log");', 'system("bash", "-c", "$command &> $log");')
open(f, 'w').write(s)
EOF
echo "patched Chromosome3D (CNS path: $TOOLS/ext/cns_solve_1.3)"

if [ ! -d "$TOOLS/ext/cns_solve_1.3" ]; then
  cat <<EOF

CNS 1.3 is still missing (needed only for Chromosome3D). It requires a free academic licence:
  1. request it at http://cns-online.org and download the Linux x86-64 build of CNS 1.3
  2. unpack it to   $TOOLS/ext/cns_solve_1.3
  3. in cns_solve_env and cns_solve_env.sh inside that folder, set CNS_SOLVE to that same path
CNS 1.3 is an x86-64 Linux program; on other systems run the pipeline on an x86-64 Linux machine.
EOF
fi
