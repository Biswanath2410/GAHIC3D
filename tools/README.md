# Competing methods

`install_linux.sh` sets all of them up. Details:

| folder | method | licence | how it gets here |
|---|---|---|---|
| `pastis-0.4.0/` | PASTIS (Varoquaux et al. 2014) | BSD | bundled; bundled iced extension disabled |
| `miniMDS/` | miniMDS (Rieber & Mahony 2017) | MIT (`miniMDS/LICENSE`) | bundled |
| `ext/HSA/` | HSA (Zou et al. 2016), R package of the authors' scripts | GPL ≥ 2 | bundled; `R CMD INSTALL tools/ext/HSA` |
| `ext/ChromSDE/` | ChromSDE (Zhang et al. 2013) with SDPT3, PPA and YALMIP as distributed by its authors | GPL v2 (ChromSDE, SDPT3); YALMIP © J. Löfberg | bundled; mex files built by `install_linux.sh` |
| `LorDG/` | LorDG (Trieu & Cheng 2017) | none stated upstream | `fetch_external_tools.sh` (BDM-Lab/LorDG @ f40f1f7) |
| `ext/GEM/` | GEM (Zhu et al. 2018) | none stated upstream | `fetch_external_tools.sh` (GuangxiangZhu/GEM @ f464ef1) |
| `ext/Chromosome3D/` | Chromosome3D (Adhikari et al. 2016) | none stated upstream | `fetch_external_tools.sh` (multicom-toolbox/Chromosome3D @ c6e30d6) + 3 patches |
| `ext/cns_solve_1.3/` | CNS 1.3 (Brünger et al. 1998), used by Chromosome3D | academic licence | download yourself from cns-online.org |

Wrappers used by `run_baselines.py`: `ext/hsa_run.R`, `ext/chromsde_run.m`.

Chromosome3D patches applied by `fetch_external_tools.sh`: the CNS location, output-folder creation, and running the
`&>` redirection through bash (it fails under `/bin/sh` on Ubuntu).
