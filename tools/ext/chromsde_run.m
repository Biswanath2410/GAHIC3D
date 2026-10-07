function chromsde_run(infile, outfile)
  pkg load statistics;
  here = pwd;
  cd(fullfile(fileparts(mfilename('fullpath')), 'ChromSDE', 'program'));
  addpath(genpath(fullfile(pwd, 'helperfunctions')));
  H = load(infile);
  [err, P, ci] = ChromSDE_knownAlpha(H, 0.5, 1);
  cd(here);
  dlmwrite(outfile, P', '\t');
end
