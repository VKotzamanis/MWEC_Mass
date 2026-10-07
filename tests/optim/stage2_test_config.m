function config = stage2_test_config(mode)
%STAGE2_TEST_CONFIG Geometry-only config of the stand-in cylinder deck, made evaluable by
%   properties_3d with one cached draft and arbitrary added-mass inputs (a formulation test needs
%   the code path, not the hydrodynamic values). Test helper, not a test.
  repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
  deck = fullfile(repo_root, 'tests', 'standins', 'fixtures', 'cylinder.ms2');
  addpath(fullfile(repo_root, 'tests', 'driver'));
  cleanup = onCleanup(@() rmpath(fullfile(repo_root, 'tests', 'driver')));
  in = deck_input(deck, mode, '');
  evalc('config = mwecmass.driver.build_config(in, [], struct());');
  config.hydro_drafts = 0;
  config.added_mass_diagonal = [1000, 2000, 3000];
  config.hydro_z_cg = 0;
  config.stage2_algorithm = 'sqp';
end
