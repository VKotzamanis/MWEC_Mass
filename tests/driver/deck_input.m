function in = deck_input(deck_file, mode, cache_dir, repo_root)
%DECK_INPUT Author inputs for a small stand-in deck: geometry-only build_config on it takes seconds.
%   deck_file is an absolute path; build_config resolves it as Input/<in.files.ms2_file>, so the
%   path is given relative to Input/. The BEM cache and HAMS paths are not used (hams_dir empty).
%   repo_root (default: this repository) is the folder whose Input/ build_config resolves against;
%   absolute POSIX paths only. Test helper, not a test.
  if nargin < 4
    repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
  end
  input_dir = fullfile(repo_root, 'Input');
  depth = numel(strsplit(strrep(input_dir, '\', '/'), '/')) - 1;
  in = WEC_User_Input();
  in.files.ms2_file = [repmat('../', 1, depth) strrep(deck_file(2:end), '\', '/')];
  in.files.geometry_cache_dir = cache_dir;
  in.materials.realisation_type = mode;
  in.materials.modular_precast.wall_height = 1.0;
  in.materials.modular_precast.n_sub = 101;
  in.geometry.n_z_levels = 20;
  in.geometry.num_ballast_sections = 4;
  in.bem.hams_dir = '';
end
