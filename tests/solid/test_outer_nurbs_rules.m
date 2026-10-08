function test_outer_nurbs_rules()
%TEST_OUTER_NURBS_RULES  F1 decisions on decks written here: no exact patch, a decimal arc, and the
%   contract F1 order of the rules (1a, 1b) on a box whose top is three constant-z rectangles with a
%   T-junction. A test-only recording stub of fit_z_faces (fit_z_faces_recorder) shows which
%   entries F1 hands to the general path; the deck coordinates are the independent reference.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
fix = fullfile(root, 'tests', 'solid', 'fixtures');
sk = fullfile(root, 'tests', 'standins', 'fixtures');
tmp = tempname();
mkdir(tmp);
stub = fit_z_faces_recorder();
addpath(stub);
cleanup = onCleanup(@() leave(stub, tmp));
opts = struct('t_min', 0.0254);

% ---------------------------------------------------------------- no exact patch: reaches fit_z_faces
% capped_cylinder (T1 fixture) with its axis from (0.3, 0, -1) to T: a revolution about an axis
% that is not vertical, so no patch converts exactly
txt = fileread(fullfile(fix, 'capped_cylinder.ms2'));
txt = strrep(txt, 'Line AX 6 -1 8x4 / * K T ;', ...
    sprintf('FramePoint KT 14 -1 / 0 * * 0.3 0.0 -1.0 ;\nLine AX 6 -1 8x4 / * KT T ;'));
f = fullfile(tmp, 'tilted.ms2');
write_text(f, txt);
model = parse(f);
cache = mwecmass.geometry.precompute_boundary_cache(model, 21);
rec = run_to_stub(@() mwecmass.solid.outer_nurbs(model, cache, opts));
check(numel(rec.faces) == numel(model.visible_surfs) && all(rec.general) && strcmp(rec.opts.kind, 'outer') && ...
    rec.opts.t_min == opts.t_min, 'tilted axis: every patch handed to fit_z_faces');
fprintf('tilted axis (capped_cylinder, axis (0.3, 0, -1) to T): %d of %d entries handed to fit_z_faces, kind %s\n', ...
    nnz(rec.general), numel(rec.general), rec.opts.kind);

% ---------------------------------------------------------------- a decimal arc converts exactly
deck = {
    'MultiSurf 1.44'
    'Units: m kg'
    'Symmetry:  x y'
    'BeginModel;'
    'FramePoint K 14 -1 / 0 * * 0.0 0.0 -3.0 ;'
    'FramePoint F1 14 -1 / 0 * * 1.45 0.0 -3.0 ;'
    'FramePoint CF 14 -1 / 0 * * 1.45 0.0 -2.95 ;'
    'FramePoint F2 14 -1 / 0 * * 1.5 0.0 -2.95 ;'
    'FramePoint F3 14 -1 / 0 * * 1.5 0.0 1.0 ;'
    'FramePoint T 14 -1 / 0 * * 0.0 0.0 1.0 ;'
    'Line keel_line 11 -1 8x4 / * K F1 ;'
    'Arc fillet 11 -1 8x4 / * 2 F1 CF F2 ;'
    'Line wall_line 11 -1 8x4 / * F2 F3 ;'
    'Line top_line 11 -1 8x4 / * F3 T ;'
    'PolyCurve2 profile 11 -1 8x4 / * { keel_line fillet wall_line top_line } ;'
    'Line axis 6 -1 8x4 / * K T ;'
    'RevSurf hull 2 11 36x4 12x4 0 / * profile axis 0.0 90.0 ;'
    'EndModel;'};
f = fullfile(tmp, 'decimal_fillet.ms2');
write_text(f, sprintf('%s\n', deck{:}));
model = parse(f);
geo = mwecmass.solid.outer_nurbs(model);
P = geo.outer;
S = model.eval_point('F1');
Cn = model.eval_point('CF');
E = model.eval_point('F2');
check(numel(P) == 4 && all([P.exact]) && all(strcmp({P.offset_kind}, 'rev_z')), 'decimal fillet: 4 exact rev_z entries');
col = reshape(P(1).surf.ctrl(:, 1, :), [], 3);
check(any(all(col == E, 2)) && any(all(col == S, 2)), 'decimal fillet: the arc ends at the deck points F1 and F2');
fprintf('decimal fillet: 4 entries exact, rev_z; |F2 - CF| - |F1 - CF| = %.3g m, arc ends at F1 and F2 bitwise\n', ...
    norm(E - Cn) - norm(S - Cn));

% ---------------------------------------------------------------- order of the rules 1a, 1b
% SK box with its top drawn as three rectangles: [-1, 0] x [-0.75, 0.75], [0, 1] x [-0.75, 0] and
% [0, 1] x [0, 0.75]; the corner (0, 0, 0.5) of the last two lies inside a boundary of the first
txt = fileread(fullfile(sk, 'box.ms2'));
pts = sprintf(['FramePoint M0 14 -1 / 0 * * 0.0 0.0 0.5 ;\nFramePoint MS 14 -1 / 0 * * 0.0 -0.75 0.5 ;\n' ...
    'FramePoint MN 14 -1 / 0 * * 0.0 0.75 0.5 ;\nFramePoint ME 14 -1 / 0 * * 1.0 0.0 0.5 ;\n' ...
    'Line TA1 11 -1 8x4 / * T4 MS ;\nLine TA2 11 -1 8x4 / * T3 MN ;\n' ...
    'Line TB1 11 -1 8x4 / * MS T1 ;\nLine TB2 11 -1 8x4 / * M0 ME ;\nLine TC2 11 -1 8x4 / * MN T2 ;\n' ...
    'RuledSurf top_a 2 11 36x4 12x4 0 / * TA1 TA2 ;\nRuledSurf top_b 2 11 36x4 12x4 0 / * TB1 TB2 ;\n' ...
    'RuledSurf top_c 2 11 36x4 12x4 0 / * TB2 TC2 ;\nEndModel;']);
txt = regexprep(txt, 'RuledSurf box_top [^\n]*\n', '');
txt = strrep(txt, 'EndModel;', pts);
f = fullfile(tmp, 'top_t_junction.ms2');
write_text(f, txt);
model = parse(f);
cache = mwecmass.geometry.precompute_boundary_cache(model, 21);
rec = run_to_stub(@() mwecmass.solid.outer_nurbs(model, cache, opts));
names = {rec.faces.name};
tops = ismember(names, {'top_a', 'top_b', 'top_c'});
check(numel(names) == 8 && nnz(tops) == 3, 'top T-junction: entries');
check(all(rec.general(tops)) && ~any(rec.general(~tops)), ...
    'top T-junction: general mask %s for %s', mat2str(rec.general), strjoin(names, ' '));
check(all([rec.faces(~tops).exact]), 'top T-junction: walls and bottom exact');
for k = 1:numel(names)
    fprintf('  %-11s general %d\n', names{k}, rec.general(k));
end
fprintf('top T-junction: the three tops join one flat region (general), walls and bottom stay exact\n');
end

function rec = run_to_stub(f)
global FIT_Z_FACES_RECORD
FIT_Z_FACES_RECORD = [];
try
    f();
catch err
    if ~strcmp(err.identifier, 'fit_z_faces_stub:reached')
        error('test_outer_nurbs_rules:fail', 'expected the fit_z_faces stub, got %s (%s)', err.identifier, err.message);
    end
    rec = FIT_Z_FACES_RECORD;
    return
end
error('test_outer_nurbs_rules:fail', 'fit_z_faces was not called');
end

function leave(stub, tmp)
rmpath(stub);
rmdir(stub, 's');
if exist(tmp, 'dir')
    rmdir(tmp, 's');
end
end

function m = parse(f)
evalc('m = mwecmass.geometry.MS2Parser.parse(f);');
end

function write_text(f, txt)
fid = fopen(f, 'w');
fprintf(fid, '%s', txt);
fclose(fid);
end

function setup(root)
if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:shadowed-function');
    addpath(fullfile(root, 'tests', 'octave_shims'));
end
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'tests', 'solid'));
end

function check(cond, varargin)
if ~cond
    error('test_outer_nurbs_rules:fail', varargin{:});
end
end
