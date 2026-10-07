function test_precast_inner_ranges()
%TEST_PRECAST_INNER_RANGES  z ranges of the inner sets realise_modules requests, and their reuse.
%   Box fixture on the SK stand-ins (sti_inner_box ignores the range; the provider records it).
%   A set of t below ctx.t_max_hollow covers the whole hollow range and serves every module of
%   that t; a thicker set covers only the span of the modules that use it (their own bound
%   t_max,i allows it, contract F2b and the T6 bound t_max,i); a kept set is reused only when its
%   range covers the span. Asserted exactly: the requested ranges and the number of provider calls.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
config = sti_config('box');
geo = mwecmass.solid.outer_nurbs(config.ms2_model);
e = config.strip_edges;
t_min = 0.0762;
record('reset');
ctx = struct('geo', geo, 'rho', struct('uhpc', 2500, 'air', 1.2), 'z_hollow', [e(1) e(end)], ...
    't_max_hollow', 0.2, 'knots_from', [], 'sets', [], 'set_ranges', zeros(0, 2), ...
    'inner_fn', @(t, knots_from, zr) record(geo, t, t_min, zr));
design = struct('mode', 'modular_precast', 'edges', e, 'vs', 0, 't', [0.1; 0.1; 0.3], ...
    'z_ballast', -2.2, 'solid_modules', zeros(1, 0));

[ev, ctx] = mwecmass.realise.modular_precast.realise_modules(ctx, design);
calls = record('log');
fprintf('first build: %d provider calls, ranges %s\n', size(calls, 1), mat2str(calls, 6));
check(isequal(calls, [0.1, e(1), e(end); 0.3, e(3), e(4)]), 'ranges: hollow range below t_max_hollow, span above');
check(numel(ev.inner) == 2 && isequal(sort([ev.inner.t]), [0.1 0.3]), 'one set per distinct t');

design.t = [0.1; 0.3; 0.3];
[~, ctx] = mwecmass.realise.modular_precast.realise_modules(ctx, design);
calls = record('log');
fprintf('second build: %d provider calls in total, last range %s\n', size(calls, 1), mat2str(calls(end, :), 6));
check(size(calls, 1) == 3 && isequal(calls(3, :), [0.3, e(2), e(4)]), ...
    'a kept thick set that does not cover the new span is not reused');

design.t = [0.1; 0.1; 0.3];
[~, ctx] = mwecmass.realise.modular_precast.realise_modules(ctx, design);
check(size(record('log'), 1) == 3, 'kept sets that cover the spans are reused');
check(size(ctx.set_ranges, 1) == numel(ctx.sets), 'one range per kept set');
end

function out = record(geo, t, t_min, zr)
persistent calls
if ischar(geo)
    if strcmp(geo, 'reset')
        calls = zeros(0, 3);
    end
    out = calls;
    return
end
calls(end + 1, :) = [t, zr(:)'];
out = sti_inner_box(geo, t, t_min);
end

function setup(root)
if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:shadowed-function');
    addpath(fullfile(root, 'tests', 'octave_shims'));
end
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'tests', 'standins'), '-end');
addpath(fullfile(root, 'tests', 'standins', 'fixtures'), '-end');
end

function check(cond, msg)
if ~cond
    error('test_precast_inner_ranges:fail', '%s', msg);
end
end
