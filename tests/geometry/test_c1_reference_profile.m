function test_c1_reference_profile()
%TEST_C1_REFERENCE_PROFILE Self-consistency of the Python reference's normals and normal offset.
%   Reference: tests/reference/c1_reference.py (--profile, --offset), written from the deck text.
%   No MATLAB kernel code is involved; later tasks compare the kernel against these values.
%   Asserted, all exact by construction, at 1e-14 (the reference evaluates B-spline and arc
%   derivatives analytically, so the residuals are a few hundred roundoffs at most):
%     - outward normals are unit vectors orthogonal to the unit tangents;
%     - the normal turned about the axis equals unit(S_s x S_theta) of the revolved patch;
%     - C1 oracle (deck knowledge, used only as an independent check): the neck is the vertical
%       line x = -0.1, so the offset half-width at t = 0.0762 m is 0.1 - 0.0762 m for z in the neck;
%       a wall of t = 0.12 m, thicker than the half-width, leaves no valid void in the neck.
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    script = fullfile(root, 'tests', 'reference', 'c1_reference.py');
    deck = fullfile(root, 'Input', 'C1.ms2');
    tol = 1e-14;

    prof = run(script, deck, '--profile 201').profile;
    gap = [prof.revolved_normal_gap];  % null (point on the axis) decodes to [] and drops out
    n_unit = max(abs(hypot([prof.nx], [prof.nz]) - 1));
    n_orth = max(abs([prof.nx] .* [prof.tx] + [prof.nz] .* [prof.tz]));
    fprintf('profile points %d: |n|-1 max %.3e, n.t max %.3e, turned normal vs S_s x S_theta max %.3e (%d off-axis points)\n', ...
            numel(prof), n_unit, n_orth, max(gap), numel(gap));
    require(n_unit < tol && n_orth < tol, 'normals unit and orthogonal to tangents');
    require(max(gap) < tol, 'turned normal equals S_s x S_theta');

    z_neck = [-3.0 -2.0 -1.0 -0.5 0.0 0.5 0.9];
    z_neck = z_neck(z_neck > -0.3671875 & z_neck < 0.95);
    t = 0.0762;
    off = run(script, deck, sprintf('--offset %.17g %s', t, sprintf(' %.17g', z_neck)));
    got = arrayfun(@(r) r.x_half_width, off);
    fprintf('neck offset half-width at t = %.4f: %s (expected %.6f)\n', t, mat2str(got, 10), 0.1 - t);
    require(max(abs(got - (0.1 - t))) < tol, 'neck offset half-width equals 0.1 - t');

    closed = run(script, deck, sprintf('--offset 0.12 %.17g', 0.5));
    fprintf('t = 0.12 at z = 0.5: valid crossings %d, raw crossings %d\n', ...
            numel(closed.x_half_width_valid), numel(closed.x_half_width_raw));
    require(isempty(closed.x_half_width_valid), 'a wall thicker than the neck half-width leaves no void');

    swept = run(script, deck, sprintf('--offset %.17g %s', t, sprintf(' %.17g', linspace(-3.1, 1.02, 60))));
    n_valid = arrayfun(@(r) numel(r.x_half_width_valid), swept);
    n_raw = arrayfun(@(r) numel(r.x_half_width_raw), swept);
    fprintf('t = %.4f, 60 heights from -3.1 to 1.02: heights with a fold or closing (raw ~= valid) %d\n', t, sum(n_raw ~= n_valid));
end

function out = run(script, deck, args)
    [status, txt] = system(sprintf('python3 -I "%s" --deck "%s" %s', script, deck, args));
    if status ~= 0
        error('test_c1_reference_profile:python', 'reference evaluator failed:\n%s', txt);
    end
    out = jsondecode(txt);
end

function require(cond, what)
    if ~cond, error('test_c1_reference_profile:assert', 'failed: %s', what); end
end
