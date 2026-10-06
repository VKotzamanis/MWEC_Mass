function test_c1_reference_sections()
%TEST_C1_REFERENCE_SECTIONS Compare C1 sections from the MATLAB kernel with the Python deck reference.
%   Reference: tests/reference/c1_reference.py, written from the deck text. The MATLAB side is
%   extract_isocurve_at_z (n_u = 100, as the driver uses) followed by waterplane_properties.
%   Asserted: properties that hold by construction (see below). Printed: all differences.
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(fullfile(root, 'tests', 'octave_shims'), fullfile(root, 'src'));

    deck = fullfile(root, 'Input', 'C1.ms2');
    z = linspace(-3.20, 1.08, 24);
    ref = run_reference(root, deck, z);

    parser = mwecmass.geometry.MS2Parser.parse(deck);
    n_u = 100;
    cache = mwecmass.geometry.precompute_boundary_cache(parser, n_u);

    d_area = zeros(size(z)); d_ixx = d_area; d_iyy = d_area; d_hw = d_area;
    mirror_gap = d_area; oracle_gap = d_area;
    fprintf('%8s %12s %12s %11s %11s %11s %11s\n', 'z [m]', 'A_ref [m2]', 'hw_ref [m]', ...
            'dA/A', 'dIxx/Ixx', 'dIyy/Iyy', 'dhw [m]');
    for k = 1:numel(z)
        pts = mwecmass.geometry.extract_isocurve_at_z(parser, z(k), n_u, cache);
        [A, Ixx, Iyy] = mwecmass.hydrostatics.waterplane_properties(pts);
        hw = max(abs(pts(:, 1)));
        d_area(k) = (A - ref(k).area) / ref(k).area;
        d_ixx(k) = (Ixx - ref(k).Ixx) / ref(k).Ixx;
        d_iyy(k) = (Iyy - ref(k).Iyy) / ref(k).Iyy;
        d_hw(k) = hw - ref(k).x_half_width;
        fprintf('%8.4f %12.6f %12.6f %11.3e %11.3e %11.3e %11.3e\n', z(k), ref(k).area, ...
                ref(k).x_half_width, d_area(k), d_ixx(k), d_iyy(k), d_hw(k));

        % Exact by construction: the contour is built from a source quadrant plus coordinate
        % sign flips, so its point set is identical under x -> -x and y -> -y.
        xy = sortrows(pts(:, 1:2));
        mirror_gap(k) = max([max(abs(sortrows([-xy(:, 1), xy(:, 2)]) - xy)(:)), ...
                             max(abs(sortrows([xy(:, 1), -xy(:, 2)]) - xy)(:))]);

        % Independent oracle using C1-specific knowledge: the C1 section is a stadium (two
        % semicircles of radius w about (0, +-1 m) joined by straight sides), whose area is
        % 4 w + pi w^2 with w the reference half-width; checks the reference against itself.
        w = ref(k).x_half_width;
        oracle_gap(k) = abs(ref(k).area - (4 * w + pi * w^2));
    end

    fprintf('max |dA/A|      = %.3e\n', max(abs(d_area)));
    fprintf('max |dIxx/Ixx|  = %.3e\n', max(abs(d_ixx)));
    fprintf('max |dIyy/Iyy|  = %.3e\n', max(abs(d_iyy)));
    fprintf('max |dhw|       = %.3e m\n', max(abs(d_hw)));
    fprintf('max mirror gap  = %.3e m\n', max(mirror_gap));
    fprintf('max stadium-oracle gap of the reference = %.3e m2\n', max(oracle_gap));

    % Bound 0: the point set and its mirror images are the same floating-point numbers.
    if max(mirror_gap) ~= 0
        error('test_c1_reference_sections:mirror', ...
              'contour is not exactly mirror symmetric (gap %.3e m)', max(mirror_gap));
    end
    % Bound 1e-13 m2: the reference sums 32-point Gauss-Legendre pieces of arcs and lines
    % (roundoff only); the stadium area is at most 4*0.6 + pi*0.6^2 = 3.6 m2 here.
    if max(oracle_gap) > 1e-13
        error('test_c1_reference_sections:oracle', ...
              'reference area departs from the stadium formula by %.3e m2', max(oracle_gap));
    end
end

function ref = run_reference(root, deck, z)
    args = sprintf(' %.17g', z);
    cmd = sprintf('python3 -I "%s" --deck "%s"%s', ...
                  fullfile(root, 'tests', 'reference', 'c1_reference.py'), deck, args);
    [status, txt] = system(cmd);
    if status ~= 0
        error('test_c1_reference_sections:python', 'reference evaluator failed:\n%s', txt);
    end
    ref = jsondecode(txt);
end
