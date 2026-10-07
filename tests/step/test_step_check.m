function test_step_check()
%TEST_STEP_CHECK tests/step_check.py on gmsh-written solids whose volumes are known in closed form.
%   Closed-form oracles: unit cube 1 m3, cylinder pi*0.5^2*2 m3, box 4^3 minus a sphere of radius 1
%   (64 - 4*pi/3 m3). The exact OpenCASCADE volume of each is analytic to double roundoff, so the
%   gate is a relative 1e-12 (a few thousand roundoffs); the cube's planar facets make its mesh
%   volume exact too. The bounding box gate 2e-7 m is twice the 1e-7 m enlargement that OpenCASCADE
%   adds to every box. The cylinder and void mesh volumes carry chord error and are printed only.
    here = fileparts(mfilename('fullpath'));
    root = fileparts(fileparts(here));
    work = tempname();
    mkdir(work);
    cleanup = onCleanup(@() confirm_and_remove(work));
    run_ok(sprintf('python3 -I "%s" "%s"', fullfile(here, 'make_step_fixtures.py'), work));
    check_py = fullfile(root, 'tests', 'step_check.py');

    cube = report(check_py, fullfile(work, 'cube.step'), '--expect-solids 1 --expect-units METRE --expect-open-edges 0 --expect-volume 1 --volume-rtol 1e-12');
    require(cube.n_solids == 1 && cube.n_open_edges == 0 && isequal(cube.declared_length_units, {'METRE'}), 'cube counts');
    require(abs(cube.solids.volume_mesh - 1) < 1e-12, 'cube mesh volume');
    require(all(abs(cube.bbox(:)' - [0 0 0 1 1 1]) < 2e-7), 'cube bounding box');
    fprintf('cube: volume %.17g, mesh volume %.17g, bbox %s\n', cube.volume_occ_total, cube.solids.volume_mesh, mat2str(cube.bbox, 8));

    cyl = report(check_py, fullfile(work, 'cylinder.step'), sprintf('--expect-solids 1 --expect-open-edges 0 --expect-volume %.17g --volume-rtol 1e-12', pi * 0.25 * 2));
    fprintf('cylinder: volume %.17g (analytic %.17g), mesh volume %.6g\n', cyl.volume_occ_total, pi * 0.5, cyl.solids.volume_mesh);

    vol = 64 - 4 * pi / 3;
    voided = report(check_py, fullfile(work, 'void.step'), sprintf('--expect-solids 1 --expect-open-edges 0 --expect-volume %.17g --volume-rtol 1e-12', vol));
    fprintf('box with spherical void: volume %.17g (analytic %.17g), mesh volume %.6g, faces %d\n', ...
            voided.volume_occ_total, vol, voided.solids.volume_mesh, voided.n_faces);

    sheet = report(check_py, fullfile(work, 'sheet.step'), '--expect-solids 0 --expect-open-edges 4');
    fprintf('sheet: solids %d, open edges %d\n', sheet.n_solids, sheet.n_open_edges);

    % The checker must fail, not pass silently, when an expectation is wrong.
    mm_file = fullfile(work, 'cube_mm.step');
    [status, ~] = system(sprintf('python3 -I "%s" "%s" --expect-units METRE', check_py, mm_file));
    require(status ~= 0, 'a MILLI.METRE file must fail --expect-units METRE');
    [status, ~] = system(sprintf('python3 -I "%s" "%s" --expect-solids 2', check_py, fullfile(work, 'cube.step')));
    require(status ~= 0, 'a wrong solid count must fail');
    fprintf('negative checks: wrong units and wrong solid count both exit non-zero\n');
end

function r = report(check_py, file, flags)
    [status, text] = system(sprintf('python3 -I "%s" "%s" %s', check_py, file, flags));
    if status ~= 0
        error('test_step_check:fail', 'step_check.py failed for %s:\n%s', file, text);
    end
    r = jsondecode(text(find(text == '{', 1):end));
end

function run_ok(cmd)
    [status, text] = system(cmd);
    if status ~= 0, error('test_step_check:fixtures', 'fixture writer failed:\n%s', text); end
end

function require(cond, what)
    if ~cond, error('test_step_check:assert', 'failed: %s', what); end
end

function confirm_and_remove(folder)
    confirm_recursive_rmdir(false, 'local');
    rmdir(folder, 's');
end
