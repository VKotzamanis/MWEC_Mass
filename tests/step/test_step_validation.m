function test_step_validation()
%TEST_STEP_VALIDATION  The writer rejects an open solid shell, a broken loop and a repeated body name.

addpath(fileparts(mfilename('fullpath')));
b = stp_new();
[b, shell] = stp_box(b, [0 0 0], [1 1 1]);
good = b;
good.bodies(1) = struct('name', 'cube', 'kind', 'solid', 'shells', {{shell}});
file = [tempname() '.step'];
cleanup = onCleanup(@() remove_if_present(file));

open_solid = good;
open_solid.bodies(1).shells = {shell(1:5)};
broken_loop = good;
broken_loop.faces(1).loops{1} = broken_loop.faces(1).loops{1}([1 3 2 4]);
twice = good;
twice.bodies(2) = twice.bodies(1);

cases = {open_solid, 'open solid shell'; broken_loop, 'broken loop'; twice, 'repeated body name'};
for k = 1:size(cases, 1)
    id = '';
    try
        mwecmass.output.step.write_step(cases{k, 1}, file);
    catch err
        id = err.identifier;
    end
    fprintf('validation: %s -> %s\n', cases{k, 2}, id);
    if ~strcmp(id, 'mwecmass:step:invalid')
        error('test_step_validation: %s was not rejected', cases{k, 2});
    end
end
end

function remove_if_present(file)
if exist(file, 'file')
    delete(file);
end
end
