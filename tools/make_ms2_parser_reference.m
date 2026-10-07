% Write tests/fixtures/ms2_parser_reference.mat from the MS2Parser on the current path.
%   octave --no-gui --quiet tools/make_ms2_parser_reference.m
% The committed fixture was written by the parser before the entity references were resolved
% once (T0c); do not regenerate it from a changed parser.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'tests', 'octave_shims'), fullfile(root, 'src'), ...
        fullfile(root, 'tests', 'geometry'));
cases = ms2_parser_cases();
reference = struct();
for k = 1:numel(cases)
    t0 = tic;
    model = mwecmass.geometry.MS2Parser.parse(fullfile(root, cases(k).deck));
    reference.(cases(k).name) = ms2_parser_record(model, cases(k).opts);
    fprintf('%s: %d calls recorded in %.1f s\n', cases(k).name, ...
            numel(reference.(cases(k).name).labels), toc(t0));
end
save('-v7', fullfile(root, 'tests', 'fixtures', 'ms2_parser_reference.mat'), 'reference');
