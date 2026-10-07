function test_ms2_parser_identity()
%TEST_MS2_PARSER_IDENTITY The parser returns the arrays it returned before entity references were
% resolved once. The fixture holds the outputs of the unchanged parser for every curve, snake,
% point and visible surface of Input/C1.ms2 and of a deck with every entity type
% (tests/fixtures/ms2_all_entity_types.ms2), the boundary caches, and extract_isocurve_at_z over
% the hull's z range. Identity is isequaln on each output: bit for bit, errors included.
    here = fileparts(mfilename('fullpath'));
    root = fileparts(fileparts(here));
    addpath(here);
    fx = load(fullfile(root, 'tests', 'fixtures', 'ms2_parser_reference.mat'));
    cases = ms2_parser_cases();
    n_total = 0;
    n_err = 0;
    n_iso = 0;
    for k = 1:numel(cases)
        ref = fx.reference.(cases(k).name);
        t0 = tic;
        model = mwecmass.geometry.MS2Parser.parse(fullfile(root, cases(k).deck));
        got = ms2_parser_record(model, cases(k).opts);
        dt = toc(t0);

        same_model = isequaln(got.model, ref.model);
        if ~same_model
            error('%s: parsed entities or header fields differ from the reference', cases(k).name);
        end
        if ~isequal(got.labels, ref.labels)
            error('%s: the recorded calls differ from the reference', cases(k).name);
        end
        for c = 1:numel(ref.calls)
            a = got.calls{c};
            b = ref.calls{c};
            if ~isequaln(a, b)
                error('%s: "%s" differs from the reference (reference ok=%d, now ok=%d; %s | %s)', ...
                      cases(k).name, ref.labels{c}, b.ok, a.ok, a.msg, b.msg);
            end
            n_err = n_err + ~b.ok;
            n_iso = n_iso + strncmp(ref.labels{c}, 'extract_isocurve_at_z', 21) * b.ok;
        end
        if ~isequaln(got.heights, ref.heights)
            error('%s: iso-z heights differ from the reference', cases(k).name);
        end
        n_total = n_total + numel(ref.calls);
        fprintf('%s: %d calls identical to the reference (%d of them errors, same id and message), %.1f s\n', ...
                cases(k).name, numel(ref.calls), sum(cellfun(@(r) ~r.ok, ref.calls)), dt);
    end
    fprintf('%d calls, %d iso-z contours, %d expected errors; all identical\n', n_total, n_iso, n_err);
    if n_iso < 50
        error('the fixture must hold at least 50 iso-z heights, found %d', n_iso);
    end
end
