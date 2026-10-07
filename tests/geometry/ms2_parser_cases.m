function cases = ms2_parser_cases()
%MS2_PARSER_CASES Decks and parameter grids of the MS2Parser identity test and its fixture.
%   cases(k).name, .deck (path from the repository root) and .opts (see ms2_parser_record).
%   The fixture tests/fixtures/ms2_parser_reference.mat holds the outputs of the parser as it stood
%   before the entity references were resolved once; tools/make_ms2_parser_reference.m wrote it.

    edge_t = [1/3, 2/3, 1 - 1e-9, 1 - 2e-10, 1 - 1e-11];
    cases = struct('name', {}, 'deck', {}, 'opts', {});

    o = struct();
    o.sg  = unique([linspace(0, 1, 21), 1/3, 2/3, 1 - 1e-11]);
    o.tg  = unique([linspace(0, 1, 401), edge_t]);
    o.n_u = [100, 37];
    o.n_z = [96, 52];
    cases(end + 1) = struct('name', 'C1', 'deck', fullfile('Input', 'C1.ms2'), 'opts', o);

    o = struct();
    o.sg  = unique([linspace(0, 1, 9), 1/3, 1 - 1e-11]);
    o.tg  = unique([linspace(0, 1, 101), edge_t]);
    o.n_u = 40;
    o.n_z = 30;
    cases(end + 1) = struct('name', 'all_entity_types', ...
                            'deck', fullfile('tests', 'fixtures', 'ms2_all_entity_types.ms2'), 'opts', o);
end
