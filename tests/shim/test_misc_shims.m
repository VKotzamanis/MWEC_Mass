function test_misc_shims()
%TEST_MISC_SHIMS Known-answer checks of the small Octave shims (all exact, no tolerance).
%   SHA-256 vectors are the FIPS 180-4 examples for "abc" and the empty message.
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(fullfile(root, 'tests', 'octave_shims'));

    expect('contains char', contains('abc_180.ms2', '_180'), true);
    expect('contains absent', contains('abc.ms2', '_180'), false);
    expect('contains cellstr', contains({'a1', 'b2', 'c3'}, '2'), [false true false]);
    expect('contains pattern list', contains('intel/oneapi', {'xx', 'oneapi'}), true);
    expect('contains IgnoreCase', contains('Y=0', 'y=0', 'IgnoreCase', true), true);
    expect('contains case-sensitive', contains('Y=0', 'y=0'), false);

    edges = [0 1 2 4];
    expect('discretize bins', discretize([0 0.5 1 3.9 4 -1 5], edges), [1 1 2 3 3 NaN NaN]);

    expect('issorted strictascend true', issorted([1 2 3], 'strictascend'), true);
    expect('issorted strictascend dup', issorted([1 2 2], 'strictascend'), false);
    expect('issorted default', issorted([1 2 2]), true);

    expect('startsWith single space', startsWith({' a', 'a '}, ' '), [true false]);
    expect('endsWith single space', endsWith({' a', 'a '}, ' '), [false true]);

    stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    if isempty(regexp(stamp, '^\d{8}_\d{6}$', 'once'))
        error('test_misc_shims:datetime', 'datetime format gave "%s"', stamp);
    end
    fprintf('  %-28s %s\n', 'datetime yyyyMMdd_HHmmss', stamp);

    d = java.security.MessageDigest.getInstance('SHA-256');
    d.update(uint8('abc')');
    expect('sha256 abc', lower(sprintf('%02x', typecast(d.digest(), 'uint8'))), ...
           'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad');
    d = java.security.MessageDigest.getInstance('SHA-256');
    expect('sha256 empty', lower(sprintf('%02x', typecast(d.digest(), 'uint8'))), ...
           'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855');
    d = java.security.MessageDigest.getInstance('SHA-256');
    d.update(uint8('a')'); d.update(uint8('bc')');
    expect('sha256 two updates', lower(sprintf('%02x', typecast(d.digest(), 'uint8'))), ...
           'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad');

    a = double.empty(0, 0);
    b = logical.empty(0, 0);
    expect('double.empty class and size', [strcmp(class(a), 'double') size(a)], [1 0 0]);
    expect('logical.empty class and size', [strcmp(class(b), 'logical') size(b)], [1 0 0]);
    expect('double(x) is the built-in', class(double(int8(3))), 'double');
    expect('logical(x) is the built-in', logical([0 2]), [false true]);
end

function expect(name, value, expected)
    ok = isequaln(value, expected) && strcmp(class(value), class(expected));
    fprintf('  %-28s %s\n', name, mat2str_safe(value));
    if ~ok
        error('test_misc_shims:value', '%s: got %s, expected %s', name, mat2str_safe(value), mat2str_safe(expected));
    end
end

function s = mat2str_safe(v)
    if ischar(v), s = ['"' v '"']; else, s = mat2str(v); end
end
