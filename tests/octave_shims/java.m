function obj = java()
%JAVA shim for the only Java use in the pipeline: java.security.MessageDigest SHA-256 (export_results.m).
%   Octave 8.4 here has no Java. Supports java.security.MessageDigest.getInstance('SHA-256'),
%   then d.update(uint8 bytes) any number of times and d.digest(), which returns the 32 digest
%   bytes as int8 like Java's byte[]. The hash is Octave's hash(). Nothing else of Java exists.
  obj = struct('security', struct('MessageDigest', struct('getInstance', @get_instance)));
end

function d = get_instance(algorithm)
  if ~strcmpi(algorithm, 'SHA-256')
    error('java:unsupported', 'java shim supports only the SHA-256 MessageDigest.');
  end
  state = containers.Map({'bytes'}, {uint8([])});
  d = struct('update', @(b) update_bytes(state, b), 'digest', @() finish(state));
end

function update_bytes(state, b)
  state('bytes') = [state('bytes'); uint8(b(:))];
end

function out = finish(state)
  hex = hash('sha256', char(state('bytes')'));
  out = typecast(uint8(hex2dec(reshape(hex, 2, [])'))', 'int8');
end
