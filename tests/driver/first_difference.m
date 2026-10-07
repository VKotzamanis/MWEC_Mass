function where = first_difference(a, b, path)
%FIRST_DIFFERENCE Path of the first difference between two configs, '' when they are identical.
%   Exact comparison (NaN equals NaN, class and size must match). A containers.Map is compared by
%   keys and values, an MS2Parser by its public properties. Test helper, not a test.
  if nargin < 3, path = 'config'; end
  where = '';
  if isa(a, 'containers.Map') || isa(b, 'containers.Map')
    if ~(isa(a, 'containers.Map') && isa(b, 'containers.Map')) || ~isequal(keys(a), keys(b))
      where = [path '.keys'];
      return;
    end
    names = keys(a);
    for k = 1:numel(names)
      where = first_difference(a(names{k}), b(names{k}), [path '{' names{k} '}']);
      if ~isempty(where), return; end
    end
  elseif isobject(a) || isobject(b)
    if ~strcmp(class(a), class(b))
      where = [path ':class'];
      return;
    end
    props = {'entities', 'visible_surfs', 'units', 'extents', 'filename', 'file_symmetry'};
    for k = 1:numel(props)
      where = first_difference(a.(props{k}), b.(props{k}), [path '.' props{k}]);
      if ~isempty(where), return; end
    end
  elseif isstruct(a)
    if ~isstruct(b) || ~isequal(size(a), size(b)) || ~isequal(sort(fieldnames(a)), sort(fieldnames(b)))
      where = [path ':struct'];
      return;
    end
    names = fieldnames(a);
    for e = 1:numel(a)
      for k = 1:numel(names)
        where = first_difference(a(e).(names{k}), b(e).(names{k}), [path '.' names{k}]);
        if ~isempty(where), return; end
      end
    end
  elseif iscell(a)
    if ~iscell(b) || ~isequal(size(a), size(b))
      where = [path ':cell'];
      return;
    end
    for k = 1:numel(a)
      where = first_difference(a{k}, b{k}, sprintf('%s{%d}', path, k));
      if ~isempty(where), return; end
    end
  elseif ~(isequaln(a, b) && strcmp(class(a), class(b)))
    where = path;
  end
end
