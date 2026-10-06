function emit(fids, fmt, varargin)
%EMIT Write one formatted report message to every requested file identifier.
    for k = 1:numel(fids)
        fprintf(fids(k), fmt, varargin{:});
    end
end
