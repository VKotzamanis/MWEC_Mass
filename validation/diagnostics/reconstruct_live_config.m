function cfg = reconstruct_live_config(cfg)
%RECONSTRUCT_LIVE_CONFIG Restore the two class-instance fields the export strips from a config.
% Inputs: cfg (results.config from a .mat file or a live config). Reconstructs cfg.ms2_model
% from cfg.ms2_file and verifies it using SHA-256 digest to catch deck edits since the run.
% Rewraps cfg.boundary_cache.data from exported struct form to containers.Map. Passes through
% unchanged if cfg already carries live ms2_model. Outputs: cfg with restored class objects ready for re-solve or geometry re-draw.
    if ~isfield(cfg, 'ms2_model') || isempty(cfg.ms2_model)
        if ~isfile(cfg.ms2_file)
            error('reconstruct_live_config:MissingDeck', 'Hull-deck file not found: %s', cfg.ms2_file);
        end
        digest = local_sha256(cfg.ms2_file);
        if isfield(cfg, 'ms2_deck_sha256') && ~isempty(cfg.ms2_deck_sha256) ...
                && ~strcmpi(digest, cfg.ms2_deck_sha256)
            warning('reconstruct_live_config:DeckDigestMismatch', ...
                ['%s has changed since the run: SHA-256 %s now vs %s recorded in the result ' ...
                 'file. This reads the deck as it stands today, not the deck the run actually ' ...
                 'used.'], cfg.ms2_file, digest, cfg.ms2_deck_sha256);
        end
        cfg.ms2_model = mwecmass.geometry.MS2Parser.parse(cfg.ms2_file);
    end

    if isfield(cfg, 'boundary_cache') && isstruct(cfg.boundary_cache) && isscalar(cfg.boundary_cache) ...
            && isfield(cfg.boundary_cache, 'data') && ~isa(cfg.boundary_cache.data, 'containers.Map')
        bc_data = cfg.boundary_cache.data;
        if isstruct(bc_data) && isscalar(bc_data) && isequal(fieldnames(bc_data), {'keys'; 'values'})
            % the fallback export form: a surface name that was not a valid MATLAB identifier
            cfg.boundary_cache.data = containers.Map(bc_data.keys{1}, bc_data.values{1});
        else
            % the common export form: surface names as the struct's own field names
            surface_names = fieldnames(bc_data);
            surface_data  = struct2cell(bc_data);
            cfg.boundary_cache.data = containers.Map(surface_names, surface_data);
        end
    end
end

function hex = local_sha256(file)
%LOCAL_SHA256 SHA-256 hex digest of a file's bytes, matching mwecmass.output.export_results's own digest.
% Reads raw bytes with fread (*uint8) so the digest is of the file as stored, independent of text
% encoding or line endings -- the same method export_results.m used to compute
% config.ms2_deck_sha256 at export time, so the two digests are comparable byte for byte.
    fid = fopen(file, 'r');
    assert(fid >= 0, 'reconstruct_live_config:local_sha256:CannotOpen', 'Cannot open %s.', file);
    cleaner = onCleanup(@() fclose(fid));
    bytes = fread(fid, Inf, '*uint8');
    digest = java.security.MessageDigest.getInstance('SHA-256');
    if ~isempty(bytes)
        digest.update(bytes);
    end
    hex = lower(sprintf('%02x', typecast(digest.digest(), 'uint8')));
end
