function close_log(fid)
%CLOSE_LOG Close one report log opened by mwecmass.output.open_log.
    if ~isempty(fid) && fid > 0
        fclose(fid);
    end
end
