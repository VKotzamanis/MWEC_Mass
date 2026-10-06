function save_hydro_cache(path, hydro_table)
%SAVE_HYDRO_CACHE Write a hydro_table struct to a .mat cache file.
    save(path, 'hydro_table', '-v7.3');
end
