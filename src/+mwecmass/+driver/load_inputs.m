function in = load_inputs()
%LOAD_INPUTS Return the author-defined input struct through the driver namespace.
% Syntax: in = mwecmass.driver.load_inputs(). No inputs; output in is WEC_User_Input()'s struct.
% WEC_User_Input must be on the MATLAB path; this wrapper performs no validation or transformation.

    in = WEC_User_Input();
end
