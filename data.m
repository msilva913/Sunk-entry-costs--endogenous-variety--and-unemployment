%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE_DUMMY_DATA - Generates placeholder data for Dynare estimation
%
% Purpose:
%   Dynare requires a data file to begin estimation, even when using
%   impulse response matching (IRF) estimation methods. This script creates
%   a minimal dummy dataset to satisfy Dynare's file requirement without
%   affecting the IRF matching estimation.
%
% Output:
%   Saves dummy_data.mat containing random data that will be ignored during
%   IRF matching estimation.
%
% Note:
%   The actual estimation uses impulse response matching, so this data is
%   purely procedural and doesn't affect results.
%
% During IRF matching estimation, Dynare will:
%   1. Require this file to start estimation
%   2. Ignore the actual data values
%   3. Use only the specified impulse responses for estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create minimal dataset that satisfies Dynare's file format requirements
T = 10;          % Number of fake time periods (arbitrary small number)
Y = rand(T, 1);  % Generate random data column vector (values unimportant)

% Save in Dynare-compatible format
save('dummy_data.mat', 'Y');  % Single variable MAT-file for Dynare