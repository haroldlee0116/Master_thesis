function y = lscpds_krt_eval(U, X, N)
%LSCPDS_KRT_EVAL Compute RHS of a KRT- and symmetric CPD structured system
% lscpds_krt_eval(U, X, N) computes the right hand side (RHS) y of a linear 
% system Ax = y with a Khatri-Rao transposed (KRT) structured coefficient 
% matrix A = krt(X, ..., X) and a symmetric Nth order rank R CPD 
% solution x = vec([|B, ..., B, c|]). A and x respectively contain N
% repetitions of X and B. B and c can be specified using a cell array
% U = {B, c} or a column vector U = [vec(B); c^T].

%   Authors:    Stijn Hendrikx          (Stijn.Hendrikx@kuleuven.be)
%               Martijn Bousse          (Martijn.Bousse@kuleuven.be)
%               Nico Vervliet           (Nico.Vervliet@kuleuven.be)              
%               Lieven De Lathauwer     (Lieven.DeLathauwer@kuleuven.be) 
%
%   References:
%   [1] Stijn Hendrikx, Martijn Bousse, Nico Vervliet, Lieven De Lathauwer, 
%   "Algebraic and Optimization Based Algorithms for Multivariate 
%   Regression Using Symmetric Tensor Decomposition", 2019 IEEE 8th 
%   International Workshop on Computational Advances in Multi-Sensor 
%   Adaptive Processing (CAMSAP 2019).
%
%   Version History:
%       - 11/12/2020   SH      Initial version.

    % Check inputs.
    if ~iscell(U)
        I = size(X, 2);
        R = length(U) / (I + 1);
        if R ~= round(R)
            error('lscpds_krt_eval:U', 'If U is a vector, it should be of length (I+1)R.');
        end
        B = reshape(U(1:I*R), [I R]);
        c = reshape(U(I*R + 1:end), [1 R]);
    elseif length(U) ~= 2
        error('lscpds_krt_eval:U', 'If U is a cell, it should contain two elements.');
    else
        B = U{1}; c = U{2};
    end
    
    N = round(N);    

    % Compute right-hand-side.
    y = (((X*B).^N)*c.');
    
end