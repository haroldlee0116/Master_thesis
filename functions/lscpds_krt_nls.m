function [U,output] = lscpds_krt_nls(X,N,y,U0,varargin)
%LSCPDS_KRT_NLS KRT-structured LSCPDS by nonlinear least squares.
%   U = LSCPDS_KRT_NLS(X,N,y,U0) computes the factor matrix U{1} and the
%   weight vector U{2} belonging to the symmetric canonical polyadic
%   decomposition constrained solution of the row-wise Khatri-Rao
%   structured linear system of equations defined by the generator matrix
%   X, the order N, and the right-hand-side y, by minimizing
%   0.5*frob(krt(X,X,...,X)*tens2vec(cpdgen({U0{1},...,U0{1},U0{2}}))-y)^2.
%
%   [U,output] = LSCPDS_KRT_NLS(X,N,y,U0) additionally returns the
%   structure output containing the following information:
%
%       output.Name  - The name of the selected algorithm.
%       output.<...> - The output of the selected algorithm.
%
%   LSCPDS_KRT_NLS(X,N,y,U0,options), in which options is a struct, or
%   LSCPDS_KRT_NLS(X,N,y,U0,'key',value) may be used to set the following
%   options:
%
%      options.Algorithm =         - The desired optimization method.
%      [@nls_gncgs| ...            
%       {@nls_gndl}|@nls_lm]       
%      options.LargeScale          - If true, the Gauss-Newton or Levenberg-
%      = [true|false|{'auto'}]       Marquardt steps are computed using a
%                                    preconditioned conjugate gradient algorithm.
%                                    Otherwise, a direct solver is used. When
%                                    'auto' is selected, the problem is
%                                    large-scale when the number of variables
%                                    exceeds 100, i.e., sum(cellfun(@numel,U0))
%                                    > 100.
%      options.M =                 - The preconditioner to use when
%      [{'block-Jacobi'}|...         options.LargeScale is true.
%       false]                     
%      options.<...>               - Parameters passed to the selected method,
%                                    e.g., options.TolFun, options.TolX and
%                                    options.CGMaxIter. See also help
%                                    [options.Algorithm].
%
%   See also lscpds_krt_minf.

%   Authors:    Martijn Bousse          (Martijn.Bousse@kuleuven.be)
%               Stijn Hendrikx          (Stijn.Hendrikx@kuleuven.be)
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
%       - 14/12/2020   MB      Initial version.

    % Check options structure.
    p = inputParser;
    p.addOptional('Algorithm', @nls_gndl);
    p.addOptional('MaxIter', 200);
    p.addOptional('TolFun', 1e-12);
    p.addOptional('TolX', 1e-8);
    p.addOptional('M', 'block-Jacobi');
    p.addOptional('CGMaxIter', 25);
    p.addOptional('LargeScale', false);
    p.addOptional('CheckInputs', true);
    p.addOptional('Display', false);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;
    
    fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
    data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
    options = cell2struct(data, fn);
    
    % Large-scale test
    if ischar(options.LargeScale) || isnan(options.LargeScale) 
        options.LargeScale = sum(cellfun(@numel, U0)) > 100;
    end 
    
    % Construct kernel options
    if ~isfield(options, 'UsePreconditioner')
        options.UsePreconditioner = options.LargeScale && ...
            (~islogical(options.M) || options.M);
    end

    % Construct kernel and check for errors.        
    newoptions = [fieldnames(options)'; struct2cell(options)'];
    kernel = LSCPDS_KRT_Kernel(X, N, y, newoptions{:});
    if options.CheckInputs, kernel.validate(U0); end
    
    % Initialize kernel.
    kernel.initialize(U0);

    % Construct optimization object.
    fval = @kernel.objfun;
    dF = struct;
    dF.JHF = @kernel.grad;
    if options.LargeScale
        dF.JHJx = @kernel.JHDJx;
        if kernel.usePreconditioner
            dF.M = @kernel.M_blockJacobi;
        end
    else
        dF.JHJ = @kernel.JHDJ;
    end
    
    % Call optimization routine.
    [U,output] = options.Algorithm(fval,dF,U0(:).',newoptions{:});    
    output.kernel = kernel;
end