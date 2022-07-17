function [U,U_initial,t_initial,t_optimization,output] = lscpds_krt(X, N, y, R, varargin)
%LSCPDS_KRT Symmetric CPD solution of a KRT-structured system
%   [U,output] = LSCPDS_KRT(X,N,y,R) computes the factor matrix U{1} and
%   the weight vector U{2} belonging to a rank-R symmetric canonical 
%   polyadic decomposition constrained solution of the row-wise Khatri-Rao
%   structured linear system of equations defined by the generator matrix
%   X, the order N, and the right-hand-side y.
%
%   [U,output] = LSCPDS_KRT(X,N,y,U0) allows the user to provide initial
%   factor matrices U0, to be used by the main algorithm. 
%
%   The structure output contains the output of the selected algorithms:
%
%       output.Initialization   - The output of the initialization step.
%       output.Algorithm        - the output of the main algorithm.
%
%   LSCPDS_KRT(X,N,y,R,options) and LSCPDS_KRT(X,N,y,U0,options) allows the
%   user to choose which algorithms will be used in the different steps and
%   also set parameters for those algorithms:
%
%       Use options.Initialization = [{'auto'}|@lscpds_krt_algebraic|
%       @cpd_rnd] to choose the initialization method. The structure
%       options.InitializationOptions will be passed to the chosen
%       initialization method. On 'auto', @lscpds_krt_algebraic is used 
%       when possible, otherwise @cpd_rnd.
%
%
%       Use options.Algorithm = [{@lscpds_krt_nls}|@lscpds_krt_minf] to 
%       choose the main algorithm. The structure options.AlgorithmOptions 
%       will be passed to the chosen algorithm.
%
%   Further options are:
%       options.Display = false   - Set to true to enable printing output
%                                   information to the command line. If
%                                   options.Display > 1, it is passed on to
%                                   AlgorithmOptions and RefinementOptions,
%                                   unless these structs define Display.
%
%   The following options are passed to AlgorithmOptions and
%   RefinementOptions unless these structs define these options:
%       options.Tolx          - Tolerance for step length
%       options.TolFun        - Tolerance for function value
%       options.MaxIter       - Maximum number of iterations
%       options.CGMaxIter     - Maximum number of CG iterations for inexact
%                                NLS type algorithms.
%
%   Example:
%       I = 10; N = 3; R = 2;
%       U = cpd_rnd([I 1],R);
%       T = cpdgen([repmat(U(1),1,N),U{2}]);% 
%       M = 2*nchoosek(I+N-1,N);
%       X = [ones(M,1), randn(M,I-1)];
%       y = lscpds_krt_eval(U,X,N);
%       [Uest,output] = lscpds_krt(X,N,y,R,'Display',true);
%       cpderr(U,Uest)
%
%   See also ccpd_rnd, ccpdgen.
    
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


    % Check inputs.    
    if size(X,1) ~= length(y)
        error('lscpds_krt:X','size(X,1) should equal size(y,1).');
    end    
    [M,I] = size(X);
        
    % Check R.
    if iscell(R)
        U = R;
        R = size(U{1},2);
        
        if ~iscell(U) || ~ismatrix(U{1}) || ~isvector(U{2})
            error('lscpds_krt:U0','U0 should be a cell with a factor matrix and a vector.')
        end
        if size(U{1},2) ~= length(U{2})
            error('lscpds_krt:U0','size(U0{1},2) should be equal to length(U0{2}).')
        end
        
    end
    
%     % Check if method can be applied.
%     if M < I*(R+1)
%         error('lscpds_krt:X',['The number of equations I=size(X,1)'...
%             'should be at least equal to the number of variables'...
%             'of the symmetric CPD I*(R+1).']);
%     end
    
    % Check options structure
    isfunc = @(f)isa(f,'function_handle');
    xsfunc = @(f)isfunc(f)&&exist(func2str(f),'file');
    p = inputParser;
    p.KeepUnmatched = true;    
    p.addOptional('Initialization','auto');
    p.addOptional('InitializationOptions',struct());
    p.addOptional('Algorithm', @lscpds_krt_nls); %@lscpds_krt_nls, @lscpds_krt_minf
    p.addOptional('AlgorithmOptions',struct());
    p.addOptional('Display', true);
    p.addOptional('TolX', nan);
    p.addOptional('TolFun', nan);
    p.addOptional('TolAbs', nan);    
    p.addOptional('MaxIter', nan);
    p.addOptional('CGMaxIter', nan);    
    p.parse(varargin{:});
    options = p.Results;
    
    fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
    data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
    options = cell2struct(data, fn); 
    
    if ~options.Display, print = @(varargin)true; else print = @fprintf; end

    % STEP 1: initialize the factor matrices unless they were all provided
    % by the user.
    tic;
    print('Step 1: Initialization ');
    if ~exist('U','var'), U = cell(1,N); end
    Uempty = cellfun(@isempty,U);
    if any(Uempty)
        if ischar(options.Initialization) && strcmpi(options.Initialization,'auto')           
            if M >= nchoosek(I+N-1,N)
                options.Initialization = @lscpds_krt_algebraic;
                print('is lscpds_krt_algebraic.\n');
            else
                options.Initialization = @cpd_rnd;
                print('is cpd_rnd (default).\n');
            end
        elseif isfunc(options.Initialization)
            print('is %s (requested by user).\n',func2str(options.Initialization));
        end
        if ~xsfunc(options.Initialization)
            error('lscpds_krt:Initialization','Not a valid initialization.')
        end
    else
        options.Initialization = false;
        print('is manual... \n');
    end

    if isfunc(options.Initialization)        
        % Generate initial factor matrices.
        output.Initialization.Name = func2str(options.Initialization); 
        if strcmpi(output.Initialization.Name,'cpd_rnd')
            [U,output.Initialization] = options.Initialization([I 1],R);
        else        
            [U,~,output.Initialization] = options.Initialization(X,N,y,R,...
                                            options.InitializationOptions);
        end
        % Name gets overwritten
        output.Initialization.Name = func2str(options.Initialization);
        
    else
        output.Initialization.Name = 'manual';        
    end
%     output.Algorithm.relerr = frob(y-lscpds_krt_eval(U,X,N))/frob(y);
%     print('Relative error on y = %.6g.\n',output.Algorithm.fval(end));
%     output.relativeerror_1=output.Algorithm.fval(end);
    U_initial=U;
    toc;
    t_initial=toc;
    % Step 2: run the selected algorithm.
    tic;
    print('Step 2: Algorithm is %s... ',func2str(options.Algorithm));
    if xsfunc(options.Algorithm)
        if ~isfield(options,'AlgorithmOptions')
            options.AlgorithmOptions = struct;
        end
        if ~isfield(options.AlgorithmOptions, 'TolX') && ~isnan(options.TolX)
            options.AlgorithmOptions.TolX = options.TolX;
        end
        if ~isfield(options.AlgorithmOptions, 'TolFun') && ~isnan(options.TolFun)
            options.AlgorithmOptions.TolFun = options.TolFun;
        end
        if ~isfield(options.AlgorithmOptions, 'MaxIter') && ~isnan(options.MaxIter)
            options.AlgorithmOptions.MaxIter = options.MaxIter;
        end
        if ~isfield(options.AlgorithmOptions, 'CGMaxIter') && ~isnan(options.CGMaxIter)
            options.AlgorithmOptions.CGMaxIter = options.CGMaxIter;
        end
        if ~isfield(options.AlgorithmOptions, 'Display') && options.Display > 1
            options.AlgorithmOptions.Display = options.Display;
        end
        
        % Run algorithm.
        [U,output.Algorithm] = options.Algorithm(X,N,y,U,options.AlgorithmOptions);
        output.Algorithm.Name = func2str(options.Algorithm);        
    else
        error('lscpds_krt:Algorithm','Not a valid algorithm.');
    end
    if isfield(output.Algorithm,'iterations')
        print('iterations = %i. ',output.Algorithm.iterations);
    end  
    toc;
    t_optimization=toc;
    % Format output.
    fn = fieldnames(output);
    for f = 1:length(fn)
        output.(fn{f}) = orderfields(output.(fn{f}));
    end
end