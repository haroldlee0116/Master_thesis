function [U,T,output] = lscpds_krt_algebraic(X,N,y,R,varargin)
%LSCPDS_KRT_ALGEBRAIC Algebraic algorithm for KRT-structured LSCPDS.
%   [T,output] = LSCPDS_KRT_ALGEBRAIC(X,N,y,R) computes the factor matrix
%   U{1} and the weight vector U{2} belonging to the rank-R symmetric 
%   canonical polyadic decomposition of T which is the solution to row-wise
%   Khatri-Rao structured linear system of equations defined by the
%   generator matrix X, the order N, and the right-hand-side y. Note that T
%   is obtained algebraically while U is computed via the cpd routine which
%   also entails optimization.
%
%   LSCPDS_KRT_ALGEBRAIC(X,n,b,options) allows the user to set following
%   options:
%
%       options.Solver =        - Select the solver for the internal linear
%       ['qr',{'chol'}]           system (either 'qr' or 'chol'). If not
%                                 given, 'chol' is used.
%   
%   Example:
%     I = 10; N = 3; R = 3;
%     U = cpd_rnd([I 1],R);
%     T = cpdgen({U{1},U{1},U{1},U{2}});
% 
%     M = 500;
%     X = [ones(M,1), randn(M,I-1)];
%     y = lscpds_krt_eval(U,X,N);
% 
%     [Uest,Test] = lscpds_krt_algebraic(X,N,y,R);
%     disp(cpderr(U,Uest))
%     disp(frob(T-Test)/frob(T))
%
%   See also lscpds_krt_nls, lscpds_krt_minf, lscpds_krt.

%   Author(s):  Nico Vervliet           (Nico.Vervliet@kuleuven.be)
%               Martijn Bousse          (Martijn.Bousse@kuleuven.be)
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
%       - 2019/07/18   NV      Initial version.
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('Solver', 'chol');
    p.addOptional('SymmetryTolerance', 1e-10, @isscalar);
    p.addOptional('UseSDF', false);
    p.addOptional('CheckInputs', true);
    p.KeepUnmatched = true;
    p.parse(varargin{:});

    fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
    data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
    options = cell2struct(data, fn);
    
    % Check inputs.
    if options.CheckInputs
        % Check arguments.
        % X
        if ~isnumeric(X) || ~ismatrix(X)
            error('lscpds_krt_algebraic:X','X should be a matrix.')
        end
        % N
        if ~isscalar(N) || N ~= abs(round(N)) || N == 0
            error('lscpds_krt_algebraic:N','N should be a positive integer, excluding zero.')
        end
        % R
        if ~isscalar(R) || R ~= abs(round(R)) || R == 0
            error('lscpds_krt_algebraic:R','R should be a positive integer, excluding zero.')
        end
        % y
        if ~isnumeric(y) || ~isvector(y)
            error('lscpds_krt_algebraic:y','y should be a vector.')
        end
        if size(X,1) ~= length(y)
            error('lscpds_krt_algebraic:y','length(y) should be equal to size(X,2).')
        end  
        % Check if method can be applied.
        [M,I] = size(X);
        if M < nchoosek(I+N-1,N)
            error('lscpds_krt_algebraic:X',['The number of equations' ...
                'I=size(X,1) should be at least nchoosek(I+N-1,N).']);
        end 
    end
    
    % Parameters
    I = size(X, 2);

    % Construct indices of unique columns and how often they appear
    ind = repmat({1:I},1,N);
    [ind{:}] = ndgrid(ind{:});
    ind = cellfun(@(i) i(:), ind, 'UniformOutput', false);
    ind = horzcat(ind{:});
    ind = sort(ind, 2); % sort to find unique entries
    ind = sum((ind - 1) .* (I.^(N-1:-1:0)),2) + 1; % sub2ind
    [~,ib,ic] = unique(ind);
    
    % Subscripts of unique columns of A. 
    sub = cell(1, N);
    [sub{:}] = ind2sub(ones(1,N)*I, ib);
    
    % Construct only unique columns, i.e., compress A.
    K = X(:,sub{1});
    for n = 2:N
        K = K .* X(:,sub{n});
    end 

    % Solve compressed linear system
    if nargin < 4 || strcmpi(options.Solver,'chol')
        C = chol(K'*K);
        T = C\(C'\(K'*y));
    elseif strcmpi(options.Solver,'qr')
        [Q,R] = qr(K,0);
        T = R\(Q'*y);
    else
        warning('lscpds_krt_algebraic:Solver', ...
            'Unknown solver specified, falling back on default (''chol'').');
    end
    
    % Undo compression and return tensor
    cnt = histc(ic,1:numel(ib));
    T = reshape(1./cnt(ic) .* T(ic), ones(1,N)*I);
    
    % Compute symmetric CPD
    if N == 2
        [U,s,V] = svds(T,R);
        s = sqrt(diag(s)).';
        Ucpd = {U.*s, V.*s};
    else
        [Ucpd,output] = cpd(T,R);
    end
    
    % Check symmetry
    Ucpd1 = repmat(Ucpd(1), 1, N-1);
    UcpdN_1 = Ucpd(2:end);
    % Symmetric up to flipped signs of vectors in factor matrices.
    if ~options.UseSDF || ...
       (max(cpderr(Ucpd1, UcpdN_1)) < options.SymmetryTolerance) 
        % Account for flipped signs.
        pos_signs = sign(Ucpd{1}(1,:));
        signs = ones(N-1,R);
        for n = 2:N
            ind = pos_signs ~= sign(Ucpd{n}(1,:));
            signs(n-1,ind) = -1;
        end
        flip_ind = prod(signs,1) == -1;

        % Determine B and c.
        coeff = vecnorm(Ucpd{1});
        c = coeff.^N;    
        B = Ucpd{1}./coeff; % Normalize columns.
        c(flip_ind) = -c(flip_ind); % Flip signs.
        U = {B, c};
    else    % Not symmetric, enforce symmetry using sdf.
        model = struct;
        model.variables = cpd_rnd([I 1], R); 
        model.factors = 1:2;
        model.factorizations{1}.data = T;
        model.factorizations{1}.cpd  = [ones(1,N) 2];
        model.factorizations{1}.issymmetric = true;

        U = ccpd_nls(model);
    end
        
end