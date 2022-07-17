classdef LSCPDS_KRT_Kernel 
%LSCPDS_KRT_Kernel Computational routiens for KRT-structured LSCPDS.
%   LSCPDS_KRT_Kernel(X,N,y) is a class collecting computation routines for
%   solving a row-wise Khatri-Rao structured linear system of equations
%   with an Nth-order rank-R symmetric polyadic decomposition constrained 
%   solution, see [1]. The output of the algorithm are the factor matrix 
%   U{1} and the weight vector U{2}.
%
%   LSCPDS_KRT_Kernel is an internal class. lscpds_krt_nls and
%   lscpds_krt_minf are interfaces to the user. Using LSCPDS_KRT_Kernel is
%   only recommended for advanced users who want to have fine-grained
%   control over all settings and functions. To access more advanced
%   information, type
%
%       edit LSCPD_Kernel
%
%   See also lscpds_krt_nls, lscpds_krt_minf.

%   LSCPDS_KRT_Kernel defines the following functions. Let kernel be an
%   instance of the LSCPDS_KRT_Kernel class, constructed for X, N, and y
%   (or alternatively, an initial guess U0):
%
%   The following objective function routines can be used:
%
%       kernel.objfun(z)
%
%   The following gradient functions
%
%       kernel.grad(z)
%
%   The following Gramian functions can be used:
%
%       kernel.JHDJ(z)
%
%   The following Gramian-vector product functions can be used:
%   
%       kernel.JHDJx(z,x)
%
%   The following preconditioners can be used:
%
%       kernel.M_blockJacobi(~, v)
%
%   For large-scale problems, the Gramian-vector products are used, and a
%   preconditioner is used to improve convergence of the CG algorithm. The
%   following options are available:
%
%       usePreconditioner   % Disable preconditioner if false
%
%   See doc LSCPDS_KRT_Kernel for more advanced parameters.
%
%   To use LSCPDS_KRT_Kernel directly, a kernel is initialized, options are 
%   set, initialize is called and an optimization algorithm is run. For 
%   example:
%
%       % Construct kernel with options
%       kernel = LSCPDS_KRT_Kernel(X, N, y, 'usePreconditioner', false);
%       % Initialize kernel
%       kernel.initialize(U0);
%       % Call optimization routine
%       fval    = @(z) kernel.objfun(z);
%       dF.JHF  = @(z) kernel.grad(z);
%       dF.JHJX = @(z) kernel.JHJx(z,x);
%       dF.M    = @(z,b) kernel.M_blockJacobi(z,v);
%       [Ures, output] = nls_gndl(fval, dF, U0);

%   Author(s):  Stijn Hendrikx      (Stijn.Hendrikx@kuleuven.be)
%               Martijn Boussé      (Martijn.Bousse@kuleuven.be)
%               Nico Vervliet       (Nico.Vervliet@kuleuven.be)
%               Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven.be)
%
%   References:
%   [1] Stijn Hendrikx, Martijn Bousse, Nico Vervliet, Lieven De Lathauwer, 
%   "Algebraic and Optimization Based Algorithms for Multivariate 
%   Regression Using Symmetric Tensor Decomposition", 2019 IEEE 8th 
%   International Workshop on Computational Advances in Multi-Sensor 
%   Adaptive Processing (CAMSAP 2019).
%
%   Version History:
%       - 11/12/2020    SH      Initial implementation.
%       - 14/12/2020    MB      Improved functionality. 
    
    properties (SetAccess = protected)
        % Properties of the LSCPDS_KRT_Kernel that cannot be changed after
        % construction, i.e., are read-only. A new kernel should be created
        % for other coefficient matrices, right-hand-sides, or ranks, etc.
        
        % System variables     
        X;                 % matrix that defines coefficient matrix
        y;                 % linear system right hand side
        M;                 % number of equations in system
        
        % Tensor parameters
        R;                 % rank of polyadic decomposition
        N;                 % tensor order
        I;                 % tensor dimension in each mode
    end
    
    properties
        %% Behavioral settings for KRT-structured LSCPDS.
        % Several implementation specific properties can be selected.
        
        useGramian              % Use Gramian information.
        usePreconditioner       % Use preconditioner
                
        %% Advance settings.
        % /
        
        %% Cached variables
        residual            % residual
        J                   % Jacobian
        invJHJ              % (pseudo) inverse of Gramian of Jacobian
        offset              % offsets for vectorized variables
        XTB                 % inner products of X and U{1}
        XTB_N               % inner products of X and U{1} to the power N
        XTB_N_1             % inner products of X and U{1} to the power N-1
        
        isInitialized       % true if initialize has been run.
                
    end
    
    methods
        
        % Constructor
        function this = LSCPDS_KRT_Kernel(X, N, y, varargin)
            % Matrix X determines the coefficient matrix A = krt(X, ..., X)
            % , with N repetitions of X, of the structured linear system.
            % y contains the right hand side of this system. If the system
            % contains M equations, then X is of dimensions M x I and y of
            % dimensions M x 1.
            
            % Check arguments
            % X
            if ~isnumeric(X) || ~ismatrix(X)
                error('LSCPDS_KRT_Kernel:X','X should be a matrix.')
            end
            % N
            if ~isscalar(N) || N ~= abs(round(N)) || N == 0
                error('LSCPDS_KRT_Kernel:N','N should be a positive integer, excluding zero.')
            end
            % y
            if ~isnumeric(y) || ~isvector(y)
                error('LSCPDS_KRT_Kernel:y','y should be a vector.')
            end
            if size(X,1) ~= length(y)
                error('LSCPDS_KRT_Kernel:y','length(y) should be equal to size(X,2).')
            end
   
            % Variables
            this.X = X; 
            this.y = y;
            this.N = N; 
            [this.M, this.I] = size(this.X);
            
            % Parse varargin
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('UseGramian', false);
            p.addOptional('UsePreconditioner', false);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            options = p.Results;
            
            % Set behavior
            this.useGramian             = options.UseGramian;
            this.usePreconditioner      = options.UsePreconditioner;
            
            this.isInitialized = false;
        end
        
        
        function this = initialize(this, z)
        %INITIALIZE Initialize kernel before computations.
        %   INITIALIZE(z) initializes the kernel before iterations start
        %   and should be called before calling the optimization routine.
        %   Results independent of the variables are precomputed here. z is
        %   a cell containing initial factor matrices and is used to derive
        %   properties of the data, e.g., wheter the factors are complex or
        %   have higher order than the tensor. INITIALIZE should be called
        %   every time a setting is altered to ensure a correct state.

            this.offset = cumsum([1 cellfun(@numel,z)]);
            this.R = size(z{1},2);
            
            this.isInitialized = true;            
        end

        function this = state(this, z)
        %STATE Update the state of the kernel.
        %   STATE(z) updates the state for a regular iteration of the
        %   algorithm, which includes caching results that change every
        %   iteration but are used by multiple functions.
            this = this.initialize(z);
            % Cache some variables.
            this.XTB = this.X*z{1};
            this.XTB_N = this.XTB.^this.N;
            this.XTB_N_1 = this.XTB.^(this.N-1);
            
            % Compute sub-Jacobians.
            this.J{1} = krt(z{2} .* (this.N*this.XTB_N_1), this.X);
            %this.J{1} = kr((z{2} .* (this.N*this.XTB_N_1)),this.X(:,1:size( (z{2} .* (this.N*this.XTB_N_1)),2)));
            %this.J{1} = (z{2} .* (this.N*this.XTB_N_1)).* this.X(:,1:size( (z{2} .* (this.N*this.XTB_N_1)),2));
            this.J{2} = this.XTB_N;         
            
            % Cache variables for preconditioner vector product
            if this.usePreconditioner && ~this.useGramian
                this.invJHJ{1} = pinv(this.J{1}.'*conj(this.J{1}));
                this.invJHJ{2} = pinv(this.J{2}.'*conj(this.J{2}));
            end
        end  
        
        
        function fval = objfun(this, z)
        %OBJFUN Objective function value for KRT-structured LSCPDS.  
            this = this.state(z);       
            this.residual = this.XTB_N*z{2}.' - this.y;            
            fval = 0.5*this.residual'*this.residual; % conj
        end
                
        function grad = grad(this, z) 
        %GRAD Gradient for KRT-structured LSCPDS.    
            this = this.state(z);       
            this.residual = this.XTB_N*z{2}.' - this.y;           
            grad = [conj(this.J{1}), conj(this.J{2})].'*this.residual;            
        end       
        
        % Multiplies preconditioner matrix with vector x.
        function y = M_blockJacobi(this, ~, x)         
        %M_BLOCKJACOBI Block-Jacobi preconditioner for KRT-structured LSCPDS.                          
            this = this.state(z);
            y = nan(size(x));
            idxb = this.offset(1):this.offset(2)-1;
            idxc = this.offset(2):this.offset(3)-1;
            y(idxb) = this.invJHJ{1}*x(idxb);
            y(idxc) = this.invJHJ{2}*x(idxc);
        end
                                     
        function gram = JHDJ(this, z)
        %JHDJ Gramian of KRT-structured LSCPDS.
            this = this.state(z);
            jac = [this.J{1}, this.J{2}];
            gram = jac'*jac;
        end
                
        function y = JHDJx(this, z, x)
        %JHDJx Gramian-vector product for KRT-structured LSCPDS.
            this = this.state(z);
            XB = reshape(x(1:this.I*this.R), [this.I this.R]);
            xc = x(this.I*this.R+1:end);
            
            % JB*XB: krt product times vectorized matrix
            JBXB = -this.N*(((this.X*XB) .* this.XTB_N_1)*z{2}.');
            % Add Jc*xc
            Jx = JBXB - this.XTB_N*xc;
            % JB^T*(J*x): kr product times vector --> vectorized CPD
            JBTJx = -this.N*cpdgen({this.X', z{2}' .* this.XTB_N_1', Jx.'});
            % Add Jc^T*(J*x)
            y = [JBTJx(:); -this.XTB_N' * Jx(:)];
        end
        
        function y = JHDJx_jac(this, z, x)
        %JHDJx_jac Gramian-vector product for KRT-structured LSCPDS using
        % explicit Jacobians.
            this = this.state(z);
             
            jac = [this.J{1}, this.J{2}];
            y = jac'*(jac*x); 
        end
       
        % Validate
        function isvalid = validate(this, z)
        %VALIDATE Check whether this kernel is in a valid state.
        %   VALIDATE(z) checks whether the combination of data and initial 
        %   variables z is ivalid. This is a deep check, and may take
        %   some time, but it ensures that the computation will succeed.
        
            if ~iscell(z) || length(z) ~= 2 || ~isnumeric(z{1}) || ...
               ~isnumeric(z{2}) || ~ismatrix(z{1}) || ~isvector(z{2}) 
                error('LSCPDS_KRT_Kernel:U0','z should be a cell containing a factor matrix and a vector.');
            end
            if size(z{1},2) ~= length(z{2})
                error('LSCPDS_KRT_Kernel:U0','size(z{1},2) should be equal to length(z{2}).')
            end
            if size(z{1},1) ~= this.I
                error('LSCPDS_KRT_Kernel:U0','size(z{1},1) should be equal to size(X,2).')
            end
            
            % Check if method can be applied.
            if size(this.X,1) < size(this.X,2)*(this.R+1)
                error('LSCPDS_KRT_Kernel:X','The number of equations I = size(X,1) should be at least I*(R+1) with R = size(U0{1},2).');
            end
            
            isvalid = true;            
        end        
    end
  
end