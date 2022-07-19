%%  Experiment 1:
%   In this tutorial we illustrate how to use the KRT-structured LSCPDS [1]
%   to fit the sythenic MISO Volterra datasets [2] in this 
%   framework.

%   Authors:    Zhehan Li                  (z.li-50@student.tudelft.nl) 
%               Kim Batselier [2]          (K.Batselier@tudelft.nl) 
%               Martijn Bousse [1]         (Martijn.Bousse@kuleuven.be)
%               Stijn Hendrikx [1]         (Stijn.Hendrikx@kuleuven.be)
%               Nico Vervliet [1]          (Nico.Vervliet@kuleuven.be)              
%               Lieven De Lathauwer [1]    (Lieven.DeLathauwer@kuleuven.be) 
%
%   References:
%   [1] Stijn Hendrikx, Martijn Bousse, Nico Vervliet, Lieven De Lathauwer, 
%   "Algebraic and Optimization Based Algorithms for Multivariate 
%   Regression Using Symmetric Tensor Decomposition", 2019 IEEE 8th 
%   International Workshop on Computational Advances in Multi-Sensor 
%   Adaptive Processing (CAMSAP 2019).
%   [2] Kim Batselier. “Enforcing symmetry in tensor network MIMO Volterra 
%   identification”. In: IFAC-PapersOnLine 54.7 (2021), pp. 469–474.


close all; clear variables; clc;
%  
%% In order to generate a random KRT-structured LSCPDS problem, we first
% define the size, order, and rank of the symmetric CPD:
p = 2;
M = 3;
I = p*M+1; 
l = 1;
d = 3;
R = 7;
N = 50*nchoosek(I+d-1,I-1);
%% Check whether we are working on under-determined case or determined case
N_samples = N*l;
N_free = (I+1)*R*l;
disp('Step 1: Ground Truth System Parameters: ');
fprintf('First we build the ground truth system with p = %g inputs and l = %g outputs. \n',p,l);
fprintf('This MIMO system has truncated memory M = %g, system order d = %g, and the rank of CPD R = %g. \n',M,d,R);
if N_samples < N_free;
    disp('Under-determined case: num(Y) < ( num(B) + num(c) )')
else N_samples >= N_free;
    disp('Determined case: num(Y) >= ( num(B) + num(c) )')
end
%% Next, we generate a random factor matrix B and a weight vector c as well 
% as the corresponding symmetric tensor.
for i=1:R
    B(:,i)=abs(1e0*randn(1))*exp(-randi(10,1).*[1:I]);
end
%B = rand(I,R)*sin(I);
c = rand(1,R);
V = {B,c};
T = cpdgen([repmat(V(1),1,d),V{2}]);
issym = norm(symcheck(reshape(T,I*ones(1,d))));
strcat("Symmetry coefficient of the ground truth tensor: ",num2str(issym))
% Now we define the linear system of equations by defining the number of
% equations, the generator for the coefficient matrix, and the right-hand
% side.
u = rand(N,p);
for t=1:N
        for k=1:M
            k=k-1;
           if t-k<1
                U(t,p*(k+1):p*(k+2)-1)=1;
            else
                U(t,p*(k+1):p*(k+2)-1)=u(t-k,1:p);
            end
        end
end
U(1:N,1)=ones(N,1);

y = lscpds_krt_eval(V,U,d);
% Compensation rates
alpha1 = (0.5*N*I*(nchoosek(I+d-1,I-1)^2)+(7/6)*(nchoosek(I+d-1,I-1)^3)+N*l*nchoosek(I+d-1,I-1))/(15*(N*(R^2)*(I^2)+(R^3)*(I^3)));
alpha2 = (nchoosek(I+d-1,I-1)*(I^d))/(R*(I+1)*(I+3));
strcat("The compensation rate alpha1= ",num2str(alpha1))
strcat("The compensation rate alpha2= ",num2str(alpha2))
%%  Algebraic initialization
% For the computation of the KRT-structured LSCPDS problem, you can use the
% lscpds_krt routine which entails two steps: (algebraic) initialization
% and optimization.
disp('--------------------------New Experiment-------------------------------');
disp('--------------------------New Experiment-------------------------------');
disp('First problem:');
options.Initialization = @lscpds_krt_algebraic;
[Vest1,Vest1_initial,t_initial1,t_optimization1,recording1] = lscpds_krt(U,d,y,R,'Display',true,options);
disp('--------------------------Result Discussion--------------------------------');
disp('Totoal computational time spent in identification (s) ');
disp(t_optimization1+t_initial1);
e_y_ini=frob(y-lscpds_krt_eval(Vest1_initial,U,d))/frob(y);
disp('S1: Relative error on y after algebraic initialization');
%disp(recording1.Algorithm.fval(1));
disp(e_y_ini);
Tini = cpdgen([repmat(Vest1_initial(1),1,d),Vest1_initial{2}]);
disp('    R2 (relative error on tensor) after algebraic initialization');
disp(frob(T-Tini)/frob(T))
disp('    Computational time spent in initialization (s)');
disp(t_initial1);
e_y_est=frob(y-lscpds_krt_eval(Vest1,U,d))/frob(y);
disp('S2: Relative error on y after optimization');
%disp(recording1.Algorithm.fval(end));
disp(e_y_est);
Test = cpdgen([repmat(Vest1(1),1,d),Vest1{2}]);
disp('    R2 (relative error on tensor) after optimization');
disp(frob(T-Test)/frob(T))
disp('    Computational time spent in optimization (s) ');
disp(t_optimization1);
issym = norm(symcheck(reshape(Test,I*ones(1,d))));
disp('S3: Symmetry coefficient of identified tensor obtained via LS-CPD');
disp(issym);
U_val=randn(N*l,I)*sin(N);
y_val=lscpds_krt_eval(V,U_val,d);
e_valy=frob(y_val-lscpds_krt_eval(Vest1,U_val,d))/frob(y_val);
disp('S3: Validation error on y');
disp(e_valy);
%%  Random initialization
% You can inspect the relative error between the true and estimated factor
% matrix after optimal scaling and permutation using cpderr.
disp('--------------------------New Experiment-------------------------------');
disp('--------------------------New Experiment-------------------------------');
disp('Random initialization:');
options.Initialization = @cpd_rnd;
[Vest2,Vest2_initial,t_initial2,t_optimization2,recording2] = lscpds_krt(U,d,y,R,'Display',true,options);
disp('--------------------------Result Discussion--------------------------------');
disp('Totoal computational time spent in identification (s) ');
disp(t_optimization2+t_initial2);
disp('S1: Relative error on y after random initialization');
e_y_ini=frob(y-lscpds_krt_eval(Vest2_initial,U,d))/frob(y);
disp(e_y_ini);
Tini = cpdgen([repmat(Vest2_initial(1),1,d),Vest2_initial{2}]);
disp('    R2 (relative error on tensor) after random initialization');
disp(frob(T-Tini)/frob(T))
disp('    Computational time spent in initialization (s)');
disp(t_initial2);

e_y_est=frob(y-lscpds_krt_eval(Vest2,U,d))/frob(y);
disp('S2: Relative error on y after optimization');
%disp(recording2.Algorithm.fval(end));
disp(e_y_est)
Test = cpdgen([repmat(Vest2(1),1,d),Vest2{2}]);
disp('    R2 (relative error on tensor) after optimization');
disp(frob(T-Test)/frob(T))
disp('    Computational time spent in optimization (s) ');
disp(t_optimization2);
T = cpdgen([repmat(Vest2(1),1,d),Vest2{2}]);
issym = norm(symcheck(reshape(T,I*ones(1,d))));
disp('S3: Symmetry coefficient of solution obtained via LS-CPD');
disp(issym);
U_val=randn(N*l,I)*sin(N);
y_val=lscpds_krt_eval(V,U_val,d);
e_valy=frob(y_val-lscpds_krt_eval(Vest2,U_val,d))/frob(y_val);
disp('S3: Validation error on y');
disp(e_valy);
%% Add noise on the output samples
disp('--------------------------New Experiment-------------------------------');
disp('--------------------------New Experiment-------------------------------');
disp('Same problem with noise');
disp('Random initialization');
options = struct;
options.Display = false;
tic;
B0 = cpd_rnd([I 1],R);
toc;
t_initial = toc;
tic;
[Best_nls,out_nls] = lscpds_krt_nls(U,d,noisy(y,20),B0,options);
toc;
t_optimization = toc;
disp('--------------------------Result Discussion--------------------------------');
disp('Totoal computational time spent identification (s) ');
disp(t_initial+t_optimization)
disp('S1: Relative error on y after initialization');
disp(out_nls.fval(:,1))
disp('    Computational time spent in initialization (s)');
disp(t_initial);
disp('S2: Relative error on y after optimization');
disp(out_nls.fval(end))
disp('    Computational time spent in optimization (s)');
disp(t_optimization);
disp('    Number of CG iterations:');
disp(out_nls.iterations)
T = cpdgen([repmat(Best_nls(1),1,d),Best_nls{2}]);
issym = norm(symcheck(reshape(T,I*ones(1,d))));
disp('S3: Symmetry coefficient of solution obtained via LS-CPD');
disp(issym);
% If we could use a better initialization (than random), we might obtain
% the solution in fewer iterations. For this purpose, we have also 
% implemented an algebraic method, which can be called using 
% lscpds_krt_algebraic. We observe fewer iterations for NLS starting from 
% this solution.
disp('--------------------------New Experiment-------------------------------');
disp('--------------------------New Experiment-------------------------------');
disp('Same problem with noise');
disp('Algebraic initialization');
tic;
B0 = lscpds_krt_algebraic(U,d,y,R);
toc;
t_initial = toc;
tic;
[Best_nls,out_nls] = lscpds_krt_nls(U,d,noisy(y,20),B0,options);
toc;
t_optimization=toc;
disp('--------------------------Result Discussion--------------------------------');
disp('Totoal computational time spent identification (s) ');
disp(t_initial+t_optimization)
disp('S1: Relative error on y after initialization');
disp(out_nls.fval(:,1))
disp('    Computational time spent in initialization (s)');
disp(t_initial);
disp('S2: Relative error on y after optimization');
disp(out_nls.fval(end))
disp('    Computational time spent in optimization (s)');
disp(t_optimization);
disp('    Number of CG iterations:');
disp(out_nls.iterations)
T = cpdgen([repmat(Best_nls(1),1,d),Best_nls{2}]);
issym = norm(symcheck(reshape(T,I*ones(1,d))));
disp('S3: Symmetry coefficient of solution obtained via LS-CPD');
disp(issym);