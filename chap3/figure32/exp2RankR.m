%%  Experiment 4:
%   In this tutorial we illustrate how the CPD-rank influnce the identification
%   performance for different size of tensor and noise, when use the 
%   KRT-structured LSCPDS [1] routines lscpds_krt_eval, lscpds_krt_algebraic,
%   lscpds_krt_nls, lscpds_krt, and fit the sythenic MIMO Volterra datasets 
%   [2] in this framework.

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

close all; clear variables; clc; rng(654);
%% In order to generate a random KRT-structured LSCPDS problem, we first
% define the size, order, and rank of the symmetric CPD:
p = 2;
M = 3;
I = p*M+1;
l = 1;
d = 3;
R = 5;
N = 50*nchoosek(I+d-1,I-1);
B = rand(I,R);
c = rand(1,R);
V = {B,c};
T = cpdgen([repmat(V(1),1,d),V{2}]);
issym = norm(symcheck(reshape(T,I*ones(1,d))));
strcat("Symmetry coefficient of the ground truth tensor: ",num2str(issym))
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
options.MinRelErr = 10^(-10);
rankest(T,options);
%% The outer loop changes I, and the inner loop changes R
index = 1;
index_w = 1;
Runs=20;
for M = 2:2:8 
    % The outer loop changes I, and update the U and y
    I = p*M+1;
    N = 3000;
    B = rand(I,R);
    c = rand(1,R);
    V = {B,c};
    clear U
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
    
    for R_exp= 1:8
        % The inner loop changes R
        index_w = 1;
        for run=1:Runs
            [Vest1,Vest1_initial,t_initial1,t_optimization1,recording1] = lscpds_krt(U(1:0.7*N,:),d,y(1:0.7*N),R_exp);
            e_y_opt=frob(y(1:0.7*N)-lscpds_krt_eval(Vest1,U(1:0.7*N,:),d))/frob(y(1:0.7*N));
            data.exp1_R.t_total(index_w,index) = t_optimization1+t_initial1;
            data.exp1_R.relerr_optimization(index_w,index) = e_y_opt;
            data.exp1_R.t_ini(index_w,index) = t_initial1;
            data.exp1_R.t_gn(index_w,index) = t_optimization1;
            T = cpdgen([repmat(Vest1(1),1,d),Vest1{2}]);
            data.exp1_R.s(index_w,index) = norm(symcheck(reshape(T,I*ones(1,d))));
            e_y_val=frob(y(0.7*N+1:end)-lscpds_krt_eval(Vest1,U(0.7*N+1:end,:),d))/frob(y(0.7*N+1:end));
            data.exp1_R.relerr_validation(index_w,index) = e_y_val;
            index_w = index_w + 1;
        end
        index = index + 1;
    end
    
end
data.exp1_R.R_exp = 1:8;
data.exp1_R.I_exp = 5:5:20;
