%%  Experiment 3: PCR under different tensor size

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
%% Performance compensation rate experiment
p = 2;
a = 1;
R = 7;
for d = [3,4,5,6]
    for M= [4,5,6]
        I = p*M+1;
        N = nchoosek(I+d-1,I-1);
        clear B;
        for i=1:R
            B(:,i)=abs(1e0*randn(1))*exp(-randi(10,1).*[1:I]);
        end
        c = rand(1,R);
        V = {B,c};
        u = rand(N,p);
        clear U;
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
        options.Initialization = @lscpds_krt_algebraic;
        [Vest1,Vest1_initial,t_initial1,t_optimization1,recording1] = lscpds_krt(U,d,y,R,'Display',true,options);
        for trial=1:30
            options.Initialization = @cpd_rnd;
            [Vest2,Vest2_initial,t_initial2,t_optimization2,recording2] = lscpds_krt(U,d,y,R,'Display',true,options);
            r=log10(recording2.Algorithm.fval(1)/recording1.Algorithm.fval(1));
            t=log10(t_initial2/t_initial1);
            pcr(a)=100*(r+t)/r;
            a=a+1;
        end
    end
end
pcr=pcr';