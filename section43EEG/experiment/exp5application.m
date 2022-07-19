%% Experiment 5: Modeling EEG dataset

%   In this experiment, we model the evoked cortical reponses data set [1] using
%   the proposed framework that uses the KRT-structured LSCPDS [2] routines
%   lscpds_krt_eval, lscpds_krt_algebraic, lscpds_krt_nls, lscpds_krt.  

%   line 63-129: record the performance with d=[1,2,3],M=25,R=20
%   line 131-194: record the performance with M=[10,20,30,40], d=2,R=20
%   line 196-261: record the performance with R=[10,15,20,30,40,70], M=30,d=2

%   Authors:    Zhehan Li                  (z.li-50@student.tudelft.nl) 

%   References:
%   [1] Modeling the nonlinear cortical response in EEG evoked by wrist 
%   joint manipulation, by Vlaar et al., IEEE Trans Neural Syst Rehabil Eng 
%   26:205-305, 2018
%   [2] Stijn Hendrikx, Martijn Bousse, Nico Vervliet, Lieven De Lathauwer, 
%   "Algebraic and Optimization Based Algorithms for Multivariate 
%   Regression Using Symmetric Tensor Decomposition", 2019 IEEE 8th 
%   International Workshop on Computational Advances in Multi-Sensor 
%   Adaptive Processing (CAMSAP 2019).

%% Reading EEG dataset
close all
clear 
clc

load('Benchmark_EEG_medium')
% download 'Benchmark_EEG_medium' data [1] from:
% https://drive.google.com/file/d/1XrkV43ZKq-vlcwVz5OUiW2A2N51Bz5Lt/view

S = size(EEGdata,1);          % #participants 
M = size(EEGdata{2}.angle,1); % #realizations

% reorganize data to 'small' version
ds = 8;             % downsample rate (from 2048 to 256 Hz)

% average over periods and downsample, see Section II-B
for Sidx=1:S
    data.u(Sidx,:,:)=squeeze(mean(EEGdata{Sidx}.angle(:,:,1:ds:end),2));
    data.y(Sidx,:,:)=squeeze(mean(EEGdata{Sidx}.comp(:,:,1:ds:end),2));
end
N = size(data.u,3); % #samples

% zero-mean and scale per participant, see Section II-G
data.u=(data.u-repmat(mean(data.u,3),1,1,N))./repmat(mean(std(data.u,[],3),2),1,M,N);
data.y=(data.y-repmat(mean(data.y,3),1,1,N))./repmat(mean(std(data.y,[],3),2),1,M,N);

% response is delayed, so also delay the input, see Section II-G
%  delay with 5 samples, i.e. 19.5 ms
data.u=circshift(data.u,5,3);

fs = 2048/ds;       % sample frequency
T = N/fs;           % segment length [s]
t=(0:N-1)'/fs;      % time vector

u=data.u;           % input, handle angle (normalized)
y=data.y;           % output, ICA component with highest SNR (normalized)
%% Individual Trainning: Different d
i=1;   % Participant 1
val=1; % First realization as validation data
r=20;
M=25;
options.Initialization=@cpd_rnd;

% validation input and output
K=[];
for k=1:7
    if k~=val
        K=[K,k];
    end
end
y_val = reshape(y(i,val,M+1:end),[N-M,1]);
u_val = reshape(u(i,val,:),[N,1]);
for t=1:N
    for k=1:M
        k=k-1;
        if t-k<1
            u_val(t,k+1)=1;
        else
            u_val(t,k+1)=u_val(t-k);
        end
    end
end
U_val=u_val(M+1:end,:);

%tranning input and output
for real=K
    u_est = reshape(u(i,real,:),[N,1]);
    y_est = reshape(y(i,real,M+1:end),[N-M,1]);
    for t=1:N
        for k=1:M
            k=k-1;
            if t-k<1
                u_est(t,k+1)=1;
            else
                u_est(t,k+1)=u_est(t-k);
            end
        end
    end
    % the input data is stored seperately for each realization
    
    if real==K(1)
        U_est=u_est(M+1:end,:);
        Y_est=y_est;
        
    else
        U_est=[U_est;u_est(M+1:end,:)];
        Y_est=[Y_est;y_est];
    end
end

% LS-CPD optimization:
for d=1:3
    [Vest1,Vest1_initial,t_initial1,t_optimization1,recording1] = lscpds_krt(U_est,d,Y_est,r,options);
    y_hat_est(:,d) = lscpds_krt_eval(Vest1,U_est,d)-ones(6*N-6*M,1)*mean(lscpds_krt_eval(Vest1,U_est,d));
    y_hat_val(:,d) = lscpds_krt_eval(Vest1,U_val,d)-ones(N-M,1)*mean(lscpds_krt_eval(Vest1,U_val,d));
    vaf_est(d) = 1-var(Y_est-y_hat_est(:,d))/var(Y_est);
    vaf_val(d) = 1-var(y_val-y_hat_val(:,d))/var(y_val);
    
    recording.y_est(1:6*(N-M),1)=Y_est;
    recording.y_val(1:N-M,1)=y_val;
    recording.y_hat_est(:,d)=y_hat_est(:,d);
    recording.y_hat_val(:,d)=y_hat_val(:,d);
end
%% Individual Tranning: Different M
i=1;
val=1; % First realization as validation data
r=80;
d=2;
options.Initialization=@cpd_rnd;
y_hat_est=zeros(6*(N-10),4);
y_hat_val=zeros(N-10,4);
index=1;
for M=10:10:40
    % validation input and output
    K=[];
    for k=1:7
        if k~=val
            K=[K,k];
        end
    end
    y_val = reshape(y(i,val,M+1:end),[N-M,1]);
    u_val = reshape(u(i,val,:),[N,1]);
    for t=1:N
        for k=1:M
            k=k-1;
            if t-k<1
                u_val(t,k+1)=1;
            else
                u_val(t,k+1)=u_val(t-k);
            end
        end
    end
    U_val=u_val(M+1:end,:);
    
    %tranning input and output
    for real=K
        u_est = reshape(u(i,real,:),[N,1]);
        y_est = reshape(y(i,real,M+1:end),[N-M,1]);
        for t=1:N
            for k=1:M
                k=k-1;
                if t-k<1
                    u_est(t,k+1)=1;
                else
                    u_est(t,k+1)=u_est(t-k);
                end
            end
        end
        % the input data is stored seperately for each realization
        
        if real==K(1)
            U_est=u_est(M+1:end,:);
            Y_est=y_est;
            
        else
            U_est=[U_est;u_est(M+1:end,:)];
            Y_est=[Y_est;y_est];
        end
    end
    % LS-CPD
    [Vest1,Vest1_initial,t_initial1,t_optimization1,recording1] = lscpds_krt(U_est,d,Y_est,r,options);
    y_hat_est(1:6*(N-M),index) = lscpds_krt_eval(Vest1,U_est,d)-ones(6*N-6*M,1)*mean(lscpds_krt_eval(Vest1,U_est,d));
    y_hat_val(1:N-M,index) = lscpds_krt_eval(Vest1,U_val,d)-ones(N-M,1)*mean(lscpds_krt_eval(Vest1,U_val,d));
    index=index+1;
    recoring2.y_hat_val(1:N-50,index) = y_hat_val(51-M:N-M,index);
end
recording2.y_val=y_val;
%% Individual Tranning: Different R
i=1;   % Participant 1
val=1; % First realization as validation data
M=30;
options.Initialization=@cpd_rnd;
d=2;

% validation input and output
K=[];
for k=1:7
    if k~=val
        K=[K,k];
    end
end
y_val = reshape(y(i,val,M+1:end),[N-M,1]);
u_val = reshape(u(i,val,:),[N,1]);
for t=1:N
    for k=1:M
        k=k-1;
        if t-k<1
            u_val(t,k+1)=1;
        else
            u_val(t,k+1)=u_val(t-k);
        end
    end
end
U_val=u_val(M+1:end,:);

%tranning input and output
for real=K
    u_est = reshape(u(i,real,:),[N,1]);
    y_est = reshape(y(i,real,M+1:end),[N-M,1]);
    for t=1:N
        for k=1:M
            k=k-1;
            if t-k<1
                u_est(t,k+1)=1;
            else
                u_est(t,k+1)=u_est(t-k);
            end
        end
    end
    % the input data is stored seperately for each realization
    
    if real==K(1)
        U_est=u_est(M+1:end,:);
        Y_est=y_est;
        
    else
        U_est=[U_est;u_est(M+1:end,:)];
        Y_est=[Y_est;y_est];
    end
end

% LS-CPD optimization:
index=1;
for r=[10,15,20,30,40,70]
    [Vest1,Vest1_initial,t_initial1,t_optimization1,recording1] = lscpds_krt(U_est,d,Y_est,r,options);
    y_hat_est(:,index) = lscpds_krt_eval(Vest1,U_est,d)-ones(6*N-6*M,1)*mean(lscpds_krt_eval(Vest1,U_est,d));
    y_hat_val(:,index) = lscpds_krt_eval(Vest1,U_val,d)-ones(N-M,1)*mean(lscpds_krt_eval(Vest1,U_val,d));
    recording4.y_val(1:N-M,1)=y_val;
    recording4.y_hat_val(:,d)=y_hat_val(:,d);
    
    index=index+1;
end


