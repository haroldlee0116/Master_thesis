% Please acknowledge the use of this data in any publications by refering
% to the article: "Modeling the nonlinear cortical response in EEG evoked
% by wrist joint manipulation" by 
%  Vlaar et al., IEEE Trans Neural Syst Rehabil Eng 26:205-305, 2018
%
% For the 'small' dataset, the data has been averaged, downsampled (Section
% II-B), scaled and time-delayed (Section II-G). 
% For this 'medium' dataset, the data is not averaged, not downsampled, not
% scaled and not time-delayed
%
% The data is organized in a struct (EEGdata), every cell in the struct
% contains the data of a participant. The input (EEGdata{S}.angle) and
% output (EEGdata{S}.comp) are matrices with dimensions:
%  [#realization[M], #periods[P], #samples[N]
% where S indicates the participant number (1-10).
%
% The data may be used, copied, or redistributed as long as it is not sold
% and this copyright notice is reproduced on each copy made. The data is
% provided as is, without any express or implied warranties whatsoever. 
%
% On behalf of all authors, 
% Alfred C. Schouten
% Delft Laboratory for NeuroMechanics and Motor Control (NMC lab)
% Department of Biomechanical Engineering
% Delft University of Technology
%
% Mekelweg 2
% 2628 CD Delft
% The Netherlands
% e-mail: A.C.Schouten@tudelft.nl
% url:    www.3me.tudelft.nl/nmc
%
% February 20, 2019

close all
clear 
clc
% download 'Benchmark_EEG_medium' data from:
% https://drive.google.com/file/d/1XrkV43ZKq-vlcwVz5OUiW2A2N51Bz5Lt/view
load('Benchmark_EEG_medium')

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
r=80;
M=20;
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

% input signal pe condition
rank(U_est)
nchoosek(1+M-1,M-1)

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
y_hay_val=zeros(N-10,4);
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
%     recording4.y_est(1:6*(N-M),1)=Y_est;
      recording4.y_val(1:N-M,1)=y_val;
%     recording4.y_hat_est(:,d)=y_hat_est(:,d);
      recording4.y_hat_val(:,d)=y_hat_val(:,d);
    
    index=index+1;
end


