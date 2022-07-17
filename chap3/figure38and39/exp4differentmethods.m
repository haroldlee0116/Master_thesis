clear all, close all
clc
NOISE_ON = 1;
%% System generation
d=3;
p=1;
l=1;
M=5;
n=p*M+1;
Ntrain=1000;
Nval=300;
r=20;

method_index=1;
for SNR=[inf,40,20,-20]
    for index=1:30

B = rand(n,r);
c = rand(1,r);
V = {B,c};

u=randn(Ntrain+Nval,p);
U=makeU(u,M,1);
y = lscpds_krt_eval(V,U,d);
%% Enforcing Method
tic;
% % rank vector construction
if mod(d,2)==1
    % odd degree
    ranks=zeros(1,(d-1)/2);
    for i=1:length(ranks)
        ranks(i)=nchoosek(n+i-1,n-1);
    end
    ranks=[ranks,fliplr(ranks)];
%
%     %even degree
else
    ranks=zeros(1,d/2);
    for i=1:length(ranks)
        ranks(i)=nchoosek(n+i-1,n-1);
    end
    ranks=[ranks,fliplr(ranks(1:d/2-1))];
end
% % noise generation 
if SNR~=inf
    norme=norm(y(1:Ntrain))/10^(SNR/20);
else
    norme=0;
end
e=randn(Ntrain,1);
e=e/norm(e)*norme;

% % identification via pseudo inverse
[TN,tolerance]=rkh2tn(U(1:Ntrain,:)',d); % Algorithm 1 but the sequence is flipped
[Q,S,V]=svd(reshape(TN.core{d},[TN.n(d,1)*n,Ntrain]),'econ');
PE_rank=nchoosek(n+d-1,n-1);
PE_rank=12;
Q=Q(:,1:PE_rank);
S=S(1:PE_rank,1:PE_rank);
V=V(:,1:PE_rank);
for i=1:d-1
   TN2.core{i}= TN.core{i};
   TN2.n(i,:)=TN.n(i,:);
end
TN2.core{d}=reshape((y(1:Ntrain)+e)'*V*inv(S)*Q',[TN.n(d,1),n,1,1]);
% %
TN2.n(d,:)=[TN.n(d,1),n,1,1];
TN3=roundTN(TN2,1e-9);  % uniform rank-3 Volterra kernels
h_hat3 = contract(TN3);
issym = norm(symcheck(reshape(h_hat3,n*ones(1,d))));
v1=toc;
if SNR==inf
    y_noise=y;
else
    y_noise=noisy(y,SNR);
end
% strcat("Symmetry coefficient of solution obtained via pseudoinverse: ",num2str(issym))
yhat_val=sim_volterraTN(u(Ntrain+1:end),transposeTN(TN3));
yhat_train=sim_volterraTN(u(1:Ntrain),transposeTN(TN3));
validation_error=frob(yhat_val(M+1:end)-y_noise(Ntrain+1+M:end))/frob(y_noise(Ntrain+1+M:end));
tranning_error=frob(yhat_train-y(1:Ntrain)+e)/frob(y(1:Ntrain)+e);
storage = Ntrain*l+d*(p*M+1)*(PE_rank^2)*(Ntrain*l+1);
% strcat("Training error of solution obtained via pseudoinverse: ",num2str(validation_error))

recording.method1.train_error(method_index,index)=tranning_error;
recording.method1.val_error(method_index,index)=validation_error;
recording.method1.time(method_index,index)=v1;
recording.method1.symmetric(method_index,index)=issym;
recording.method1.storage(method_index,index)=storage;
%% ALS / MALS method
MAXITR=30;
tic;
ranks=[4,4];
% % rank vector construction
[TN1,etrain1,issym1,sweep,rmax]=mvals(y(1:Ntrain)+e,u(1:Ntrain),M,ranks,1e-20,MAXITR);
v2=toc;
% strcat("Symmetry coefficient of solution obtained via ALS: ",num2str(issym1(1,end)))
% strcat("Trainning error of solution obtained via ALS: ",num2str(etrain1(1,end)))
yhat_val=sim_volterraTN(u(Ntrain+1:end),TN1);
yhat_train=sim_volterraTN(u(1:Ntrain),TN1);
validation_error=frob(yhat_val(M+1:end)-y_noise(Ntrain+1+M:end))/frob(y_noise(Ntrain+1+M:end));
tranning_error=frob(yhat_train-y(1:Ntrain)+e)/frob(y(1:Ntrain)+e);
r=max(ranks);
storage = Ntrain*l+(d-1)*(r^2)*((p*M+1)^2)*(Ntrain*l+1);
recording.method2.train_error(method_index,index)=tranning_error;
recording.method2.val_error(method_index,index)=validation_error;
recording.method2.time(method_index,index)=v2;
recording.method2.symmetric(method_index,index)=issym1(1,end);
recording.method2.storage(method_index,index)=storage;

tic;
[TN,etrain,issym,sweep,rmax]=mvmals2(y(1:Ntrain)+e,u(1:Ntrain),M,d,1e-20,MAXITR);
v3=toc;
% strcat("Symmetry coefficient of solution obtained via MALS: ",num2str(issym(1,end)))
% strcat("Trainning error of solution obtained via MALS: ",num2str(etrain(1,end)))
r=max(max(TN.n));
storage = Ntrain*l+Ntrain*(p*M+1)*r+(p*M+1)*r*l+2*(d-1)*r*r*(p*M+1);
yhat_val=sim_volterraTN(u(Ntrain+1:end),TN);
yhat_train=sim_volterraTN(u(1:Ntrain),TN);
validation_error=frob(yhat_val(M+1:end)-y_noise(Ntrain+1+M:end))/frob(y_noise(Ntrain+1+M:end));
tranning_error=frob(yhat_train-y(1:Ntrain)+e)/frob(y(1:Ntrain)+e);



recording.method3.train_error(method_index,index)=tranning_error;
recording.method3.val_error(method_index,index)=validation_error;
recording.method3.time(method_index,index)=v3;
recording.method3.symmetric(method_index,index)=issym(1,end);
recording.method3.storage(method_index,index)=storage;
%% LS-CPD method
y_train = noisy(y,SNR);
r = 11;
[Vest1,Vest1_initial,t_initial1,t_optimization1,recording1] = lscpds_krt(U(1:Ntrain,:),d,y_train(1:Ntrain),r);
disp(t_optimization1+t_initial1);
e_y_est=frob(y_train(1:Ntrain)-lscpds_krt_eval(Vest1,U(1:Ntrain,:),d))/frob(y_train(1:Ntrain));
e_y_val=frob(y_train(Ntrain+1:end)-lscpds_krt_eval(Vest1,U(Ntrain+1:end,:),d))/frob(y_train(Ntrain+1:end));
% strcat("Trainning error of solution obtained via LS-CPD: ",num2str(e_y_est))
T = cpdgen([repmat(Vest1(1),1,d),Vest1{2}]);
issym = norm(symcheck(reshape(T,n*ones(1,d))));
% strcat("Symmetry coefficient of solution obtained via LS-CPD: ",num2str(issym))
storage = r*((p*M+1)+1)*(3+Ntrain)+Ntrain*((p*M+1)+2)+r*r*((p*M+1)+1)^2;

recording.method5.train_error(method_index,index)=e_y_est;
recording.method5.val_error(method_index,index)=e_y_val;
recording.method5.time(method_index,index)=t_optimization1+t_initial1;
recording.method5.symmetric(method_index,index)=issym;
recording.method5.storage(method_index,index)=storage;
    end
    method_index=method_index+1;
    disp('SNR Update')
end
save recording
