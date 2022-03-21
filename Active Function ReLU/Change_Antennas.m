clc
clear 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Uplink
%Consider a square are of DxD m^2
%M distributed APs serves K terminals, they all randomly located in the area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Inital parameters
M=200; %number of access points
K=60; %number of terminals
N_Antenna_vec = 1:6; %number of AP equiped antenna

D=1; %in kilometer
tau=10;%training length
[U,S,V]=svd(randn(tau,tau));%U includes tau orthogonal sequences 

T=200;%coherence interval
B=20; %Mhz
Hb = 15; % Base station height in m
Hm = 1.65; % Mobile height in m
f = 1900; % Frequency in MHz
aL = (1.1*log10(f)-0.7)*Hm-(1.56*log10(f)-0.8);
L = 46.3+33.9*log10(f)-13.82*log10(Hb)-aL;

power_f=0.1; %uplink power: 100 mW
noise_p = 10^((-203.975+10*log10(20*10^6)+9)/10); %noise power
Pu = power_f/noise_p;%nomalized receive SNR
Pp=Pu;%pilot power

d0=0.01;%km
d1=0.05;%km

N_Loop=100 ;
AP_Threshold = 30;
R_cf_1_min=zeros(1,N_Loop);%min rate, cell-free, without power allocation

Req = 2.0;

Req_step = 1;


range = 5;

rp_T1_satis = zeros(1,range);
rp_T2_satis = zeros(1,range);
rp_TWC_satis = zeros(1,range);
rp_H_satis = zeros(1,range);
rp_T2_opt_satis = zeros(1,range);
rp_H_opt_satis = zeros(1,range);
rp_idx = 1;

parfor rp_idx=1:6
    N_Antenna = N_Antenna_vec(rp_idx);
    
num_T1 = zeros(N_Loop,K,Req_step);
rate_user_T1 = zeros(N_Loop,K,Req_step);
status_T1 = zeros(N_Loop,K,Req_step);
RReq = zeros(N_Loop,K);
A_T1   = zeros(M,K,N_Loop);
 
num_T2 = zeros(N_Loop,K,Req_step);
rate_user_T2 = zeros(N_Loop,K,Req_step);
status_T2 = zeros(N_Loop,K,Req_step);
A_T2   = zeros(M,K,N_Loop);
AP_used = zeros(N_Loop,1);

num_H = zeros(N_Loop,K,Req_step);
rate_user_H = zeros(N_Loop,K,Req_step);
status_H = zeros(N_Loop,K,Req_step);

num_TWC = zeros(N_Loop,K,Req_step);
rate_user_TWC = zeros(N_Loop,K,Req_step);
status_TWC = zeros(N_Loop,K,Req_step);
A_TWC   = zeros(M,K,N_Loop);

num = zeros(N_Loop,K,Req_step);
rate_user     = zeros(N_Loop,K);
rate_user_opt = zeros(N_Loop,K);
A_opt = zeros(M,K,N_Loop);

status = zeros(N_Loop,Req_step);
status_SVM = zeros(N_Loop,Req_step);
cvx_status = zeros(N_Loop,Req_step);
cvx_status_SVM = zeros(N_Loop,Req_step);
cvx_fully_status = zeros(N_Loop,Req_step);
RReq = zeros(N_Loop,K);
A   = zeros(M,K,N_Loop);
unSatisfy = zeros(Req_step,N_Loop);
num_Satisfy = zeros(Req_step,K,N_Loop);
n = 1;
check =1;

notSatis = zeros(1,N_Loop);
Satis = zeros(1,N_Loop);
mean_H_rate = zeros(1,N_Loop);
min_H_rate = zeros(1,N_Loop);
    
while n<=N_Loop 
% for (n = 1:N)
    clc
    rp_idx
    n
%%%%%Randomly locations of M APs%%%%
AP=unifrnd(-D/2,D/2,M,2);

%Wrapped around (8 neighbor cells)
D1=zeros(M,2);
D1(:,1)=D1(:,1)+ D*ones(M,1);
AP(:,:,2)=AP(:,:,1)+D1;

D2=zeros(M,2);
D2(:,2)=D2(:,2)+ D*ones(M,1);
AP(:,:,3)=AP(:,:,1)+D2;

D3=zeros(M,2);
D3(:,1)=D3(:,1)- D*ones(M,1);
AP(:,:,4)=AP(:,:,1)+D3;

D4=zeros(M,2);
D4(:,2)=D4(:,2)- D*ones(M,1);
AP(:,:,5)=AP(:,:,1)+D4;

D5=zeros(M,2);
D5(:,1)=D5(:,1)+ D*ones(M,1);
D5(:,2)=D5(:,2)- D*ones(M,1);
AP(:,:,6)=AP(:,:,1)+D5;

D6=zeros(M,2);
D6(:,1)=D6(:,1)- D*ones(M,1);
D6(:,2)=D6(:,2)+ D*ones(M,1);
AP(:,:,7)=AP(:,:,1)+D6;

D7=zeros(M,2);
D7=D7+ D*ones(M,2);
AP(:,:,8)=AP(:,:,1)+D7;

D8=zeros(M,2);
D8=D8- D*ones(M,2);
AP(:,:,9)=AP(:,:,1)+D8;

%Randomly locations of K terminals:
Ter=zeros(K,2,9);

Ter(:,:,1) = unifrnd(-D/2, D/2, K, 2);

%Wrapped around (8 neighbor cells)
D1=zeros(K,2);
D1(:,1)=D1(:,1)+ D*ones(K,1);
Ter(:,:,2)=Ter(:,:,1)+D1;

D2=zeros(K,2);
D2(:,2)=D2(:,2)+ D*ones(K,1);
Ter(:,:,3)=Ter(:,:,1)+D2;

D3=zeros(K,2);
D3(:,1)=D3(:,1)- D*ones(K,1);
Ter(:,:,4)=Ter(:,:,1)+D3;

D4=zeros(K,2);
D4(:,2)=D4(:,2)- D*ones(K,1);
Ter(:,:,5)=Ter(:,:,1)+D4;

D5=zeros(K,2);
D5(:,1)=D5(:,1)+ D*ones(K,1);
D5(:,2)=D5(:,2)- D*ones(K,1);
Ter(:,:,6)=Ter(:,:,1)+D5;

D6=zeros(K,2);
D6(:,1)=D6(:,1)- D*ones(K,1);
D6(:,2)=D6(:,2)+ D*ones(K,1);
Ter(:,:,7)=Ter(:,:,1)+D6;

D7=zeros(K,2);
D7=D7+ D*ones(K,2);
Ter(:,:,8)=Ter(:,:,1)+D7;

D8=zeros(K,2);
D8=D8- D*ones(K,2);
Ter(:,:,9)=Ter(:,:,1)+D8;

sigma_shd=8; %in dB
D_cor=0.1;

%%%%%%Create the MxK correlated shadowing matrix %%%%%%%

    %%%%M correlated shadowing cofficients of M APs:
    Dist=zeros(K,K);%distance matrix
    Cor=zeros(K,K);%correlation matrix

    for m1=1:M
        for m2=1:M
            Dist(m1,m2) = min([norm(AP(m1,:,1)-AP(m2,:,1)), norm(AP(m1,:,1)-AP(m2,:,2)),norm(AP(m1,:,1)-AP(m2,:,3)),norm(AP(m1,:,1)-AP(m2,:,4)),norm(AP(m1,:,1)-AP(m2,:,5)),norm(AP(m1,:,1)-AP(m2,:,6)),norm(AP(m1,:,1)-AP(m2,:,7)),norm(AP(m1,:,1)-AP(m2,:,8)),norm(AP(m1,:,1)-AP(m2,:,9)) ]); %distance between AP m1 and AP m2
            %Dist(m1,m2)=norm(AP(m1,:,1)-AP(m2,:,1));
            Cor(m1,m2)=exp(-log(2)*Dist(m1,m2)/D_cor);
        end
    end
    A1 = chol(Cor,'lower');
    x1 = randn(M,1);
    sh_AP = A1*x1;
    for m=1:M
        sh_AP(m)=(1/sqrt(2))*sigma_shd*sh_AP(m)/norm(A1(m,:));
    end

    %%%%K correlated shadowing matrix of K terminal:
    Dist=zeros(K,K);%distance matrix
    Cor=zeros(K,K);%correlation matrix

    for k1=1:K
        for k2=1:K
            Dist(k1,k2)=min([norm(Ter(k1,:,1)-Ter(k2,:,1)), norm(Ter(k1,:,1)-Ter(k2,:,2)),norm(Ter(k1,:,1)-Ter(k2,:,3)),norm(Ter(k1,:,1)-Ter(k2,:,4)),norm(Ter(k1,:,1)-Ter(k2,:,5)),norm(Ter(k1,:,1)-Ter(k2,:,6)),norm(Ter(k1,:,1)-Ter(k2,:,7)),norm(Ter(k1,:,1)-Ter(k2,:,8)),norm(Ter(k1,:,1)-Ter(k2,:,9)) ]); %distance between Terminal k1 and Terminal k2
            Cor(k1,k2)=exp(-log(2)*Dist(k1,k2)/D_cor);
        end
    end
    A2 = chol(Cor,'lower');
    x2 = randn(K,1);
    sh_Ter = A2*x2;
    for k=1:K
        sh_Ter(k)=(1/sqrt(2))*sigma_shd*sh_Ter(k)/norm(A2(k,:));
    end

%%% The shadowing matrix
Z_shd=zeros(M,K);
for m=1:M
    for k=1:K
        Z_shd(m,k)= sh_AP(m)+ sh_Ter(k);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create an MxK large-scale coefficients beta_mk
BETAA = zeros(M,K);
dist=zeros(M,K);
for m=1:M  
    for k=1:K
    [dist(m,k),index] = min([norm(AP(m,:,1)-Ter(k,:,1)), norm(AP(m,:,2)-Ter(k,:,1)),norm(AP(m,:,3)-Ter(k,:,1)),norm(AP(m,:,4)-Ter(k,:,1)),norm(AP(m,:,5)-Ter(k,:,1)),norm(AP(m,:,6)-Ter(k,:,1)),norm(AP(m,:,7)-Ter(k,:,1)),norm(AP(m,:,8)-Ter(k,:,1)),norm(AP(m,:,9)-Ter(k,:,1)) ]); %distance between Terminal k and AP m
    if dist(m,k)<d0
         betadB=-L - 35*log10(d1) + 20*log10(d1) - 20*log10(d0);
    elseif ((dist(m,k)>=d0) && (dist(m,k)<=d1))
         betadB= -L - 35*log10(d1) + 20*log10(d1) - 20*log10(dist(m,k));
    else
    betadB = -L - 35*log10(dist(m,k)) + Z_shd(m,k); %large-scale in dB
    end
    
    BETAA(m,k)=10^(betadB/10); 
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                       K-MEANS                          %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pilot Asignment: (random choice)
Ter_k = Ter(:,:,1);
if K >1 
    [Ter(:,:,1), Phii_3, ~, ~] = Kmeans_Bjonson_new(U, Ter_k, AP(:,:,1));
else
    Phii_3 = U;
end    

%% Create Gamma matrix
Gammaa_2 = zeros(M,K);
mau_2=zeros(M,K);
Phii_cf_2 = Phii_3;


for m=1:M
    for k=1:K
        mau_2(m,k)=norm( (BETAA(m,:).^(1/2)).*(Phii_cf_2(:,k)'*Phii_cf_2))^2;
    end
end

for m=1:M
    for k=1:K
        Gammaa_2(m,k)=tau*Pp*BETAA(m,k)^2/(tau*Pp*mau_2(m,k) + 1);
    end
end

Rate_index = zeros(M,K);

for step =1:Req_step
    RReq_co = 1*step;
    
    %Te_test(Gammaa_2,BETAA,Phii_3,Pu,N_Antenna);
    %[rate_user_H(n,:,step),status_H(n,:,step),num_H(n,:,step),~]= test_connect(Gammaa_2,BETAA,Phii_3,Pu,RReq_co,Req,N_Antenna);
    %[Rate_index] = Rate_compueting(BETAA,Gammaa_2,Phii_3,Pu,N_Antenna);
    [Check,status_T2(n,:,step),status(n,step),cvx_status(n,step),cvx_fully_status(n,step)] = Rate_single_computing(BETAA,Gammaa_2,Phii_3,Pu,N_Antenna,Req,RReq_co);
    [status_SVM(n,step),cvx_status_SVM(n,step)] = SVM_computing(BETAA,Gammaa_2,Phii_3,Pu,N_Antenna,Req,RReq_co);
    %[rate_user_TWC(n,:,step),A_TWC(:,:,n),status_TWC(n,:,step),num_TWC(n,:,step),~] = Rate_computing_TWC(Gammaa_2,BETAA,Phii_3,Pu,RReq_co,Req,N_Antenna);    
    [rate_user_H(n,:,step),status_H(n,:,step),num_H(n,:,step),~]                = Rate_computing_H(Gammaa_2,BETAA,Phii_3,Pu,RReq_co,Req,N_Antenna);

    
end
n = n+1;
    
end

rp_T1_satis(rp_idx) = (mean( sum( status_T2(:,:,Req_step),2 ) ))*100/K;
rp_T2_satis(rp_idx) = mean(status)*100;
rp_SVM_satis(rp_idx) = mean(status_SVM)*100;
rp_TWC_satis(rp_idx) = (K-mean( sum( status_TWC(:,:,Req_step),2 ) ))*100/K;
rp_H_satis(rp_idx) = (K-mean( sum( status_H(:,:,Req_step),2 ) ))*100/K;
rp_T2_opt_satis(rp_idx) = mean(cvx_status)*100;
rp_NN_opt_satis(rp_idx) = mean(cvx_status_SVM)*100;
rp_H_opt_satis(rp_idx) = mean(cvx_fully_status)*100;


end
   
figure(5)
grid on;
hold on;

plot(N_Antenna_vec,rp_T2_satis,'r');
plot(N_Antenna_vec,rp_T1_satis,'k');
plot(N_Antenna_vec,rp_H_satis,'b');
plot(N_Antenna_vec,rp_SVM_satis,'m');
% plot(number_user,rp_T2_opt_satis,'k');
% plot(number_user,rp_TWC_satis,'k');
% plot(number_user,rp_TWC_satis,'g');

xlabel('Number of Antennas');
ylabel('Number of user satisfy QoS(%)');
legend('Lemma Method', 'T Method', 'Fully Connected Method','Neural Network method');

% legend('T Method', 'H Method', 'TWC Method');

figure(6)

grid on;
hold on;

plot(N_Antenna_vec,rp_T2_opt_satis,'r');
plot(N_Antenna_vec,rp_H_opt_satis,'k');
plot(N_Antenna_vec,rp_NN_opt_satis,'m');

xlabel('Number of Antennas');
ylabel('Number of user satisfy QoS(%)');

legend('Lemma Method', 'H Method','Neural Network Method');



save change_Antennas

