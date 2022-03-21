function [status,cvx_status] = SVM_computing(BETAA,Gamma,Phi,Pu,N,Req,RReq_co)
%This function is used to compute AP single connected uplink rate and
%maximum uplink rate based AP connected.
%
%This function returns status (satis RReq or not) with input net connecting variable C
%   Detailed explanation goes here
[M,K] = size(Gamma);
C = ones(1,K);
Check = 0;         %If there exists a rate value that does not satisfy RReq, check=0
reject_index = [];
PhiPhi = zeros(K,K);
cvx_check = 0;
cvx_status = 0;


for ii=1:K
     for k=1:K
         PhiPhi(ii,k) = round( Phi(:,k)'*Phi(:,ii) );
     end
end
    
count = 1;

while Check == 0
    Rate_single_user = zeros(M,K);
    Rate_user = zeros(1,K);
    PC_2 = zeros(K,K);
    IN   = zeros(1,K);


%Computing AP single connected uplink rate and sorting them 

for m_iter = 1:M
    A = zeros(M,K);
    A(m_iter,:) = ones(1,K);
    Gamma_new = Gamma;
    for ii=1:K
         for k=1:K
             PC_2(ii,k) = sum( (C(ii)*A(m_iter,k)*Gamma_new(m_iter,k)/BETAA(m_iter,k))*BETAA(m_iter,ii)  )*PhiPhi(ii,k);
         end
    end
    PC22=N^2*(abs(PC_2)).^2;


    for k=1:K
        deno1=0;
        for m=1:M
            deno1=deno1 + C(k)*A(m,k)*Gamma_new(m,k)*sum(BETAA(m,:));
        end
        IN(k) = (N*sum(Gamma_new(:,k).*A(:,k)) + Pu*deno1*N + Pu*sum(PC22(:,k)) - Pu*PC22(k,k));
        

        if C(k)==1
            Rate_single_user(m_iter,k) = log2( 1+Pu*N^2*( sum(C(k)*A(:,k).*Gamma_new(:,k)) )^2/IN(k) );
        end
    end 

end

[~,Rate_index] = sort(-1*Rate_single_user);
RReq = log2(1+Req);

%Check max rate of user 

% [value,index,Rate_connect] = test_connect(Gamma,BETAA,Phi,Pu,RReq_co,Req,N,C,Rate_index);
% [value,index,~,A] = test_connect_2(Gamma,BETAA,Phi,Pu,RReq_co,Req,N,C,Rate_index);
if count == 1
    [A_max,value_1] = test_connect_new(Gamma,BETAA,Phi,Pu,RReq_co,Req,N,C,Rate_index,Rate_single_user);
    [~,index] = sort(value_1);
    
    reject_threshold = round(0.05*K);
    for rej_count = 1:reject_threshold
        C(index(rej_count)) = 0;
    end
    
%     A_NN = NN_model_new(W1,W2,B1,B2,M,K,Gamma,BETAA,Rate_single_user,C);
    count = count+1;
    continue
end

   if count == 2
       [W1,W2,B1,B2] = Neural_network(Gamma,BETAA,Rate_single_user,A_max,Pu,N,C);
   end

% Use SVM model

[value,A_NN] = NN_rate(W1,W2,B1,B2,Rate_single_user,Gamma,BETAA,Phi,Pu,N,C);
% [sort_value,idx_value] = sort(value);
% %CVX solving

if cvx_check == 0
    [cvx_check,cvx_status] = PowerControl_TMethod(Gamma,BETAA,Pu,Phi,N,A_NN,C);
end

% 
% %Compare max rate with QoS
% 
if sum(value>=RReq) == sum(C)
    Check = 1;
else
    [sort_value,idx_value] = sort(value);
    for k=1:K
        if sort_value(k) == 0
            continue
        else
            C( idx_value(k) ) = 0;
            break
        end
    end
end





%Compute connected matrix A if Check == 1



end


status = sum(C)/K;



    
end

