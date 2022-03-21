function [rate_user,A,status,num,RReq,test_reject,stop_threshold,AP_used_per_user] = Rate_computing_T2(Gamma,BETAA,Phi,Pu,RReq_co,Req,A_T1_check,N)
%M: number of AP
%K: Number of user
[M,K] = size(Gamma);
rate_user = zeros(1,K);
PC_2 = zeros(K,K);
IN   = zeros(1,K);
A    = zeros(M,K);
connect = ones(1,K);
status = zeros(1,K);
check  = 0;
Gamma_frac = zeros(1,K);
count = 1;
tau   = 10;


while check == 0
    Gamma_new = Gamma.*connect;
    BETAA     = BETAA.*connect;
    check = 1;
    for ii=1:K
         for k=1:K
             PC_2(ii,k) = sum( (Gamma_new(:,k)./BETAA(:,k)).*BETAA(:,ii)  )*Phi(:,k)'*Phi(:,ii);
         end
    end
    PC22=N^2*(abs(PC_2)).^2;


    for k=1:K
        deno1=0;
        for m=1:M
            deno1=deno1 + Gamma_new(m,k)*sum(BETAA(m,:));
        end
        IN(k) = N*(sum(Gamma_new(:,k)) + Pu*deno1*N + Pu*sum(PC22(:,k)) - Pu*PC22(k,k));
    end 
    
    %computing Gamma require
    RReq = RReq_co*log2(1+Req);
    Gamma_req = zeros(1,K);
    for k = 1:K
        if connect(k) == 1
            SINR_req = 2^RReq-1;
            Gamma_req(k) = sqrt(SINR_req*IN(k)/(Pu*N^2));
        end
    end    

    for k=1:K
        Gamma_check = 0;
        [tmp_Gamma, idx_Gamma] = sort(-1*Gamma_new(:,k));
        tmp_Gamma = -1*tmp_Gamma;
        
        for m=1:M
            Gamma_check = Gamma_check+tmp_Gamma(m);
            A( idx_Gamma(m) ,k) = 1;
            if Gamma_check>=Gamma_req(k)
                break
            end
        end
        
%         if Gamma_check < Gamma_req
%             connect(k) = 0;
%         end
        
        Gamma_frac(k) = Gamma_check/Gamma_req(k);
        
    end
    
    [connect_check, connect_idx] = sort(Gamma_frac);

    for k=1:K
        if connect_check(k)<1 && connect_check(k)>0
            connect( connect_idx(k) ) = 0;
            test_reject(count) = connect_idx(k);
            check = 0;
            A(:,k) = zeros(M,1); 
            break
        end
    end
   
    count = count+1;

        
    
end

%Number of APs connect to user
num = sum(A);
status = -( connect -  ones(1,K));

Gamma_new = Gamma_new.*A;

for k = 1:K
    if connect(k) == 1
        rate_user(k) = log2( 1 + Pu*N^2*sum(Gamma_new(:,k))^2/IN(k) );
    end
end

idx = find(connect == 1);

AP_used = 0;

for k = 1:sum(connect)
    AP_used = AP_used + sum(A(:,k));
end

AP_used_per_user = AP_used/sum(connect);

%Test CVX

 C   = ones(1,K);
 count = 1;        %count CVX step
 count_1 = 1;
 check = 0;
 stop = sum(status);
 disp(stop);
 rate = zeros(1,K);
 tau = 10;
  
 stop_threshold = length(test_reject);
 

 check_reject = find(connect == 0);
end

