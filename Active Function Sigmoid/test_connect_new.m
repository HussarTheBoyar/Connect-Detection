function [A_max,value] = test_connect_new(Gamma,BETAA,Phi,Pu,RReq_co,Req,N,C,Rate_index,Rate_value)
%M: number of AP
%K: Number of user
[M,K] = size(Gamma);
PhiPhi = zeros(K,K);
for ii=1:K
    for k=1:K
        PhiPhi(ii,k) = round( Phi(:,k)'*Phi(:,ii) );
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%     THEORY METHOD     %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rate_user_new = zeros(M,K);
A    = zeros(M,K);
Check = 0;
DS = zeros(1,K);
IN_1 = zeros(1,K);
rate_single = Rate_value;
Rate_combine=zeros(1,K);
A_max = zeros(M,K);
RReq = RReq_co*log2(1+Req);
theory_status = ones(1,K);
tmp_PC = zeros(K,K);
PC2 = zeros(K,K);
itial_check = zeros(1,K);
IN_intial = zeros(1,K);


%Set up IN_1itial parameter 
for k=1:K
    if C(k) == 1
        A_max(Rate_index(1,k),k) = 1;
        DS(k) = Gamma(Rate_index(1,k),k);
        IN_1(k) = ( Pu*N^2*DS(k)^2 )/( 2^(rate_single( Rate_index(1,k),k )) -1 );
        Rate_combine(k) = rate_single(Rate_index(1,k),k);
    end
end


for ii=1:K
    for k=1:K
        PC2(ii,k) = sum( (C(ii)*A_max( Rate_index(1,k),k )*Gamma( Rate_index(1,k),k )/BETAA( Rate_index(1,k),k )).*BETAA( Rate_index(1,k),ii )  )*PhiPhi(ii,k);
    end
end



for m=2:M
    tmp_IN = zeros(1,K);
    tmp_A  = zeros(M,K);
    
%Set up single connection to check    
    for k = 1:K
        if C(k) == 1
            tmp_A( Rate_index(m,k),k ) = 1;
        end
    end
    
  
%Computing single connection IN_1     
    for k =1:K
        if C(k) == 1
            tmp_IN(k) = Pu*N^2*Gamma( Rate_index(m,k),k )^2 /(2^rate_single( Rate_index(m,k),k )-1);
        end
    end

%Computing single PC2    

    for ii=1:K
        for k=1:K
            tmp_PC(ii,k) = sum( (C(ii)*tmp_A( Rate_index(m,k),k )*Gamma( Rate_index(m,k),k )/BETAA( Rate_index(m,k),k )).*BETAA( Rate_index(m,k),ii )  )*PhiPhi(ii,k);
        end
    end
    
    PC_check = tmp_PC.*PC2;
    
    for k=1:K
        PC_check(k,k) = 0;
    end
    
%Connecting condition

    connection_check = zeros(1,K);
    tmp_Gamma = zeros(1,K);
    tmp_rate  = zeros(1,K);
    
    for k=1:K
        if C(k) == 1
            tmp_Gamma(k) = Gamma( Rate_index(m,k),k );
            tmp_rate(k) = rate_single( Rate_index(m,k), k );
            connection_check(k) = ( tmp_Gamma(k)*DS(k) - 1/2*tmp_Gamma(k)^2*( Rate_combine(k)/tmp_rate(k) -1 ) -sum(PC_check(:,k)) > 0 );
        end
    end

%Update parameter: IN_1, PC2, Rate, connection matrix    
    for k = 1:K
        if connection_check(k) == 1
            A_max( Rate_index(m,k),k ) = 1;     
            tmp_A(k) = 0;
        end
    end
    
    
    PC2 = PC2+tmp_PC.*connection_check; 
    DS = DS + connection_check.*tmp_Gamma;
    PC_check = connection_check.*tmp_PC.*PC2;
    
    
    for k=1:K
        PC_check(k,k) = 0;
    end
    
    IN_1 = IN_1+tmp_IN.*connection_check+2*N^2*Pu*sum(PC_check);
    for k = 1:K
        Rate_combine(k) = log2(1+N^2*Pu*DS(k)^2/IN_1(k));
    end
        
    
    for k = 1:K
        if  C(k) == 0
            Rate_combine(k) = 0;
        end
    end
    
end
    
value = Rate_combine;



if sum(value>=RReq) == sum(C)
    Check = 1;
end

if Check == 1
    for k = 1:K
        if C(k) == 0
            continue
        end
     
      
        for m=1:M
            if Rate_combine(m,k)<RReq
                A( Rate_index(m,k),k ) = A_max(Rate_index(m,k),k);
            elseif Rate_combine(m,k)>=RReq
                A( Rate_index(m,k),k ) = A_max(Rate_index(m,k),k);
                break
            end
        end
    end
end

end

