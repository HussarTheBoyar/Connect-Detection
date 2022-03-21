function [cvx_status,T1_status] = PowerControl_FullyConnect(Gamma,BETAA,Pu,Phi,N)
%This function uses CVX to solve power control coefficient satisfying
%Uplink rate condition

[M,K] = size(Gamma);
cvx_check = 0;
C = ones(1,K);

while cvx_check == 0


A_max = ones(M,K);
[PC22,deno1,~,rate,~] = PC_comp(Gamma,BETAA,Pu,Phi,N,A_max,C);
tmp_rate = sort(rate);
tmp_K = sum(C);
idx = find(C==1);

deno1 = deno1(idx);
PC22 = PC22(idx,idx);
RReq = log2(3);

cvx_quiet true
cvx_begin %sdp
variables x(1,tmp_K) 
    minimize 0
        subject to
            for k=1:tmp_K
                x(k)*deno1(k)*N*Pu + Pu*[x(1:k-1), x(k+1:tmp_K)]*[PC22(1:k-1,k); PC22(k+1:tmp_K,k)] + N*sum( A_max(:,idx(k)).*Gamma(:,idx(k)) ) <= N^2*Pu*sum( A_max(:,idx(k)).*Gamma(:,idx(k)) )^2*x(k)/(2^RReq-1) ;
            end   
                
            for k=1:tmp_K
                x(k)<=1;
                x(k)>=0;
            end

                            
cvx_end

cvx_check = contains(cvx_status, 'Solve');

if cvx_check == 0

    [value,index] = sort(C.*rate);
    for k = 1:K
        if value(k) == 0
            continue
        else  
            C( index(k) ) = 0;
            break
        end
    end
    [~,~,~,rate,~] = PC_comp_new(Gamma,BETAA,Pu,Phi,N,A_max,C);
end

end

if cvx_check == 1
    cvx_status = tmp_K/K;
end

while sum(rate>=log2(3)) < sum(C)
    [value,index] = sort(C.*rate);
    for k = 1:K
        if value(k) == 0
            continue
        else  
            C( index(k) ) = 0;
            break
        end
    end
    
end

T1_status = sum(C)/K;

end

