function [cvx_check,cvx_status_num] = PowerControl_TMethod(Gamma,BETAA,Pu,Phi,N,A_max,C)
%This function uses CVX to solve power control coefficient satisfying
%Uplink rate condition

[PC22,deno1,~,rate_PC,~] = PC_comp_new(Gamma,BETAA,Pu,Phi,N,A_max,C);
tmp_K = sum(C);
[~,K] = size(Gamma);
idx = find(C==1);
cvx_status_num = 0;

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
if cvx_check == 1
    cvx_status_num = tmp_K/K;
end

end

