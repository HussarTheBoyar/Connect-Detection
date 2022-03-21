function [unSatisfy, num_Satisfy,check,rate,A,power_co] = solvePowerControl(Gammaa,BETAA,Phii,Pu,RReq_co,Req,A,sort_order,num,reject,status,test_reject,stop_threshold,N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   C: Considered vector, =1 if we considering it in cvx prob, =0 if not
%   A: Choosing matrix

[M,K]  = size(Gammaa);
 C   = ones(1,K);
 count = 1;        %count CVX step
 count_1 = 1;
 check = 0;
 stop = sum(status);
%  stop_threshold = stop;
 disp(stop_threshold);
 rate = zeros(1,K);
 tau = 10;
 
 test_reject = [test_reject zeros(1, K-length(test_reject) )];

% test_reject = reject;

Gamma_H = mean(Gammaa);

for k = 1:K
    for m=1:M
        if Gammaa(m,k) >= Gamma_H(k)
            A(m,k) = 1;
        end
    end
end 
 
 
 RReq = RReq_co*log2(1+Req);
 

 
%CVX solve
cvx_status = 'Infeasible';
cvx_quiet true

while ( contains( cvx_status,'Infeasible' ) ) 
    %Parameter Setting: tmp_K, tmp_index, Gamma, BETA, Phi, A
    tmp_K = sum(C);
    tmp_index = zeros(1,tmp_K);
    tmp_Phi = zeros(tau,tmp_K);
%     tmp_A = zeros(M,tmp_K);
    tmp_A = ones(M,tmp_K);
    j = 1;
        
    for k = 1:K
        if C(k) == 1
            tmp_index(j) = k;
            j = j+1;
        end
    end
    
%     for k = 1:tmp_K
%         tmp_A(:,k) = A(:, tmp_index(k) );
%     end
    
    
    for k = 1:tmp_K
        tmp_Phi(:,k) = Phii(:, tmp_index(k) );
    end
    
    BETAAn = zeros(M,tmp_K);
    Gammaan_1 = zeros(M,tmp_K);
    Gammaan_2 = zeros(M,tmp_K);
    
    for k = 1:tmp_K
        BETAAn(:,k) = Pu*BETAA(:, tmp_index(k) );
        Gammaan_1(:,k) = Pu*Gammaa(:, tmp_index(k) );
        Gammaan_2(:,k) = Pu*Gammaa(:, tmp_index(k) );
    end
    
    BETAA_1 = BETAAn/Pu;
    %Interference and Noise Computing
    PhiPhi = zeros(tmp_K,tmp_K);
    
    for ii=1:tmp_K
        for k=1:tmp_K
            PhiPhi(ii,k)=norm(tmp_Phi(:,ii)'*tmp_Phi(:,k));
        end
    end
    
    Te1 =zeros(tmp_K,tmp_K);
    Te2 =zeros(tmp_K,tmp_K);   
    disp(count_1);

    
    for ii=1:tmp_K
        for k=1:tmp_K
            Te1(ii,k)=N*   sum(BETAAn(:,ii).*Gammaan_1(:,k));
            Te2(ii,k)=N^2*(sum((Gammaan_1(:,k)./BETAA_1(:,k)).*BETAA_1(:,ii)) )^2*PhiPhi(k,ii)^2;   
        end
    end
    
    if count_1 >= stop_threshold+1
        check = 1;
        break
    end
    
    
%     while contains(cvx_status, 'Infeasible') && count <= M - min(num)+1 
        if ~isequal( size(Gammaan_2), size(tmp_A) )
            check = 1;
            break
        end
        Gammaan_2 = Gammaan_2.*tmp_A;
        
        cvx_begin %sdp
              variables x(tmp_K,1) 
              minimize(sum(x))
              subject to
                for k=1:tmp_K
                    Te1(:,k)'*x + [Te2(1:(k-1),k); Te2((k+1):tmp_K,k) ]'*[x(1:(k-1)); x((k+1):tmp_K)] + N*sum(Gammaan_1(:,k)) <= N^2*(sum(Gammaan_2(:,k)))^2*x(k)/(2^RReq-1) ;
                end   
                
                for k=1:tmp_K
                    x(k)<=1;
                    x(k)>=0;
                end

                            
         cvx_end
         
         %Update tmp_A
         
%          for k= 1:tmp_K
%              if sum(tmp_A(:,k)) < M
%                  if  num( tmp_index(k) )+count-1 > M
%                      check = 1;
%                      break
%                  end
%                  tmp_A( sort_order( num( tmp_index(k) )+count-1, tmp_index(k) ), tmp_index(k) ) = 1;
%              end
%          end
         
%          count = count+1;
%          tmp_A_check = sum(tmp_A);
%     end
            
         %Update C
         
        if count_1<=stop_threshold
            tmp_x = ones(tmp_K,1);
            
            check_SINR = zeros(1,tmp_K);
            check_status = zeros(1,tmp_K);
            
            for k=1:tmp_K
                check_SINR(k) = tmp_x(k)*sum(Gammaan_2(:,k))^2/( Te1(:,k)'*tmp_x + [Te2(1:(k-1),k); Te2((k+1):tmp_K,k) ]'*[tmp_x(1:(k-1)); tmp_x((k+1):tmp_K)] + sum(Gammaan_1(:,k)) ) ;
                check_status(k) = Te1(:,k)'*tmp_x + [Te2(1:(k-1),k); Te2((k+1):tmp_K,k) ]'*[tmp_x(1:(k-1)); tmp_x((k+1):tmp_K)] + sum(Gammaan_1(:,k)) <= (sum(Gammaan_2(:,k)))^2*tmp_x(k)/(2^RReq-1) ;
            end 
            
            [~,idx_reject] = sort(check_SINR);
            
            C( tmp_index(idx_reject(1)) ) = 0;
            
        end

        
        count_1 = count_1+1;
        count = 1;
        

    
end

rate_buff = zeros(1,tmp_K);

if check == 1
    x = ones(tmp_K,1);
end

if contains(cvx_status, 'Solve')
    check = 0;
end

for k = 1:tmp_K
    rate_buff(k) = log2( 1+x(k)*N^2*sum(Gammaan_2(:,k))^2/( Te1(:,k)'*x + [Te2(1:(k-1),k); Te2((k+1):tmp_K,k) ]'*[x(1:(k-1)); x((k+1):tmp_K)] + sum(Gammaan_1(:,k)) ) );
end

for k=1:tmp_K
    rate( tmp_index(k) ) = rate_buff(k);
end

unSatisfy = K-sum(C);


for k = 1:tmp_K
    A(:, tmp_index(k) ) = tmp_A(:,k);
end
    
num_Satisfy = sum(A);
power_co = mean(x);

end

