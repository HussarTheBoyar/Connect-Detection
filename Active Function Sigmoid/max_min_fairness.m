function [notSatis,Satis,average_rate,min_rate] = max_min_fairness(Gammaa,BETAA,Pu,Phi,Req_co,Req,N)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[M,K] = size(Gammaa);
tau = 10;
RReq = Req_co*log2(1+Req);
rate = zeros(1,K);
Pp = Pu;

%Pilot contamination
Phii_cf = Phi;
PC = zeros(K,K);
for ii=1:K
    for k=1:K
        PC(ii,k) = sum( (Gammaa(:,k)./BETAA(:,k)).*BETAA(:,ii)  )*Phii_cf(:,k)'*Phii_cf(:,ii);
    end
end
PC1=N^2*(abs(PC)).^2;

for k=1:K
    deno1=0;
    for m=1:M
        deno1=deno1 + Gammaa(m,k)*sum(BETAA(m,:));
    end
    SINR(k) = Pu*(sum(Gammaa(:,k)))^2/(sum(Gammaa(:,k)) + Pu*deno1*N + Pu*sum(PC1(:,k)) - Pu*PC1(k,k));
    %Rate:
    R_cf(k) = log2(1+ SINR(k));
end

R_cf_min = min(R_cf);

% stepp=5;
% Ratestep=zeros(stepp,K);
% Ratestep(1,:)=R_cf;
% 
% %% Find the pilot sequences using greedy method
% for st=2:stepp
%     [minvalue minindex]=min(Ratestep(st-1,:));
%     
%         Mat=zeros(tau,tau)-Pu*sum(BETAA(:,minindex))*Phii_cf(:,minindex)*Phii_cf(:,minindex)';
%         for kk=1:K
%                 Mat = Mat + Pu*sum(BETAA(:,kk))*Phii_cf(:,kk)*Phii_cf(:,kk)';
%         end
%         [U1,S1,V1] = svd(Mat);
%         Phii_cf(:,minindex) = U1(:,tau);
%     
%  
%   %% Create Gamma matrix
% Gammaa = zeros(M,K);
% mau=zeros(M,K);
% for m=1:M
%     for k=1:K
%         mau(m,k)=norm( (BETAA(m,:).^(1/2)).*(Phii_cf(:,k)'*Phii_cf))^2;
%     end
% end
% 
% for m=1:M
%     for k=1:K
%         Gammaa(m,k)=tau*Pp*BETAA(m,k)^2/(tau*Pp*mau(m,k) + 1);
%     end
% end
% 
% %% Compute Rate
% 
% SINR=zeros(1,K);
% R_cf=zeros(1,K);
% 
% %Pilot contamination
% PC = zeros(K,K);
% for ii=1:K
%     for k=1:K
%         PC(ii,k) = sum( (Gammaa(:,k)./BETAA(:,k)).*BETAA(:,ii)  )*Phii_cf(:,k)'*Phii_cf(:,ii);
%     end
% end
% PC1=N^2*(abs(PC)).^2;
% 
% for k=1:K
%     deno1=0;
%     for m=1:M
%         deno1=deno1 + Gammaa(m,k)*sum(BETAA(m,:));
%     end
%     SINR(k) = N^2*Pu*(sum(Gammaa(:,k)))^2/(sum(Gammaa(:,k)) + Pu*deno1*N + Pu*sum(PC1(:,k)) - Pu*PC1(k,k));
%     %Rate:
%     Ratestep(st,k) = log2(1+ SINR(k));
% end
% 
%        
% end
% 
% R_cf_min=min(Ratestep(stepp,:));
% R_cf_user=Ratestep(stepp,:);

%% 2) Max-Min power allocations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmin=2^R_cf_min-1;
tmax=2^(2*R_cf_min+1.2)-1;
epsi=max(tmin/5,0.01);

BETAAn=BETAA*Pu;
Gammaan=Gammaa*Pu;
PhiPhi=zeros(K,K);
Te1 =zeros(K,K);
Te2 =zeros(K,K);
for ii=1:K
    for k=1:K
        PhiPhi(ii,k)=norm(Phii_cf(:,ii)'*Phii_cf(:,k));
    end
end

for ii=1:K
    for k=1:K
        Te1(ii,k)=N*    sum(BETAAn(:,ii).*Gammaan(:,k));
        Te2(ii,k)=N^2* (sum((Gammaan(:,k)./BETAA(:,k)).*BETAA(:,ii)) )^2*PhiPhi(k,ii)^2;   
    end
end

loop_count = 0;

%cvx_solver sedumi
cvx_quiet true
            while( tmax - tmin > epsi)
                loop_count = loop_count+1;

            tnext = (tmax+tmin)/2; 
           cvx_begin %sdp
              variables x(K,1) 
              minimize(0)
              subject to
                for k=1:K
                    Te1(:,k)'*x + [Te2(1:(k-1),k); Te2((k+1):K,k) ]'*[x(1:(k-1)); x((k+1):K)] + N*sum(Gammaan(:,k)) <= N^2*(1/tnext)*(sum(Gammaan(:,k)))^2*x(k) ;
                end              
                for k=1:K
                    x(k)>=0;
                    x(k)<=1;
                end
                            
            cvx_end
            
            if contains(cvx_status,'Solved')
                tmp_x = x;
            end
            
            % bisection
            if strfind(cvx_status,'Solved') % feasible
            fprintf(1,'Problem is feasible ',tnext);
            tmin = tnext;
            else % not feasible
            fprintf(1,'Problem not feasible ',tnext);
            tmax = tnext;   
            end

            end

 for k = 1:K
     rate(k) = N^2*(sum(Gammaan(:,k)))^2*tmp_x(k)/( Te1(:,k)'*tmp_x + [Te2(1:(k-1),k); Te2((k+1):K,k) ]'*[tmp_x(1:(k-1)); tmp_x((k+1):K)] + sum(Gammaan(:,k)) );
 end
 
 notSatis = sum( rate<RReq );
 Satis = K-notSatis;
 average_rate = mean(rate);
 min_rate = min(rate);
end

