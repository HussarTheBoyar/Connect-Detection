function [A] = choosing(rate,RReq)
%M : number of APs
%K : Number of users
%N : Connect Threshold
%RReq: if rate of single AP less than this value, disconnect to the AP
[M,K] = size(rate);
%Choosing matrix A: A(m,k)=1 if connect, =0 if disconnect
A = zeros(M,K);
k =1;
tmp_rate = sort(rate(:,k));
tmp_index = find(tmp_rate>=RReq);

if sum(size(tmp_index) == [1 0]) == 2
    

end

