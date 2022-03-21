function [res, Phii, ind, np, loop] = Kmeans_Bjonson_new(U, Ter, AP)

%%Training phase

D = 1;
[K,~] = size(Ter);
Kp = unifrnd(-D/2, D/2, K, 2);
tau = 10;
Phii = zeros(tau,0);

if mod(K,tau) ~= 0
    num = floor(K/tau)+1;
else
    num = floor(K/tau);
end

L = length(AP);
muy_setup = unifrnd(-D/2, D/2, num, 2);
muy_num = zeros(1,floor(K/tau));
index   = zeros(1,K);
epsilon = 0.01;
error   = 1;
muy     = zeros(L,K);
res = zeros(0,2);
loop = 1;

%Set up muy
muy = distance(muy_setup, AP);
dp  = distance(Kp, AP);
    


while error>= epsilon             %convergence conditions
    muy_num = zeros(1,num);
    tmp_muy = muy;
    for i = 1:K
        test = check(dp(:,i), muy);
        [~,idx] = min(test);
        muy_num(idx) = muy_num(idx)+1;
        index(i)  = idx;
    end
    
    %update muy
    for z = 1:num
        a = find(index==z);
        S = 0;
        for j = 1:muy_num(z)
            S = S+dp(:,a(j));
        end
        muy(:,z) = S/muy_num(z);
    end
    error = max(norm(muy-tmp_muy)^2);
    loop = loop+1;
end

%%Checking phase

muy_num = zeros(1,num);
index   = zeros(1,K);
dk = distance(Ter,AP);

for i=1:num
    j = 1;
    sz = tau;
    tmp = dk;
    while j<=tau && sz>=0
        [~,indx] = min(check(muy(:,i),tmp));
        index(indx) = i;
        muy_num(i) = muy_num(i)+1;
        j = j+1;
        sz = sz-1;
        tmp(:,i) = [];
    end
end
    
ind = index;
np = muy_num;
tmp_dk = dk;

%%Reorder users and distance vector 
cen_Ter = zeros(tau,2,num);
for i = 1:num
    a = find(index == i);
    for j = 1:length(a)
        cen_Ter(j,:,i) = Ter(a(j),:);
    end
    res = [res; cen_Ter(:,:,i)];
end



%Pilot assignment
cen_dk = zeros(L,tau,num);
cen_phii = zeros(tau,tau,num);
cen_D = zeros(tau,tau,num);

for i = 1:num
    cen_dk(:,:,i) = dk(:,(i-1)*tau+1:i*tau);
end

cen_phii(:,:,1) = U;

%Create distance matrix to choose used sharing pilot sequence 
for i = 1:num
    for j = 1:tau
        cen_D(j,:,i) = Dis(cen_dk(j,:,1),cen_dk(:,:,i));
    end
end



%Choosing users sharing pilot sequence
for i = 2:num
   cen_d = cen_D(:,:,i);
   [row,~] = size(cen_d);
   cen_tmp = cen_d;
   while row>0
       
       
       
       op = max(max(cen_tmp));
       [op_row, op_col] = find(cen_tmp == op);
       if sum(size(op_row) == [1 1])>2||sum(size(op_col) == [1 1])>2
           [op_row_fix,~] = find(cen_D(:,:,i) == cen_tmp(op_row,:));
       else
           op_row = op_row(1);
           op_col = op_col(1);
           [op_row_fix,~] = find(cen_D(:,:,i) == cen_tmp(op_row,:));
       end
       
       if sum(cen_phii(:,op_col,i) == zeros(tau,1)) == tau
           cen_phii(:,op_col,i) = U(:,op_row_fix(1));
           cen_d(op_row,:) = [];
           row = row-1;
           cen_tmp = cen_d;
       else
           cen_tmp(op_row, op_col) = 0;
       end
       
       
   end
end

for i = 1:num
    Phii = [Phii cen_phii(:,:,i)];
end

end
        
        
        

