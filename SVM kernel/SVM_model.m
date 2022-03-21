function [Mdl] = SVM_model(Rate,Gamma,BETAA,A_max)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[M,K] = size(A_max);
class = [];
data = [];

Gamma = Gamma./max(Gamma);
Rate  = Rate./max(Rate);
BETAA = BETAA./max(BETAA);

for k=1:K
    class = [class; A_max(:,k)];
    for m = 1:M
        tmp_data = [Rate(m,k) Gamma(m,k) BETAA(m,k)];
        data = [data; tmp_data];
    end
end

class = 2*class - ones( size(class) );

%Train the SVM Classifier

Mdl = fitcsvm(data,class,'KernelFunction','rbf',...
    'BoxConstraint',Inf,'ClassNames',[-1,1]);


end

