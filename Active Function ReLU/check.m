function res = check(dp, muy)
[~,M] = size(muy);
res = zeros(1,M);
for i= 1:M
    res(i) = norm(dp-muy(:,i))^2;
end
end