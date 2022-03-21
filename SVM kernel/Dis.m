function res = Dis(d1,dk)
[~,tau] = size(dk);
res = zeros(1,tau);
for i = 1:tau
   res(i) = norm(d1-dk(:,i))^2;
end

end
