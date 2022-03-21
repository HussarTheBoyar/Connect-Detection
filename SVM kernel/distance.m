function res = distance(matrix, AP)
[X,~] = size(matrix);
L = length(AP);
res = zeros(L,X);
for i = 1:L
   for j = 1:X
       res(i,j) = norm(matrix(j,:)-AP(i,:));
   end
end

end