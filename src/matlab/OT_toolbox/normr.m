function out = normr(X)
% normalize each row of matrix X

nrows = length(X(:,1));
ncols = length(X(1,:));
out = zeros(nrows,ncols);
for i = 1:nrows
    row = X(i,:);
    normfactor = sqrt(sum(X(i,:).^2));
    out(i,:) = row/normfactor;
end

end
