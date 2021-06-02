h = 1/10;
x=0:h:1; M=length(x);

% 2-lines
col = repmat(x, 1, M);
row = reshape(repmat(x, M, 1), [numel(repmat(x, M, 1)), 1]);
plot(row, col, ".")

% A matrix
A = sparse(M, M);

B = diag(ones(1,M));
C = diag(-ones(1,M-1),1);
D = diag(-ones(1,M-1),-1);

A = (2*B+C+D)/h;