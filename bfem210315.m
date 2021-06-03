h = 1/10;
x=0:h:1; M=length(x);

% Double for-loops
for i =1:length(x)
    for j = 1:length(x)
        N(length(x)*(i-1)+j, 1:2)=[(i-1)*h (j-1)*h];
    end
end
plot(N(:,1),N(:,2),".")

% Single for-loops
for i =1:M
    N1(M*(i-1)+1: M*i, 1:2)=[(i-1)*h*ones(M,1) x'];
end
plot(N1(:,1),N1(:,2),".")

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