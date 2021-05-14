h = 1/10;
x=0:h:1; M=length(x);

% 2-lines
col = repmat(x, 1, M);
row = reshape(repmat(x, M, 1), [numel(col), 1]);
plot(row, col, ".")

N = [row, col'];
T=delaunay(N);
trimesh(T, N(:,1),N(:,2))

% stiffness matrix
A = sparse(11);
for m = 1:length(T)
    k=T(m,:);
    M=[N(k(1),1) N(k(1),2) 1;
       N(k(2),1) N(k(2),2) 1;
       N(k(3),1) N(k(3),2) 1];
    p1 = M\([1 0 0]');
    p2 = M\([0 1 0]');
end
