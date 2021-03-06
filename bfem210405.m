% Construct the stiffness matrix A on 2D
clear
h=1/4;
x=0:h:1;
M = length(x);

% mesh
col = repmat(x',M,1);
row = reshape(repmat(x,M,1),[numel(col),1]);
N=[row, col];

% triangluation
T = [];
for j = 1:M-1
    for i = 1:M-1
        T = [T;
            i+M*(j-1) i+M*j i+1+M*(j-1);
            i+M*j i+1+M*j i+1+M*(j-1)];
    end
end

% % Visualize mesh
% trimesh(T,N(:,1),N(:,2))

% stiffness matrix A
A =sparse(M^2,M^2);
for n = 1:length(T)
    k = T(n,:);
    % gradient
    m=[N(k(1),1) N(k(1),2) 1;
       N(k(2),1) N(k(2),2) 1;
       N(k(3),1) N(k(3),2) 1];
    p1 = m\([1 0 0]');
    p2 = m\([0 1 0]');
    p3 = m\([0 0 1]');
    % element stiff matrix
    p11=p1(1:2)'*p1(1:2)*abs(det(m))/2;
    p22=p2(1:2)'*p2(1:2)*abs(det(m))/2;
    p33=p3(1:2)'*p3(1:2)*abs(det(m))/2;
    p12=p1(1:2)'*p2(1:2)*abs(det(m))/2;
    p23=p2(1:2)'*p3(1:2)*abs(det(m))/2;
    p13=p1(1:2)'*p3(1:2)*abs(det(m))/2;
    p21=p12;
    p32=p23;
    p31=p13;
    % Assembly
    A(k(1),k(1))=A(k(1),k(1))+p11;
    A(k(2),k(2))=A(k(2),k(2))+p22;
    A(k(3),k(3))=A(k(3),k(3))+p33;
    A(k(1),k(2))=A(k(1),k(2))+p12;
    A(k(2),k(3))=A(k(2),k(3))+p23;
    A(k(1),k(3))=A(k(1),k(3))+p13;
    A(k(2),k(1))=A(k(2),k(1))+p21;
    A(k(3),k(2))=A(k(3),k(2))+p32;
    A(k(3),k(1))=A(k(3),k(1))+p31;
end

spy(A)