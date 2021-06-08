% Solve the inhomogeneous Dirichlet BVP by seperating final solution to boundary value and interior
% value.
clear
% Source
f = @(x,y) -2*pi*pi*cos(pi*x)*cos(pi*y);
% Quadrature points and weights
xw=[0.33333333333333    0.33333333333333   -0.56250000000000
    0.20000000000000    0.20000000000000    0.52083333333333
    0.60000000000000    0.20000000000000    0.52083333333333
    0.20000000000000    0.60000000000000    0.52083333333333];

% Generate mesh
h = 1/8;
x = 0:h:1;
M = length(x);
col = repmat(x',M,1);
row = reshape(repmat(x,M,1),[numel(col),1]);
N=[row, col];

% triangluation
T = [];
for j = 1:M-1
    for i = 1:M-1
        % counter-clockwise
        T = [T;
            i+M*(j-1) i+M*j i+1+M*(j-1); 
            i+M*j i+1+M*j i+1+M*(j-1)];
    end
end

%stiffness matrix A and RHS b
A =sparse(M^2,M^2);
b=zeros(length(N),1);

for n = 1:length(T)
    k = T(n,:);
    m=[N(k(1),1) N(k(1),2) 1;
       N(k(2),1) N(k(2),2) 1;
       N(k(3),1) N(k(3),2) 1];
    % basis functions   
    p1 = m\([1 0 0]');
    p2 = m\([0 1 0]');
    p3 = m\([0 0 1]');
    % Area of a triangle
    S = abs(det(m))*0.5; 

    % Contruct stiff matrix
    % element stiff matrix can be computed b.c. basis ftns are 1st order.
    p11 = p1(1:2)'*p1(1:2)*S;
    p22 = p2(1:2)'*p2(1:2)*S;
    p33 = p3(1:2)'*p3(1:2)*S;
    p12 = p1(1:2)'*p2(1:2)*S;
    p23 = p2(1:2)'*p3(1:2)*S;
    p13 = p1(1:2)'*p3(1:2)*S;
    p21 = p12;
    p32 = p23;
    p31 = p13;
    % Assembly
    A(k(1),k(1)) = A(k(1),k(1)) + p11;
    A(k(2),k(2)) = A(k(2),k(2)) + p22;
    A(k(3),k(3)) = A(k(3),k(3)) + p33;
    A(k(1),k(2)) = A(k(1),k(2)) + p12;
    A(k(2),k(3)) = A(k(2),k(3)) + p23;
    A(k(1),k(3)) = A(k(1),k(3)) + p13;
    A(k(2),k(1)) = A(k(2),k(1)) + p21;
    A(k(3),k(2)) = A(k(3),k(2)) + p32;
    A(k(3),k(1)) = A(k(3),k(1)) + p31;

    % Construct RHS
    % Original Triangle
    x1 = m(1,1); 
    y1 = m(1,2);
    x2 = m(2,1); 
    y2 = m(2,2);
    x3 = m(3,1); 
    y3 = m(3,2);  
    % Integration w/ Gaussian quadrature on reference triangle
    qd1 =0;
    qd2 =0;
    qd3 =0;
    for i = 1:4
        % calculate quadrature pts on ref triangle
        s = x1 + (x2-x1)*xw(i,1) + (x3-x1)*xw(i,2);
        t = y1 + (y2-y1)*xw(i,1) + (y3-y1)*xw(i,2);
        % Gaussian quadrature
        qd1 = qd1 + S*xw(i,3)*f(s,t)*(p1(1)*s + p1(2)*t + p1(3));
        qd2 = qd2 + S*xw(i,3)*f(s,t)*(p2(1)*s + p2(2)*t + p2(3));
        qd3 = qd3 + S*xw(i,3)*f(s,t)*(p3(1)*s + p3(2)*t + p3(3));
    end  
    % Assembly
    b(k(1)) = b(k(1)) + qd1;
    b(k(2)) = b(k(2)) + qd2;
    b(k(3)) = b(k(3)) + qd3;  
end

% Find indice of bdr points and evaluate values of theirs.
BDR = [];
ub=zeros(length(N),1);

for k = 1:length(N)
    x = N(k,1);
    y = N(k,2);
    
    if x==0
        BDR=[BDR k];
        ub(k) = cos(pi*y);
    elseif x == 1
        BDR=[BDR k];
        ub(k) = -cos(pi*y);
    elseif y == 0
        BDR=[BDR k];
        ub(k) = cos(pi*x);
    elseif y == 1
        BDR=[BDR k];
        ub(k) = -cos(pi*x);
    end
end
BDR = unique(BDR, 'rows');

% Subtract RHS(A*ub) on b to find the interior solution
b = b-A*ub;
final_sol = zeros(length(N),1);
% Find indice of the interior points
all = 1:length(N);
INT = setdiff(all, BDR);
% Solve the equation wrt only interior pts
A(BDR,:)=[];
A(:,BDR)=[];
b(BDR)=[];
z=A\b;
final_sol(INT)=z;
final_sol = final_sol + ub;

% Error Analysis
u = @(x,y) cos(pi*x).*cos(pi*y);
u_ext = u(N(:,1),N(:,2));

err_l1 = max(final_sol - u_ext);

err_l2 = sqrt(sum((final_sol - u_ext).^2));



























































