% Solve the Neuman BVP by the mean information
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

% expand the stiff matrix with mean information
A = [A;(ones(1,length(N)))];
A = [A (ones(length(N)+1,1))];
% reduce the condition number
A(end,end) = 1e-3;
% expand RHS with mean information
b = [b;0];
x = A\b;
condest(A)

