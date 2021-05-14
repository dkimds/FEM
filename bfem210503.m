% Generate right hand side matrix b(=Ax) using gaussian quadrature  
clear
% Source
f = @(x,y) 2*pi*pi*sin(pi*x)*sin(pi*y);
% Quadrature points and weights
xw=[0.33333333333333    0.33333333333333   -0.56250000000000
    0.20000000000000    0.20000000000000    0.52083333333333
    0.60000000000000    0.20000000000000    0.52083333333333
    0.20000000000000    0.60000000000000    0.52083333333333];

% Generate mesh
h=1/2^3;
x=0:h:1;
M = length(x);
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

% RHS b
b=zeros(length(N),1);
for n = 1:length(T)
    k = T(n,:);
    m=[N(k(1),1) N(k(1),2) 1;
       N(k(2),1) N(k(2),2) 1;
       N(k(3),1) N(k(3),2) 1];
    % Original Triangle
    x1 = m(1,1); 
    y1 = m(1,2);
    x2 = m(2,1); 
    y2 = m(2,2);
    x3 = m(3,1); 
    y3 = m(3,2);
    % basis functions
    p1 = m\([1 0 0]');
    p2 = m\([0 1 0]');
    p3 = m\([0 0 1]');
    % Integration
    area_tri = abs(det(m))*0.5;
    qd1 =0;
    qd2 =0;
    qd3 =0;
    for i = 1:4
        s = x1 + (x2-x1)*xw(i,1) + (x3-x1)*xw(i,2);
        t = y1 + (y2-y1)*xw(i,1) + (y3-y1)*xw(i,2);
        qd1 = qd1 + area_tri*xw(i,3)*f(s,t)*(p1(1)*s + p1(2)*t + p1(3));
        qd2 = qd2 + area_tri*xw(i,3)*f(s,t)*(p2(1)*s + p2(2)*t + p2(3));
        qd3 = qd3 + area_tri*xw(i,3)*f(s,t)*(p3(1)*s + p3(2)*t + p3(3));
    end  
    % Assembly
    b(k(1)) = b(k(1)) + qd1;
    b(k(2)) = b(k(2)) + qd2;
    b(k(3)) = b(k(3)) + qd3;
end

% Remove basis ftns on boundary
B = [];
for k = 1:length(N)
    x = N(k,1);
    y = N(k,2);
    
    if x==0
        B=[B k];
    elseif x == 1
        B=[B k];
    elseif y == 0
        B=[B k];
    elseif y == 1
        B=[B k];
    end
end
B = unique(B, 'rows');
b(B) = [];





















