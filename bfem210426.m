% Integrate sin(pi*x)*sin(pi*y) on [0,1]x[0,1]
clear
% Integrand
g = @(x,y) sin(pi*x)*sin(pi*y);
% Quadrature points and weights
xw=[0.33333333333333    0.33333333333333   -0.56250000000000
    0.20000000000000    0.20000000000000    0.52083333333333
    0.60000000000000    0.20000000000000    0.52083333333333
    0.20000000000000    0.60000000000000    0.52083333333333];
qd = 0; % Initial value

% Generate mesh
h=1/2^5;
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

% Summation on the whole domain
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
    % Summation on a triangle
    for i=1:4
        s = x1 + (x2-x1)*xw(i,1) + (x3-x1)*xw(i,2);
        t = y1 + (y2-y1)*xw(i,1) + (y3-y1)*xw(i,2);
        qd = qd + abs(det(m))*0.5*xw(i,3)*g(s,t);
    end
end

% Display error
err = abs((2/pi)^2-qd)