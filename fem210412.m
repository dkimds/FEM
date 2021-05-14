h = 1/4;
x=0:h:1; M=length(x);

% 2-lines
col = repmat(x, 1, M);
row = reshape(repmat(x, M, 1), [numel(col), 1]);
plot(row, col, ".")

N = [row, col'];

% Save indice of boundary points
B= [];

for k = 1:length(N)
    x = N(k,1);
    y= N(k,2);
    
    if x==0
        B=[B k];
    elseif x ==1
        B=[B k];
    elseif y == 0
        B=[B k];
    elseif y==1
        B=[B k];
    end
end

B = unique(B, 'rows');

% matrix 
f = @(x, y) 2* pi * pi * sin(pi*x) * sin(pi*x*y);
RHS = f(node
