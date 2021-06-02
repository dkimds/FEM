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

B=[];
for n = 1:length(T)
    k = T(n,:);
    x1 = N(k(1),1);
    y1 = N(k(1),2);
    x2 = N(k(2),1);
    y2 = N(k(2),2);
    x3 = N(k(3),1);
    y3 = N(k(3),2);
    
    if x1 > 0.5 && y1 > 0.5
        B=[B n];
    end
    if x2 > 0.5 && y2 > 0.5
        B=[B n];
    end
    if x3 > 0.5 && y3 > 0.5
        B=[B n];
    end
end

B = unique(B);

T(B,:)=[];
% Visualize mesh
trimesh(T,N(:,1),N(:,2))
