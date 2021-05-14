h = 1/10;
x=(0:h:1)';
for i =1:length(x)
    for j = 1:length(x)
        N(length(x)*(i-1)+j, 1:2)=[(i-1)*h (j-1)*h];
    end
end

plot(N(:,1),N(:,2),".")
