h = 1/10;
x=0:h:1; M=length(x);
for i =1:M
      N(M*(i-1)+1: M*i, 1:2)=[(i-1)*h*ones(M,1) x'];
end

plot(N(:,1),N(:,2),".")
