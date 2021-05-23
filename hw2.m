h = 1/10;
x=0:h:1; M=length(x);
for i =1:M
      N(M*(i-1)+1: M*i, 1:2)=[(i-1)*h*ones(M,1) x'];
end

%%
nx = length(x);
T=[];
nx = length(x);
T=[];
for m = 1:nx-1
    for k =1:nx-1
        T = [T; nx*(m-1)+k nx*(m-1)*k+1 nx*m+k+1];
        T = [T; nx*(m-1)+k nx*m+k+1 nx*m+k];
    end
end
%%
trimesh(T, N(:,1), N(:,2))