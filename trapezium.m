a1=[1 0]; a2=[0.5 sqrt(3)/2];
x =zeros(100);
y=zeros(100);
for i = 1:100
    for j = 1:100
        x(i,j)=a1(1)*i/100+a2(1)*j/100;
        y(i,j)=a1(2)*i/100+a2(2)*j/100;
    end
end

daspect([1 1 1])
plot(x,y,'o')
axis equal
z=cos(x*3)+cos(y*3);
daspect([1 1 1])
figure
daspect([1 1 1])
surf(x,y,z)
daspect([1 1 1])