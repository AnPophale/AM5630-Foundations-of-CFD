clc; clear all; close all; format shortg; format compact;

%Domain dimensions
L = 1;
H = 0.5;

%Number of cells
Nx = 20;
Ny = 20;

%Delta x and delta y
deltax = zeros(1,Nx);
deltay = zeros(1,Ny);

deltax(1) = L/Nx;
deltay(1) = H/Ny;

%Stretch factors
Sx = 1;
Sy = 1;

for i = 2:Nx
    deltax(i) = deltax(i-1)*Sx; 
end
for j = 2:Ny
    deltay(j) = deltay(j-1)*Sy;
end

%Resizing
deltax = deltax/sum(deltax);
deltay = deltay/(2*sum(deltay));

%Cell Face Locations
xface = zeros(1,Nx+1);
yface = zeros(1,Ny+1);

xface(1) = 0; yface(1) = 0;

for i = 2:Nx+1
    xface(i) = xface(i-1) + deltax(i-1);
end
for j = 2:Ny+1
    yface(j) = yface(j-1) + deltay(j-1);
end

%Node Loactions
xnode = zeros(1,Nx+2);
ynode = zeros(1,Ny+2);

xnode(1) = 0; xnode(end) = L;
ynode(1) = 0; ynode(end) = H;

for i = 2:Nx+1
    xnode(i) = (xface(i-1) + xface(i))/2;
end
for j = 2:Ny+1
    ynode(j) = (yface(j-1) + yface(j))/2;
end

%Delta xe,xw,ye,yw
deltaxe = zeros(1,Nx);
deltaxw = zeros(1,Ny);

for i = 1:Nx
    deltaxe(i) = xnode(i+2) - xnode(i+1);
    deltaxw(i) = xnode(i+1) - xnode(i);
end

deltayn = zeros(1,Ny);
deltays = zeros(1,Ny);

for j = 1:Ny
    deltayn(j) = ynode(j+2) - ynode(j+1);
    deltays(j) = ynode(j+1) - ynode(j);
end

%k at nodes
knode = zeros(Ny+2,Nx+2);
for i = Ny+2:-1:1
    knode(i,:) = 16*(1+ynode(Ny+3-i)/H);
end

%k at faces using interpolation function
ke = zeros(Ny,Nx);
kw = zeros(Ny,Nx);
kn = zeros(Ny,Nx);
ks = zeros(Ny,Nx);

for i = Ny:-1:1
    for j = 1:Nx

        if j == Nx
            ke(i,j) = knode(i+1,j+2);
        else
            fxe = deltax(j)*0.5/deltaxe(j);
            ke(i,j) = knode(i+1,j+2)*fxe + (1-fxe)*knode(i+1,j+1);
        end

        if j == 1
            kw(i,j) = knode(i+1,j);
        else    
            fxw = deltax(j)*0.5/deltaxw(j);
            kw(i,j) = knode(i+1,j)*fxw + (1-fxw)*knode(i+1,j+1);
        end

        if i == 1
            kn(i,j) = knode(i,j);
        else
            fyn = deltay(Ny+1-i)*0.5/deltayn(Ny+1-i);
            kn(i,j) = fyn*knode(i,j) + (1-fyn)*knode(i+1,j);
        end

        if i == Ny
            ks(i,j) = knode(i+2,j);
        else
            fys = deltay(Ny+1-i)*0.5/deltays(Ny+1-i);
            ks(i,j) = fys*knode(i+2,j) + (1-fys)*knode(i+1,j);
        end
            
    end
end

%Dirichilet Boundary Conditions
Te = zeros(Ny,1);
for i = Ny:-1:1
    Te(i) = 5*(1-ynode(Ny+2-i)/H) + 15*sin(pi*ynode(Ny+2-i)/H);
end

Tn =10*ones(Nx+2,1);
Ts = 15*ones(Nx+2,1);

%Calculating aE,aW,aN,aS,aP,Su using Implicit BC implementation
aE = zeros(Ny,Nx);
aW = zeros(Ny,Nx);
aN = zeros(Ny,Nx);
aS = zeros(Ny,Nx);
aP = zeros(Ny,Nx);
Su = zeros(Ny,Nx);

for i = Ny:-1:1
    for j = 1:Nx

        if j == Nx
            aE(i,j) = 0;
            aP(i,j) = aP(i,j) + ke(i,j)*deltay(Ny+1-i)/deltaxe(j);
            Su(i,j) = Su(i,j) + Te(i)*ke(i,j)*deltay(Ny+1-i)/deltaxe(j);
        else
            aE(i,j) = ke(i,j)*deltay(Ny+1-i)/deltaxe(j);
        end

        if j == 1
            aW(i,j) = 0;
        else
            aW(i,j) = kw(i,j)*deltay(Ny+1-i)/deltaxw(j);
        end

        if i == 1
            aN(i,j) = 0;
            aP(i,j) = aP(i,j) + kn(i,j)*deltax(j)/deltayn(Ny+1-i);
            Su(i,j) = Su(i,j) + Tn(j)*kn(i,j)*deltax(j)/deltayn(Ny+1-i);
        else
            aN(i,j) = kn(i,j)*deltax(j)/deltayn(Ny+1-i);
        end

        if i == Ny
            aS(i,j) = 0;
            aP(i,j) = aP(i,j) + ks(i,j)*deltax(j)/deltays(Ny+1-i);
            Su(i,j) = Su(i,j) + Ts(j)*ks(i,j)*deltax(j)/deltays(Ny+1-i);
        else
            aS(i,j) = ks(i,j)*deltax(j)/deltays(Ny+1-i);
        end

        aP(i,j) = aP(i,j) + aE(i,j) + aW(i,j) + aN(i,j) + aS(i,j);
    end
end

%Gauss Sidel Iteration loop
T = ones(Ny+2,Nx+2);
T(2:end-1,end) = Te;
T(1,:) = Tn;
T(end,:) = Ts;

S = -1.5;
err = 1;
tol = 10^-1;
iter = 0;

while err>tol
    err = 0;
    iter = iter + 1;
    Told = T;
    for i = Ny:-1:1
        for j = 1:1:Nx
            T(i+1,j+1) = (aE(i,j)*T(i+1,j+2) + aW(i,j)*T(i+1,j) + aS(i,j)*T(i+2,j+1) + aN(i,j)*T(i,j+1) + Su(i,j))/(aP(i,j) - S*deltax(j)*deltay(Ny+1-i)/Told(i+1,j+1));
        end
    end
    for i = Ny:-1:1
        for j = 1:1:Nx
            %Calculating the residual at each iteration
            err = err + abs(aE(i,j)*T(i+1,j+2) + aW(i,j)*T(i+1,j) + aS(i,j)*T(i+2,j+1) + aN(i,j)*T(i,j+1) + Su(i,j) - (aP(i,j) - S*deltax(j)*deltay(Ny+1-i)/Told(i+1,j+1))*T(i+1,j+1));
        end
    end
    res(iter) = err;
    %Neumann Boundary Condition
    T(:,1) = T(:,2);
end

%Plotting the mesh
figure(1)
[xface_grid,yface_grid] = meshgrid(xface,yface);
[xnode_grid,ynode_grid] = meshgrid(xnode,ynode);

plot(xface_grid,yface','-r')
hold on
plot(xface,yface_grid','-r')
hold on

plot(xnode_grid,ynode','.k')
hold on
plot(xnode,ynode_grid','.k')
xlabel("X")
ylabel("Y")
title("Mesh")

%Plotting Temperature profile
figure(2);
pcolor(xnode,flip(ynode),T)
set(gca, 'YDir', 'normal');
shading('interp')
colormap('jet')
colorbar
xlabel('X')
ylabel('Y')
title('Temperature Distribution (20x20 mesh)')
hcb=colorbar;
title(hcb,'Temp')

%Calculating and plotting temperatures at x = 0.5, y = 0.25
V1 = interp2(xnode,flip(ynode),T,xnode,0.25);
V2 = interp2(xnode,flip(ynode),T,0.5,ynode);

figure(3);
plot(xnode,V1,'-')
hold on
xlabel("X")
ylabel("Temperature")
title("Temperature vs X at y = 0.25")

figure(4);
plot(ynode,V2,'-')
hold on
xlabel("Y")
ylabel("Temperature")
title("Temperature vs Y at x = 0.5")

%Residual vs Number of Iterations
figure(5)
semilogy(1:1:iter,res)
xlabel("Number of Iterations")
ylabel("Residual")
title("Residual vs Number of Iterations (Mesh size 20x20)")
hold on

%Error Tolerance vs Number of Iterations
figure(6)
semilogy(iter,tol,'-o')
xlabel("Error Tolerance")
ylabel("Number of iterations")
title("Error Tolerance vs Number of Iterations (Mesh Size 20x20)")
hold on

%Plotting the heat flux
deltax2 = deltaxw;
deltax2(end+1) = deltaxe(end);

kxface = kw;
kxface(:,end+1) = ke(:,end);

for i = Ny:-1:1
    for j = 1:Nx+1
        qx(i,j) = kxface(i,j)*(T(i+1,j+1) - T(i+1,j))/deltax2(j);
    end
end

deltay2 = deltayn;
deltay2(end+1) = deltays(end);

kyface = kn;
kyface(end+1,:) = ks(end,:);

for i = Ny+1:-1:1
    for j = 1:Nx
        qy(i,j) = kyface(i,j)*(T(i,j+1) - T(i+1,j+1))/deltay2(i);
    end
end

figure(7)
plot(xface_grid,yface','-k')
hold on
plot(xface,yface_grid','-k')
hold on
quiver(xface,ynode(2:end-1),qx,zeros(size(qx)),0.5,'-r',"LineWidth",0.75);
hold on
quiver(xnode(2:end-1),yface,zeros(size(qy)),qy,0.5,'-b',"LineWidth",0.75);
xlim([0,1.05]);
ylim([-0.05,0.5]);
xlabel("X")
ylabel("Y")
title("Heat flux through cell faces in X and Y direction ")
