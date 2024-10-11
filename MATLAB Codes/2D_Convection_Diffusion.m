clear all; clc; format compact; format shortg;

%Physical Constants
rho = 1;
gamma = 1/50;

%Loading the cell face locations
load yc.dat;
load xc.dat;

%Number of cells
Nx = numel(xc) - 1;
Ny = numel(yc) - 1;

%Calculating node positions
xnode = zeros(Nx+2,1);
ynode = zeros(Ny+2,1);

xnode(1) = xc(1); ynode(1) = yc(1);
xnode(end) = xc(end); ynode(end) = yc(end);

for i = 1:Nx
    xnode(i+1) = (xc(i+1) + xc(i))/2;
end

for j = 1:Ny
    ynode(j+1) = (yc(j+1) + yc(j))/2;
end

%deltax and deltay
deltax = zeros(Nx,1);
deltay = zeros(Ny,1);

for i = 1:Nx
    deltax(i) = xc(i+1) - xc(i);
end

for j = 1:Ny
    deltay(j) = yc(j+1) - yc(j);
end

%Delta xe,xw,yn,ys
deltaxe = zeros(Nx,1);
deltaxw = zeros(Ny,1);

for i = 1:Nx
    deltaxe(i) = xnode(i+2) - xnode(i+1);
    deltaxw(i) = xnode(i+1) - xnode(i);
end

deltayn = zeros(Ny,1);
deltays = zeros(Ny,1);

for j = 1:Ny
    deltayn(j) = ynode(j+2) - ynode(j+1);
    deltays(j) = ynode(j+1) - ynode(j);
end

%Loading velocity profile
load u.dat;
load v.dat;

unode = flip(reshape(u,Nx+2,Ny+2));
vnode = flip(reshape(v,Nx+2,Ny+2));

%Calculating the velocities at the faces
uface = zeros(Ny,Nx+1);
vface = zeros(Ny+1,Nx);

uface(:,1) = unode(2:end-1,1);
uface(:,end) = unode(2:end-1,end);

for i = 1:Ny
    for j = 2:Nx
        fx = 0.5*deltax(j)/(xnode(j+1) - xnode(j));
        uface(i,j) = fx*unode(i,j) + (1-fx)*unode(i+1,j);
    end
end

vface(1,:) = vnode(1,2:end-1);
vface(end,:) = vnode(end,2:end-1);

for j = 1:Nx
    for i = 2:Ny
        fy = 0.5*deltay(i)/(ynode(i+1)-ynode(i));
        vface(i,j) = fy*vnode(i,j) + (1-fy)*vnode(i,j+1);
    end
end

ue = uface(:,2:end);
uw = uface(:,1:end-1);

vn = vface(1:end-1,:);
vs = vface(2:end,:);

%Calculating the Peclet Number
Pe_e = zeros(Ny,Nx);
Pe_w = zeros(Ny,Nx);
Pe_n = zeros(Ny,Nx);
Pe_s = zeros(Ny,Nx);

for i = Ny:-1:1
    for j = 1:Nx
        Pe_e(i,j) = rho*ue(i,j)*deltaxe(j)/gamma;
        Pe_w(i,j) = rho*uw(i,j)*deltaxw(j)/gamma;
        Pe_n(i,j) = rho*vn(i,j)*deltayn(Ny+1-i)/gamma;
        Pe_s(i,j) = rho*vs(i,j)*deltays(Ny+1-i)/gamma;
    end
end

%Calculating aE,aW,aN,aS,aP
aE = zeros(Ny,Nx);
aW = zeros(Ny,Nx);
aN = zeros(Ny,Nx);
aS = zeros(Ny,Nx);
aP = zeros(Ny,Nx);

for i = Ny:-1:1
    for j = 1:Nx
        
        %Diffusion Terms
        De = gamma*deltay(Ny+1-i)/deltaxe(j);
        Dw = gamma*deltay(Ny+1-i)/deltaxw(j);
        Dn = gamma*deltax(j)/deltayn(Ny+1-i);
        Ds = gamma*deltax(j)/deltays(Ny+1-i);

        %Convection Terms
        Fe = rho*ue(i,j)*deltay(Ny+1-i);
        Fw = rho*uw(i,j)*deltay(Ny+1-i);
        Fn = rho*vn(i,j)*deltax(j);
        Fs = rho*vs(i,j)*deltax(j);

        %Interpolation functions for stretch mesh
        fxe = deltax(j)*0.5/deltaxe(j);
        fxw = deltax(j)*0.5/deltaxw(j);
        fyn = deltay(Ny+1-i)*0.5/deltayn(Ny+1-i);
        fys = deltay(Ny+1-i)*0.5/deltays(Ny+1-i);

        %Hybrid Scheme
        %We use 1/fxe instead to 2 for Peclet number condition as it is not an equidistant mesh
        %Explicit Peclet number condition is not required
        aE(i,j) = max([0, -Fe, De - fxe*Fe]);
        aW(i,j) = max([0, Fw, Dw + fxw*Fw]);
        aN(i,j) = max([0, -Fn, Dn - fyn*Fn]);
        aS(i,j) = max([0, Fs, Ds + fys*Fs]);
       
        deltaf = Fe - Fw + Fn - Fs;
        aP(i,j) = aE(i,j) + aW(i,j) + aN(i,j) + aS(i,j);
    end
end

%Main Iteration Loop
T = 10*ones(Ny+2,Nx+2);

%Dirichilet Boundary Conditions
T(1:6,1) = 20;
T(1:end-6,end) = 50;

%Neumann Boundary Conditions
T(7:end,1) = T(7:end,2);
T(1,:) = T(2,:);
T(end,:) = T(end-1,:);
T(end-5:end,end) = T(end-5:end,end-1);

err = 1;
tol = 10^-4;
iter = 0;

while err>tol
    err = 0;
    iter = iter + 1;

    % % TDMA Sweep in X direction
    % for i = Ny:-1:1
    %         a = aP(i,1);
    %         b = aE(i,1);
    %         c = aW(i,1);
    %         d = aN(i,1)*T(i,2) + aS(i,1)*T(i+2,2);
    % 
    %         P(1) = b/a;
    %         Q(1) = (c*T(i+1,1) + d)/a;
    % 
    %     for j = 2:1:Nx
    %         a = aP(i,j);
    %         b = aE(i,j);
    %         c = aW(i,j);
    %         d = aN(i,j)*T(i,j+1) + aS(i,j)*T(i+2,j+1);
    % 
    %         P(j) = b/(a-c*P(j-1));
    %         Q(j) = (c*Q(j-1) + d)/(a-c*P(j-1));
    %     end
    % 
    %     for j = Nx:-1:1
    %         T(i+1,j+1) = P(j)*T(i+1,j+2) + Q(j);
    %     end
    % end

    % % TDMA Sweep in Y direction
    % for j = 1:1:Nx
    %         a = aP(1,j);
    %         b = aS(1,j);
    %         c = aN(1,j);
    %         d = aE(1,j)*T(2,j+2) + aW(1,j)*T(2,j);
    % 
    %         P(1) = b/a;
    %         Q(1) = (c*T(1,j+1) + d)/a;
    % 
    %     for i = 2:1:Ny
    %         a = aP(i,j);
    %         b = aS(i,j);
    %         c = aN(i,j);
    %         d = aE(i,j)*T(i+1,j+2) + aW(i,j)*T(i+1,j);
    % 
    %         P(i) = b/(a-c*P(i-1));
    %         Q(i) = (c*Q(i-1) + d)/(a-c*P(i-1));
    %     end
    % 
    %     for i = Ny:-1:1
    %         T(i+1,j+1) = P(i)*T(i+2,j+1) + Q(i);
    %     end
    % end

    %Gauss Sidel Method
    for i = Ny:-1:1
        for j = 1:1:Nx
            T(i+1,j+1) = (aE(i,j)*T(i+1,j+2) + aW(i,j)*T(i+1,j) + aS(i,j)*T(i+2,j+1) + aN(i,j)*T(i,j+1))/aP(i,j);
        end
    end

    for i = Ny:-1:1
        for j = 1:1:Nx
            %Calculating the residual at each iteration
            err = err + abs(aE(i,j)*T(i+1,j+2) + aW(i,j)*T(i+1,j) + aS(i,j)*T(i+2,j+1) + aN(i,j)*T(i,j+1) - aP(i,j)*T(i+1,j+1));
        end
    end
    res(iter) = err;

    %Neumann Boundary Condition
    T(7:end,1) = T(7:end,2);
    T(1,:) = T(2,:);
    T(end,:) = T(end-1,:);
    T(end-5:end,end) = T(end-5:end,end-1);
end

%Plotting Temperature profile
figure(1);
pcolor(xnode,flip(ynode),T)
set(gca, 'YDir', 'normal');
shading('interp')
colormap('jet')
colorbar
xlabel('X')
ylabel('Y')
title('Temperature Distribution')
hcb=colorbar;
title(hcb,'Temp')

% Plotting the mesh
figure(2);
[xc_grid,yc_grid] = meshgrid(xc,yc);
plot(xc_grid,yc','-k')
hold on
plot(xc,yc_grid','-k')

[xnode_grid,ynode_grid] = meshgrid(xnode,ynode);
plot(xnode_grid,ynode','.r')
hold on
plot(xnode,ynode_grid','.r')

title('Mesh')
xlabel('X')
ylabel('Y')

% Plotting the velocity profile
figure(3);
quiver(xnode,flip(ynode),unode,vnode,5);
axis('equal');
title('Velocity Profile')
xlabel('X')
ylabel('Y')

%Calculating and plotting temperatures at x = 1, y = 1
V1 = interp2(xnode,flip(ynode),T,xnode,1);
V2 = interp2(xnode,flip(ynode),T,1,ynode);

figure(4);
plot(xnode,V1,'--r')
hold on
xlabel("X")
ylabel("Temperature")
title("Temperature vs X at y = 1")

figure(5);
plot(ynode,V2,'--r')
hold on
xlabel("Y")
ylabel("Temperature")
title("Temperature vs Y at x  = 1")
