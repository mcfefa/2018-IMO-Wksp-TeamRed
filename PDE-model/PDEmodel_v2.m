function [u] = PDEmodel_v2()
close all;

%% Initialization - Boundary conditions, time step and # of time steps (usually large)
Nx = [51;51];                    % Number of grid points
Lx = [1;1];                      % Length of domain
dt = 1e-2;                   % Time step
Nt = 1e4;                    % # of time steps
BC = {'DN','NN'};                    % Boundary condition
tol = 1e-6;
k   = 1;                  % Impact of high mutation rate rho = r(1 - exp(-k(1-m)))
v   = 0.0;
RADtreat = 0.0;

ICval  = [];                 % Set to control heights of IC
ICtype = 'domain wall';       % Pick Initial condition type.

u = 100*ones(Nx(1)*Nx(2),1);

% Discritize domain
r = linspace(0,Lx(1),Nx(1))'; m = linspace(0,Lx(1),Nx(2))';
dr = r(2) - r(1); dm = m(2) - m(1);
dx = [dr;dm];

[R,M] = meshgrid(r,m);

% Mutation (diffusion rate)
prolifDiff = 0.001;
mutDiff = 0.001;

% Make Laplacian and Gradient
Laplacian = make_Laplacian(Nx,dx,BC,2,prolifDiff,mutDiff);
gradient = makeGradient(Nx,dx,BC,2);
I = speye(Nx(1)*Nx(2));

% Make proliferation term
rho = R.*(1-exp(-k*(1-M)));
% rho = 4*R.*M.*(1-M);
% rho = R.*(1-M);

% Kill rate
rcrit = 0.2;
fmin = 1.0;
fmax = 0.3;
a = fmin^2*rcrit^2/(4*fmax*(fmin + fmax));
delta = rho - fmin*rho.*(rcrit-rho)./(a+rho.^2).*(1-M);

% delta = 0.5*ones(size(delta));
% rho = ones(size(rho));

% % surf(R,M,rho-delta)
% surf(R,M,rho-delta);
% xlabel('Proliferation'); ylabel('Mutation'); colorbar;
% title('Treatment on');
% figure(2);
% surf(R,M,rho);
% xlabel('Proliferation'); ylabel('Mutation'); colorbar;
% title('Treatment off');

% We create a matrix solve using LHS*x_{n+1} = RHS*x_n, but we can save
% some time by initializing the part that stays fixed in time
CN = 0.5;
LHS_null= I - dt*CN*Laplacian;
RHS_null = I + dt*(1-CN)*Laplacian;

nskip = 50;

plotting = true;

BCGstart = [];
BCGoff = [];

% treatmentstart = [1;(100:500:10000)'];
% treatmentoff = [50;(400:500:10000)'];

BCGstart = [1;5000];
BCGoff = [2000;7000];

RADstart = [1000000;1000000];
RADoff = [200; 2100];

BCGon = false;
RADon = false;

moviecount = 1;
vidObj = VideoWriter('testmovie.avi');
open(vidObj);

% Movie(50) = struct('cdata',[],'colormap',[]);

%% Time loop using C-N scheme
for n = 1 : Nt
    
    % treatment time
    if sum(n == BCGstart)
        BCGon = true;
    elseif sum(n == BCGoff)
        BCGon = false;
    end
    if sum(n == RADstart)
        RADon = true;
    elseif sum(n == RADoff)
        RADon = false;
    end
    
    if BCGon
        LHS = spdiags(spdiags(LHS_null,0) - dt*CN*(rho(:) - delta(:)),0,LHS_null);
        RHS = spdiags(spdiags(RHS_null,0) + dt*(1-CN)*(rho(:) - delta(:)),0,RHS_null);
    else
        LHS = spdiags(spdiags(LHS_null,0) - dt*CN*rho(:),0,LHS_null);
        RHS = spdiags(spdiags(RHS_null,0) + dt*(1-CN)*rho(:),0,RHS_null);
    end
    if RADon
        LHS = LHS + dt*CN*(RADtreat+v).*gradient;
        RHS = RHS - dt*(1-CN)*(RADtreat+v).*gradient;
    else
        LHS = LHS + dt*CN*v*gradient;
        RHS = RHS - dt*(1-CN)*v*gradient;
    end
    
%     unew = LHS\(RHS*u)
%     
    [unew, ~, ~, ~, ~] = pcg(LHS,RHS*u,1e-10,1000);
    
%     if norm(unew - u) < tol
%         u = unew;
%         break;
%     end
    
    u = unew;
    TumorVolume(n,1) = sum(u);
    
    if plotting && mod(n,nskip)==0
        axis tight manual
        if BCGon && RADon
            txt = 'BCG and Radiation therapy';
        elseif BCGon
            txt = 'BCG on';
        elseif RADon
            txt = 'Radiation therapy applied';
        else
            txt = 'Treatment off';
        end
        drawnow;
        subplot(2,1,1)
        surf(R,M,reshape(u,Nx(1),Nx(2)));
        xticks([]);
        yticks([]);
        xlabel('Proliferation');
        ylabel('Mutation');
        shading interp
        view(0,90);
%         colorbar;
        title(txt);
        
        subplot(2,1,2)
        semilogy((1:n),TumorVolume)
        xlabel('time'); ylabel('Population size');
        ylim([1e5 1e11])
        xlim([0 10000])
        
        % Capture the plot as an image
       currFrame = getframe(gcf);
       writeVideo(vidObj,currFrame);
        
%         subplot(3,1,3);
%         if BCGon
%             surf(R,M,rho-delta);
%         else
%             surf(R,M,rho);
%         end
%         xlabel('Proliferation'); ylabel('Mutation'); colorbar;
%         view(0,90)
        

    end
    
    
end

% Close the file.
close(vidObj);

% % Output the movie as an avi file
% VideoWriter(Movie,'testMovie.avi');

end


function Laplacian = make_Laplacian(Nxvec,dxvec,BC,dim,mut,D)

if dim >= 1
    
    Nx = Nxvec(1);
    dx = dxvec(1);

    Ix = speye(Nx);
    ex = ones(Nx,1);
    Dxx = spdiags([ex -2*ex ex],[-1 0 1],Nx,Nx);
    
    if strcmp(BC{1},'NN')
        Dxx(1,2) = 2;
        Dxx(end,end-1) = 2;
    elseif strcmp(BC{1},'P')
        Dxx(1,end) = 1;
        Dxx(end,1) = 1;
    elseif strcmp(BC{1},'QP')
        Dxx(1,end) = -1;
        Dxx(end,1) = -1;
    elseif strcmp(BC{1},'ND')
        Dxx(1,2) = 2;
        Dxx(end,end-1) = 0;
    elseif strcmp(BC{1},'DN')
        Dxx(1,2) = 0;
        Dxx(end,end-1) = 2;
    end
    
    Dxx = mut*Dxx/dx^2;
    
    if dim == 1
        Laplacian = Dxx;
    end
    
end

if dim >= 2
    
    Ny = Nxvec(2);
    dy = dxvec(2);
    
    Iy = speye(Ny);
    ey = ones(Ny,1);
    Dyy = spdiags([ey -2*ey ey],[-1 0 1],Ny,Ny);
    
    if strcmp(BC{2},'NN')
        Dyy(1,2) = 2;
        Dyy(end,end-1) = 2;
    elseif strcmp(BC{2},'P')
        Dyy(1,end) = 1;
        Dyy(end,1) = 1;
    elseif strcmp(BC{2},'QP')
        Dyy(1,end) = -1;
        Dyy(end,1) = -1;
    elseif strcmp(BC{2},'ND')
        Dyy(1,2) = 2;
        Dyy(end,end-1) = 0;
    elseif strcmp(BC{2},'DN')
        Dyy(1,2) = 0;
        Dyy(end,end-1) = 2;
    end
    
    Dyy = D*Dyy/dy^2;
    
    if dim == 2
        Laplacian = kron(Iy,Dxx) + kron(Dyy,Ix);
    end
    
end

if dim == 3
    
    Nz = Nxvec(3);
    dz = dxvec(3);
    
    Iz = speye(Nz);
    ez = ones(Nz,1);
    Dzz = spdiags([ez -2*ez ez],[-1 0 1],Nz,Nz);
    dz = dx_v(3);
    
    if strcmp(BC{3},'NN') == 1
        Dzz(1,2) = 2;
        Dzz(end,end-1) = 2;
    end
    if strcmp(BC{3},'P') == 1
        Dzz(1,end) = 1;
        Dzz(end,1) = 1;
    end
    
    Dzz = Dzz/dz^2;
    
    Laplacian = kron(Iz,kron(Iy,Dxx)) + kron(Iz,kron(Dyy,Ix)) + ...
        kron(Dzz,kron(Iy,Ix));
    
end

end

function gradient = makeGradient(Nxvec,dxvec,BC,dim)

    Nx = Nxvec(1); Ny = Nxvec(2);
    dx = dxvec(1); dy = dxvec(2);

    Ix = speye(Nx);
    ex = ones(Nx,1);
    Dx = spdiags([-ex ex],[-1 1],Nx,Nx);

    Iy = speye(Ny);
    ey = ones(Ny,1);
    Dy = spdiags([-ey ey],[-1 1],Ny,Ny);

    if strcmp(BC{2},'NN')
        Dy(1,2) = 2;
        Dy(end,end-1) = 2;
    elseif strcmp(BC{2},'P')
        Dy(1,end) = 1;
        Dy(end,1) = 1;
    elseif strcmp(BC{2},'QP')
        Dy(1,end) = -1;
        Dy(end,1) = -1;
    elseif strcmp(BC{2},'ND')
        Dy(1,2) = 2;
        Dy(end,end-1) = 0;
    elseif strcmp(BC{2},'DN')
        Dy(1,2) = 0;
        Dy(end,end-1) = 2;
    end

    Dy = Dy/2/dy;

    if dim == 2
        gradient = kron(Iy,0*Dx) + kron(Dy,Ix);
    end

end