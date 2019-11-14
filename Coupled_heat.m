clc
clear all
close all
gf_workspace('clear all')

PLOTTO=1;%vuoi che plotto?
plot_slices=0;

% Code equal to "Coupled_pressure_velocity". 
% See section %% Build FEM over the mesh
%% Omega 
%creazione ambiente
%dimensioni Omega

lx=8;   ly=8;   lz=8;
x0=0;   xf=1;   dx=xf-x0;
y0=0;   yf=1;   dy=yf-y0;
z0=0;   zf=1;   dz=zf-z0;

vectX= x0:(dx/lx):xf ;
vectY= y0: dy/ly :yf ;
vectZ= z0: dz/lz :zf ;

Omegaref=zeros(lx+1,ly+1,lz+1);
OmegaM=gf_mesh('regular simplices',vectX, vectY, vectZ,'degree', 1); %degree x = grado x dei polinomi sui segmenti, e quindi n+1 elementi su ciascuno
%Dimensioni Lambda

nbsk=ceil(lx*2);
r=(min([dz/lz,dx/lx,dy/ly]))/20;
khat=8;
Lambda=zeros(nbsk,3);
Lambda(:,1)= linspace(0+eps,1-eps, nbsk );%    
Lambda(:,2)= linspace(1-eps,0+eps, nbsk );%   
Lambda(:,3)= linspace(1-eps,0+eps, nbsk );%    

Lambdaref= 0:1/(nbsk-1):1;
LambdaM=gf_mesh('cartesian', Lambdaref);


%% Boundary Regions

fb1 = gf_mesh_get(OmegaM, 'outer faces with direction', [ 1  0  0], 0.01);
fb2 = gf_mesh_get(OmegaM, 'outer faces with direction', [-1  0  0], 0.01);
fb3 = gf_mesh_get(OmegaM, 'outer faces with direction', [ 0  1  0], 0.01);
fb4 = gf_mesh_get(OmegaM, 'outer faces with direction', [ 0 -1  0], 0.01);
fb5 = gf_mesh_get(OmegaM, 'outer faces with direction', [ 0  0  1], 0.01);
fb6 = gf_mesh_get(OmegaM, 'outer faces with direction', [ 0  0 -1], 0.01);
fb7 = gf_mesh_get(OmegaM, 'outer faces');

RIGHT_BOUND=1; LEFT_BOUND=2; BACK_BOUND=3; FRONT_BOUND=4; TOP_BOUND=5; BOTTOM_BOUND=6; VESSEL_BOUND=7; DE_OMEGA=8;

gf_mesh_set(OmegaM, 'region',  RIGHT_BOUND, fb1);
gf_mesh_set(OmegaM, 'region',   LEFT_BOUND, fb2);
gf_mesh_set(OmegaM, 'region',   BACK_BOUND, fb3);
gf_mesh_set(OmegaM, 'region',  FRONT_BOUND, fb4);
gf_mesh_set(OmegaM, 'region',    TOP_BOUND, fb5);
gf_mesh_set(OmegaM, 'region', BOTTOM_BOUND, fb6);
gf_mesh_set(OmegaM, 'region',     DE_OMEGA, fb7);
%DoF_Front=length(gf_mesh_fem_get(mf3, 'basic dof on region',FRONT_BOUND));

%gf_mesh_set(OmegaM,'region merge', FRONT_BOUND,  RIGHT_BOUND);
%gf_mesh_set(OmegaM,'region merge',  BACK_BOUND,   LEFT_BOUND);
%gf_mesh_set(OmegaM,'region merge',   TOP_BOUND, BOTTOM_BOUND);

% detect the border of the mesh
borderL = gf_mesh_get(LambdaM,'outer faces');       %bordi di Lmabda, entrambi. per ora è inutile
gf_mesh_set(LambdaM, 'region',42, borderL);

BorderOut = gf_mesh_get(LambdaM, 'outer faces with direction', +1, 0);
BorderIn  = gf_mesh_get(LambdaM, 'outer faces with direction', -1, 0);
INFLOW_BOUNDARY=11;  OUTFLOW_BOUNDARY=22;
gf_mesh_set(LambdaM, 'region', OUTFLOW_BOUNDARY, BorderOut);
gf_mesh_set(LambdaM, 'region',  INFLOW_BOUNDARY, BorderIn);

%% Tubocavo
tubocavo=zeros(nbsk,khat,3); %Per ogni foglio hai le coordinate dei punti dell'anello
normali=zeros(nbsk,3);
theta=0:2*pi/(khat-1):2*pi;

for i=1:nbsk
    if i~=nbsk
        normali(i,:)=Lambda(i+1,:)-Lambda(i,:);
    else
        normali(i,:)=Lambda(i-1,:)-Lambda(i,:);
    end
    v=null(normali(i,:));
    tubocavo(i,:,:)=(repmat(Lambda(i,:)',1,size(theta,2))+r*(v(:,1)*cos(theta)+v(:,2)*sin(theta)))';
end
if(0)
    plot3(Lambda(:,1),Lambda(:,2),Lambda(:,3),'*');
    hold on
    plot3(tubocavo(:,:,1),tubocavo(:,:,2),tubocavo(:,:,3),'o')
    xlim([0 1]); ylim([0 1]); zlim([0 1]);
end

%% Build FEM over the mesh

mfT = gf_mesh_fem(OmegaM, 1);
gf_mesh_fem_set(mfT, 'fem', gf_fem('FEM_PK(3,1)') ); %Continuous Pk of degree 1
DoF_OmT= gf_mesh_fem_get(mfT, 'nbdof');

%% Intergration methods

IT = gf_integ('IM_TETRAHEDRON(1)');   %Integrazione Tessuto, che va bene sia per 1 che per 3        
mimT = gf_mesh_im( mfT, IT);

%% Deltalambda without interpolation matrix 
%See "Coupled_pressure_velocity"

% Deltalambda
CVnumber=gf_mesh_get(OmegaM, 'max cvid'); %numero di convessi nella mesh
deltalambda=[];
radius=4*max(gf_mesh_get(OmegaM,'convex radius'))+r;

%grandeciclo
for i=1:CVnumber
    PIDs = gf_mesh_get(OmegaM, 'pid in cvids', i);
    vertici= gf_mesh_get(OmegaM, 'pts', PIDs);
    for j=1:nbsk
        if pdist([Lambda(j,:);vertici(:,1)']) < radius
            for k=1:khat
                p=squeeze(tubocavo(j,k,:))';
                if ptintetra(vertici(:,1)',vertici(:,2)',vertici(:,3)',vertici(:,4)',p)
                    if not(ismember(i,deltalambda))
                        deltalambda=[deltalambda,i];
                    end                                        
                end
            end
        end
    end
end
deltalambda(2,:)=0;

%% Boundary Regions

fb1 = gf_mesh_get(OmegaM, 'outer faces with direction', [ 1  0  0], 0.01);
fb2 = gf_mesh_get(OmegaM, 'outer faces with direction', [-1  0  0], 0.01);
fb3 = gf_mesh_get(OmegaM, 'outer faces with direction', [ 0  1  0], 0.01);
fb4 = gf_mesh_get(OmegaM, 'outer faces with direction', [ 0 -1  0], 0.01);
fb5 = gf_mesh_get(OmegaM, 'outer faces with direction', [ 0  0  1], 0.01);
fb6 = gf_mesh_get(OmegaM, 'outer faces with direction', [ 0  0 -1], 0.01);

RIGHT_BOUND=1; LEFT_BOUND=2; BACK_BOUND=3; FRONT_BOUND=4; TOP_BOUND=5; BOTTOM_BOUND=6; VESSEL_BOUND=7;

gf_mesh_set(OmegaM, 'region',  RIGHT_BOUND, fb1);
gf_mesh_set(OmegaM, 'region',   LEFT_BOUND, fb2);
gf_mesh_set(OmegaM, 'region',   BACK_BOUND, fb3);
gf_mesh_set(OmegaM, 'region',  FRONT_BOUND, fb4);
gf_mesh_set(OmegaM, 'region',    TOP_BOUND, fb5);
gf_mesh_set(OmegaM, 'region', BOTTOM_BOUND, fb6);

DELTALAMBDA=17;
gf_mesh_set(OmegaM,'region', DELTALAMBDA, deltalambda);

%% Parameters

kappa = 1;
rho   = 1;
gamma = 1; 
BetaT = 15; 
R     = 10;
cV    = 1;
Lpsv  = 0.5; 
pL    = 0.5;
SAR   = 1.3; 
Tbl   = 37;
ct    = 1; 

pt = ones(1, DoF_OmT)*1.1;
D  = Lpsv*(pt-pL);          % D=L_p^{SF}*s/v*(p_t - p_L)
N  = BetaT/(rho*gamma);     %N= \beta_t / (rho*gamma)
A  = pi*R^2*cV;             % A=pi*R^2*c_v
ut = ones(DoF_OmT, 3)*0.1;
P  = 2*pi*R*BetaT;          %P=2*pi*R*\beta_t
vecFd = rho*gamma*Tbl*D;
 
%% Model

md=gf_model('real'); %modello  vuoto, reale (non complesso)

gf_model_set(md, 'add fem variable', 'T', mfT);
% Initialize data
gf_model_set(md, 'add initialized data', 'k', kappa); 
gf_model_set(md, 'add initialized data', 'rho', rho);
gf_model_set(md, 'add initialized data', 'gamma', gamma);
gf_model_set(md, 'add initialized data', 'N', N);
gf_model_set(md, 'add initialized data', 'A', A);
gf_model_set(md, 'add initialized data', 'SAR', SAR);
gf_model_set(md, 'add initialized data', 'ct', ct);
gf_model_set(md, 'add initialized data', 'Tbl', Tbl);
gf_model_set(md, 'add initialized data', 'P', P);

AA=gf_asm('mass matrix', mimT, mfT); AA=rho*gamma*(AA.*D);  
gf_model_set(md, 'add explicit matrix', 'T', 'T', AA, 1, 0);                    % rho*gamma* D(T,\xi)_Omega  -> CONTIENE Pt !!!!!

gf_model_set(md, 'add mass brick', mimT, 'T', 'rho*gamma*N', DE_OMEGA );        % rho*gamma* N(T,\xi)_De_Omega
gf_model_set(md, 'add mass brick', mimT, 'T', 'P', DELTALAMBDA );               %(\delta_\Lambda)*P*(T,\xi)_Omega

gf_model_set(md, 'add nonlinear generic assembly brick', mimT, 'k*Grad_T.Grad_Test_T', -1, 0, 0); %k(\Nabla T, \Nabla \xi)_Omega

gf_model_set(md, 'add initialized fem data', 'ut', mfT, ut); 
gf_model_set(md, 'add nonlinear generic assembly brick', mimT, '-rho*gamma*(ut*T).Grad_Test_T', -1);  %-rho*gamma* (T*u,\Nabla \xi)_Omega

%Questo è il theta method, theta = 1, first-order implicit/backward Euler Scheme, assolut/. incondizionat/. st
gf_model_set(md, 'add theta method for first order', 'T', 1);               %stiamo dicendo "esiste Dot_T"
gf_model_set(md, 'add mass brick', mimT, 'Dot_T', 'rho*gamma', -1 );        %rho*gamma*([\de/\de(t)](T),\xi)_Omega

%Source Terms
gf_model_set(md, 'add source term brick', mimT, 'T', 'P*Tbl', DELTALAMBDA );        %\delta_\Lambda*P(T_bl,\xi)_Omega
gf_model_set(md, 'add source term brick', mimT, 'T', 'A*SAR', DELTALAMBDA );        %\delta_\Lambda*(A,\xi)_Omega
gf_model_set(md, 'add source term brick', mimT, 'T', 'SAR*ct', -1 );                %SAR(c_t,\xi)_Omega
gf_model_set(md, 'add source term brick', mimT, 'T', 'rho*gamma*Tbl', DE_OMEGA );   %rho*gamma(T_bl,\xi)_DEOmega
V= gf_asm('volumic source', mimT, mfT, mfT, vecFd ) ;   
gf_model_set(md, 'add explicit rhs', 'T', V);                                       %rho*gamma*D(T_bl,\xi)_Omega  -> CONTIENE Pt !!!
%il giochino del tempo

%% Solve
% Interpolate the initial data
T0=(30)*ones(1,DoF_OmT);                %Initial Temperature
DotT0=T0*0;

Tempo=30; NumIntervalli=10;
dt=Tempo/NumIntervalli; ddt=dt/20.;

gf_model_set(md, 'set time', 0);
gf_model_set(md, 'set time step', dt);

% Initial data.
gf_model_set(md, 'variable', 'Previous_T',  T0);
gf_model_set(md, 'variable', 'Previous_Dot_T',  DotT0);

gf_model_set(md, 'perform init time derivative', dt/20.);
gf_model_get(md, 'solve'); 

gf_model_get(md, 'assembly', 'build_all');
%A =gf_model_get(md, 'assembly', 'build_matrix');
AA = gf_model_get(md, 'tangent_matrix');
rhs= gf_model_get(md, 'rhs');
TT=gf_linsolve('superlu',AA,rhs);
minn=zeros(1,NumIntervalli+1); ii=1;
maxx=minn(:);
minn(ii)=min(TT);
maxx(ii)=max(TT); ii=ii+1;
s2=gf_slice({'planar', 0, [0.5;0.5;0.51], [0;0;1]}, OmegaM, 1);
nbslices=4; sliceh=linspace(0,1,nbslices);
tic,
for t=0:dt:Tempo                                               
    
    gf_model_get(md, 'solve');
    
    DotT=gf_model_get(md, 'variable', 'Dot_T');
    T = gf_model_get(md, 'variable', 'T');
    minn(ii)=min(T);
    maxx(ii)=max(T); ii=ii+1;
    if(0)
        figure(1), colorbar, caxis([0,40]);
        
        for i=1:nbslices-1            
            sl = gf_slice({'planar', 0, [0.5;0.5;sliceh(i)], [0;0;1]}, OmegaM, 1);
            PlotT=gf_compute(mfT,T,'interpolate on', sl);
             subplot(1,2,1),colorbar, caxis([0,40]);
            gf_plot_slice(sl,'mesh_faces','off','mesh','off','data', PlotT,'mesh_slice_edges','off'); hold on,
        end
        subplot(1,2,2),colorbar, caxis([0,40]);
        %gf_plot(mfT, T, 'mesh', 'off', 'cvlst', gf_mesh_get(OmegaM, 'outer faces')); hold on,
        
        
        PlotT2=gf_compute(mfT, T, 'interpolate on', s2);
        gf_plot_slice(s2,'mesh_faces','off','mesh','off','data', PlotT2,'mesh_slice_edges','off');  hold on,
        
        pause(0.1);
        title(['Tempo: ' num2str(t) '/' num2str(Tempo)])
        %gf_plot(mfT, PlotT, 'mesh', 'off','cvlst', gf__mesh_get(OmegaM, 'outer faces')); hold on;
        %xlim=[0,1]; ylim=[0,1]; zlim=[0.1];
    end
    gf_model_set(md, 'shift variables for time integration');
end
toc,

if (1)
    figure(3), colorbar;
    gf_plot(mfT, T, 'mesh', 'off', 'cvlst', gf_mesh_get(OmegaM, 'inner faces')); hold on,
    title('Temperatura')
    
end

if (1)
    figure(2),
    plot(linspace(0,Tempo,NumIntervalli+2), minn); hold on,
    plot(linspace(0,Tempo,NumIntervalli+2), maxx); hold off
    %ylim([0,40]), 
    ylabel('Temperature'), xlabel('Time');
    legend('Temp Min','Temp Max');
    title('Max Temp & Min Temp');
end


%% Functions

function samesid = sameside(v1,v2,v3,v4,p)
    normal = cross(v2-v1,v3-v1);
    dotV4 = dot(normal, v4-v1);
    dotP = dot(normal,p-v1);
    samesid = sign(dotV4) == sign(dotP);
end

function ptintetr = ptintetra(v1,v2,v3,v4,p)
    ptintetr = sameside(v1, v2, v3, v4, p) && sameside(v2, v3, v4, v1, p) && sameside(v3, v4, v1, v2, p) && sameside(v4, v1, v2, v3, p); 
end




