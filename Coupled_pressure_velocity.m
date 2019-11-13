clc
clear all
close all
gf_workspace('clear all')

%%
%We wrote the code being consistent with the notation used in Notaro,
%Cattaneo, Formaggia, Scotti, Zunino "A mixed finite element method for modeling 
%the fluid exchange between microcirculation and tissue interstitium" 2016

%% Omega e Lambda

PLOTTO=0; 

%% Omega 

%lx, ly and lz are the number of intervals on each direction of the mesh.
%There will be lx+1, ly+1 and lz+1 nodes/points on each direction
%It is a parallelepiped whose dimensions are defined by x0,y0 and z0,
%xf,yf,zf.

lx=5;   ly=5;   lz=5;
x0=0;   xf=1;   dx=xf-x0;
y0=0;   yf=1;   dy=yf-y0;
z0=0;   zf=1;   dz=zf-z0;

vectX= x0:(dx/lx):xf ;              
vectY= y0: (dy/ly) :yf ;            
vectZ= z0: dz/lz :zf ;
Omegaref=zeros(lx+1,ly+1,lz+1);
OmegaM=gf_mesh('regular simplices',vectX, vectY, vectZ,'degree', 1);    %Mesh over the cube DON'T TOUCH

nbsk=ceil(lx*2); %nbsk is the number of nodes of our vessel. 
                 %We defined it depending on the size of our cube, useful
                 %for quickly testing the code 

r=(min([dz/lz,dx/lx,dy/ly]))/20;    %Geometric radius of the vessel used to create the interpolation matrix
                                    %It is used to define the region
                                    %DELTALAMBDA (see later)
                                    
khat=16;                            %Number of elements of the discrete ring defining the vessel/pipe ( \gamma(s_k))

Lambda=zeros(nbsk,3);               
Lambda(:,1)= linspace(0+eps, 1-eps, nbsk ); %Coordinates of the vessel points embedded in the cube.                                            
Lambda(:,2)= 0.41;      
Lambda(:,3)= 0.41;          

Lambdaref= 0:1/(nbsk-1):1;
LambdaM=gf_mesh('cartesian', Lambdaref);    %Mesh over the vessel defined as the segment [0 1]. When considering the bifurcations, 
                                            %we feel that here is our first problem  


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

borderL = gf_mesh_get(LambdaM,'outer faces');      
gf_mesh_set(LambdaM, 'region',42, borderL);

BorderOut = gf_mesh_get(LambdaM, 'outer faces with direction', +1, 0);
BorderIn  = gf_mesh_get(LambdaM, 'outer faces with direction', -1, 0);
INFLOW_BOUNDARY=11;  OUTFLOW_BOUNDARY=22;
gf_mesh_set(LambdaM, 'region', OUTFLOW_BOUNDARY, BorderOut);
gf_mesh_set(LambdaM, 'region',  INFLOW_BOUNDARY, BorderIn);

%% Tubocavo 
tubocavo=zeros(nbsk,khat,3);            %Coordinates of \gamma(s_k) points for every s_k
normali=zeros(nbsk,3);
theta=0:2*pi/(khat-1):2*pi;             %The result will be a discretized pipe describing the vessel in Omega
                                        %The pipe will be surrounding our 1D domain
                                        
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


%% Building a finite element method on the mesh 
%Omega
mf3 = gf_mesh_fem(OmegaM, 3);       %3D Finite element method
mf1 = gf_mesh_fem(OmegaM, 1);       %1D                                    

gf_mesh_fem_set( mf3, 'fem', gf_fem('FEM_RT0(3)'));                 %Raviart-Thomas degree 1
gf_mesh_fem_set( mf1, 'fem', gf_fem('FEM_PK_DISCONTINUOUS(3,0)') ); %Pk disc. degree 0

DoF_Om3= gf_mesh_fem_get(mf3, 'nbdof');    %Degrees of Freedom over Omega 3D for velocity
DoF_Om1= gf_mesh_fem_get(mf1, 'nbdof');    %Same but for pressure

%Lambda
mfLv = gf_mesh_fem(LambdaM,1);       %1D FEM over Lambda for velocity
mfLp = gf_mesh_fem(LambdaM,1);       %1D FEM over Lambda for pressure
gf_mesh_fem_set(mfLv, 'fem', gf_fem('FEM_PK(1,1)'));                %Pk cont. degree 1 for velocity over Lambda                   
gf_mesh_fem_set(mfLp, 'fem', gf_fem('FEM_PK_DISCONTINUOUS(1,0)'));  %Pk disc. degree 0 for pressure over Lambda
DoF_Lv = gf_mesh_fem_get(mfLv, 'nbdof');
DoF_Lp = gf_mesh_fem_get(mfLp, 'nbdof');
DoF_TOT=DoF_Lp+DoF_Lv+DoF_Om1+DoF_Om3;

%% Intergration methods
% Omega
IT = gf_integ('IM_TETRAHEDRON(1)');   

mim3 = gf_mesh_im( mf3, IT);
mim1 = gf_mesh_im( mf1, IT);

%Lambda
IL = gf_integ('IM_GAUSS1D(1)');         

mimLv = gf_mesh_im(mfLv, IL);
mimLp = gf_mesh_im(mfLp, IL);

%% Interpolation Matrix 

WW=gf_mesh_im_get(mim1, 'im nodes', 1);         %Weights of quadrature formula 
W=WW(4,1);                                      
Pi_gamma = zeros(nbsk, khat, DoF_Om1);          %Local dicrete interpolation matrix
Bar_Pi = zeros(DoF_Lp,DoF_Om1);                 %Interpolation matrix

% Deltalambda
CVnumber=gf_mesh_get(OmegaM, 'max cvid');               %Number of convexes
deltalambda=[];                                         %Will contain the global id number of convexes belonging to the region of the vessel aka DELTALAMBDA
radius=4*max(gf_mesh_get(OmegaM,'convex radius'))+r;    %Radius used to filter if a node is close enough to a given convex, see later

%grandeciclo                                          
%long story short, find the convexes that contains at least one point of
%'tubocavo'. Absolutely not computationally efficient. 
% And creates the interpolation matrix
for i=1:CVnumber                                        %For every convex 
    PIDs = gf_mesh_get(OmegaM, 'pid in cvids', i);      %Obtain the coordinates of its vertices
    vertici= gf_mesh_get(OmegaM, 'pts', PIDs);
    for j=1:nbsk                                        %For every node (s_k) of the vessel 
        if pdist([Lambda(j,:);vertici(:,1)']) < radius  %Filter the nodes close enough to our convex
            for k=1:khat                                %For every point on the ring \gamma(s_k)
                p=squeeze(tubocavo(j,k,:))';            
                if ptintetra(vertici(:,1)',vertici(:,2)',vertici(:,3)',vertici(:,4)',p)     %Check if the k-th point of \gamma(s_k) is 
                    if not(ismember(i,deltalambda))                                         %inside the i-th convex (tetrahedron) CHECK FUNCTIONS BELOW
                        deltalambda=[deltalambda,i];                                        %If yes, we add the convex id to our vector deltalambda
                    end   
                    DOFs=gf_mesh_fem_get(mf1, 'basic dof from cvid', i);                    %Get the DoF of pressure in tissue of our i-th convex belonging to deltalambda   
                    if j==1
                        Bar_Pi(j,DOFs)=1;
                    elseif j==nbsk
                        Bar_Pi(j-1,DOFs)=1;
                    else
                        Bar_Pi([j-1,j],DOFs)=1;                                             %Link this DoF to the DoF of the pressure in the vessel
                    end
                    
                end
            end
        end
    end
end
deltalambda(2,:)=0;
DELTALAMBDA=17;                             %Create the region DELTALAMBDA from the vector of ids deltalambda
gf_mesh_set(OmegaM,'region', DELTALAMBDA, deltalambda);

Bar_PiT= sparse(Bar_Pi'*W);                 %Interpolation matrix and
Bar_Pi = sparse(Bar_Pi *W);                 %its transpose 

%% Parameters 
% Our cube is considered to be 2cmx2cmx2cm. So one space unit equals 2cm.
% We had to rescale all the values in Nabil et al. 2016 keeping that in mind. 

mu=4;                   % Blood viscosity
kappa=2.5e-13;          % Hydraulic permeability, interstitial volume
B=0.5/3600;             % Effective permeability, lymphatic vessels
R=7.64e-2;              % Capillary radius
Lp=2.5e-5;              % Hydraulic permeability, capillary wall

%Source terms parameters 
A=27.9;  %A=sigmap*(pi*v_v-pi*t_v); 
pL=1;           % Hydrostatic pressure within lymphatic channels
gvin=35;        % Pressure in the inflow boundary of the vessel   
gvout=0;        % Pressure in the outflow boundary of the vessel 
gt=.2;          % Pressure in (part of) the boundary of the cube 
gt_front=repmat(gt,3,3); 

Q=2*pi*R*Lp;        
BV = Q/(pi*R^2);    
musukappa=mu/kappa; 
C=8*mu /R^2;          

%% Model 

md=gf_model('real');  %Model built with the brick structure of GetFem++
gf_model_set(md, 'add fem variable', 'ut', mf3);
gf_model_set(md, 'add fem variable', 'pt', mf1);
gf_model_set(md, 'add fem variable', 'uv', mfLv);
gf_model_set(md, 'add fem variable', 'pv', mfLp);

gf_model_set(md, 'add initialized data', 'C', C);
gf_model_set(md, 'add initialized data', 'Q', Q);
gf_model_set(md, 'add initialized data', 'B', B);
gf_model_set(md, 'add initialized data', 'musukappa', musukappa); 
gf_model_set(md, 'add initialized data', 'BV' , BV ); 
gf_model_set(md, 'add initialized data', 'gvin',  gvin);
gf_model_set(md, 'add initialized data', 'gvout', gvout);
gf_model_set(md, 'add initialized data', 'pL', pL);
gf_model_set(md, 'add initialized data', 'A', A);
gf_model_set(md, 'add initialized data', 'R', R);
gf_model_set(md, 'add initialized data', 'Lp', Lp);
gf_model_set(md, 'add initialized data', 'gt', gt);
gf_model_set(md, 'add initialized data', 'gt_front', gt_front);

% Tangent Matrix
gf_model_set(md, 'add mass brick', mim3, 'ut', 'musukappa',-1); %Mtt
gf_model_set(md, 'add mass brick', mim1, 'pt', 'B', -1);        %Btt solo una parte
gf_model_set(md, 'add nonlinear generic assembly brick', mim3, '-Div_Test_ut.pt', -1, 0, 0); %-DttT
gf_model_set(md, 'add nonlinear generic assembly brick', mim1, 'Div_ut.Test_pt' , -1, 0, 0); % Dtt
gf_model_set(md, 'add mass brick', mimLv, 'uv', 'C', -1); %Mvv
gf_model_set(md, 'add nonlinear generic assembly brick', mimLp,  'Grad_uv.Test_pv', -1, 0, 0); %Dvv
gf_model_set(md, 'add nonlinear generic assembly brick', mimLv, '-pv.Grad_Test_uv', -1, 0, 0); %-DvvT

% Coupling bricks, construction though sparse matrices (assembly command of GetFem++)
MvvP = gf_asm('mass matrix', mimLp, mfLp);

Btt = gf_spmat('mult',Bar_PiT, MvvP); Btt = gf_spmat('mult', Btt, Bar_Pi);  gf_spmat_set(Btt, 'scale',  Q); % Btt su Lambda
Btv = gf_spmat('mult',Bar_PiT, MvvP);                                       gf_spmat_set(Btv, 'scale', -Q);
Bvt = gf_spmat('mult', MvvP, Bar_Pi);                                       gf_spmat_set(Bvt, 'scale', -BV);
Bvv = gf_spmat('copy', MvvP);                                               gf_spmat_set(Bvv, 'scale',  BV);
% Coupling bricks, adding to the model using again the brick structure of
% GetFem++
gf_model_set(md, 'add explicit matrix', 'pt', 'pt', Btt, 1, 0); %Btt
gf_model_set(md, 'add explicit matrix', 'pt', 'pv', Btv, 0, 0); %Btv
gf_model_set(md, 'add explicit matrix', 'pv', 'pt', Bvt, 0, 0); %Bvt
gf_model_set(md, 'add explicit matrix', 'pv', 'pv', Bvv, 1, 0); %Bvv

% Buondaries, adding to te model
gf_model_set(md, 'add source term brick', mimLv, 'uv', 'gvin',   INFLOW_BOUNDARY);      %Fv
gf_model_set(md, 'add source term brick', mimLv, 'uv', 'gvout', OUTFLOW_BOUNDARY);      %Fv
gf_model_set(md, 'add source term brick', mim1,  'pt', 'B*pL',   -1);                   %S1 su tutto Omega
gf_model_set(md, 'add source term brick', mim1,  'pt', 'A*Q', DELTALAMBDA);             %S1 su DELTALAMBDA
gf_model_set(md, 'add source term brick', mimLp, 'pv', 'BV*A',   -1);                   %S2
gf_model_set(md, 'add normal source term brick', mim3, 'ut', 'gt_front', FRONT_BOUND);  %Ft

%% Solve

gf_model_get(md, 'assembly', 'build_all');
AA = gf_model_get(md, 'tangent_matrix');
rhs= gf_model_get(md, 'rhs');

%% Preconditioning of the matrix
[PP,RR,CC] = equilibrate(AA);
BB=RR*PP*AA*CC;
d=RR*PP*rhs';
Y = gf_linsolve('lu', BB, d );

%% Final part

X = CC*Y;

Pt=X(1:DoF_Om1)';                               %Pressure field on the tissue
pv=X(DoF_Om1+1:DoF_Om1+DoF_Lp)';                %Pressure vector on the vessel
Ut=X(DoF_Om1+DoF_Lp+1:DoF_Om1+DoF_Lp+DoF_Om3)'; %Velocity field on the tissue
uv=X(end-DoF_Lv+1:end)';                        %Velocity vector on the vessel

if (PLOTTO)
    figure(1)
    gf_plot(mfLp, pv);
    title('Pressione nel vaso');
    figure(2)
    gf_plot(mfLv, uv);
    title('Velocità nel vaso');
    figure(3)
    gf_plot(mf1, Pt, 'mesh', 'on', 'cvlst', gf_mesh_get(OmegaM, 'outer faces') ); hold on;
    gf_plot(mf1, Pt, 'mesh', 'on', 'cvlst', gf_mesh_get(OmegaM, 'inner faces') ); hold off;
    title('Mappa pressione nel tessuto');
    colorbar;
    figure(4)
    gf_plot(mf3, Ut, 'mesh', 'off', 'cvlst', gf_mesh_get(OmegaM, 'outer faces') ); hold on;
    gf_plot(mf3, Ut, 'mesh', 'off', 'cvlst', gf_mesh_get(OmegaM, 'inner faces') ); hold off;
    title('Mappa velocità nel tessuto');

end

%% Functions
%Functions used to check if a given point p is inside the tetrahedron of
%vertices v1 v2 v3 v4.
function samesid = sameside(v1,v2,v3,v4,p)
    normal = cross(v2-v1,v3-v1);
    dotV4 = dot(normal, v4-v1);
    dotP = dot(normal,p-v1);
    samesid = sign(dotV4) == sign(dotP);
end

function ptintetr = ptintetra(v1,v2,v3,v4,p)
    ptintetr = sameside(v1, v2, v3, v4, p) && sameside(v2, v3, v4, v1, p) && sameside(v3, v4, v1, v2, p) && sameside(v4, v1, v2, v3, p); 
end




