 clc
 clear all
 close all
 gf_workspace('clear all')

%% Omega e Lambda

PLOTTO=0; 

%% Omega 

lx=5;   ly=5;   lz=5;
x0=0;   xf=1;   dx=xf-x0;
y0=0;   yf=1;   dy=yf-y0;
z0=0;   zf=1;   dz=zf-z0;

vectX= x0:(dx/lx):xf ;             
vectY= y0: (dy/ly) :yf ;            
vectZ= z0: dz/lz :zf ;
Omegaref=zeros(lx+1,ly+1,lz+1);
OmegaM=gf_mesh('regular simplices',vectX, vectY, vectZ,'degree', 3);    %Mesh over the cube DON'T TOUCH

nbsk=5;
nbranches=3;
r=(min([dz/lz,dx/lx,dy/ly]))/20;   

khat=16;

Lambda=zeros(nbranches,nbsk,3);                 %For simplicity, we have 3 branches, each one with the same number of nodes         
Lambda(1,:,1)= linspace(0+eps,0.5,nbsk);        %Lambda=Lambda(#branch, #nodes, x-y-z coordinates)                                    
Lambda(1,:,2)= repmat(0.5,1,nbsk);       
Lambda(1,:,3)= repmat(0.5,1,nbsk);   


Lambda(2,:,1)= linspace(0.5,1,nbsk);                                            
Lambda(2,:,2)= linspace(0.5,1,nbsk);         
Lambda(2,:,3)= repmat(0.5,1,nbsk);


Lambda(3,:,1)= linspace(0.5,1,nbsk);                                            
Lambda(3,:,2)= linspace(0.5,0,nbsk);        
Lambda(3,:,3)= repmat(0.5,1,nbsk);



LambdaM=[];
for i=1:nbranches
    LambdaM =[LambdaM,gf_mesh('regular simplices',linspace(0,norm(squeeze(Lambda(i,end,:)-Lambda(i,1,:))),nbsk),[0 0 0],[0 0 0] ,'degree',1)]; % canot define as 3D (*)
    %LambdaM=[LambdaM,gf_mesh('regular simplices',linspace(0,norm(squeeze(Lambda(i,end,:)-Lambda(i,1,:))),nbsk), 'degree',1)];                   nor as 1D, see later
   
    P = gf_mesh_get(LambdaM(i), 'pts');
    PIDtotali=gf_mesh_get(LambdaM(i), 'pid');
    PIDs = gf_mesh_get(LambdaM(i), 'pid from coords', P);
    PIDdatenere=unique(PIDs);
    delete=[];
    for j=1:length(PIDtotali)
        if not(ismember(PIDtotali(j),PIDdatenere))
            delete=[delete,PIDtotali(j)];
        end
    end
%     CVIDdelete=gf_mesh_get(LambdaM(i), 'cvid from pid', delete, 1);
%     gf_mesh_set(LambdaM(i), 'del convex', CVIDdelete);
%     %gf_mesh_set(LambdaM(i), 'del point', delete);
%     

end
%Our idea is to define 3 simple separate branches over wich we define 3 independent 1D
%meshes in the space. Then, we want to translate and rotate our branches in
%the 3D space according to the coordinates of the matrix Lambda. We checked
%that using GetFem++ we can translate and rotate a mesh and we can merge
%different meshes.
% (*) We had problems with the very definition of a 1D mesh embedded in a 3D
%domain. In fact we get that our convexes are degenerate tetrahedrons if we
%try to define a 1D mesh in the space.
%If we try to define the vessel in a 1D space, we cannot move it in a 3D
%space.

%Our idea was actually to position and then merge the previously single
%defined independent branches

%Compute R=rotation matrix given the two branches
p0=squeeze(Lambda(1,2,:)-Lambda(1,1,:));
p1=squeeze(Lambda(2,2,:)-Lambda(2,1,:));

C = cross(p0, p1); 
D = dot(p0, p1);
NP0 = norm(p0,3) ;  % used for scaling
if ~all(C==0)       % check for colinearity    
    Z = [0 -C(3) C(2); C(3) 0 -C(1); -C(2) C(1) 0] ; 
    R = (eye(3) + Z + Z^2 * (1-D)/(norm(C)^2)) / NP0^2 ; % rotation matrix
else
    R = sign(D) * (norm(p1) / NP0) ; % orientation and scaling
end

gf_mesh_set(LambdaM(2), 'transform', R)
gf_mesh_set(LambdaM(2), 'translate', squeeze(Lambda(2,1,:)));

gf_mesh_set(LambdaM(1), 'merge', LambdaM(2));
gf_mesh_set(LambdaM(1), 'merge', LambdaM(1));
figure()
gf_plot_mesh(LambdaM(1),'vertices','on'); 
gf_mesh_get(LambdaM(2),'nbpts')




