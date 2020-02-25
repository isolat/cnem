clc;
clear all;
close all;

addpath('../bin');

%% load domain (nodes + boundary(s)) and creat dummy  then plot

load('data_interpol_1');

%% creat dummy point location for interpolation evaluation

XYZ_Point=rand(100,3)*1.5-0.75;

%% plot

figure;
hold on;
tri_out_handle=trimesh(IN_Tri_Ini,XYZ_Noeud(:,1),XYZ_Noeud(:,2),XYZ_Noeud(:,3),'edgecolor','black');
alpha(tri_out_handle,0.5);
plot3(XYZ_Noeud(:,1),XYZ_Noeud(:,2),XYZ_Noeud(:,3),'.','color','green');
plot3(XYZ_Point(:,1),XYZ_Point(:,2),XYZ_Point(:,3),'o','color','blue');
axis vis3d
axis equal
hold off;

%% compute raw Sibson shape function and gradiant

Var = XYZ_Noeud; % the var field are the nodes coordonate

Fxyz=naturalInterpolant(XYZ_Noeud,Var,'SibsonRaw'); % with convex hull as constrain
%Fxyz=naturalInterpolant(XYZ_Noeud,Var,'SibsonRaw',IN_Tri_Ini,'tetgen'); %
%with boundary faces as constrain

%% check that gradiant at XYZ_Point is [1,0,0;0,1,0;0,0,1]

Var_Int=Fxyz.eval(XYZ_Point);

grad_x_node=Fxyz.Mat_GradX*XYZ_Noeud;
grad_y_node=Fxyz.Mat_GradY*XYZ_Noeud;
grad_z_node=Fxyz.Mat_GradZ*XYZ_Noeud;

nb_Pnt=size(XYZ_Point,1);

j=0;
for i=1:nb_Pnt
    if Fxyz.Nb_V(i)>3 
        % ie :
        %   Nb_V == 0, the point is out of the domain
        %   Nb_V == 1, the point is coincident with a node
        %   Nb_V == 3, the point is on boundary triangle or the convex hull
        %   if no bondary face are supplied
        %   Nb_V > 3, the point is in the domain and not coincident with a
        %             node
        
        j=j+1;
    end
end

if j~=0 
    check_grad=zeros(j,9);
    j=1;
    for i=1:nb_Pnt
        if Fxyz.Nb_V(i)>3
            check_grad(j,:)=[grad_x_node(i,:),grad_y_node(i,:),grad_z_node(i,:)]-[1,0,0,0,1,0,0,0,1];
            j=j+1;
        end
    end 

    err=max(max(abs(check_grad)))
end

%%






