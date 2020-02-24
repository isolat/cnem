clc;
clear all;
close all;

addpath('../bin');
%addpath('Z:/DATA/1/dev/old/cnem/cpp/win/test/Debug')
addpath('Z:/DATA/1/dev/colab/nico/thermo3d')
%%
with_matlab=false;

%% load domain (nodes + boundary(s)) and creat dummy  then plot

load('data_interpol_1');
tic

%% initialize the interpolator

[X,Y,Z] = meshgrid(-0.5:0.01:0.5,-0.5:0.01:0.5,-0.5:0.01:0.5);
XYZ_Noeud=[reshape(X,[],1),reshape(Y,[],1),reshape(Z,[],1)];
Var=XYZ_Noeud;

if ~with_matlab

F_xyz=naturalInterpolant(XYZ_Noeud,Var,'Sibson');
%F_xyz=naturalInterpolant(XYZ_Noeud,Var,'Sibson',IN_Tri_Ini);% with
%Constraints (domaine boundary,can be non convex)
else
F_x = scatteredInterpolant(XYZ_Noeud(:,1),XYZ_Noeud(:,2),XYZ_Noeud(:,3),Var(:,1),'natural');
F_y = scatteredInterpolant(XYZ_Noeud(:,1),XYZ_Noeud(:,2),XYZ_Noeud(:,3),Var(:,2),'natural');
F_z = scatteredInterpolant(XYZ_Noeud(:,1),XYZ_Noeud(:,2),XYZ_Noeud(:,3),Var(:,3),'natural');
end

for I=1:1

%% creat dummy point location for interpolation evaluation

XYZ_Point=rand(1000,3)*1.5-0.75;

%% plot

% figure;
% hold on;
% tri_out_handle=trimesh(IN_Tri_Ini,XYZ_Noeud(:,1),XYZ_Noeud(:,2),XYZ_Noeud(:,3),'edgecolor','black');
% alpha(tri_out_handle,0.5);
% plot3(XYZ_Noeud(:,1),XYZ_Noeud(:,2),XYZ_Noeud(:,3),'.','color','green');
% plot3(XYZ_Point(:,1),XYZ_Point(:,2),XYZ_Point(:,3),'o','color','blue');
% axis vis3d
% axis equal
% hold off;

if ~with_matlab

%% interpolate a filds
% test : Var = XYZ_Noeud ==> interpolated Var on XYZ_Point = XYZ_Point

Var_Int=F_xyz.eval(XYZ_Point);
%Var_Int=F_xyz.mat_interpol_glob*Var;

%% cal error
 
dif=Var_Int-XYZ_Point;
%err=max(max(abs(dif)))

j=0;
ind_p_in=zeros(sum(F_xyz.In_Out),1);
for i=1:size(F_xyz.In_Out,1)
    if F_xyz.In_Out(i)
        j=j+1;
        ind_p_in(j)=i;
    end
end

int_p_out=setdiff(1:size(F_xyz.In_Out,1),ind_p_in);

err=max(max(abs(dif(ind_p_in,:))))

else

%% interpolate a filds
% test : Var = XYZ_Noeud ==> interpolated Var on XYZ_Point = XYZ_Point

Var_Int=[F_x(XYZ_Point),F_y(XYZ_Point),F_z(XYZ_Point)];

%% cal error
 
dif=Var_Int-XYZ_Point;

err=max(max(abs(dif)))

end

%% plot in out point

% figure;
% hold on;
% tri_out_handle=trimesh(IN_Tri_Ini,XYZ_Noeud(:,1),XYZ_Noeud(:,2),XYZ_Noeud(:,3),'edgecolor','black');
% alpha(tri_out_handle,0.5);
% plot3(XYZ_Point(ind_p_in,1),XYZ_Point(ind_p_in,2),XYZ_Point(ind_p_in,3),'*','color','red');
% plot3(XYZ_Point(int_p_out,1),XYZ_Point(int_p_out,2),XYZ_Point(int_p_out,3),'*','color','blue');
% axis vis3d
% axis equal
% hold off;


%%
end

toc

%%
if ~with_matlab
%delete(F_xyz);
%clear('F_xyz');
end