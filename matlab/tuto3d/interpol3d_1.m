clc;
clear all;
close all;

addpath('../bin');

%% load domain (nodes + boundary(s)) and creat dummy  then plot

load('data_interpol_1');

%% creat dummy point location for interpolation evaluation

XYZ_Point=rand(1000,3)*1.5-0.75;

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

%% initialize the interpolator and interpolate a filds
% test : Var = XYZ_Noeud ==> interpolated Var on XYZ_Point = XYZ_Point

Fxyz=naturalInterpolant(XYZ_Noeud,XYZ_Noeud,'Sibson',IN_Tri_Ini,'tetgen');% tetgen : convex + non convex delaunay mesher
%Fxyz=naturalInterpolant(XYZ_Noeud,XYZ_Noeud,'Sibson');% matlab : only convex delaunay mesher

Var_Int=Fxyz.eval(XYZ_Point);

%% cal error

dif=Var_Int-XYZ_Point;
j=0;
ind_p_in=zeros(sum(Fxyz.In_Out),1);
for i=1:size(Fxyz.In_Out,1)
    if Fxyz.In_Out(i)
        j=j+1;
        ind_p_in(j)=i;
    end
end

int_p_out=setdiff(1:size(Fxyz.In_Out,1),ind_p_in);

err=max(max(abs(dif(ind_p_in,:))))

%% plot in out point

figure;
hold on;
tri_out_handle=trimesh(IN_Tri_Ini,XYZ_Noeud(:,1),XYZ_Noeud(:,2),XYZ_Noeud(:,3),'edgecolor','black');
alpha(tri_out_handle,0.5);
plot3(XYZ_Point(ind_p_in,1),XYZ_Point(ind_p_in,2),XYZ_Point(ind_p_in,3),'*','color','red');
plot3(XYZ_Point(int_p_out,1),XYZ_Point(int_p_out,2),XYZ_Point(int_p_out,3),'*','color','blue');
axis vis3d
axis equal
hold off;
