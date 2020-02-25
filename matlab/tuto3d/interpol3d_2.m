clc;
clear all;
close all;

addpath('../bin');

%%
with_matlab=true;% comparison between matlab natural interpolation and the supplied library
%%
tic
%% initialize the interpolator

[X,Y,Z] = meshgrid(-0.5:0.1:0.5,-0.5:0.1:0.5,-0.5:0.1:0.5);
XYZ_Noeud=[reshape(X,[],1),reshape(Y,[],1),reshape(Z,[],1)];
Var=XYZ_Noeud;

%%

if ~with_matlab

F_xyz=naturalInterpolant(XYZ_Noeud,Var,'Sibson');

else
F_x = scatteredInterpolant(XYZ_Noeud(:,1),XYZ_Noeud(:,2),XYZ_Noeud(:,3),Var(:,1),'natural');
F_y = scatteredInterpolant(XYZ_Noeud(:,1),XYZ_Noeud(:,2),XYZ_Noeud(:,3),Var(:,2),'natural');
F_z = scatteredInterpolant(XYZ_Noeud(:,1),XYZ_Noeud(:,2),XYZ_Noeud(:,3),Var(:,3),'natural');
end

for I=1:1000

%% creat dummy point location for interpolation evaluation

XYZ_Point=rand(1000,3)*1.5-0.75;

if ~with_matlab

%% interpolate a filds
% test : Var = XYZ_Noeud ==> interpolated Var on XYZ_Point = XYZ_Point

Var_Int=F_xyz.eval(XYZ_Point);
%Var_Int=F_xyz.Mat_INT*Var;

%% cal error
 
dif=Var_Int-XYZ_Point;

j=0;
ind_p_in=zeros(sum(F_xyz.In_Out),1);
for i=1:size(F_xyz.In_Out,1)
    if F_xyz.In_Out(i)
        j=j+1;
        ind_p_in(j)=i;
    end
end

int_p_out=setdiff(1:size(F_xyz.In_Out,1),ind_p_in);

err=max(max(abs(dif(ind_p_in,:))));

else

%% interpolate a filds
% test : Var = XYZ_Noeud ==> interpolated Var on XYZ_Point = XYZ_Point

Var_Int=[F_x(XYZ_Point),F_y(XYZ_Point),F_z(XYZ_Point)];

%% cal error
 
dif=Var_Int-XYZ_Point;

err=max(max(abs(dif)));

end

%%
end

toc
