
classdef naturalInterpolant < handle
   properties (Access=public) 
      Points(:,3) double {mustBeReal,mustBeNonNan,mustBeNonsparse}
      Values(:,3) double {mustBeReal,mustBeNonNan,mustBeNonsparse}
      Locations(:,3) double {mustBeReal,mustBeNonNan,mustBeNonsparse}
      Constraints(:,3) double {mustBeReal,mustBeInteger,mustBeNonNan,mustBeNonsparse}
      Method {mustBeMember(Method,{'Sibson','Laplace','SibsonRaw'})}='Sibson'
      Mesher {mustBeMember(Mesher,{'tetgen','matlab'})}='matlab'
   end 
   properties (Access=private)
      Points_changed=false
      Values_changed=false
      Constraints_changed=false
      Locations_changed=false
      Method_changed=false
      Mesher_changed=false
   end
   methods
       function obj = set.Points(obj,Points)
           obj.Points_changed=true;
           obj.Points=Points;
       end
       function obj = set.Values(obj,Values)
           obj.Values_changed=true;
           obj.Values=Values;
       end
       function obj = set.Locations(obj,Locations)
           obj.Locations_changed=true;
           obj.Locations=Locations;
       end
       function obj = set.Method(obj,Method)
           obj.Method_changed=true;
           obj.Method=Method;
       end
       function obj = set.Mesher(obj,Mesher)
           obj.Mesher_changed=true;
           obj.Mesher=Mesher;
       end
       function obj = set.Constraints(obj,Constraints)
           obj.Constraints_changed=true;
           obj.Constraints=Constraints;
       end
       
       function obj=naturalInterpolant(Points,Values,Method,Constraints,Mesher)
           if nargin >= 3
               obj.Points=Points;
               obj.Values=Values;           
               obj.Method=Method;
           end
           
           if nargin >=4 
               obj.Constraints=Constraints;
           end
           
           if nargin >=5
               obj.Mesher=Mesher;
           end
       end
       
       function Evaluation=eval(obj,Locations)
           
           if nargin == 2
               obj.Locations=Locations;
           end
               
           remesh=false;
           if obj.Points_changed
              remesh=true;
              obj.Points_changed=false;
           end
           if obj.Mesher_changed
               remesh=true;
               obj.Mesher_changed=false;
           end
           if obj.Constraints_changed
               remesh=true;
               obj.Constraints_changed=false;
           end
           if remesh
               if size(obj.Points,1) == 0
                   warning('Points not set !')
                   return 
               end
               obj.init();
           end
           
           reinterpol=false;
           if obj.Locations_changed
               reinterpol=true;
               obj.Locations_changed=false;
           end
           if obj.Method_changed
               reinterpol=true;
               obj.Method_changed=false;
           end
           if reinterpol
               if size(obj.Locations,1) == 0
                   warning('Locations not set !')
                   return 
               end
               obj.set_point();
           end
           
           obj.interpolate();
           Evaluation=obj.Var_Int;
       end
       
%        function ret=subsref(obj,Locations)
%            if Locations.type()=='()'
%                obj.Locations=Locations.subs{1}
%                return obj.eval()
%            elseif Locations.type()=='.'
%                 % ???
%            end
%        end
	end
%    methods (Attributes)
%        function obj = methodName(obj,arg2,...)
%           ...
%        end
%     end
%     events (Attributes) 
%        EventName
%     end

    properties(Access=public)
        In_Out
        Nb_V
        Var_Int
        Mat_INT
        Mat_INT_NN
        Mat_GradX
        Mat_GradY
        Mat_GradZ
    end
    properties(Access=private)
        
        nb_noeud_ini
        
        IN_New_Old
        IN_Old_New
        INV_NN
        PNV_NN
        
        my_id
        
        point_set
        
        Type_Call_Tet=0
        
        Type_FF=containers.Map({'Sibson','Laplace','TetLinear','SibsonRaw'},{0,1,2,3})
        
    end
    methods (Access=private)
        function init(obj)
            
            XYZ_Noeud=obj.Points;
            IN_Tri_Ini=obj.Constraints;
            
            obj.nb_noeud_ini=size(XYZ_Noeud,1);
            obj.point_set=false;
            obj.my_id=0;
            
            if obj.Mesher=='tetgen' 
                
                if size(IN_Tri_Ini,1)~=0
                    fid=fopen('bbox_03','r');
                    nb_node_bbox=fread(fid,1,'uint32');
                    nb_tri_bbox=double(fread(fid,1,'uint32'));
                    xyz_node_bbox=fread(fid,nb_node_bbox*3,'float32');
                    in_tri_bbox=fread(fid,nb_tri_bbox*3,'uint32');
                    fclose(fid);
                    xyz_node_bbox=reshape(xyz_node_bbox,3,nb_node_bbox)';
                    in_tri_bbox=reshape(in_tri_bbox,3,nb_tri_bbox)'+1;

                    max_XYZ_Noeud=max(XYZ_Noeud,[],1);
                    min_XYZ_Noeud=min(XYZ_Noeud,[],1);
                    O_dbox=(max_XYZ_Noeud+min_XYZ_Noeud)/2;
                    L_dbox=max(max_XYZ_Noeud-min_XYZ_Noeud);
                    R_dbox=sqrt(3)*L_dbox/2.;
                    coef=1.5;

                    xyz_node_bbox=xyz_node_bbox*(coef*R_dbox)+repmat(O_dbox,size(xyz_node_bbox,1),1);
                    IN_Tri_Ini=[IN_Tri_Ini;in_tri_bbox+size(XYZ_Noeud,1)];
                    XYZ_Noeud=[XYZ_Noeud;xyz_node_bbox];
                end

                %a=sparse(1,1);s=whos('a');
                %if s.bytes==20 %32bit
                %   IN_Tri_Ini=uint32(IN_Tri_Ini);
                %else %64bit
                %   IN_Tri_Ini=uint64(IN_Tri_Ini);
                %end
                
                %S = dbstack('-completenames');
                %my_path=S(1).file;
                [folder, name, ext] = fileparts(which('cnem3d'));
                my_dir=[folder filesep];

                [IN_New_Old,IN_Old_New,INV_NN,PNV_NN,my_id]=...
                cnem3d(2,XYZ_Noeud',IN_Tri_Ini'-1,obj.Type_Call_Tet,my_dir);

                obj.IN_New_Old=double(IN_New_Old)';
                obj.IN_Old_New=double(IN_Old_New)';
                obj.INV_NN=double(INV_NN)'+1;
                obj.PNV_NN=PNV_NN';
                obj.my_id=my_id;
            else
                DT = delaunayTriangulation(XYZ_Noeud(:,1),XYZ_Noeud(:,2),XYZ_Noeud(:,3));
                T=DT.ConnectivityList;
                %T=circshift(T,1,2);
                N=DT.neighbors;
                %N=circshift(N,2,2);
                N(isnan(N)) = 0;
                F=DT.convexHull;
                %DT.Points
                
                [my_id]=...
                cnem3d(3,XYZ_Noeud',F'-1,T'-1,N'-1);
                
                obj.my_id=my_id;
                obj.IN_New_Old=(1:obj.nb_noeud_ini)';
                obj.IN_Old_New=(1:obj.nb_noeud_ini)';
                obj.INV_NN=zeros(0,3);
                obj.PNV_NN=zeros(0,3);
            end
        end
        
        function set_point(obj)
                
            XYZ_Point=obj.Locations;
            Type_FF=obj.Type_FF(obj.Method);
            
            obj.point_set=true;
            
            [Ind_Point,Nb_V,Ind_V,FF,Grad]=...
            cnem3d(4,obj.my_id,XYZ_Point',Type_FF);

            Ind_Point=double(Ind_Point)'+1;
            Nb_V=double(Nb_V)';
            Ind_V=double(Ind_V)'+1;

            nb_noeud_ini_new=size(obj.IN_Old_New,1);

            new_ind_new_noeud=zeros(size(obj.INV_NN,1),1);
            nb_new_noeud=size(obj.INV_NN,1);
            INV_NN_new=zeros(nb_new_noeud,3);
            PNV_NN_new=zeros(nb_new_noeud,3);
            I=0;
            for i=1:size(Ind_V,1)
                if Ind_V(i)<=nb_noeud_ini_new
                    Ind_V(i)=obj.IN_Old_New(Ind_V(i));
                else
                    j=Ind_V(i)-nb_noeud_ini_new;
                    if new_ind_new_noeud(j)==0
                        I=I+1;
                        new_ind_new_noeud(j)=obj.nb_noeud_ini+I;
                        INV_NN_new(I,:)=obj.IN_Old_New(obj.INV_NN(j,:))';
                        PNV_NN_new(I,:)=obj.PNV_NN(j,:);
                    end
                    Ind_V(i)=new_ind_new_noeud(j);
                end
            end
            INV_NN=INV_NN_new(1:I,:);
            PNV_NN=PNV_NN_new(1:I,:);

            INT_NN=[];
            if I~=0
                i_NN=reshape(repmat(1:I,3,1),3*I,1);
                INV_NN=reshape(INV_NN',3*I,1);
                PNV_NN=reshape(PNV_NN',3*I,1);
                INT_NN=sparse(i_NN,INV_NN,PNV_NN,I,obj.nb_noeud_ini);
            end

            Inv_Ind_Point=(1:size(Ind_Point,1));
            Inv_Ind_Point(Ind_Point)=Inv_Ind_Point;

            i_Point=zeros(1,size(Ind_V,1));
            k=1;
            for i=1:size(Nb_V,1)
                for j=1:Nb_V(i)
                    i_Point(k)=i;
                    k=k+1;
                end
            end

            Mat_INT=sparse(i_Point',Ind_V,FF',size(XYZ_Point,1),obj.nb_noeud_ini+I);
            obj.Mat_INT=Mat_INT(Inv_Ind_Point,:);
            
            if size(Grad,2)~=0
                Mat_GradX=sparse(i_Point',Ind_V,Grad(1,:)',size(XYZ_Point,1),obj.nb_noeud_ini+I);
                obj.Mat_GradX=Mat_GradX(Inv_Ind_Point,:);

                Mat_GradY=sparse(i_Point',Ind_V,Grad(2,:)',size(XYZ_Point,1),obj.nb_noeud_ini+I);
                obj.Mat_GradY=Mat_GradY(Inv_Ind_Point,:);

                Mat_GradZ=sparse(i_Point',Ind_V,Grad(3,:)',size(XYZ_Point,1),obj.nb_noeud_ini+I);
                obj.Mat_GradZ=Mat_GradZ(Inv_Ind_Point,:);
            end

            Nb_V(Ind_Point)=Nb_V;
            obj.Nb_V=double(Nb_V);
            obj.In_Out=logical(Nb_V);
            
            obj.Mat_INT_NN=INT_NN;
        end
        
        function delete(obj)
            if obj.my_id ~= 0   
                cnem3d(5,obj.my_id);
            end
        end
        
        function interpolate(obj)
            Var=obj.Values;
            Var_Add_Node=[];
            if size(obj.Mat_INT_NN,1)~=0
                Var_Add_Node=obj.Mat_INT_NN*Var;
            end 
            obj.Var_Int=obj.Mat_INT*[Var;Var_Add_Node];
        end
        
        function Mat=mat_interpol_glob(obj)
            if size(obj.Mat_INT_NN,1)~=0
                n=size(obj.Mat_INT_NN,2);
                m=size(obj.Mat_INT,2);
                Mat=obj.Mat_INT(:,1:n)+obj.Mat_INT(:,n+1:m)*obj.Mat_INT_NN;
            else
                Mat=obj.Mat_INT;
            end
        end
    end
    methods(Static)
        function in_use()
            cnem3d(6);
        end 
    end
end

       

       
       