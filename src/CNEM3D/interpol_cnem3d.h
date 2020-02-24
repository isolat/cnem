/* This file is part of CNEMLIB.
 
Copyright (C) 2003-2011
Lounes ILLOUL (illoul_lounes@yahoo.fr)
Philippe LORONG (philippe.lorong@ensam.eu)
Arts et Metiers ParisTech, Paris, France
 
CNEMLIB is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CNEMLIB is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

Please report bugs to illoul_lounes@yahoo.fr */

#pragma once

#include "C_Meshless_3d.h"

#include "InterpolParal.h"

#include<vector>
using namespace std;

long interpol_cnem3d
(//IN
size_t Nb_Noeud,
double* Tab_Noeud,
size_t Nb_Tri_Front,
size_t* Tab_Ind_Noeud_Tri_Front,
size_t Nb_Point,
double* Tab_Point,
size_t Type_Appel_Tetgen,
size_t Type_FF,
/*size_t nb_core_for_ff_cal,*/
//OUT
vector<size_t>* P_Vec_Ind_Noeud_New_Old,
vector<size_t>* P_Vec_Ind_Noeud_Old_New,
vector<size_t>* P_Vec_INVNN,
vector<double>* P_Vec_PNVNN,
vector<size_t>* P_Vec_Ind_Point,
vector<size_t>* P_Vec_Nb_Contrib,
vector<size_t>* P_Vec_INV,
vector<double>* P_Vec_Phi,
vector<double>* P_Vec_Gard);

//---------------------------------------------------------------------------//

class C_Interpolator3d {

public:

	C_Interpolator3d(size_t Nb_Noeud, double* Tab_Noeud, size_t Nb_Tri_Front, double* Tab_Ind_Noeud_Tri_Front, size_t Type_Appel_Tetgen, char* bin_path_);
	C_Interpolator3d(size_t Nb_Noeud, double* Tab_Noeud, size_t nb_tri_dm, double* tri_dm, size_t nb_tet_dm, double* tetra_dm, double* neighbor_dm);
	~C_Interpolator3d();

	void Interpolate(size_t Nb_Point, double* Tab_Point, size_t Type_FF);

public:
	vector<size_t> Vec_Ind_Noeud_New_Old;
	vector<size_t> Vec_Ind_Noeud_Old_New;
	vector<size_t> Vec_INVNN;
	vector<double> Vec_PNVNN;
	//
	vector<size_t> Vec_Ind_Point;
	vector<size_t> Vec_Nb_Contrib;
	vector<size_t> Vec_INV;
	vector<double> Vec_Phi;
	vector<double> Vec_Grad;
	bool initialised;
	size_t nb_core_for_ff_cal;
	char* bin_path;

private:

	double* Copy_Tab_Noeud;
	long* Copy_Tab_Ind_Noeud_Tri_Front;
	C_Meshless_3d* PML;
};

class C_int_man
{
public:
	C_int_man() { free_id = 1; }
	~C_int_man() {}
	int add(C_Interpolator3d* p_int) { dict[free_id] = p_int; free_id++; return free_id - 1; }
	C_Interpolator3d* get(int id)
	{
		map<int, C_Interpolator3d*>::iterator it = dict.find(id);
		if (it != dict.end())return (C_Interpolator3d*)(*it).second;
		else NULL;
	}
	int del(int id)
	{
		map<int, C_Interpolator3d*>::iterator it = dict.find(id);
		if (it != dict.end()) { dict.erase(it); return 0; }
		else return 1;
	}
	void print_in_use()
	{
		printf("\n%zd interpolator id in use :\n", dict.size());
		map<int, C_Interpolator3d*>::iterator i;
		for (i = dict.begin(); i != dict.end(); i++)printf("    %d\n", (*i).first);
	}

private:
	int free_id;
	map<int, C_Interpolator3d*> dict;
};
