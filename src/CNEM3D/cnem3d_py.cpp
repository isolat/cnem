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

#include <Python.h>
#include "scni_cnem3d.h"
#include "interpol_cnem3d.h"
#include "mesh_cnem3d.h"
#include "utile_py.h"

PyObject* Py_SCNI_CNEM3D(PyObject* self,PyObject* args)
{
    PyObject* Py_Noeud;
    PyObject* Py_Ind_Noeud_Tri_Front;
    PyObject* Py_Type_Appel_Tetgen;
    PyObject* Py_Type_FF;
    PyObject* Py_Sup_NN_GS;
	char *bin_path;
    //PyObject* Py_nb_core_for_gs_cal;
    
    if(!PyArg_ParseTuple(args,"OOOOOs",&Py_Noeud,&Py_Ind_Noeud_Tri_Front,&Py_Type_Appel_Tetgen,&Py_Type_FF,&Py_Sup_NN_GS, &bin_path/*,&Py_nb_core_for_gs_cal*/))
        return (PyObject*)Py_BuildValue("");
    if(!PySequence_Check(Py_Noeud))return (PyObject*)Py_BuildValue("");
    if(!PySequence_Check(Py_Ind_Noeud_Tri_Front))return (PyObject*)Py_BuildValue("");
    if(!PyNumber_Check(Py_Type_Appel_Tetgen))return (PyObject*)Py_BuildValue("");
    if(!PyNumber_Check(Py_Type_FF))return (PyObject*)Py_BuildValue("");
    if(!PyNumber_Check(Py_Sup_NN_GS))return (PyObject*)Py_BuildValue("");
    //if(!PyNumber_Check(Py_nb_core_for_gs_cal))return (PyObject*)Py_BuildValue("");
    
    size_t Size;
    Size=PySequence_Size(Py_Noeud);
    size_t Nb_Noeud=Size/3;
    double* Tab_Noeud=Tuple_F_2_Tab_D(Py_Noeud,Size);

    Size=PySequence_Size(Py_Ind_Noeud_Tri_Front);
    size_t Nb_Tri_Front=Size/3;
    size_t* Tab_Ind_Noeud_Tri_Front=Tuple_L_2_Tab_UL(Py_Ind_Noeud_Tri_Front,Size);
    
    size_t Type_Appel_Tetgen=Num_2_UL(Py_Type_Appel_Tetgen);
    size_t Type_FF=Num_2_UL(Py_Type_FF);
    size_t Sup_NN_GS=Num_2_UL(Py_Sup_NN_GS);
    //size_t nb_core_for_gs_cal=Num_2_UL(Py_nb_core_for_gs_cal);

    size_t Type_Int=0;

    //-----------------------------------------------------------------------//

    printf("\n\nCalling cnem3d...\n\n");

    vector<size_t> Vec_Ind_Noeud_New_Old;
    vector<size_t> Vec_Ind_Noeud_Old_New;
    vector<double> Vec_New_Noeud;
    vector<size_t> Vec_INVNN;
    vector<double> Vec_PNVNN;
    vector<double> Vec_Vol_Cel;
    vector<size_t> Vec_Nb_Contrib;
    vector<size_t> Vec_INV;
    vector<double> Vec_Grad;
    vector<size_t> Vec_Ind_Noeud_New_Tri;
    vector<size_t> Vec_Ind_Noeud_Tet;

    if(scni_cnem3d
       (Nb_Noeud,Tab_Noeud,Nb_Tri_Front,Tab_Ind_Noeud_Tri_Front,
        Type_Appel_Tetgen,Sup_NN_GS,Type_Int,Type_FF, bin_path,/*nb_core_for_gs_cal,*/
        &Vec_Ind_Noeud_New_Old,&Vec_Ind_Noeud_Old_New,&Vec_New_Noeud,&Vec_INVNN,&Vec_PNVNN,
        &Vec_Vol_Cel,&Vec_Nb_Contrib,&Vec_INV,&Vec_Grad,
        &Vec_Ind_Noeud_New_Tri,&Vec_Ind_Noeud_Tet))
        return (PyObject*)Py_BuildValue("");

    //-----------------------------------------------------------------------//
    
    PyObject* Py_Res=PyTuple_New(11);

    PyTuple_SetItem(Py_Res,0,Vec_UL_2_Tuple_L(&Vec_Ind_Noeud_New_Old));
    PyTuple_SetItem(Py_Res,1,Vec_UL_2_Tuple_L(&Vec_Ind_Noeud_Old_New));
    PyTuple_SetItem(Py_Res,2,Vec_D_2_Tuple_F(&Vec_New_Noeud));
    PyTuple_SetItem(Py_Res,3,Vec_UL_2_Tuple_L(&Vec_INVNN));
    PyTuple_SetItem(Py_Res,4,Vec_D_2_Tuple_F(&Vec_PNVNN));
    PyTuple_SetItem(Py_Res,5,Vec_D_2_Tuple_F(&Vec_Vol_Cel));
    PyTuple_SetItem(Py_Res,6,Vec_UL_2_Tuple_L(&Vec_Nb_Contrib));
    PyTuple_SetItem(Py_Res,7,Vec_UL_2_Tuple_L(&Vec_INV));
    PyTuple_SetItem(Py_Res,8,Vec_D_2_Tuple_F(&Vec_Grad));
    PyTuple_SetItem(Py_Res,9,Vec_UL_2_Tuple_L(&Vec_Ind_Noeud_New_Tri));
    PyTuple_SetItem(Py_Res,10,Vec_UL_2_Tuple_L(&Vec_Ind_Noeud_Tet));
    
    //-----------------------------------------------------------------------//

    return Py_Res;
}

//---------------------------------------------------------------------------//

PyObject* Py_MESH_CNEM3D(PyObject* self, PyObject* args)
{
	PyObject* Py_Noeud;
	PyObject* Py_Ind_Noeud_Tri_Front;
	PyObject* Py_Type_Appel_Tetgen;
	char *bin_path;

	if (!PyArg_ParseTuple(args, "OOOs", &Py_Noeud, &Py_Ind_Noeud_Tri_Front, &Py_Type_Appel_Tetgen, &bin_path))
		return (PyObject*)Py_BuildValue("");
	if (!PySequence_Check(Py_Noeud))return (PyObject*)Py_BuildValue("");
	if (!PySequence_Check(Py_Ind_Noeud_Tri_Front))return (PyObject*)Py_BuildValue("");
	if (!PyNumber_Check(Py_Type_Appel_Tetgen))return (PyObject*)Py_BuildValue("");

	size_t Size;
	Size = PySequence_Size(Py_Noeud);
	size_t Nb_Noeud = Size / 3;
	double* Tab_Noeud = Tuple_F_2_Tab_D(Py_Noeud, Size);

	Size = PySequence_Size(Py_Ind_Noeud_Tri_Front);
	size_t Nb_Tri_Front = Size / 3;
	size_t* Tab_Ind_Noeud_Tri_Front = Tuple_L_2_Tab_UL(Py_Ind_Noeud_Tri_Front, Size);

	size_t Type_Appel_Tetgen = Num_2_UL(Py_Type_Appel_Tetgen);

	//-----------------------------------------------------------------------//

	printf("\n\nCalling cnem3d...\n\n");

	vector<size_t> Vec_Ind_Noeud_New_Old;
	vector<size_t> Vec_Ind_Noeud_Old_New;
	vector<double> Vec_New_Noeud;
	vector<size_t> Vec_Ind_Noeud_New_Tri;
	vector<size_t> Vec_Ind_Noeud_Tet;

	if (mesh_cnem3d
	(Nb_Noeud, Tab_Noeud, Nb_Tri_Front, Tab_Ind_Noeud_Tri_Front,
		Type_Appel_Tetgen, bin_path,
		&Vec_Ind_Noeud_New_Old, &Vec_Ind_Noeud_Old_New, &Vec_New_Noeud,
		&Vec_Ind_Noeud_New_Tri, &Vec_Ind_Noeud_Tet))
		return (PyObject*)Py_BuildValue("");

	//-----------------------------------------------------------------------//

	PyObject* Py_Res = PyTuple_New(5);

	PyTuple_SetItem(Py_Res, 0, Vec_UL_2_Tuple_L(&Vec_Ind_Noeud_New_Old));
	PyTuple_SetItem(Py_Res, 1, Vec_UL_2_Tuple_L(&Vec_Ind_Noeud_Old_New));
	PyTuple_SetItem(Py_Res, 2, Vec_D_2_Tuple_F(&Vec_New_Noeud));
	PyTuple_SetItem(Py_Res, 3, Vec_UL_2_Tuple_L(&Vec_Ind_Noeud_New_Tri));
	PyTuple_SetItem(Py_Res, 4, Vec_UL_2_Tuple_L(&Vec_Ind_Noeud_Tet));

	//-----------------------------------------------------------------------//

	return Py_Res;
}

//---------------------------------------------------------------------------//

C_int_man int_man;
//Py_INCREF(Py_None); return Py_None;
PyObject* Py_INTERPOL_CNEM3D_ini_tg(PyObject* self, PyObject* args)
{
	PyObject* Py_Noeud;
	PyObject* Py_Ind_Noeud_Tri_Front;
	PyObject* Py_Type_Appel_Tetgen;
	char *bin_path;
	//PyObject* Py_bin_path;
	//PyObject* Py_nb_core_for_ff_cal;

	if (!PyArg_ParseTuple(args, "OOOs", &Py_Noeud, &Py_Ind_Noeud_Tri_Front, &Py_Type_Appel_Tetgen, &bin_path /*,&Py_nb_core_for_ff_cal*/))
		return (PyObject*)Py_BuildValue("");
	if (!PySequence_Check(Py_Noeud))return (PyObject*)Py_BuildValue("");
	if (!PySequence_Check(Py_Ind_Noeud_Tri_Front))return (PyObject*)Py_BuildValue("");
	if (!PyNumber_Check(Py_Type_Appel_Tetgen))return (PyObject*)Py_BuildValue("");
	//if (!PyString_Check(Py_bin_path))return (PyObject*)Py_BuildValue("");
	//if(!PyNumber_Check(Py_nb_core_for_ff_cal))return (PyObject*)Py_BuildValue("");

	printf("\n\bin_path : %s \n\n",bin_path);

	size_t Size;
	Size = PySequence_Size(Py_Noeud);
	size_t Nb_Noeud = Size / 3;
	double* Tab_Noeud = Tuple_F_2_Tab_D(Py_Noeud, Size);

	Size = PySequence_Size(Py_Ind_Noeud_Tri_Front);
	size_t Nb_Tri_Front = Size / 3;
	double* Tab_Ind_Noeud_Tri_Front = Tuple_F_2_Tab_D(Py_Ind_Noeud_Tri_Front, Size);

	size_t Type_Appel_Tetgen = Num_2_UL(Py_Type_Appel_Tetgen);
	//size_t nb_core_for_ff_cal=Num_2_UL(Py_nb_core_for_ff_cal);

	//-----------------------------------------------------------------------//

	C_Interpolator3d* P_Interpolator = new C_Interpolator3d(Nb_Noeud, Tab_Noeud, Nb_Tri_Front, Tab_Ind_Noeud_Tri_Front, Type_Appel_Tetgen, bin_path);
	if (!P_Interpolator->initialised) { printf("\nfaill to build the vc !\n"); delete P_Interpolator; return (PyObject*)Py_BuildValue(""); }

	printf("\n\ninterpolator initialised \n\n");

	int id = int_man.add(P_Interpolator);

	//-----------------------------------------------------------------------//

	PyObject* Py_Res = PyTuple_New(5);

	PyTuple_SetItem(Py_Res, 0, Vec_UL_2_Tuple_L(&P_Interpolator->Vec_Ind_Noeud_New_Old));
	PyTuple_SetItem(Py_Res, 1, Vec_UL_2_Tuple_L(&P_Interpolator->Vec_Ind_Noeud_Old_New));
	PyTuple_SetItem(Py_Res, 2, Vec_UL_2_Tuple_L(&P_Interpolator->Vec_INVNN));
	PyTuple_SetItem(Py_Res, 3, Vec_D_2_Tuple_F(&P_Interpolator->Vec_PNVNN));
	PyTuple_SetItem(Py_Res, 4, Py_BuildValue("i",id));

	//-----------------------------------------------------------------------//

	free(Tab_Noeud);
	free(Tab_Ind_Noeud_Tri_Front);
	
	//-----------------------------------------------------------------------//

	return Py_Res;
}

PyObject* Py_INTERPOL_CNEM3D_ini_py(PyObject* self, PyObject* args)
{
	PyObject* Py_Noeud;
	PyObject* Py_Tri_DM;
	PyObject* Py_Tet_DM;
	PyObject* Py_Net_DM;

	//PyObject* Py_nb_core_for_ff_cal;

	if (!PyArg_ParseTuple(args, "OOOO", &Py_Noeud, &Py_Tri_DM, &Py_Tet_DM, &Py_Net_DM /*,&Py_nb_core_for_ff_cal*/))
		return (PyObject*)Py_BuildValue("");
	if (!PySequence_Check(Py_Noeud))return (PyObject*)Py_BuildValue("");
	if (!PySequence_Check(Py_Tri_DM))return (PyObject*)Py_BuildValue("");
	if (!PySequence_Check(Py_Tet_DM))return (PyObject*)Py_BuildValue("");
	if (!PySequence_Check(Py_Net_DM))return (PyObject*)Py_BuildValue("");
	//if(!PyNumber_Check(Py_nb_core_for_ff_cal))return (PyObject*)Py_BuildValue("");

	size_t Size;
	Size = PySequence_Size(Py_Noeud);
	size_t Nb_Noeud = Size / 3;
	double* Tab_Noeud = Tuple_F_2_Tab_D(Py_Noeud, Size);

	Size = PySequence_Size(Py_Tri_DM);
	size_t Nb_Tri_DM = Size / 3;
	double* Tab_Tri_DM = Tuple_F_2_Tab_D(Py_Tri_DM, Size);

	Size = PySequence_Size(Py_Tet_DM);
	size_t Nb_Tet_DM = Size / 4;
	double* Tab_Tet_DM = Tuple_F_2_Tab_D(Py_Tet_DM, Size);

	Size = PySequence_Size(Py_Net_DM);
	size_t Nb_Net_DM = Size / 4;
	double* Tab_Net_DM = Tuple_F_2_Tab_D(Py_Net_DM, Size);

	//-----------------------------------------------------------------------//

	C_Interpolator3d* P_Interpolator = new C_Interpolator3d(Nb_Noeud, Tab_Noeud, Nb_Tri_DM, Tab_Tri_DM, Nb_Tet_DM, Tab_Tet_DM, Tab_Net_DM);
	if (!P_Interpolator->initialised) { printf("\nfaill to build the vc !\n"); delete P_Interpolator; return (PyObject*)Py_BuildValue(""); }

	printf("\n\ninterpolator initialised \n\n");

	int id = int_man.add(P_Interpolator);

	//-----------------------------------------------------------------------//

	PyObject* Py_Res = PyTuple_New(1);

	PyTuple_SetItem(Py_Res, 0, Py_BuildValue("i", id));

	//-----------------------------------------------------------------------//

	free(Tab_Noeud);
	free(Tab_Tri_DM);
	free(Tab_Tet_DM);
	free(Tab_Net_DM);

	//-----------------------------------------------------------------------//

	return Py_Res;
}

PyObject* Py_INTERPOL_CNEM3D_del(PyObject* self, PyObject* args)
{
	int id;
	if (!PyArg_ParseTuple(args,"i",&id))return (PyObject*)Py_BuildValue("");
	C_Interpolator3d* P_Interpolator = int_man.get(id);

	if (P_Interpolator == NULL)
	{
		printf("\ninterpolator not initialised !\n");
		Py_INCREF(Py_None); return Py_None;
	}

	printf("\ndeleting interpolator \n");

	delete P_Interpolator;
	int_man.del(id);
	Py_INCREF(Py_None); return Py_None;
}

PyObject* Py_INTERPOL_CNEM3D_int(PyObject* self, PyObject* args)
{
	PyObject* Py_Point;
	PyObject* Py_Type_FF;
	int id;

	if (!PyArg_ParseTuple(args, "iOO", &id, &Py_Point, &Py_Type_FF/*,&Py_nb_core_for_ff_cal*/))
		return (PyObject*)Py_BuildValue("");
	if (!PySequence_Check(Py_Point))return (PyObject*)Py_BuildValue("");
	if (!PyNumber_Check(Py_Type_FF))return (PyObject*)Py_BuildValue("");
	//if(!PyNumber_Check(Py_nb_core_for_ff_cal))return (PyObject*)Py_BuildValue("");

	C_Interpolator3d* P_Interpolator = int_man.get(id);
	if (P_Interpolator == NULL)
	{
		printf("\ninterpolator not initialised !\n");
		Py_INCREF(Py_None); return Py_None;
	}

	printf("\ncalling interpolator... \n");

	size_t Size;
	Size = PySequence_Size(Py_Point);
	size_t Nb_Point = Size / 3;
	double* Tab_Point = Tuple_F_2_Tab_D(Py_Point, Size);

	size_t Type_FF = Num_2_UL(Py_Type_FF);
	//size_t nb_core_for_ff_cal=Num_2_UL(Py_nb_core_for_ff_cal);

	//-----------------------------------------------------------------------//

	P_Interpolator->Interpolate(Nb_Point, Tab_Point, Type_FF/*,nb_core_for_ff_cal,*/);

	//-----------------------------------------------------------------------//

	PyObject* Py_Res = PyTuple_New(5);

	PyTuple_SetItem(Py_Res, 0, Vec_UL_2_Tuple_L(&P_Interpolator->Vec_Ind_Point));
	PyTuple_SetItem(Py_Res, 1, Vec_UL_2_Tuple_L(&P_Interpolator->Vec_Nb_Contrib));
	PyTuple_SetItem(Py_Res, 2, Vec_UL_2_Tuple_L(&P_Interpolator->Vec_INV));
	PyTuple_SetItem(Py_Res, 3, Vec_D_2_Tuple_F(&P_Interpolator->Vec_Phi));
	PyTuple_SetItem(Py_Res, 4, Vec_D_2_Tuple_F(&P_Interpolator->Vec_Grad));

	//-----------------------------------------------------------------------//

	return Py_Res;
}

PyObject* Py_INTERPOL_CNEM3D_use(PyObject* self, PyObject* args)
{
	int_man.print_in_use();
	return (PyObject*)Py_BuildValue("");
}

static PyMethodDef Methods_Module_CNEM3D[]={
    {"SCNI_CNEM3D",Py_SCNI_CNEM3D,METH_VARARGS},
    {"MESH_CNEM3D",Py_MESH_CNEM3D,METH_VARARGS},
	{"INTERPOL_CNEM3D_ini_tg",Py_INTERPOL_CNEM3D_ini_tg,METH_VARARGS},
	{"INTERPOL_CNEM3D_ini_py",Py_INTERPOL_CNEM3D_ini_py,METH_VARARGS},
	{"INTERPOL_CNEM3D_del",Py_INTERPOL_CNEM3D_del,METH_VARARGS},
	{"INTERPOL_CNEM3D_int",Py_INTERPOL_CNEM3D_int,METH_VARARGS},
	{"INTERPOL_CNEM3D_use",Py_INTERPOL_CNEM3D_use,METH_VARARGS},
    {NULL,NULL}};

#if PYTHON_VERSION == 2

extern "C" void initCNEM3D()
{
    Py_InitModule("CNEM3D",Methods_Module_CNEM3D);
}

#elif PYTHON_VERSION == 3

static struct PyModuleDef Module_CNEM3D = {
   PyModuleDef_HEAD_INIT,
   "CNEM3D",   /* name of module */
   NULL, /* module documentation, may be NULL */
   -1,   /* size of per-interpreter state of the module,
            or -1 if the module keeps state in global variables. */
   Methods_Module_CNEM3D
};

extern "C" PyMODINIT_FUNC PyInit_CNEM3D(void)
{
    return PyModule_Create(&Module_CNEM3D);
}

#endif
