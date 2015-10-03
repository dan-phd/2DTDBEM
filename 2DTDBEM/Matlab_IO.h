// My Matlab file control functions

#pragma once

#include "Define.h"
#include "matio.h"

int CreateMatFile(mat_t **file_ptr, const char *filename)
{
	mat_t *file_p = Mat_CreateVer(filename, NULL, MAT_FT_MAT73);
	*file_ptr = file_p;
	if (NULL == *file_ptr) {
		fprintf(stderr, "Error opening MAT file %s\n", filename);
		return -1;
	}
	return 1;
}

void InsertVar(mat_t **file_ptr, const char *name, double *D)
{
	//write double D
	size_t dims[2] = { 1, 1 };			// default size for a Matlab scalar

	matvar_t *matvar;
	matvar = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, D, 0);

	if (NULL == matvar) {
		fprintf(stderr, "Error creating variable for '%s'\n", name);
	}
	else {
		Mat_VarWrite(*file_ptr, matvar, MAT_COMPRESSION_ZLIB);
		Mat_VarFree(matvar);
	}
}

void InsertCube(mat_t **file_ptr, const char *name, cube &D)
{
	//write cube D
	size_t dims[3] = { D.n_rows, D.n_cols, D.n_slices };
	int tot = D.n_elem;
	double *tmp = new double[tot];
	for (int i = 0; i<tot; i++)
		tmp[i] = D(i);

	matvar_t *field;
	field = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 3, dims, tmp, 0);

	if (NULL == field) {
		fprintf(stderr, "Error creating variable for '%s'\n", name);
	}
	else {
		Mat_VarWrite(*file_ptr, field, MAT_COMPRESSION_ZLIB);
		Mat_VarFree(field);
	}
	delete tmp;
}

void FinishMatFile(mat_t **file_ptr)
{
	Mat_Close(*file_ptr);
}

void CreateStruct(matvar_t **struct_in, const char *mat_name,
	const char **fieldnames, unsigned int nfields)
{
	size_t struct_dims[2] = { 1, 1 };		// default size for a Matlab struct

	*struct_in = Mat_VarCreateStruct(mat_name, 2, struct_dims, fieldnames, nfields);
	if (NULL == *struct_in) {
		fprintf(stderr, "Error creating %s struct\n", mat_name);
		return;
	}
}

void InsertCubeIntoStruct(matvar_t **struct_in, const char *name, cube &D)
{
	//write cube D
	size_t dims[3] = { D.n_rows, D.n_cols, D.n_slices };
	int tot = D.n_elem;
	double *tmp = new double[tot];
	for (int i = 0; i<tot; i++)
		tmp[i] = D(i);

	matvar_t *field;
	field = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 3, dims, tmp, 0);
	Mat_VarSetStructFieldByName(*struct_in, name, 0, field);
	delete tmp;
}

void InsertMatrixIntoStruct(matvar_t **struct_in, const char *name, MATRIX &D)
{
	//size_t dims[2] = { D.n_cols, D.n_rows };		// row major
	size_t dims[2] = { D.n_rows, D.n_cols };		// column major

	int tot = D.n_elem;
	double *tmp = new double[tot];

	for (int i = 0; i<tot; i++)
		tmp[i] = D(i);

	matvar_t *field;
	field = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, tmp, 0);
	Mat_VarSetStructFieldByName(*struct_in, name, 0, field);
	delete tmp;
}

void FinishStruct(mat_t **file_ptr, matvar_t **struct_in)
{
	Mat_VarWrite(*file_ptr, *struct_in, MAT_COMPRESSION_ZLIB);
	Mat_VarFree(*struct_in);
}

bool ReadGeometry(GEOMETRY& geometry, double& dt, double& mu, double& eps, double& num_shapes, const char* fname)
{
	mat_t *matfp = NULL;
	matvar_t *matvar = NULL, *tmp = NULL;
	matfp = Mat_Open(fname, MAT_ACC_RDONLY);

	if (NULL == matfp)
	{
		fprintf(stderr, "Error opening %s file %s\n", fname);
		return false;
	}

	// Get dt
	matvar = Mat_VarRead(matfp, "dt");
	memcpy(&dt, matvar->data, matvar->data_size);

	// Get mu
	matvar = Mat_VarRead(matfp, "mu");
	memcpy(&mu, matvar->data, matvar->data_size);

	// Get eps
	matvar = Mat_VarRead(matfp, "eps");
	memcpy(&eps, matvar->data, matvar->data_size);

	// Get number of shapes in this file
	matvar = Mat_VarRead(matfp, "num_shapes");
	memcpy(&num_shapes, matvar->data, matvar->data_size);

	// Get geometry
	matvar = Mat_VarRead(matfp, "geometry");
	if (matvar == NULL)	return false;
	const int num_edges = (int)matvar->dims[0];
	matvar_t* tmpp;
	EDGE edge;
	int index = 0;
	while (index<num_edges)
	{
		tmp = Mat_VarGetStructFieldByName(matvar, "he_idx", index);
		memcpy(&edge.he_idx, tmp->data, tmp->nbytes);
		tmp = Mat_VarGetStructFieldByName(matvar, "a", index);
		memcpy(&edge.a, tmp->data, tmp->nbytes);
		tmp = Mat_VarGetStructFieldByName(matvar, "b", index);
		memcpy(&edge.b, tmp->data, tmp->nbytes);
		tmp = Mat_VarGetStructFieldByName(matvar, "l", index);
		memcpy(&edge.l, tmp->data, tmp->nbytes);
		tmp = Mat_VarGetStructFieldByName(matvar, "t", index);
		memcpy(&edge.t, tmp->data, tmp->nbytes);
		tmp = Mat_VarGetStructFieldByName(matvar, "n", index);
		memcpy(&edge.n, tmp->data, tmp->nbytes);
		tmp = Mat_VarGetStructFieldByName(matvar, "shape", index);
		memcpy(&edge.shape, tmp->data, tmp->nbytes);

		geometry.push_back(edge);
		tmpp = matvar;
		index++;
	}

	Mat_VarFree(tmpp);
	Mat_Close(matfp);
	return true;
}
