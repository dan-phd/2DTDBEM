// My Matlab file control functions

#pragma once

#include "Define.h"
#include "matio.h"

int CreateMatFile(mat_t **file_ptr, const char *filename)
{
	// Create .mat file container, and decide which standard to use -
	// choose from MAT_FT_DEFAULT (chooses for you), MAT_FT_MAT5, or MAT_FT_MAT73
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

void InsertMatrix(mat_t **file_ptr, const char *name, MATRIX &D)
{
	//write MATRIX D
	size_t dims[2] = { D.n_rows, D.n_cols };
	int tot = D.n_elem;
	double *tmp = new double[tot];
	for (int i = 0; i < tot; i++)
		tmp[i] = D(i);

	matvar_t *field;
	field = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, tmp, 0);

	if (NULL == field) {
		fprintf(stderr, "Error creating variable for '%s'\n", name);
	}
	else {
		Mat_VarWrite(*file_ptr, field, MAT_COMPRESSION_ZLIB);
		Mat_VarFree(field);
	}
	delete tmp;
}

void InsertCube(mat_t **file_ptr, const char *name, cube &D)
{
	//write cube D
	size_t dims[3] = { D.n_rows, D.n_cols, D.n_slices };
	int tot = D.n_elem;
	double *tmp = new double[tot];
	for (int i = 0; i < tot; i++)
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
	for (int i = 0; i < tot; i++)
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

	for (int i = 0; i < tot; i++)
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

bool ReadGeometry(GEOMETRY& geometry, double& dt, double& c,
	double& num_shapes, double& dual, const char* fname)
{
	mat_t *matfp = NULL;
	matvar_t *matvar = NULL, *tmp = NULL;
	matfp = Mat_Open(fname, MAT_ACC_RDONLY);

	if (NULL == matfp)
	{
		fprintf(stderr, "Error opening file %s\n", fname);
		return false;
	}

	// Get dt
	matvar = Mat_VarRead(matfp, "dt");
	memcpy(&dt, matvar->data, matvar->data_size);

	// Get c
	matvar = Mat_VarRead(matfp, "c");
	memcpy(&c, matvar->data, matvar->data_size);

	// Get number of shapes in this file
	matvar = Mat_VarRead(matfp, "num_shapes");
	memcpy(&num_shapes, matvar->data, matvar->data_size);

	// Get the "dual" flag
	matvar = Mat_VarRead(matfp, "dual");
	memcpy(&dual, matvar->data, matvar->data_size);

	// Get geometry
	matvar = Mat_VarRead(matfp, "boundary");
	if (matvar == NULL)	return false;
	const int num_edges = (int)matvar->dims[0];
	matvar_t* tmpp;
	EDGE edge;
	int index = 0;
	while (index < num_edges)
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


bool ReadScatteredFieldGeometry(MATRIX& M, MATRIX& J, GRID& rho,
	double& material_parameter, UINT N_T_set_by_user, const char* fname)
{

	mat_t *matfp = NULL;
	matvar_t *matvar = NULL, *tmp = NULL;
	matfp = Mat_Open(fname, MAT_ACC_RDONLY);

	// Get c
	matvar = Mat_VarRead(matfp, "material_parameter");
	memcpy(&material_parameter, matvar->data, matvar->data_size);

	// Get M
	matvar = Mat_VarRead(matfp, "M");
	if (matvar == NULL)	return false;
	const int N_V = (int)matvar->dims[0];
	const int N_T = (int)matvar->dims[1];
	if ((int)N_T_set_by_user>N_T)
	{
		fprintf(stderr, "\n\nSpecified number of timesteps exceeds that of input file. \n");
		return false;
	}
	MATRIX tmp_mat(N_V, N_T, fill::zeros);
	int num_elements = N_T * N_V;
	const double *M_data = static_cast<const double*>(matvar->data);
	for (int i = 0; i < num_elements; ++i)
	{
		tmp_mat(i) = M_data[i];
	}
	M = tmp_mat;
	tmp_mat.zeros();

	// Get J
	matvar = Mat_VarRead(matfp, "J");
	if (matvar == NULL)	return false;
	const double *J_data = static_cast<const double*>(matvar->data);
	for (int i = 0; i < num_elements; ++i)
	{
		tmp_mat(i) = J_data[i];
	}
	J = tmp_mat;
	tmp_mat.clear();

	// Get rho
	matvar = Mat_VarRead(matfp, "rho");
	if (matvar == NULL)	return false;
	const int num_points = (int)matvar->dims[0];
	matvar_t* tmpp;
	POINT2D point;
	int index = 0;
	while (index < num_points)
	{
		tmp = Mat_VarGetStructFieldByName(matvar, "x", index);
		memcpy(&point.x, tmp->data, tmp->nbytes);
		tmp = Mat_VarGetStructFieldByName(matvar, "y", index);
		memcpy(&point.y, tmp->data, tmp->nbytes);

		rho.push_back(point);
		tmpp = matvar;
		index++;
	}

	Mat_VarFree(tmpp);
	Mat_Close(matfp);
	return true;
}
