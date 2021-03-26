/*
 * Copyright 2020 Uriel Sandoval
 *
 * This file is part of CVXOPT.
 *
 * CVXOPT is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * CVXOPT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "osqp.h"
#include "cs.h"

#include "cvxopt.h"
#include "misc.h"


PyDoc_STRVAR(osqp__doc__,
     "Interface to OSQP LP and QP solver");

 static PyObject *osqp_module;


 void print_csc_matrix2(csc *M, const char *name)
 {
   c_int j, i, row_start, row_stop;
   c_int k = 0;

   // Print name
   c_print("%s :\n", name);

   for (j = 0; j < M->n; j++) {
     row_start = M->p[j];
     row_stop  = M->p[j + 1];

     if (row_start == row_stop) continue;
     else {
       for (i = row_start; i < row_stop; i++) {
         c_print("\t[%3u,%3u] = %.3g\n", (int)M->i[i], (int)j, M->x[k++]);
       }
     }
   }
 }

 void print_dns_matrix2(c_float *M, c_int m, c_int n, const char *name)
 {
   c_int i, j;

   c_print("%s : \n\t", name);

   for (i = 0; i < m; i++) {   // Cycle over rows
     for (j = 0; j < n; j++) { // Cycle over columns
       if (j < n - 1)
         // c_print("% 14.12e,  ", M[j*m+i]);
         c_print("% .3g,  ", M[j * m + i]);

       else
         // c_print("% 14.12e;  ", M[j*m+i]);
         c_print("% .3g;  ", M[j * m + i]);
     }

     if (i < m - 1) {
       c_print("\n\t");
     }
   }
   c_print("\n");
 }

 void print_vec2(c_float *v, c_int n, const char *name) {
   print_dns_matrix2(v, 1, n, name);
 }

static char doc_qp[] =
   "        minimize    (1/2)*x'*P*x + q'*x \n"
   "        subject to  G*x <= h \n"
   "                     A*x = b \n"
   "                                \n"
   "minimize        0.5 x' P x + q' x \n"
   "subject to      l <= A x <= u";

static PyObject *qp(PyObject *self, PyObject *args, PyObject *kwargs){

    matrix *q, *h, *b = NULL, *x = NULL, *y = NULL, *z = NULL;
    spmatrix *P = NULL, *G, *A = NULL, *Pt = NULL;
    PyObject *opts = NULL, *res = NULL;
    int_t i, j, k, m, n, p, exitflag;
    csc *A_o;
    double *l, *u;

    OSQPWorkspace *work;
    OSQPSettings  *settings;
    OSQPData      *data;
    char *kwlist[] = {"q", "G", "h", "A", "b", "P", "options", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOO|OOOO", kwlist, &q,
        &G, &h, &A, &b, &P, &opts)) return NULL;



    if (!(SpMatrix_Check(G) && SP_ID(G) == DOUBLE)){
            PyErr_SetString(PyExc_TypeError, "G must be a sparse 'd' matrix");
            return NULL;
    }
    if ((m = SP_NROWS(G)) <= 0)
        err_p_int("m");
    if ((n = SP_NCOLS(G)) <= 0)
        err_p_int("n");

    if (!Matrix_Check(h) || h->id != DOUBLE)
        err_dbl_mtrx("h");
    if (h->nrows != m || h->ncols != 1){
        PyErr_SetString(PyExc_ValueError, "incompatible dimensions");
        return NULL;
    }

    if (!Matrix_Check(q) || q->id != DOUBLE)
        err_dbl_mtrx("q");
    if (q->nrows != n || q->ncols != 1){
        PyErr_SetString(PyExc_ValueError, "incompatible dimensions");
        return NULL;
    }

    if ((PyObject *) A == Py_None)
        A = NULL;
    if (A){
        if (!(SpMatrix_Check(A) && SP_ID(A) == DOUBLE)){
                PyErr_SetString(PyExc_ValueError, "A must be a sparse 'd' matrix");
                return NULL;
        }
        if ((p = SP_NROWS(A)) < 0)
            err_p_int("p");
        if (SP_NCOLS(A) != n){
            PyErr_SetString(PyExc_ValueError, "incompatible dimensions");
            return NULL;
        }
    }
    else
        p = 0;

    if ((PyObject *) b == Py_None) b = NULL;
    if (b && (!Matrix_Check(b) || b->id != DOUBLE))
        err_dbl_mtrx("b");
    if ((b && (b->nrows != p || b->ncols != 1)) || (!b && p !=0 )){
        PyErr_SetString(PyExc_ValueError, "incompatible dimensions");
        return NULL;
    }

    if (P == Py_None)
        P = NULL;
    if (P) {
        if (!(SpMatrix_Check(P) && SP_ID(P) == DOUBLE)){
                PyErr_SetString(PyExc_ValueError, "P must be a sparse 'd' matrix");
                return NULL;
        }

        if (SP_NCOLS(P) != n){
            PyErr_SetString(PyExc_ValueError, "incompatible dimensions");
            return NULL;
        }

        if (SP_NROWS(P) != n){
            PyErr_SetString(PyExc_ValueError, "incompatible dimensions");
            return NULL;
        }
        /* We transpose the matrix since CVXOPT formulation receives the
         * matrix in a lower triangular form. OSQP needs it in upper triangular */
        Pt = SpMatrix_Trans(P);

    }



    if (!(settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings))))
        return PyErr_NoMemory();


    if (!(data = (OSQPData *)c_malloc(sizeof(OSQPData)))){
        c_free(settings);
        return PyErr_NoMemory();
    }


    /*
    Transform from CVXOPT formulation:

    minimize    (1/2)*x'*P*x + q'*x
    subject to  G*x <= h
                A*x = b

    To OSQP

    minimize     0.5 x' P x + q' x

    subject to      l <= A x <= u


    */

    data->n = n;
    /*  Rows(G) + Rows(A) */
    data->m = m + p;

    if (p > 0) {
        /*
        Transform from CVXOPT formulation:

        minimize    (1/2)*x'*P*x + q'*x
        subject to  G*x <= h
                    A*x = b

        To OSQP form

        minimize     0.5 x' P x + q' x

        subject to      l <= A x <= u

        */

        if (!(A_o = csc_spalloc(data->m, data->n, SP_NNZ(G) + SP_NNZ(A), 1, 0))){
            c_free(settings);
            c_free(data);
            return PyErr_NoMemory();
        }

        k = 0;
        for (j = 0; j < SP_NCOLS(G); j++){
            for (i = SP_COL(G)[j]; i < SP_COL(G)[j + 1]; k++, i++){
                A_o->i[k] = SP_ROW(G)[i];
                A_o->x[k] = SP_VALD(G)[i];
            }
            A_o->p[j] = SP_COL(G)[j];

            for (i = SP_COL(A)[j]; i < SP_COL(A)[j + 1]; k++, i++){
                A_o->i[k] = SP_ROW(A)[i] + SP_NROWS(G);
                A_o->x[k] = SP_VALD(A)[i];
            }
            A_o->p[j] += SP_COL(A)[j];
        }
        A_o->p[j] = SP_NNZ(G) + SP_NNZ(A);

        data->A = A_o;
    }
    else{
        data->A = csc_matrix(m, n, SP_NNZ(G), SP_VALD(G), (c_int *)SP_ROW(G),
                             (c_int *)SP_COL(G));
    }


    l = malloc(data->m * sizeof(double));
    u = malloc(data->m * sizeof(double));

    if (!l || !u){
        c_free(settings);
        c_free(data);
        if (p > 0) csc_spfree(A_o);
        else c_free(data->A);
        free(l);
        free(u);
        return PyErr_NoMemory();
    }

    for (i = 0; i < m; i ++) {
        l[i] = -OSQP_INFTY;
        u[i] = MAT_BUFD(h)[i];
    }

    if (b)
        for (j = 0, i = m; i < m + p; i++, j++) {
            l[i] = u[i] = MAT_BUFD(b)[j];
        }

    data->l = l;
    data->u = u;

    if (P)
        data->P = csc_matrix(data->n, data->n, SP_NNZ(Pt), SP_VALD(Pt),
                             (c_int *)SP_ROW(Pt), (c_int *)SP_COL(Pt));
    else {
        data->P = csc_spalloc(data->n, data->n, 0, 0, 0);
        for (i = 0; i < data->n + 1; i++)
            data->P->p[i] = 0;
    }

    data->q = MAT_BUFD(q);


    osqp_set_default_settings(settings);
    settings->verbose= 1;

    print_csc_matrix2(data->P, "P");
    print_vec2(data->q, data->n, "q");

    print_csc_matrix2(data->A, "A");
    print_vec2(data->l, data->m, "l");
    print_vec2(data->u, data->m, "u");


    // Setup workspace
    exitflag = osqp_setup(&work, data, settings);

    // Solve Problem
    exitflag = osqp_solve(work);


    // We can free the following variables since the solution is stored in work
    if (p > 0) csc_spfree(A_o);
    else c_free(data->A);
    if (P){
        c_free(data->P);
        Py_XDECREF(Pt);
    }
    else
        csc_spfree(data->P);
    c_free(data);
    free(l);
    free(u);
    c_free(settings);

    if (!(res = PyTuple_New(4))){
        osqp_cleanup(work);
        return PyErr_NoMemory();
    }

    x = (matrix *) Matrix_New(n, 1, DOUBLE);
    z = (matrix *) Matrix_New(m, 1, DOUBLE);
    y = (matrix *) Matrix_New(p, 1, DOUBLE);

    if (!x || !z || !y){
        Py_XDECREF(x);
        Py_XDECREF(z);
        Py_XDECREF(y);
        Py_XDECREF(res);
        osqp_cleanup(work);
        return PyErr_NoMemory();
    }

    PyTuple_SET_ITEM(res, 0, (PyObject *) PYSTRING_FROMSTRING(work->info->status));


    if (work->info->status_val == OSQP_SOLVED ||
        work->info->status_val == OSQP_SOLVED_INACCURATE){

        memcpy(MAT_BUFD(x), (double *) work->solution->x, n * sizeof(double));
        memcpy(MAT_BUFD(z), (double *) work->solution->y, m * sizeof(double));
        if (p > 0)
            memcpy(MAT_BUFD(y), (double *) &work->solution->y[m],
                   p * sizeof(double));

    }
    else if (work->info->status_val == OSQP_PRIMAL_INFEASIBLE ||
             work->info->status_val == OSQP_PRIMAL_INFEASIBLE_INACCURATE){
        /* Return the primal certificate */
        memcpy(MAT_BUFD(z), (double *) work->delta_y, m * sizeof(double));
        if (p > 0)
            memcpy(MAT_BUFD(y), (double *) &work->delta_y[m], p * sizeof(double));

    }
    else if (work->info->status_val == OSQP_DUAL_INFEASIBLE ||
             work->info->status_val == OSQP_DUAL_INFEASIBLE_INACCURATE){

        memcpy(MAT_BUFD(x), (double *) work->delta_x, n * sizeof(double));

    }


    PyTuple_SET_ITEM(res, 1, (PyObject *) x);
    PyTuple_SET_ITEM(res, 2, (PyObject *) z);
    PyTuple_SET_ITEM(res, 3, (PyObject *) y);


    osqp_cleanup(work);
    return res;


}




static PyMethodDef osqp_functions[] = {
    {"qp", (PyCFunction) qp, METH_VARARGS|METH_KEYWORDS, doc_qp},
    {NULL}  /* Sentinel */
};

#if PY_MAJOR_VERSION >= 3

static PyModuleDef osqp_module_def = {
    PyModuleDef_HEAD_INIT,
    "osqp",
    osqp__doc__,
    -1,
    osqp_functions,
    NULL, NULL, NULL, NULL
};

PyMODINIT_FUNC PyInit_osqp(void)
{
  if (!(osqp_module = PyModule_Create(&osqp_module_def))) return NULL;
  PyModule_AddObject(osqp_module, "options", PyDict_New());
  if (import_kvxopt() < 0) return NULL;
  return osqp_module;
}

#else

PyMODINIT_FUNC initosqp(void)
{
    osqp_module = Py_InitModule3("kvxopt.osqp", osqp_functions,
        osqp__doc__);
    PyModule_AddObject(osqp_module, "options", PyDict_New());
    if (import_kvxopt() < 0) return;
}

#endif
