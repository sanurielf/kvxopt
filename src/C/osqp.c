/*
 * Copyright 2020 Uriel Sandoval
 *
 * This file is part of KVXOPT.
 *
 * KVXOPT is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * KVXOPT is distributed in the hope that it will be useful,
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


void print_csc_matrix2(csc *M, const char *name) {
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
                c_print("\t[%3u,%3u] = % .3f\n", (int)M->i[i], (int)j, M->x[k++]);
            }
        }
    }
}

void print_dns_matrix2(c_float *M, c_int m, c_int n, const char *name) {
    c_int i, j;

    c_print("%s : \n\t", name);

    for (i = 0; i < m; i++) {   // Cycle over rows
        for (j = 0; j < n; j++) { // Cycle over columns
            if (j < n - 1)
                // c_print("% 14.12e,  ", M[j*m+i]);
                c_print("% .3f,  ", M[j * m + i]);

            else
                // c_print("% 14.12e;  ", M[j*m+i]);
                c_print("% .3f;  ", M[j * m + i]);
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


static PyObject *resize_problem(spmatrix *G, matrix *h, spmatrix *A,
                                matrix *b) {


    /*
    Transform from CVXOPT formulation:

    minimize    (1/2)*x'*P*x + q'*x
    subject to  G*x <= h
                A*x = b

    To OSQP

    minimize     0.5 x' P x + q' x

    subject to      l <= A x <= u


    */
    PyObject *res;
    spmatrix *Anew;
    matrix *l, *u;
    int_t k, i, j, m, n, nnz, p;


    if (A)
        p = SP_NROWS(A);
    else
        p = 0;


    n = SP_NCOLS(G);
    m = SP_NROWS(G) + (p > 0 ? SP_NROWS(A) : 0);
    nnz = SP_NNZ(G) + (p > 0 ? SP_NNZ(A) : 0);
    if (p > 0) {
        if (!(Anew = SpMatrix_New(m, n, nnz, DOUBLE)))
            return NULL;
    }

    k = 0;

    if (p > 0) {
        for (j = 0; j < SP_NCOLS(G); j++) {
            for (i = SP_COL(G)[j]; i < SP_COL(G)[j + 1]; k++, i++) {
                SP_ROW(Anew)[k] = SP_ROW(G)[i];
                SP_VALD(Anew)[k] = SP_VALD(G)[i];
            }
            SP_COL(Anew)[j] = SP_COL(G)[j];

            for (i = SP_COL(A)[j]; i < SP_COL(A)[j + 1]; k++, i++) {
                SP_ROW(Anew)[k] = SP_ROW(A)[i] + SP_NROWS(G);
                SP_VALD(Anew)[k] = SP_VALD(A)[i];
            }
            SP_COL(Anew)[j] += SP_COL(A)[j];
        }

        SP_COL(Anew)[j] = SP_NNZ(G) + SP_NNZ(A);
    }

    l = Matrix_New(m, 1, DOUBLE);
    u = Matrix_New(m, 1, DOUBLE);

    if (!l || !u) {
        if (p > 0)
            Py_DECREF(Anew);
        return NULL;
    }

    for (i = 0; i < m; i ++) {
        MAT_BUFD(l)[i] = -OSQP_INFTY;
        MAT_BUFD(u)[i] = MAT_BUFD(h)[i];
    }

    for (j = 0, i = SP_NROWS(G); i < m; i++, j++)
        MAT_BUFD(l)[i] = MAT_BUFD(u)[i] = MAT_BUFD(b)[j];

    if (!(res = PyTuple_New(3))) {
        Py_DECREF(Anew);
        Py_DECREF(l);
        Py_DECREF(u);
        return NULL;
    }

    if (p > 0)
        PyTuple_SET_ITEM(res, 0, (PyObject *) Anew);
    else {
        Py_INCREF(Py_None);
        PyTuple_SET_ITEM(res, 0, (PyObject *) Py_None);
    }

    PyTuple_SET_ITEM(res, 1, (PyObject *) l);
    PyTuple_SET_ITEM(res, 2, (PyObject *) u);


    return res;

}

static PyObject *solve_problem(spmatrix *P, matrix *q, spmatrix *A, matrix *l,
                               matrix *u) {

    PyObject *res;
    matrix *x, *z;
    int_t i, exitflag;

    OSQPWorkspace *work;
    OSQPSettings  *settings;
    OSQPData      *data;


    if (!(settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings))))
        return NULL;


    if (!(data = (OSQPData *)c_malloc(sizeof(OSQPData)))) {
        c_free(settings);
        return NULL;
    }

    data->m = SP_NROWS(A);
    data->n = SP_NCOLS(A);
    data->A = csc_matrix(data->m, data->n, SP_NNZ(A), (c_float *) SP_VALD(A),
                         (c_int *)SP_ROW(A), (c_int *)SP_COL(A));
    data->l = MAT_BUFD(l);
    data->u = MAT_BUFD(u);


    if (P)
        data->P = csc_matrix(data->n, data->n, SP_NNZ(P), (c_float *) SP_VALD(P),
                             (c_int *)SP_ROW(P), (c_int *)SP_COL(P));
    else {
        data->P = csc_spalloc(data->n, data->n, 0, 0, 0);
        for (i = 0; i < data->n + 1; i++)
            data->P->p[i] = 0;
    }

    data->q = MAT_BUFD(q);

    osqp_set_default_settings(settings);
    settings->verbose = 1;


    //print_csc_matrix2(data->P, "P");
    //print_vec2(data->q, data->n, "q");

    //print_csc_matrix2(data->A, "A");
    //print_vec2(data->l, data->m, "l");
    //print_vec2(data->u, data->m, "u");

    // Setup workspace
    exitflag = osqp_setup(&work, data, settings);

    // Solve Problem
    Py_BEGIN_ALLOW_THREADS;
    exitflag = osqp_solve(work);
    Py_END_ALLOW_THREADS;

    if (P)
        c_free(data->P);
    else
        csc_spfree(data->P);




    if (!(res = PyTuple_New(3))) {

        c_free(data);
        c_free(settings);
        osqp_cleanup(work);
        return NULL;
    }


    x = Matrix_New(data->n, 1, DOUBLE);
    z = Matrix_New(data->m, 1, DOUBLE);

    if (!x || !z) {
        Py_XDECREF(x);
        Py_XDECREF(z);
        Py_XDECREF(res);

        c_free(data);
        c_free(settings);
        osqp_cleanup(work);
        return NULL;
    }



    PyTuple_SET_ITEM(res, 0, (PyObject *) PYSTRING_FROMSTRING(work->info->status));



    if (work->info->status_val == OSQP_SOLVED ||
            work->info->status_val == OSQP_SOLVED_INACCURATE) {

        memcpy(MAT_BUFD(x), (double *) work->solution->x, data->n * sizeof(double));
        memcpy(MAT_BUFD(z), (double *) work->solution->y, data->m * sizeof(double));

        //if (p > 0)
        //    memcpy(MAT_BUFD(y), (double *) &work->solution->y[m],
        //           p * sizeof(double));

    }
    else if (work->info->status_val == OSQP_PRIMAL_INFEASIBLE ||
             work->info->status_val == OSQP_PRIMAL_INFEASIBLE_INACCURATE) {
        /* Return the primal certificate */
        memcpy(MAT_BUFD(z), (double *) work->delta_y, data->m * sizeof(double));
        //if (p > 0)
        //    memcpy(MAT_BUFD(y), (double *) &work->delta_y[m], p * sizeof(double));

    }
    else if (work->info->status_val == OSQP_DUAL_INFEASIBLE ||
             work->info->status_val == OSQP_DUAL_INFEASIBLE_INACCURATE) {

        memcpy(MAT_BUFD(x), (double *) work->delta_x, data->n * sizeof(double));

    }


    PyTuple_SET_ITEM(res, 1, (PyObject *) x);
    PyTuple_SET_ITEM(res, 2, (PyObject *) z);


    c_free(data);
    c_free(settings);

    osqp_cleanup(work);
    return res;

}

static char doc_solve[] =
    "minimize        0.5 x' P x + q' x \n"
    "subject to      l <= A x <= u";

static PyObject *solve(PyObject *self, PyObject *args, PyObject *kwargs) {
    matrix *q, *u, *l;
    spmatrix *P = NULL, *A;
    PyObject *opts = NULL, *res = NULL;
    int_t m, n;

    char *kwlist[] = {"q", "A", "u", "l", "P", "options", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOO|OO", kwlist, &q, &A,
                                     &u, &l, &P, &opts)) return NULL;

    if (!(SpMatrix_Check(A) && SP_ID(A) == DOUBLE)) {
        PyErr_SetString(PyExc_TypeError, "A must be a sparse 'd' matrix");
        return NULL;
    }
    if ((m = SP_NROWS(A)) <= 0)
        err_p_int("m");
    if ((n = SP_NCOLS(A)) <= 0)
        err_p_int("n");

    if (!Matrix_Check(q) || q->id != DOUBLE)
        err_dbl_mtrx("q");
    if (q->nrows != n || q->ncols != 1) {
        PyErr_SetString(PyExc_ValueError, "incompatible dimensions");
        return NULL;
    }

    if (!Matrix_Check(u) || u->id != DOUBLE)
        err_dbl_mtrx("u");
    if (u->nrows != m || u->ncols != 1) {
        PyErr_SetString(PyExc_ValueError, "incompatible dimensions");
        return NULL;
    }

    if (!Matrix_Check(l) || l->id != DOUBLE)
        err_dbl_mtrx("l");
    if (l->nrows != m || l->ncols != 1) {
        PyErr_SetString(PyExc_ValueError, "incompatible dimensions");
        return NULL;
    }

    if ((PyObject *) P == Py_None)
        P = NULL;
    if (P) {
        if (!(SpMatrix_Check(P) && SP_ID(P) == DOUBLE)) {
            PyErr_SetString(PyExc_ValueError, "P must be a sparse 'd' matrix");
            return NULL;
        }

        if (SP_NCOLS(P) != n || SP_NROWS(P) != n) {
            PyErr_SetString(PyExc_ValueError, "incompatible dimensions");
            return NULL;
        }

    }


    res = solve_problem(P, q, A, l, u);

    if (!res)
        return PyErr_NoMemory();

    return res;

}

static char doc_qp[] =
    "        minimize    (1/2)*x'*P*x + q'*x \n"
    "        subject to  G*x <= h \n"
    "                     A*x = b \n"
    "                                \n"
    "minimize        0.5 x' P x + q' x \n"
    "subject to      l <= A x <= u";

static PyObject *qp(PyObject *self, PyObject *args, PyObject *kwargs) {

    matrix *q, *h, *b = NULL, *l = NULL, *u = NULL, *x, *z, *y, *z1;
    spmatrix *P = NULL, *G, *A = NULL, *Pt = NULL, *Anew = NULL;
    PyObject *opts = NULL, *res = NULL, *res_osqp = NULL, *resized = NULL,
              *status;
    int_t m, n, p;


    char *kwlist[] = {"q", "G", "h", "A", "b", "P", "options", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOO|OOOO", kwlist, &q,
                                     &G, &h, &A, &b, &P, &opts)) return NULL;



    if (!(SpMatrix_Check(G) && SP_ID(G) == DOUBLE)) {
        PyErr_SetString(PyExc_TypeError, "G must be a sparse 'd' matrix");
        return NULL;
    }
    if ((m = SP_NROWS(G)) <= 0)
        err_p_int("m");
    if ((n = SP_NCOLS(G)) <= 0)
        err_p_int("n");

    if (!Matrix_Check(h) || h->id != DOUBLE)
        err_dbl_mtrx("h");
    if (h->nrows != m || h->ncols != 1) {
        PyErr_SetString(PyExc_ValueError, "incompatible dimensions");
        return NULL;
    }

    if (!Matrix_Check(q) || q->id != DOUBLE)
        err_dbl_mtrx("q");
    if (q->nrows != n || q->ncols != 1) {
        PyErr_SetString(PyExc_ValueError, "incompatible dimensions");
        return NULL;
    }

    if ((PyObject *) A == Py_None)
        A = NULL;
    if (A) {
        if (!(SpMatrix_Check(A) && SP_ID(A) == DOUBLE)) {
            PyErr_SetString(PyExc_ValueError, "A must be a sparse 'd' matrix");
            return NULL;
        }
        if ((p = SP_NROWS(A)) < 0)
            err_p_int("p");
        if (SP_NCOLS(A) != n) {
            PyErr_SetString(PyExc_ValueError, "incompatible dimensions");
            return NULL;
        }
    }
    else
        p = 0;

    if ((PyObject *) b == Py_None) b = NULL;
    if (b && (!Matrix_Check(b) || b->id != DOUBLE))
        err_dbl_mtrx("b");
    if ((b && (b->nrows != p || b->ncols != 1)) || (!b && p != 0 )) {
        PyErr_SetString(PyExc_ValueError, "incompatible dimensions");
        return NULL;
    }

    if ((PyObject *) P == Py_None)
        P = NULL;
    if (P) {
        if (!(SpMatrix_Check(P) && SP_ID(P) == DOUBLE)) {
            PyErr_SetString(PyExc_ValueError, "P must be a sparse 'd' matrix");
            return NULL;
        }

        if (SP_NCOLS(P) != n) {
            PyErr_SetString(PyExc_ValueError, "incompatible dimensions");
            return NULL;
        }

        if (SP_NROWS(P) != n) {
            PyErr_SetString(PyExc_ValueError, "incompatible dimensions");
            return NULL;
        }
        /* We transpose the matrix since CVXOPT formulation receives the
         * matrix in a lower triangular form. OSQP needs it in upper triangular */
        Pt = SpMatrix_Trans(P);

    }



    if (!(resized = resize_problem(G, h, A, b))) {
        if (P)
            Py_DECREF(Pt);
        return PyErr_NoMemory();
    }

    Anew = (spmatrix *) PyTuple_GET_ITEM(resized, 0);
    l = (matrix *) PyTuple_GET_ITEM(resized, 1);
    u = (matrix *) PyTuple_GET_ITEM(resized, 2);


    res_osqp = solve_problem(P ? Pt : NULL, q, p > 0 ? Anew : G, l, u);


    Py_DECREF(resized);


    if (P)
        Py_DECREF(Pt);

    if (!res_osqp)
        return PyErr_NoMemory();

    if (!(res = PyTuple_New(4)))
        return PyErr_NoMemory();

    status = PyTuple_GET_ITEM(res_osqp, 0);
    x = (matrix *) PyTuple_GET_ITEM(res_osqp, 1);
    z = (matrix *) PyTuple_GET_ITEM(res_osqp, 2);
    Py_INCREF(status);
    Py_INCREF(x);
    Py_INCREF(z);

    Py_DECREF(res_osqp);



    y = (matrix *) Matrix_New(p, 1, DOUBLE);

    PyTuple_SET_ITEM(res, 0, status);
    PyTuple_SET_ITEM(res, 1, (PyObject *) x);

    if (p > 0) {
        z1 = (matrix *) Matrix_New(m, 1, DOUBLE);
        memcpy(MAT_BUFD(z1), MAT_BUFD(z) , m * sizeof(double));
        memcpy(MAT_BUFD(y), &MAT_BUFD(z)[m] , p * sizeof(double));
        Py_DECREF(z);
        PyTuple_SET_ITEM(res, 2, (PyObject *) z1);
    }
    else {
        PyTuple_SET_ITEM(res, 2, (PyObject *) z);
    }

    PyTuple_SET_ITEM(res, 3, (PyObject *) y);


    return res;

}




static PyMethodDef osqp_functions[] = {
    {"qp", (PyCFunction) qp, METH_VARARGS | METH_KEYWORDS, doc_qp},
    {"solve", (PyCFunction) solve, METH_VARARGS | METH_KEYWORDS, doc_solve},
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
