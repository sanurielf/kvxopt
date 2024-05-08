/*
 * Copyright 2021 Uriel Sandoval
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
#include "kvxopt.h"
#include "gurobi_c.h"
#include "misc.h"

// Status codes from
// https://www.gurobi.com/documentation/9.1/refman/optimization_status_codes.html#sec:StatusCodes
const char *STATUS_CODES[] = {"",
                              "loaded",
                              "optimal",
                              "infeasible",
                              "inf_or_unbd",
                              "unbounded",
                              "cutoff",
                              "iteration_limit",
                              "node_limit",
                              "time_limit",
                              "solution_limit",
                              "interrupted",
                              "numeric",
                              "suboptimal",
                              "inprogress",
                              "user_obj_limit"};

PyDoc_STRVAR(gurobi__doc__, "Interface to Gurobi LP and QP solver");

static PyObject *gurobi_module;

static int solve_problem(spmatrix *P, matrix *_q, matrix *G_l, spmatrix *G, 
                         matrix *G_u, spmatrix *A, matrix *b, matrix *x_l, 
                         matrix *x_u, PyObject *opts, PyObject **res) {
    /* Solve a QP/LP problem in the following form:

    minimize    (1/2)*x'*P*x + q'*x
    subject to  G_l <= G*x <= G_u
                A*x = b
                x_l <= x <= x_u
    */

    PyObject *key, *value;
    matrix *x, *z;
    int_t error = 0, n, i, j, p, q, nz, n_const = 0, n_nnz = 0, pos = 0;
    double  *cval = NULL, *P_v = NULL;
    int *cbeg = NULL, *cind = NULL, *buf = NULL, status, *P_i = NULL, *P_j = NULL;
    char *c_sense = NULL;

    GRBenv *env = NULL;
    GRBmodel *mod = NULL;

    if ((error = GRBemptyenv(&env))) goto CLEAN;

    if ((error = GRBsetstrparam(env, "LogFile", "gurobi_kvxopt.log")))
        goto CLEAN;

    if ((error = GRBstartenv(env))) goto CLEAN;

    /* Here we detect if any user defined options are available through
     * the module options or a dictionary. Otherwise, use standard
     * settings
     */
    if (!(opts && PyDict_Check(opts)))
        opts = PyObject_GetAttrString(gurobi_module, "options");
    if (!opts || !PyDict_Check(opts)) {
        PyErr_SetString(PyExc_AttributeError,
                        "missing gurobi.options dictionary");   
        error = 1;
        goto CLEAN;
    }

    while (PyDict_Next(opts, &pos, &key, &value)) {
        if (PYINT_CHECK(value)) {
            error =
                GRBsetintparam(env, PyStr_AsString(key), PYINT_AS_LONG(value));
        } else if (PyFloat_Check(value)) {
            error = GRBsetdblparam(env, PyStr_AsString(key), PyFloat_AsDouble(value));
        } else {
            PyErr_SetString(PyExc_ValueError,
                            "Option value must be integer or float number");
            error = 1;
            goto CLEAN;
        }
        if (error == 10007) {
            PyErr_SetString(PyExc_ValueError, "Invalid option");
            goto CLEAN;
        } else if (error) {
            PyErr_SetString(PyExc_ValueError,
                            "Unkown error when setting option");
            goto CLEAN;
        }
    }

    n = MAT_NROWS(x_l);

    if (_q) {
        if ((error = GRBnewmodel(env, &mod, "KVXOPT", n, MAT_BUFD(_q), 
                                MAT_BUFD(x_l), MAT_BUFD(x_u),
                                NULL, NULL)))
            goto CLEAN;
    }
    else { 
        if ((error = GRBnewmodel(env, &mod, "KVXOPT", n, NULL,
                             MAT_BUFD(x_l), MAT_BUFD(x_u), NULL, NULL)))
            goto CLEAN;
    }

    if ((error = GRBupdatemodel(mod))) goto CLEAN;
    if ((error = GRBsetintattr(mod, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE)))
        goto CLEAN;

    // Add constraints
    // We add first the inequality constraints:
    // G*x <= h
    // but alloc the workspace to the maximum number of constraints
    // to be added
    if (G) {
        n_const = SP_NROWS(G);
        n_nnz = SP_NNZ(G);
        
        cbeg = calloc(n_const, sizeof(int_t));
        c_sense = malloc((n_const) * sizeof(char));
        cind = calloc(n_nnz, sizeof(int_t));
        cval = calloc(n_nnz, sizeof(double));
        buf = calloc(n_const, sizeof(int_t));  // workspace


        if (!cbeg || !c_sense || !cind || !cval || !buf) {
            error = 100;
            goto CLEAN;
        }

        // Get compressed sparse row (CSR) format of G and
        // add constraints
        // Run through matrix and count number of elms in each row
        for (i = 0; i < SP_NNZ(G); i++) buf[SP_ROW(G)[i]]++;

        // Cummulative sum to get row pointers (cbeg)
        // Also build the sense matrix
        
        nz = 0;
        for (i = 0; i < SP_NROWS(G); i++) {
            cbeg[i] = nz;
            nz += buf[i];
            buf[i] = cbeg[i];
            c_sense[i] = GRB_LESS_EQUAL;
        }

        for (j = 0; j < SP_NCOLS(G); j++) {
            for (p = SP_COL(G)[j]; p < SP_COL(G)[j + 1]; p++) {
                cind[q = buf[SP_ROW(G)[p]]++] = j;
                cval[q] = SP_VALD(G)[p];
            }
        }

        if (G_u){
            // G * x  <= G_u
            for (i = 0; i < SP_NROWS(G); i++) c_sense[i] = GRB_LESS_EQUAL;

            if ((error = GRBaddconstrs(mod, SP_NROWS(G), SP_NNZ(G), cbeg, cind, cval,
                                    c_sense, MAT_BUFD(G_u), NULL)))
                goto CLEAN;
        }

        if (G_l) {
            // G * x >= G_l

            for (i = 0; i < SP_NROWS(G); i++) c_sense[i] = GRB_GREATER_EQUAL;

            if ((error = GRBaddconstrs(mod, SP_NROWS(G), SP_NNZ(G), cbeg, cind,
                                    cval, c_sense, MAT_BUFD(G_l), NULL)))
                goto CLEAN;
        }

        if ((error = GRBupdatemodel(mod))) goto CLEAN;
    }

    if (A) {
        // Previously allocated workspace
        n_const = SP_NROWS(A);
        n_nnz = SP_NNZ(A);
        if (G) {
            if (SP_NROWS(A) > SP_NROWS(G)) {
                // Need to increase size of workspace
                cbeg = realloc(cbeg, n_const * sizeof(int_t));
                c_sense = realloc(c_sense, n_const * sizeof(char));
                buf = realloc(buf, n_const * sizeof(int_t));
                if (!cbeg || !c_sense || !buf) {
                    error = 100;
                    goto CLEAN;
                }
            }

            if (SP_NNZ(A) > SP_NNZ(G)) {
                // Need to increase size of workspace
                cind = realloc(cind, n_nnz * sizeof(int_t));
                cval = realloc(cval, n_nnz * sizeof(double));
                if (!cind || !cval) {
                    error = 100;
                    goto CLEAN;
                }
            }
            for (i = 0; i < n_const; i++) {
                buf[i] = 0;
                cbeg[i] = 0;
            }
            for (i = 0; i < n_nnz; i++) {
                cind[i] = 0;
                cval[i] = 0;
            }
        }
        else{
            // First time we need this workspace
            cbeg = calloc(n_const, sizeof(int_t));
            c_sense = malloc((n_const) * sizeof(char));
            cind = calloc(n_nnz, sizeof(int_t));
            cval = calloc(n_nnz, sizeof(double));
            buf = calloc(n_const, sizeof(int_t));
        }


        // Get compressed sparse row (CSR) format of A and
        // add constraints
        // Run through matrix and count number of elms in each row
        for (i = 0; i < SP_NNZ(A); i++) buf[SP_ROW(A)[i]]++;

        // Cummulative sum to get row pointers (cbeg)
        // Also build the sense matrix
        nz = 0;
        for (i = 0; i < SP_NROWS(A); i++) {
            cbeg[i] = nz;
            nz += buf[i];
            buf[i] = cbeg[i];
            c_sense[i] = GRB_EQUAL;
        }

        for (j = 0; j < SP_NCOLS(A); j++) {
            for (p = SP_COL(A)[j]; p < SP_COL(A)[j + 1]; p++) {
                cind[q = buf[SP_ROW(A)[p]]++] = j;
                cval[q] = SP_VALD(A)[p];
            }
        }

        if ((error = GRBaddconstrs(mod, n_const, n_nnz, cbeg, cind,
                                cval, c_sense, MAT_BUFD(b), NULL)))
            goto CLEAN;

        if ((error = GRBupdatemodel(mod))) goto CLEAN;
    }

    if (P) {
        P_i = calloc(SP_NNZ(P), sizeof(int));
        P_j = calloc(SP_NNZ(P), sizeof(int));
        P_v = calloc(SP_NNZ(P), sizeof(double));

        for (i = 0; i < SP_NNZ(P); i++) {
            P_i[i] = SP_ROW(P)[i];
        }

        for (i = 0; i < SP_NCOLS(P); i++) {
            for (j = SP_COL(P)[i]; j < SP_COL(P)[i + 1]; j++)
                P_j[j] = i;
        }

        for (i = 0; i < SP_NNZ(P); i++) {
            P_v[i] = SP_VALD(P)[i] / 2;

            // Since CVXOPT formulation only takes into account one triangular
            // part of P we duplicate it and set to zero the other part. 
            //  we assume that P[i,j] == P[j,i]
            if (P_i[i] < P_j[i])
                P_v[i] = 0.0;
            if (P_i[i] > P_j[i]) 
                P_v[i] *= 2.0;
        }

        if ((error =
                 GRBaddqpterms(mod, SP_NNZ(P), P_i, P_j, P_v)))
            goto CLEAN;
    }

    GRBwrite(mod, "output.lp");

    if ((error = GRBoptimize(mod))) goto CLEAN;

    if ((error = GRBgetintattr(mod, "Status", &status))) goto CLEAN;

    n_const = 0;
    if (G) n_const += SP_NROWS(G);
    if (A) n_const += SP_NROWS(A);
    x = Matrix_New(n, 1, DOUBLE);
    z = Matrix_New(n_const, 1, DOUBLE);

    if (!x || !z) {
        Py_XDECREF(x);
        Py_XDECREF(z);
        error = 100;
        goto CLEAN;
    }

    if (!(*res = PyTuple_New(3))) {
        Py_XDECREF(x);
        Py_XDECREF(z);
        error = 100;
        goto CLEAN;
    }



    PyTuple_SET_ITEM(*res, 0,
                     (PyObject *)PYSTRING_FROMSTRING(STATUS_CODES[status]));

    GRBgetdblattrarray(mod, "X", 0, n, MAT_BUFD(x));
    PyTuple_SET_ITEM(*res, 1, (PyObject *)x);

    GRBgetdblattrarray(mod, "Pi", 0, n_const, MAT_BUFD(z));
    for (i = 0; i < n_const; i++) MAT_BUFD(z)[i] *= -1;
    PyTuple_SET_ITEM(*res, 2, (PyObject *)z);

CLEAN:
    GRBfreemodel(mod);
    GRBfreeenv(env);

    free(cbeg);
    free(c_sense);
    free(cind);
    free(cval);
    free(buf);
    if (P) {
        free(P_i);
        free(P_j);
        free(P_v);
    }

    if (error) {
        printf("Error %ld\n", error);
        if (error > 1000) printf("%s\n", GRBgeterrormsg(env));
        return error;
    }

    return 0;
}

static char doc_solve[] =
    "minimize        0.5 x' P x + q' x \n"
    "subject to      G_l <= G x <= G_u \n"
    "                       A x = b \n"
    "               x_l <=  x <= x_u\n";

static PyObject *solve(PyObject *self, PyObject *args, PyObject *kwargs) {
    matrix *q = NULL, *G_l = NULL, *G_u = NULL, *x_l = NULL, *x_u = NULL, 
            *b = NULL;
    spmatrix *G = NULL, *A = NULL, *P = NULL;
    PyObject *opts = NULL, *res = NULL;
    int_t i, m, n = -1, n1, p,  error;

    char *kwlist[] = {"q", "G_l", "G", "G_u", "A", "b", "P", "x_l", "x_u", 
                      "options",  NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|OOOOOOOOOO", kwlist, &q, 
                                     &G_l, &G, &G_u, &A, &b, &P, &x_l, &x_u,
                                     &opts))
        return NULL;

    // Since all the arguments are optional we need, user can use positional 
    // arguments and pass None. Then we do Py_None -> NULL

    // G_l <= G *x <= G_u
    if ((PyObject *) G == Py_None) G = NULL;
    if ((PyObject *) G_l == Py_None) G_l = NULL;
    if ((PyObject *) G_u == Py_None) G_u = NULL;
    if (G){
        if (!(SpMatrix_Check(G) && SP_ID(G) == DOUBLE)) {
            PyErr_SetString(PyExc_TypeError, "G must be a sparse 'd' matrix");
            return NULL;
        }


        if (!G_l || !G_u){
            PyErr_SetString(PyExc_TypeError, 
                            "At least one bound matrix must be provided for G");
            return NULL;
        }
        // Check G constraints dimenssions
        if ((m = SP_NROWS(G)) <= 0) err_p_int("m");
        if ((n = SP_NCOLS(G)) <= 0) err_p_int("n");

        if (G_l){
            if (!Matrix_Check(G_l) || G_l->id != DOUBLE) err_dbl_mtrx("G_l");
            if (G_l->nrows != m || G_l->ncols != 1) {
                PyErr_SetString(PyExc_ValueError,
                                "incompatible dimensions: nrows(G) != "
                                "nrows(G_l) or ncols(G) !=1 ");
                return NULL;
            }
        }
        if (G_u) {
            if (!Matrix_Check(G_u) || G_u->id != DOUBLE) err_dbl_mtrx("G_u");
            if (G_u->nrows != m || G_u->ncols != 1) {
                PyErr_SetString(PyExc_ValueError,
                                "incompatible dimensions: nrows(G) != "
                                "nrows(G_u) or ncols(G) !=1 ");
                return NULL;
            }
        }
    }
    // A *x = b
    if ((PyObject *) A == Py_None) A = NULL;
    if ((PyObject *) b == Py_None) b = NULL;
    if (A) {
        if (!(SpMatrix_Check(A) && SP_ID(A) == DOUBLE)) {
            PyErr_SetString(PyExc_TypeError, "A must be a sparse 'd' matrix");
            return NULL;
        }
        if (!b) {
            PyErr_SetString(PyExc_TypeError,
                            "b vector missing");
            return NULL;
        }
        // Check G constraints dimenssions
        if ((p = SP_NROWS(A)) <= 0) err_p_int("p");
        if ((n1 = SP_NCOLS(A)) <= 0) err_p_int("n");

        if (n >= 0 && n !=  n1) {
            PyErr_SetString(PyExc_ValueError,
                            "incompatible dimensions: cols(G) != cols(A)");
            return NULL;
        }
        if (n < 0) n = n1;
        if (!Matrix_Check(b) || b->id != DOUBLE) err_dbl_mtrx("b");
        if (b->nrows != p || b->ncols != 1) {
            PyErr_SetString(PyExc_ValueError,
                            "incompatible dimensions: nrows(b) != nrows(A) or "
                            "ncols(b) != 1");
            return NULL;
        }
    }

    if ((PyObject *)q == Py_None) q = NULL;
    if (q) {
        if (!Matrix_Check(q) || q->id != DOUBLE) err_dbl_mtrx("q");
        
        // We may reach the case where we don't have set n
        if (n < 0) n = MAT_NROWS(q);

        if (q->nrows != n || q->ncols != 1) {
            PyErr_SetString(
                PyExc_ValueError,
                "incompatible dimensions: nrows(q) != n or ncols(q) != 1");
            return NULL;
        }
    }

    if ((PyObject *)P == Py_None) P = NULL;
    if (P) {
        if (!(SpMatrix_Check(P) && SP_ID(P) == DOUBLE)) {
            PyErr_SetString(PyExc_ValueError, "P must be a sparse 'd' matrix");
            return NULL;
        }
        // We may reach the case where we don't have set n
        if (n < 0) n = SP_NROWS(P);

        if (SP_NCOLS(P) != n || SP_NROWS(P) != n) {
            PyErr_SetString(PyExc_ValueError,
                            "P must be square matrix of n x n");
            return NULL;
        }
    }


    if ((PyObject *) x_l == Py_None) x_l = NULL;
    if ((PyObject *)x_u == Py_None) x_u = NULL;
    
    if (!q && !P){
        PyErr_SetString(PyExc_ValueError, "At least q or P must be defined");
        return NULL;
    }
    if (!G && !A && !x_l && !x_u) {
        PyErr_SetString(PyExc_ValueError, "Unbounded problem");
        return NULL;
    }

    if (x_l){
        if (!Matrix_Check(x_l) || x_l->id != DOUBLE) err_dbl_mtrx("x_l");

        if (x_l->nrows != n || x_l->ncols != 1) {
            PyErr_SetString(
                PyExc_ValueError,
                "incompatible dimensions: nrows(x_l) != n or ncols(x_l) != 1");
            return NULL;
        }
        Py_INCREF(x_l);
    } else{
        if (!(x_l = Matrix_New(n, 1, DOUBLE))){
            PyErr_NoMemory();
            return NULL;
        }
        for (i = 0; i < n; i++)
            MAT_BUFD(x_l)[i] = -GRB_INFINITY;
    }
    if (x_u) {
        if (!Matrix_Check(x_u) || x_u->id != DOUBLE) err_dbl_mtrx("x_u");

        if (x_u->nrows != n || x_u->ncols != 1) {
            PyErr_SetString(
                PyExc_ValueError,
                "incompatible dimensions: nrows(x_u) != n or ncols(x_u) != 1");
            return NULL;
        }
    } else {
        if (!(x_u = Matrix_New(n, 1, DOUBLE))) {
            PyErr_NoMemory();
            return NULL;
        }
        for (i = 0; i < n; i++) MAT_BUFD(x_u)[i] = GRB_INFINITY;
    }

    error = solve_problem(P, q, G_l, G, G_u, A, b, x_l, x_u, opts, &res);
    Py_DECREF(x_l);
    Py_DECREF(x_u);

    if (error == 100) {
        PyErr_NoMemory();
        return NULL;
    }
    else if (error) {
        return NULL;
    }
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
    matrix *q, *h, *b = NULL, *x, *z, *y, *z1, *x_l, *x_u;
    spmatrix *P = NULL, *G, *A = NULL;
    PyObject *opts = NULL, *res = NULL, *res_grb = NULL, *status;
    int_t i, m, n, p = 0;
    int error;

    char *kwlist[] = {"q", "G", "h", "A", "b", "P", "options", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOO|OOOO", kwlist, &q, &G,
                                     &h, &A, &b, &P, &opts))
        return NULL;

    if (!(SpMatrix_Check(G) && SP_ID(G) == DOUBLE)) {
        PyErr_SetString(PyExc_TypeError, "G must be a sparse 'd' matrix");
        return NULL;
    }
    if ((m = SP_NROWS(G)) <= 0) err_p_int("m");
    if ((n = SP_NCOLS(G)) <= 0) err_p_int("n");

    if (!Matrix_Check(h) || h->id != DOUBLE) err_dbl_mtrx("h");
    if (h->nrows != m || h->ncols != 1) {
        PyErr_SetString(
            PyExc_ValueError,
            "incompatible dimensions: nrows(h) != m or ncols(h) != 1");
        return NULL;
    }

    if (!Matrix_Check(q) || q->id != DOUBLE) err_dbl_mtrx("q");
    if (q->nrows != n || q->ncols != 1) {
        PyErr_SetString(
            PyExc_ValueError,
            "incompatible dimensions: nrows(q) != n or ncols(h) != 1");
        return NULL;
    }

    if ((PyObject *)A == Py_None ) A = NULL;
    if (A) {
        if (!(SpMatrix_Check(A) && SP_ID(A) == DOUBLE)) {
            PyErr_SetString(PyExc_ValueError, "A must be a sparse 'd' matrix");
            return NULL;
        }
        if ((p = SP_NROWS(A)) < 0) err_p_int("p");
        if (SP_NCOLS(A) != n) {
            PyErr_SetString(
                PyExc_ValueError,
                "incompatible dimensions: ncols(A) != n");
            return NULL;
        }
    }

    if ((PyObject *)b == Py_None ) b = NULL;
    if (b){
    
        if (!Matrix_Check(b) || b->id != DOUBLE) err_dbl_mtrx("b");
        if (b->nrows != p || b->ncols != 1) {
            PyErr_SetString(
                PyExc_ValueError,
                "incompatible dimensions: nrows(b) != p or ncols(b) != 1");
            return NULL;
        }
    }

    if ((PyObject *)P == Py_None) P = NULL; 
    if (P){
        if (!(SpMatrix_Check(P) && SP_ID(P) == DOUBLE)) {
            PyErr_SetString(PyExc_ValueError, "P must be a sparse 'd' matrix");
            return NULL;
        }

        if (SP_NCOLS(P) != n || SP_NROWS(P) != n) {
            PyErr_SetString(PyExc_ValueError,
                            "P must be square matrix of n x n");
            return NULL;
        }

    }

    // CVXOPT forlumation does not consider the variable limits. We create
    // ficticius bounds. To use variable bounds use gurobi.solve

    x_l = Matrix_New(n, 1, DOUBLE);
    x_u = Matrix_New(n, 1, DOUBLE);
    if (!x_l || !x_u) {
        PyErr_NoMemory();
        return NULL;
    }
    for (i = 0; i < n; i++) MAT_BUFD(x_l)[i] = -(MAT_BUFD(x_u)[i] = GRB_INFINITY);

    error = solve_problem(P, q, NULL, G, h, A, b, x_l, x_u, opts, &res_grb);
    Py_DECREF(x_l);
    Py_DECREF(x_u);

    if (error == 100){
        PyErr_NoMemory();
        return NULL;
    }
    else if (error)
        return NULL;

    if (!(res = PyTuple_New(4))) {
        PyErr_NoMemory();
        return NULL;
    } 
    if (!(y = (matrix *)Matrix_New(p, 1, DOUBLE))) {
        PyErr_NoMemory();
        return NULL;
    }

    status = PyTuple_GET_ITEM(res_grb, 0);
    x = (matrix *)PyTuple_GET_ITEM(res_grb, 1);
    z = (matrix *)PyTuple_GET_ITEM(res_grb, 2);
    Py_INCREF(status);
    Py_INCREF(x);
    Py_INCREF(z);

    Py_DECREF(res_grb);


    PyTuple_SET_ITEM(res, 0, status);
    PyTuple_SET_ITEM(res, 1, (PyObject *)x);

    if (p > 0) {

        z1 = (matrix *)Matrix_New(m, 1, DOUBLE);
        memcpy(MAT_BUFD(z1), MAT_BUFD(z), m * sizeof(double));
        memcpy(MAT_BUFD(y), &MAT_BUFD(z)[m], p * sizeof(double));
        PyTuple_SET_ITEM(res, 2, (PyObject *)z1);
        Py_DECREF(z);

    } else {
        PyTuple_SET_ITEM(res, 2, (PyObject *)z);
    }

    PyTuple_SET_ITEM(res, 3, (PyObject *)y);

    return res;
}

static PyMethodDef gurobi_functions[] = {
    {"qp", (PyCFunction)qp, METH_VARARGS | METH_KEYWORDS, doc_qp},
    {"solve", (PyCFunction)solve, METH_VARARGS | METH_KEYWORDS, doc_solve},
    {NULL} /* Sentinel */
};

static PyModuleDef gurobi_module_def = {PyModuleDef_HEAD_INIT,
                                        "gurobi",
                                        gurobi__doc__,
                                        -1,
                                        gurobi_functions,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL};

PyMODINIT_FUNC PyInit_gurobi(void) {
    GRBenv *env = NULL;
    int_t error = 0;

    if (!(gurobi_module = PyModule_Create(&gurobi_module_def))) return NULL;
    PyModule_AddObject(gurobi_module, "options", PyDict_New());
    if (import_kvxopt() < 0) return NULL;

    // Check for Gurobi license
    error = GRBemptyenv(&env);
    if (error) {
        GRBfreeenv(env);
        PY_ERR(PyExc_ImportError, "Gurobi license not found.")
    }
    error = GRBstartenv(env);
    GRBfreeenv(env);
    if (error) {
        PY_ERR(PyExc_ImportError, "Gurobi license not found.")
    }

    return gurobi_module;
}
