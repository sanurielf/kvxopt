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
#include "cvxopt.h"
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

static PyObject *solve_problem(spmatrix *P, matrix *_q, spmatrix *G, matrix *h,
                               spmatrix *A, matrix *b) {
  /* Solve a QP/LP problem in the following form:

  minimize    (1/2)*x'*P*x + q'*x
  subject to  G*x <= h
              A*x = b
  */

  PyObject *res;
  matrix *x, *z;
  int_t error = 0, n, i, j, p, q, nz, n_const, n_nnz;
  double *l_x = NULL, *u_x = NULL, *cval = NULL;
  int *cbeg = NULL, *cind = NULL, *buf = NULL, status;
  char *c_sense = NULL;

  GRBenv *env = NULL;
  GRBmodel *mod = NULL;

  if ((error = GRBemptyenv(&env))) goto CLEAN;

  if ((error = GRBsetstrparam(env, "LogFile", "gurobi_kvxopt.log"))) goto CLEAN;

  if ((error = GRBstartenv(env))) goto CLEAN;

  n = _q->nrows;

  // Allocate variables
  l_x = (double *)malloc(n * sizeof(double));
  u_x = (double *)malloc(n * sizeof(double));
  if (!l_x || !u_x) {
    error = 100;
    goto CLEAN;
  }
  for (i = 0; i < n; i++) l_x[i] = -(u_x[i] = GRB_INFINITY);

  if ((error = GRBnewmodel(env, &mod, "KVXOPT", n, MAT_BUFD(_q), l_x, u_x, NULL,
                           NULL)))
    goto CLEAN;

  if ((error = GRBupdatemodel(mod))) goto CLEAN;
  if ((error = GRBsetintattr(mod, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE)))
    goto CLEAN;

  // Add constraints
  // We add first the inequality constraints:
  // G*x <= h
  // but alloc the workspace to the maximum number of constraints
  // to be added
  n_const = SP_NROWS(G) > SP_NROWS(A) ? SP_NROWS(G) : SP_NROWS(A);
  n_nnz = SP_NNZ(G) > SP_NNZ(A) ? SP_NNZ(G) : SP_NNZ(A);

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

  if ((error = GRBaddconstrs(mod, SP_NROWS(G), SP_NNZ(G), cbeg, cind, cval,
                             c_sense, MAT_BUFD(h), NULL)))
    goto CLEAN;

  if ((error = GRBupdatemodel(mod))) goto CLEAN;

  if (A) {
    for (i = 0; i < n_const; i++) {
      buf[i] = 0;
      cbeg[i] = 0;
    }
    for (i = 0; i < n_nnz; i++) {
      cind[i] = 0;
      cval[i] = 0;
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

    if ((error = GRBaddconstrs(mod, SP_NROWS(A), SP_NNZ(A), cbeg, cind, cval,
                               c_sense, MAT_BUFD(b), NULL)))
      goto CLEAN;

    if ((error = GRBupdatemodel(mod))) goto CLEAN;
  }

  GRBwrite(mod, "output.lp");

  if ((error = GRBoptimize(mod))) goto CLEAN;

  if ((error = GRBgetintattr(mod, "Status", &status))) goto CLEAN;

  n_const = SP_NROWS(G) + SP_NROWS(A);
  x = Matrix_New(n, 1, DOUBLE);
  z = Matrix_New(n_const, 1, DOUBLE);

  if (!x || !z) {
    Py_XDECREF(x);
    Py_XDECREF(z);
    error = 100;
    goto CLEAN;
  }

  if (!(res = PyTuple_New(3))) {
    Py_XDECREF(x);
    Py_XDECREF(z);
    error = 100;
    goto CLEAN;
  }

  PyTuple_SET_ITEM(res, 0,
                   (PyObject *)PYSTRING_FROMSTRING(STATUS_CODES[status]));

  GRBgetdblattrarray(mod, "X", 0, n, MAT_BUFD(x));
  PyTuple_SET_ITEM(res, 1, (PyObject *)x);

  GRBgetdblattrarray(mod, "Pi", 0, n_const, MAT_BUFD(z));
  for (i = 0; i < n_const; i++)
    MAT_BUFD(z)[i] *= -1;
  PyTuple_SET_ITEM(res, 2, (PyObject *)z);

CLEAN:
  GRBfreemodel(mod);
  GRBfreeenv(env);
  free(l_x);
  free(u_x);

  free(cbeg);
  free(c_sense);
  free(cind);
  free(cval);
  free(buf);

  if (error) {
    printf("Error %ld\n", error);
    if (error > 1000) printf("%s\n", GRBgeterrormsg(env));
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
  matrix *q, *h, *b = NULL, *x, *z, *y, *z1;
  spmatrix *P = NULL, *G, *A = NULL;
  PyObject *opts = NULL, *res = NULL, *res_grb = NULL, *status;
  int_t m, n, p;

  char *kwlist[] = {"q", "G", "h", "A", "b", "P", "options", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOO|OOOO", kwlist, &q, &G, &h,
                                   &A, &b, &P, &opts))
    return NULL;

  if (!(SpMatrix_Check(G) && SP_ID(G) == DOUBLE)) {
    PyErr_SetString(PyExc_TypeError, "G must be a sparse 'd' matrix");
    return NULL;
  }
  if ((m = SP_NROWS(G)) <= 0) err_p_int("m");
  if ((n = SP_NCOLS(G)) <= 0) err_p_int("n");

  if (!Matrix_Check(h) || h->id != DOUBLE) err_dbl_mtrx("h");
  if (h->nrows != m || h->ncols != 1) {
    PyErr_SetString(PyExc_ValueError, "incompatible dimensions");
    return NULL;
  }

  if (!Matrix_Check(q) || q->id != DOUBLE) err_dbl_mtrx("q");
  if (q->nrows != n || q->ncols != 1) {
    PyErr_SetString(PyExc_ValueError, "incompatible dimensions");
    return NULL;
  }

  if ((PyObject *)A == Py_None) A = NULL;
  if (A) {
    if (!(SpMatrix_Check(A) && SP_ID(A) == DOUBLE)) {
      PyErr_SetString(PyExc_ValueError, "A must be a sparse 'd' matrix");
      return NULL;
    }
    if ((p = SP_NROWS(A)) < 0) err_p_int("p");
    if (SP_NCOLS(A) != n) {
      PyErr_SetString(PyExc_ValueError, "incompatible dimensions");
      return NULL;
    }
  } else
    p = 0;

  if ((PyObject *)b == Py_None) b = NULL;
  if (b && (!Matrix_Check(b) || b->id != DOUBLE)) err_dbl_mtrx("b");
  if ((b && (b->nrows != p || b->ncols != 1)) || (!b && p != 0)) {
    PyErr_SetString(PyExc_ValueError, "incompatible dimensions");
    return NULL;
  }

  if ((PyObject *)P == Py_None) P = NULL;
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
  }

  res_grb = solve_problem(P, q, G, h, A, b);

  if (!res_grb) return PyErr_NoMemory();

  if (!(res = PyTuple_New(4))) return PyErr_NoMemory();

  status = PyTuple_GET_ITEM(res_grb, 0);
  x = (matrix *)PyTuple_GET_ITEM(res_grb, 1);
  z = (matrix *)PyTuple_GET_ITEM(res_grb, 2);
  Py_INCREF(status);
  Py_INCREF(x);
  Py_INCREF(z);

  Py_DECREF(res_grb);

  y = (matrix *)Matrix_New(p, 1, DOUBLE);

  PyTuple_SET_ITEM(res, 0, status);
  PyTuple_SET_ITEM(res, 1, (PyObject *)x);

  if (p > 0) {
    z1 = (matrix *)Matrix_New(m, 1, DOUBLE);
    memcpy(MAT_BUFD(z1), MAT_BUFD(z), m * sizeof(double));
    memcpy(MAT_BUFD(y), &MAT_BUFD(z)[m], p * sizeof(double));
    Py_DECREF(z);
    PyTuple_SET_ITEM(res, 2, (PyObject *)z1);
  } else {
    PyTuple_SET_ITEM(res, 2, (PyObject *)z);
  }

  PyTuple_SET_ITEM(res, 3, (PyObject *)y);

  return res;
}

static PyMethodDef gurobi_functions[] = {
    {"qp", (PyCFunction)qp, METH_VARARGS | METH_KEYWORDS, doc_qp},
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
  if (error){
    PY_ERR(PyExc_ImportError, "Gurobi license not found.")
  }

  return gurobi_module;
}
