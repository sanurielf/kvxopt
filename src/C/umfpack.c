/*
 * Copyright 2012-2023 M. Andersen and L. Vandenberghe.
 * Copyright 2010-2011 L. Vandenberghe.
 * Copyright 2004-2009 J. Dahl and L. Vandenberghe.
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

#include "cvxopt.h"
#include "umfpack.h"
#include "misc.h"

#if (SIZEOF_INT < SIZEOF_SIZE_T)
#define UMFD(name) umfpack_dl_ ## name
#define UMFZ(name) umfpack_zl_ ## name
#else
#define UMFD(name) umfpack_di_ ## name
#define UMFZ(name) umfpack_zi_ ## name
#endif

#ifndef _MSC_VER
const int E_SIZE[] = {sizeof(int_t), sizeof(double), sizeof(double complex)};
#else
const int E_SIZE[] = {sizeof(int_t), sizeof(double), sizeof(_Dcomplex)};
#endif

static char umfpack_error[40];

const char *descrdFs = "UMFPACK SYM D FACTOR";
const char *descrzFs = "UMFPACK SYM Z FACTOR";

const char *descrdFn = "UMFPACK NUM D FACTOR";
const char *descrzFn = "UMFPACK NUM Z FACTOR";

PyDoc_STRVAR(umfpack__doc__,"Interface to the UMFPACK library.\n\n"
    "Routines for symbolic and numeric LU factorization of sparse\n"
    "matrices and for solving sparse sets of linear equations.\n"
    "The default control settings of UMPFACK are used.\n\n"
    "See also www.suitesparse.com.");

static void free_umfpack_d_symbolic(void *F)
{
    void *Fptr = PyCapsule_GetPointer(F, PyCapsule_GetName(F));
    UMFD(free_symbolic)(&Fptr);
}

static void free_umfpack_z_symbolic(void *F)
{
    void *Fptr = PyCapsule_GetPointer(F, PyCapsule_GetName(F));
    UMFZ(free_symbolic)(&Fptr);
}

static void free_umfpack_d_numeric(void *F)
{
    void *Fptr = PyCapsule_GetPointer(F, PyCapsule_GetName(F));
    UMFD(free_numeric)(&Fptr);
}

static void free_umfpack_z_numeric(void *F)
{
    void *Fptr = PyCapsule_GetPointer(F, PyCapsule_GetName(F));
    UMFZ(free_numeric)(&Fptr);
}

static char doc_linsolve[] =
    "Solves a sparse set of linear equations.\n\n"
    "linsolve(A, B, trans='N', nrhs=B.size[1], ldB=max(1,B.size[0]),\n"
    "         offsetB=0)\n\n"
    "PURPOSE\n"
    "If trans is 'N', solves A*X = B.\n"
    "If trans is 'T', solves A^T*X = B.\n"
    "If trans is 'C', solves A^H*X = B.\n"
    "A is a sparse n by n matrix, and B is n by nrhs.\n"
    "On exit B is replaced by the solution.  A is not modified.\n\n"
    "ARGUMENTS\n"
    "A         square sparse matrix\n\n"
    "B         dense matrix of the same type as A, stored following \n"
    "          the BLAS conventions\n\n"
    "trans     'N', 'T' or 'C'\n\n"
    "nrhs      integer.  If negative, the default value is used.\n\n"
    "ldB       nonnegative integer.  ldB >= max(1,n).  If zero, the\n"
    "          default value is used.\n\n"
    "offsetB   nonnegative integer";

static PyObject* linsolve(PyObject *self, PyObject *args,
    PyObject *kwrds)
{
    spmatrix *A;
    matrix *B;
    int trans_ = 'N';
    char trans='N';
    double info[UMFPACK_INFO];
    int oB=0, n, nrhs=-1, ldB=0, k;
    void *symbolic, *numeric, *x;
    char *kwlist[] = {"A", "B", "trans", "nrhs", "ldB", "offsetB",
        NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwrds, "OO|Ciii", kwlist,
        &A, &B, &trans_, &nrhs, &ldB, &oB)) return NULL;
    trans = (char) trans_;

    if (!SpMatrix_Check(A) || SP_NROWS(A) != SP_NCOLS(A))
        PY_ERR_TYPE("A must be a square sparse matrix");
    n = SP_NROWS(A);
    if (!Matrix_Check(B) || MAT_ID(B) != SP_ID(A))
        PY_ERR_TYPE("B must a dense matrix of the same numeric type "
            "as A");

    if (nrhs < 0) nrhs = B->ncols;
    if (n == 0 || nrhs == 0) return Py_BuildValue("i", 0);
    if (ldB == 0) ldB = MAX(1,B->nrows);
    if (ldB < MAX(1,n)) err_ld("ldB");
    if (oB < 0) err_nn_int("offsetB");
    if (oB + (nrhs-1)*ldB + n > MAT_LGT(B)) err_buf_len("B");

    if (trans != 'N' && trans != 'T' && trans != 'C')
        err_char("trans", "'N', 'T', 'C'");

    if (SP_ID(A) == DOUBLE)
        UMFD(symbolic)(n, n, SP_COL(A), SP_ROW(A), SP_VAL(A), &symbolic,
            NULL, info);
    else
        UMFZ(symbolic)(n, n, SP_COL(A), SP_ROW(A), SP_VAL(A), NULL,
            &symbolic, NULL, info);

    if (info[UMFPACK_STATUS] != UMFPACK_OK){
        if (SP_ID(A) == DOUBLE)
            UMFD(free_symbolic)(&symbolic);
        else
            UMFZ(free_symbolic)(&symbolic);
        if (info[UMFPACK_STATUS] == UMFPACK_ERROR_out_of_memory)
            return PyErr_NoMemory();
        else {
            snprintf(umfpack_error, 40,"%s %i","UMFPACK ERROR",
                (int) info[UMFPACK_STATUS]);
            PyErr_SetString(PyExc_ValueError, umfpack_error);
            return NULL;
        }
    }

    if (SP_ID(A) == DOUBLE) {
        UMFD(numeric)(SP_COL(A), SP_ROW(A), SP_VAL(A), symbolic,
            &numeric, NULL, info);
        UMFD(free_symbolic)(&symbolic);
    } else {
        UMFZ(numeric)(SP_COL(A), SP_ROW(A), SP_VAL(A), NULL, symbolic,
            &numeric, NULL, info);
        UMFZ(free_symbolic)(&symbolic);
    }
    if (info[UMFPACK_STATUS] != UMFPACK_OK){
        if (SP_ID(A) == DOUBLE)
            UMFD(free_numeric)(&numeric);
        else
            UMFZ(free_numeric)(&numeric);
        if (info[UMFPACK_STATUS] == UMFPACK_ERROR_out_of_memory)
            return PyErr_NoMemory();
        else {
            if (info[UMFPACK_STATUS] == UMFPACK_WARNING_singular_matrix)
                PyErr_SetString(PyExc_ArithmeticError, "singular "
                    "matrix");
            else {
                snprintf(umfpack_error, 40,"%s %i","UMFPACK ERROR",
                    (int) info[UMFPACK_STATUS]);
                PyErr_SetString(PyExc_ValueError, umfpack_error);
            }
            return NULL;
        }
    }

    if (!(x = malloc(n*E_SIZE[SP_ID(A)]))) {
        if (SP_ID(A) == DOUBLE)
            UMFD(free_numeric)(&numeric);
        else
            UMFZ(free_numeric)(&numeric);
        return PyErr_NoMemory();
    }
    for (k=0; k<nrhs; k++){
        if (SP_ID(A) == DOUBLE)
            UMFD(solve)(trans == 'N' ? UMFPACK_A: UMFPACK_Aat,
                SP_COL(A), SP_ROW(A), SP_VAL(A), x,
                MAT_BUFD(B) + k*ldB + oB, numeric, NULL, info);
        else
            UMFZ(solve)(trans == 'N' ? UMFPACK_A: trans == 'C' ?
                UMFPACK_At : UMFPACK_Aat, SP_COL(A), SP_ROW(A),
                SP_VAL(A), NULL, x, NULL,
                (double *)(MAT_BUFZ(B) + k*ldB + oB), NULL, numeric,
                NULL, info);
        if (info[UMFPACK_STATUS] == UMFPACK_OK)
	  memcpy((unsigned char*)B->buffer + (k*ldB + oB)*E_SIZE[SP_ID(A)], x,
                n*E_SIZE[SP_ID(A)]);
        else
	    break;
    }
    free(x);
    if (SP_ID(A) == DOUBLE)
        UMFD(free_numeric)(&numeric);
    else
        UMFZ(free_numeric)(&numeric);

    if (info[UMFPACK_STATUS] != UMFPACK_OK){
        if (info[UMFPACK_STATUS] == UMFPACK_ERROR_out_of_memory)
            return PyErr_NoMemory();
        else {
            if (info[UMFPACK_STATUS] == UMFPACK_WARNING_singular_matrix)
                PyErr_SetString(PyExc_ArithmeticError, "singular "
                    "matrix");
            else {
                snprintf(umfpack_error, 40,"%s %i","UMFPACK ERROR",
                    (int) info[UMFPACK_STATUS]);
                PyErr_SetString(PyExc_ValueError, umfpack_error);
            }
        return NULL;
        }
    }
    return Py_BuildValue("");
}


static char doc_symbolic[] =
    "Symbolic LU factorization of a sparse matrix.\n\n"
    "F = symbolic(A)\n\n"
    "ARGUMENTS\n"
    "A         sparse matrix with at least one row and at least one\n"
    "          column.  A may be rectangular.\n\n"
    "F         the symbolic factorization as an opaque C object";

static PyObject* symbolic(PyObject *self, PyObject *args)
{
    spmatrix *A;
    double info[UMFPACK_INFO];
    void *symbolic;

    if (!PyArg_ParseTuple(args, "O", &A)) return NULL;

    if (!SpMatrix_Check(A)) PY_ERR_TYPE("A must be a sparse matrix");
    if (SP_NCOLS(A) == 0 || SP_NROWS(A) == 0) {
        PyErr_SetString(PyExc_ValueError, "A must have at least one "
            "row and column");
        return NULL;
    }

    switch (SP_ID(A)){
        case DOUBLE:
            UMFD(symbolic)(SP_NROWS(A), SP_NCOLS(A), SP_COL(A),
                SP_ROW(A), SP_VAL(A), &symbolic, NULL, info);
            if (info[UMFPACK_STATUS] == UMFPACK_OK)
                return (PyObject *) PyCapsule_New( (void *) symbolic,
                    descrdFs,
                    (PyCapsule_Destructor) &free_umfpack_d_symbolic);

            else
                UMFD(free_symbolic)(&symbolic);
            break;

        case COMPLEX:
            UMFZ(symbolic)(SP_NROWS(A), SP_NCOLS(A), SP_COL(A),
                SP_ROW(A), SP_VAL(A), NULL, &symbolic, NULL, info);
            if (info[UMFPACK_STATUS] == UMFPACK_OK)
                return (PyObject *) PyCapsule_New(
                    (void *) symbolic, descrzFs,
                    (PyCapsule_Destructor) &free_umfpack_z_symbolic);

            else
                UMFZ(free_symbolic)(&symbolic);
            break;
    }

    if (info[UMFPACK_STATUS] == UMFPACK_ERROR_out_of_memory)
        return PyErr_NoMemory();
    else {
        snprintf(umfpack_error, 40,"%s %i","UMFPACK ERROR",
            (int) info[UMFPACK_STATUS]);
        PyErr_SetString(PyExc_ValueError, umfpack_error);
        return NULL;
    }
}


static char doc_numeric[] =
    "Numeric LU factorization of a sparse matrix, given a symbolic\n"
    "factorization computed by umfpack.symbolic.  Raises an\n"
    "ArithmeticError if A is singular.\n\n"
    "Fn = numeric(A, Fs)\n\n"
    "ARGUMENTS\n"
    "A         sparse matrix; may be rectangular\n\n"
    "Fs        symbolic factorization of A, or a matrix with the same\n"
    "          sparsity pattern, dimensions, and typecode as A, \n"
    "          created by umfpack.symbolic\n\n"
    "Fn        the numeric factorization, as an opaque C object";

static PyObject* numeric(PyObject *self, PyObject *args)
{
    spmatrix *A;
    PyObject *Fs;
    double info[UMFPACK_INFO];
    void *numeric;
    void *Fsptr;

    if (!PyArg_ParseTuple(args, "OO", &A, &Fs)) return NULL;

    if (!SpMatrix_Check(A)) PY_ERR_TYPE("A must be a sparse matrix");

    if (!PyCapsule_CheckExact(Fs)) err_CO("Fs");


    switch (SP_ID(A)) {
	case DOUBLE:
            TypeCheck_Capsule(Fs, descrdFs, "Fs is not the UMFPACK symbolic "
                "factor of a 'd' matrix");
            if (!(Fsptr = (void *) PyCapsule_GetPointer(Fs, descrdFs)))
                err_CO("Fs");
            UMFD(numeric)(SP_COL(A), SP_ROW(A), SP_VAL(A), Fsptr, &numeric,
                NULL, info);

            if (info[UMFPACK_STATUS] == UMFPACK_OK)
                return (PyObject *) PyCapsule_New(
                    (void *) numeric, descrdFn,
                    (PyCapsule_Destructor) &free_umfpack_d_numeric);

            else
                UMFD(free_numeric)(&numeric);
	    break;

        case COMPLEX:
            TypeCheck_Capsule(Fs, descrzFs, "Fs is not the UMFPACK symbolic "
                "factor of a 'z' matrix");
            if (!(Fsptr = (void *) PyCapsule_GetPointer(Fs, descrzFs)))
                err_CO("Fs");
            UMFZ(numeric)(SP_COL(A), SP_ROW(A), SP_VAL(A), NULL, Fsptr,
                &numeric, NULL, info);

            if (info[UMFPACK_STATUS] == UMFPACK_OK)
                return (PyObject *) PyCapsule_New(
                    (void *) numeric, descrzFn,
                    (PyCapsule_Destructor) &free_umfpack_z_numeric);

	    else
                 UMFZ(free_numeric)(&numeric);
	    break;
    }

    if (info[UMFPACK_STATUS] == UMFPACK_ERROR_out_of_memory)
        return PyErr_NoMemory();
    else {
        if (info[UMFPACK_STATUS] == UMFPACK_WARNING_singular_matrix)
	    PyErr_SetString(PyExc_ArithmeticError, "singular matrix");
        else {
	    snprintf(umfpack_error, 40,"%s %i","UMFPACK ERROR",
	        (int) info[UMFPACK_STATUS]);
	    PyErr_SetString(PyExc_ValueError, umfpack_error);
        }
        return NULL;
    }
}

static char doc_get_numeric[] =
    "This routine copies the LU factors and permutation vectors from the \n"
    "Numeric object into user-accessible arrays.  This routine is not \n"
    "needed to solve a linear system..\n\n"
    "L, U, P, Q, R = get_numeric(A, Fn)\n\n"
    "ARGUMENTS\n"
    "A         sparse matrix; must be square\n\n"
    "Fn        numeric factorization, as returned by umfpack.numeric\n\n";

static PyObject* get_numeric(PyObject *self, PyObject *args, PyObject *kwrds)
{
    spmatrix *A, *L, *U, *P, *Q, *R;
    PyObject *F, *res;
    void *numeric;
    int trans_ = 'N';
    char trans='N';
    int_t i, n_row, n_col, nn, n_inner, lnz, unz, status, ignore1, ignore2, ignore3,
          *Ltp, *Ltj, *Up, *Ui, *Pt, *Qt, *Rp, *Ri, do_recip;
    double *Ltx, *Ux, *Rs;
    char *kwlist[] = {"A", "F", "trans", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwrds, "OO|C", kwlist,
        &A, &F, &trans_)) return NULL;
    trans = (char) trans_;

    if (!SpMatrix_Check(A)) PY_ERR_TYPE("A must be a sparse matrix");


    n_row = SP_NROWS(A);
    n_col = SP_NCOLS(A);
    nn = MAX(n_row, n_col);
    n_inner = MIN(n_row, n_col);

    if (!PyCapsule_CheckExact(F)) err_CO("F");
    if (SP_ID(A) == DOUBLE) {
        TypeCheck_Capsule(F, descrdFn, "F is not the UMFPACK numeric factor "
            "of a 'd' matrix");
    }
    else  {
        TypeCheck_Capsule(F, descrzFn, "F is not the UMFPACK numeric factor "
            "of a 'z' matrix");
    }

    switch (SP_ID(A)) {
        case DOUBLE:
            numeric = (void *) PyCapsule_GetPointer(F, descrdFn);
            status = UMFD(get_lunz)(&lnz, &unz, &ignore1, &ignore2, &ignore3,
                                    numeric);
            break;

        case COMPLEX:
            numeric = (void *) PyCapsule_GetPointer(F, descrzFn);
            status = UMFZ(get_lunz)(&lnz, &unz, &ignore1, &ignore2, &ignore3,
                                    numeric);
            break;

    }

    if (status < 0){
        snprintf(umfpack_error, 40,"%s %i","Extracting LUnz factors failed",
        (int) status);
        PyErr_SetString(PyExc_ValueError, umfpack_error);
        return NULL;
    }
    Ltp = (int_t *) malloc ((n_row+1) * sizeof (int_t));
    Ltj = (int_t *) malloc (lnz * sizeof (int_t));

    if (SP_ID(A) == COMPLEX){
        Ltx = (double *) malloc (2 * lnz * sizeof (double)) ;
    }
    else{
        Ltx = (double *) malloc (lnz * sizeof (double));
    }

    /* Matrix U */
    U = SpMatrix_New(n_inner, n_col, unz, SP_ID(A));
    Up = SP_COL(U);
    Ui = SP_ROW(U);
    Ux = SP_VAL(U);

    /* temporary space for the integer permutation vectors */
    Pt = (int_t *) malloc(n_row * sizeof(int_t));
    Qt = (int_t *) malloc(n_col * sizeof(int_t));

    /* create a diagonal sparse matrix for the scale factors R*/
    R = SpMatrix_New(n_row, n_row, n_row, DOUBLE);
    Rp = SP_COL(R);
    Ri = SP_ROW(R);
    for (i = 0; i < n_row; i++){
        Rp[i] = i;
        Ri[i] = i;
    }
    Rp[n_row] = n_row;
    Rs = SP_VALD(R);

    /* get Lt, U, P, Q, and Rs from the numeric object */
    switch (SP_ID(A)) {
    case DOUBLE:

        status = UMFD(get_numeric)(Ltp, Ltj, Ltx, Up, Ui, Ux,
                      Pt, Qt, NULL, &do_recip, Rs, numeric);
        break;

    case COMPLEX:

        status = UMFZ(get_numeric)(Ltp, Ltj, Ltx, NULL, Up, Ui, Ux, NULL, Pt,
                                   Qt, NULL, NULL, &do_recip, Rs, numeric);
        break;
    }

    if (status < 0){
        free(Ltp);
        free(Ltj);
        free(Ltx);
        free(Pt);
        free(Qt);
        // TODO: Free R, U??
        snprintf(umfpack_error, 40,"%s %i","Extracting LU factors failed",
        (int) status);
        PyErr_SetString(PyExc_ValueError, umfpack_error);
        return NULL;
    }

    /* We do the reciprocal here, instead of passing the flag to user */
    if (!do_recip){
        for (i = 0; i < n_row; i++)
            Rs[i] = 1.0 / Rs[i];
    }

    /* create sparse permutation matrix for P */
    P = SpMatrix_New(n_row, n_row, n_row, DOUBLE);
    for (i = 0; i < n_row; i++){
        SP_COL(P)[i] = i;
        SP_ROW(P)[Pt[i]] = i;
        SP_VALD(P)[i] = 1;
    }
    SP_COL(P)[n_row] = n_row;

    /* create sparse permutation matrix for Q */
    Q = SpMatrix_New(n_col, n_col, n_col, DOUBLE);
    for (i = 0; i < n_col; i++){
        SP_COL(Q)[i] = i;
        SP_ROW(Q)[i] = Qt[i];
        SP_VALD(Q)[i] = 1;
    }
    SP_COL(Q)[n_col] = n_col;

    /* Matrix L */
    L = SpMatrix_New(n_row, n_inner, lnz, SP_ID(A));

    /* convert L from row form to column form */
    switch (SP_ID(A)) {
    case DOUBLE:
        status = UMFD(transpose)(n_inner, n_row, Ltp, Ltj, Ltx, NULL, NULL,
                                 SP_COL(L), SP_ROW(L), SP_VALD(L));
        break;
    case COMPLEX:
        status = UMFZ(transpose)(n_inner, n_row, Ltp, Ltj, Ltx, NULL, NULL,
                                 NULL, SP_COL(L), SP_ROW(L), SP_VALD(L), NULL, 0);
        break;
    }


    /* Free workspace */
    free(Ltp);
    free(Ltj);
    free(Ltx);
    free(Pt);
    free(Qt);

    if (status < 0){
        // TODO: Free L, U, P, Q, R?
        snprintf(umfpack_error, 40,"%s %i","Transposing L failed",
        (int) status);
        PyErr_SetString(PyExc_ValueError, umfpack_error);
        return NULL;
    }

    if (!(res = PyTuple_New(5))) return PyErr_NoMemory();
    PyTuple_SET_ITEM(res, 0, (PyObject *) L);
    PyTuple_SET_ITEM(res, 1, (PyObject *) U);
    PyTuple_SET_ITEM(res, 2, (PyObject *) P);
    PyTuple_SET_ITEM(res, 3, (PyObject *) Q);
    PyTuple_SET_ITEM(res, 4, (PyObject *) R);

    return res;


}

static char doc_solve[] =
    "Solves a factored set of linear equations.\n\n"
    "solve(A, F, B, trans='N', nrhs=B.size[1], ldB=max(1,B.size[0]),\n"
    "      offsetB=0)\n\n"
    "PURPOSE\n"
    "If trans is 'N', solves A*X = B.\n"
    "If trans is 'T', solves A^T*X = B.\n"
    "If trans is 'C', solves A^H*X = B.\n"
    "A is a sparse n by n matrix, and B is n by nrhs.  F is the\n"
    "numeric factorization of A, computed by umfpack.numeric.\n"
    "On exit B is replaced by the solution.  A is not modified.\n\n"
    "ARGUMENTS\n"
    "A         square sparse matrix\n\n"
    "F         numeric factorization, as returned by umfpack.numeric\n"
    "\n"
    "B         dense matrix of the same type as A, stored following \n"
    "          the BLAS conventions\n\n"
    "trans     'N', 'T' or 'C'\n\n"
    "nrhs      integer.  If negative, the default value is used.\n\n"
    "ldB       nonnegative integer.  ldB >= max(1,n).  If zero, the\n"
    "          default value is used.\n\n"
    "offsetB   nonnegative integer";

static PyObject* solve(PyObject *self, PyObject *args, PyObject *kwrds)
{
    spmatrix *A;
    PyObject *F;
    matrix *B;
    int trans_ = 'N';
    char trans='N';
    double *x, info[UMFPACK_INFO];
    int oB=0, n, ldB=0, nrhs=-1, k;
    char *kwlist[] = {"A", "F", "B", "trans", "nrhs", "ldB", "offsetB",
        NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwrds, "OOO|Ciii", kwlist,
        &A, &F, &B, &trans_, &nrhs, &ldB, &oB)) return NULL;
    trans = (char) trans_;

    if (!SpMatrix_Check(A) || SP_NROWS(A) != SP_NCOLS(A))
        PY_ERR_TYPE("A must a square sparse matrix");
    n = SP_NROWS(A);

    if (!PyCapsule_CheckExact(F)) err_CO("F");
    if (SP_ID(A) == DOUBLE) {
        TypeCheck_Capsule(F, descrdFn, "F is not the UMFPACK numeric factor "
            "of a 'd' matrix");
    }
    else  {
        TypeCheck_Capsule(F, descrzFn, "F is not the UMFPACK numeric factor "
            "of a 'z' matrix");
    }


    if (!Matrix_Check(B) || MAT_ID(B) != SP_ID(A))
        PY_ERR_TYPE("B must a dense matrix of the same numeric type "
            "as A");
    if (nrhs < 0) nrhs = B->ncols;
    if (n == 0 || nrhs == 0) return Py_BuildValue("");
    if (ldB == 0) ldB = MAX(1,B->nrows);
    if (ldB < MAX(1,n)) err_ld("ldB");
    if (oB < 0) err_nn_int("offsetB");
    if (oB + (nrhs-1)*ldB + n > MAT_LGT(B)) err_buf_len("B");

    if (trans != 'N' && trans != 'T' && trans != 'C')
        err_char("trans", "'N', 'T', 'C'");

    if (!(x = malloc(n*E_SIZE[SP_ID(A)]))) return PyErr_NoMemory();

    for (k=0; k<nrhs; k++) {
        if (SP_ID(A) == DOUBLE)
            UMFD(solve)(trans == 'N' ? UMFPACK_A : UMFPACK_Aat,
                SP_COL(A), SP_ROW(A), SP_VAL(A), x,
                MAT_BUFD(B) + k*ldB + oB,
                (void *) PyCapsule_GetPointer(F, descrdFn), NULL, info);

        else
            UMFZ(solve)(trans == 'N' ? UMFPACK_A : trans == 'C' ?
                UMFPACK_At : UMFPACK_Aat, SP_COL(A), SP_ROW(A),
                SP_VAL(A), NULL, x, NULL,
                (double *)(MAT_BUFZ(B) + k*ldB + oB), NULL,
                (void *) PyCapsule_GetPointer(F, descrzFn), NULL, info);

        if (info[UMFPACK_STATUS] == UMFPACK_OK)
	    memcpy((unsigned char*)B->buffer + (k*ldB + oB)*E_SIZE[SP_ID(A)], x,
                n*E_SIZE[SP_ID(A)]);
        else
	    break;
    }
    free(x);

    if (info[UMFPACK_STATUS] != UMFPACK_OK){
        if (info[UMFPACK_STATUS] == UMFPACK_ERROR_out_of_memory)
            return PyErr_NoMemory();
        else {
            if (info[UMFPACK_STATUS] == UMFPACK_WARNING_singular_matrix)
                PyErr_SetString(PyExc_ArithmeticError,
                    "singular matrix");
            else {
                snprintf(umfpack_error, 40,"%s %i","UMFPACK ERROR",
                    (int) info[UMFPACK_STATUS]);
                PyErr_SetString(PyExc_ValueError, umfpack_error);
            }
            return NULL;
        }
    }

    return Py_BuildValue("");
}



static char doc_get_det[] =
    "Returns determinant of a UMFPACK symbolic/numeric object\n"
    "d = get_det(A, Fs, Fn)\n\n"
    "PURPOSE\n"
    "A is a real sparse n by n matrix, F and N its symbolic and numeric \n"
    "factorizations respectively.\n"
    "On exit a double/complex is returned with the value of the determinant.\n\n"
    "ARGUMENTS\n"
    "A         square sparse matrix\n\n"
    "Fs        the symbolic factorization as an opaque C object\n\n"
    "Fn        the numeric factorization as an opaque C object\n\n"
    "d         the determinant value of the matrix\n\n";

static PyObject* get_det(PyObject *self, PyObject *args, PyObject *kwrds) {
    spmatrix *A;
    PyObject *Fs, *Fn;

    double dx, dz, info[UMFPACK_INFO];
    int_t status;

    if (!PyArg_ParseTuple(args, "OOO", &A, &Fs, &Fn))
        return NULL;

    if (!SpMatrix_Check(A)) PY_ERR_TYPE("A must be a sparse matrix");
    if (!PyCapsule_CheckExact(Fn)) err_CO("Fn");
    if (!PyCapsule_CheckExact(Fs)) err_CO("Fs");

    if (SP_ID(A) == DOUBLE){
        TypeCheck_Capsule(Fn, descrdFn, "F is not the UMFPACK numeric factor "
                          "of a 'd' matrix");
        status = UMFD(get_determinant)(&dx, NULL,
                              (void *) PyCapsule_GetPointer(Fn, descrdFn),
                              info);
    }
    else{
        TypeCheck_Capsule(Fn, descrzFn, "F is not the UMFPACK numeric factor "
                          "of a 'z' matrix");
        status = UMFZ(get_determinant)(&dx, &dz, NULL,
                              (void *) PyCapsule_GetPointer(Fn, descrzFn),
                              info);
    }

    if (status < 0){
        snprintf(umfpack_error, 40,"%s %i","UMFPACK ERROR",
                (int) info[UMFPACK_STATUS]);
        PyErr_SetString(PyExc_ValueError, umfpack_error);
        return NULL;
    }


    if (SP_ID(A) == DOUBLE)
        return Py_BuildValue("d", dx);
    else
        return PyComplex_FromDoubles(dx, dz);

}

static PyMethodDef umfpack_functions[] = {
    {"linsolve", (PyCFunction) linsolve, METH_VARARGS|METH_KEYWORDS,
        doc_linsolve},
    {"symbolic", (PyCFunction) symbolic, METH_VARARGS, doc_symbolic},
    {"numeric", (PyCFunction) numeric, METH_VARARGS, doc_numeric},
    {"get_numeric", (PyCFunction) get_numeric, METH_VARARGS|METH_KEYWORDS, doc_get_numeric},
    {"solve", (PyCFunction) solve, METH_VARARGS|METH_KEYWORDS, doc_solve},
    {"get_det", (PyCFunction) get_det, METH_VARARGS, doc_get_det},
    {NULL}  /* Sentinel */
};

static PyModuleDef umfpack_module = {
    PyModuleDef_HEAD_INIT,
    "umfpack",
    umfpack__doc__,
    -1,
    umfpack_functions,
    NULL, NULL, NULL, NULL
};

PyMODINIT_FUNC PyInit_umfpack(void)
{
  PyObject *m;
  if (!(m = PyModule_Create(&umfpack_module))) return NULL;
  if (import_kvxopt() < 0) return NULL;
  return m;
}
