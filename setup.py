from setuptools import setup, Extension
from glob import glob
import os, sys

# Modifiy this if BLAS and LAPACK libraries are not in /usr/lib.
BLAS_LIB_DIR = '/usr/lib'

# Default names of BLAS and LAPACK libraries
BLAS_LIB = ['blas']
LAPACK_LIB = ['lapack']
FFTW_LIB = ['fftw3']
BLAS_EXTRA_LINK_ARGS = []

# Set environment variable BLAS_NOUNDERSCORES=1 if your BLAS/LAPACK do
# not use trailing underscores
BLAS_NOUNDERSCORES = False

# Set to 1 if you are using the random number generators in the GNU
# Scientific Library.
BUILD_GSL = 0

# Directory containing libgsl (used only when BUILD_GSL = 1).
GSL_LIB_DIR = '/usr/lib'

# Directory containing the GSL header files (used only when BUILD_GSL = 1).
GSL_INC_DIR = '/usr/include/gsl'

# Set to 1 if you are installing the fftw module.
BUILD_FFTW = 0

# Directory containing libfftw3 (used only when BUILD_FFTW = 1).
FFTW_LIB_DIR = '/usr/lib'

# Directory containing fftw.h (used only when BUILD_FFTW = 1).
FFTW_INC_DIR = '/usr/include'

# Set to 1 if you are installing the glpk module.
BUILD_GLPK = 0

# Directory containing libglpk (used only when BUILD_GLPK = 1).
GLPK_LIB_DIR = '/usr/lib'

# Directory containing glpk.h (used only when BUILD_GLPK = 1).
GLPK_INC_DIR = '/usr/include'

# Set to 1 if you are installing the osqp module.
BUILD_OSQP = 0

# Directory containing libglpk (used only when BUILD_OSQP = 1).
OSQP_LIB_DIR = '/usr/local/lib'

# Directory containing glpk.h (used only when BUILD_OSQP = 1).
OSQP_INC_DIR = '/usr/local/include/osqp'

# Set to 1 if you want to compile the Gurobi C extension for DCOPF
BUILD_GRB = 0

# Directory containing libgurobi (used only when BUILD_GRB = 1).
GRB_LIB_DIR = '/Library/gurobi911/mac64/lib'

# Directory containing the Gurobi header files (used only when BUILD_GRB = 1).
GRB_INC_DIR = '/Library/gurobi911/mac64/include'

# Name of the Gurobi lib
GRB_LIB = 'gurobi91'

# Set to 1 if you are installing the DSDP module.
BUILD_DSDP = 0

# Directory containing libdsdp (used only when BUILD_DSDP = 1).
DSDP_LIB_DIR = '/usr/lib'

# Directory containing dsdp5.h (used only when BUILD_DSDP = 1).
DSDP_INC_DIR = '/usr/include/dsdp'

# Guess SUITESPARSE_LIB_DIR and SUITESPARSE_INC_DIR
if sys.platform.startswith("darwin"):
    # macOS
    SUITESPARSE_LIB_DIR = '/usr/local/lib'
    SUITESPARSE_INC_DIR = '/usr/local/include'
else:
    if glob("/usr/lib/x86_64-linux-gnu/libsuitesparse*"):
        # Ubuntu/Debian
        SUITESPARSE_LIB_DIR = "/usr/lib/x86_64-linux-gnu"
        SUITESPARSE_INC_DIR = "/usr/include/suitesparse"
    elif glob("/usr/lib64/libsuitesparse*"):
        # CentOS/Fedora/RedHat
        SUITESPARSE_LIB_DIR = "/usr/lib64"
        SUITESPARSE_INC_DIR = "/usr/include/suitesparse"
    else:
        # Default
        SUITESPARSE_LIB_DIR = '/usr/lib'
        SUITESPARSE_INC_DIR = '/usr/include'

if sys.platform.startswith("win"):
    GSL_MACROS = [('GSL_DLL',''),('WIN32','')]
    FFTW_MACROS = [('FFTW_DLL',''),('FFTW_NO_Complex','')]
else:
    GSL_MACROS = []
    FFTW_MACROS = []

# Directory containing SuiteSparse source
SUITESPARSE_SRC_DIR = ''

# For SuiteSparse Versions before to 4.0.0 SuiteSparse_config does not exist
# We can avoid the search and link with this flag
SUITESPARSE_CONFIG = 1

# Set to 1 if compiling with MSVC 14 or later
MSVC=0

# No modifications should be needed below this line.

BLAS_NOUNDERSCORES = int(os.environ.get("KVXOPT_BLAS_NOUNDERSCORES",BLAS_NOUNDERSCORES)) == True
BLAS_LIB = os.environ.get("KVXOPT_BLAS_LIB",BLAS_LIB)
LAPACK_LIB = os.environ.get("KVXOPT_LAPACK_LIB",LAPACK_LIB)
BLAS_LIB_DIR = os.environ.get("KVXOPT_BLAS_LIB_DIR",BLAS_LIB_DIR)
BLAS_EXTRA_LINK_ARGS = os.environ.get("KVXOPT_BLAS_EXTRA_LINK_ARGS",BLAS_EXTRA_LINK_ARGS)
if type(BLAS_LIB) is str: BLAS_LIB = BLAS_LIB.strip().split(';')
if type(LAPACK_LIB) is str: LAPACK_LIB = LAPACK_LIB.strip().split(';')
if type(BLAS_EXTRA_LINK_ARGS) is str: BLAS_EXTRA_LINK_ARGS = BLAS_EXTRA_LINK_ARGS.strip().split(';')
BUILD_GSL = int(os.environ.get("KVXOPT_BUILD_GSL",BUILD_GSL))
GSL_LIB_DIR = os.environ.get("KVXOPT_GSL_LIB_DIR",GSL_LIB_DIR)
GSL_INC_DIR = os.environ.get("KVXOPT_GSL_INC_DIR",GSL_INC_DIR)
BUILD_FFTW = int(os.environ.get("KVXOPT_BUILD_FFTW",BUILD_FFTW))
FFTW_LIB = os.environ.get("KVXOPT_FFTW_LIB",FFTW_LIB)
FFTW_LIB_DIR = os.environ.get("KVXOPT_FFTW_LIB_DIR",FFTW_LIB_DIR)
FFTW_INC_DIR = os.environ.get("KVXOPT_FFTW_INC_DIR",FFTW_INC_DIR)
if type(FFTW_LIB) is str: FFTW_LIB = FFTW_LIB.strip().split(';')
BUILD_GLPK = int(os.environ.get("KVXOPT_BUILD_GLPK",BUILD_GLPK))
GLPK_LIB_DIR = os.environ.get("KVXOPT_GLPK_LIB_DIR",GLPK_LIB_DIR)
GLPK_INC_DIR = os.environ.get("KVXOPT_GLPK_INC_DIR",GLPK_INC_DIR)
BUILD_DSDP = int(os.environ.get("KVXOPT_BUILD_DSDP",BUILD_DSDP))
DSDP_LIB_DIR = os.environ.get("KVXOPT_DSDP_LIB_DIR",DSDP_LIB_DIR)
DSDP_INC_DIR = os.environ.get("KVXOPT_DSDP_INC_DIR",DSDP_INC_DIR)
SUITESPARSE_LIB_DIR = os.environ.get("KVXOPT_SUITESPARSE_LIB_DIR",SUITESPARSE_LIB_DIR)
SUITESPARSE_INC_DIR = os.environ.get("KVXOPT_SUITESPARSE_INC_DIR",SUITESPARSE_INC_DIR)
SUITESPARSE_SRC_DIR = os.environ.get("KVXOPT_SUITESPARSE_SRC_DIR",SUITESPARSE_SRC_DIR)
SUITESPARSE_CONFIG = os.environ.get("KVXOPT_SUITESPARSE_CONFIG",SUITESPARSE_CONFIG) == True
BUILD_OSQP = int(os.environ.get("KVXOPT_BUILD_OSQP",BUILD_OSQP))
OSQP_LIB_DIR = os.environ.get("KVXOPT_OSQP_LIB_DIR",OSQP_LIB_DIR)
OSQP_INC_DIR = os.environ.get("KVXOPT_OSQP_INC_DIR",OSQP_INC_DIR)
BUILD_GRB = int(os.environ.get("KVXOPT_BUILD_GRB", BUILD_GRB))
GRB_LIB_DIR = os.environ.get("KVXOPT_GRB_LIB_DIR", GRB_LIB_DIR)
GRB_INC_DIR = os.environ.get("KVXOPT_GRB_INC_DIR", GRB_INC_DIR)
GRB_LIB = os.environ.get("KVXOPT_GRB_LIB",GRB_LIB)
if type(GRB_LIB) is str: GRB_LIB = GRB_LIB.strip().split(';')
MSVC = int(os.environ.get("KVXOPT_MSVC",MSVC)) == True
PYTHON_REQUIRES = (
    '>=3.5, !=3.0.*, !=3.1.*, '
    '!=3.2.*, !=3.3.*, !=3.4.*')
INSTALL_REQUIRES = os.environ.get("CVXOPT_INSTALL_REQUIRES",[])
if type(INSTALL_REQUIRES) is str: INSTALL_REQUIRES = INSTALL_REQUIRES.strip().split(';')

RT_LIB = ["rt"] if sys.platform.startswith("linux") else []
M_LIB = ["m"] if not MSVC else []
UMFPACK_EXTRA_COMPILE_ARGS = ["-Wno-unknown-pragmas"] if not MSVC else []

extmods = []


SUITESPARSE_CONF_LIB = ([], ['suitesparseconfig'])[SUITESPARSE_CONFIG]

if sys.maxsize > 2**31:
    DLONG = True
else:
    DLONG = False


# Macros
MACROS = []
if BLAS_NOUNDERSCORES: MACROS.append(('BLAS_NO_UNDERSCORE',''))

# optional modules

if BUILD_GSL:
    gsl = Extension('gsl', libraries = M_LIB + ['gsl'] + BLAS_LIB,
        include_dirs = [ GSL_INC_DIR ],
        library_dirs = [ GSL_LIB_DIR, BLAS_LIB_DIR ],
        define_macros = GSL_MACROS,
        extra_link_args = BLAS_EXTRA_LINK_ARGS,
        sources = ['src/C/gsl.c'] )
    extmods += [gsl];

if BUILD_FFTW:
    fftw = Extension('fftw', libraries = FFTW_LIB + BLAS_LIB,
        include_dirs = [ FFTW_INC_DIR ],
        library_dirs = [ FFTW_LIB_DIR, BLAS_LIB_DIR ],
        define_macros = FFTW_MACROS,
        extra_link_args = BLAS_EXTRA_LINK_ARGS,
        sources = ['src/C/fftw.c'] )
    extmods += [fftw];

if BUILD_GLPK:
    glpk = Extension('glpk', libraries = ['glpk'],
        include_dirs = [ GLPK_INC_DIR ],
        library_dirs = [ GLPK_LIB_DIR ],
        sources = ['src/C/glpk.c'] )
    extmods += [glpk];

if BUILD_OSQP:
    osqp = Extension('osqp', libraries = ['osqp'],
        include_dirs = [ OSQP_INC_DIR ],
        library_dirs = [ OSQP_LIB_DIR ],
        define_macros = MACROS + [('DDEBUG', ''), ('PRINTING', '')],
        sources = ['src/C/osqp.c'] )
    extmods += [osqp];

if BUILD_GRB:
    grb = Extension('gurobi', libraries=GRB_LIB,
                    include_dirs=[GRB_INC_DIR],
                    library_dirs=[GRB_LIB_DIR],
                    define_macros=MACROS,
                    sources=['src/C/gurobi.c'])
    extmods += [grb]

if BUILD_DSDP:
    dsdp = Extension('dsdp', libraries = ['dsdp'] + LAPACK_LIB + BLAS_LIB,
        include_dirs = [ DSDP_INC_DIR ],
        library_dirs = [ DSDP_LIB_DIR, BLAS_LIB_DIR ],
        extra_link_args = BLAS_EXTRA_LINK_ARGS,
        sources = ['src/C/dsdp.c'] )
    extmods += [dsdp];

# Required modules

base = Extension('base', libraries = M_LIB + LAPACK_LIB + BLAS_LIB,
    library_dirs = [ BLAS_LIB_DIR ],
    define_macros = MACROS,
    extra_link_args = BLAS_EXTRA_LINK_ARGS,
    sources = ['src/C/base.c','src/C/dense.c','src/C/sparse.c'])

blas = Extension('blas', libraries = BLAS_LIB,
    library_dirs = [ BLAS_LIB_DIR ],
    define_macros = MACROS,
    extra_link_args = BLAS_EXTRA_LINK_ARGS,
    sources = ['src/C/blas.c'] )

lapack = Extension('lapack', libraries = LAPACK_LIB + BLAS_LIB,
    library_dirs = [ BLAS_LIB_DIR ],
    define_macros = MACROS,
    extra_link_args = BLAS_EXTRA_LINK_ARGS,
    sources = ['src/C/lapack.c'] )

if not SUITESPARSE_SRC_DIR:
    umfpack = Extension('umfpack',
        libraries = ['umfpack','cholmod','amd','colamd'] + SUITESPARSE_CONF_LIB + LAPACK_LIB + BLAS_LIB + RT_LIB,
        include_dirs = [SUITESPARSE_INC_DIR],
        library_dirs = [SUITESPARSE_LIB_DIR, BLAS_LIB_DIR],
        sources = ['src/C/umfpack.c'])
else:
    umf_sources =  [ 'src/C/umfpack.c',
            SUITESPARSE_SRC_DIR + '/UMFPACK/Source/umfpack_tictoc.c',
            SUITESPARSE_SRC_DIR + '/SuiteSparse_config/SuiteSparse_config.c']

    if DLONG:
        umf_sources += \
        glob(SUITESPARSE_SRC_DIR + '/UMFPACK/Source2/*_l_*.c') +\
        glob(SUITESPARSE_SRC_DIR + '/UMFPACK/Source2/*_dl_*.c') +\
        glob(SUITESPARSE_SRC_DIR + '/UMFPACK/Source2/*_zl_*.c') +\
        glob(SUITESPARSE_SRC_DIR + '/AMD/Source/*_l*.c')
    else:
        umf_sources += \
        glob(SUITESPARSE_SRC_DIR + '/UMFPACK/Source2/*_i_*.c') +\
        glob(SUITESPARSE_SRC_DIR + '/UMFPACK/Source2/*_di_*.c') +\
        glob(SUITESPARSE_SRC_DIR + '/UMFPACK/Source2/*_zi_*.c') +\
        glob(SUITESPARSE_SRC_DIR + '/AMD/Source/*[!_l]*.c')

    umfpack = Extension('umfpack',
        include_dirs = [ SUITESPARSE_SRC_DIR + '/UMFPACK/Include',
            SUITESPARSE_SRC_DIR + '/AMD/Include',
            SUITESPARSE_SRC_DIR + '/UMFPACK/Source',
            SUITESPARSE_SRC_DIR + '/AMD/Source',
            SUITESPARSE_SRC_DIR + '/SuiteSparse_config' ],
        library_dirs = [ BLAS_LIB_DIR ],
        define_macros = MACROS + [('NTIMER', '1'), ('NCHOLMOD', '1')],
        libraries = LAPACK_LIB + BLAS_LIB,
        extra_compile_args = UMFPACK_EXTRA_COMPILE_ARGS,
        extra_link_args = BLAS_EXTRA_LINK_ARGS,
        sources = umf_sources)

if not SUITESPARSE_SRC_DIR:
    klu = Extension('klu',
    libraries=['klu', 'amd', 'colamd', 'btf',] + SUITESPARSE_CONF_LIB + LAPACK_LIB + BLAS_LIB + RT_LIB,
    include_dirs = [SUITESPARSE_INC_DIR],
    define_macros = MACROS,
    library_dirs = [SUITESPARSE_LIB_DIR, BLAS_LIB_DIR],
    sources = ['src/C/klu.c'])
else:

    klu_sources = [ 'src/C/klu.c',
                  SUITESPARSE_SRC_DIR + '/SuiteSparse_config/SuiteSparse_config.c']

    if DLONG:
        klu_sources += (
            glob(SUITESPARSE_SRC_DIR + "/AMD/Source/*_l*.c")
            + glob(SUITESPARSE_SRC_DIR + "/BTF/Source/*_l_*.c")
            + [SUITESPARSE_SRC_DIR + "/COLAMD/Source/colamd_l.c"]
            + [
                SUITESPARSE_SRC_DIR + "/KLU/Source/" + x + ".c"
                for x in [
                    "klu_l",
                    "klu_zl",
                    "klu_l_memory",
                    "klu_l_analyze_given",
                    "klu_l_scale",
                    "klu_zl_scale",
                    "klu_l_kernel",
                    "klu_zl_kernel",
                    "klu_zl_factor",
                    "klu_zl_solve",
                    "klu_l_free_numeric",
                    "klu_zl_free_numeric",
                    "klu_l_analyze",
                    "klu_l_defaults",
                    "klu_l_solve",
                    "klu_l_factor",
                    "klu_zl_extract",
                    "klu_zl_tsolve",
                    "klu_l_tsolve",
                    "klu_l_free_symbolic",
                    "klu_l_extract",
                ]
            ]
        )
    else:
        klu_sources += (
            glob(SUITESPARSE_SRC_DIR + "/AMD/Source/*[!_l]*.c")
            + glob(SUITESPARSE_SRC_DIR + "/BTF/Source/*[!_l_]*.c")
            + [SUITESPARSE_SRC_DIR + "/COLAMD/Source/colamd.c"]
            + [
                SUITESPARSE_SRC_DIR + "/KLU/Source/" + x + ".c"
                for x in [
                    "klu",
                    "klu_z",
                    "klu_memory",
                    "klu_analyze_given",
                    "klu_scale",
                    "klu_z_scale",
                    "klu_analyze",
                    "klu_kernel",
                    "klu_z_kernel",
                    "klu_free_numeric",
                    "klu_z_free_numeric",
                    "klu_defaults",
                    "klu_solve",
                    "klu_tsolve",
                    "klu_z_solve",
                    "klu_z_tsolve",
                    "klu_factor",
                    "klu_z_factor",
                    "klu_extract",
                    "klu_z_extract",
                    "klu_free_symbolic",
                ]
            ]
        )

    klu = Extension('klu',
        include_dirs = [ SUITESPARSE_SRC_DIR + '/KLU/Include',
            SUITESPARSE_SRC_DIR + '/KLU/Source',
            SUITESPARSE_SRC_DIR + '/AMD/Include',
            SUITESPARSE_SRC_DIR + '/AMD/Source',
            SUITESPARSE_SRC_DIR + '/COLAMD/Include',
            SUITESPARSE_SRC_DIR + '/COLAMD/Source',
            SUITESPARSE_SRC_DIR + '/BTF/Include',
            SUITESPARSE_SRC_DIR + '/BTF/Source',
            SUITESPARSE_SRC_DIR + '/SuiteSparse_config' ],
        library_dirs = [ BLAS_LIB_DIR ],
        define_macros = MACROS + [('NTIMER', '1'), ('NCHOLMOD', '1')],
        libraries = LAPACK_LIB + BLAS_LIB,
        extra_link_args = BLAS_EXTRA_LINK_ARGS,
        sources = klu_sources)

# Build for int or long?
if DLONG:
    MACROS += [('DLONG',None)]

if not SUITESPARSE_SRC_DIR:
    cholmod = Extension('cholmod',
        libraries = ['cholmod','colamd','amd'] + SUITESPARSE_CONF_LIB + LAPACK_LIB + BLAS_LIB + RT_LIB,
        include_dirs = [SUITESPARSE_INC_DIR],
        library_dirs = [SUITESPARSE_LIB_DIR, BLAS_LIB_DIR],
        sources = [ 'src/C/cholmod.c' ])
else:
    cholmod_sources = (
        ["src/C/cholmod.c"]
        + [SUITESPARSE_SRC_DIR + "/SuiteSparse_config/SuiteSparse_config.c"]
        + glob(SUITESPARSE_SRC_DIR + "/CHOLMOD/Check/cholmod_*.c")
        + glob(SUITESPARSE_SRC_DIR + "/CHOLMOD/Utility/cholmod_*.c")
        + glob(SUITESPARSE_SRC_DIR + "/CHOLMOD/MatrixOps/cholmod_*.c")
    )

    if DLONG:
        cholmod_sources += \
        [SUITESPARSE_SRC_DIR + '/AMD/Source/amd_l' + s for s in ['_postorder.c', '_post_tree.c', '2.c']] +\
        [SUITESPARSE_SRC_DIR + '/COLAMD/Source/colamd_l.c'] +\
        glob(SUITESPARSE_SRC_DIR + '/CHOLMOD/Core/cholmod_l_*.c') +\
        glob(SUITESPARSE_SRC_DIR + '/CHOLMOD/Cholesky/cholmod_l_*.c') +\
        glob(SUITESPARSE_SRC_DIR + '/CHOLMOD/Supernodal/c*_l_*.c')
    else:
        cholmod_sources += \
        [SUITESPARSE_SRC_DIR + '/AMD/Source/amd_' + s for s in ['postorder.c', 'post_tree.c', '2.c']] +\
        [SUITESPARSE_SRC_DIR + '/COLAMD/Source/colamd.c'] +\
        glob(SUITESPARSE_SRC_DIR + '/CHOLMOD/Core/cholmod_[!l_]*.c') +\
        glob(SUITESPARSE_SRC_DIR + '/CHOLMOD/Cholesky/cholmod_[!l_]*.c') +\
        glob(SUITESPARSE_SRC_DIR + '/CHOLMOD/Supernodal/cholmod_[!l_]*.c')

    cholmod = Extension('cholmod',
        library_dirs = [ BLAS_LIB_DIR ],
        libraries = LAPACK_LIB + BLAS_LIB,
        include_dirs = [ SUITESPARSE_SRC_DIR + '/CHOLMOD/Include',
            SUITESPARSE_SRC_DIR + '/COLAMD',
            SUITESPARSE_SRC_DIR + '/AMD/Include',
            SUITESPARSE_SRC_DIR + '/COLAMD/Include',
            SUITESPARSE_SRC_DIR + '/SuiteSparse_config' ],
        define_macros = MACROS + [('NPARTITION', '1'), ('NTIMER', '1')],
        extra_link_args = BLAS_EXTRA_LINK_ARGS,
        sources = cholmod_sources)

if not SUITESPARSE_SRC_DIR:
    amd = Extension('amd',
        libraries = ['amd'] + SUITESPARSE_CONF_LIB + RT_LIB,
        include_dirs = [SUITESPARSE_INC_DIR],
        library_dirs = [SUITESPARSE_LIB_DIR],
        sources = ['src/C/amd.c'])
else:
    amd_sources = [ 'src/C/amd.c', SUITESPARSE_SRC_DIR + '/SuiteSparse_config/SuiteSparse_config.c']

    if DLONG:
        amd_sources += glob(SUITESPARSE_SRC_DIR + '/AMD/Source/*_l*.c')
    else:
        amd_sources += glob(SUITESPARSE_SRC_DIR + '/AMD/Source/*[!_l]*.c')


    amd = Extension('amd',
        include_dirs = [SUITESPARSE_SRC_DIR + '/AMD/Include',
            SUITESPARSE_SRC_DIR + '/SuiteSparse_config' ],
        define_macros = MACROS + [('NTIMER', '1')],
        sources = amd_sources )

misc_solvers = Extension('misc_solvers',
    libraries = LAPACK_LIB + BLAS_LIB,
    library_dirs = [ BLAS_LIB_DIR ],
    define_macros = MACROS,
    extra_link_args = BLAS_EXTRA_LINK_ARGS,
    sources = ['src/C/misc_solvers.c'] )

extmods += [base, blas, lapack, umfpack, klu, cholmod, amd, misc_solvers]


setup (name = 'kvxopt',
    description = 'Convex optimization package and Suite Sparse interface',
    long_description = '''
KVXOPT is a fork from CVXOPT wich contains more functions and
wrappers to Suite Sparse library.

CVXOPT is a free software package for convex optimization based on the
Python programming language. It can be used with the interactive Python
interpreter, on the command line by executing Python scripts, or
integrated in other software via Python extension modules. Its main
purpose is to make the development of software for convex optimization
applications straightforward by building on Python's extensive standard
library and on the strengths of Python as a high-level programming
language.
''',
    author = 'M. Andersen, J. Dahl, L. Vandenberghe, and U. Sandoval',
    author_email = 'martin.skovgaard.andersen@gmail.com, dahl.joachim@gmail.com, vandenbe@ee.ucla.edu, sanurielf@gmail.com',
    url = '',
    project_urls = {'Source': 'https://github.com/sanurielf/kvxopt'},
    license = 'GNU GPL version 3',
    ext_package = "kvxopt",
    ext_modules = extmods,
    package_dir = {"kvxopt": "src/python"},
    package_data = {'': ['.libs/*.dll', '*LICENSE']},
    packages = ["kvxopt"],
    python_requires=PYTHON_REQUIRES,
    install_requires = INSTALL_REQUIRES,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        'Programming Language :: C',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
        ],
    zip_safe=False
    )
