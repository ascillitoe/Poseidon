source ../libs/OPS/scripts/source_pgilaptop
make clean
../libs/OPS/translator/python/fortran/ops_fortran.py main.F90
make
 mv *.optrpt reports/
