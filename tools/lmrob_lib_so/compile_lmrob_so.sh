#!/bin/bash
pushd tools/lmrob_lib_so/src

gcc -std=gnu99 -I../inc -DNDEBUG -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-3.3.3=. -fstack-protector-strong -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 -g  -c R-rng4ftn.c -o R-rng4ftn.o
gcc -std=gnu99 -I../inc -DNDEBUG -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-3.3.3=. -fstack-protector-strong -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 -g  -c init.c -o init.o
gcc -std=gnu99 -I../inc -DNDEBUG -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-3.3.3=. -fstack-protector-strong -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 -g  -c lmrob.c -o lmrob.o
gcc -std=gnu99 -I../inc -DNDEBUG -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-3.3.3=. -fstack-protector-strong -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 -g  -c mc.c -o mc.o
gcc -std=gnu99 -I../inc -DNDEBUG -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-3.3.3=. -fstack-protector-strong -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 -g  -c monitor.c -o monitor.o
gcc -std=gnu99 -I../inc -DNDEBUG -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-3.3.3=. -fstack-protector-strong -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 -g  -c qn_sn.c -o qn_sn.o
gcc -std=gnu99 -I../inc -DNDEBUG -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-3.3.3=. -fstack-protector-strong -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 -g  -c rob-utils.c -o rob-utils.o
gcc -std=gnu99 -I../inc -DNDEBUG -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-3.3.3=. -fstack-protector-strong -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 -g  -c rowMedians.c -o rowMedians.o
gcc -std=gnu99 -I../inc -DNDEBUG -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-3.3.3=. -fstack-protector-strong -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 -g  -c wgt_himed.c -o wgt_himed.o

gfortran -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-3.3.3=. -fstack-protector-strong  -c eigen.f -o eigen.o
gfortran -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-3.3.3=. -fstack-protector-strong  -c rf-common.f -o rf-common.o
gfortran -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-3.3.3=. -fstack-protector-strong  -c rffastmcd.f -o rffastmcd.o
gfortran -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-3.3.3=. -fstack-protector-strong  -c rfltsreg.f -o rfltsreg.o
gfortran -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-3.3.3=. -fstack-protector-strong  -c rllarsbi.f -o rllarsbi.o

gfortran -Wall -g -ffixed-form -fno-second-underscore -Wall -g -fno-second-underscore -fPIC -O3 -funroll-loops -c blas.f -o blas.o
gfortran -Wall -g -ffixed-form -fno-second-underscore -Wall -g -fno-second-underscore -fPIC -O3 -funroll-loops -c cmplxblas.f -o cmplxblas.o
gfortran -Wall -g -ffixed-form -fno-second-underscore -Wall -g -fno-second-underscore -fPIC -O3 -funroll-loops -c dchdc.f -o dchdc.o
gfortran -Wall -g -ffixed-form -fno-second-underscore -Wall -g -fno-second-underscore -fPIC -O3 -funroll-loops -c dpbfa.f -o dpbfa.o
gfortran -Wall -g -ffixed-form -fno-second-underscore -Wall -g -fno-second-underscore -fPIC -O3 -funroll-loops -c dpbsl.f -o dpbsl.o
gfortran -Wall -g -ffixed-form -fno-second-underscore -Wall -g -fno-second-underscore -fPIC -O3 -funroll-loops -c dpoco.f -o dpoco.o
gfortran -Wall -g -ffixed-form -fno-second-underscore -Wall -g -fno-second-underscore -fPIC -O3 -funroll-loops -c dpodi.f -o dpodi.o
gfortran -Wall -g -ffixed-form -fno-second-underscore -Wall -g -fno-second-underscore -fPIC -O3 -funroll-loops -c dpofa.f -o dpofa.o
gfortran -Wall -g -ffixed-form -fno-second-underscore -Wall -g -fno-second-underscore -fPIC -O3 -funroll-loops -c dposl.f -o dposl.o
gfortran -Wall -g -ffixed-form -fno-second-underscore -Wall -g -fno-second-underscore -fPIC -O3 -funroll-loops -c dqrdc2.f -o dqrdc2.o
gfortran -Wall -g -ffixed-form -fno-second-underscore -Wall -g -fno-second-underscore -fPIC -O3 -funroll-loops -c dqrdc.f -o dqrdc.o
gfortran -Wall -g -ffixed-form -fno-second-underscore -Wall -g -fno-second-underscore -fPIC -O3 -funroll-loops -c dqrls.f -o dqrls.o
gfortran -Wall -g -ffixed-form -fno-second-underscore -Wall -g -fno-second-underscore -fPIC -O3 -funroll-loops -c dqrsl.f -o dqrsl.o
gfortran -Wall -g -ffixed-form -fno-second-underscore -Wall -g -fno-second-underscore -fPIC -O3 -funroll-loops -c dqrutl.f -o dqrutl.o
gfortran -Wall -g -ffixed-form -fno-second-underscore -Wall -g -fno-second-underscore -fPIC -O3 -funroll-loops -c dsvdc.f -o dsvdc.o
gfortran -Wall -g -ffixed-form -fno-second-underscore -Wall -g -fno-second-underscore -fPIC -O3 -funroll-loops -c dtrco.f -o dtrco.o
gfortran -Wall -g -ffixed-form -fno-second-underscore -Wall -g -fno-second-underscore -fPIC -O3 -funroll-loops -c dtrsl.f -o dtrsl.o
gfortran -Wall -g -ffixed-form -fno-second-underscore -Wall -g -fno-second-underscore -fPIC -O3 -funroll-loops -c xerbla.f -o xerbla.o

gcc -std=gnu99 -shared blas.o cmplxblas.o dchdc.o dpbfa.o dpbsl.o dpoco.o dpodi.o dpofa.o dposl.o dqrdc2.o dqrdc.o dqrls.o dqrsl.o dqrutl.o dsvdc.o dtrco.o dtrsl.o xerbla.o -L/usr/lib64/R/lib -Wl,-z,relro -o lmrob.so R-rng4ftn.o eigen.o init.o lmrob.o mc.o monitor.o qn_sn.o rf-common.o rffastmcd.o rfltsreg.o rllarsbi.o rob-utils.o rowMedians.o wgt_himed.o -llapack -lblas -lgfortran -lm -lquadmath -lgfortran -lm -lquadmath -L/usr/lib/R/lib -lR
popd
cp tools/lmrob_lib_so/src/lmrob.so lmrob/liblmrob.so
rm tools/lmrob_lib_so/src/*.o -v
rm tools/lmrob_lib_so/src/*.so -v
