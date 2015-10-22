#.RECIPEPREFIX = >
.ONESHELL:

Riemann_solver/Riemann_solver.a: Riemann_solver/Riemann_solver_exact.c \
                                 Riemann_solver/Riemann_solver_Toro.c \
                                 Riemann_solver/linear_GRP_solver.c \
                                 Riemann_solver/linear_GRP_solver_Li.c \
                                 Riemann_solver/GRPSUB.c
	cd ./Riemann_solver/; \
         gcc -c ./Riemann_solver_exact.c -g; \
         gcc -c ./Riemann_solver_Toro.c -I ../ -g; \
         gcc -c ./linear_GRP_solver.c -I ../ -g; \
         gcc -c ./linear_GRP_solver_Li.c -I ../ -g; \
         gcc -c ./GRPSUB.c -I ../ -g; \
         ar crv Riemann_solver.a Riemann_solver_exact.o Riemann_solver_Toro.o linear_GRP_solver.o linear_GRP_solver_Li.o GRPSUB.o; \
         ranlib Riemann_solver.a; \
         cd ..

file_io/file_io.a: file_io/runhist.c \
                   file_io/str_opt.c \
                   file_io/file_i.c \
                   file_io/file_o.c
	cd ./file_io/; \
         gcc -c ./runhist.c -I ../ -g; \
         gcc -c ./str_opt.c        -g; \
         gcc -c ./file_i.c -I ../  -g; \
         gcc -c ./file_o.c -I ../  -g; \
         ar crv file_io.a str_opt.o file_i.o file_o.o runhist.o; \
         ranlib file_io.a; \
         cd ..


.PHONY: clean

clean:
	rm -f ./file_io/*.[oa]; \
        rm -f ./Riemann_solver/*.[oa]; \
        rm -f ./reconstruction/*.[oa]; \
        rm -f ./FD_WENO/*.[oa]; \
        rm -f ./GRP_fix/*.[oa]; \
        rm -f ./GRP_HWENO_fix/*.[oa]; \
        rm -f *.[oa]; \
        rm -f MainEngine;

