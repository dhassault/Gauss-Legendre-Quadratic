# TP1 quadrature de Gauss-Legendre
FC = gfortran
FFLAGS =
LDFLAGS =

SRCS = gauss_legendre.f90
OBJS = $(SRCS:.f90=.o)
EXEC = $(SRCS:.f90=)

all: $(EXEC)

%.o: %.f90
					$(FC) $(FFLAGS) -c $<

%: %.o
					$(FC) $(FFLAGS) -o $@ $^

gauss_legendre: functions.o

clean :
					rm -f $(OBJS) $(EXEC) functions.mod functions.o
