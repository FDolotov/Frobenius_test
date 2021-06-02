#ifndef FROBENIUS_H
#define FROBENIUS_H

#include <gcrypt.h>

extern gcry_mpi_t zero;
extern gcry_mpi_t one;
extern gcry_mpi_t two;

struct params 
{
	gcry_mpi_t b;
	gcry_mpi_t c;
};

extern void set_nums();
extern void release_memory();
extern void set_params(struct params*, const gcry_mpi_t);
extern void square_root(gcry_mpi_t, const gcry_mpi_t);
extern void step_1(const gcry_mpi_t);
extern void step_2(const gcry_mpi_t);
extern int jacobi(const gcry_mpi_t,const gcry_mpi_t);

//void mult_mod();

#endif
