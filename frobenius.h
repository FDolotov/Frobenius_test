#ifndef FROBENIUS_H
#define FROBENIUS_H

#include <gcrypt.h>

extern gcry_mpi_t zero;
extern gcry_mpi_t one;
extern gcry_mpi_t two;

extern gcry_mpi_t n;
extern gcry_mpi_t b;
extern gcry_mpi_t c;

extern void set_params();
extern void release_memory();
//extern int is_square(gcry_mpi_t);
extern void step_1(const gcry_mpi_t);
extern void step_2(const gcry_mpi_t);

//void mult_mod();

#endif
