#ifndef FROBENIUS_H
#define FROBENIUS_H

#include <gcrypt.h>

extern gcry_mpi_t zero;
extern gcry_mpi_t one;
extern gcry_mpi_t two;

extern void set_params();
extern void release_memory();

//void mult_mod();

#endif
