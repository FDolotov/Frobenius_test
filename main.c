#include <stdio.h>
#include "frobenius.h"
#include "primes.h"

gcry_mpi_t zero;
gcry_mpi_t one;
gcry_mpi_t two;

int main()
{
	gcry_mpi_t buff = gcry_mpi_new(0);
	
	u_int16_t i;
	
	for(i = 0; i < 300; i++)
	{
		//set_params();
		gcry_mpi_set_ui(buff, prime_list[i]);
		gcry_mpi_dump(buff);
		printf("\n");
		//release_memory();
	};
	gcry_mpi_release(buff);
};
