#include <stdio.h>
#include "frobenius.h"

gcry_mpi_t zero;
gcry_mpi_t one;
gcry_mpi_t two;

int main()
{
	set_nums();

	gcry_mpi_t buff = gcry_mpi_new(0);
	gcry_mpi_scan(&buff, GCRYMPI_FMT_HEX, "B", 0, 0);

	gcry_mpi_t buff2 = gcry_mpi_new(0);
	gcry_mpi_scan(&buff2, GCRYMPI_FMT_HEX, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD97", 0, 0);

	gcry_mpi_t buff3 = gcry_mpi_new(0);
	gcry_mpi_scan(&buff3, GCRYMPI_FMT_HEX, "605F6B7C183FA81578BC39CFAD518132B9DF62897009AF7E522C32D6DC7BFFB", 0, 0); 
	
	struct params p;
	set_params(&p, buff);

	printf("\n");
	gcry_mpi_dump(buff);
	printf("\n");
	QFT(buff, &p);
	
	set_params(&p, buff2);
	printf("\n");
	gcry_mpi_dump(buff2);
	printf("\n");
	QFT(buff2, &p);

	set_params(&p, buff3);
	printf("\n");
	gcry_mpi_dump(buff3);
	printf("\n");
	RQFT(buff3);

	gcry_mpi_release(buff);
	gcry_mpi_release(buff2);
	gcry_mpi_release(buff3);
	release_params(&p);
	release_memory();
};
