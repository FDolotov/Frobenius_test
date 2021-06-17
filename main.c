#include <stdio.h>
#include "frobenius.h"

gcry_mpi_t zero;
gcry_mpi_t one;
gcry_mpi_t two;

int main()
{
	set_nums();
	size_t scanned;

	gcry_mpi_t buff = gcry_mpi_new(0);
	gcry_mpi_scan(&buff, GCRYMPI_FMT_HEX, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF333333333333", 0, &scanned);

	gcry_mpi_t buff2 = gcry_mpi_new(0);
	gcry_mpi_scan(&buff2, GCRYMPI_FMT_HEX, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD97", 0, &scanned);

	gcry_mpi_t buff3 = gcry_mpi_new(0);
	gcry_mpi_scan(&buff3, GCRYMPI_FMT_HEX, "688589CC0E9505E2F2FEE557FFFFFFF", 0, &scanned); 
	
	struct params p;
	set_params(&p, buff);

	printf("\n");
	gcry_mpi_dump(buff);
	printf("\nAnticipated result: composite\nResult: ");
	QFT(buff, &p);
	release_params(&p);
	
	set_params(&p, buff2);
	printf("\n");
	gcry_mpi_dump(buff2);
	printf("\nAnticipated result: prime\nResult: ");
	QFT(buff2, &p);
	release_params(&p);

	printf("\n");
	gcry_mpi_dump(buff3);
	printf("\nAnticipated result: prime\nResult: ");
	RQFT(buff3);

	gcry_mpi_release(buff);
	gcry_mpi_release(buff2);
	gcry_mpi_release(buff3);
	release_memory();
};
