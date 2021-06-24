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

	gcry_mpi_t buff4 = gcry_mpi_new(0);
    	gcry_mpi_scan(&buff4, GCRYMPI_FMT_HEX, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDC7", 0, &scanned);

	gcry_mpi_t buff5 = gcry_mpi_new(0);
   	gcry_mpi_scan(&buff5, GCRYMPI_FMT_HEX, "9E4F5D8C017D8D9F13A5CF3CDF5BFE4DAB402D54198E31EBDE28A0621050439CA6B39E0A515C06B304E2CE43E79E369E91A0CFC2BC2A22B4CA302DBB33EE7551", 0, &scanned);    
    
    	gcry_mpi_t buff6 = gcry_mpi_new(0);
    	gcry_mpi_scan(&buff6, GCRYMPI_FMT_HEX, "3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC98CDBA46506AB004C33A9FF5147502CC8EDA9E7A769A12694623CEF47F023ED", 0, &scanned);

	gcry_mpi_t buff7 = gcry_mpi_new(0);
    	gcry_mpi_scan(&buff7, GCRYMPI_FMT_HEX, "469AF79D1FB1F5E16B99592B77A01E2A0FDFB0D01794368D9A56117F7B38669522DD4B650CF789EEBF068C5D139732F0905622C04B2BAAE7600303EE73001A3D", 0, &scanned);
	
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

	printf("\n");
	gcry_mpi_dump(buff4);
	printf("\nAnticipated result: prime\nResult: ");
	RQFT(buff4);

	printf("\n");
	gcry_mpi_dump(buff5);
	printf("\nAnticipated result: composite\nResult: ");
	RQFT(buff5);

	printf("\n");
	gcry_mpi_dump(buff6);
	printf("\nAnticipated result: prime\nResult: ");
	RQFT(buff6);

	printf("\n");
	gcry_mpi_dump(buff7);
	printf("\nAnticipated result: composite\nResult: ");
	RQFT(buff7);

	gcry_mpi_release(buff);
	gcry_mpi_release(buff2);
	gcry_mpi_release(buff3);
	gcry_mpi_release(buff4);
	gcry_mpi_release(buff5);
	gcry_mpi_release(buff6);
	gcry_mpi_release(buff7);
	release_memory();
};
