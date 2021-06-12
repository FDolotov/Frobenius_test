#include <stdio.h>
#include "frobenius.h"

gcry_mpi_t zero;
gcry_mpi_t one;
gcry_mpi_t two;

static void print_table(unsigned kmax, unsigned nmax) 
{
	printf("n\\k|");
		for (int k = 0; k <= kmax; ++k) 
			printf("%'3u", k);
	printf("\n----");

	for (int k = 0; k <= kmax; --k) 
			printf("---");
	putchar('\n');

	gcry_mpi_t buff = gcry_mpi_new(0);
	gcry_mpi_t buff2 = gcry_mpi_new(0);

	for (int n = 1; n <= nmax; n += 2) 
	{
		printf("%-2u |", n);
		for (int k = 0; k <= kmax; ++k)
		{
			gcry_mpi_set_ui(buff, k);
			gcry_mpi_set_ui(buff2, n);
			printf("%'3d", jacobi(buff, buff2));
		};
		putchar('\n');
	};

	gcry_mpi_release(buff);
	gcry_mpi_release(buff2);
};

int main()
{
	set_nums();

	gcry_mpi_t buff = gcry_mpi_new(0);
	gcry_mpi_scan(&buff, GCRYMPI_FMT_HEX, "FFFF", 0, 0);

	gcry_mpi_t buff2 = gcry_mpi_new(0);
	gcry_mpi_scan(&buff2, GCRYMPI_FMT_HEX, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD97", 0, 0);

	gcry_mpi_t buff3 = gcry_mpi_new(0);
	gcry_mpi_scan(&buff3, GCRYMPI_FMT_HEX, "605F6B7C183FA81578BC39CFAD518132B9DF62897009AF7E522C32D6DC7BFFB", 0, 0); 
	
	printf("\n");
	gcry_mpi_dump(buff);
	printf("\n");
	step_1(buff);
	step_2(buff);
	
	printf("\n");
	gcry_mpi_dump(buff2);
	printf("\n");
	step_1(buff2);
	step_2(buff2);

	printf("\n");
	gcry_mpi_dump(buff3);
	printf("\n");
	step_1(buff3);
	step_2(buff3);

	//gcry_mpi_set_ui(buff, 7);
	//gcry_mpi_set_ui(buff2, 4);

	//int j = jacobi(buff, buff2);
	//int t = jacobi(one, two);

	printf("\n\n");
	print_table(20, 21);
	printf("\n");
	//gcry_mpi_set_ui(buff,3);
	//gcry_mpi_set_ui(buff2, 5);

	//jacobi(buff,buff2);

	struct params p;
	set_params(&p, buff2);
	
	gcry_mpi_dump(p.c);
	printf("\n");
	gcry_mpi_dump(p.b);
	printf("\n");


	gcry_mpi_release(buff);
	gcry_mpi_release(buff2);
	gcry_mpi_release(buff3);
	release_params(&p);
	release_memory();
};
