#include "frobenius.h"
#include "primes.h"

u_int64_t multipl;

void set_params()
{
	size_t scanned;

	zero = gcry_mpi_new(0);
	gcry_mpi_scan(&zero, GCRYMPI_FMT_HEX, "0", 0, &scanned);

	one = gcry_mpi_new(0);
	gcry_mpi_scan(&one, GCRYMPI_FMT_HEX, "1", 0, &scanned);

	two = gcry_mpi_new(0);
	gcry_mpi_scan(&two, GCRYMPI_FMT_HEX, "2", 0, &scanned);
};

void release_memory()
{
	gcry_mpi_release(zero);
    	gcry_mpi_release(one);
    	gcry_mpi_release(two);
};

void square_root(gcry_mpi_t res, const gcry_mpi_t num)
{
	int bit = gcry_mpi_get_nbits(num);
	bit /= 2;
	
	gcry_mpi_t buff = gcry_mpi_new(0);
	gcry_mpi_t bit_t = gcry_mpi_new(0);
	gcry_mpi_t current_t = gcry_mpi_new(0);
	gcry_mpi_set_ui(bit_t, bit);
	gcry_mpi_set(current_t, bit_t);

	while(1)
	{
		gcry_mpi_div(buff, NULL, num, bit_t, 0);
		gcry_mpi_add(buff, bit_t, buff);
		gcry_mpi_rshift(buff, buff, 1);

		if(gcry_mpi_cmp(buff, bit_t) == 0 || gcry_mpi_cmp(buff, current_t) == 0)
		{
			gcry_mpi_set(res, buff);
			return;
		}

		gcry_mpi_set(current_t, bit_t);
		gcry_mpi_set(bit_t, buff);
	}

	gcry_mpi_release(buff);
	gcry_mpi_release(bit_t);
	gcry_mpi_release(current_t);
};


//Step 1: test n for divisibility by primes <= 1 + n/2
void step_1(const gcry_mpi_t num)
{
	if(gcry_mpi_cmp(zero, num) == 0 || gcry_mpi_cmp(one, num) == 0 || gcry_mpi_cmp(two, num) == 0)
	{
		printf("Number is prime, ");
		return;
	}

	gcry_mpi_t buff = gcry_mpi_new(0);
	square_root(buff, num);

	gcry_mpi_t prime = gcry_mpi_new(0);
	gcry_mpi_set_ui(prime, prime_list[0]);

	gcry_mpi_t mod = gcry_mpi_new(0);

	for(int i = 1; gcry_mpi_cmp(num, prime) >= 0; i++)
	{
		if(i > 5133)
		{
			break;
		}
	
		if(gcry_mpi_cmp(prime, num) == 0)
		{
			printf("Number is prime, ");
			return;
		}

		gcry_mpi_set_ui(prime, prime_list[i]);
	};

	gcry_mpi_set_ui(prime, prime_list[0]);

	//i will start from 1, because number should be odd
	//prime_list[1] is 3
	
	for(int i = 0; gcry_mpi_cmp(buff, prime) >= 0; i++)
	{
		if(i > 5133)
		{
			break;
		}

		gcry_mpi_mod(mod, num, prime);

		if(gcry_mpi_cmp(zero, mod) == 0)
		{
			printf("Number is composite, ");
			return;
		}
		else
		{
			gcry_mpi_set_ui(prime, prime_list[i]);
		}
	}

	printf("Number is probably prime, ");

	gcry_mpi_release(buff);
	gcry_mpi_release(prime);
	gcry_mpi_release(mod);
};

//Step 2: test whether sqrt(n) is integer
void step_2 (const gcry_mpi_t num)
{
	gcry_mpi_t buff = gcry_mpi_new(0);

	square_root(buff, num);
	gcry_mpi_mul(buff, buff, buff);
	
	if(gcry_mpi_cmp(buff, num) == 0)
	{
		printf("Number is not prime\n");
	}
	else
	{
		printf("Number is probably prime\n");
	}

	gcry_mpi_release(buff);
};


/*Calculate f(x) * g(x) mod (n, x^2 - b*x - c), f(x) = f_x*x + f_1, g(x) = g_x*x + g_1
void mult_mod(gcry_mpi_t rez_x, gcry_mpi_t rez_1, const gcry_mpi_t f_x, const gcry_mpi_t f_1, const gcry_mpi_t g_x, const gcry_mpi_t g_1)
{
	if(f_x == zero)
	{
		gcry_mpi_mulm(rez_x, f_1, g_x, n);
		gcry_mpi_mulm(rez_1, f_1, g_x, n);

		multipl =+ 2;

		return;
	};

	gcry_mpi_t buff = gcry_mpi_new(0);
	gcry_mpi_t buff2 = gcry_mpi_new(0);
	
	//rez_x = (f_x * g_x * b + f_x * g_1 + f_1 * g_x) % n

	gcry_mpi_mulm(buff, f_x, g_x, n); // f_x * g_x
	gcry_mpi_mulm(buff, buff, b, n); // f_x * g_x * b
	gcry_mpi_mulm(buff2, f_x, g_1, n); // f_x * g_1
	gcry_mpi_addm(buff, buff, buff2, n); // f_x * g_x * b + f_x *g_1
	gcry_mpi_mulm(buff2, f_x, g_x, n); // f_1 * g_x
	gcry_mpi_addm(rez_x, buff, buff2, n);

	//rez_1 = (f_x * g_x * c + f_1 * g_1) % n
	
	gcry_mpi_mulm(buff, f_x, g_x, n); // f_x * g_x
	gcry_mpi_mulm(buff, buff, c, n); // f_X * g_x * c
	gcry_mpi_mulm(buff2, f_1, g_1, n); //f_1 * g_1
	gcry_mpi_addm(rez_1, buff, buff2, n); //
};*/
