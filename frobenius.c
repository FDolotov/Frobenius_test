#include "frobenius.h"
#include "primes.h"
#include <assert.h>

void set_nums()
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


//Step 1: test n for divisibility by primes <= square root of num
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

int jacobi(const gcry_mpi_t q, const gcry_mpi_t p) 
{
	gcry_mpi_t n = gcry_mpi_new(0);
	gcry_mpi_set(n, q);
	gcry_mpi_t k = gcry_mpi_new(0);
	gcry_mpi_set(k, p);

	gcry_mpi_t buff = gcry_mpi_new(0);
	gcry_mpi_t buff2 = gcry_mpi_new(0);
	gcry_mpi_t r = gcry_mpi_new(0);
	gcry_mpi_t four = gcry_mpi_new(0);
	gcry_mpi_t eight = gcry_mpi_new(0);
	gcry_mpi_set_ui(four, 4);
	gcry_mpi_set_ui(eight, 8);

	gcry_mpi_mod(buff, k, two);

	assert(gcry_mpi_cmp_ui(buff, 1) == 0);

	gcry_mpi_mod(n, n, k);
    
	int t = 1;

	while (gcry_mpi_cmp_ui(n, 0) != 0) 
	{
		gcry_mpi_mod(buff2, n, two);
        	while (gcry_mpi_cmp_ui(buff2, 0) == 0) 
		{
			gcry_mpi_div(n, NULL, n, two, 0);           
			gcry_mpi_mod(r, k, eight);
        
			if (gcry_mpi_cmp_ui(r, 3) == 0 || gcry_mpi_cmp_ui(r, 5) == 0)
			{
				t = -t;
			}
			gcry_mpi_mod(buff2, n, two);
		}

		gcry_mpi_set(buff, n);
		gcry_mpi_set(n, k);
		gcry_mpi_set(k, buff);

		gcry_mpi_mod(buff, n, four);
		gcry_mpi_mod(buff2, k, four);
	
		if (gcry_mpi_cmp_ui(buff, 3) == 0 && gcry_mpi_cmp_ui(buff2, 3) == 0)
		{
			t = -t;		
		}
		
		gcry_mpi_mod(n, n, k);
	}

	t = (gcry_mpi_cmp(k, one) == 0) ? t : 0;
	return t;
	
	gcry_mpi_release(n);
	gcry_mpi_release(k);
	gcry_mpi_release(buff);
	gcry_mpi_release(buff2);
	gcry_mpi_release(r);
	gcry_mpi_release(four);
	gcry_mpi_release(eight);
}

void set_params(struct params *p, const gcry_mpi_t n)
{
	p->b = gcry_mpi_new(0);
	p->c = gcry_mpi_new(0);

	gcry_mpi_set_ui(p->b, 0);
	gcry_mpi_set_ui(p->c, 0);

	int size = gcry_mpi_get_nbits(n);
	
	do
	{
		gcry_mpi_randomize(p->c, size, GCRY_STRONG_RANDOM);
		gcry_mpi_mod(p->c, p->c, n); //if c would be greater than n

		gcry_mpi_mulm(p->c, p->c, p->c, n);
		gcry_mpi_sub(p->c, n, p->c);
	}
	while(gcry_mpi_cmp_ui(p->c, 3) < 0);

	gcry_mpi_t bc = gcry_mpi_new(0); //for b^2 + 4c
	gcry_mpi_t buff = gcry_mpi_new(0);

	do
	{
		gcry_mpi_randomize(p->b, size, GCRY_STRONG_RANDOM);
		gcry_mpi_mod(p->b, p->b, n);

		gcry_mpi_mul(bc, p->b, p->b);
		gcry_mpi_mul_ui(buff, p->c, 4);
		gcry_mpi_add(bc, bc, buff);

	}
	while(jacobi(bc, n) != -1);

	gcry_mpi_release(buff);
	gcry_mpi_release(bc);
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
	gcry_mpi_mulm(buff2, f_1, g_1, n); // f_1 * g_1
	gcry_mpi_addm(rez_1, buff, buff2, n);
};*/

//Calculate square of f(x) = (f_x*x + f_1)^2 mod (n, x^2 - b*x - c)

/*void square_mod(gcry_mpi_t rez_x, gcry_mpi_t res_1, const gcry_mpi_t f_x, const gcry_mpi_t f_1)
{
	if(gcry_mpi_cmp_ui(f_x, 0) == 0)
	{
		gcry_mpi_set_ui(res_x, 0);
		gcry_mpi_mulm(res_1, , f_1, f_1, n);

		return;
	}	

	//res_x = b * f_x^2 + 2 * f_x * f_1 % n
	gcry_mpi_t buff = gcry_mpi_new(0);
	gcry_mpi_powm(buff, f_x, two, n); // f_x^2
	gcry_mpi_mulm(buff, buff, b, n); // b * f_x^2
	gcry_mpi_mulm(rez_x, two, f_x, n); // 2 * f_x
	gcry_mpi_mulm(rez_x, rez_x, f_1, n); // 2 * f_x * f_1
	gcry_mpi_addm(rez_x, fez_x, buff, n);

	//rez_1 = c * f_x^2 + f_1 ^ 2
	gcry_mpi_powm(buff, f_x, two, n); // f_x^2
	gcry_mpi_mulm(buff, buff, c, n); // c * f_x^2
       	gcry_mpi_powm(res_1, f_1, f_1, n); // f_1^2
	gcry_mpi_addm(rez_1, buff, rez_1, n);
};*/
