#include "frobenius.h"
#include "primes.h"
#include <assert.h>

void set_nums()
{
	zero = gcry_mpi_new(0);
	gcry_mpi_set_ui(zero, 0);

	one = gcry_mpi_new(0);
	gcry_mpi_set_ui(one, 1);

	two = gcry_mpi_new(0);
	gcry_mpi_set_ui(two, 2);
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
			break;
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
	int primality = -1; //1 - prime, -1 - composite

	gcry_mpi_t buff = gcry_mpi_new(0);
	square_root(buff, num);

	gcry_mpi_t prime = gcry_mpi_new(0);
	gcry_mpi_set_ui(prime, prime_list[0]);

	gcry_mpi_t mod = gcry_mpi_new(0);

	for(int i = 1; gcry_mpi_cmp(num, prime) >= 0; i++)
	{
		if(i > 5133)
			break;
	
		if(gcry_mpi_cmp(prime, num) == 0)
		{
			primality = 1;
			break;
		}

		gcry_mpi_set_ui(prime, prime_list[i]);
	};

	gcry_mpi_set_ui(prime, prime_list[0]);

	//i will start from 1, because number should be odd
	//prime_list[1] is 3
	
	for(int i = 0; gcry_mpi_cmp(buff, prime) >= 0; i++)
	{
		if(i > 5133)
			break;

		gcry_mpi_mod(mod, num, prime);

		if(gcry_mpi_cmp(zero, mod) == 0)
		{
			primality = -1;
			break;
		}
		else
		{
			primality = 1;
			gcry_mpi_set_ui(prime, prime_list[i]);
		}
	}

	if(primality == 1)
		printf("Number is probably prime, ");
	else
		printf("Number is composite, ");

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
	
	gcry_mpi_release(n);
	gcry_mpi_release(k);
	gcry_mpi_release(buff);
	gcry_mpi_release(buff2);
	gcry_mpi_release(r);
	gcry_mpi_release(four);
	gcry_mpi_release(eight);
	
	return t;
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

void release_params(struct params* p)
{
	gcry_mpi_release(p->c);
	gcry_mpi_release(p->b);
};


//Proposion 3.2
//Calculate f(x) * g(x) mod (n, x^2 - b*x - c), f(x) = f*x + g, g(x) = d*x + e

static void mult_mod(gcry_mpi_t rez_x, gcry_mpi_t rez_1, const gcry_mpi_t f, const gcry_mpi_t g, const gcry_mpi_t d, const gcry_mpi_t e, const gcry_mpi_t n, struct params *p)
{
	if(f == zero)
	{
		gcry_mpi_mulm(rez_x, g, d, n);
		gcry_mpi_mulm(rez_1, g, d, n);
		return;
	};

	gcry_mpi_t buff = gcry_mpi_new(0);
	gcry_mpi_t buff2 = gcry_mpi_new(0);
	
	//rez_x = (dg + ef + bdf) % n

	gcry_mpi_mulm(buff, d, g, n); // d * g
	gcry_mpi_mulm(buff2, e, f, n); // e * f
	gcry_mpi_addm(buff, buff, buff2, n); // d * g + e * f
	gcry_mpi_mulm(buff2, p->b, d, n); // b * d
	gcry_mpi_mulm(buff2, buff2, f, n); // b * d * f
	gcry_mpi_addm(rez_x, buff, buff2, n);

	//rez_1 = (eg + cdf) % n
	
	gcry_mpi_mulm(buff, e, g, n); // e * g
	gcry_mpi_mulm(buff2, p->c, d, n); // c * d
	gcry_mpi_mulm(buff2, buff2, f, n); // c * d * f
	gcry_mpi_addm(rez_1, buff, buff2, n);
};

//From proposion 3.2
//Calculate square of f(x) = (f*x + g)^2 mod (n, x^2 - b*x - c)

static void square_mod(gcry_mpi_t rez_x, gcry_mpi_t rez_1, const gcry_mpi_t f, const gcry_mpi_t g, const gcry_mpi_t n, struct params *p)
{
	if(gcry_mpi_cmp_ui(f, 0) == 0)
	{
		gcry_mpi_set_ui(rez_x, 0);
		gcry_mpi_mulm(rez_1, f, f, n);

		return;
	}	

	//res_x = (2fg + b(f^2)) % n
	gcry_mpi_t buff = gcry_mpi_new(0);

	gcry_mpi_mulm(buff, f, g, n); // f * g
	gcry_mpi_mulm(buff, buff, two, n); // 2 * f * g
	gcry_mpi_powm(rez_x, f, two, n); // f^2
	gcry_mpi_mulm(rez_x, rez_x, p->b, n); // b * f^2
	gcry_mpi_addm(rez_x, rez_x, buff, n);

	//rez_1 = (g^2 + c(f^2)) % n
	gcry_mpi_powm(buff, g, two, n); // g^2
	gcry_mpi_powm(rez_1, f, two, n); // f^2
       	gcry_mpi_mulm(rez_1, rez_1, p->c, n); // c * f^2
	gcry_mpi_addm(rez_1, buff, rez_1, n);
};
