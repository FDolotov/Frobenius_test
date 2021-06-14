#include "frobenius.h"
#include "primes.h"
#include <stdbool.h>

//Step 1: test n for divisibility by primes <= square root of num
bool step_1(const gcry_mpi_t n)
{
	printf("Step 1:\n");
	bool primality = false; //1 - prime, -1 - composite

	gcry_mpi_t buff = gcry_mpi_new(0);
	square_root(buff, n);

	gcry_mpi_t prime = gcry_mpi_new(0);
	gcry_mpi_set_ui(prime, prime_list[0]);

	gcry_mpi_t mod = gcry_mpi_new(0);

	for(int i = 1; gcry_mpi_cmp(n, prime) >= 0; i++)
	{
		if(i > 5133)
			break;
	
		if(gcry_mpi_cmp(prime, n) == 0)
		{
			primality = true;
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

		gcry_mpi_mod(mod, n, prime);

		if(gcry_mpi_cmp(zero, mod) == 0)
		{
			primality = false;
			break;
		}
		else
		{
			primality = true;
			gcry_mpi_set_ui(prime, prime_list[i]);
		}
	};

	gcry_mpi_release(buff);
	gcry_mpi_release(prime);
	gcry_mpi_release(mod);

	return primality;
};

//Step 2: test whether sqrt(n) is integer
bool step_2 (const gcry_mpi_t n)
{
	printf("Step 2:\n");

	bool primality;
	gcry_mpi_t buff = gcry_mpi_new(0);

	square_root(buff, n);
	gcry_mpi_mul(buff, buff, buff);
	
	if(gcry_mpi_cmp(buff, n) == 0)
		primality = false;
	else
		primality = true;

	gcry_mpi_release(buff);

	return primality;
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

//x * (fx + g) mod(n, x^2 - bx - c) = (fb + g)x + fc
static void mult_x_mod(gcry_mpi_t rez_x, gcry_mpi_t rez_1, const gcry_mpi_t f, const gcry_mpi_t g, const gcry_mpi_t n, struct params *p)
{
	gcry_mpi_t buff = gcry_mpi_new(0);
	gcry_mpi_t buff2 = gcry_mpi_new(0);

	gcry_mpi_set(buff, f);
	gcry_mpi_set(buff2, g);

	gcry_mpi_mulm(rez_1, p->c, f, n);

	gcry_mpi_mulm(rez_x, p->b, buff, n);
	gcry_mpi_addm(rez_x, rez_x, buff, n);
}

//From proposion 3.2
//Calculate square of f(x) = (f*x + g)^2 mod (n, x^2 - b*x - c)

static void square_mod(gcry_mpi_t rez_x, gcry_mpi_t rez_1, const gcry_mpi_t f, const gcry_mpi_t g, const gcry_mpi_t n, struct params *p)
{
	if(gcry_mpi_cmp(f, zero) == 0)
	{
		gcry_mpi_set(rez_x, zero);
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


//From proposion 3.3 and theorem 3.4
static void power_mod(gcry_mpi_t rez_x, gcry_mpi_t rez_1, const gcry_mpi_t power, const gcry_mpi_t n, struct params *p)
{
	gcry_mpi_t Aj = gcry_mpi_new(0);
	gcry_mpi_t Bj = gcry_mpi_new(0);
	gcry_mpi_t Cj = gcry_mpi_new(0);
	gcry_mpi_t buff = gcry_mpi_new(0);
	gcry_mpi_t buff2 = gcry_mpi_new(0);

	//start computing Aj, Bj & Cj
	//A1 = x^1 + (b - x)^1 = b
	gcry_mpi_set(Aj, p->b);

	//B1 = (x^1 - (b - x)^1)/(2x - b) = 1
	gcry_mpi_set(Bj, one);

	//C1 = c^1 = c
	gcry_mpi_set(Cj, p->c);

	int j = -1;

	for(int base = gcry_mpi_get_nbits(power) -1; base > 0; base--)
	{
		//Doubling
		//Compute B2j = Bj*Aj
		gcry_mpi_mulm(Bj, Bj, Aj, n);

		//Compute A2j = Aj^2 (-1)^j * 2Cj
		gcry_mpi_mulm(Aj, Aj, Aj, n);
		gcry_mpi_addm(buff, Cj, Cj, n);

		if(j == 1)
			gcry_mpi_sub(Aj, Aj, buff);
		else
			gcry_mpi_addm(Aj, Aj, buff, n);

		//Compute C2j = Cj^2
		gcry_mpi_mulm(Cj, Cj, Cj, n);

		j = 1;

		if(gcry_mpi_test_bit(power, base) == 1)
		{
			//Chain addition
			//Compute A(j+1) = 2^-1(AjA1 + (b^2 + 4c)BjB1)
			gcry_mpi_mulm(buff, p->b, p->b, n);
			gcry_mpi_mul_ui(buff2, p->c, 4);
			gcry_mpi_addm(buff, buff, buff2, n);

			gcry_mpi_mulm(buff, buff, Bj, n); //As B1 = 1, we can ignore it
			gcry_mpi_mulm(buff2, Aj, p->b, n); //A1 = b
			gcry_mpi_addm(buff, buff, buff2, n);
			gcry_mpi_div(buff, NULL, buff, two, 0);

			//Compute B(j+1) = 2^-1(AjB1 + A1Bj)
			gcry_mpi_mulm(buff2, p->b, Bj, n); // A1 = b
			gcry_mpi_addm(buff2, buff2, Aj, n); //Use Aj, not A(j+1), B1 = 1
			gcry_mpi_div(Bj, NULL, buff2, two, 0);

			gcry_mpi_set(Aj, buff); //Aj -> A(j+1)

			//Compute C(j+1) = CjC1
			gcry_mpi_mulm(Cj, Cj, p->c, n); //C1 = c

			j = -1;
		};

	};

	//Compute x^j = Bjx + 2^-1(Aj - bBj)
	gcry_mpi_set(rez_x, Bj);

	gcry_mpi_mulm(buff, p->b, Bj, n);
	gcry_mpi_set(rez_1, Aj);
	gcry_mpi_sub(rez_1, rez_1, buff);
	gcry_mpi_div(rez_1, NULL, rez_1, two, 0);

	gcry_mpi_release(Aj);
	gcry_mpi_release(Bj);
	gcry_mpi_release(Cj);
	gcry_mpi_release(buff);
	gcry_mpi_release(buff2);
};

static void sigma(gcry_mpi_t rez_x, gcry_mpi_t rez_1, const gcry_mpi_t f, const gcry_mpi_t g, const gcry_mpi_t n, struct params *p)
{
	gcry_mpi_sub(rez_x, n, f);
	
	gcry_mpi_mulm(rez_1, f, p->b, n);
	gcry_mpi_addm(rez_1, rez_1, g, n);
};

//set f*x + g mod (n, x^2 - bx - c) to multiplicative inverse
static int invert_mod(gcry_mpi_t rez_x, gcry_mpi_t rez_1, const gcry_mpi_t f, const gcry_mpi_t g, const gcry_mpi_t n, struct params *p)
{
	gcry_mpi_t buff = gcry_mpi_new(0);

	if(gcry_mpi_cmp(g, zero) == 0)
	{
		gcry_mpi_mulm(buff, f, p->c, n);
		gcry_mpi_invm(rez_x, buff, n);
		gcry_mpi_sub(rez_1, n, rez_x);
		gcry_mpi_mulm(rez_1, rez_1, p->b, n);
		gcry_mpi_release(buff);
		return 1;
	}

	//g^(-1)
	//if g != 0 is not invertable, gcd(f_1, n) is non-trivial, than n is composite
	if(!gcry_mpi_invm(buff, g, n))
	{
		gcry_mpi_release(buff);
		return 0;
	}

	// f * b * g + g^2 - f^2 * c
	gcry_mpi_mulm(buff, f, f, n); //f^2
	gcry_mpi_mulm(buff, buff, p->c, n); //f^2 *c

	gcry_mpi_t buff1 = gcry_mpi_new(0);

	gcry_mpi_mulm(buff1, f, p->b, n); //f *b

	gcry_mpi_t buff2 = gcry_mpi_new(0);

	gcry_mpi_mulm(buff1, buff1, g, n); //f * b * g
	gcry_mpi_mulm(buff2, g, g, n); //f_1^2
	gcry_mpi_addm(buff1, buff1, buff2, n); //f * b * g + g^2
	gcry_mpi_subm(buff1, buff1, buff, n); // - f^2 * c

	//f * b * g + g^2 - f^2 * c != 0
	//if its not invertable, n should be composite
	if(!gcry_mpi_invm(buff1, buff1, n))
	{
		gcry_mpi_release(buff);
		gcry_mpi_release(buff1);
		gcry_mpi_release(buff2);
		return 0;
	}

	//rez_x = -f/(f * b * g + g^2 - f^2 * c)
	gcry_mpi_mulm(rez_x, f, buff1, n);
	gcry_mpi_sub(rez_x, n, rez_x);

	//rez_1 = (f * b + g)/(f * b * g + g^2 - f^2 * c)
	gcry_mpi_mulm(rez_1, f, p->b, n);
	gcry_mpi_addm(rez_1, rez_1, g, n);
	gcry_mpi_mulm(rez_1, rez_1, buff1, n);

	gcry_mpi_release(buff);
	gcry_mpi_release(buff1);
	gcry_mpi_release(buff2);

	return 1;
};


bool step3_4_5(const gcry_mpi_t n, struct params *p)
{
	u_int64_t r;

	bool primality = false;
	gcry_mpi_t s = gcry_mpi_new(0);
	gcry_mpi_t t = gcry_mpi_new(0); //if 2^r*s = n +- 1, t = (s-1)/2

	//args of poly x^t mod (n, x^2 - bx - c)
	gcry_mpi_t x_t_x = gcry_mpi_new(0);
	gcry_mpi_t x_t_1 = gcry_mpi_new(0);

	//args of temporary polynom
	gcry_mpi_t f = gcry_mpi_new(0);
	gcry_mpi_t g = gcry_mpi_new(0);

	//Theorem 3.4: we should check if n = 1 mod 4 or n = 3 mod 4
	//Step 3: check if x^(n+1)/2 is in Z
	printf("Step 3:\n");
	
	gcry_mpi_t buff = gcry_mpi_new(0);
	gcry_mpi_set_ui(buff, 4);
	gcry_mpi_mod(buff, n, buff);
	

	bool n_1_mod_4 = gcry_mpi_cmp_ui(buff, 1) == 0;
	if(n_1_mod_4)
		gcry_mpi_sub(buff, n, one);
	else
		gcry_mpi_add(buff, n, one);

	split(&r, s, buff);
	gcry_mpi_div(t, NULL, s, two, 0);//t = (s-1)/2

	power_mod(x_t_x, x_t_1, t, n, p); //(x_t_x * x + x_t_1) = x^t mod(n, x^2 - bx - c)
	square_mod(f, g, x_t_x, x_t_1, n, p); //(x^t) ^2 = x^(s-1)
	mult_x_mod(f, g, f, g, n, p); //x * x^(s-1) = x^s

	//fx + g = x^s
	//Theorem 3.4:
	//To calculate x^(n-1)/2 or x^(n+1)/2 we need to square x^s r-1 times
	for(u_int64_t i = 0; i < r - 1; i++)
		square_mod(f, g, f, g, n, p);

	if(n_1_mod_4)
		//in this case fx + g = x^(n-1)/2
		//as we need co calculate x^(n+1)/2, we just multiply (fx + g) by x
		mult_x_mod(f, g, f, g, n, p);

	//Params to store x^(n+1)/2
	gcry_mpi_t power_f = gcry_mpi_new(0);
	gcry_mpi_t power_g = gcry_mpi_new(0);

	gcry_mpi_set(power_f, f);
	gcry_mpi_set(power_g, g);

	if(gcry_mpi_cmp(f, zero) == 0)
		primality = true;

	//Step 4: check if x^(n+1) == -c
	printf("Step 4:\n");
	gcry_mpi_mul(f, f, f);
	gcry_mpi_sub(buff, n, p->c);

	if(gcry_mpi_cmp(f, buff) != 0)
		primality = false;
	else
		primality = true;

	//Step 5:
	printf("Step 5:\n");
	gcry_mpi_mul(buff, n, n);
	gcry_mpi_sub(buff, buff, one);
	
	//calculate r, s, that 2^r*s == n^2 - 1
	split(&r, s, buff);
	if(n_1_mod_4)
	{
		sigma(f, g, x_t_x, x_t_1, n, p);
		mult_mod(f, g, f, g, x_t_x, x_t_1, n, p);
		mult_mod(f, g, f, g, power_f, power_g, n, p);
	}
	else
	{
		sigma(f, g, x_t_x, x_t_1, n, p);
		mult_mod(f, g, f, g, power_f, power_g, n, p);
		mult_x_mod(x_t_x, x_t_1, x_t_x, x_t_1, n, p);
		invert_mod(x_t_x, x_t_1, x_t_x, x_t_1, n, p);
		mult_mod(f, g, f, g, x_t_x, x_t_1, n, p);
	}

	gcry_mpi_sub(buff, n, one);

	if(gcry_mpi_cmp(f, zero) == 0 && gcry_mpi_cmp(g, one) == 0)
		primality = true;

	for(u_int64_t i = 0; i < r - 1; i++)
	{
		if(gcry_mpi_cmp(f, zero) == 0 && gcry_mpi_cmp(g, buff) == 0)
			primality = true;

		square_mod(f, g, f, g, n, p);	
	};
	
	gcry_mpi_release(t);
	gcry_mpi_release(s);
	gcry_mpi_release(f);
	gcry_mpi_release(g);
	gcry_mpi_release(x_t_x);
	gcry_mpi_release(x_t_1);
	gcry_mpi_release(power_f);
	gcry_mpi_release(power_g);
	gcry_mpi_release(buff);

	return primality;
};

extern void QFT(const gcry_mpi_t n, struct params *p)
{
	bool prime = step_1(n);
	prime = step_2(n);
	if(!prime)
	{
		printf("Number is composite\n");
		return;
	}

	if(step3_4_5(n, p))
		printf("Number is probably prime\n");
	else
		printf("number is composite\n");
};

extern void RQFT(const gcry_mpi_t n)
{
	struct params parameters;
	set_params(&parameters, n);
	
	QFT(n, &parameters);

	release_params(&parameters);
};
