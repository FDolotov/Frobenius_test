#include "help.h"
#include <assert.h>
#include <string.h>

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

void square_root(gcry_mpi_t res, const gcry_mpi_t n)
{
	int bit = gcry_mpi_get_nbits(n);
	bit /= 2;
	
	gcry_mpi_t buff = gcry_mpi_new(0);
	gcry_mpi_t bit_t = gcry_mpi_new(0);
	gcry_mpi_t current_t = gcry_mpi_new(0);
	gcry_mpi_set_ui(bit_t, bit);
	gcry_mpi_set(current_t, bit_t);

	while(1)
	{
		gcry_mpi_div(buff, NULL, n, bit_t, 0);
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
	assert(gcry_mpi_cmp(buff, one) == 0);

	gcry_mpi_mod(n, n, k);
    
	int t = 1;

	while (gcry_mpi_cmp(n, zero) != 0) 
	{
		gcry_mpi_mod(buff2, n, two);
        	while (gcry_mpi_cmp(buff2, zero) == 0) 
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
};

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

void split(u_int64_t *s, gcry_mpi_t d, const gcry_mpi_t num)
{
	u_int64_t bit = gcry_mpi_get_nbits(num);
	
	for(u_int64_t i = 1; i <= bit; i++)
	{
		if(gcry_mpi_test_bit(num, i))
		{
			*s = i;
			break;
		}
	};

	gcry_mpi_mul_2exp(d, one, *s);
	gcry_mpi_div(d, NULL, num, d, 0);
};

//Base could be 2 or 16
int number_length(const gcry_mpi_t n, const int base)
{	
	int l = 0;
	gcry_mpi_t del = gcry_mpi_new(0);

	gcry_mpi_set_ui(del, 16);

	gcry_mpi_t buff = gcry_mpi_new(0);
	gcry_mpi_set(buff, n);

	while(gcry_mpi_cmp(buff, zero) != 0)
	{
		gcry_mpi_div(buff, NULL, buff, del, 0);
		//gcry_mpi_rshift(buff, buff, bit);
		l++;
	}

	gcry_mpi_release(buff);
	return l;
};

void hex_to_bin(gcry_mpi_t rez, const gcry_mpi_t n)
{
	gcry_mpi_t number = gcry_mpi_new(0);
	gcry_mpi_set(number, n);
	
	gcry_mpi_t buff = gcry_mpi_new(0);
	
	int len = number_length(n, 16);

	char num[len*4];
	strcpy(num, "0");

	char bin[16][5] = {"0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111", "1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111"};

	gcry_mpi_t buff2 = gcry_mpi_new(0);

	int k;
	for(int i = len - 1; i >= 0; i--)
	{
		gcry_mpi_rshift(buff, number, 4 * i);
		gcry_mpi_set(buff2, buff);
		for(int m = 0; m < i; m++)
			gcry_mpi_mul_ui(buff, buff, 16);

		gcry_mpi_sub(number, number, buff);

		k = 0;
		while(gcry_mpi_cmp_ui(buff2, k) != 0) 
			k++;

		strcat(num, bin + k);
	};

	gcry_mpi_t buff3 = gcry_mpi_new(0);
	size_t scanned;
	gcry_mpi_scan(&buff3, GCRYMPI_FMT_HEX, num, 0, &scanned);
	gcry_mpi_set(rez, buff3);

	
	gcry_mpi_release(buff);
	gcry_mpi_release(buff2);
	gcry_mpi_release(buff3);
	gcry_mpi_release(number);
};
