#include "stdafx.h"
#include <float.h>
#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <intrin.h>

constexpr uint32_t INDEX_BITS = 7;
constexpr int MAX_NEWTON_RAPHSON_ITER = 2;

double min_error[MAX_NEWTON_RAPHSON_ITER] = {};
double max_error[MAX_NEWTON_RAPHSON_ITER] = {};
uint64_t worst_divisor_min[MAX_NEWTON_RAPHSON_ITER] = {};
uint64_t worst_divisor_max[MAX_NEWTON_RAPHSON_ITER] = {};
uint64_t smallest_divisor[MAX_NEWTON_RAPHSON_ITER] = {};

#pragma section("udiv128", read, execute)
__declspec(allocate("udiv128"))
static const unsigned char MachineCode[] =
{
	0x48, 0x89, 0xD0, // mov rax,rdx
	0x48, 0x89, 0xCA, // mov rdx,rcx
	0x49, 0xF7, 0xF0, // div r8
	0x49, 0x89, 0x11, // mov [r9],rdx
	0xC3,             // ret
};

uint64_t(*udiv128)(uint64_t numhi, uint64_t numlo, uint64_t den, uint64_t* rem) = (uint64_t(*)(uint64_t, uint64_t, uint64_t, uint64_t*))((const unsigned char*)MachineCode);

struct SReciprocal
{
	uint32_t C0;
	uint16_t C1;
	uint16_t C2;
} reciprocals[1UL << INDEX_BITS];
static_assert(sizeof(SReciprocal) == sizeof(uint64_t), "SReciprocal size is incorrect");

void first_order_approximation_init()
{
	_control87(RC_NEAR, MCW_RC);

	for (uint64_t i = 0; i < (1UL << INDEX_BITS); ++i)
	{
		const uint64_t divisor1 = (1ULL << 31) + (i << (31 - INDEX_BITS)) + 1;
		const uint64_t divisor2 = (1ULL << 31) + ((i + 1) << (31 - INDEX_BITS)) - 1;
		const uint64_t divisor_mid = (divisor1 + divisor2) >> 1;

		uint64_t rem;

		uint64_t r1 = udiv128(1, 0, divisor1, &rem);
		if (rem * 2 > divisor1) ++r1;

		uint64_t r2 = udiv128(1, 0, divisor2, &rem);
		if (rem * 2 > divisor2) ++r1;

		uint64_t r_mid = udiv128(1, 0, divisor_mid, &rem);
		if (rem * 2 > divisor_mid) ++r_mid;

		uint64_t c1 = udiv128(r1 - r_mid, 0, (divisor_mid - divisor1) << 2, &rem);
		if (rem * 2 > ((divisor_mid - divisor1) << 2)) ++c1;

		uint64_t c2 = udiv128(r_mid - r2, 0, (divisor2 - divisor_mid) << 2, &rem);
		if (rem * 2 > ((divisor2 - divisor_mid) << 2)) ++c2;

		c1 = (c1 >> 48) + ((c1 & (1ULL << 47)) ? 1 : 0);
		c2 = (c2 >> 48) + ((c2 & (1ULL << 47)) ? 1 : 0);

		reciprocals[i].C0 = static_cast<uint32_t>(r_mid);
		reciprocals[i].C1 = static_cast<uint16_t>(c1);
		reciprocals[i].C2 = static_cast<uint16_t>(c2);
	}

	std::cout << "Lookup table size is " << (1UL << INDEX_BITS) << " entries" << std::endl;
	std::cout << "Calculating piecewise linear approximation of 1/x with 64-bit precision (ULP = 2^-64)" << std::endl;
	for (uint64_t i = 0; i < (1UL << INDEX_BITS); ++i)
	{
		const uint64_t divisor1 = 0x80000000ULL + (i << (31 - INDEX_BITS)) + 1;
		const uint64_t divisor2 = 0x80000000ULL + ((i + 1) << (31 - INDEX_BITS)) - 1;

		double min_error = 0.0;
		double max_error = 0.0;
		for (uint64_t a = divisor1; a <= divisor2; ++a)
		{
			uint64_t r;

			double exact_r = static_cast<double>(udiv128(1, 0, a, &r));
			exact_r += static_cast<double>(r) / a;

			const int32_t index2 = static_cast<int32_t>(a & ((1UL << (31 - INDEX_BITS)) - 1)) - (1 << (30 - INDEX_BITS)) + 1;
			const int64_t r1 = reciprocals[i].C0 | 0x100000000LL;
			int64_t r2 = ((index2 < 0) ? reciprocals[i].C1 : reciprocals[i].C2);
			r = r1 - ((r2 * index2) >> 14);

			const double err = r - exact_r;
			if (err < min_error)
			{
				min_error = err;
			}
			if (err > max_error)
			{
				max_error = err;
			}
		}
		const int32_t fix = static_cast<int32_t>(round((max_error + min_error) * 0.5));
		reciprocals[i].C0 -= fix;
		min_error -= fix;
		max_error -= fix;
		std::cout << "err[" << i << "] = (" << min_error << ", " << max_error << ") ULP" << std::endl;
	}
	std::cout << std::endl;
}

void first_order_approximation_check()
{
	std::cout << "Calculating maximum error after up to " << MAX_NEWTON_RAPHSON_ITER << " Newton-Raphson iteration(s)" << std::endl << std::endl;

	uint64_t prev_index = 0;
	for (uint64_t a = 0xFFFFFFFFUL; a >= 3; --a)
	{
		if ((a >> 20) != prev_index)
		{
			std::cout << std::setw(4) << (a >> 20) << '\r';
			prev_index = (a >> 20);
		}

		uint64_t r;

		uint64_t exact_r = udiv128(1, 0, a, &r);
		double exact_r_fraction = static_cast<double>(r) / a;

		unsigned long shift;
		_BitScanReverse64(&shift, a);
		shift = 31 - shift;
		const uint64_t a1 = a << shift;

		const uint32_t index1 = (static_cast<uint32_t>(a1 << 1) >> (32 - INDEX_BITS));
		const int32_t index2 = static_cast<int32_t>(a1 & ((1UL << (31 - INDEX_BITS)) - 1)) - (1 << (30 - INDEX_BITS)) + 1;
		const int64_t r1 = reciprocals[index1].C0 | 0x100000000LL;
		int64_t r2 = ((index2 < 0) ? reciprocals[index1].C1 : reciprocals[index1].C2);
		r = (r1 - ((r2 * index2) >> 14)) << shift;

		for (int i = 0; i < MAX_NEWTON_RAPHSON_ITER; ++i)
		{
			uint64_t lo, hi, hi2;
			lo = _umul128(r, a, &hi);

			// Pre-condition: r*a ~ 2^64
			// Instead of doing r*(2 - a*r) we do
			// r*(2 + a/2^65 - a*r) =
			// r*(2 - a*r) + r*a/2^65 ~=
			// r*(2 - a*r) + 0.5
			// We only have to do truncation after that and get almost perfect rounding
			// Truncation allows to skip computing low 64 bits of the last 128-bit multiplication
			_subborrow_u64(_subborrow_u64(0, a >> 1, lo, &lo), 2, hi, &hi);
			_umul128(r, lo, &hi2);
			r = hi2 + (hi ? r : 0);

			const double err = static_cast<int64_t>(r - exact_r) - exact_r_fraction;
			if (err < min_error[i])
			{
				min_error[i] = err;
				worst_divisor_min[i] = a;
			}
			if (err > max_error[i])
			{
				max_error[i] = err;
				worst_divisor_max[i] = a;
			}
			if (((min_error[i] <= -0.999999) || (max_error[i] >= 0.999999)) && (smallest_divisor[i] == 0))
			{
				smallest_divisor[i] = a;
				std::cout << "Smallest incorrect divisor after " << (i + 1) << " Newton-Raphson iteration(s):  " << smallest_divisor[i] <<
							", err = (" << min_error[i] << ", " << max_error[i] << ") ULP" << std::endl << std::endl;
			}
		}
	}

	std::cout.precision(15);
	for (int i = 0; i < MAX_NEWTON_RAPHSON_ITER; ++i)
	{
		std::cout << (i + 1) << " Newton-Raphson iteration(s):" << std::endl;
		std::cout << "err = (" << min_error[i] << ", " << max_error[i] << ") ULP" << std::endl;
		std::cout << "Worst divisor (undershoot): " << worst_divisor_min[i] << std::endl;
		std::cout << "Worst divisor (overshoot):  " << worst_divisor_max[i] << std::endl << std::endl;
	}
}

extern uint64_t get_reciprocal(uint32_t a);
extern void fast_div(uint64_t a, uint32_t b, uint64_t &q, uint64_t &r);

void check_generated_code()
{
	std::cout << "Checking generated code" << std::endl << std::endl;

	double min_err = 0.0;
	double max_err = 0.0;
	uint64_t min_err_divisor = 0;
	uint64_t max_err_divisor = 0;
	uint64_t prev_index = 0;
	for (uint32_t a = 0xFFFFFFFFUL; a >= 0x80000001; a -= 2)
	{
		if ((a >> 20) != prev_index)
		{
			std::cout << std::setw(4) << (a >> 20) << '\r';
			prev_index = (a >> 20);
		}

		uint64_t r;
		uint64_t exact_r = udiv128(1, 0, a, &r);
		double exact_r_fraction = static_cast<double>(r) / a;
		r = get_reciprocal(a);

		const double err = static_cast<int64_t>(r - exact_r) - exact_r_fraction;
		if (err < min_err)
		{
			min_err = err;
			min_err_divisor = a;
		}
		if (err > max_err)
		{
			max_err = err;
			max_err_divisor = a;
		}
	}

	std::cout << "err = (" << min_err << ", " << max_err << ") ULP" << std::endl;
	std::cout << "Worst divisor (undershoot): " << min_err_divisor << std::endl;
	std::cout << "Worst divisor (overshoot):  " << max_err_divisor << std::endl << std::endl;
}

void check_and_benchmark_generated_code()
{
	std::vector<std::pair<uint64_t, uint64_t>> test_numbers;
	test_numbers.resize(100000000);
	std::mt19937_64 rnd;
	rnd.seed(123);

	std::cout << "Generating " << test_numbers.size() << " random numbers" << std::endl;
	for (std::pair<uint64_t, uint64_t>& n : test_numbers)
	{
		n.first = rnd();
		n.second = rnd();
	}

	std::cout << "Comparing DIV instruction and fast_div" << std::endl;
	for (const std::pair<uint64_t, uint64_t>& n : test_numbers)
	{
		const uint64_t dividend = n.first;
		const uint32_t divisor1 = static_cast<uint32_t>(n.second) | 0x80000001UL;
		const uint32_t divisor2 = static_cast<uint32_t>(n.second >> 32) | 0x80000001UL;

		uint64_t q1, r1;
		fast_div(dividend, divisor1, q1, r1);
		const uint64_t q1_check = dividend / divisor1;
		const uint64_t r1_check = dividend % divisor1;
		if ((q1 != q1_check) || (r1 != r1_check))
		{
			std::cout << "Incorrect result for " << dividend << '/' << divisor1 << std::endl;
			return;
		}

		uint64_t q2, r2;
		fast_div(dividend, divisor2, q2, r2);
		const uint64_t q2_check = dividend / divisor2;
		const uint64_t r2_check = dividend % divisor2;
		if ((q2 != q2_check) || (r2 != r2_check))
		{
			std::cout << "Incorrect result for " << dividend << '/' << divisor2 << std::endl;
			return;
		}
	}

	std::cout << "Benchmarking DIV instruction and fast_div" << std::endl;

	using namespace std::chrono;

	uint64_t sum0 = 0;
	uint64_t sum1 = 0;
	uint64_t sum2 = 0;

	const uint64_t t0 = time_point_cast<milliseconds>(high_resolution_clock::now()).time_since_epoch().count();

	for (const std::pair<uint64_t, uint64_t>& n : test_numbers)
	{
		sum0 += n.first + n.second;
	}

	const uint64_t t1 = time_point_cast<milliseconds>(high_resolution_clock::now()).time_since_epoch().count();

	for (const std::pair<uint64_t, uint64_t>& n : test_numbers)
	{
		const uint64_t dividend = n.first;
		const uint32_t divisor1 = static_cast<uint32_t>(n.second) | 0x80000001UL;
		const uint32_t divisor2 = static_cast<uint32_t>(n.second >> 32) | 0x80000001UL;
		const uint64_t q1 = dividend / divisor1;
		const uint64_t r1 = dividend % divisor1;
		const uint64_t q2 = dividend / divisor2;
		const uint64_t r2 = dividend % divisor2;
		sum1 += q1 + r1 + q2 + r2;
	}

	const uint64_t t2 = time_point_cast<milliseconds>(high_resolution_clock::now()).time_since_epoch().count();

	for (const std::pair<uint64_t, uint64_t>& n : test_numbers)
	{
		const uint64_t dividend = n.first;
		const uint32_t divisor1 = static_cast<uint32_t>(n.second) | 0x80000001UL;
		const uint32_t divisor2 = static_cast<uint32_t>(n.second >> 32) | 0x80000001UL;
		uint64_t q1, r1;
		fast_div(dividend, divisor1, q1, r1);
		uint64_t q2, r2;
		fast_div(dividend, divisor2, q2, r2);
		sum2 += q1 + r1 + q2 + r2;
	}

	const uint64_t t3 = time_point_cast<milliseconds>(high_resolution_clock::now()).time_since_epoch().count();

	if (sum1 != sum2)
	{
		std::cout << "Incorrect overall result" << std::endl;
		return;
	}

	std::cout << "Empty loop: " << (t1 - t0) << " ms, sum0 = " << sum0 << std::endl;
	std::cout << "DIV instruction: " << (t2 - t1 - (t1 - t0)) << " ms" << std::endl;
	std::cout << "fast_div: " << (t3 - t2 - (t1 - t0)) << " ms" << std::endl;
}


int main()
{
	std::cout.precision(15);

	check_generated_code();
	check_and_benchmark_generated_code();
	first_order_approximation_init();
	first_order_approximation_check();

	std::ofstream f("reciprocal64.cpp");

	f << "#include <stdint.h>\n";
	f << "#include <intrin.h>\n\n";
	f << "static const uint32_t R[" << (1 << (INDEX_BITS + 1)) << "] = {\n";
	f << std::hex;
	for (int i = 0; i < (1 << (INDEX_BITS + 1)); i += 8)
	{
		f << '\t';
		for (int j = 0; j < 8; ++j)
		{
			f << "0x" << reinterpret_cast<uint32_t*>(reciprocals)[i + j] << "u,";
		}
		f << '\n';
	}
	f << std::dec;
	f << "};\n\n";
	f << "struct SReciprocal\n";
	f << "{\n";
	f << "	uint32_t C0;\n";
	f << "	uint16_t C1;\n";
	f << "	uint16_t C2;\n";
	f << "};\n\n";
	f << "#define RCP ((SReciprocal*)R)\n\n";
	f << "uint64_t get_reciprocal(uint32_t a)\n";
	f << "{\n";
	f << "	unsigned long shift;\n";
	f << "	_BitScanReverse64(&shift, a);\n";
	f << "	shift = 31 - shift;\n";
	f << "	const uint64_t a1 = a << shift;\n";
	f << "\n";
	f << "	const uint32_t index1 = (static_cast<uint32_t>(a1 << 1) >> " << (32 - INDEX_BITS) << ");\n";
	f << "	const int32_t index2 = " << "static_cast<int32_t>(a1 & " << ((1UL << (31 - INDEX_BITS)) - 1) << ") - " << ((1 << (30 - INDEX_BITS)) - 1) << ";\n";
	f << "	const int64_t r1 = RCP[index1].C0 | 0x100000000LL;\n";
	f << "	const int32_t r2_1 = RCP[index1].C1;\n";
	f << "	const int32_t r2_2 = RCP[index1].C2;\n";
	f << "	const int32_t r2 = (index2 < 0) ? r2_1 : r2_2;\n";
	f << "	uint64_t r = (r1 - ((static_cast<int64_t>(r2) * index2) >> 14)) << shift;\n";
	f << "\n";
	f << "	uint64_t lo, hi, hi2;\n";
	f << "\n";
	f << "	lo = _umul128(r, a, &hi);\n";
	f << "	_subborrow_u64(_subborrow_u64(0, a >> 1, lo, &lo), 2, hi, &hi);\n";
	f << "	_umul128(r, lo, &hi2);\n";
	f << "	r = hi2 + (hi ? r : 0);\n";
	f << "\n";
	f << "	if (a <= " << smallest_divisor[0] << ")\n";
	f << "	{\n";
	f << "		lo = _umul128(r, a, &hi);\n";
	f << "		_subborrow_u64(_subborrow_u64(0, a >> 1, lo, &lo), 2, hi, &hi);\n";
	f << "		_umul128(r, lo, &hi2);\n";
	f << "		r = hi2 + (hi ? r : 0);\n";
	f << "	}\n";
	f << "\n";
	f << "	return r;\n";
	f << "}\n\n";
	f << "void fast_div(uint64_t a, uint32_t b, uint64_t &q, uint64_t &r)\n";
	f << "{\n";
	f << "	_umul128(a, get_reciprocal(b), &q);\n";
	f << "	int64_t tmp = a - q * b;\n";
	f << "	if (tmp < 0)\n";
	f << "	{\n";
	f << "		--q;\n";
	f << "		tmp += b;\n";
	f << "	}\n";
	f << "	else if (tmp >= b)\n";
	f << "	{\n";
	f << "		++q;\n";
	f << "		tmp -= b;\n";
	f << "	}\n";
	f << "	r = static_cast<uint64_t>(tmp);\n";
	f << "}\n";

	f.close();

	return 0;
}
