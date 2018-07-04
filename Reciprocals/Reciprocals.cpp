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

// 31st bit is always 1 (divisor is shifted left until this bit becomes 1)
// index1: bits 30-24, unsigned, range=[0,127]
// index2: bits 23-0, signed, range=[-8388608,8388607]
constexpr uint32_t INDEX1_BITS = 7;
constexpr uint32_t INDEX2_BITS = 31 - INDEX1_BITS;

// 1 ULP = 2^-64
constexpr int C0_BITS = 37; // 38 bits (37 bits stored) = 0.015625 ULP precision
constexpr int C1_BITS = 27; // 27 bits = 0.015625 ULP precision for C1*(x-x0) if |x-x0| <= 4194304
constexpr int C2_BITS = 20; // 20 bits = 0.015625 ULP precision for C2*(x-x0)^2 if |x-x0| <= 4194304
constexpr int C3_BITS = 12; // 12 bits = 0.0078125 ULP precision for C3*(x-x0)^3 if |x-x0| <= 4194304

// Total error = error from coefficients + error from the skipped remainder of Taylor series + error from rounding
//
// Error from coefficients <= 0.015625 * 3 + 0.0078125 = 0.0546875 ULP
// Error from the skipped remainder of Taylor series <= 2^64/((2^31)^5) * (2^22)^4 = 2^(64-155+88) = 2^-3 = 0.125 ULP
// Error from rounding <= 0.5 ULP
//
// Upper bound for total error: 0.0546875 + 0.125 + 0.5 = 0.6796875 ULP
// Actual bounds for total error: [-0.601, 0.602] ULP

constexpr int MAX_NEWTON_RAPHSON_ITER = 1;

double min_error[MAX_NEWTON_RAPHSON_ITER + 1] = {};
double max_error[MAX_NEWTON_RAPHSON_ITER + 1] = {};
uint64_t worst_divisor_min[MAX_NEWTON_RAPHSON_ITER + 1] = {};
uint64_t worst_divisor_max[MAX_NEWTON_RAPHSON_ITER + 1] = {};
uint64_t smallest_divisor[MAX_NEWTON_RAPHSON_ITER + 1] = {};

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

#pragma pack(push, 1)
struct SReciprocal
{
	uint64_t C0 : C0_BITS; 
	uint64_t C1 : C1_BITS;
	uint32_t C2 : C2_BITS;
	uint32_t C3 : C3_BITS;
} reciprocals[1UL << INDEX1_BITS];
#pragma pack(pop)

static_assert(sizeof(SReciprocal) == sizeof(uint32_t) * 3, "SReciprocal size is incorrect");

void cubic_approximation_init()
{
	_control87(RC_NEAR, MCW_RC);

	for (uint64_t i = 0; i < (1UL << INDEX1_BITS); ++i)
	{
		const uint64_t divisor1 = (1ULL << 31) + (i << (31 - INDEX1_BITS)) + 1;
		const uint64_t divisor2 = (1ULL << 31) + ((i + 1) << (31 - INDEX1_BITS)) - 1;
		const uint64_t divisor_mid = (divisor1 + divisor2) >> 1;

		uint64_t rem;

		uint64_t c0 = udiv128(1ULL << (C0_BITS - 32), 0, divisor_mid, &rem);
		if (rem > divisor_mid - rem) ++c0;

		const uint64_t divisor_mid2 = uint64_t(divisor_mid) * divisor_mid;
		const double x1 = (18446744073709551616.0 / 4.0 * (1ULL << C1_BITS)) / divisor_mid2;

		uint64_t divisor_mid3_hi;
		const uint64_t divisor_mid3_lo = _umul128(divisor_mid2, divisor_mid, &divisor_mid3_hi);
		const double divisor_mid3 = divisor_mid3_hi * 18446744073709551616.0 + divisor_mid3_lo;
		const double x2 = (18446744073709551616.0 * 4294967296.0 / 8.0 * (1ULL << C2_BITS)) / divisor_mid3;

		const double divisor_mid4 = divisor_mid3 * divisor_mid;
		const double x3 = (18446744073709551616.0 * 18446744073709551616.0 / 16.0 * (1ULL << C3_BITS)) / divisor_mid4;

		reciprocals[i].C0 = c0;
		reciprocals[i].C1 = static_cast<uint64_t>(round(x1));
		reciprocals[i].C2 = static_cast<uint32_t>(round(x2));
		reciprocals[i].C3 = static_cast<uint32_t>(round(x3));
	}

	std::cout << "Lookup table size is " << (1UL << INDEX1_BITS) << " entries" << std::endl;
	std::cout << "Calculating piecewise cubic polynomial approximation of 1/x with 64-bit precision (ULP = 2^-64)" << std::endl;
	std::cout << std::setprecision(3);
	double min_error_global = 0.0;
	double max_error_global = 0.0;
	for (uint64_t i = 0; i < (1UL << INDEX1_BITS); ++i)
	{
		const uint64_t divisor1 = 0x80000000ULL + (i << (31 - INDEX1_BITS));
		const uint64_t divisor2 = 0x80000000ULL + ((i + 1) << (31 - INDEX1_BITS)) - 1;

		double min_error = 0.0;
		double max_error = 0.0;
		for (uint64_t a = divisor1; a <= divisor2; ++a)
		{
			uint64_t rem;
			int64_t exact_r = udiv128(1, 0, a, &rem);
			double exact_r_fract = static_cast<double>(rem) / a;

			const int32_t index2 = static_cast<int32_t>(a & ((1UL << INDEX2_BITS) - 1)) - (1 << (INDEX2_BITS - 1));
			const int64_t r0 = int64_t(reciprocals[i].C0 | (1LL << C0_BITS)) << 22;
			const int64_t r1 = int64_t(reciprocals[i].C1) << 12;
			const int64_t r2 = (int64_t(reciprocals[i].C2) << 19) + (reciprocals[i].C2 << 3);
			const int32_t r3 = reciprocals[i].C3;

#if 0
			double r = -r3;

			r *= (index2 >> 4);
			if (fabs(r) >= 2147483648.0) __debugbreak();

			r = (r + r2) * index2;
			if (fabs(r) >= 9223372036854775808.0) __debugbreak();
			r /= static_cast<double>(1LL << 31);

			r = (r - r1) * index2;
			if (fabs(r) >= 9223372036854775808.0) __debugbreak();
			r /= static_cast<double>(1LL << 10);

			r = r + r0;
			if (fabs(r) >= 9223372036854775808.0) __debugbreak();
			r /= static_cast<double>(1LL << 27);
#else
			int64_t r;
			r = static_cast<int32_t>(-r3 * (index2 >> 4)); // 32x32 -> 32 bit mul
			r = ((r + r2) * index2) >> 31; // 64x32 -> 64 bit mul
			r = ((r - r1) * index2) >> 10; // 64x32 -> 64 bit mul
			r = (r + r0) >> 27;
#endif

			const double err = r - exact_r - exact_r_fract;
			if (err < min_error)
			{
				min_error = err;
			}
			if (err > max_error)
			{
				max_error = err;
			}
		}
		const int32_t fix = static_cast<int32_t>(round((max_error + min_error) * 0.5 * (1LL << (C0_BITS - 32))));
		reciprocals[i].C0 -= fix;
		min_error -= static_cast<double>(fix) / (1LL << (C0_BITS - 32));
		max_error -= static_cast<double>(fix) / (1LL << (C0_BITS - 32));
		std::cout << "err[" << std::setw(4) << i << "] = (" << std::setw(5) << min_error << ", " << std::setw(5) << max_error << ") ULP" << std::endl;
		if (min_error < min_error_global)
		{
			min_error_global = min_error;
		}
		if (max_error > max_error_global)
		{
			max_error_global = max_error;
		}
	}
	std::cout << "Global error = (" << std::setw(5) << min_error_global << ", " << std::setw(5) << max_error_global << ") ULP" << std::endl;
	std::cout << std::endl;
}

void cubic_approximation_check()
{
	std::cout << "Calculating maximum error after up to " << MAX_NEWTON_RAPHSON_ITER << " Newton-Raphson iteration(s)" << std::endl << std::endl;

	uint64_t prev_index = 0;
	for (uint32_t a = 0xFFFFFFFFUL; a >= 3; --a)
	{
		if ((a >> 20) != prev_index)
		{
			std::cout << std::setw(4) << (a >> 20) << '\r';
			prev_index = (a >> 20);
		}

		uint64_t rem;
		int64_t exact_r = udiv128(1, 0, a, &rem);
		double exact_r_fraction = static_cast<double>(rem) / a;

		unsigned long shift;
		_BitScanReverse(&shift, a);
		shift = 31 - shift;
		const uint32_t a1 = a << shift;

		const uint32_t index1 = ((a1 << 1) >> (INDEX2_BITS + 1));
		const int32_t index2 = static_cast<int32_t>(a1 & ((1UL << INDEX2_BITS) - 1)) - (1 << (INDEX2_BITS - 1));
		const int64_t r0 = int64_t(reciprocals[index1].C0 | (1LL << C0_BITS)) << 22;
		const int64_t r1 = int64_t(reciprocals[index1].C1) << 12;
		const int64_t r2 = (int64_t(reciprocals[index1].C2) << 19) + (reciprocals[index1].C2 << 3);
		const int32_t r3 = reciprocals[index1].C3;

		int64_t rr;
		rr = static_cast<int32_t>(-r3 * (index2 >> 4)); // 32x32 -> 32 bit mul
		rr = ((rr + r2) * index2) >> 31; // 64x32 -> 64 bit mul
		rr = ((rr - r1) * index2) >> 10; // 64x32 -> 64 bit mul
		uint64_t r = static_cast<uint64_t>((rr + r0) >> 27) << shift;

		for (int i = 0; i <= MAX_NEWTON_RAPHSON_ITER; ++i)
		{
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
				std::cout << "Smallest incorrect divisor after " << i << " Newton-Raphson iteration(s):  " << smallest_divisor[i] <<
					", err = (" << min_error[i] << ", " << max_error[i] << ") ULP" << std::endl;
			}

			if (i == MAX_NEWTON_RAPHSON_ITER)
			{
				break;
			}

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
		}
	}
	std::cout << "    " << std::endl;

	std::cout.precision(15);
	for (int i = 0; i <= MAX_NEWTON_RAPHSON_ITER; ++i)
	{
		std::cout << i << " Newton-Raphson iteration(s):" << std::endl;
		std::cout << "err = (" << min_error[i] << ", " << max_error[i] << ") ULP" << std::endl;
		std::cout << "Worst divisor (undershoot): " << worst_divisor_min[i] << std::endl;
		std::cout << "Worst divisor (overshoot):  " << worst_divisor_max[i] << std::endl << std::endl;
	}
}

extern uint64_t get_reciprocal(uint32_t a);
extern void fast_div(uint64_t a, uint32_t b, uint64_t &q, uint64_t &r);

constexpr uint32_t DIVISOR_OR_MASK = 0x80000001UL;

void check_generated_code()
{
	std::cout << "Checking generated code" << std::endl << std::endl;

	double min_err = 0.0;
	double max_err = 0.0;
	uint64_t min_err_divisor = 0;
	uint64_t max_err_divisor = 0;
	uint64_t prev_index = 0;
	for (uint32_t a = 0xFFFFFFFFUL; a >= DIVISOR_OR_MASK; --a)
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
	test_numbers.reserve(100000000);
	std::mt19937_64 rnd;
	rnd.seed(123);

	std::cout << "Generating " << test_numbers.capacity() << " random numbers" << std::endl;
	for (size_t i = 0, n = test_numbers.capacity(); i < n; ++i)
	{
		test_numbers.emplace_back(rnd(), rnd());
	}

	std::cout << "Comparing DIV instruction and fast_div" << std::endl;
	for (const std::pair<uint64_t, uint64_t>& n : test_numbers)
	{
		const uint64_t dividend = n.first;
		const uint32_t divisor1 = static_cast<uint32_t>(n.second) | DIVISOR_OR_MASK;
		const uint32_t divisor2 = static_cast<uint32_t>(n.second >> 32) | DIVISOR_OR_MASK;

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
		const uint32_t divisor1 = static_cast<uint32_t>(n.second) | DIVISOR_OR_MASK;
		const uint32_t divisor2 = static_cast<uint32_t>(n.second >> 32) | DIVISOR_OR_MASK;
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
		const uint32_t divisor1 = static_cast<uint32_t>(n.second) | DIVISOR_OR_MASK;
		const uint32_t divisor2 = static_cast<uint32_t>(n.second >> 32) | DIVISOR_OR_MASK;
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
	check_generated_code();
	check_and_benchmark_generated_code();
	cubic_approximation_init();
	cubic_approximation_check();

	std::ofstream f("reciprocal64.cpp");

	const size_t R_size = (1 << INDEX1_BITS) * sizeof(SReciprocal) / sizeof(uint32_t);

	f << "#include <stdint.h>\n";
	f << "#include <intrin.h>\n\n";
	f << "// Set it to 1 to correctly process divisors < 2^31\n";
	f << "#define SMALL_DIVISORS_ENABLED 0\n\n";
	f << "static const uint32_t R[" << R_size << "] = {\n";
	f << std::hex;
	for (size_t i = 0; i < R_size; i += 8)
	{
		f << '\t';
		for (size_t j = 0; j < 8; ++j)
		{
			f << "0x" << reinterpret_cast<uint32_t*>(reciprocals)[i + j] << "u,";
		}
		f << '\n';
	}
	f << std::dec;
	f << "};\n\n";
	f << "#pragma pack(push, 1)\n";
	f << "struct SReciprocal\n";
	f << "{\n";
	f << "	uint64_t C0 : " << C0_BITS << ";\n";
	f << "	uint64_t C1 : " << C1_BITS << ";\n";
	f << "	uint32_t C2 : " << C2_BITS << ";\n";
	f << "	uint32_t C3 : " << C3_BITS << ";\n";
	f << "};\n";
	f << "#pragma pack(pop)\n\n";
	f << "#define RCP ((SReciprocal*)R)\n\n";
	f << "uint64_t get_reciprocal(uint32_t a)\n";
	f << "{\n";
	f << "#if SMALL_DIVISORS_ENABLED\n";
	f << "	unsigned long shift;\n";
	f << "	_BitScanReverse(&shift, a);\n";
	f << "	shift = 31 - shift;\n";
	f << "	const uint32_t a1 = a << shift;\n\n";
	f << "	const uint32_t index1 = (a1 << 1) >> " << (INDEX2_BITS + 1) << ";\n";
	f << "	const int32_t index2 = (a1 & " << ((1UL << INDEX2_BITS) - 1) << ") - " << (1 << (INDEX2_BITS - 1)) << ";\n";
	f << "#else\n";
	f << "	const uint32_t index1 = (a << 1) >> " << (INDEX2_BITS + 1) << ";\n";
	f << "	const int32_t index2 = (a & " << ((1UL << INDEX2_BITS) - 1) << ") - " << (1 << (INDEX2_BITS - 1)) << ";\n";
	f << "#endif\n\n";
	f << "	const int64_t r0 = int64_t(RCP[index1].C0 | (1LL << " << C0_BITS << ")) << 22;\n";
	f << "	const int64_t r1 = int64_t(RCP[index1].C1) << 12;\n";
	f << "	const int64_t r2 = (int64_t(RCP[index1].C2) << 19) + (RCP[index1].C2 << 3);\n";
	f << "	const int32_t r3 = RCP[index1].C3;\n\n";
	f << "	int64_t rr;\n";
	f << "	rr = static_cast<int32_t>(-r3 * (index2 >> 4));\n";
	f << "	rr = ((rr + r2) * index2) >> 31;\n";
	f << "	rr = ((rr - r1) * index2) >> 10;\n";
	f << "	uint64_t r = static_cast<uint64_t>((rr + r0) >> 27);\n\n";
	f << "#if SMALL_DIVISORS_ENABLED\n";
	f << "	r <<= shift;\n";
	f << "	if (a <= " << smallest_divisor[0] << ")\n";
	f << "	{\n";
	f << "		uint64_t lo, hi, hi2;\n";
	f << "		lo = _umul128(r, a, &hi);\n";
	f << "		_subborrow_u64(_subborrow_u64(0, a >> 1, lo, &lo), 2, hi, &hi);\n";
	f << "		_umul128(r, lo, &hi2);\n";
	f << "		r = hi2 + (hi ? r : 0);\n";
	f << "	}\n";
	f << "#endif\n\n";
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
