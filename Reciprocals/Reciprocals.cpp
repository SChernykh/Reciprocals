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

constexpr uint32_t INDEX_BITS = 6;
constexpr int MAX_NEWTON_RAPHSON_ITER = 1;
constexpr int POLYNOMIAL_ORDER = 3;

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

template<int order>
void RemezApproximationStep(double a, double b, double (&coeff)[order + 1])
{
	constexpr int N = order + 2;

	double x[N], y[N];
	double m[N][N + 1];
	double x_base = (a + b) * 0.5;
	double y_base = 18446744073709551616.0 / x_base;

	for (int i = 0; i < N; ++i)
	{
		x[i] = (a - b) * 0.5 * cos(static_cast<double>(i) / (N - 1) * 3.1415926535897932384626433832795);
		y[i] = 18446744073709551616.0 / (x[i] + x_base) - y_base;

		m[i][0] = 1.0;
		for (int j = 1; j < N - 1; ++j)
		{
			m[i][j] = m[i][j - 1] * x[i];
		}
		m[i][N - 1] = ((i & 1) == 0) ? 1.0 : -1.0;
		m[i][N] = y[i];
	}

	for (int i = 0; i < N; ++i)
	{
		double t = m[i][i];
		if (t != 1.0)
		{
			m[i][i] = 1.0;
			for (int j = i + 1; j < N + 1; ++j)
			{
				m[i][j] /= t;
			}
		}
		for (int k = 0; k < N; ++k)
		{
			if (k != i)
			{
				t = m[k][i];
				for (int j = i; j < N + 1; ++j)
				{
					m[k][j] -= m[i][j] * t;
				}
			}
		}
	}

	coeff[0] = y_base + m[0][N];
	for (int i = 1; i <= order; ++i)
	{
		coeff[i] = m[i][N];
	}

	double err[N];
	double min_abs_err = 1e10;
	double max_abs_err = -1e10;
	for (int i = 0; i < N; ++i)
	{
		double value = 0.0;
		for (int j = N - 2; j >= 0; --j)
		{
			value = value * x[i] + coeff[j];
		}
		err[i] = value - y_base - y[i];
		if (fabs(err[i]) < min_abs_err) min_abs_err = fabs(err[i]);
		if (fabs(err[i]) > max_abs_err) max_abs_err = fabs(err[i]);
	}

	if (max_abs_err - min_abs_err > 1e-6)
	{
		printf("Sample points are suboptimal, approximation errors: ");
		for (int i = 0; i < N; ++i)
		{
			printf("%.06f ", err[i]);
		}
		printf("\n");
	}
}

struct
{
	uint32_t C0;
	float C[POLYNOMIAL_ORDER];

	int64_t value(float x) const
	{
		float delta = C[POLYNOMIAL_ORDER - 1];
		for (int i = POLYNOMIAL_ORDER - 1; i > 0; --i)
		{
			delta = delta * x + C[i - 1];
		}
		delta *= x;
		return C0 + _mm_cvtss_si32(_mm_set_ss(delta));
	}
} reciprocals[1UL << INDEX_BITS];

void approximation_init()
{
	for (uint64_t i = 0; i < (1UL << INDEX_BITS); ++i)
	{
		double RemezCoeffs[POLYNOMIAL_ORDER + 1];
		RemezApproximationStep<POLYNOMIAL_ORDER>(0x80000000ULL + (i << (31 - INDEX_BITS)), 0x80000000ULL + ((i + 1) << (31 - INDEX_BITS)), RemezCoeffs);
		reciprocals[i].C0 = static_cast<uint32_t>(_mm_cvttsd_si64(_mm_set_sd(RemezCoeffs[0] * 0.5)));
		for (int j = 0; j < POLYNOMIAL_ORDER; ++j)
		{
			reciprocals[i].C[j] = static_cast<float>(RemezCoeffs[j + 1] * 0.5);
		}
	}

	std::cout << "Lookup table size is " << (1UL << INDEX_BITS) << " entries" << std::endl;
	std::cout << "Calculating piecewise cubic approximation of 1/x with 64-bit precision (ULP = 2^-64)" << std::endl;
	for (uint64_t i = 0; i < (1UL << INDEX_BITS); ++i)
	{
		const uint64_t divisor1 = 0x80000000ULL + (i << (31 - INDEX_BITS));
		const uint64_t divisor2 = 0x80000000ULL + ((i + 1) << (31 - INDEX_BITS)) - 1;

		double min_error = 1e10;
		double max_error = -1e10;
		int64_t min_error_divisor = 0;
		int64_t max_error_divisor = 0;
		for (uint64_t a = divisor1; a <= divisor2; ++a)
		{
			uint64_t rem;
			uint64_t exact_r = udiv128(1, 0, a, &rem);
			double exact_r_fract = static_cast<double>(rem) / a;

			const float index1 = static_cast<float>(static_cast<int32_t>(a & ((1UL << (31 - INDEX_BITS)) - 1)) - (1 << (30 - INDEX_BITS)));
			int64_t r = ((reciprocals[i].value(index1) + 1) >> 1) << 2;

			const double err = static_cast<int64_t>(r - exact_r) - exact_r_fract;
			if (err < min_error)
			{
				min_error = err;
				min_error_divisor = a - (divisor1 + divisor2) / 2;
			}
			if (err > max_error)
			{
				max_error = err;
				max_error_divisor = a - (divisor1 + divisor2) / 2;
			}
		}
		int32_t fix = static_cast<int32_t>(round((min_error + max_error) * 0.5) + 2) >> 2;
		reciprocals[i].C0 -= fix << 1;
		min_error -= fix * 4;
		max_error -= fix * 4;
		std::cout << "err[" << i << "] = (" << min_error << ", " << max_error << ") ULP at (" << min_error_divisor << ", " << max_error_divisor << ')' << std::endl;
	}
	std::cout << std::endl;
}

void approximation_check()
{
	std::cout << "Calculating maximum error after up to " << MAX_NEWTON_RAPHSON_ITER << " Newton-Raphson iteration(s)" << std::endl << std::endl;

	for (int i = 0; i < MAX_NEWTON_RAPHSON_ITER; ++i)
	{
		min_error[i] = 1e10;
		max_error[i] = -1e10;
	}

	uint64_t prev_index = 0;
	for (uint64_t a = 0xFFFFFFFFUL; a >= 3; a -= 2)
	{
		if ((a >> 20) != prev_index)
		{
			std::cout << std::setw(4) << (a >> 20) << '\r';
			prev_index = (a >> 20);
		}

		int64_t r;

		uint64_t exact_r = udiv128(1, 0, a, (uint64_t*) &r);
		double exact_r_fraction = static_cast<double>(r) / a;

		unsigned long shift;
		_BitScanReverse64(&shift, a);
		shift = 31 - shift;
		const uint64_t a1 = a << shift;

		const uint32_t index0 = (static_cast<uint32_t>(a1 << 1) >> (32 - INDEX_BITS));
		const float index1 = static_cast<float>(static_cast<int32_t>(a1 & ((1UL << (31 - INDEX_BITS)) - 1)) - (1 << (30 - INDEX_BITS)));
		r = (reciprocals[index0].value(index1) + 1) >> 1;

		for (int i = 0; i < MAX_NEWTON_RAPHSON_ITER; ++i)
		{
			int64_t h = int64_t((1LL << (62 - shift)) - (r * a));
			const int h_shift = (shift < 2) ? (2 - shift) : 0;
			h >>= h_shift;

			const int64_t r1 = (r * h) >> (64 - ((shift + 2) * 2 + h_shift));
			r = (r << (shift + 2)) + r1 + 1;

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
			if ((r != exact_r) && (r != exact_r + 1) && (smallest_divisor[i] == 0))
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

extern uint64_t get_reciprocal(uint64_t a);
extern void fast_div(uint64_t a, uint32_t b, uint64_t &q, uint64_t &r);
extern int64_t min_h, max_h;

void check_generated_code()
{
	std::cout << "Checking generated code" << std::endl << std::endl;

	double min_err = 1e10;
	double max_err = -1e10;
	uint64_t min_err_divisor = 0;
	uint64_t max_err_divisor = 0;
	uint64_t prev_index = 0;
	for (uint32_t a = 0xFFFFFFFFUL; a >= 3; a -= 2)
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
		const uint32_t divisor1 = static_cast<uint32_t>(n.second) | 5;
		const uint32_t divisor2 = static_cast<uint32_t>(n.second >> 32) | 5;

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
		const uint32_t divisor1 = static_cast<uint32_t>(n.second) | 5;
		const uint32_t divisor2 = static_cast<uint32_t>(n.second >> 32) | 5;
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
		const uint32_t divisor1 = static_cast<uint32_t>(n.second) | 5;
		const uint32_t divisor2 = static_cast<uint32_t>(n.second >> 32) | 5;
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

//#define PLATFORM_32_BIT 1

__forceinline void fast_div_v2(const uint64_t a, const uint32_t b, uint32_t& q, uint32_t& r)
{
	const float rcp = 1.0f / static_cast<float>(b);

	float rcp2;
	*((uint32_t*)&rcp2) = *((uint32_t*)&rcp) + (32U << 23);

	float q1f = static_cast<float>(((uint32_t*)&a)[1]) * rcp2;
	uint32_t q1 = static_cast<uint32_t>(round(q1f));
	int64_t a1 = static_cast<int64_t>(a - uint64_t(q1) * b);

	((int32_t*)&a1)[1] -= (q1f >= 4294967296.0f) ? b : 0;

	int32_t q2 = static_cast<int32_t>(round(static_cast<float>(a1) * rcp));
	int32_t sign = q2 >> 31;
	int64_t a2 = static_cast<uint64_t>(static_cast<uint32_t>((q2 ^ sign) - sign)) * b;
	a1 = (q2 < 0) ? (a1 + a2) : (a1 - a2);

	sign = ((int32_t*)&a1)[1] >> 31;
	q = q1 + q2 + sign;
	r = a1 + (static_cast<uint32_t>(sign) & b);
}

__forceinline int64_t fast_div_heavy(int64_t a, const int32_t b)
{
	const int32_t sign_a = ((int32_t*)&a)[1] >> 31;
	const int32_t sign_b = b >> 31;
	const int32_t sign2 = ((sign_a ^ sign_b) << 1) + 1;

	const float rcp = 1.0f / static_cast<float>(b);

	float rcp2;
	*((uint32_t*)&rcp2) = *((uint32_t*)&rcp) + (32U << 23);

	int64_t q1 = static_cast<int64_t>(round(static_cast<float>(((int32_t*)&a)[1]) * rcp2));
	a -= q1 * b;

	*((uint32_t*)&rcp2) = *((uint32_t*)&rcp) + (12U << 23);
	int64_t q2 = static_cast<int64_t>(round(static_cast<float>(static_cast<int32_t>(a >> 12)) * rcp2));
	a -= q2 * b;

	int32_t q3 = static_cast<int32_t>(round(static_cast<float>(static_cast<int32_t>(a)) * rcp));
	int32_t r = static_cast<int32_t>(a) - q3 * b;
	r = (r ^ sign_a) - sign_a;
	q3 -= sign2 & (r >> 31);

	return q1 + q2 + q3;
}

#include <Windows.h>
#include <thread>

__declspec(noinline) void test_div_v2(uint32_t d1, uint32_t d2, uint64_t seed, uint64_t id)
{
	SetThreadAffinityMask(GetCurrentThread(), 1ULL << id);
	SetThreadPriority(GetCurrentThread(), THREAD_MODE_BACKGROUND_BEGIN);
	SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_IDLE);

	printf("v2: id = %llu, d1 = %u, d2 = %u\n", id, d1, d2);

	for (uint32_t d = d1; d < d2; d += 2)
	{
		if ((d1 == 0) && ((d % 65536) == 0) && (d > 0))
		{
			std::cout << "v2: " << d << '\r';
		}
		for (uint32_t i = 0; i < 64; ++i)
		{
			seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
			uint64_t a = seed;
			uint32_t b = d | 0x80000001UL;
			uint32_t q, r;
			fast_div_v2(a, b, q, r);
			if ((q != static_cast<uint32_t>(a / b)) || (r != a % b))
			{
				printf("%llu / %u != %u\n%u\n%u\n%d\n\n", a, b, q, static_cast<uint32_t>(a / b), q, static_cast<uint32_t>(a / b) - q);
				__debugbreak();
			}
		}
	}
}

__declspec(noinline) void test_div_heavy(uint32_t d1, uint32_t d2, uint64_t seed, uint64_t id)
{
	SetThreadAffinityMask(GetCurrentThread(), 1ULL << id);
	SetThreadPriority(GetCurrentThread(), THREAD_MODE_BACKGROUND_BEGIN);
	SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_IDLE);

	d1 |= 1;
	if (d1 < 3) d1 = 3;

	printf("heavy: id = %llu, d1 = %u, d2 = %u\n", id, d1, d2);

	for (uint32_t d = d1; d < d2; d += 2)
	{
		if ((d1 == 3) && ((d % 65536) == 3) && (d > 3))
		{
			std::cout << "heavy: " << d << '\r';
		}
		for (uint32_t i = 0; i < 64; ++i)
		{
			seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
			int64_t a = static_cast<int64_t>(seed);
			int32_t b = d;
			int64_t q = fast_div_heavy(a, b);
			if (q != a / b)
			{
				printf("%lld / %d != %lld\n%lld\n%lld\n%lld\n\n", a, b, q, a / b, q, (a / b) - q);
				__debugbreak();
			}
			q = fast_div_heavy(a, -b);
			if (q != a / -b)
			{
				printf("%lld / %d != %lld\n%lld\n%lld\n%lld\n\n", a, -b, q, a / -b, q, (a / -b) - q);
				__debugbreak();
			}
		}
	}
}

int main()
{
	_control87(RC_NEAR, MCW_RC);
	{
		SetPriorityClass(GetCurrentProcess(), IDLE_PRIORITY_CLASS);

		uint32_t q, r;
		fast_div_v2(9468199478243335150, 2147483649, q, r);
		fast_div_heavy(-2614408371619423461, -134217729);

		constexpr int N = 16;
		std::vector<std::thread> threads;
		threads.reserve(N);

		//for (uint64_t i = 0; i < N; ++i)
		//{
		//	const uint32_t d1 = (i << 31) / N;
		//	const uint32_t d2 = ((i + 1) << 31) / N;
		//	const uint64_t seed = i + 123;
		//	threads.emplace_back(test_div_v2, d1, d2, seed, i);
		//}
		//for (uint64_t i = 0; i < N; ++i)
		//{
		//	threads[i].join();
		//}
		//threads.clear();

		for (uint64_t i = 0; i < N; ++i)
		{
			const uint32_t d1 = (i << 31) / N;
			const uint32_t d2 = ((i + 1) << 31) / N;
			const uint64_t seed = i + 123;
			threads.emplace_back(test_div_heavy, d1, d2, seed, i);
		}
		for (uint64_t i = 0; i < N; ++i)
		{
			threads[i].join();
		}
		threads.clear();
	}

	std::cout.precision(15);

	check_generated_code();
	check_and_benchmark_generated_code();
	approximation_init();
	approximation_check();

	std::ofstream f("reciprocal64.cpp");

	f << "#include <stdint.h>\n";
	f << "#include <intrin.h>\n\n";
	f << "static const uint32_t R[" << (sizeof(reciprocals) / sizeof(uint32_t)) << "] = {\n";
	f << std::hex;
	for (int i = 0; i < sizeof(reciprocals) / sizeof(uint32_t); i += 8)
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
	f << "	float C[3];\n\n";
	f << "	int64_t value(float x) const\n";
	f << "	{\n";
	f << "		float delta = C[" << (POLYNOMIAL_ORDER - 1) << "];\n";
	for (int i = POLYNOMIAL_ORDER - 2; i >= 0; --i)
	{
		f << "		delta = delta * x + C[" << i << "];\n";
	}
	f << "		delta *= x;\n";
	f << "		return C0 + _mm_cvtss_si32(_mm_set_ss(delta));\n";
	f << "	}\n";
	f << "};\n\n";
	f << "#define RCP ((SReciprocal*)R)\n\n";
	f << "uint64_t get_reciprocal(uint64_t a)\n";
	f << "{\n";
	f << "	unsigned long shift;\n";
	f << "	_BitScanReverse(&shift, static_cast<uint32_t>(a));\n";
	f << "	shift = 31 - shift;\n";
	f << "	const uint64_t a1 = a << shift;\n";
	f << "\n";
	f << "	const uint32_t index0 = static_cast<uint32_t>(a1 << 1) >> " << (32 - INDEX_BITS) << ";\n";
	f << "	const float index1 = " << "static_cast<float>(static_cast<int32_t>(a1 & " << ((1UL << (31 - INDEX_BITS)) - 1) << ") - " << (1 << (30 - INDEX_BITS)) << ");\n";
	f << "	const int64_t r = (RCP[index0].value(index1) + 1) >> 1;\n\n";
	f << "	int64_t h = int64_t((1LL << (62 - shift)) - (r * a));\n";
	f << "	const int h_shift = (shift < 2) ? (2 - shift) : 0;\n";
	f << "	h >>= h_shift;\n\n";
	f << "	const int32_t r1 = (r * h) >> (64 - ((shift + 2) * 2 + h_shift));\n";
	f << "	return (r << (shift + 2)) + r1 + 1;\n";
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
