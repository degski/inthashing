
#include <cstdint>
#include <ciso646>
#include <cmath>
#include <intrin.h>

#include <iostream>
#include <utility> // std::move...
#include <functional> // std::greater
#include <vector>
#include <algorithm>
#include <numeric>

#include <boost/random/xoroshiro.hpp>
#include <boost/random/seed_seq_fe.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>

#include <integer_utils.hpp>

#include <autotimer.hpp>

/*
https://crypto.stackexchange.com/questions/34269/calculation-of-the-avalanche-effect-coefficient#35172

I've done it this way for SHA 512:-

1. Generate a random number 512 bits long, called O

2. Randomly flip one bit to generate number F

3. Compute X = SHA(O) xor SHA(F)

4. Calculate no. of set bits in X

5. Coefficient Ksac = X / 512 (note)

Rinse and repeat a million times to get a mean for Ksac

Standard deviation should be 0.5 * sqrt ( N = 512 ) = 11.3137...

*/

namespace inthashing {

	template<typename T>
	T hash ( T x_, const T m_ ) {

		x_ = ( ( x_ >> ( ( sizeof ( T ) / 2 ) * 8 ) ) ^ x_ ) * m_;
		x_ = ( ( x_ >> ( ( sizeof ( T ) / 2 ) * 8 ) ) ^ x_ ) * m_;
		//x_ = ( ( x_ >> ( ( sizeof ( T ) / 2 ) * 8 ) ) ^ x_ ) * m_;

		return ( x_ >> ( ( sizeof ( T ) / 2 ) * 8 ) ) ^ x_;
	}

	template<typename T>
	T mul_m ( T x_, const T m_ ) {

		return x_ * m_;
	}
};


template < typename N >
N set_mask_0 ( const std::uint32_t bit_ ) noexcept {

	return N ( 1 ) << bit_;
}


// boost::random::seed_seq_fe256 g_seq { 36534427, 234237, 144433, 2144421 };
boost::random::xoroshiro128plus g_rng ( iu::seed<typename boost::random::xoroshiro128plus::result_type> ( ) );


template<typename T>
T getRandom ( ) noexcept {

	return boost::random::uniform_int_distribution<T> ( ) ( g_rng );
}


template<typename T>
T flipBit ( const T i_, const std::uint32_t n_ ) noexcept {

	return i_ ^ set_mask_0<T> ( n_ );
}

template<typename T>
T flipRandomBit ( const T i_ ) noexcept {

	static boost::random::uniform_int_distribution<std::uint32_t> uid ( 0, sizeof ( T ) * 8 - 1 );

	return flipBit ( i_, uid ( g_rng ) );
}


// strict avalanche criterion
// https://en.wikipedia.org/wiki/Avalanche_effect

template<typename T>
double getKsac ( const T x_, const T y_, const T m_ ) noexcept {

	return ( double ) iu::popCount ( iu::bit_xor ( inthashing::hash ( x_, m_ ), inthashing::hash ( y_, m_ ) ) ) / double ( sizeof ( T ) * 8 );
}

template<typename T>
double getRandomKsac ( const T m_ ) noexcept {

	const T x = getRandom<T> ( );

	return getKsac ( x, flipRandomBit ( x ), m_ );
}

template<typename T>
double getLowEntropyKsac ( const T m_ ) noexcept {

	static T x = getRandom<T> ( );

	++x;

	return getKsac ( x, flipRandomBit ( x ), m_ );
}


template<typename T>
double getFsacError ( const T x_, const T m_ ) noexcept {

	double error = getKsac ( x_, flipBit ( x_, 0 ), m_ ) - 0.5, mean_squared_error = error * error;

	for ( std::uint32_t i = 1; i < sizeof ( T ) * 8; ++i ) {

		error = getKsac ( x_, flipBit ( x_, i ), m_ ) - 0.5;

		mean_squared_error += ( error * error - mean_squared_error ) / i;
	}

	return mean_squared_error;
}

template<typename T>
double getTotalFsacError ( const T m_ ) noexcept {

	double mean_error = getFsacError ( ( T ) 0, m_ );

	for ( std::uint32_t i = 1; i < std::numeric_limits<T>::max ( ); ++i ) {

		mean_error += ( getFsacError ( ( T ) i, m_ ) - mean_error ) / i;
	}

	return mean_error;
}

template<typename T>
T getBestM ( ) noexcept {

	T best_m = 3;
	double best_error = getTotalFsacError ( best_m );

	//std::cout << "i " << ( std::uint32_t ) 0 << " " << " error " << best_error << '\n';

	for ( std::uint32_t i = 5; i <std::numeric_limits<T>::max ( ); i += 2 ) {

		const double error = getTotalFsacError ( ( T ) i );

		if ( error < best_error ) {

			best_m = ( T ) i;
			best_error = error;

			///std::cout << "i " << ( std::uint32_t ) i << " " << " best error " << best_error << '\n';

			continue;
		}

		std::cout << "i " << ( std::uint32_t ) i << " " << " error " << error << " best m " << best_m << " best error " << best_error << '\n';
	}

	return best_m;
}

int main3345678 ( ) {

	std::cout << getBestM<std::uint16_t> ( ) << '\n';

	return 0;
}

template<typename T>
double getRandomFsacError ( const T m_ ) noexcept {

	return getFsacError ( getRandom<T> ( ), m_ );
}

template<typename T>
double getLowEntropyFsacError ( const T m_ ) noexcept {

	static T x = getRandom<T> ( );

	return getFsacError ( x++, m_ );
}

template<typename T>
double getCombinedFsacError ( const T m_, const size_t i_ ) noexcept {

	double mean_error = 0.0, error;

	std::size_t i = 1, last = i_ + 1;

	for ( ; i < last; ++i ) {

		error = getRandomFsacError ( m_ );

		mean_error += ( error - mean_error ) / i;
	}

	last = i_ * 2 + 1;

	for ( ; i < last; ++i ) {

		error = getLowEntropyFsacError ( m_ );

		mean_error += ( error - mean_error ) / i;
	}

	return mean_error;
}

template<typename T>
double getCombinedKsacMSR ( const T m_, const size_t i_ ) noexcept {

	double mean_squared_error = 0.0, error;

	size_t i = 1, last = i_ + 1;

	for ( ; i < last; ++i ) {

		error = getRandomKsac<T> ( m_ ) - 0.5;

		mean_squared_error += ( error * error - mean_squared_error ) / i;
	}

	last = i_ * 2 + 1;

	for ( ; i < last; ++i ) {

		error = getLowEntropyKsac<T> ( m_ ) - 0.5;

		mean_squared_error += ( error * error - mean_squared_error ) / i;
	}

	return mean_squared_error;
}


template<typename T>
bool isInverse ( const T m1_, const T m2_ ) noexcept {

	for ( uint32_t i = 0; i < 1'000'000'000; ++i ) {

		const T x = iu::make_odd ( getRandom<std::uint64_t> ( ) );

		if ( iu::hash ( iu::hash ( x, m1_ ), m2_ ) != x ) {

			return false;
		}
	}

	return true;
}


int32_t main455435 ( ) {

	const char format_string [ ] { " %10u - 0x%016llX 0x%016llX - %.10f\n" };

	uint64_t best_m = iu::make_odd ( getRandom<std::uint64_t> ( ) );

	double avg_ksacd = getCombinedKsacMSR ( best_m, 500 );
	uint64_t a_cnt = 1;

	double best_ksacd = getCombinedKsacMSR ( best_m, 50'000 );

	printf ( format_string, 0u, best_m, iu::mod_mul_inv ( best_m ), best_ksacd );

	for ( uint32_t i = 1; i < UINT32_MAX; ++i ) {

		const uint64_t m = iu::make_odd ( getRandom<std::uint64_t> ( ) );

		const double a_ksacd = getCombinedKsacMSR ( m, 500 );

		if ( a_ksacd < avg_ksacd ) {

			avg_ksacd += ( a_ksacd - avg_ksacd ) / ++a_cnt;

			const double b_ksacd = getCombinedKsacMSR ( m, 50'000 );

			if ( b_ksacd < best_ksacd ) {

				best_m = m;
				best_ksacd = b_ksacd;

				printf ( format_string, i, best_m, iu::mod_mul_inv ( m ), best_ksacd );
			}
		}
	}

	return 0;
}

template<typename T>
struct candidate {

	T value;
	std::size_t eval_unit, evaluations, ctr = 1;
	double score;

	candidate ( ) : value ( iu::make_odd ( getRandom<std::uint64_t> ( ) ) ), eval_unit ( 6 * 1'024 ), evaluations ( 0 ), score ( getCombinedKsacMSR ( value, eval_unit ) ) { }
	candidate ( const T v_, const std::size_t e_ ) : value ( v_ ), evaluations ( e_ ), score ( getCombinedKsacMSR ( value, evaluations ) ) { }

	void evaluate ( ) noexcept {

		evaluations += eval_unit;
		score += ( getCombinedKsacMSR ( value, eval_unit ) - score ) / ++ctr;
	}

	bool operator < ( candidate & rhs_ ) noexcept {

		return score < rhs_.score;
	}

	static constexpr std::size_t max ( ) { return 1024 * 1024 * 4; }
};

int32_t main ( ) {

	/*
	const std::uint64_t tested = 0x39FF916554A55555, initial = 0x0CF3FD1B9997F637;

	std::uint64_t cnt = ~0ULL - 1ULL, value = initial;

	while ( cnt-- ) {

		value *= tested;

		if ( value == initial ) {

			std::cout << "cycle after " << ( ~0ULL - cnt ) << '\n';

			exit ( 0 );
		}
	}

	std::cout << "no cycles\n";

	exit ( 0 );
	*/

	using int_type = std::uint64_t;

	const char format_string [ ] { " %10u - 0x%016llX 0x%016llX - %.10f - %20llu\n" };

	constexpr std::size_t pop_size = 16 * 1024;

	std::vector<candidate<int_type>> population ( pop_size );
	population.reserve ( pop_size );

	std::sort ( std::begin ( population ), std::end ( population ) );

	for ( std::uint32_t i = 0; ; ++i ) {

		for ( auto & p : population ) {

			p.evaluate ( );
		}

		std::sort ( std::begin ( population ), std::end ( population ) );

		// printf ( format_string, i, population [ 0 ].value, 0ULL, population [ 0 ].score, population [ 0 ].evaluations );

		printf ( format_string, i, ( std::uint64_t ) population [ 0 ].value, ( std::uint64_t ) iu::mod_mul_inv ( population [ 0 ].value ), population [ 0 ].score, population [ 0 ].evaluations );
		printf ( format_string, i, ( std::uint64_t ) population [ 1 ].value, ( std::uint64_t ) iu::mod_mul_inv ( population [ 1 ].value ), population [ 1 ].score, population [ 1 ].evaluations );
		printf ( format_string, i, ( std::uint64_t ) population [ 2 ].value, ( std::uint64_t ) iu::mod_mul_inv ( population [ 2 ].value ), population [ 2 ].score, population [ 2 ].evaluations );

		std::cout << std::endl;

		std::size_t m = pop_size * 0.6f;

		while ( m < pop_size ) {

			population [ m++ ] = candidate<int_type> { };
		}
	}

	return 0;
}




// Release:
// -Xclang -fcxx-exceptions -Xclang -std=c++14 -Qunused-arguments -Wno-unused-variable -Xclang -O3 -Xclang -ffast-math -mmmx  -msse  -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -mavx2 -Xclang -Wno-deprecated-declarations -Xclang -Wno-unknown-pragmas -Xclang -Wno-ignored-pragmas -Xclang -Wno-unused-private-field



template<typename T>
constexpr T gcd1 ( const T a_, const T b_ ) {

	return a_ < b_ ? gcd1 ( b_, a_ ) : a_ % b_ == 0 ? b_ : gcd1 ( b_, a_ % b_ );
}


template<typename T>
T gcd2 ( T a_, T b_ ) {

	repeat:

	if ( ! a_ ) {

		return b_;
	}

	b_ %= a_;

	if ( ! b_ ) {

		return a_;
	}

	a_ %= b_;

	goto repeat;
}


int interp_cgoto ( unsigned char* code, int initval ) {

	// The indices of labels in the dispatch_table are the relevant opcodes...

	static void* dispatch_table [ ] = { &&do_halt, &&do_inc, &&do_dec, &&do_mul2, &&do_div2, &&do_add7, &&do_neg };

#define DISPATCH() goto *dispatch_table [ code [ pc++ ] ]

	int pc = 0;
	int val = initval;

	DISPATCH ( );
	while ( 1 ) {
	do_halt:
		return val;
	do_inc:
		val++;
		DISPATCH ( );
	do_dec:
		val--;
		DISPATCH ( );
	do_mul2:
		val *= 2;
		DISPATCH ( );
	do_div2:
		val /= 2;
		DISPATCH ( );
	do_add7:
		val += 7;
		DISPATCH ( );
	do_neg:
		val = -val;
		DISPATCH ( );
	}
}




int main9878970 ( ) {

	return 0;
}


/*

it 0 18430510334115963153 10335394720245952497 0.004008422851562498667732370449812151491641998291015625000000
it 1 7477738296582267083 7148602372516906211 0.003883107800292810544984245524346988531760871410369873046875
it 4 5420600479766853331 17806009454484328795 0.003866896215820473346747032650227993144653737545013427734375
it 19 12077015678325397305 3272976820258023177 0.003852416577147856849988594163392008340451866388320922851563
it 33 5895393256363349695 3348953932214376767 0.003846883691406476936391900522949072183109819889068603515625
it 632 12673909164712089953 8876396653587663521 0.003822637756348148431462252716528382734395563602447509765625
it 917 5553275399470672577 2427780294986310977 0.003792151782226740961562194698331040854100137948989868164063
it 33394 14049795677467858241 8526623269437389505 0.003787331982421852753784019540717054042033851146697998046875
it 45138 18410416285450594625 1605937283254986433 0.003785130444336448109210468970786678255535662174224853515625
it 46527 6699509707676226881 2365871517788338881 0.003780415039062876274983393543038800999056547880172729492188
it 66566 15902845056476064385 2464940388702033281 0.003755318481446032023718384351695931400172412395477294921875
it 5628951 6783152649535924737 12664843139132773889 0.003752792846679828264733203013747697696089744567871093750000
it 7548686 16003472649440565889 18247563199220589953 0.003744663146972564398556881926083406142424792051315307617188
it 75930845 15270417933577923073 10310843644990019073 0.003742728759764692419292897440641354478430002927780151367188
it 82109243 10064833879070834689 16371722055553024001 0.003729221057129303035920786513202074274886399507522583007813
it 150768790 14831393149721859073 1360584342912019457 0.003719900988770050154208490766905015334486961364746093750000
it 266552271 2064167726340199425 2993319394596793345 0.003703780920410166752065883599698281614109873771667480468750
0xB25880B9598CFE8D 0x606EB83A7025F445 - 0.0034293757
0x1AEC805299990163 0xCDFB859A3DD0884B - 0.0033731416
0x329FF22C54C063B9 0x7E75D92ADF54B289 - 0.0034275592
0xDBBF59B09980D163 0x7177B03E4D47384B - 0.0033654531
0x0CF3FD1B9997F637 0xAFC1530680179F87 - 0.0033298184
*/

struct Snitch {   // Note: All methods have side effects
	Snitch ( ) { std::cout << "c'tor" << std::endl; }
	~Snitch ( ) { std::cout << "d'tor" << std::endl; }

	Snitch ( const Snitch& ) { std::cout << "copy c'tor" << std::endl; }
	Snitch ( Snitch&& ) { std::cout << "move c'tor" << std::endl; }

	Snitch& operator=( const Snitch& ) {
		std::cout << "copy assignment" << std::endl;
		return *this;
	}

	Snitch& operator=( Snitch&& ) {
		std::cout << "move assignment" << std::endl;
		return *this;
	}
};



/*

__uint128_t FASTMUL128 ( const __uint128_t TA, const __uint128_t TB ) {
union {
__uint128_t WHOLE;
struct {
unsigned long long int LWORDS [ 2 ];
} SPLIT;
} KEY;
register unsigned long long int __RAX, __RDX, __RSI, __RDI;
__uint128_t RESULT;

KEY.WHOLE = TA;
__RAX = KEY.SPLIT.LWORDS [ 0 ];
__RDX = KEY.SPLIT.LWORDS [ 1 ];
KEY.WHOLE = TB;
__RSI = KEY.SPLIT.LWORDS [ 0 ];
__RDI = KEY.SPLIT.LWORDS [ 1 ];
__asm__ __volatile__ (
"movq           %0,                             %%rax                   \n\t"
"movq           %1,                             %%rdx                   \n\t"
"movq           %2,                             %%rsi                   \n\t"
"movq           %3,                             %%rdi                   \n\t"
"movq           %%rsi,                          %%rbx                   \n\t"
"movq           %%rdi,                          %%rcx                   \n\t"
"movq           %%rax,                          %%rsi                   \n\t"
"movq           %%rdx,                          %%rdi                   \n\t"
"xorq           %%rax,                          %%rax                   \n\t"
"xorq           %%rdx,                          %%rdx                   \n\t"
"movq           %%rdi,                          %%rax                   \n\t"
"mulq           %%rbx                                                   \n\t"
"xchgq          %%rbx,                          %%rax                   \n\t"
"mulq           %%rsi                                                   \n\t"
"xchgq          %%rax,                          %%rsi                   \n\t"
"addq           %%rdx,                          %%rbx                   \n\t"
"mulq           %%rcx                                                   \n\t"
"addq           %%rax,                          %%rbx                   \n\t"
"movq           %%rsi,                          %%rax                   \n\t"
"movq           %%rbx,                          %%rdx                   \n\t"
"movq           %%rax,                          %0                      \n\t"
"movq           %%rdx,                          %1                      \n\t"
"movq           %%rsi,                          %2                      \n\t"
"movq           %%rdi,                          %3                      \n\t"
: "=m"( __RAX ), "=m"( __RDX ), "=m"( __RSI ), "=m"( __RDI )
: "m"( __RAX ), "m"( __RDX ), "m"( __RSI ), "m"( __RDI )
: "rax", "rbx", "ecx", "rdx", "rsi", "rdi"
);
KEY.SPLIT.LWORDS [ 0 ] = __RAX;
KEY.SPLIT.LWORDS [ 1 ] = __RDX;
RESULT = KEY.WHOLE;
return RESULT;
}


{
uint64_t hi, lo;
// hi,lo = 64bit x 64bit multiply of c[0] and b[0]

__asm__ ( "mulq %3\n\t"
: "=d" ( hi ),
"=a" ( lo )
: "%a" ( c [ 0 ] ),
"rm" ( b [ 0 ] )
: "cc" );

a [ 0 ] += hi;
a [ 1 ] += lo;
}

*/
