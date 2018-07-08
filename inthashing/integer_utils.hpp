
#pragma once

#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <immintrin.h>
#include <cstdint>

#include <type_traits>

#include <boost/random/uniform_int_distribution.hpp>


#ifdef NDEBUG
#pragma comment ( lib, "integer_utils-s.lib" )
#else
#pragma comment ( lib, "integer_utils-s-d.lib" )
#endif


//
//                       _oo0oo_
//                      o8888888o
//                      88" . "88
//                      (| -_- |)
//                      0\  =  /0
//                    ___/`---'\___
//                  .' \\|     |// '.
//                 / \\|||  :  |||// \
//                / _||||| -:- |||||- \
//               |   | \\\  -  /// |   |
//               | \_|  ''\---/''  |_/ |
//               \  .-\__  '-'  ___/-. /
//             ___'. .'  /--.--\  `. .'___
//          ."" '<  `.___\_<|>_/___.' >' "".
//         | | :  `- \`.;`\ _ /`;.`/ - ` : | |
//         \  \ `_.   \_ __\ /__ _/   .-` /  /
//     =====`-.____`.___ \_____/___.-`___.-'=====
//                       `=---='
//


namespace jimi {


	// Greatest Common Denominator...

	template < typename T, typename = std::enable_if_t < std::is_unsigned < T >::value && std::is_integral < T >::value, T > >
	inline constexpr T gcd (  T a_, T b_ ) {

		while ( true ) {

			if ( !a_ ) {

				return b_;
			}

			b_ %= a_;

			if ( !b_ ) {

				return a_;
			}

			a_ %= b_;
		}
	}


	// Least Common Multiple...

	template < typename T, typename = std::enable_if_t < std::is_signed < T >::value && std::is_integral < T >::value, T > >
	inline constexpr T lcm ( const T a_, const T b_ ) {

		const T t = gcd < T > ( a_, b_ );

		return t ? a_ / t * b_ : T ( 0 );
	}


	// In number theory, two integers a and b are said to be
	// relatively prime, mutually prime, or coprime (also
	// spelled co-prime)[1] if the only positive integer that
	// divides both of them is 1...

	template < typename T, typename = std::enable_if_t < std::is_signed < T >::value && std::is_integral < T >::value, T > >
	inline constexpr bool areCoPrime ( const T a_, const T b_ ) {

		// Assumes a_ and b_ > 0...

		return gcd < T > ( a_, b_ ) == T ( 1 );
	}


	// Integer Log2...

	template < typename T, typename = std::enable_if_t < std::is_unsigned < T >::value && std::is_integral < T >::value, T > >
	constexpr T iLog2 ( const T n_, const T p_ = T ( 0 ) ) {

		return n_ <= 1 ? p_ : iLog2 < T > ( n_ / 2, p_ + 1 );
	}


	template < typename T, typename = std::enable_if_t < std::is_unsigned < T >::value && std::is_integral < T >::value, T > >
	constexpr T nextPowerOfTwo ( const T n_ ) {

		return n_ > 2 ? T ( 1 ) << ( iLog2 < T > ( n_ - 1 ) + 1 ) : n_;
	}


	template < typename T, typename = std::enable_if_t < std::is_unsigned < T >::value && std::is_integral < T >::value, T > >
	constexpr bool isPowerOfTwo ( const T n_ ) {

		return n_ && ! ( n_ & ( n_ - 1 ) );
	}


	template < typename T, typename = std::enable_if_t < std::is_unsigned < T >::value && std::is_integral < T >::value, T > >
	constexpr T sumToN ( const T n_ ) {

		return ( n_ * ( n_ + 1 ) ) / T ( 2 );
	}


	template < typename T, typename = std::enable_if_t < std::is_unsigned < T >::value && std::is_integral < T >::value, T > >
	constexpr T sumMToN ( const T m_, const T n_ ) {

		return ( ( n_ * ( n_ + 1 ) ) - ( m_ * ( m_ + 1 ) ) ) / T ( 2 );
	}


	// Pointer Alignment...

	inline uint32_t pointerAlignment ( const void * p_ ) {

		return ( uint32_t ) ( ( uintptr_t ) p_ & ( uintptr_t ) -( ( intptr_t ) p_ ) );
	}


	// Gray Coding...

	template < typename T, typename = std::enable_if_t < std::is_unsigned < T >::value && std::is_integral < T >::value, T > >
	inline constexpr T decimalToGray ( const T i_ ) {

		return i_ ^ ( i_ >> 1 );
	}


	inline uint8_t grayToDecimal ( uint8_t g_ ) {

		g_ ^= g_ >> 4;
		g_ ^= g_ >> 2;
		g_ ^= g_ >> 1;

		return g_;
	}


	inline uint16_t grayToDecimal ( uint16_t g_ ) {

		g_ ^= g_ >> 8;
		g_ ^= g_ >> 4;
		g_ ^= g_ >> 2;
		g_ ^= g_ >> 1;

		return g_;
	}


	inline uint32_t grayToDecimal ( uint32_t g_ ) {

		g_ ^= g_ >> 16;
		g_ ^= g_ >> 8;
		g_ ^= g_ >> 4;
		g_ ^= g_ >> 2;
		g_ ^= g_ >> 1;

		return g_;
	}


	inline uint64_t grayToDecimal ( uint64_t g_ ) {

		g_ ^= g_ >> 32;
		g_ ^= g_ >> 16;
		g_ ^= g_ >>  8;
		g_ ^= g_ >>  4;
		g_ ^= g_ >>  2;
		g_ ^= g_ >>  1;

		return g_;
	}


	uint32_t modularMultiplicativeInverse ( const uint32_t a_ );
	uint64_t modularMultiplicativeInverse ( const uint64_t a_ );


	bool isPrime ( const uint32_t n_ ); // Odd only...
	bool isPrime ( const uint64_t n_ ); // Odd only...


	// Integer Hashing...

	inline uint32_t hash ( uint32_t x ) {

		x = ( ( x >> 16 ) ^ x ) * 0x45d9f3b;
		x = ( ( x >> 16 ) ^ x ) * 0x45d9f3b;
		x = ( ( x >> 16 ) ^ x );

		return x;
	}


	inline uint32_t unHash ( uint32_t x ) {

		x = ( ( x >> 16 ) ^ x ) * 0x119de1f3;
		x = ( ( x >> 16 ) ^ x ) * 0x119de1f3;
		x = ( ( x >> 16 ) ^ x );

		return x;
	}

	// 0x0CF3FD1B9997F637 0xAFC1530680179F87 - 0.0033298184

	inline uint64_t hash ( uint64_t x ) {

		x = ( ( x >> 32 ) ^ x ) * 0x0CF3FD1B9997F637;
		x = ( ( x >> 32 ) ^ x ) * 0x0CF3FD1B9997F637;
		x = ( ( x >> 32 ) ^ x );

		return x;
	}


	inline uint64_t unHash ( uint64_t x ) {

		x = ( ( x >> 32 ) ^ x ) * 0xAFC1530680179F87;
		x = ( ( x >> 32 ) ^ x ) * 0xAFC1530680179F87;
		x = ( ( x >> 32 ) ^ x );

		return x;
	}


	uint32_t popCount ( const uint8_t x_ );
	uint32_t popCount ( const uint16_t x_ );
	uint32_t popCount ( const uint32_t x_ );
	uint32_t popCount ( const uint64_t x_ );


	template < typename T, typename = std::enable_if_t < std::is_integral < T >::value, T > >
	inline T makeOdd ( const T i_ ) {

		return i_ | T ( 1 );
	}


	template < typename T, typename = std::enable_if_t < std::is_integral < T >::value, T > >
	inline T makeEven ( const T i_ ) {

		return i_ & ~T ( 1 );
	}


	// Random...

	// Seeding, from Intel Broadwell CPU onwards...

	void seed ( uint8_t  & s_ );
	void seed ( uint16_t & s_ );
	void seed ( uint32_t & s_ );
	void seed ( uint64_t & s_ );


	// http://xoroshiro.di.unimi.it/

	// XoRoShiRo128Plus: This is the successor to xorshift128+. It is the
	// fastest full-period generator passing BigCrush without systematic
	// failures, but due to the relatively short period it is acceptable only
	// for applications with a mild amount of parallelism; otherwise, use a
	// xorshift1024* generator.
	//
	// Beside passing BigCrush, this generator passes the PractRand test suite
	// up to (and included) 16TB, with the exception of binary rank tests,
	// which fail due to the lowest bit being an LFSR; all other bits pass all
	// tests. We suggest to use a sign test to extract a random Boolean value.
	//
	// Note that the generator uses a simulated rotate operation, which most C
	// compilers will turn into a single instruction. In Java, you can use
	// Long.rotateLeft(). In languages that do not make low-level rotation
	// instructions accessible xorshift128+ could be faster.
	//
	// The state must be seeded so that it is not everywhere zero. If you have
	// a 64-bit seed, we suggest to seed a splitmix64 generator and use its
	// output to fill s.
	//
	// Period 2^128 - 1.

	class XoRoShiRo128Plus {

		mutable uint64_t m_s0, m_s1;

	public:

		typedef uint64_t result_type;

		constexpr static result_type min ( ) {

			return result_type ( 0 );
		}

		constexpr static result_type max ( ) {

			return ~result_type ( 0 );
		}

		XoRoShiRo128Plus ( );
		XoRoShiRo128Plus ( const uint64_t s_ );

		void seed ( const uint64_t s_ );

		result_type operator ( ) ( ) const {

			const uint64_t r = m_s0 + m_s1;

			m_s1 ^= m_s0;

			m_s0 = _rotl64 ( m_s0, 55 );

			m_s0 ^= m_s1;
			m_s0 ^= m_s1 << 14; // a, b

			m_s1 = _rotl64 ( m_s1, 36 ); // c

			return r;
		}

		// This is the jump function for the generator. It is equivalent
		// to 2^64 calls to next(); it can be used to generate 2^64
		// non-overlapping subsequences for parallel computations.

		void jump ( ) const;
	};

#ifdef __AVX2__

	class XoRoShiRo128PlusAvx {

		mutable __declspec ( align ( 32 ) ) __m256i m_s0, m_s1, m_r;
		mutable uint32_t m_i = 3;

	public:

		typedef uint64_t result_type;

		static result_type min ( ) {

			return result_type ( 0 );
		}

		static result_type max ( ) {

			return ~result_type ( 0 );
		}

		XoRoShiRo128PlusAvx ( );
		XoRoShiRo128PlusAvx ( const uint64_t s_ );

		void seed ( const uint64_t s_ );

		result_type operator ( ) ( ) const {

			if ( m_i < 3UL ) {

				return ( ( result_type * ) & m_r ) [ ++m_i ];
			}

			m_r = _mm256_add_epi64 ( m_s0, m_s1 );

			m_s1 = _mm256_xor_si256 ( m_s1, m_s0 ); // is s here
// 			m_s0 = _mm256_xor_si256 ( _mm256_xor_si256 ( m_s1, _mm256_xor_si256 ( _mm256_slli_epi64 ( m_s1, 55 ), _mm256_srli_epi64 ( m_s1, 9 ) ) ), _mm256_slli_epi64 ( m_s1, 14 ) );
			m_s0 = _mm256_xor_si256 ( _mm256_xor_si256 ( m_s1, _mm256_xor_si256 ( _mm256_slli_epi64 ( m_s1, 211 ), _mm256_srli_epi64 ( m_s1, 31 ) ) ), _mm256_slli_epi64 ( m_s1, 56 ) );

			m_s1 = _mm256_xor_si256 ( _mm256_slli_epi64 ( m_s1, 36 ), _mm256_srli_epi64 ( m_s1, 28 ) );

			return ( ( result_type * ) & m_r ) [ m_i = 0 ];
		}
	};

#endif

	template < typename T, typename = std::enable_if_t < std::is_unsigned < T >::value && std::is_integral < T >::value, T > >
	void printBits ( const T n ) {

		T i = T ( 1 ) << ( sizeof ( T ) * 8 - 1 );

		while ( i ) {

			putchar ( int ( ( n & i ) > 0 ) + int ( 48 ) );

			i >>= 1;
		}

		putchar ( '\n' );
	}
}
