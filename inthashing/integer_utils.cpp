
#include <intrin.h>
#include <immintrin.h>
#define _CRT_RAND_S
#include <math.h>

#include <boost/random/taus88.hpp>

#include "../inthashing/sprp32.h"
#include "../inthashing/sprp64.h"


#include "integer_utils.hpp"


namespace jimi {

	uint32_t modularMultiplicativeInverse ( const uint32_t a_ ) {

		// Given odd a, compute x such that a * x = 1 over 32 bits...

		const uint16_t b = ( uint16_t ) a_;				// low 16 bits of a
		uint16_t x = ( ( ( b + 2u ) & 4u ) << 1 ) + b;	// low  4 bits of inverse

		x = ( 2u - b * x ) * x;							// low  8 bits of inverse
		x = ( 2u - b * x ) * x;							// low 16 bits of inverse

		return ( 2u - a_ * x ) * x;						//     32 bits of inverse
	}


	uint64_t modularMultiplicativeInverse ( const uint64_t a_ ) {

		// Given odd a, compute x such that a * x = 1 over 64 bits...

		const uint32_t b = ( uint32_t ) a_;				// low 32 bits of a
		uint32_t x = ( ( ( b + 2u ) & 4u ) << 1 ) + b;	// low  4 bits of inverse

		x = ( 2u - b * x ) * x;							// low  8 bits of inverse
		x = ( 2u - b * x ) * x;							// low 16 bits of inverse
		x = ( 2u - b * x ) * x;							// low 32 bits of inverse

		return ( 2u - a_ * x ) * x;						//     64 bits of inverse
	}


	bool isPrime ( const uint32_t n_ ) {

		// Assumes n = odd...

		static const uint32_t bases [ 3 ] = { 2UL, 7UL, 61UL };

		return efficient_mr32 ( bases, 3, n_ ) == 1;
	}


	bool isPrime ( const uint64_t n_ ) {

		// Assumes n = odd...

		static const uint64_t bases [ 7 ] = { 2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL };

		return efficient_mr64 ( bases, 7, n_ ) == 1;
	}


	uint32_t popCount ( const uint8_t x_ ) {

		return ( uint32_t ) __popcnt ( ( uint32_t ) x_ );
	}


	uint32_t popCount ( const uint16_t x_ ) {

		return ( uint32_t ) __popcnt ( ( uint32_t ) x_ );
	}


	uint32_t popCount ( const uint32_t x_ ) {

		return ( uint32_t ) __popcnt ( x_ );
	}


	uint32_t popCount ( const uint64_t x_ ) {

		return ( uint32_t ) __popcnt64 ( x_ );
	}


	// Random...

	// Seeding, from Intel Broadwell CPU onwards...

	static inline uint32_t get_seedu32_ ( ) {

		uint32_t seed;

		rand_s ( & seed );

		return seed;
	}

	extern "C" uint64_t get_seedu64 ( ) {

		return ( ( ( uint64_t ) get_seedu32_ ( ) << 32 ) | ( uint64_t ) get_seedu32_ ( ) );
	}



	void seed ( uint8_t & s_ ) {

		uint32_t random_val;

		rand_s ( & random_val );

		s_ = ( uint8_t ) ( random_val >> 24 ) ^ ( uint8_t ) ( random_val >> 16 ) ^ ( uint8_t ) ( random_val >> 8 ) ^ ( uint8_t ) random_val;
	}


	void seed ( uint16_t & s_ ) {

		uint32_t random_val;

		rand_s ( &random_val );

		s_ = ( uint16_t ) ( random_val >> 16 ) ^ ( uint16_t ) random_val;

	}


	void seed ( uint32_t & s_ ) {

		rand_s ( & s_ );
	}


	void seed ( uint64_t & s_ ) {

		rand_s ( ( ( uint32_t * ) & s_ ) + 0 );
		rand_s ( ( ( uint32_t * ) & s_ ) + 1 );
	}



	void seed_bw ( uint16_t & s_ ) {

		_rdseed16_step ( &s_ );
	}


	void seed_bw ( uint32_t & s_ ) {

		_rdseed32_step ( & s_ );
	}


	void seed_bw ( uint64_t & s_ ) {

		_rdseed64_step ( & s_ );
	}


	XoRoShiRo128Plus::XoRoShiRo128Plus ( ) {

		m_s0 = [ ] ( ) -> result_type { result_type s; jimi::seed ( s ); return s; } ( );
		m_s1 = [ ] ( ) -> result_type { result_type s; jimi::seed ( s ); return s; } ( );
	}


	XoRoShiRo128Plus::XoRoShiRo128Plus ( const uint64_t s_ ) {

		seed ( s_ );
	}


	void XoRoShiRo128Plus::seed ( const uint64_t s_ ) {

		boost::random::taus88 rng ( s_ );

		m_s0 = boost::random::uniform_int_distribution < uint64_t > ( ) ( rng );
		m_s1 = boost::random::uniform_int_distribution < uint64_t > ( ) ( rng );
	}


	void XoRoShiRo128Plus::jump ( ) const {

		static const uint64_t j0 = 0xbeac0467eba5facb, j1 = 0xd86b048b86aa9922;

		uint64_t s0 = 0, s1 = 0;

		for ( uint32_t b = 0; b < 64UL; ++b ) {

			if ( j0 & 1ULL << b ) {

				s0 ^= m_s0;
				s1 ^= m_s1;
			}

			m_s1 = m_s0 ^ m_s1;
			m_s0 = _rotl64 ( m_s0, 55 ) ^ m_s1 ^ ( m_s1 << 14 ); // a, b
			m_s1 = _rotl64 ( m_s1, 36 ); // c
		}

		for ( uint32_t b = 0; b < 64UL; ++b ) {

			if ( j1 & 1ULL << b ) {

				s0 ^= m_s0;
				s1 ^= m_s1;
			}

			m_s1 = m_s0 ^ m_s1;
			m_s0 = _rotl64 ( m_s0, 55 ) ^ m_s1 ^ ( m_s1 << 14 ); // a, b
			m_s1 = _rotl64 ( m_s1, 36 ); // c
		}

		m_s0 = s0;
		m_s1 = s1;
	}


#ifdef __AVX2__


	XoRoShiRo128PlusAvx::XoRoShiRo128PlusAvx ( ) {

		m_s0 = _mm256_set_epi64x ( [ ] ( ) -> result_type { result_type s; jimi::seed ( s ); return s; } ( ), [ ] ( ) -> result_type { result_type s; jimi::seed ( s ); return s; } ( ), [ ] ( ) -> result_type { result_type s; jimi::seed ( s ); return s; } ( ), [ ] ( ) -> result_type { result_type s; jimi::seed ( s ); return s; } ( ) );
		m_s1 = _mm256_set_epi64x ( [ ] ( ) -> result_type { result_type s; jimi::seed ( s ); return s; } ( ), [ ] ( ) -> result_type { result_type s; jimi::seed ( s ); return s; } ( ), [ ] ( ) -> result_type { result_type s; jimi::seed ( s ); return s; } ( ), [ ] ( ) -> result_type { result_type s; jimi::seed ( s ); return s; } ( ) );
	}


	XoRoShiRo128PlusAvx::XoRoShiRo128PlusAvx ( const uint64_t s_ ) {

		seed ( s_ );
	}


	void XoRoShiRo128PlusAvx::seed ( const uint64_t s_ ) {

		boost::random::taus88 rng ( s_ );

		m_s0 = _mm256_set_epi64x ( boost::random::uniform_int_distribution < int64_t > ( ) ( rng ), boost::random::uniform_int_distribution < int64_t > ( ) ( rng ), boost::random::uniform_int_distribution < int64_t > ( ) ( rng ), boost::random::uniform_int_distribution < int64_t > ( ) ( rng ) );
		m_s1 = _mm256_set_epi64x ( boost::random::uniform_int_distribution < int64_t > ( ) ( rng ), boost::random::uniform_int_distribution < int64_t > ( ) ( rng ), boost::random::uniform_int_distribution < int64_t > ( ) ( rng ), boost::random::uniform_int_distribution < int64_t > ( ) ( rng ) );
	}

#endif

}
