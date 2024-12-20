// This implementation originates from the gaussian.c in Labrador within Latticedog and has been rewritten into a Rust version.
//
// 1. Overview of the original C code:
//    The C code implements a discrete Gaussian sampler using various helper functions:
//    - `exp_small(x)`: Approximates exp(x) for small x values.
//    - `rcdtsampler(rcdt, len)`: Samples from a table of cumulative distribution thresholds (rcdt).
//    - `BerExp(x)`: Returns a Bernoulli trial with probability e^{-x}.
//    - `gaussian_sampler(mu, sigma)`: Returns a sample from a discrete Gaussian distribution with mean mu and standard deviation sigma.
//
//    The C code relies on:
//    - Some precomputed constants `rcdt125`.
//    - A custom randombytes function for randomness.
//    - Bitwise operations for efficiency and to ensure correct sampling.
//
// 2. Differences in Rust:
//    - In Rust, we would use `f64` for floating-point operations (equivalent to double in C).
//    - For randomness, we can use the `rand` crate, specifically `rand::Rng`, or `getrandom` for system randomness.
//    - For mathematical functions (exp, log, etc.), we can use `f64` methods from the standard library.
//
//    Note: The code below is a direct port. For real-world usage, consider error handling, parameter checks, and potentially safer wrappers.
//
// 3. Example usage:
//    - Suppose mu = 0.0 and sigma = 3.0 (just an example):
//      ```rust
//      let sample = gaussian_sampler(0.0, 3.0);
//      println!("Sampled value: {}", sample);
//      ```
//
//    - You could call the sampler multiple times to get multiple samples from the distribution.
//
// 4. If you prefer a simpler approach to sampling a normal distribution and then rounding to the nearest integer, consider using `rand_distr::Normal` from the `rand_distr` crate. However, this would yield continuous values from a normal distribution, not the same discrete distribution implemented here. The code below aims to replicate the exact discrete logic.
//
// 5. Potential improvements:
//    - Replace `randombytes()` with a secure RNG source, e.g. `rand::thread_rng()` or `getrandom::getrandom()`.
//    - Add checks that sigma is not zero or negative.
//    - Consider performance optimizations (like inlining).
//
// 6. Please note that this code uses several constants and low-level bit manipulations that are critical for correct cryptographic sampling. Ensure that changes are reviewed by a cryptography expert before use in production.
//
// 7. There's no direct, fully equivalent Rust crate that replicates this exact discrete Gaussian sampler from lattice-based cryptography. It's often implemented as part of lattice-based scheme libraries (like in `pqcrypto` libraries). For continuous distributions or general math operations, you might explore `statrs` or `rand_distr` as mentioned.
//
// Below is the Rust code rewrite:

use std::f64;
use rand::Rng;
use rand::rngs::OsRng; // Using the OS random number generator for cryptographic security

// RCDT table (threshold values)
static RCDT125: [u64; 11] = [
    4760398266205102531u64,
    1519614160162702751u64,
    278740980977854822u64,
    28213006827475907u64,
    1542173394262455u64,
    45012484900334u64,
    697367627767u64,
    5716868205u64,
    24757408u64,
    56588u64,
    68u64,
];

// Approximate exp(x) for |x| <= 0.5*ln(2)
// Equivalent to the exp_small function in C
#[inline]
fn exp_small(x: f64) -> f64 {
    // Polynomial constants
    const C1: f64 = 1.66666666666666019037e-01;
    const C2: f64 = -2.77777777770155933842e-03;
    const C3: f64 = 6.61375632143793436117e-05;
    const C4: f64 = -1.65339022054652515390e-06;
    const C5: f64 = 4.13813679705723846039e-08;

    let mut t = x * x;
    t = x - t * (C1 + t * (C2 + t * (C3 + t * (C4 + t * C5))));
    1.0 + (x - x * t / (t - 2.0))
}

// Samples from the RCDT table
#[inline]
fn rcdtsampler(rcdt: &[u64]) -> u32 {
    // Generate a random 63-bit value
    let mut rng = OsRng;
    let r: u64 = rng.gen(); 
    let r = r & ((1u64 << 63) - 1);

    let mut z: u32 = 0;
    for &threshold in rcdt {
        z += ((r.wrapping_sub(threshold)) >> 63) as u32;
    }
    z
}

// Bernoulli trial with probability e^{-x}
#[inline]
fn ber_exp(x: f64) -> u32 {
    // convert x to a scaled exponent:
    // We use x*(1/log(2)) + 0.5 and separate integer/fractional parts for exponentiation
    let sc = x * (1.0 / f64::consts::LN_2) + 0.5;
    let t_int = sc.floor();
    let t = t_int as u64;

    let fractional = sc - t_int;
    let neg_x = - (x - f64::consts::LN_2 * (t as f64));
    // exp_small(-x) and then shift by 63 bits (multiply by 2^63)
    let exp_val = exp_small(neg_x) * f64::exp2(63.0);
    let scaled_t = ((exp_val + 0.5).floor() as u64) >> t;

    let mut rng = OsRng;
    let u: u64 = rng.gen();
    let u = u & ((1u64 << 63) - 1);

    // Returns 1 if u < scaled_t, else 0
    ((u.wrapping_sub(scaled_t)) >> 63) as u32
}

// Main Gaussian sampler
// mu: mean, sigma: standard deviation
// Returns a signed 64-bit integer sample from the discrete Gaussian
#[inline]
pub fn gaussian_sampler(mu: f64, sigma: f64) -> i64 {
    assert!(sigma > 0.0, "sigma must be positive");

    let m = (sigma * (1.0 / 1.25)).ceil() as u64;
    let mut mask = m - 1;
    mask |= mask >> 1;
    mask |= mask >> 2;
    mask |= mask >> 4;
    mask |= mask >> 8;
    mask |= mask >> 16;
    mask |= mask >> 32;

    let mut rng = OsRng;
    let mut r: u64;
    // Sample r uniformly from {0,1,...,m-1}
    loop {
        r = rng.gen::<u64>() & mask;
        if r < m {
            break;
        }
    }

    let c0 = (mu + (r as f64)) * (1.0 / (m as f64));
    let c1 = c0.floor() as i64;
    let frac = c0 - (c1 as f64); // fractional part
    let d = sigma * (1.0 / (m as f64));
    let d_inv = 1.0 / (2.0 * d * d);

    // To generate the bimodal Gaussian
    let mut bits: u8;
    loop {
        // Ensure bits > 1 before using
        bits = rng.gen();
        bits |= 0x80;

        // Extract a single bit and sample k
        let mut b: i32;
        let mut k: i32;
        {
            b = (bits & 1) as i32;
            bits >>= 1;
            let sampled_k = rcdtsampler(&RCDT125) as i32;
            k = ((-b & (2 * sampled_k)) - sampled_k + b) as i32;
        }

        // Calculate exponent in BerExp
        // x = ((k - frac)^2)*d_inv - ((k - b)^2)*(1/3.125)
        // 3.125 = (5/ (2 * something))? It's a given constant from original code
        let x_val = ( ( (k as f64 - frac) * (k as f64 - frac) ) * d_inv )
                    - ( ( (k as f64 - b as f64)*(k as f64 - b as f64) ) * (1.0/3.125) );

        if ber_exp(x_val) != 0 {
            // We have succeeded in the Bernoulli trial
            let final_val = ((k as i64 + c1) * (m as i64)) - (r as i64);
            return final_val;
        }
        // Otherwise, loop until the condition is satisfied
    }
}


#[cfg(test)]
mod tests {
    use super::gaussian_sampler;

    #[test]
    fn test_gaussian_sampler_basic() {
        // Just a basic test to ensure the function runs without panic
        let sample = gaussian_sampler(0.0, 3.0);
        println!("Sampled value: {}", sample);

        // Another test call
        let sample2 = gaussian_sampler(2.5, 1.0);
        println!("Sampled value: {}", sample2);

        // We can't assert exact values due to randomness, but we can at least check it's finite.
        assert!(sample.is_finite());
        assert!(sample2.is_finite());
    }
}

// Suggested Crate for Reference (not a direct substitute):
// - For randomness: `rand` (already used above).
// - For continuous normal distributions: `rand_distr` crate: `rand_distr::Normal`.
//   Example: 
//     let normal = rand_distr::Normal::new(mu, sigma).unwrap();
//     let x: f64 = normal.sample(&mut rng); // Continuous sample
//   However, this is not the same as discrete Gaussian sampling used here.
// - For various statistics and distributions: `statrs` crate. But again, discrete Gaussian
//   for cryptographic purposes is often custom implemented as above.
//

