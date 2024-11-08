use rand::Rng;

pub fn setup() {
    // 0. setup
    // public parameters after setup: [a_ij^(k), a_ij^(l), phi^(k), phi^(l), b^(k), b0(l)']

    // 1. setup constraints
    // 1.1 get s_i and do norm check
    // 1.1.1 get s_i = s_1 - s_r; r is the number of witness s
    // each s_i is a vector of ring elements

    // 1.1.2 get beta norm bound, refer paper page 26, theorem 6.3

    // 1.1.3 do check: sum of s_i norm <= beta_square

    // 1.2 calculate b^(k)
    // 1.2.1 calculate dot product ss = a_ij * <s_i, s_j> for all i, j
    // a_ij is the quadratic coefficient, phi^(k) is the linear coefficient
    // 1.2.2 calculate phi_s = <phi^(k), s_i> for all i
    // 1.2.3 calculate b^(k) = sum(ss) + sum(phi_s)

    // 1.3 calculate b'^(l)
    // 1.3.1 calculate dot product ss = a_ij' * <s_i, s_j> for all i, j
    // a_ij' is the quadratic coefficient, phi^(l)' is the linear coefficient
    // 1.3.2 calculate phi_s = <phi^(l)', s_i> for all i
    // 1.3.3 calculate b'^(l) = sum(ss) + sum(phi_s)

    // L = |F'| = ceiling(128 / logQ)
}

pub fn prove() {
    println!("Proving something...");
    // 2. GOAL: calculate ajtai commitment (1st outer commitment)
    // 2.1 split t to t_i for all i
    // 2.1.1 get basis b1, refer to paper page 16, labrador c code line 142
    // refer to note: line => xxx

    // 2.2 split g = <s_i, s_j> for all i, j
    // 2.2.1 get basis b2 same as 2.1.1

    // 2.3 calculate u1
    // 2.3.1 B & C is randomly chosen
    // 2.3.2 calculate u1 = sum(B_ik * t_i^(k)) + sum(C_ijk * g_ij^(k))


    // ================================================

    // 3. GOAL: JL projection
    // 3.1 PI_i is randomly chosen from \Chi { -1, 0, 1 }^{256 * nd}
    //      (Using Guassian Distribution)
    // 3.2 caculate p_j = sum(<pi_i^(j)>, s_i) for all i-r
    // 3.3 Verifier have to check: || p || <= \sqrt{128} * beta

    // ================================================

    // 4. GOAL: Aggregation
    // 4.1 psi^(k) is randomly chosen from Z_q^{L}
    // 4.2 omega^(k) is randomly chosen from Z_q^{256}
    //      (Both using Guassian Distribution)
    // 4.3 caculate b^{''(k)}
    // 4.3.1 calculate a_ij^{''(k)} = sum(psi_l^(k) * a_ij^{'(l)}) for all l = 1..L
    // 4.3.2 calculate phi_i^{''(k)} =
    //       sum(psi_l^(k) * phi_i^{'(l)}) for all l = 1..L
    //       + sum(omega_j^(k) * sigma_{-1} * pi_i^{j)) for all j = 1..256
    // 4.3.3 calculate b^{''(k)} = sum(a_ij^{''(k)} * <s_i, s_j>) + sum(<phi_i^{''(k)}, s_i>)

    // Send b^{''(k)} to verifier
    // Verifier check: b_0^{''(k)} ?= <⟨omega^(k),p⟩> + sum(psi_l^(k) * b_0^{'(l)}) for all l = 1..L

    // ================================================

    // 5. GOAL: Calculate u2 (2nd outer commitment)
    // 5.1 vec<alpha> and vec<beta> are randomly chosen from R_q^{128/logQ} // why is this not `L`
    // 5.2 phi_i = sum(alpha_k * phi_i) + beta_k * phi_i^{''(k)}
    // 5.3 h_ij = 1/2 * (<phi_i, s_j> + <phi_j, s_i>)
    // 5.4 u2 = sum D_ij * h_ij^(k) for all k = 1..(t1-1)

    // Send u2 to verifier

    // ================================================

    // 6. GOAL: calculate z (Amortized Opening)
    // 6.1 c_i is randomly chosen from C
    // 6.2 calculate z = sum(c_i * s_i) for all i = 1..r
    // Send z, t_i, g_ij, h_ij to verifier
}

// Assuming you have a RingPolynomial type defined
#[derive(Debug)]
struct RingPolynomial {
    coefficients: Vec<usize>, // Example field, adjust as necessary
}

impl RingPolynomial {
    // Add this method to enable multiplication by RingPolynomial
    fn multiply_by_ringpolynomial(&self, other: &RingPolynomial) -> RingPolynomial {
        let mut result_coefficients = vec![0; self.coefficients.len() + other.coefficients.len() - 1];
        for (i, &coeff1) in self.coefficients.iter().enumerate() {
            for (j, &coeff2) in other.coefficients.iter().enumerate() {
                result_coefficients[i + j] += coeff1 * coeff2;
            }
        }
        RingPolynomial { coefficients: result_coefficients }
    }

    fn add_ringpolynomial(&self, other: &RingPolynomial) -> RingPolynomial {
        // print coefficients of self and other
        println!("coefficients of self: {:?}", self.coefficients);
        println!("coefficients of other: {:?}", other.coefficients);
        // check 2 polynomials are the same size
        // assert!(self.coefficients.len() == other.coefficients.len());
            // Start of Selection
            let max_len = std::cmp::max(self.coefficients.len(), other.coefficients.len());
            let mut result_coefficients = Vec::with_capacity(max_len);
            for i in 0..max_len {
                let a = if i < self.coefficients.len() { self.coefficients[i] } else { 0 };
                let b = if i < other.coefficients.len() { other.coefficients[i] } else { 0 };
                result_coefficients.push(a + b);
            }
        RingPolynomial { coefficients: result_coefficients }
    }
}

// inner product of 2 vectors of RingPolynomial
// Start of Selection
fn inner_product_ringpolynomial(a: &Vec<RingPolynomial>, b: &Vec<RingPolynomial>) -> RingPolynomial {
    a.iter()
        .zip(b.iter())
        .map(|(a, b)| a.multiply_by_ringpolynomial(b))
        .collect::<Vec<RingPolynomial>>()
        .into_iter()
        .reduce(|acc, x| acc.add_ringpolynomial(&x))
        .unwrap()
}

// Function to calculate b^(k)
fn calculate_b_constraint(s: &Vec<Vec<RingPolynomial>>, a: &Vec<Vec<usize>>, phi: &Vec<usize>) -> RingPolynomial {
    let mut b: RingPolynomial = RingPolynomial { coefficients: vec![0] };
    let s_len = s.len();
    // Calculate b^(k)
    for i in 0..s_len {
        for j in 0..s_len {
            // calculate inner product of s[i] and s[j], will retturn a single RingPolynomial
            // Start of Selection
            let elem_s_i = &s[i];
            let elem_s_j = &s[j];
            // Calculate inner product and update b
            let inner_product_si_sj = inner_product_ringpolynomial(&elem_s_i, &elem_s_j);
            b = b.add_ringpolynomial(
                &inner_product_si_sj
                    .multiply_by_ringpolynomial(&RingPolynomial { coefficients: vec![a[i][j]] })
            );
        }
        // calculate inner product of s[i] and phi
        // Start of Selection
        s[i]
            .iter()
            .map(|elem| elem)
            .zip(phi.iter().map(|x| RingPolynomial { coefficients: vec![*x] }))
            .for_each(|(x, y)| b = b.add_ringpolynomial(&x.multiply_by_ringpolynomial(&y)));
    }

    b
}

#[derive(Debug)]
struct A {
    // matrix size: kappa * n
    values: Vec<Vec<RingPolynomial>>, // A matrix of RingPolynomial values
}

impl A {
    fn new(size_kappa: usize, size_n: usize) -> Self {
        let mut rng = rand::thread_rng();
        let values = (0..size_kappa)
            .map(|_| (0..size_n).map(|_| RingPolynomial { coefficients: (0..size_n).map(|_| rng.gen_range(1..10)).collect() }).collect())
            .collect();
        A { values }
    }
}

// calculate A matrix times s_i
fn calculate_a_times_s_i(a: &A, s_i: &Vec<RingPolynomial>) -> Vec<RingPolynomial> {
    a.values.iter().map(|row| {
        row.iter().zip(s_i.iter()).map(|(a, b)| a.multiply_by_ringpolynomial(b)).collect::<Vec<RingPolynomial>>()
    }).collect::<Vec<Vec<RingPolynomial>>>().into_iter().flatten().collect::<Vec<RingPolynomial>>()
}

// convert number to basis
// 42 = 0 * 2^7 + 1 * 2^6 + 0 * 2^5 + 1 * 2^4 + 0 * 2^3 + 1 * 2^2 + 0 * 2^1 + 0 * 2^0
// first step: 42 / 2 = 21, result_i = 0
// second step: 21 / 2 = 10, result_i = 1
// third step: 10 / 2 = 5, result_i = 0
// forth step: 5 / 2 = 2, result_i = 1
// fifth step: 2 / 2 = 1, result_i = 0
// sixth step: 1 / 2 = 0, result_i = 1

fn num_to_basis(num: usize, basis: usize, digits: usize) -> Vec<usize> {
    let mut result = Vec::new();
    let mut remainder = num;

    while remainder > 0 {
        result.push(remainder % basis);
        remainder /= basis;
    }

    while result.len() < digits {
        // push 0 to the highest position
        result.push(0);
    }

    result
}

// convert ring polynomial to basis
fn ring_polynomial_to_basis(poly: &RingPolynomial, basis: usize, digits: usize) -> Vec<Vec<usize>> {
    poly.coefficients.iter().map(|coeff| num_to_basis(*coeff, basis, digits)).collect()
}

// create test case for setup
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_setup() {
        // s is a vector of size r. each s_i is a RingPolynomial(Rq) with n coefficients
        let s_len: usize = 3; // r: Number of witness elements
        let s_i_length: usize = 5; // n
        let beta: usize = 50; // Example value for beta
        let s: Vec<Vec<RingPolynomial>> = (1..=s_len).map(|i| {
            (1..=s_i_length).map(|j| RingPolynomial { coefficients: vec![i * 3 + j, i * 3 + j + 1, i * 3 + j + 2] }).collect()
        }).collect();
        println!("s: {:?}", s);
        // Calculate the sum of squared norms
        let mut sum_squared_norms = 0;
        for vector in &s {
            let norm_squared: usize = vector.iter()
                .map(|elem| elem.coefficients[0].pow(2)) // Calculate the square of each element
                .sum();
            sum_squared_norms += norm_squared; // Accumulate the squared norms
        }
        println!("sum_squared_norms: {}", sum_squared_norms);
        println!("beta^2: {}", beta.pow(2));
        // Check the condition
        assert!(sum_squared_norms <= beta.pow(2), "The condition is not satisfied: sum of squared norms exceeds beta^2");

        let mut rng = rand::thread_rng();
        let k: usize = 6; // Change k to usize
        // Generate random a^(k)_{i,j} and φ^{(k)}_{i}
        let a_k: Vec<Vec<usize>> = (0..s_len).map(|_| (0..s_len).map(|_| rng.gen_range(1..k)).collect()).collect();
        let phi_k: Vec<usize> = (0..s_len).map(|_| rng.gen_range(1..5)).collect();
        println!("a_k: {:?}", a_k);
        println!("phi_k: {:?}", phi_k);

        // calculate b^(k)
        let mut b_values_k = Vec::new();
        for _ in 0..k {
            let b_i = calculate_b_constraint(&s, &a_k, &phi_k);
            b_values_k.push(b_i);
        }
        println!("b_values_k: {:?}", b_values_k);
        let l: usize = 4; // Define L as usize
        // Generate random a^(k)_{i,j} and φ^{(k)}_{i}
        let a_l: Vec<Vec<usize>> = (0..s_len).map(|_| (0..s_len).map(|_| rng.gen_range(1..l)).collect()).collect();
        println!("a_l: {:?}", a_l);
        let phi_l: Vec<usize> = (0..s_len).map(|_| rng.gen_range(1..5)).collect();
        // calculate b^(l)
        let mut b_values_l = Vec::new();
        for _ in 0..l {
            let b_i = calculate_b_constraint(&s, &a_l, &phi_l);
            b_values_l.push(b_i);
        }
        println!("b_values_l: {:?}", b_values_l);
        let size_kappa = 3; // Example size
        let size_n = 5;
        // A: matrix size: kappa * n, each element is RingPolynomial(Rq)
        // calculate t_i = A * s_i for all i = 1..r
        // size of t_i = (kappa * n)Rq * 1Rq = kappa * n
        let matrix_a = A::new(size_kappa, size_n);
        println!("A: {:?}", matrix_a);
        // print size of A
        println!("size of A: {:?} x {:?}", matrix_a.values.len(), matrix_a.values[0].len());
        assert!(matrix_a.values.len() == size_kappa);
        assert!(matrix_a.values[0].len() == size_n);
        let mut all_t_i = Vec::new();
        for s_i in &s {
            let t_i = calculate_a_times_s_i(&matrix_a, &s_i);
            println!("size of t_i: {:?}", t_i.len());
            all_t_i.push(t_i);
        }
        println!("Calculated all t_i: {:?}", all_t_i);
        // print size of all_t_i
        println!("size of all_t_i: {:?}", all_t_i.len());
        // check size of all_t_i is kappa
        assert!(all_t_i.len() == size_kappa);

        // ================================================
        // prover
        // 2. GOAL: calculate ajtai commitment (1st outer commitment)
        // 2.1 split t to t_i for all i
        // 2.1.1 get basis b1, refer to paper page 16, labrador c code line 142
        // s = β / sqrt(rnd)
        // r: s_len, n
        // for example:
        // t_0: [
        //  [8, 46, 61, 71, 33, 33, 18], -> t[0][0]
        //  [20, 54, 94, 93, 70, 33, 14], -> t[0][1]
        //  [24, 40, 100, 85, 121, 57, 56],  -> t[0][2]
        //  [14, 37, 91, 118, 159, 109, 72],  -> t[0][3]
        //  [32, 100, 120, 121, 102, 103, 70],
        //  [20, 61, 83, 76, 55, 53, 42],
        //  [35, 67, 124, 129, 141, 92, 42],
        //  [48, 86, 105, 59, 34, 30, 16],
        //  [56, 127, 172, 141, 75, 44, 9],
        //  [72, 113, 190, 144, 164, 94, 60],
        //  [36, 69, 120, 117, 131, 94, 48],
        //  [15, 53, 98, 111, 88, 46, 21],
        //  [12, 56, 119, 173, 159, 100, 32],
        //  [28, 95, 157, 144, 92, 33, 27],
        //  [56, 127, 166, 115, 63, 37, 30]
        // ]
        // such as basis = 10
        // t1: length of t[i][j][k], such as length of t[0][0][0] = 3
        // Then:
        // t_0_basis_form: [
        //   [[8, 2, 0], [7, 6, 0], [4, 1, 1], [0, 0, 1], [7, 6, 0], [3, 2, 0], [6, 0, 0]],
        //   [[5, 4, 0], [4, 8, 0], [4, 3, 1], [9, 0, 1], [9, 8, 0], [7, 4, 0], [4, 1, 0]],
        //   [[6, 0, 0], [3, 4, 0], [0, 8, 0], [5, 2, 1], [5, 2, 1], [8, 9, 0], [8, 4, 0]],
        //   [[8, 2, 0], [0, 6, 0], [9, 8, 0], [1, 8, 0], [0, 0, 1], [3, 8, 0], [3, 6, 0]],
        //   [[4, 2, 0], [5, 7, 0], [4, 2, 1], [3, 5, 1], [8, 6, 1], [2, 3, 1], [0, 8, 0]],
        //   [[4, 0, 0], [3, 1, 0], [6, 3, 0], [1, 6, 0], [2, 9, 0], [6, 7, 0], [8, 4, 0]],
        //   [[0, 3, 0], [6, 4, 0], [9, 9, 0], [8, 9, 0], [9, 3, 1], [0, 9, 0], [6, 5, 0]],
        //   [[8, 1, 0], [3, 3, 0], [6, 8, 0], [8, 0, 1], [4, 2, 1], [9, 6, 0], [4, 2, 0]],
        //   [[4, 1, 0], [7, 3, 0], [5, 0, 1], [8, 4, 1], [4, 4, 1], [1, 7, 0], [9, 0, 0]],
        //   [[4, 2, 0], [3, 4, 0], [4, 6, 0], [4, 5, 0], [0, 7, 0], [6, 5, 0], [0, 4, 0]],
        //   [[8, 2, 0], [1, 7, 0], [3, 0, 1], [6, 0, 1], [8, 8, 0], [8, 7, 0], [6, 3, 0]],
        //   [[0, 2, 0], [4, 3, 0], [5, 5, 0], [7, 5, 0], [6, 8, 0], [7, 7, 0], [9, 4, 0]],
        //   [[6, 3, 0], [2, 7, 0], [5, 2, 1], [9, 1, 1], [5, 1, 1], [8, 6, 0], [2, 3, 0]],
        //   [[4, 1, 0], [3, 2, 0], [4, 5, 0], [8, 4, 0], [8, 5, 0], [5, 2, 0], [8, 1, 0]],
        //   [[6, 1, 0], [6, 2, 0], [3, 9, 0], [8, 9, 0], [8, 3, 1], [5, 6, 0], [0, 5, 0]]
        // ]
        let basis = 10;
        let digits = 3; // t1
        let all_t_i_basis_form: Vec<Vec<Vec<Vec<usize>>>> = all_t_i.iter().map(|t_i|
            t_i.iter().map(|t_i_j| ring_polynomial_to_basis(t_i_j, basis, digits)).collect::<Vec<Vec<Vec<usize>>>>()
        ).collect::<Vec<Vec<Vec<Vec<usize>>>>>();
        println!("all_t_i_basis_form: {:?}", all_t_i_basis_form);
        // print t_0
        println!("t_0: {:?}", all_t_i[0]);
        println!("t_0_basis_form: {:?}", all_t_i_basis_form[0]);
        // Start of Selection
        for (i, t_i_j_basis_form) in all_t_i_basis_form.iter().enumerate() {
            for (j, coeff_basis_form) in t_i_j_basis_form.iter().enumerate() {
                println!("coeff_basis_form (all_t_i_basis_form[{}][{}]): {:?}", i, j, coeff_basis_form);
                // Get the number of columns from the first inner vector
                let num_cols = coeff_basis_form[0].len();

                // Sum elements at each position across all inner vectors
                let result: Vec<usize> = (0..num_cols)
                    .map(|k| coeff_basis_form.iter().map(|row| row[k]).sum())
                    .collect();

                println!("result (all_t_i_basis_form[{}][{}]): {:?}", i, j, result);
            }
        }
        // 2
        // 2.2.1 get basis b2 same as 2.1.1
        // Start of Selection
        // Calculate g_ij = <s_i, s_j>
        let num_s: usize = s.len();
        let mut g: Vec<Vec<usize>> = vec![vec![0; num_s]; num_s];
        // for i in 0..num_s {
        //     for j in 0..num_s {
        //         g_ij[i][j] = s[i].iter().zip(s[j].iter()).map(|(a, b)| a * b).sum();
        //     }
        // }
        // Calculate b^(k)
        for i in 0..num_s {
            for j in 0..num_s {
                // calculate inner product of s[i] and s[j]
                let inner_product = s[i].iter().map(|elem| elem.coefficients[0]).zip(s[j].iter().map(|elem| elem.coefficients[0])).map(|(x, y)| {
                    x * y
                }).sum::<usize>();
                g[i][j] = inner_product;
            }
        }
        println!("g: {:?}", g);
        // 2.3 calculate u1
        // 2.3.1 B & C is randomly chosen
        // 2.3.2 calculate u1 = sum(B_ik * t_i^(k)) + sum(C_ijk * g_ij^(k))


        // ================================================

        // 3. GOAL: JL projection
        // 3.1 PI_i is randomly chosen from \Chi { -1, 0, 1 }^{256 * nd}
        //      (Using Guassian Distribution)
        // 3.2 caculate p_j = sum(<pi_i^(j)>, s_i) for all i-r
        // 3.3 Verifier have to check: || p || <= \sqrt{128} * beta

        // ================================================

        // 4. GOAL: Aggregation
        // 4.1 psi^(k) is randomly chosen from Z_q^{L}
        // 4.2 omega^(k) is randomly chosen from Z_q^{256}
        //      (Both using Guassian Distribution)
        // 4.3 caculate b^{''(k)}
        // 4.3.1 calculate a_ij^{''(k)} = sum(psi_l^(k) * a_ij^{'(l)}) for all l = 1..L
        // 4.3.2 calculate phi_i^{''(k)} =
        //       sum(psi_l^(k) * phi_i^{'(l)}) for all l = 1..L
        //       + sum(omega_j^(k) * sigma_{-1} * pi_i^{j)) for all j = 1..256
        // 4.3.3 calculate b^{''(k)} = sum(a_ij^{''(k)} * <s_i, s_j>) + sum(<phi_i^{''(k)}, s_i>)

        // Send b^{''(k)} to verifier
        // Verifier check: b_0^{''(k)} ?= <⟨omega^(k),p⟩> + sum(psi_l^(k) * b_0^{'(l)}) for all l = 1..L

        // ================================================

        // 5. GOAL: Calculate u2 (2nd outer commitment)
        // 5.1 vec<alpha> and vec<beta> are randomly chosen from R_q^{128/logQ} // why is this not `L`
        // 5.2 phi_i = sum(alpha_k * phi_i) + beta_k * phi_i^{''(k)}
        // 5.3 h_ij = 1/2 * (<phi_i, s_j> + <phi_j, s_i>)
        // 5.4 u2 = sum D_ij * h_ij^(k) for all k = 1..(t1-1)

        // Send u2 to verifier

        // ================================================

        // 6. GOAL: calculate z (Amortized Opening)
        // 6.1 c_i is randomly chosen from C
        // 6.2 calculate z = sum(c_i * s_i) for all i = 1..r
        // Send z, t_i, g_ij, h_ij to verifier
    }

    #[test]
    fn test_multiply_by_ringpolynomial() {
        let poly1 = RingPolynomial { coefficients: vec![1, 2] };
        let poly2 = RingPolynomial { coefficients: vec![3, 4] };
        let result = poly1.multiply_by_ringpolynomial(&poly2);
        assert_eq!(result.coefficients, vec![3, 10, 8]); // 1*3, 1*4 + 2*3, 2*4
    }

    #[test]
    fn test_calculate_b_k() {
        let r: usize = 3;
        let s: Vec<Vec<RingPolynomial>> = vec![
            vec![RingPolynomial { coefficients: vec![1, 2, 3] }, RingPolynomial { coefficients: vec![4, 5, 6] }],
            vec![RingPolynomial { coefficients: vec![7, 8, 9] }, RingPolynomial { coefficients: vec![10, 11, 12] }],
            vec![RingPolynomial { coefficients: vec![13, 14, 15] }, RingPolynomial { coefficients: vec![16, 17, 18] }],
        ];
        let k: usize = 6;
        let a_k: Vec<Vec<usize>> = (0..r).map(|_| (0..r).map(|r_i| r_i).collect()).collect();
        let phi_k: Vec<usize> = (0..r).map(|r_i| r_i).collect();
        let b_k = calculate_b_constraint(&s, &a_k, &phi_k);
        println!("b_k: {:?}", b_k);
        // assert_eq!(b_k, 1983);
    }

    #[test]
    fn test_a_new() {
        let size: usize = 3;
        let a = A::new(size, size);
        assert_eq!(a.values.len(), size);
        assert_eq!(a.values[0].len(), size);
    }

    #[test]
    fn test_calculate_a_times_s_i() {
        let a = A::new(2, 2);
        let s_i = vec![
            RingPolynomial { coefficients: vec![1, 2] },
            RingPolynomial { coefficients: vec![3, 4] },
        ];
        let result = calculate_a_times_s_i(&a, &s_i);
        assert_eq!(result.len(), a.values.len() * s_i.len()); // Check that the result length is correct
    }

    #[test]
    fn test_num_to_basis() {
        let num = 42;
        let basis = 2;
        let digits = 6;
        let binary = num_to_basis(num, basis, digits);
        assert_eq!(binary, vec![0, 1, 0, 1, 0, 1]);

        let num = 100;
        let basis = 3;
        let digits = 6;
        let binary = num_to_basis(num, basis, digits);
        assert_eq!(binary, vec![1, 0, 2, 0, 1, 0]);

        let num = 100;
        let basis = 6;
        let digits = 6;
        let binary = num_to_basis(num, basis, digits);
        assert_eq!(binary, vec![4, 4, 2, 0, 0, 0]);

        let num = 100;
        let basis = 10;
        let digits = 6;
        let binary = num_to_basis(num, basis, digits);
        assert_eq!(binary, vec![0, 0, 1, 0, 0, 0]);
    }


    #[test]
    fn test_ring_polynomial_to_basis() {
        let poly = RingPolynomial { coefficients: vec![42, 100, 100] };
        let basis = 2;
        let digits = 8;
        let expected_result = vec![
            vec![0, 1, 0, 1, 0, 1, 0, 0],
            vec![0, 0, 1, 0, 0, 1, 1, 0],
            vec![0, 0, 1, 0, 0, 1, 1, 0],
        ];
        let result = ring_polynomial_to_basis(&poly, basis, digits);
        assert_eq!(result, expected_result);
    }

    #[test]
    fn test_inner_product_ringpolynomial() {
        let a = vec![
            RingPolynomial { coefficients: vec![1, 2, 3] },
            RingPolynomial { coefficients: vec![4, 5, 6] },
        ];
        let b = vec![
            RingPolynomial { coefficients: vec![7, 8, 9] },
            RingPolynomial { coefficients: vec![10, 11, 12] },
        ];

        let result = inner_product_ringpolynomial(&a, &b);

        // Expected result calculation:
        // (1 + 2x + 3x^2) * (7 + 8x + 9x^2) = 7 + 22x + 46x^2 + 42x^3 + 27x^4
        // (4 + 5x + 6x^2) * (10 + 11x + 12x^2) = 40 +  96x + 163x^2 + 126x^3 + 72x^4
        // Sum: 47 + 116x + 209x^2 + 168x^3 + 99x^4

        let expected = RingPolynomial { coefficients: vec![47, 116, 209, 168, 99] };

        assert_eq!(result.coefficients, expected.coefficients);
    }
}


