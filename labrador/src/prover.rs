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
}

// Function to calculate b^(k)
fn calculate_b_constraint(s: &Vec<Vec<RingPolynomial>>, a: &Vec<Vec<usize>>, phi: &Vec<usize>) -> usize {
    let mut b = 0;
    let s_len = s.len();
    // Calculate b^(k)
    for i in 0..s_len {
        for j in 0..s_len {
            // calculate inner product of s[i] and s[j]
            let inner_product = s[i].iter().map(|elem| elem.coefficients[0]).zip(s[j].iter().map(|elem| elem.coefficients[0])).map(|(x, y)| {
                // println!("x: {}, y: {}", x, y);
                x * y
            }).sum::<usize>();
            b += a[i][j] * inner_product;
        }
        // calculate inner product of s[i] and phi
        let inner_product_phi = s[i].iter().map(|elem| elem.coefficients[0]).zip(phi.iter()).map(|(x, y)| x * y).sum::<usize>();
        b += inner_product_phi;
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
        for i in 0..k {
            let b_i = calculate_b_constraint(&s, &a_k, &phi_k);
            b_values_k.push(b_i);
            println!("b^({}) = {}", i, b_i);
        }

        let l: usize = 4; // Define L as usize
        // Generate random a^(k)_{i,j} and φ^{(k)}_{i}
        let a_l: Vec<Vec<usize>> = (0..s_len).map(|_| (0..s_len).map(|_| rng.gen_range(1..l)).collect()).collect();
        println!("a_l: {:?}", a_l);
        let phi_l: Vec<usize> = (0..s_len).map(|_| rng.gen_range(1..5)).collect();
        // calculate b^(l)
        let mut b_values_l = Vec::new();
        for i in 0..l {
            let b_i = calculate_b_constraint(&s, &a_l, &phi_l);
            b_values_l.push(b_i);
            println!("b^({}) = {}", i, b_i);
        }

        let size_kappa = 3; // Example size
        let size_n = 5;
        // A: matrix size: kappa * n, each element is RingPolynomial(Rq)
        // calculate t_i = A * s_i for all i = 1..r
        // size of t_i = (kappa * n)Rq * 1Rq = kappa * n
        let matrix_a = A::new(size_kappa, size_n);
        println!("A: {:?}", matrix_a);
        // print size of A
        println!("size of A: {:?} x {:?}", matrix_a.values.len(), matrix_a.values[0].len());
        let mut all_t_i = Vec::new();
        for s_i in &s {
            let t_i = calculate_a_times_s_i(&matrix_a, &s_i);
            println!("size of t_i: {:?}", t_i.len());
            all_t_i.push(t_i);
        }
        println!("Calculated all t_i: {:?}", all_t_i);
        // print size of t_i
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
        println!("b_k: {}", b_k);
        assert_eq!(b_k, 1983);
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
}
