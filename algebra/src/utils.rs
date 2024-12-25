use crate::polynomial_ring::PolynomialRing;
use crate::zq::Zq;
use rand::Rng;

pub fn zero_poly() -> PolynomialRing {
    PolynomialRing {
        coefficients: vec![Zq::from(0); 1],
    }
}

// a: Vec<PolynomialRing>, b: PolynomialRing
// calculate c = a * b, c_i = a_i * b
// c: Vec<PolynomialRing>
pub fn poly_vec_times_poly(a: &[PolynomialRing], b: &PolynomialRing) -> Vec<PolynomialRing> {
    a.iter().map(|a_i| a_i * b).collect()
}
// a: Vec<PolynomialRing>, b: Vec<PolynomialRing>
// calculate c = a + b, c_i = a_i + b_i
// c: Vec<PolynomialRing>
pub fn poly_vec_add_poly_vec(a: &[PolynomialRing], b: &[PolynomialRing]) -> Vec<PolynomialRing> {
    a.iter().zip(b.iter()).map(|(a_i, b_i)| a_i + b_i).collect()
}

#[allow(clippy::ptr_arg)]
// inner product of 2 vectors of PolynomialRing
pub fn inner_product_polynomial_ring_vector(
    a: &Vec<PolynomialRing>,
    b: &Vec<PolynomialRing>,
) -> PolynomialRing {
    assert_eq!(
        a.len(),
        b.len(),
        "inner_product_polynomial_ring_vector: a and b must have the same length"
    );
    a.iter()
        .zip(b.iter())
        .map(|(a, b)| a * b)
        .collect::<Vec<PolynomialRing>>()
        .into_iter()
        .reduce(|acc, x| acc + x)
        .unwrap()
}

pub fn inner_product_zq_vector(a: &[Zq], b: &[Zq]) -> Zq {
    a.iter().zip(b.iter()).map(|(a, b)| *a * *b).sum()
}

// a: Vec<Vec<PolynomialRing>>, b: Vec<PolynomialRing>
// calculate c = sum(c_i), c_i = poly_vec_times_poly(a_i, b_i)
// c: Vec<PolynomialRing>
pub fn inner_product_poly_matrix_and_poly_vector(
    poly_matrix: &[Vec<PolynomialRing>],
    poly_vector: &[PolynomialRing],
) -> Vec<PolynomialRing> {
    assert_eq!(poly_matrix.len(), poly_vector.len(), "inner_product_poly_matrix_and_poly_vector: poly_matrix and poly_vector must have the same length");
    poly_matrix
        .iter()
        .zip(poly_vector.iter())
        .map(|(poly_matrix_row, poly_vector_element)| {
            poly_vec_times_poly(poly_matrix_row, poly_vector_element)
        })
        .fold(
            vec![zero_poly(); poly_matrix[0].len()],
            |acc, x: Vec<PolynomialRing>| poly_vec_add_poly_vec(&acc, &x),
        )
}

pub fn generate_random_polynomial_ring(deg_bound_d: usize) -> PolynomialRing {
    let mut rng = rand::thread_rng();
    PolynomialRing {
        coefficients: (0..deg_bound_d)
            .map(|_| Zq::from(rng.gen_range(1..5)))
            .collect(),
    }
}

// calculate matrix times vector of PolynomialRing
pub fn matrix_poly_times_poly_vector(
    poly_matrix: &[Vec<PolynomialRing>],
    poly_vec: &Vec<PolynomialRing>,
) -> Vec<PolynomialRing> {
    poly_matrix
        .iter()
        .map(|row| inner_product_polynomial_ring_vector(row, poly_vec))
        .collect::<Vec<PolynomialRing>>()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::vec;

    #[test]
    fn test_poly_vec_times_poly() {
        // Arrange
        let a = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(1), Zq::from(2), Zq::from(3)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(4), Zq::from(5), Zq::from(6)],
            },
        ];
        let b = PolynomialRing {
            coefficients: vec![Zq::from(2), Zq::from(3), Zq::from(4)],
        };

        // Act
        let result = poly_vec_times_poly(&a, &b);
        // expected[0] = a[0] * b = (1 + 2x + 3x^2) * (2 + 3x + 4x^2) = 2 + 7x + 16x^2 + 17x^3 + 12x^4
        // expected[1] = a[1] * b = (4 + 5x + 6x^2) * (2 + 3x + 4x^2) = 8 + 22x + 43x^2 + 38x^3 + 24x^4
        // Assert
        let expected = vec![
            PolynomialRing {
                coefficients: vec![
                    Zq::from(2),
                    Zq::from(7),
                    Zq::from(16),
                    Zq::from(17),
                    Zq::from(12),
                ],
            },
            PolynomialRing {
                coefficients: vec![
                    Zq::from(8),
                    Zq::from(22),
                    Zq::from(43),
                    Zq::from(38),
                    Zq::from(24),
                ],
            },
        ];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_poly_vec_add_poly_vec() {
        // Arrange
        let a = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(1), Zq::from(2), Zq::from(3)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(4), Zq::from(5), Zq::from(6)],
            },
        ];
        let b = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(7), Zq::from(8), Zq::from(9)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(10), Zq::from(11), Zq::from(12)],
            },
        ];

        // Act
        let result = poly_vec_add_poly_vec(&a, &b);

        // Assert
        let expected = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(8), Zq::from(10), Zq::from(12)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(14), Zq::from(16), Zq::from(18)],
            },
        ];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_inner_product_polynomial_ring() {
        let a = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(1), Zq::from(2), Zq::from(3)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(4), Zq::from(5), Zq::from(6)],
            },
        ];
        let b = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(7), Zq::from(8), Zq::from(9)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(10), Zq::from(11), Zq::from(12)],
            },
        ];

        let result = inner_product_polynomial_ring_vector(&a, &b);

        // Expected result calculation:
        // (1 + 2x + 3x^2) * (7 + 8x + 9x^2) = 7 + 22x + 46x^2 + 42x^3 + 27x^4
        // (4 + 5x + 6x^2) * (10 + 11x + 12x^2) = 40 +  96x + 163x^2 + 126x^3 + 72x^4
        // Sum: 47 + 116x + 209x^2 + 168x^3 + 99x^4

        let expected = PolynomialRing {
            coefficients: vec![
                Zq::from(47),
                Zq::from(116),
                Zq::from(209),
                Zq::from(168),
                Zq::from(99),
            ],
        };

        assert_eq!(result.coefficients, expected.coefficients);
    }

    #[test]
    fn test_inner_product_zq_vector() {
        // Arrange
        let a = vec![Zq::from(1), Zq::from(2), Zq::from(3)];
        let b = vec![Zq::from(4), Zq::from(5), Zq::from(6)];

        // Act
        let result = inner_product_zq_vector(&a, &b);

        // Assert
        let expected =
            (Zq::from(1) * Zq::from(4)) + (Zq::from(2) * Zq::from(5)) + (Zq::from(3) * Zq::from(6));
        assert_eq!(
            result, expected,
            "inner_product_zq_vector did not return the correct result"
        );
    }

    #[test]
    fn test_inner_product_poly_matrix_and_poly_vector() {
        // Arrange
        let poly_matrix = vec![
            vec![
                PolynomialRing {
                    coefficients: vec![Zq::from(1), Zq::from(2)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(3), Zq::from(4)],
                },
            ],
            vec![
                PolynomialRing {
                    coefficients: vec![Zq::from(5), Zq::from(6)],
                },
                PolynomialRing {
                    coefficients: vec![Zq::from(7), Zq::from(8)],
                },
            ],
        ];
        let poly_vector = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(9), Zq::from(10)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(11), Zq::from(12)],
            },
        ];

        // Expected Calculation:
        // u = sum D_ij * h_ij^(k) for all k = 1..(t1-1)
        // For this test case:
        // Result[0] = (1+2x) * (9+10x) + (5+6x) * (11+12x) = 64 + 154x + 92x^2
        // Result[1] = (3+4x) * (9+10x) + (7+8x) * (11+12x) = 104 + 238x + 136x^2
        let expected = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(64), Zq::from(154), Zq::from(92)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(104), Zq::from(238), Zq::from(136)],
            },
        ];

        // Act
        let result = inner_product_poly_matrix_and_poly_vector(&poly_matrix, &poly_vector);

        // Assert
        assert_eq!(result, expected, "The inner product of the polynomial matrix and vector did not produce the expected result.");
    }

    #[test]
    fn test_inner_product_polynomial_ring_vector() {
        // Define sample PolynomialRing vectors
        let a = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(1), Zq::from(2), Zq::from(3)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(4), Zq::from(5), Zq::from(6)],
            },
        ];
        let b = vec![
            PolynomialRing {
                coefficients: vec![Zq::from(7), Zq::from(8), Zq::from(9)],
            },
            PolynomialRing {
                coefficients: vec![Zq::from(10), Zq::from(11), Zq::from(12)],
            },
        ];
        // (1 + 2x + 3x^2) * (7 + 8x + 9x^2) + (4 + 5x + 6x^2) * (10 + 11x + 12x^2)
        // = (7 + 22x + 46x^2 + 42x^3 + 27x^4) + (40 + 94x + 163x^2 + 126x^3 + 72x^4)
        // = 47 + 116x + 209x^2 + 168x^3 + 99x^4
        let result = inner_product_polynomial_ring_vector(&a, &b);
        let expected = PolynomialRing {
            coefficients: vec![
                Zq::from(47),  // 47
                Zq::from(116), // 116
                Zq::from(209), // 209
                Zq::from(168), // 168
                Zq::from(99),  // 99
            ],
        };

        assert_eq!(result, expected);
    }
}
