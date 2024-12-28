use crate::polynomial_ring::PolynomialRing;
use crate::zq::Zq;

#[derive(Debug)]
pub struct RqMatrix {
    pub values: Vec<Vec<PolynomialRing>>,
}

impl RqMatrix {
    pub fn new(kappa: Zq, size_n: Zq) -> Self {
        let size_kappa_usize: usize = kappa.value();
        let size_n_usize: usize = size_n.value();
        let values: Vec<Vec<PolynomialRing>> = (0..size_kappa_usize)
            .map(|_| {
                (0..size_n_usize)
                    .map(|_| PolynomialRing {
                        coefficients: (0..size_n_usize)
                            .map(|_| Zq::from(2)) // we just use 2 as random number to facilitate test
                            .collect(),
                    })
                    .collect()
            })
            .collect();
        assert_eq!(
            values.len(),
            size_kappa_usize,
            "values must have the same length as size_kappa"
        );
        assert_eq!(
            values[0].len(),
            size_n_usize,
            "values[0] must have the same length as size_n"
        );
        RqMatrix { values }
    }
}
