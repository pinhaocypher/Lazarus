use algebra::{polynomial_ring::PolynomialRing, zq::Zq};
// statement
pub struct St {
    pub a_constraint: Vec<Vec<Vec<PolynomialRing>>>,
    pub phi_constraint: Vec<Vec<Vec<PolynomialRing>>>,
    pub b_constraint: Vec<PolynomialRing>,
    pub a_constraint_ct: Vec<Vec<Vec<PolynomialRing>>>,
    pub phi_constraint_ct: Vec<Vec<Vec<PolynomialRing>>>,
    pub b_constraint_ct: Vec<Zq>,
}
