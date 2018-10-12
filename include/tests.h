#ifndef TESTS_H
#define TESTS_H

#include <chrono>
#include "PolynomialFit.h"
#include "argon_cross.h"

void test_polynomial_fit (void);
void test_phase_shift_fit (void);
void test_factor_helper (void);
void test_legendre_polynomial(void);
void test_legendre_intregral (void);
void test_diff_tot_cross (void);
void test_backward_scatter_prob (void);
void test_TM_forward (void);
void test_TM_backward (void);
void test_data_table (void);
void test_resonance_cross (void);
void test_all(void);

#endif
