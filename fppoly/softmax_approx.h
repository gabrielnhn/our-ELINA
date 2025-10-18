/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2021 Department of Computer Science, ETH Zurich
 *  This software is distributed under GNU Lesser General Public License Version 3.0.
 *  For more information, see the ELINA project website at:
 *  http://elina.ethz.ch
 *
 *  THE SOFTWARE IS PROVIDED "AS-IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER
 *  EXPRESS, IMPLIED OR STATUTORY, INCLUDING BUT NOT LIMITED TO ANY WARRANTY
 *  THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS OR BE ERROR-FREE AND ANY
 *  IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
 *  TITLE, OR NON-INFRINGEMENT.  IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY     
 *  DAMAGES, INCLUDING BUT NOT LIMITED TO DIRECT, INDIRECT,
 *  SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN
 *  ANY WAY CONNECTED WITH THIS SOFTWARE (WHETHER OR NOT BASED UPON WARRANTY,
 *  CONTRACT, TORT OR OTHERWISE).
 *
 */

#ifndef __SOFTMAX_APPROX_H_INCLUDED__
#define __SOFTMAX_APPROX_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif

#include "backsubstitute.h"

/**
 * Handle softmax layer with LSE-based approximation
 * 
 * Uses linear LSE upper bound to create tight softmax approximations
 * where all input variables appear in the bound coefficients.
 * 
 * @param man ELINA manager
 * @param element Abstract element
 * @param num_neurons Number of neurons (softmax dimension)
 * @param predecessors Predecessor layer indices
 * @param num_predecessors Number of predecessor layers
 * @param use_default_heuristic Use default heuristic (temperature = 1.0)
 */
void handle_softmax_layer(elina_manager_t *man, elina_abstract0_t* element, 
                          size_t num_neurons, size_t *predecessors, 
                          size_t num_predecessors, bool use_default_heuristic);

/**
 * Apply softmax lower bound to expression
 * 
 * @param pr FPPOLY internal structure
 * @param lexpr_p Pointer to lower expression
 * @param neurons Array of input neurons
 * @param num_neurons Number of neurons
 * @param output_idx Which softmax output component
 * @param temperature Temperature parameter
 * @return Upper bound on softmax output
 */
double apply_softmax_lexpr(fppoly_internal_t *pr, expr_t **lexpr_p, 
                           neuron_t **neurons, size_t num_neurons,
                           size_t output_idx, double temperature);

/**
 * Apply softmax upper bound to expression
 * 
 * @param pr FPPOLY internal structure
 * @param uexpr_p Pointer to upper expression
 * @param neurons Array of input neurons
 * @param num_neurons Number of neurons
 * @param output_idx Which softmax output component
 * @param temperature Temperature parameter
 * @return Upper bound on softmax output
 */
double apply_softmax_uexpr(fppoly_internal_t *pr, expr_t **uexpr_p,
                           neuron_t **neurons, size_t num_neurons,
                           size_t output_idx, double temperature);

#ifdef __cplusplus
 }
#endif

#endif