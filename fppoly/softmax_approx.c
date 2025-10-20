#include "softmax_approx.h"
#include <math.h>
#include <glpk.h>
#include <lapacke.h>

// -----------------------------------------------------------
// Utility: LogSumExp computation at a given point
// -----------------------------------------------------------
static double compute_lse_at_point(const double *z, size_t dim, double temperature) {
    double max_val = z[0];
    for (size_t i = 1; i < dim; i++) {
        if (z[i] > max_val) max_val = z[i];
    }
    double sum = 0.0;
    for (size_t i = 0; i < dim; i++) {
        sum += exp((z[i] - max_val) / temperature);
    }
    return max_val + temperature * log(sum);
}

// -----------------------------------------------------------
// Utility: Gradient of LogSumExp (softmax at point z)
// -----------------------------------------------------------
static void compute_lse_gradient(double *grad, const double *z, size_t dim, double temperature) {
    double max_val = z[0];
    for (size_t i = 1; i < dim; i++) {
        if (z[i] > max_val) max_val = z[i];
    }
    double sum = 0.0;
    for (size_t i = 0; i < dim; i++) {
        grad[i] = exp((z[i] - max_val) / temperature);
        sum += grad[i];
    }
    for (size_t i = 0; i < dim; i++) grad[i] /= sum;
}

// -----------------------------------------------------------
// Utility: Generate all 2^dim corners of the input box
// -----------------------------------------------------------
static void generate_corners(double **corners, neuron_t **neurons, size_t dim) {
    size_t num_corners = 1 << dim;
    for (size_t i = 0; i < num_corners; i++) {
        for (size_t j = 0; j < dim; j++) {
            corners[i][j] = (i & (1 << j)) ? neurons[j]->ub : neurons[j]->lb;
        }
    }
}

// -----------------------------------------------------------
// Tangent (lower) affine plane at a given point
// -----------------------------------------------------------
void compute_lse_lower_tangent(double *c_coeffs, double *c_0,
                               const double *point, size_t dim,
                               double temperature) {
    double lse_val = compute_lse_at_point(point, dim, temperature);
    compute_lse_gradient(c_coeffs, point, dim, temperature);
    *c_0 = lse_val;
    for (size_t i = 0; i < dim; i++) *c_0 -= c_coeffs[i] * point[i];
}

// -----------------------------------------------------------
// Main: Compute upper affine bound (DeepPoly convex upper hull)
// -----------------------------------------------------------
static void compute_lse_upper_bound(double *d_coeffs, double *d_0,
                                    neuron_t **neurons, size_t dim,
                                    double temperature) {

    size_t num_corners = 1 << dim;

    // ----------------- Prepare corner samples -----------------
    double **corners = (double **)malloc(num_corners * sizeof(double *));
    for (size_t i = 0; i < num_corners; i++) corners[i] = (double *)malloc(dim * sizeof(double));
    generate_corners(corners, neurons, dim);

    double *center = (double *)malloc(dim * sizeof(double));
    double *half_width = (double *)malloc(dim * sizeof(double));
    for (size_t j = 0; j < dim; j++) {
        center[j] = (neurons[j]->lb + neurons[j]->ub) / 2.0;
        half_width[j] = (neurons[j]->ub - neurons[j]->lb) / 2.0;
        if (half_width[j] < 1e-10) half_width[j] = 1.0;
    }

    double **corners_norm = (double **)malloc(num_corners * sizeof(double *));
    for (size_t i = 0; i < num_corners; i++) {
        corners_norm[i] = (double *)malloc(dim * sizeof(double));
        for (size_t j = 0; j < dim; j++)
            corners_norm[i][j] = (corners[i][j] - center[j]) / half_width[j];
    }

    double *lse_values = (double *)malloc(num_corners * sizeof(double));
    for (size_t i = 0; i < num_corners; i++)
        lse_values[i] = compute_lse_at_point(corners[i], dim, temperature);

    // ----------------- Build LP (GLPK) -----------------
    glp_prob *lp = glp_create_prob();
    glp_set_prob_name(lp, "LSE_UpperBound_LP");
    glp_set_obj_dir(lp, GLP_MAX); // FIXED: we want to maximize intercept b

    glp_add_cols(lp, dim + 1);
    for (size_t j = 0; j < dim + 1; j++) {
        char col_name[16];
        sprintf(col_name, "d%zu", j);
        glp_set_col_name(lp, (int)(j + 1), col_name);
        glp_set_col_bnds(lp, (int)(j + 1), GLP_FR, 0.0, 0.0);
        glp_set_obj_coef(lp, (int)(j + 1), (j == dim) ? 1.0 : 0.0);
    }

    glp_add_rows(lp, (int)num_corners + 1);
    int max_nnz = (int)(num_corners * (dim + 1) + dim + 2);
    int *ia = (int *)malloc((1 + max_nnz) * sizeof(int));
    int *ja = (int *)malloc((1 + max_nnz) * sizeof(int));
    double *ar = (double *)malloc((1 + max_nnz) * sizeof(double));
    int constraint_idx = 1;

    // ----------------- Corner constraints -----------------
    for (size_t i = 0; i < num_corners; i++) {
        char row_name[16];
        sprintf(row_name, "c%zu", i);
        glp_set_row_name(lp, (int)(i + 1), row_name);
        glp_set_row_bnds(lp, (int)(i + 1), GLP_LO, lse_values[i], 0.0);

        for (size_t j = 0; j < dim; j++) {
            ia[constraint_idx] = (int)(i + 1);
            ja[constraint_idx] = (int)(j + 1);
            ar[constraint_idx] = corners_norm[i][j]; // positive
            constraint_idx++;
        }

        ia[constraint_idx] = (int)(i + 1);
        ja[constraint_idx] = (int)(dim + 1);
        ar[constraint_idx] = 1.0; // bias positive
        constraint_idx++;
    }

    // ----------------- Normalization constraint (optional) -----------------
    int row_id = (int)num_corners + 1;
    glp_set_row_name(lp, row_id, "sum_constraint");
    glp_set_row_bnds(lp, row_id, GLP_DB, 0.8, 1.2);

    for (size_t j = 0; j < dim; j++) {
        ia[constraint_idx] = row_id;
        ja[constraint_idx] = (int)(j + 1);
        ar[constraint_idx] = 1.0;
        constraint_idx++;
    }

    glp_load_matrix(lp, constraint_idx - 1, ia, ja, ar);
    free(ia); free(ja); free(ar);

    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_OFF;
    int lp_status = glp_simplex(lp, &parm);

    double lp_solution_norm[dim + 1];
    bool lp_solved_optimally = false;
    if (lp_status == 0 && glp_get_status(lp) == GLP_OPT) {
        lp_solved_optimally = true;
        for (size_t j = 0; j < dim + 1; j++)
            lp_solution_norm[j] = glp_get_col_prim(lp, (int)(j + 1));
    }

    glp_delete_prob(lp);

    // ----------------- Fallback: Tangent if LP failed -----------------
    if (!lp_solved_optimally) {
        double *center_point = (double *)malloc(dim * sizeof(double));
        for (size_t j = 0; j < dim; j++)
            center_point[j] = (neurons[j]->lb + neurons[j]->ub) / 2.0;
        compute_lse_lower_tangent(d_coeffs, d_0, center_point, dim, temperature);
        *d_0 += 1e-6;
        free(center_point);
    } else {
        // Rescale coefficients to original space
        *d_0 = lp_solution_norm[dim];
        for (size_t j = 0; j < dim; j++) {
            d_coeffs[j] = lp_solution_norm[j] / half_width[j];
            *d_0 -= (lp_solution_norm[j] / half_width[j]) * center[j];
        }
    }

    // ----------------- Safety margin to ensure validity -----------------
    double max_violation = 0.0;
    for (size_t i = 0; i < num_corners; i++) {
        double approx = *d_0;
        for (size_t j = 0; j < dim; j++)
            approx += d_coeffs[j] * corners[i][j];
        if (approx < lse_values[i]) {
            double diff = lse_values[i] - approx;
            if (diff > max_violation) max_violation = diff;
        }
    }
    *d_0 += max_violation + 1e-6;

    for (size_t i = 0; i < num_corners; i++) { free(corners[i]); free(corners_norm[i]); }
    free(corners); free(corners_norm); free(center); free(half_width); free(lse_values);
}

// -----------------------------------------------------------
// Bound a linear function over neuron intervals
// -----------------------------------------------------------
static void bound_linear_function(double *phi_min, double *phi_max,
                                 const double *coeffs, double c_0,
                                 neuron_t **neurons, size_t dim) {
    *phi_min = c_0;
    *phi_max = c_0;
    for (size_t j = 0; j < dim; j++) {
        if (coeffs[j] > 0) {
            *phi_min += coeffs[j] * neurons[j]->lb;
            *phi_max += coeffs[j] * neurons[j]->ub;
        } else {
            *phi_min += coeffs[j] * neurons[j]->ub;
            *phi_max += coeffs[j] * neurons[j]->lb;
        }
    }
}

// MAIN FUNCTION TO HANDLE SOFTMAX
static expr_t *create_softmax_expr(fppoly_internal_t *pr, neuron_t *out_neuron,
                                  neuron_t **in_neurons, size_t dim,
                                  size_t output_idx, bool is_lower,
                                  double temperature) {
    
    // ========== CONCRETE PASS DETECTION ==========
    bool is_concrete = true;
    const double concrete_tol = 1e-6;
    for (size_t j = 0; j < dim; j++) {
        fprintf(stderr, "Difference for neuron %zu: %.6e\n", j, in_neurons[j]->ub - in_neurons[j]->lb);
        if (in_neurons[j]->ub - in_neurons[j]->lb > concrete_tol) {
            is_concrete = false;
            break;
        }
    }
    
    if (is_concrete) {
        fprintf(stderr, "DEBUG: Image [ID?] - Concrete pass detected.\n");
        double input_vals[dim];
        for (size_t j = 0; j < dim; j++) {
            input_vals[j] = (in_neurons[j]->lb + in_neurons[j]->ub) / 2.0;
        }
        
        double max_val = input_vals[0];
        for (size_t j = 1; j < dim; j++) {
            if (input_vals[j] > max_val) max_val = input_vals[j];
        }
        
        double sum_exp = 0.0;
        for (size_t j = 0; j < dim; j++) {
            sum_exp += exp((input_vals[j] - max_val) / temperature);
        }
        
        double softmax_val = exp((input_vals[output_idx] - max_val) / temperature) / sum_exp;
        softmax_val = (softmax_val < 0.0) ? 0.0 : ((softmax_val > 1.0) ? 1.0 : softmax_val);
        
        // assign concrete bounds!
        out_neuron->lb = softmax_val;
        out_neuron->ub = softmax_val;
        
        // assign symbolic bounds!
        expr_t *res = alloc_expr();
        res->inf_coeff = (double *)malloc(dim * sizeof(double));
        res->sup_coeff = (double *)malloc(dim * sizeof(double));
        res->inf_cst = softmax_val;
        res->sup_cst = softmax_val;

        // keep it dense, whatever this means
        res->type = DENSE;
        res->size = dim;
        for (size_t j = 0; j < dim; j++) {
            res->inf_coeff[j] = 0.0;
            res->sup_coeff[j] = 0.0;
        }
        res->dim = (size_t *)malloc(dim * sizeof(size_t));
        for (size_t j = 0; j < dim; j++) {
            res->dim[j] = j;
        }
        return res;
    }
    
    // ========== ABSTRACT CASE ==========
    // things are currently going wrong here!

    expr_t *res = alloc_expr();
    res->inf_coeff = (double *)malloc(dim * sizeof(double));
    res->sup_coeff = (double *)malloc(dim * sizeof(double));
    res->type = DENSE;
    res->size = dim;
    res->dim = (size_t *)malloc(dim * sizeof(size_t));
    for (size_t j = 0; j < dim; j++) {
        res->dim[j] = j;
    }
    
    double *d_coeffs = (double *)malloc(dim * sizeof(double));
    double *c_coeffs = (double *)malloc(dim * sizeof(double));
    double *phi_coeffs = (double *)malloc(dim * sizeof(double));
    double d_0, c_0;
    
    if (is_lower) {
        // ========== LOWER BOUND ==========
        // this method is important. look it up/. log-sum-exponente, monotonic increasing function.
        // check jupyter notebook last code cell.
        // log(sum(exp(z))) is convex therefore cool.


        compute_lse_upper_bound(d_coeffs, &d_0, in_neurons, dim, temperature);
        
        for (size_t j = 0; j < dim; j++) {
            if (j == output_idx) {
                phi_coeffs[j] = 1.0 - d_coeffs[j];
            } else {
                phi_coeffs[j] = -d_coeffs[j];
            }
        }
        
        double phi_min, phi_max;
        bound_linear_function(&phi_min, &phi_max, phi_coeffs, -d_0, in_neurons, dim);
        
        double exp_min = exp(phi_min);
        double exp_max = exp(phi_max);
        double slope, intercept;
        
        if (fabs(phi_max - phi_min) < 1e-9) {
            slope = exp_min;
            intercept = 0.0;
        } else {
            slope = (exp_max - exp_min) / (phi_max - phi_min);
            intercept = (phi_max * exp_min - phi_min * exp_max) / (phi_max - phi_min);
        }
        
        for (size_t j = 0; j < dim; j++) {
            res->inf_coeff[j] = slope * phi_coeffs[j];
            res->sup_coeff[j] = -slope * phi_coeffs[j];
        }
        res->inf_cst = slope * (-d_0) + intercept;
        res->sup_cst = -(slope * (-d_0) + intercept);
        
    } else {
        // ========== UPPER BOUND ==========
        double *center = (double *)malloc(dim * sizeof(double));
        for (size_t j = 0; j < dim; j++) {
            center[j] = (in_neurons[j]->lb + in_neurons[j]->ub) / 2.0;
        }
        
        compute_lse_lower_tangent(c_coeffs, &c_0, center, dim, temperature);
        free(center);
        
        for (size_t j = 0; j < dim; j++) {
            if (j == output_idx) {
                phi_coeffs[j] = 1.0 - c_coeffs[j];
            } else {
                phi_coeffs[j] = -c_coeffs[j];
            }
        }
        
        double phi_min, phi_max;
        bound_linear_function(&phi_min, &phi_max, phi_coeffs, -c_0, in_neurons, dim);
        
        double exp_min = exp(phi_min);
        double exp_max = exp(phi_max);
        double slope, intercept;
        
        if (fabs(phi_max - phi_min) < 1e-9) {
            slope = exp_min;
            intercept = exp_min - slope * phi_min;
        } else {
            double chord_slope = (exp_max - exp_min) / (phi_max - phi_min);
            
            if (chord_slope < 1e-9) {
                slope = exp_min;
                intercept = exp_min - slope * phi_min;
            } else {
                slope = chord_slope;
                intercept = slope * (1.0 - log(slope));
            }
        }
        
        double final_intercept = intercept + slope * (-c_0);
        
        // ========== TIER 4 SAFEGUARD ==========
        bool coefficients_safe = true;
        const double TIER4_THRESHOLD = 100.0; // Reduced from 1e4
        double max_coeff = 0.0;
        
        for (size_t j = 0; j < dim; j++) {
            double final_slope = slope * phi_coeffs[j];
            double abs_slope = fabs(final_slope);
            if (abs_slope > max_coeff) max_coeff = abs_slope;
            if (abs_slope > TIER4_THRESHOLD) {
                coefficients_safe = false;
                break;
            }
        }
        
        if (fabs(final_intercept) > TIER4_THRESHOLD) {
            coefficients_safe = false;
        }
        
        if (!coefficients_safe) {
            fprintf(stderr, "TIER4 TRIGGERED: output %zu, max_coeff=%.2e, intercept=%.2e\n",
                    output_idx, max_coeff, final_intercept);
            
            // Better fallback: center-based tangent
            double *center_point = (double *)malloc(dim * sizeof(double));
            for (size_t j = 0; j < dim; j++) {
                center_point[j] = (in_neurons[j]->lb + in_neurons[j]->ub) / 2.0;
            }
            
            double input_vals[dim];
            for (size_t j = 0; j < dim; j++) {
                input_vals[j] = center_point[j];
            }
            double max_val = input_vals[0];
            for (size_t j = 1; j < dim; j++) {
                if (input_vals[j] > max_val) max_val = input_vals[j];
            }
            double sum_exp = 0.0;
            for (size_t j = 0; j < dim; j++) {
                sum_exp += exp((input_vals[j] - max_val) / temperature);
            }
            double softmax_at_center = exp((input_vals[output_idx] - max_val) / temperature) / sum_exp;
            
            for (size_t j = 0; j < dim; j++) {
                res->inf_coeff[j] = 0.0;
                res->sup_coeff[j] = 0.0;
            }
            res->inf_cst = -fmin(1.0, softmax_at_center * 1.5);
            res->sup_cst = fmin(1.0, softmax_at_center * 1.5);
            
            free(center_point);
        } else {
            for (size_t j = 0; j < dim; j++) {
                double final_slope = slope * phi_coeffs[j];
                res->inf_coeff[j] = -final_slope;
                res->sup_coeff[j] = final_slope;
            }
            res->inf_cst = -final_intercept;
            res->sup_cst = final_intercept;
        }
    }
    
    // ========== FINAL BOUND COMPUTATION ==========
    double lb_val = res->inf_cst;
    double ub_val = res->sup_cst;
    
    for (size_t j = 0; j < dim; j++) {
        if (res->inf_coeff[j] > 0) {
            lb_val += res->inf_coeff[j] * in_neurons[j]->lb;
        } else {
            lb_val += res->inf_coeff[j] * in_neurons[j]->ub;
        }
        
        if (res->sup_coeff[j] > 0) {
            ub_val += res->sup_coeff[j] * in_neurons[j]->ub;
        } else {
            ub_val += res->sup_coeff[j] * in_neurons[j]->lb;
        }
    }
    
    if (is_lower) {
        out_neuron->lb = (lb_val < 0.0) ? 0.0 : lb_val;
    } else {
        out_neuron->ub = (ub_val > 1.0) ? 1.0 : fmax(1e-15, ub_val);
    }
    
    free(d_coeffs);
    free(c_coeffs);
    free(phi_coeffs);
    
    return res;
}

void handle_softmax_layer(elina_manager_t *man, elina_abstract0_t *element,
                         size_t num_neurons, size_t *predecessors,
                         size_t num_predecessors, bool use_default_heuristic) {
    assert(num_predecessors == 1);
    fppoly_t *fp = fppoly_of_abstract0(element);
    fppoly_internal_t *pr = fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
    size_t numlayers = fp->numlayers;
    
    fppoly_add_new_layer(fp, num_neurons, predecessors, num_predecessors, true);
    neuron_t **out_neurons = fp->layers[numlayers]->neurons;
    int k = predecessors[0] - 1;
    neuron_t **in_neurons = fp->layers[k]->neurons;
    
    // **FIXED: Compute temperature from ACTUAL interval widths**
    double max_interval_width = 0.0;
    for (size_t j = 0; j < num_neurons; j++) {
        double width = in_neurons[j]->ub - in_neurons[j]->lb;
        if (width > max_interval_width) {
            max_interval_width = width;
        }
    }
    
    // Temperature scales with interval width to stabilize approximation
    // For small intervals: T ≈ 1, for large: T increases
    double temperature = 1.0 + max_interval_width / 10.0;
    
    fprintf(stderr, "DEBUG: Adaptive temperature: %.2f (max_interval_width: %.2f)\n", 
            temperature, max_interval_width);
    
    for (size_t i = 0; i < num_neurons; i++) {
        out_neurons[i]->lexpr = create_softmax_expr(pr, out_neurons[i], in_neurons,
                                                    num_neurons, i, true, temperature);
        out_neurons[i]->uexpr = create_softmax_expr(pr, out_neurons[i], in_neurons,
                                                    num_neurons, i, false, temperature);
    }
}

static double apply_softmax_expr(fppoly_internal_t *pr, expr_t **expr_p,
                                neuron_t **neurons, size_t num_neurons,
                                size_t output_idx, bool is_lower,
                                double temperature) {
    expr_t *expr = *expr_p;
    
    neuron_t *tmp_neuron = neuron_alloc();
    expr_t *res = create_softmax_expr(pr, tmp_neuron, neurons, num_neurons,
                                     output_idx, is_lower, temperature);
    
    double *slope_inf = res->inf_coeff;
    double *slope_sup = res->sup_coeff;
    double intercept_inf = res->inf_cst;
    double intercept_sup = res->sup_cst;
    
    size_t expr_size = expr->size;
    
    for (size_t i = 0; i < expr_size; i++) {
        elina_double_interval_mul_expr_coeff(pr, &expr->inf_coeff[i], &expr->sup_coeff[i],
                                            slope_inf[i], slope_sup[i],
                                            expr->inf_coeff[i], expr->sup_coeff[i]);
    }
    
    elina_double_interval_mul_cst_coeff(pr, &expr->inf_cst, &expr->sup_cst,
                                       intercept_inf, intercept_sup,
                                       expr->inf_cst, expr->sup_cst);
    
    elina_double_interval_add_cst_coeff(pr, &expr->inf_cst, &expr->sup_cst,
                                       intercept_inf, intercept_sup,
                                       expr->inf_cst, expr->sup_cst);
    
    double result_ub = tmp_neuron->ub;
    
    free_expr(res);
    free_neuron(tmp_neuron);
    
    return result_ub;
}

double apply_softmax_lexpr(fppoly_internal_t *pr, expr_t **lexpr_p,
                          neuron_t **neurons, size_t num_neurons,
                          size_t output_idx, double temperature) {
    return apply_softmax_expr(pr, lexpr_p, neurons, num_neurons,
                             output_idx, true, temperature);
}

double apply_softmax_uexpr(fppoly_internal_t *pr, expr_t **uexpr_p,
                          neuron_t **neurons, size_t num_neurons,
                          size_t output_idx, double temperature) {
    return apply_softmax_expr(pr, uexpr_p, neurons, num_neurons,
                             output_idx, false, temperature);
}