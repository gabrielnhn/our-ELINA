// #include "expr.h"

// void elina_double_interval_add_expr_coeff(fppoly_internal_t *pr, double * res_inf, double *res_sup, double inf, double sup, double inf_expr, double sup_expr){
// 	*res_inf = inf + inf_expr;
// 	*res_sup = sup + sup_expr;
// 	double maxA = fmax(fabs(inf_expr),fabs(sup_expr));
// 	double tmp1, tmp2;
// 	elina_double_interval_mul(&tmp1,&tmp2, inf, sup, maxA*pr->ulp, maxA*pr->ulp);
// 	*res_inf += tmp1;
// 	*res_sup += tmp2;
// }


// void elina_double_interval_add_cst_coeff(fppoly_internal_t *pr, double * res_inf, double *res_sup, double inf, double sup, double inf_expr, double sup_expr){
// 	elina_double_interval_add_expr_coeff(pr, res_inf, res_sup, inf, sup, inf_expr, sup_expr);
// 	*res_inf += pr->min_denormal;
// 	*res_sup += pr->min_denormal;	
// }


// void elina_double_interval_mul_expr_coeff(fppoly_internal_t *pr, double * res_inf, double *res_sup, double inf, double sup, double inf_expr, double sup_expr){
// 	elina_double_interval_mul(res_inf,res_sup,inf,sup,inf_expr,sup_expr);
// 	double maxA = fmax(fabs(inf_expr),fabs(sup_expr));
// 	double tmp1, tmp2;
// 	elina_double_interval_mul(&tmp1,&tmp2, inf, sup, maxA*pr->ulp, maxA*pr->ulp);
// 	*res_inf += tmp1;
// 	*res_sup += tmp2;
// }

// void elina_double_interval_mul_cst_coeff(fppoly_internal_t *pr, double * res_inf, double *res_sup, double inf, double sup, double inf_expr, double sup_expr){
// 	elina_double_interval_mul_expr_coeff(pr, res_inf, res_sup, inf, sup, inf_expr, sup_expr);
// 	*res_inf += pr->min_denormal;
// 	*res_sup += pr->min_denormal;	
// }


// void expr_fprint(FILE * stream, expr_t *expr){
// 	if((expr->inf_coeff==NULL) || (expr->sup_coeff==NULL)){
// 		fprintf(stdout,"+ [%g, %g]\n",-expr->inf_cst,expr->sup_cst);
// 		return;
// 	}
// 	size_t size = expr->size;
// 	size_t i;
// 	for(i=0; i < size; i++){
// 		if(i==0){
// 			if(expr->type==DENSE){
// 				fprintf(stream, "[%g, %g]x0 ", -expr->inf_coeff[0],expr->sup_coeff[0]);
// 			}
// 			else{
// 				fprintf(stream, "[%g, %g]x%zu ", -expr->inf_coeff[0],expr->sup_coeff[0],expr->dim[0]);
// 			}
// 		}
		
// 		else{
// 			if(expr->type==DENSE){
// 				fprintf(stream,"+ [%g, %g]x%zu ",-expr->inf_coeff[i],expr->sup_coeff[i],i);
// 			}
// 			else{
// 				fprintf(stream,"+ [%g, %g]x%zu ",-expr->inf_coeff[i],expr->sup_coeff[i],expr->dim[i]);
// 			}
// 		}
// 	}
	
// 	fprintf(stdout,"+ [%g, %g]\n",-expr->inf_cst,expr->sup_cst);
	
// }


// void expr_print(expr_t * expr){
// 	expr_fprint(stdout, expr);	
// }

// expr_t * alloc_expr(void){
// 	expr_t *expr = (expr_t *)malloc(sizeof(expr_t));
// 	expr->inf_coeff = NULL;
// 	expr->sup_coeff = NULL;
// 	expr->dim = NULL;
// 	return expr;
// }

// expr_t * create_dense_expr(double *coeff, double cst, size_t size){
// 	expr_t *expr = (expr_t *)malloc(sizeof(expr_t));
// 	expr->inf_coeff = (double *)malloc(size*sizeof(double));
// 	expr->sup_coeff = (double *)malloc(size*sizeof(double));
// 	expr->dim= NULL;
// 	size_t i;
// 	expr->size = size;
// 	expr->inf_cst = -cst;
// 	expr->sup_cst = cst;
// 	expr->type = DENSE;
// 	for(i=0; i < size; i++){
// 		expr->inf_coeff[i] = -coeff[i];
// 		expr->sup_coeff[i] = coeff[i];
// 	}
// 	return expr;
// }


// expr_t * create_cst_expr(double l, double u){
// 	expr_t *expr = (expr_t*)malloc(sizeof(expr_t));
// 	expr->inf_coeff = NULL;
// 	expr->sup_coeff = NULL;
// 	expr->dim = NULL;
// 	expr->type = SPARSE;
// 	expr->size = 0;
// 	expr->inf_cst = l;
// 	expr->sup_cst = u;
// 	return expr;
// }

// expr_t * create_sparse_expr(double *coeff, double cst, size_t *dim, size_t size){
// 	expr_t *expr = (expr_t *)malloc(sizeof(expr_t));
// 	if(size>0){
// 		expr->inf_coeff = (double *)malloc(size*sizeof(double));
// 		expr->sup_coeff = (double *)malloc(size*sizeof(double));
// 		expr->dim = (size_t *)malloc(size*sizeof(size_t));
// 	}
// 	else{
// 		expr->inf_coeff = NULL;
// 		expr->sup_coeff = NULL;
// 		expr->dim = NULL;
// 	}
// 	size_t i;
// 	expr->size = size;
// 	expr->inf_cst = -cst;
// 	expr->sup_cst = cst;
// 	expr->type = SPARSE;
// 	for(i=0; i < size; i++){
// 		expr->inf_coeff[i] = -coeff[i];
// 		expr->sup_coeff[i] = coeff[i];
// 		expr->dim[i] = dim[i];
// 	}
// 	return expr;
// }


// void free_expr(expr_t *expr){
// 	if(expr->inf_coeff){
// 		free(expr->inf_coeff);
// 		expr->inf_coeff = NULL;
// 	}
// 	if(expr->sup_coeff){
// 		free(expr->sup_coeff);
// 		expr->sup_coeff = NULL;
// 	}
// 	if(expr->type==SPARSE && expr->dim){
// 		free(expr->dim);
// 	}
// 	expr->dim = NULL;
// 	free(expr);
// 	expr = NULL;  
// }

// expr_t * copy_cst_expr(expr_t *src){
// 	expr_t *dst = (expr_t *)malloc(sizeof(expr_t));
// 	dst->inf_coeff = NULL;
// 	dst->sup_coeff = NULL;
// 	dst->inf_cst = src->inf_cst;
// 	dst->sup_cst = src->sup_cst; 
// 	dst->type = src->type;
// 	dst->dim = NULL;
// 	dst->size = src->size; 
// 	return dst;
// }



// expr_t * copy_expr(expr_t *src){
// 	expr_t *dst = (expr_t *)malloc(sizeof(expr_t));
// 	dst->inf_coeff = (double *)malloc(src->size*sizeof(double));
// 	dst->sup_coeff = (double *)malloc(src->size*sizeof(double));
	
// 	size_t i;
// 	dst->inf_cst = src->inf_cst;
// 	dst->sup_cst = src->sup_cst; 
// 	dst->type = src->type;
// 	for(i=0; i < src->size; i++){
// 		dst->inf_coeff[i] = src->inf_coeff[i];
// 		dst->sup_coeff[i] = src->sup_coeff[i];
// 	}
// 	if(src->type==SPARSE){
// 		dst->dim = (size_t *)malloc(src->size*sizeof(size_t));
// 		for(i=0; i < src->size; i++){
// 			dst->dim[i] = src->dim[i];
// 		}
// 	}
// 	dst->size = src->size; 
// 	return dst;
// }

// expr_t* concretize_dense_sub_expr(fppoly_internal_t *pr, expr_t * expr, double *inf, double *sup, size_t start, size_t size){
// 	expr_t * res = (expr_t *)malloc(sizeof(expr_t));
// 	res->inf_coeff = (double *)malloc(start*sizeof(double));
// 	res->sup_coeff = (double *)malloc(start*sizeof(double));
// 	size_t i;
// 	res->inf_cst = expr->inf_cst;
// 	res->sup_cst = expr->sup_cst;
// 	res->type = expr->type;
// 	for(i=0; i < start; i++){
// 		res->inf_coeff[i] = expr->inf_coeff[i];
// 		res->sup_coeff[i] = expr->sup_coeff[i];
// 	}
// 	for(i=start; i< size;i++){
// 		double tmp1,tmp2;
// 		elina_double_interval_mul_expr_coeff(pr,&tmp1,&tmp2,inf[i-start],sup[i-start],expr->inf_coeff[i],expr->sup_coeff[i]);
// 		res->inf_cst += tmp1;
// 		res->sup_cst += tmp2;
// 	}
// 	res->size = start;
// 	return res;
// }


// void merge_sparse_expr(expr_t *expr, size_t l, size_t m, size_t r) {
//     int i, j, k;
//     int n1 = m - l + 1;
//     int n2 = r - m;

//     /* create temp arrays */
//     size_t *L = (size_t *)malloc(n1*sizeof(size_t));
//     size_t *R = (size_t *)malloc(n2*sizeof(size_t));
//     double *L2 = (double *)malloc(n1*sizeof(double));
//     double *R2 = (double *)malloc(n2*sizeof(double));
//     double *L3 = (double *)malloc(n1*sizeof(double));
//     double *R3 = (double *)malloc(n2*sizeof(double));
    
//     /* Copy data to temp arrays L[] and R[] */
//     for (i = 0; i < n1; i++) {
//         L[i] = expr->dim[l + i];
//         L2[i] = expr->inf_coeff[l + i];
// 	L3[i] = expr->sup_coeff[l + i];
//     }
//     for (j = 0; j < n2; j++) {
//         R[j] = expr->dim[m + 1 + j];
//         R2[j] = expr->inf_coeff[m + 1 + j];
// 	R3[j] = expr->sup_coeff[m + 1 + j];
//     }

//     /* Merge the temp arrays back into arr[l..r]*/
//     i = 0; // Initial index of first subarray
//     j = 0; // Initial index of second subarray
//     k = l; // Initial index of merged subarray
//     while (i < n1 && j < n2) {
//         if (L[i] <= R[j]) {
//             expr->dim[k] = L[i];
//             expr->inf_coeff[k] = L2[i];
// 	    expr->sup_coeff[k] = L3[i];
//             i++;
//         } else {
//             expr->dim[k] = R[j];
//             expr->inf_coeff[k] = R2[j];
// 	    expr->sup_coeff[k] = R3[j];
//             j++;
//         }
//         k++;
//     }

//     /* Copy the remaining elements of L[], if there
//        are any */
//     while (i < n1) {
//         expr->dim[k] = L[i];
//         expr->inf_coeff[k] = L2[i];
// 	expr->sup_coeff[k] = L3[i];
//         i++;
//         k++;
//     }

//     /* Copy the remaining elements of R[], if there
//        are any */
//     while (j < n2) {
//         expr->dim[k] = R[j];
//         expr->inf_coeff[k] = R2[j];
// 	expr->sup_coeff[k] = R3[j];
//         j++;
//         k++;
//     }
//     free(L);
//     free(R);
//     free(L2);
//     free(R2);
//     free(L3);
//     free(R3);
// }


// /* l is for left index and r is right index of the
//    sub-array of arr to be sorted */
// void merge_sort_sparse_expr(expr_t *expr, size_t l, size_t r) {
//     if (l < r) {
//         // Same as (l+r)/2, but avoids overflow for
//         // large l and h
//         size_t m = l + (r - l) / 2;

//         // Sort first and second halves
//         merge_sort_sparse_expr(expr, l, m);
//         merge_sort_sparse_expr(expr, m + 1, r);

//         merge_sparse_expr(expr, l, m, r);
//     }
// }

// void sort_sparse_expr(expr_t *expr){
// 	merge_sort_sparse_expr(expr,0,expr->size-1);
// }


// expr_t * multiply_expr(fppoly_internal_t *pr, expr_t *expr, double mul_inf, double mul_sup){
// 	expr_t * res = alloc_expr();
// 	if(expr->size > 0){
// 		res->inf_coeff = malloc(expr->size*sizeof(double));
// 		res->sup_coeff = malloc(expr->size*sizeof(double));
// 	}
// 	else{
// 		res->inf_coeff = NULL;		
// 		res->sup_coeff = NULL;
// 	}
// 	res->type = expr->type;
// 	size_t i;
// 	for(i=0; i < expr->size; i++){
// 		//res->coeff[i] = mul_coeff*expr->coeff[i];
// 		elina_double_interval_mul_expr_coeff(pr,&res->inf_coeff[i],&res->sup_coeff[i],mul_inf,mul_sup,expr->inf_coeff[i],expr->sup_coeff[i]);
		
// 	}
// 	if(expr->type==SPARSE){
// 		if(expr->size>0){
// 			res->dim = (size_t*)malloc(expr->size*sizeof(size_t));
// 			for(i=0; i < expr->size; i++){
// 				res->dim[i] = expr->dim[i];
// 			}
// 		}
// 		else{
// 			res->dim = NULL;
// 		}
// 	}
// 	res->size = expr->size;
	
// 	elina_double_interval_mul_cst_coeff(pr,&res->inf_cst,&res->sup_cst,mul_inf,mul_sup,expr->inf_cst,expr->sup_cst);
	
// 	//res->cst = mul_coeff*expr->cst;
// 	return res;
// }


// expr_t * multiply_cst_expr(fppoly_internal_t *pr, expr_t *expr, double mul_inf, double mul_sup){
// 	expr_t * res = alloc_expr();
// 	res->inf_coeff = NULL;		
// 	res->sup_coeff = NULL;
// 	res->dim = NULL;
// 	res->type = expr->type;
// 	res->size = expr->size;
// 	elina_double_interval_mul_cst_coeff(pr,&res->inf_cst,&res->sup_cst,mul_inf,mul_sup,expr->inf_cst,expr->sup_cst);
// 	//res->cst = mul_coeff*expr->cst;
// 	return res;
// }






// void add_cst_expr(fppoly_internal_t *pr, expr_t * exprA, expr_t *exprB){
// 	double maxA = fmax(fabs(exprA->inf_cst),fabs(exprA->sup_cst));
// 	double maxB = fmax(fabs(exprB->inf_cst),fabs(exprB->sup_cst));
// 	exprA->inf_cst = exprA->inf_cst + exprB->inf_cst  + (maxA + maxB)*pr->ulp + pr->min_denormal; 
// 	exprA->sup_cst = exprA->sup_cst + exprB->sup_cst + (maxA + maxB)*pr->ulp + pr->min_denormal; 
// 	return;
// }

// //A = A + B
// void add_expr(fppoly_internal_t *pr,expr_t * exprA, expr_t * exprB){
// 	//
// 	size_t sizeB = exprB->size;
// 	if(sizeB==0){
// 		double maxA = fmax(fabs(exprA->inf_cst),fabs(exprA->sup_cst));
// 		double maxB = fmax(fabs(exprB->inf_cst),fabs(exprB->sup_cst));
// 		exprA->inf_cst = exprA->inf_cst + exprB->inf_cst  + (maxA + maxB)*pr->ulp + pr->min_denormal; 
// 		exprA->sup_cst = exprA->sup_cst + exprB->sup_cst  + (maxA + maxB)*pr->ulp + pr->min_denormal; 
// 		return;
// 	}
// 	size_t i;
// 	if(exprA->size==0){
		
// 		exprA->size = exprB->size;
// 		double maxA = fmax(fabs(exprA->inf_cst),fabs(exprA->sup_cst));
// 		double maxB = fmax(fabs(exprB->inf_cst),fabs(exprB->sup_cst));
// 		exprA->inf_cst += exprB->inf_cst  + (maxA + maxB)*pr->ulp + pr->min_denormal;
// 		exprA->sup_cst += exprB->sup_cst  + (maxA + maxB)*pr->ulp + pr->min_denormal;
// 		exprA->inf_coeff = (double *)malloc(sizeB*sizeof(double));
// 		exprA->sup_coeff = (double *)malloc(sizeB*sizeof(double));
// 		for(i=0; i < sizeB; i++){
// 			exprA->inf_coeff[i] = exprB->inf_coeff[i];
// 			exprA->sup_coeff[i] = exprB->sup_coeff[i];
// 		} 
// 		exprA->type = exprB->type;
// 		if(exprA->type==SPARSE){
// 			exprA->dim = (size_t *)malloc(sizeB*sizeof(size_t));
// 			for(i=0; i < sizeB; i++){
// 				exprA->dim[i] = exprB->dim[i];
// 			}
// 		} 
		
// 		return;
// 	}
// 	else{
// 		size_t sizeA = exprA->size;
// 		assert(sizeA==sizeB);
// 		double maxA = fmax(fabs(exprA->inf_cst),fabs(exprA->sup_cst));
// 		double maxB = fmax(fabs(exprB->inf_cst),fabs(exprB->sup_cst));
// 		exprA->inf_cst += exprB->inf_cst  + (maxA + maxB)*pr->ulp + pr->min_denormal;
// 		exprA->sup_cst += exprB->sup_cst  + (maxA + maxB)*pr->ulp + pr->min_denormal;
// 		if(exprA->type==DENSE){
// 			if(exprB->type==DENSE){
// 				//printf("AA\n");
// 				//fflush(stdout);
// 				for(i=0; i < sizeB; i++){
// 					maxA = fmax(fabs(exprA->inf_coeff[i]),fabs(exprA->sup_coeff[i]));
// 					maxB = fmax(fabs(exprB->inf_coeff[i]),fabs(exprB->sup_coeff[i]));
// 					exprA->inf_coeff[i] = exprA->inf_coeff[i] + exprB->inf_coeff[i] + (maxA + maxB)*pr->ulp;
// 					exprA->sup_coeff[i] = exprA->sup_coeff[i] + exprB->sup_coeff[i] + (maxA + maxB)*pr->ulp;
// 				}
// 			}
// 			else{
// 				//printf("AB\n");	
// 				//fflush(stdout);
// 				size_t k = 0;
// 				for(i=0; i < sizeA; i++){
// 					if(k < sizeB && exprB->dim[k]==i){
// 						maxA = fmax(fabs(exprA->inf_coeff[i]),fabs(exprA->sup_coeff[i]));
// 						maxB = fmax(fabs(exprB->inf_coeff[k]),fabs(exprB->sup_coeff[k]));
// 						exprA->inf_coeff[i] = exprA->inf_coeff[i] + exprB->inf_coeff[k] + (maxA + maxB)*pr->ulp ;
// 						exprA->sup_coeff[i] = exprA->sup_coeff[i] + exprB->sup_coeff[k] + (maxA + maxB)*pr->ulp;
// 						k++;
// 					}
// 				}
// 			}
// 		}
// 		else{
// 			size_t sizeB = exprB->size;
// 			size_t k;
// 			double * new_inf_coeff;
// 			double * new_sup_coeff;
// 			if(exprB->type==DENSE){
// 				//printf("BB\n");
// 				//fflush(stdout);
// 				i=0;
// 				new_inf_coeff = (double *)malloc(sizeB*sizeof(double));
// 				new_sup_coeff = (double *)malloc(sizeB*sizeof(double));				
// 				for(k=0; k < sizeB; k++){
// 					if(i < sizeA && exprA->dim[i] == k){
// 						maxA = fmax(fabs(exprA->inf_coeff[i]),fabs(exprA->sup_coeff[i]));
// 						maxB = fmax(fabs(exprB->inf_coeff[k]),fabs(exprB->sup_coeff[k]));
// 						new_inf_coeff[k] = exprA->inf_coeff[i] + exprB->inf_coeff[k] + (maxA + maxB)*pr->ulp;
// 						new_sup_coeff[k] = exprA->sup_coeff[i] + exprB->sup_coeff[k] + (maxA + maxB)*pr->ulp;
// 						i++;
// 					}
// 					else{
// 						new_inf_coeff[k] = exprB->inf_coeff[k];
// 						new_sup_coeff[k] = exprB->sup_coeff[k];
// 					}
// 				}
// 				exprA->type = DENSE;
// 				exprA->size = sizeB;
// 				free(exprA->dim);
// 				exprA->dim = NULL;
// 			}
// 			else{
// 				//printf("BC %zu %zu\n",exprA->size, exprB->size);
// 				//fflush(stdout);
// 				i=0;
// 				k=0;
// 				size_t l = 0;
// 				//printf("before sort\n");
// 				//expr_print(exprA);
// 				//fflush(stdout);
// 				//if(exprA->size>0){
// 				//	sort_sparse_expr(exprA);
// 				//}
// 				//printf("after sort\n");
// 				//expr_print(exprA);
// 				//fflush(stdout);
				
// 				new_inf_coeff = (double *)malloc((sizeA+sizeB)*sizeof(double));
// 				new_sup_coeff = (double *)malloc((sizeA+sizeB)*sizeof(double));
// 				size_t * new_dim = (size_t *)malloc((sizeA+sizeB)*sizeof(size_t));
// 				while(i < sizeA && k < sizeB){
// 					if(exprA->dim[i] < exprB->dim[k]){
// 						new_inf_coeff[l] = exprA->inf_coeff[i];
// 						new_sup_coeff[l] = exprA->sup_coeff[i];
// 						new_dim[l] = exprA->dim[i];
// 						i++;
						
// 					}
// 					else if(exprB->dim[k] < exprA->dim[i]){
// 						new_inf_coeff[l] = exprB->inf_coeff[k];
// 						new_sup_coeff[l] = exprB->sup_coeff[k];
// 						new_dim[l] = exprB->dim[k];
// 						k++;
// 					}
// 					else{
// 						maxA = fmax(fabs(exprA->inf_coeff[i]),fabs(exprA->sup_coeff[i]));
// 						maxB = fmax(fabs(exprB->inf_coeff[k]),fabs(exprB->sup_coeff[k]));
// 						new_inf_coeff[l] = exprA->inf_coeff[i] + exprB->inf_coeff[k] + (maxA + maxB)*pr->ulp;
// 						new_sup_coeff[l] = exprA->sup_coeff[i] + exprB->sup_coeff[k] + (maxA + maxB)*pr->ulp;
// 						new_dim[l] = exprA->dim[i];
// 						i++;
// 						k++;
// 					}
// 					l++;
// 				}
// 				while(i < sizeA){
// 					new_inf_coeff[l] = exprA->inf_coeff[i];
// 					new_sup_coeff[l] = exprA->sup_coeff[i];
// 					new_dim[l] = exprA->dim[i];
// 					i++;
// 					l++;
// 				}
// 				while(k < sizeB){
// 					new_inf_coeff[l] = exprB->inf_coeff[k];
// 					new_sup_coeff[l] = exprB->sup_coeff[k];
// 					new_dim[l] = exprB->dim[k];
// 					k++;
// 					l++;
// 				}
				
// 				new_inf_coeff = (double*)realloc(new_inf_coeff,l*sizeof(double));
// 				new_sup_coeff = (double*)realloc(new_sup_coeff,l*sizeof(double));
// 				free(exprA->dim);
// 				exprA->dim = NULL;
// 				new_dim = (size_t *)realloc(new_dim,l*sizeof(size_t));
// 				exprA->dim = new_dim;
// 				exprA->size = l;
// 			}
// 			if(exprA->inf_coeff){
// 				free(exprA->inf_coeff);
// 				exprA->inf_coeff = NULL;
// 			}
// 			if(exprA->sup_coeff){
// 				free(exprA->sup_coeff);
// 				exprA->sup_coeff = NULL;
// 			}
// 			exprA->inf_coeff = new_inf_coeff;
// 			exprA->sup_coeff = new_sup_coeff;
			
// 		}
// 	}
// }

// expr_t * extract_subexpr_concatenate(expr_t * expr, size_t index, size_t* C, size_t num_neurons, size_t num_channels){
// 	size_t i, j=0, k;
// 	//size_t index_end = index_start+num_neurons;
// 	expr_t * res = alloc_expr();
// 	size_t res_size = 0;
// 	res->inf_cst = 0.0;
// 	res->sup_cst = 0.0;
	
// 	size_t hw = num_neurons/num_channels;
// 	size_t offset = 0;
// 	for(i=0; i < index; i++){
// 		offset = offset + C[i];
// 	}
// 	//printf("offset: %zu %zu %zu\n", offset, C[0],C[1]);
	
// 	if(expr->type==DENSE){
// 		res->type = DENSE;
// 		size_t num_neurons_in_layer = C[index] * hw;
// 		res->inf_coeff = (double *)malloc(num_neurons_in_layer*sizeof(double));
// 		res->sup_coeff = (double *)malloc(num_neurons_in_layer*sizeof(double));
// 		res->size = num_neurons_in_layer;
// 		for(i=0; i < hw; i++){
// 			//printf("START HERE: %zu %zu\n", i*num_channels+offset, num_neurons_in_layer);
// 			//fflush(stdout);
// 			for(k=0; k < C[index]; k++){
// 				res->inf_coeff[j] = expr->inf_coeff[i*num_channels + offset + k];
// 				res->sup_coeff[j] = expr->sup_coeff[i*num_channels + offset + k];
// 				j++;
// 			}
// 		}
// 	}
// 	else{
// 		size_t *start_indices = (size_t *)malloc(hw*sizeof(size_t));
// 		size_t *end_indices = (size_t *)malloc(hw*sizeof(size_t));
// 		for(i=0; i < hw; i++){
// 			start_indices[i] = i*num_channels + offset;
// 			end_indices[i] = i*num_channels + offset + C[index];
// 		}
// 		size_t res_size = 0;
// 		for(i = 0; i < expr->size; i++){
// 			k = expr->dim[i];
// 			size_t l;
// 			for(l=0; l < hw; l++){
// 				if(k>=start_indices[l] && k < end_indices[l]){
// 					res_size++;
// 					break;
// 				}
// 			}
			
// 		}
// 		res->inf_coeff = (double *)malloc(res_size*sizeof(double));
// 		res->sup_coeff = (double *)malloc(res_size*sizeof(double));
// 		res->dim = (size_t *)malloc(res_size*sizeof(size_t));
// 		res->size = res_size;
// 		res->type = SPARSE;
// 		for(i=0; i < expr->size; i++){
// 			k = expr->dim[i];
// 			size_t l;
// 			for(l=0; l < hw; l++){
// 				if(k>=start_indices[l] && k < end_indices[l]){
// 					res->inf_coeff[j] = expr->inf_coeff[i];
// 					res->sup_coeff[j] = expr->sup_coeff[i];
// 					res->dim[j] = l*C[index] + k-start_indices[l];
// 					j++;
// 					break;
// 				}
// 			}
// 		}
// 		free(start_indices);
// 		free(end_indices);
// 	}
// 	return res;
	
// }

// expr_t * expr_replace_bounds_affine(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons, bool is_lower){
// 	if(expr->size==0){
// 		return copy_cst_expr(expr);
// 	}	
// 	if(expr->inf_coeff==NULL || expr->sup_coeff==NULL ){
// 		return alloc_expr();
// 	}
	
// 	size_t num_neurons = expr->size;
// 	size_t i,k;
// 	expr_t * res;
// 	if(expr->type==DENSE){
// 		k = 0;
// 	}
// 	else{
// 		k = expr->dim[0];		
// 	}
		
// 	expr_t * mul_expr = NULL;
// 	neuron_t * neuron_k = neurons[k];
// 	if(is_lower){
// 		if(expr->sup_coeff[0] < 0){
// 			mul_expr = neuron_k->uexpr;
// 		}
// 		else if(expr->inf_coeff[0]<0){
// 			mul_expr = neuron_k->lexpr;
// 		}
// 	}
// 	else{
// 		if(expr->sup_coeff[0] < 0){
// 			mul_expr = neuron_k->lexpr;
// 		}
// 		else if(expr->inf_coeff[0]<0){
// 			mul_expr = neuron_k->uexpr;
// 		}
// 	}	
	
// 	if(mul_expr==NULL){
		
// 		double tmp1=0.0, tmp2=0.0;
// 		if(expr->inf_coeff[0]!=0 || expr->sup_coeff[0]!=0){
// 			elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,neuron_k->lb,neuron_k->ub,expr->inf_coeff[0],expr->sup_coeff[0]);
// 		}
// 		double coeff[1];
// 		size_t dim[1];
// 		coeff[0] = 0;
// 		dim[0] = 0;
// 		if(is_lower){
// 			res = create_sparse_expr(coeff,-tmp1,dim,1);
// 		}
// 		else{
// 			res = create_sparse_expr(coeff,tmp2,dim,1);
// 		}
// 	}
// 	else if(mul_expr->size==0){
			
// 		res = multiply_cst_expr(pr,mul_expr,expr->inf_coeff[0],expr->sup_coeff[0]);
			
// 	}
// 	else{
			
// 		res = multiply_expr(pr,mul_expr,expr->inf_coeff[0],expr->sup_coeff[0]);
// 	}
    	
// 	for(i=1; i < num_neurons; i++){
// 		if(expr->type==DENSE){
// 			k = i;
// 		}
// 		else{
// 			k = expr->dim[i];
// 		}
// 		neuron_k = neurons[k];
// 		mul_expr = NULL;
// 		if(is_lower){
// 			if(expr->sup_coeff[i] < 0){
// 				mul_expr = neuron_k->uexpr;
// 			}
// 			else if(expr->inf_coeff[i]<0){
// 				mul_expr = neuron_k->lexpr;
// 			}
// 		}
// 		else{
// 			if(expr->sup_coeff[i] < 0){
// 				mul_expr = neuron_k->lexpr;
// 			}
// 			else if(expr->inf_coeff[i]<0){
// 				mul_expr = neuron_k->uexpr;
// 			}
// 		}
// 		if(expr->sup_coeff[i]==0 && expr->inf_coeff[i]==0){
// 			continue;
// 		}
// 		expr_t * tmp_mul_expr = NULL;
// 		if(expr->sup_coeff[i] < 0 || expr->inf_coeff[i]<0){
// 			if(mul_expr->size==0){
// 				tmp_mul_expr = multiply_cst_expr(pr,mul_expr, expr->inf_coeff[i],expr->sup_coeff[i]);
// 				add_cst_expr(pr,res,tmp_mul_expr);
// 				free_expr(tmp_mul_expr);
// 			}
// 			else{
// 				tmp_mul_expr = multiply_expr(pr, mul_expr, expr->inf_coeff[i],expr->sup_coeff[i]);
				
// 				add_expr(pr,res,tmp_mul_expr);
// 				free_expr(tmp_mul_expr);
				
// 			}
// 		}
// 		else{
// 			//printf("WTF2 %g %g\n",expr->inf_coeff[i],expr->sup_coeff[i]);
// 			//fflush(stdout);
// 			double tmp1, tmp2;
// 			elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,neuron_k->lb,neuron_k->ub,expr->inf_coeff[i],expr->sup_coeff[i]);
// 			if(is_lower){
// 				res->inf_cst = res->inf_cst + tmp1;
// 				res->sup_cst = res->sup_cst - tmp1;
// 			}
// 			else{
// 				res->inf_cst = res->inf_cst - tmp2;
// 				res->sup_cst = res->sup_cst + tmp2;
// 			}
// 		}
		
// 	}
// 	res->inf_cst = res->inf_cst + expr->inf_cst; 
// 	res->sup_cst = res->sup_cst + expr->sup_cst; 
// 	return res;
// }

// expr_t * lexpr_replace_bounds_affine(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons){
// 	return expr_replace_bounds_affine(pr, expr, neurons, true);
// }

// expr_t * uexpr_replace_bounds_affine(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons){
// 	return expr_replace_bounds_affine(pr, expr, neurons, false);
// }

// expr_t * expr_replace_bounds_activation(fppoly_internal_t * pr, expr_t * expr, neuron_t ** neurons, bool is_lower){
// 	size_t num_neurons = expr->size;
// 	size_t i,k;
// 	expr_t * res = alloc_expr();  
// 	res->inf_coeff = (double *)malloc(num_neurons*sizeof(double));
// 	res->sup_coeff = (double *)malloc(num_neurons*sizeof(double));
// 	res->inf_cst = expr->inf_cst;
// 	res->sup_cst = expr->sup_cst;
// 	res->type = expr->type;
// 	res->size = num_neurons;  

// 	for(i = 0; i < num_neurons; i++){
// 		if(expr->type==DENSE){
// 			k = i;
// 		}
// 		else{
// 			k = expr->dim[i];
// 		}
// 		neuron_t *neuron_k = neurons[k];
		
// 		if((expr->sup_coeff[i]==0) && (expr->inf_coeff[i]==0)){
// 			res->inf_coeff[i] = 0.0;
// 			res->sup_coeff[i] = 0.0;
// 			continue;
// 		}
// 		expr_t * mul_expr = NULL;
// 		if(is_lower){
// 			if(expr->sup_coeff[i] < 0){
// 				mul_expr = neuron_k->uexpr;
// 			}
// 			else if(expr->inf_coeff[i]<0){
// 				mul_expr = neuron_k->lexpr;
// 			}
// 		}
// 		else{
// 			if(expr->sup_coeff[i] < 0){
// 				mul_expr = neuron_k->lexpr;
// 			}
// 			else if(expr->inf_coeff[i]<0){
// 				mul_expr = neuron_k->uexpr;
// 			}
// 		}
// 		if(expr->sup_coeff[i]<0 || expr->inf_coeff[i] < 0){
// 			double lambda_inf = mul_expr->inf_coeff[0];
// 			double lambda_sup = mul_expr->sup_coeff[0];
// 			double mu_inf = mul_expr->inf_cst;
// 			double mu_sup = mul_expr->sup_cst;
// 			//res->coeff[i] = lambda*expr->coeff[i];
// 			//res->cst = res->cst + expr->coeff[i]*mu;
// 			elina_double_interval_mul_expr_coeff(pr,&res->inf_coeff[i],&res->sup_coeff[i],lambda_inf,lambda_sup,expr->inf_coeff[i],expr->sup_coeff[i]);
// 			double tmp1, tmp2;
// 			elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,mu_inf,mu_sup,expr->inf_coeff[i],expr->sup_coeff[i]);
// 			res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
// 			res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
			
// 		}
// 		else{
// 			res->inf_coeff[i] = 0.0;
// 			res->sup_coeff[i] = 0.0;
// 			double tmp1, tmp2;
// 			elina_double_interval_mul_expr_coeff(pr, &tmp1,&tmp2, neuron_k->lb, neuron_k->ub, expr->inf_coeff[i],expr->sup_coeff[i]);
// 			if(is_lower){
// 				res->inf_cst = res->inf_cst + tmp1;
// 				res->sup_cst = res->sup_cst - tmp1;
// 			}
// 			else{
// 				res->inf_cst = res->inf_cst - tmp2;
// 				res->sup_cst = res->sup_cst + tmp2;
// 			}
			
// 		}
// 	}
// 	if(expr->type==SPARSE){
// 		res->dim = (size_t*)malloc(num_neurons*sizeof(size_t));
// 		for(i=0; i < num_neurons; i++){
// 			res->dim[i] = expr->dim[i];
// 		}
// 	}
	
// 	return res;
// }

// expr_t * lexpr_replace_bounds_activation(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons){
// 	return expr_replace_bounds_activation(pr, expr, neurons, true);
// }

// expr_t * uexpr_replace_bounds_activation(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons){
// 	return expr_replace_bounds_activation(pr, expr, neurons, false);
// }


// expr_t * lexpr_replace_bounds(fppoly_internal_t * pr, expr_t * expr, neuron_t ** neurons, bool is_activation){
// 	if(is_activation){
// 		return lexpr_replace_bounds_activation(pr, expr, neurons);
// 	}
// 	else{
// 		return lexpr_replace_bounds_affine(pr, expr, neurons);
// 	}
// }

// expr_t * uexpr_replace_bounds(fppoly_internal_t * pr, expr_t * expr, neuron_t ** neurons, bool is_activation){
// 	if(is_activation){
// 		return uexpr_replace_bounds_activation(pr, expr, neurons);
// 	}
// 	else{
// 		return uexpr_replace_bounds_affine(pr, expr, neurons);
// 	}
// }

// elina_linexpr0_t *elina_linexpr0_from_expr(expr_t *expr){
// 	//assert(expr->type==DENSE);
// 	size_t size = expr->size;
// 	size_t i, k;
// 	elina_linexpr0_t * linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_DENSE,size);
// 	for(i=0; i< size; i++){
// 		if(expr->type==SPARSE){
// 			k = expr->dim[i];
// 		}
// 		else{
// 			k = i;
// 		}
// 		elina_linexpr0_set_coeff_scalar_double(linexpr0, k, expr->sup_coeff[i]);
// 	}
// 	elina_linexpr0_set_cst_scalar_double(linexpr0, expr->sup_cst);
// 	return linexpr0;
// }

#include <stdio.h>  // For fprintf, FILE*
#include <stdlib.h> // For malloc, free, realloc
#include <string.h> // For memcpy? Not strictly needed here based on code
#include <math.h>   // For fmax, fabs, INFINITY? Not strictly needed based on code shown
#include <assert.h> // For assert()
#include <fenv.h>   // For rounding mode control

#include "expr.h"   // Assuming necessary types and function declarations are here

// Include necessary headers or forward declarations for ELINA types if not in expr.h
// e.g., elina_double_interval_mul, elina_linexpr0_t, etc.

#pragma STDC FENV_ACCESS ON // Needed for fesetround

// --- Forward declaration for the base sound multiplication ---
// --- Ensure this is defined correctly elsewhere (e.g., elina_box_meetjoin.c) ---
extern void elina_double_interval_mul(double *a_inf, double *a_sup, double b_inf, double b_sup, double c_inf, double c_sup);


/* ============================================================ */
/* Sound Interval Addition with Rounding Control               */
/* ============================================================ */
// Adds interval [L1, U1] and [L2, U2] storing result in [L_res, U_res]
// Input representation: inf1 = -L1, sup1 = U1, inf2 = -L2, sup2 = U2
// Output representation: res_inf = -L_res, res_sup = U_res
static inline void elina_double_interval_add(double *res_inf, double *res_sup,
                                        double inf1, double sup1, double inf2, double sup2) {
    // Convert to [L, U] representation
    double l1 = -inf1;
    double u1 = sup1;
    double l2 = -inf2;
    double u2 = sup2;

    int original_round = fegetround(); // Save rounding mode

    // Calculate lower bound rounded down
    fesetround(FE_DOWNWARD);
    double res_l = l1 + l2;

    // Calculate upper bound rounded up
    fesetround(FE_UPWARD);
    double res_u = u1 + u2;

    // Convert back to negated lower bound representation, rounded up
    *res_inf = -res_l; // Rounding -res_l UP
    *res_sup = res_u;  // Already rounded UP

    fesetround(original_round); // Restore rounding mode
}


// --- Removed unsound wrapper functions ---
// void elina_double_interval_add_expr_coeff(...)
// void elina_double_interval_add_cst_coeff(...)
// void elina_double_interval_mul_expr_coeff(...)
// void elina_double_interval_mul_cst_coeff(...)


void expr_fprint(FILE * stream, expr_t *expr){
    if(expr == NULL) {
        fprintf(stream, "(NULL expression)\n");
        return;
    }
    if (expr->size == 0 || expr->inf_coeff == NULL || expr->sup_coeff == NULL) {
         // Use the representation [-inf, sup] for the constant
        fprintf(stream,"+ [%g, %g]\n", -expr->inf_cst, expr->sup_cst);
        return;
    }

    size_t size = expr->size;
    size_t i;
    bool first_term = true;

    for(i=0; i < size; i++){
        // Use tolerance to avoid printing near-zero terms if desired
        if (fabs(expr->inf_coeff[i]) > 1e-12 || fabs(expr->sup_coeff[i]) > 1e-12) {
            if (!first_term) {
                fprintf(stream, "+ ");
            }
            if(expr->type==DENSE){
                 // Use the representation [-inf, sup] for coefficient
                fprintf(stream, "[%g, %g]x%zu ", -expr->inf_coeff[i], expr->sup_coeff[i], i);
            }
            else{
                 // Use the representation [-inf, sup] for coefficient
                fprintf(stream, "[%g, %g]x%zu ", -expr->inf_coeff[i], expr->sup_coeff[i], expr->dim[i]);
            }
            first_term = false;
        }
    }

    // Print constant term if it's non-zero or if no variable terms were printed
    if (first_term || fabs(expr->inf_cst) > 1e-12 || fabs(expr->sup_cst) > 1e-12) {
        if (!first_term) {
             fprintf(stream, "+ ");
        }
         // Use the representation [-inf, sup] for the constant
        fprintf(stream,"[%g, %g]\n", -expr->inf_cst, expr->sup_cst);
    } else if (first_term) {
        // If all coefficients and constant are zero
         fprintf(stream,"[0, 0]\n");
    } else {
        // If coefficients existed but constant is zero, just add newline
        fprintf(stream, "\n");
    }
}


void expr_print(expr_t * expr){
    expr_fprint(stdout, expr);
}

expr_t * alloc_expr(void){
    expr_t *expr = (expr_t *)malloc(sizeof(expr_t));
    if (!expr) return NULL; // Allocation check
    expr->inf_coeff = NULL;
    expr->sup_coeff = NULL;
    expr->dim = NULL;
    expr->inf_cst = 0.0; // Initialize constants
    expr->sup_cst = 0.0;
    expr->size = 0;
    expr->type = SPARSE; // Default type? Or should be specified?
    return expr;
}

// Creates dense expression [L, U] = coeff * x + cst
// Input coeff is scalar, cst is scalar.
// Output uses negated lower bound representation.
expr_t * create_dense_expr(double *coeff, double cst, size_t size){
    expr_t *expr = (expr_t *)malloc(sizeof(expr_t));
    if (!expr) return NULL;
    expr->inf_coeff = (double *)malloc(size*sizeof(double));
    expr->sup_coeff = (double *)malloc(size*sizeof(double));
    if (size > 0 && (!expr->inf_coeff || !expr->sup_coeff)) { // Allocation check
        free(expr->inf_coeff); free(expr->sup_coeff); free(expr);
        return NULL;
    }
    expr->dim= NULL;
    size_t i;
    expr->size = size;
    expr->inf_cst = -cst; // Negated lower constant
    expr->sup_cst = cst;  // Upper constant
    expr->type = DENSE;
    for(i=0; i < size; i++){
        expr->inf_coeff[i] = -coeff[i]; // Negated lower coeff
        expr->sup_coeff[i] = coeff[i];  // Upper coeff
    }
    return expr;
}


// Creates constant expression [l, u] (input l, u are mathematical bounds)
// Output uses negated lower bound representation.
expr_t * create_cst_expr(double l, double u){
    expr_t *expr = (expr_t*)malloc(sizeof(expr_t));
    if (!expr) return NULL;
    expr->inf_coeff = NULL;
    expr->sup_coeff = NULL;
    expr->dim = NULL;
    expr->type = SPARSE; // Type for constant expression?
    expr->size = 0;
    expr->inf_cst = -l; // Negated lower constant
    expr->sup_cst = u;  // Upper constant
    return expr;
}

// Creates sparse expression [L, U] = sum(coeff_i * x_dim[i]) + cst
// Input coeff_i are scalar, cst is scalar.
// Output uses negated lower bound representation.
expr_t * create_sparse_expr(double *coeff, double cst, size_t *dim, size_t size){
    expr_t *expr = (expr_t *)malloc(sizeof(expr_t));
    if (!expr) return NULL;
    if(size > 0){
        expr->inf_coeff = (double *)malloc(size*sizeof(double));
        expr->sup_coeff = (double *)malloc(size*sizeof(double));
        expr->dim = (size_t *)malloc(size*sizeof(size_t));
        if (!expr->inf_coeff || !expr->sup_coeff || !expr->dim) { // Allocation check
            free(expr->inf_coeff); free(expr->sup_coeff); free(expr->dim); free(expr);
            return NULL;
        }
    }
    else{
        expr->inf_coeff = NULL;
        expr->sup_coeff = NULL;
        expr->dim = NULL;
    }
    size_t i;
    expr->size = size;
    expr->inf_cst = -cst; // Negated lower constant
    expr->sup_cst = cst;  // Upper constant
    expr->type = SPARSE;
    for(i=0; i < size; i++){
        expr->inf_coeff[i] = -coeff[i]; // Negated lower coeff
        expr->sup_coeff[i] = coeff[i];  // Upper coeff
        expr->dim[i] = dim[i];
    }
    return expr;
}


void free_expr(expr_t *expr){
    if (!expr) return; // Handle NULL pointer
    free(expr->inf_coeff); // free(NULL) is safe
    free(expr->sup_coeff);
    // Only free dim if it was allocated (sparse and size > 0)
    if(expr->type == SPARSE && expr->dim != NULL){
        free(expr->dim);
    }
    expr->inf_coeff = NULL;
    expr->sup_coeff = NULL;
    expr->dim = NULL;
    free(expr);
    // expr = NULL; // Setting caller's pointer to NULL is not possible here
}

// Copies only the constant part of an expression
expr_t * copy_cst_expr(expr_t *src){
    if (!src) return NULL;
    expr_t *dst = (expr_t *)malloc(sizeof(expr_t));
    if (!dst) return NULL;
    dst->inf_coeff = NULL;
    dst->sup_coeff = NULL;
    dst->inf_cst = src->inf_cst;
    dst->sup_cst = src->sup_cst;
    dst->type = src->type; // Keep original type? Or force SPARSE?
    dst->dim = NULL;
    dst->size = 0; // Size should be 0 for constant expression
    return dst;
}


// Copies the full expression
expr_t * copy_expr(expr_t *src){
    if (!src) return NULL;
    expr_t *dst = (expr_t *)malloc(sizeof(expr_t));
    if (!dst) return NULL;

    dst->inf_cst = src->inf_cst;
    dst->sup_cst = src->sup_cst;
    dst->type = src->type;
    dst->size = src->size;
    dst->inf_coeff = NULL;
    dst->sup_coeff = NULL;
    dst->dim = NULL;

    if (src->size > 0) {
        dst->inf_coeff = (double *)malloc(src->size*sizeof(double));
        dst->sup_coeff = (double *)malloc(src->size*sizeof(double));
        if (!dst->inf_coeff || !dst->sup_coeff) { // Allocation check
            free(dst->inf_coeff); free(dst->sup_coeff); free(dst);
            return NULL;
        }
        memcpy(dst->inf_coeff, src->inf_coeff, src->size * sizeof(double));
        memcpy(dst->sup_coeff, src->sup_coeff, src->size * sizeof(double));

        if(src->type==SPARSE && src->dim != NULL){ // Check src->dim too
            dst->dim = (size_t *)malloc(src->size*sizeof(size_t));
            if (!dst->dim) { // Allocation check
                 free(dst->inf_coeff); free(dst->sup_coeff); free(dst);
                 return NULL;
            }
            memcpy(dst->dim, src->dim, src->size * sizeof(size_t));
        }
    }
    return dst;
}

// Concretizes part of a dense expression using interval arithmetic
expr_t* concretize_dense_sub_expr(fppoly_internal_t *pr, expr_t * expr, double *inf, double *sup, size_t start, size_t size){
    if (!expr || expr->type != DENSE || start > expr->size || size > expr->size) {
        // Basic error checking
        return NULL;
    }
    expr_t * res = (expr_t *)malloc(sizeof(expr_t));
    if (!res) return NULL;

    res->inf_coeff = NULL;
    res->sup_coeff = NULL;
    res->dim = NULL;
    res->inf_cst = expr->inf_cst; // Start with original constant
    res->sup_cst = expr->sup_cst;
    res->type = DENSE; // Resulting expression is dense up to 'start'
    res->size = start;

    if (start > 0) {
        res->inf_coeff = (double *)malloc(start*sizeof(double));
        res->sup_coeff = (double *)malloc(start*sizeof(double));
         if (!res->inf_coeff || !res->sup_coeff) {
             free(res->inf_coeff); free(res->sup_coeff); free(res);
             return NULL;
         }
        // Copy the coefficients that remain symbolic
        memcpy(res->inf_coeff, expr->inf_coeff, start * sizeof(double));
        memcpy(res->sup_coeff, expr->sup_coeff, start * sizeof(double));
    }

    // Concretize terms from 'start' to 'size'
    for(size_t i = start; i < size; i++){
        double term_inf, term_sup;
        // Use sound interval multiplication
        elina_double_interval_mul(&term_inf, &term_sup,
                                 expr->inf_coeff[i], expr->sup_coeff[i], // Coeff interval
                                 inf[i-start], sup[i-start]);          // Var interval

        // Use sound interval addition to add to constant term
        elina_double_interval_add(&res->inf_cst, &res->sup_cst,
                                  res->inf_cst, res->sup_cst,          // Current constant interval
                                  term_inf, term_sup);                 // Term interval
    }

    return res;
}


// --- Merge sort for sparse expressions (unchanged, seems correct) ---
void merge_sparse_expr(expr_t *expr, size_t l, size_t m, size_t r) {
    size_t i, j, k; // Use size_t for indices
    size_t n1 = m - l + 1;
    size_t n2 = r - m;

    /* create temp arrays */
    // Check malloc results
    size_t *L_dim = (size_t *)malloc(n1*sizeof(size_t));
    size_t *R_dim = (size_t *)malloc(n2*sizeof(size_t));
    double *L_inf = (double *)malloc(n1*sizeof(double));
    double *R_inf = (double *)malloc(n2*sizeof(double));
    double *L_sup = (double *)malloc(n1*sizeof(double));
    double *R_sup = (double *)malloc(n2*sizeof(double));

    if (!L_dim || !R_dim || !L_inf || !R_inf || !L_sup || !R_sup) {
        // Handle allocation failure
        free(L_dim); free(R_dim); free(L_inf); free(R_inf); free(L_sup); free(R_sup);
        // Maybe signal error? For now, just return.
        fprintf(stderr, "ERROR: Allocation failed in merge_sparse_expr\n");
        return;
    }

    /* Copy data to temp arrays */
    for (i = 0; i < n1; i++) {
        L_dim[i] = expr->dim[l + i];
        L_inf[i] = expr->inf_coeff[l + i];
        L_sup[i] = expr->sup_coeff[l + i];
    }
    for (j = 0; j < n2; j++) {
        R_dim[j] = expr->dim[m + 1 + j];
        R_inf[j] = expr->inf_coeff[m + 1 + j];
        R_sup[j] = expr->sup_coeff[m + 1 + j];
    }

    /* Merge the temp arrays back */
    i = 0; j = 0; k = l;
    while (i < n1 && j < n2) {
        if (L_dim[i] <= R_dim[j]) {
            expr->dim[k] = L_dim[i];
            expr->inf_coeff[k] = L_inf[i];
            expr->sup_coeff[k] = L_sup[i];
            i++;
        } else {
            expr->dim[k] = R_dim[j];
            expr->inf_coeff[k] = R_inf[j];
            expr->sup_coeff[k] = R_sup[j];
            j++;
        }
        k++;
    }

    /* Copy remaining elements */
    while (i < n1) {
        expr->dim[k] = L_dim[i];
        expr->inf_coeff[k] = L_inf[i];
        expr->sup_coeff[k] = L_sup[i];
        i++; k++;
    }
    while (j < n2) {
        expr->dim[k] = R_dim[j];
        expr->inf_coeff[k] = R_inf[j];
        expr->sup_coeff[k] = R_sup[j];
        j++; k++;
    }

    // Free temp arrays
    free(L_dim); free(R_dim); free(L_inf); free(R_inf); free(L_sup); free(R_sup);
}

void merge_sort_sparse_expr(expr_t *expr, size_t l, size_t r) {
    // Check if valid range and size > 0
    if (expr == NULL || expr->size == 0 || l >= r) {
        return;
    }
    if (l < r) { // Ensure l < r to avoid potential issues with m+1 later
        size_t m = l + (r - l) / 2;
        merge_sort_sparse_expr(expr, l, m);
        merge_sort_sparse_expr(expr, m + 1, r);
        merge_sparse_expr(expr, l, m, r);
    }
}

void sort_sparse_expr(expr_t *expr){
    if (expr && expr->type == SPARSE && expr->size > 1) { // Only sort if sparse and size > 1
        merge_sort_sparse_expr(expr, 0, expr->size - 1);
    }
}


// Multiplies expression by an interval [mul_l, mul_u] (negated lb rep)
expr_t * multiply_expr(fppoly_internal_t *pr, expr_t *expr, double mul_inf, double mul_sup){
    if (!expr) return NULL;
    expr_t * res = alloc_expr();
    if (!res) return NULL;

    res->type = expr->type;
    res->size = expr->size;
    res->inf_coeff = NULL;
    res->sup_coeff = NULL;
    res->dim = NULL;

    if (expr->size > 0) {
        res->inf_coeff = (double *)malloc(expr->size*sizeof(double));
        res->sup_coeff = (double *)malloc(expr->size*sizeof(double));
        if (!res->inf_coeff || !res->sup_coeff) { // Allocation check
            free(res->inf_coeff); free(res->sup_coeff); free(res);
            return NULL;
        }
        for(size_t i = 0; i < expr->size; i++){
            // Sound multiplication for coefficients
            elina_double_interval_mul(&res->inf_coeff[i], &res->sup_coeff[i],
                                     mul_inf, mul_sup,                     // Multiplier interval
                                     expr->inf_coeff[i], expr->sup_coeff[i]); // Original coeff interval
        }

        if(expr->type == SPARSE && expr->dim != NULL){ // Check expr->dim too
            res->dim = (size_t*)malloc(expr->size*sizeof(size_t));
            if (!res->dim) { // Allocation check
                 free(res->inf_coeff); free(res->sup_coeff); free(res);
                 return NULL;
            }
            memcpy(res->dim, expr->dim, expr->size * sizeof(size_t));
        }
    }

    // Sound multiplication for constant term
    elina_double_interval_mul(&res->inf_cst, &res->sup_cst,
                             mul_inf, mul_sup,         // Multiplier interval
                             expr->inf_cst, expr->sup_cst); // Original const interval

    return res;
}

// Multiplies a constant expression by an interval
expr_t * multiply_cst_expr(fppoly_internal_t *pr, expr_t *expr, double mul_inf, double mul_sup){
     if (!expr || expr->size != 0) return NULL; // Should only be called on constant expressions
    expr_t * res = alloc_expr();
    if (!res) return NULL;

    res->inf_coeff = NULL;
    res->sup_coeff = NULL;
    res->dim = NULL;
    res->type = expr->type; // Preserve type?
    res->size = 0;

    // Sound multiplication for constant term
    elina_double_interval_mul(&res->inf_cst, &res->sup_cst,
                             mul_inf, mul_sup,         // Multiplier interval
                             expr->inf_cst, expr->sup_cst); // Original const interval
    return res;
}


// Adds exprB (constant) to exprA (potentially non-constant)
void add_cst_expr(fppoly_internal_t *pr, expr_t * exprA, expr_t *exprB){
    if (!exprA || !exprB || exprB->size != 0) return; // Basic checks
    // Sound interval addition for constant terms
    elina_double_interval_add(&exprA->inf_cst, &exprA->sup_cst,
                              exprA->inf_cst, exprA->sup_cst,
                              exprB->inf_cst, exprB->sup_cst);
    return;
}

// Adds exprB to exprA (modifies exprA)
void add_expr(fppoly_internal_t *pr, expr_t * exprA, expr_t * exprB){
    if (!exprA || !exprB) return;

    // Add constant parts soundly
    elina_double_interval_add(&exprA->inf_cst, &exprA->sup_cst,
                              exprA->inf_cst, exprA->sup_cst,
                              exprB->inf_cst, exprB->sup_cst);

    // If exprB has no coefficients, we are done
    if (exprB->size == 0) {
        return;
    }

    // If exprA initially had no coefficients
    if (exprA->size == 0) {
        exprA->size = exprB->size;
        exprA->type = exprB->type;
        if (exprB->size > 0) {
            exprA->inf_coeff = (double *)malloc(exprB->size * sizeof(double));
            exprA->sup_coeff = (double *)malloc(exprB->size * sizeof(double));
            if (!exprA->inf_coeff || !exprA->sup_coeff) { /* Handle error */ return; }
            memcpy(exprA->inf_coeff, exprB->inf_coeff, exprB->size * sizeof(double));
            memcpy(exprA->sup_coeff, exprB->sup_coeff, exprB->size * sizeof(double));

            if (exprB->type == SPARSE && exprB->dim != NULL) {
                exprA->dim = (size_t *)malloc(exprB->size * sizeof(size_t));
                if (!exprA->dim) { /* Handle error */ return; }
                memcpy(exprA->dim, exprB->dim, exprB->size * sizeof(size_t));
            }
        }
        return;
    }

    // --- Both exprA and exprB have coefficients ---
    size_t sizeA = exprA->size;
    size_t sizeB = exprB->size;
    size_t i = 0, k = 0, l = 0; // i for A, k for B, l for result index

    double *new_inf_coeff = NULL;
    double *new_sup_coeff = NULL;
    size_t *new_dim = NULL;
    size_t new_size = 0;
    int new_type = DENSE; // Start assuming dense if mixing

    // Case 1: Both Dense
    if (exprA->type == DENSE && exprB->type == DENSE) {
        assert(sizeA == sizeB); // Should have same dimension for dense add
        new_size = sizeA;
        new_type = DENSE;
        new_inf_coeff = exprA->inf_coeff; // Modify A in place
        new_sup_coeff = exprA->sup_coeff;
        for (i = 0; i < new_size; i++) {
            elina_double_interval_add(&new_inf_coeff[i], &new_sup_coeff[i],
                                      exprA->inf_coeff[i], exprA->sup_coeff[i],
                                      exprB->inf_coeff[i], exprB->sup_coeff[i]);
        }
        // No need to reallocate or change dim/type/size
        return; // Done
    }
    // Case 2: A is Dense, B is Sparse
    else if (exprA->type == DENSE && exprB->type == SPARSE) {
        new_size = sizeA;
        new_type = DENSE;
        new_inf_coeff = exprA->inf_coeff; // Modify A in place
        new_sup_coeff = exprA->sup_coeff;
        for (k = 0; k < sizeB; k++) {
            size_t dimB = exprB->dim[k];
            if (dimB < sizeA) { // Check dimension bounds
                 elina_double_interval_add(&new_inf_coeff[dimB], &new_sup_coeff[dimB],
                                          exprA->inf_coeff[dimB], exprA->sup_coeff[dimB],
                                          exprB->inf_coeff[k], exprB->sup_coeff[k]);
            } else {
                fprintf(stderr, "Warning: Sparse dimension %zu out of bounds for dense add (size %zu).\n", dimB, sizeA);
            }
        }
        // No need to reallocate or change dim/type/size
        return; // Done
    }
    // Case 3: A is Sparse, B is Dense
    else if (exprA->type == SPARSE && exprB->type == DENSE) {
        new_size = sizeB;
        new_type = DENSE;
        new_inf_coeff = (double *)malloc(new_size * sizeof(double));
        new_sup_coeff = (double *)malloc(new_size * sizeof(double));
        if (!new_inf_coeff || !new_sup_coeff) { /* Handle error */ return; }

        memcpy(new_inf_coeff, exprB->inf_coeff, new_size * sizeof(double)); // Start with B
        memcpy(new_sup_coeff, exprB->sup_coeff, new_size * sizeof(double));

        for (i = 0; i < sizeA; i++) {
            size_t dimA = exprA->dim[i];
            if (dimA < new_size) { // Check bounds
                 elina_double_interval_add(&new_inf_coeff[dimA], &new_sup_coeff[dimA],
                                          new_inf_coeff[dimA], new_sup_coeff[dimA], // Current value (from B)
                                          exprA->inf_coeff[i], exprA->sup_coeff[i]); // Value from A
            } else {
                 fprintf(stderr, "Warning: Sparse dimension %zu out of bounds for dense add (size %zu).\n", dimA, new_size);
            }
        }
        // Update exprA
        free(exprA->inf_coeff); free(exprA->sup_coeff); free(exprA->dim);
        exprA->inf_coeff = new_inf_coeff;
        exprA->sup_coeff = new_sup_coeff;
        exprA->dim = NULL; // Now dense
        exprA->size = new_size;
        exprA->type = new_type;
        return; // Done
    }
    // Case 4: Both Sparse
    else { // Both SPARSE
        new_type = SPARSE;
        // Allocate worst-case size, then realloc later
        new_inf_coeff = (double *)malloc((sizeA + sizeB) * sizeof(double));
        new_sup_coeff = (double *)malloc((sizeA + sizeB) * sizeof(double));
        new_dim = (size_t *)malloc((sizeA + sizeB) * sizeof(size_t));
        if (!new_inf_coeff || !new_sup_coeff || !new_dim) { /* Handle error */ return; }

        i = 0; k = 0; l = 0; // i for A, k for B, l for result index
        while (i < sizeA && k < sizeB) {
            if (exprA->dim[i] < exprB->dim[k]) {
                new_dim[l] = exprA->dim[i];
                new_inf_coeff[l] = exprA->inf_coeff[i];
                new_sup_coeff[l] = exprA->sup_coeff[i];
                i++;
            } else if (exprB->dim[k] < exprA->dim[i]) {
                new_dim[l] = exprB->dim[k];
                new_inf_coeff[l] = exprB->inf_coeff[k];
                new_sup_coeff[l] = exprB->sup_coeff[k];
                k++;
            } else { // Dimensions match, add coefficients
                new_dim[l] = exprA->dim[i];
                elina_double_interval_add(&new_inf_coeff[l], &new_sup_coeff[l],
                                          exprA->inf_coeff[i], exprA->sup_coeff[i],
                                          exprB->inf_coeff[k], exprB->sup_coeff[k]);
                i++;
                k++;
            }
            l++;
        }
        // Copy remaining elements from A
        while (i < sizeA) {
            new_dim[l] = exprA->dim[i];
            new_inf_coeff[l] = exprA->inf_coeff[i];
            new_sup_coeff[l] = exprA->sup_coeff[i];
            i++; l++;
        }
        // Copy remaining elements from B
        while (k < sizeB) {
            new_dim[l] = exprB->dim[k];
            new_inf_coeff[l] = exprB->inf_coeff[k];
            new_sup_coeff[l] = exprB->sup_coeff[k];
            k++; l++;
        }
        new_size = l; // Actual size of the result

        // Update exprA
        free(exprA->inf_coeff); free(exprA->sup_coeff); free(exprA->dim);
        // Realloc to actual size
        exprA->inf_coeff = (double*)realloc(new_inf_coeff, new_size * sizeof(double));
        exprA->sup_coeff = (double*)realloc(new_sup_coeff, new_size * sizeof(double));
        exprA->dim = (size_t*)realloc(new_dim, new_size * sizeof(size_t));
        // Check realloc results
        if (new_size > 0 && (!exprA->inf_coeff || !exprA->sup_coeff || !exprA->dim)) {
             /* Handle realloc error - exprA is now in a bad state */
             fprintf(stderr, "ERROR: Realloc failed in add_expr (sparse+sparse)\n");
             // Maybe try to restore or signal error?
             // For now, free potentially allocated new arrays if realloc failed
             if(exprA->inf_coeff != new_inf_coeff) free(new_inf_coeff);
             if(exprA->sup_coeff != new_sup_coeff) free(new_sup_coeff);
             if(exprA->dim != new_dim) free(new_dim);
             // Set exprA to a safe (e.g., empty) state if possible
             exprA->inf_coeff = NULL; exprA->sup_coeff = NULL; exprA->dim = NULL; exprA->size = 0;
             return;
        }
        exprA->size = new_size;
        exprA->type = new_type;
        return; // Done
    }
}


// Extracts sub-expression relevant for a specific predecessor in concatenation
// Assumes expr is defined over the concatenated layer's indices
expr_t * extract_subexpr_concatenate(expr_t * expr, size_t index, size_t* C, size_t num_neurons, size_t num_channels){
    if (!expr || !C) return NULL;
    expr_t * res = alloc_expr();
    if (!res) return NULL;

    res->inf_cst = 0.0; // Sub-expression constant starts at 0
    res->sup_cst = 0.0;

    // Calculate height * width
    assert(num_channels > 0);
    size_t hw = num_neurons / num_channels;

    // Calculate the starting channel index for the desired predecessor 'index'
    size_t channel_offset = 0;
    for(size_t i = 0; i < index; i++){
        channel_offset += C[i];
    }
    size_t num_pred_channels = C[index]; // Number of channels from this predecessor

    if (expr->type == DENSE) {
        res->type = DENSE;
        size_t num_neurons_in_pred = num_pred_channels * hw;
        res->size = num_neurons_in_pred;
        if (res->size == 0) return res; // Empty result

        res->inf_coeff = (double *)malloc(res->size * sizeof(double));
        res->sup_coeff = (double *)malloc(res->size * sizeof(double));
        if (!res->inf_coeff || !res->sup_coeff) { free_expr(res); return NULL; }

        size_t res_idx = 0;
        // Iterate through spatial locations (height * width)
        for(size_t i = 0; i < hw; i++){
            // Iterate through the channels belonging to this predecessor
            for(size_t k = 0; k < num_pred_channels; k++){
                size_t original_dense_idx = i * num_channels + channel_offset + k;
                if (original_dense_idx < expr->size) { // Bounds check
                    res->inf_coeff[res_idx] = expr->inf_coeff[original_dense_idx];
                    res->sup_coeff[res_idx] = expr->sup_coeff[original_dense_idx];
                } else {
                    // Should not happen if inputs are correct, but handle defensively
                    res->inf_coeff[res_idx] = 0.0;
                    res->sup_coeff[res_idx] = 0.0;
                }
                res_idx++;
            }
        }
    } else { // expr->type == SPARSE
        res->type = SPARSE;

        // Determine the range of original dimensions for this predecessor
        // [hw_idx * num_channels + channel_offset, hw_idx * num_channels + channel_offset + num_pred_channels)
        // We need to count how many terms fall into these ranges and map their dimensions.

        // Count terms belonging to this predecessor
        size_t count = 0;
        for (size_t i = 0; i < expr->size; i++) {
            size_t original_dim = expr->dim[i];
            size_t hw_idx = original_dim / num_channels;
            size_t channel_idx = original_dim % num_channels;
            if (channel_idx >= channel_offset && channel_idx < channel_offset + num_pred_channels) {
                count++;
            }
        }

        res->size = count;
        if (res->size == 0) return res; // Empty result

        res->inf_coeff = (double *)malloc(res->size * sizeof(double));
        res->sup_coeff = (double *)malloc(res->size * sizeof(double));
        res->dim = (size_t *)malloc(res->size * sizeof(size_t));
        if (!res->inf_coeff || !res->sup_coeff || !res->dim) { free_expr(res); return NULL; }

        size_t res_idx = 0;
        for (size_t i = 0; i < expr->size; i++) {
            size_t original_dim = expr->dim[i];
            size_t hw_idx = original_dim / num_channels;
            size_t channel_idx = original_dim % num_channels;

            // Check if this term belongs to the current predecessor
            if (channel_idx >= channel_offset && channel_idx < channel_offset + num_pred_channels) {
                res->inf_coeff[res_idx] = expr->inf_coeff[i];
                res->sup_coeff[res_idx] = expr->sup_coeff[i];
                // Map original dimension to the predecessor's local dimension space
                // Local channel index: channel_idx - channel_offset
                // Local dimension: hw_idx * num_pred_channels + (channel_idx - channel_offset)
                res->dim[res_idx] = hw_idx * num_pred_channels + (channel_idx - channel_offset);
                res_idx++;
            }
        }
        assert(res_idx == res->size); // Sanity check
    }
    return res;
}

// Performs back-substitution through an affine layer
expr_t * expr_replace_bounds_affine(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons, bool is_lower){
    if (!expr) return NULL;
    if (expr->size == 0) { // If expression is constant
        return copy_expr(expr); // Just copy the constant
    }
    if (expr->inf_coeff == NULL || expr->sup_coeff == NULL) { // Invalid expression
        return NULL; // Or return an empty/error expression
    }

    size_t num_terms = expr->size;
    size_t i, k;
    expr_t * res = NULL; // Resulting expression after substitution
    bool first_term = true;

    for(i = 0; i < num_terms; i++){
        size_t neuron_idx;
        if(expr->type == DENSE){
            neuron_idx = i;
        } else { // SPARSE
            neuron_idx = expr->dim[i];
        }

        // Check bounds for neuron_idx if needed, assume valid for now
        neuron_t *neuron_k = neurons[neuron_idx];
        if (!neuron_k) { /* Handle error */ continue; } // Neuron not found?

        // Interval for the current term's coefficient: [coeff_l, coeff_u]
        // Remember inf_coeff stores -coeff_l
        double coeff_inf = expr->inf_coeff[i];
        double coeff_sup = expr->sup_coeff[i];

        // Skip if coefficient interval is zero
        if (fabs(coeff_inf) < 1e-12 && fabs(coeff_sup) < 1e-12) {
            continue;
        }

        // Get the expression for neuron_k (lexpr or uexpr)
        expr_t *sub_expr = NULL;
        // Determine which expression (lexpr/uexpr) from neuron_k to substitute
        // Based on the sign of the coefficient interval [coeff_l, coeff_u]
        // Note: inf < 0 means lower bound -inf > 0; sup < 0 means upper bound sup < 0
        if (is_lower) { // Calculating lower bound of the result
            if (coeff_sup < 0) {         // Upper coeff is negative (so coeff interval is negative) -> Use uexpr
                sub_expr = neuron_k->uexpr;
            } else if (coeff_inf < 0) {  // Lower coeff is negative (so coeff interval crosses zero or is positive)
                                         // If interval crosses zero, need careful handling. This simple check assumes lexpr if not strictly negative.
                sub_expr = neuron_k->lexpr;
            } else {                     // Coeff interval is positive -> Use lexpr
                 sub_expr = neuron_k->lexpr;
            }
        } else { // Calculating upper bound of the result
            if (coeff_sup < 0) {         // Coeff interval is negative -> Use lexpr
                sub_expr = neuron_k->lexpr;
            } else if (coeff_inf < 0) {  // Coeff interval crosses zero or is positive -> Use uexpr
                sub_expr = neuron_k->uexpr;
            } else {                     // Coeff interval is positive -> Use uexpr
                 sub_expr = neuron_k->uexpr;
            }
        }

        if (!sub_expr) { /* Handle error - neuron expression missing? */ continue; }

        expr_t * term_result = NULL;
        // Multiply the chosen sub_expr by the coefficient interval [coeff_l, coeff_u]
        term_result = multiply_expr(pr, sub_expr, coeff_inf, coeff_sup);

        if (!term_result) { /* Handle error */ continue; }

        // Add the result of this term to the overall result expression 'res'
        if (first_term) {
            res = term_result; // Assign the first term directly
            first_term = false;
        } else {
            add_expr(pr, res, term_result); // Add subsequent terms
            free_expr(term_result); // Free the temporary term result
        }
    } // End loop over terms

    // If no terms were processed (e.g., all coefficients zero), create a constant expression
    if (first_term) {
        res = create_cst_expr(-expr->inf_cst, expr->sup_cst); // Use original constant, convert representation
    } else {
        // Add the original constant term of 'expr' to the final result
         expr_t * const_expr = create_cst_expr(-expr->inf_cst, expr->sup_cst); // Convert representation
         if (const_expr) {
            add_cst_expr(pr, res, const_expr);
            free_expr(const_expr);
         }
    }
    return res;
}

expr_t * lexpr_replace_bounds_affine(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons){
    return expr_replace_bounds_affine(pr, expr, neurons, true);
}

expr_t * uexpr_replace_bounds_affine(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons){
    return expr_replace_bounds_affine(pr, expr, neurons, false);
}


// Performs back-substitution through an activation layer (like ReLU)
expr_t * expr_replace_bounds_activation(fppoly_internal_t * pr, expr_t * expr, neuron_t ** neurons, bool is_lower){
     if (!expr) return NULL;
     if (expr->size == 0) { // If expression is constant
        return copy_expr(expr);
    }
     if (expr->inf_coeff == NULL || expr->sup_coeff == NULL) { // Invalid expression
        return NULL;
    }

    size_t num_terms = expr->size;
    size_t i, k;
    // IMPORTANT: Result expression depends ONLY on the predecessor of the activation neurons
    // We assume the activation layer 'neurons' has a single predecessor layer
    // Find the predecessor layer index (should be stored or passed) - Assuming neuron_t doesn't store it directly
    // This function signature might need adjustment if predecessor info isn't easily accessible

    // Let's create the result expression based on the predecessor's size and type
    // This is complex - for now, we'll try to modify the current expression structure,
    // although ideally, the result should be defined over the predecessor's variables.
    // We will build the new expression term by term.

    expr_t * res = NULL; // Initialize result
    bool first_res_term = true;

    // Start with the constant term of the original expression
    expr_t * current_const = create_cst_expr(-expr->inf_cst, expr->sup_cst); // Convert rep
    if (!current_const) return NULL;


    for(i = 0; i < num_terms; i++){
        size_t neuron_idx;
        if(expr->type == DENSE){
            neuron_idx = i;
        } else { // SPARSE
            neuron_idx = expr->dim[i];
        }
        neuron_t *neuron_k = neurons[neuron_idx]; // This is the activation neuron
        if (!neuron_k) continue;

        // Coefficient interval [coeff_l, coeff_u] for term i
        double coeff_inf = expr->inf_coeff[i];
        double coeff_sup = expr->sup_coeff[i];

        // Skip zero coefficients
        if (fabs(coeff_inf) < 1e-12 && fabs(coeff_sup) < 1e-12) {
            continue;
        }

        // Get the appropriate symbolic expression (lexpr/uexpr) from the activation neuron
        expr_t *sub_expr = NULL; // This is the ReLU approximation: lambda * x_pre + mu
        if (is_lower) {
            if (coeff_sup < 0) sub_expr = neuron_k->uexpr; // Use upper approx (lambda_u, mu_u)
            else               sub_expr = neuron_k->lexpr; // Use lower approx (lambda_l, mu_l)
        } else { // Upper bound
            if (coeff_sup < 0) sub_expr = neuron_k->lexpr; // Use lower approx (lambda_l, mu_l)
            else               sub_expr = neuron_k->uexpr; // Use upper approx (lambda_u, mu_u)
        }

        if (!sub_expr || sub_expr->size != 1 || sub_expr->inf_coeff == NULL || sub_expr->sup_coeff == NULL) {
            // Activation approximation should be form [lambda_l, lambda_u] * x_pre + [mu_l, mu_u]
            // Represented with size 1, pointing to the pre-activation neuron
             fprintf(stderr, "Warning: Invalid activation expression format for neuron %zu.\n", neuron_idx);
             continue; // Skip this term or handle error
        }

        // Extract lambda = [lambda_l, lambda_u] and mu = [mu_l, mu_u] from sub_expr
        double lambda_inf = sub_expr->inf_coeff[0]; // Stores -lambda_l
        double lambda_sup = sub_expr->sup_coeff[0]; // Stores lambda_u
        double mu_inf = sub_expr->inf_cst;          // Stores -mu_l
        double mu_sup = sub_expr->sup_cst;          // Stores mu_u

        // Calculate term: [coeff_l, coeff_u] * ([lambda_l, lambda_u] * x_pre + [mu_l, mu_u])
        // = ([coeff_l, coeff_u] * [lambda_l, lambda_u]) * x_pre + ([coeff_l, coeff_u] * [mu_l, mu_u])

        // 1. Calculate new coefficient for x_pre: new_lambda = coeff * lambda
        double new_lambda_inf, new_lambda_sup;
        elina_double_interval_mul(&new_lambda_inf, &new_lambda_sup,
                                 coeff_inf, coeff_sup, lambda_inf, lambda_sup);

        // 2. Calculate new constant contribution: new_mu = coeff * mu
        double new_mu_inf, new_mu_sup;
        elina_double_interval_mul(&new_mu_inf, &new_mu_sup,
                                 coeff_inf, coeff_sup, mu_inf, mu_sup);

        // Add new_mu to the overall constant term
        elina_double_interval_add(&current_const->inf_cst, &current_const->sup_cst,
                                  current_const->inf_cst, current_const->sup_cst,
                                  new_mu_inf, new_mu_sup);

        // Create a temporary expression for the new_lambda * x_pre term
        // The dimension should be the index of the pre-activation neuron (sub_expr->dim[0])
        expr_t * term_expr = alloc_expr();
         if (!term_expr) continue; // Handle alloc failure
         term_expr->size = 1;
         term_expr->type = SPARSE; // Result is sparse in terms of pre-activation neuron
         term_expr->inf_coeff = (double*)malloc(sizeof(double));
         term_expr->sup_coeff = (double*)malloc(sizeof(double));
         term_expr->dim = (size_t*)malloc(sizeof(size_t));
         if (!term_expr->inf_coeff || !term_expr->sup_coeff || !term_expr->dim) {
            free_expr(term_expr); continue;
         }
         term_expr->inf_coeff[0] = new_lambda_inf;
         term_expr->sup_coeff[0] = new_lambda_sup;
         term_expr->dim[0] = sub_expr->dim[0]; // Dimension of the pre-activation neuron
         term_expr->inf_cst = 0.0;
         term_expr->sup_cst = 0.0;


        // Add this term_expr to the overall result 'res'
        if (first_res_term) {
            res = term_expr;
            first_res_term = false;
        } else {
            add_expr(pr, res, term_expr);
            free_expr(term_expr);
        }
    } // End loop over terms

    // If no terms were processed, result is just the constant
    if (first_res_term) {
        res = current_const;
    } else {
        // Add the accumulated constant term to the result expression
        add_cst_expr(pr, res, current_const);
        free_expr(current_const);
    }

    // The resulting expression 'res' is now defined over the predecessors of the activation layer
    return res;
}


expr_t * lexpr_replace_bounds_activation(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons){
    return expr_replace_bounds_activation(pr, expr, neurons, true);
}

expr_t * uexpr_replace_bounds_activation(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons){
    return expr_replace_bounds_activation(pr, expr, neurons, false);
}


expr_t * lexpr_replace_bounds(fppoly_internal_t * pr, expr_t * expr, neuron_t ** neurons, bool is_activation){
    if(is_activation){
        return lexpr_replace_bounds_activation(pr, expr, neurons);
    }
    else{
        return lexpr_replace_bounds_affine(pr, expr, neurons);
    }
}

expr_t * uexpr_replace_bounds(fppoly_internal_t * pr, expr_t * expr, neuron_t ** neurons, bool is_activation){
    if(is_activation){
        return uexpr_replace_bounds_activation(pr, expr, neurons);
    }
    else{
        return uexpr_replace_bounds_affine(pr, expr, neurons);
    }
}

// Converts expr_t (internal representation) to elina_linexpr0_t (ELINA representation)
// Assumes expr_t is defined over the primary input variables.
elina_linexpr0_t *elina_linexpr0_from_expr(expr_t *expr){
    if (!expr) return NULL;

    size_t size = expr->size;
    size_t i, k;
    // Determine max dimension needed for dense allocation
    size_t max_dim = 0;
    if (expr->type == DENSE) {
        max_dim = size;
    } else { // SPARSE
        for (i = 0; i < size; i++) {
            if (expr->dim[i] >= max_dim) { // Use >= to handle dim[i] == max_dim case correctly
                max_dim = expr->dim[i] + 1; // Need space up to max_dim index
            }
        }
    }

    elina_linexpr0_t * linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_DENSE, max_dim);
    if (!linexpr0) return NULL;

    // ELINA uses scalar/interval coefficients. We store interval [-inf, sup].
    // When converting back, we need to decide how to represent this.
    // OPTION 1: Use interval coefficients in ELINA.
    // OPTION 2: Use scalar coefficients (e.g., midpoint? upper bound?) - Depends on ELINA consumer.
    // The original code used scalar `expr->sup_coeff[i]`. Let's stick to that for now,
    // assuming the consumer wants the upper bound expression. This might be incorrect
    // depending on how ELINA uses these expressions later.

    // Initialize all dense coefficients to zero
    // for (k = 0; k < max_dim; k++) {
    //     elina_linexpr0_set_coeff_scalar_double(linexpr0, k, 0.0);
    // }

    // Set coefficients based on expr_t
    for(i = 0; i < size; i++){
        if(expr->type == SPARSE){
            k = expr->dim[i];
        } else { // DENSE
            k = i;
        }
        if (k < max_dim) { // Bounds check
             // Using sup_coeff as per original code. This represents the upper bound coeff.
             // This loses the lower bound information. Consider using interval coefficients if needed.
            elina_linexpr0_set_coeff_scalar_double(linexpr0, k, expr->sup_coeff[i]);
        }
    }

    // Set constant term. Using sup_cst as per original code.
    elina_linexpr0_set_cst_scalar_double(linexpr0, expr->sup_cst);

    return linexpr0;
}



// here

/* ====================================================================== */
/* Compatibility Wrappers (Calling Sound Base Functions)                  */
/* ====================================================================== */

// Wrapper for sound interval addition (coefficients)
void elina_double_interval_add_expr_coeff(fppoly_internal_t *pr, double * res_inf, double *res_sup,
                                        double inf1, double sup1, double inf2, double sup2) {
    // Directly call the sound addition function
    elina_double_interval_add(res_inf, res_sup, inf1, sup1, inf2, sup2);
    // NO unsound error terms added here.
}

// Wrapper for sound interval addition (constants)
void elina_double_interval_add_cst_coeff(fppoly_internal_t *pr, double * res_inf, double *res_sup,
                                       double inf1, double sup1, double inf2, double sup2) {
    // Directly call the sound addition function
    elina_double_interval_add(res_inf, res_sup, inf1, sup1, inf2, sup2);
    // NO unsound error terms added here.
    // NOTE: The original added pr->min_denormal here. If that was essential for some
    // specific constant handling logic elsewhere (beyond basic FP soundness),
    // you might need to reconsider, but generally, base sound arithmetic shouldn't need it.
    // For now, removing it aligns with relying solely on fesetround for soundness.
}

// Wrapper for sound interval multiplication (coefficients)
void elina_double_interval_mul_expr_coeff(fppoly_internal_t *pr, double * res_inf, double *res_sup,
                                        double inf1, double sup1, double inf2, double sup2) {
    // Directly call the sound multiplication function (assumed defined elsewhere)
    elina_double_interval_mul(res_inf, res_sup, inf1, sup1, inf2, sup2);
    // NO unsound error terms added here.
}

// Wrapper for sound interval multiplication (constants)
void elina_double_interval_mul_cst_coeff(fppoly_internal_t *pr, double * res_inf, double *res_sup,
                                       double inf1, double sup1, double inf2, double sup2) {
    // Directly call the sound multiplication function (assumed defined elsewhere)
    elina_double_interval_mul(res_inf, res_sup, inf1, sup1, inf2, sup2);
    // NO unsound error terms added here.
    // NOTE: See comment in elina_double_interval_add_cst_coeff regarding pr->min_denormal.
}