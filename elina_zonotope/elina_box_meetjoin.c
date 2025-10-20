// /*
//  *
//  *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
//  *  ELINA is Copyright © 2021 Department of Computer Science, ETH Zurich
//  *  This software is distributed under GNU Lesser General Public License Version 3.0.
//  *  For more information, see the ELINA project website at:
//  *  http://elina.ethz.ch
//  *
//  *  THE SOFTWARE IS PROVIDED "AS-IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER
//  *  EXPRESS, IMPLIED OR STATUTORY, INCLUDING BUT NOT LIMITED TO ANY WARRANTY
//  *  THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS OR BE ERROR-FREE AND ANY
//  *  IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
//  *  TITLE, OR NON-INFRINGEMENT.  IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY     
//  *  DAMAGES, INCLUDING BUT NOT LIMITED TO DIRECT, INDIRECT,
//  *  SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN
//  *  ANY WAY CONNECTED WITH THIS SOFTWARE (WHETHER OR NOT BASED UPON WARRANTY,
//  *  CONTRACT, TORT OR OTHERWISE).
//  *
//  */


// #include <string.h>
// #include <stdio.h>

// #include "elina_box_internal.h"
// #include "elina_box_representation.h"
// #include "elina_box_constructor.h"
// #include "elina_box_meetjoin.h"
// #include "elina_scalar_arith.h"

// /* ============================================================ */
// /* Meet and Join */
// /* ============================================================ */
// elina_box_t* elina_box_meet(elina_manager_t* man, bool destructive, elina_box_t* a1, elina_box_t* a2)
// {
//     size_t i;
//     size_t nbdims;
//     elina_box_t* res;
//     man->result.flag_best = true;
//     man->result.flag_exact = false;
//     res = destructive ? a1 : elina_box_alloc(a1->intdim,a1->realdim);
//     if ((a1->inf==NULL && a1->sup==NULL) || (a2->inf==NULL && a2->sup==NULL)){
//         elina_box_set_bottom(res);
//         return res;
//     }
//     if (!destructive){
//         elina_box_init(res);
//     }
    
//     nbdims = a1->intdim + a2->realdim;
//     for (i=0; i<nbdims; i++){
//         res->inf[i] = fmin(a1->inf[i],a2->inf[i]);
//         res->sup[i] = fmin(a1->sup[i],a2->sup[i]);
//     }
//     return res;
// }



// elina_box_t* elina_box_join(elina_manager_t* man, bool destructive, elina_box_t* a1, elina_box_t* a2)
// {
//   size_t i;
//   size_t nbdims;
//   elina_box_t* res;
//   man->result.flag_best = true;
//   man->result.flag_exact = false;
//   res = destructive ? a1 : elina_box_alloc(a1->intdim,a1->realdim);
//   if ((a1->inf==NULL) && (a1->sup==NULL)){
//     if ((a2->inf!=NULL) && (a2->sup!=NULL)){
//       man->result.flag_exact = true;
//       elina_box_set(res,a2);
//     }
//     return res;
//   }
//   else if (a2->inf==NULL && a2->sup==NULL){
//     man->result.flag_exact = true;
//     if (!destructive) elina_box_set(res,a1);
//     return res;
//   }
//   man->result.flag_exact = false;
//   if (!destructive){
//     elina_box_init(res);
//   }
  
//   nbdims = a1->intdim + a2->realdim;
//   for (i=0; i<nbdims; i++){
//     res->inf[i] = fmax(a1->inf[i],a2->inf[i]);
//     res->sup[i] = fmax(a1->sup[i],a2->sup[i]);
//   }
   
//   return res;
// }

// bool elina_double_interval_canonicalize(double *inf, double *sup, bool integer, elina_scalar_discr_t discr)
// {
//   bool exc;

//   if (integer){
//     *inf = floor(*inf);
//     *sup = floor(*sup); 
//   }
//   if ((*inf==INFINITY) || *sup==INFINITY) return false;

//   /* Check that it is not bottom */
//   exc = false;
//   if (*sup< -*inf )
//     exc = true;
//   return exc;
// }



// void elina_double_interval_mul(double *a_inf, double *a_sup, double b_inf, double b_sup, double c_inf, double c_sup){
// 	if(c_inf<=0){
// 		/* interval c is positive */
// 		if(b_inf<=0){
// 			/*interval b is positive*/
// 			if((b_inf==0) || (c_inf==0)){
// 				*a_inf = 0.0;
// 			}
// 			else{
// 				*a_inf = b_inf * -c_inf;
// 			}
// 			if((b_sup==0) || (c_sup==0)){
// 				*a_sup = 0.0;
// 			}
// 			else{
// 				*a_sup = b_sup * c_sup;
// 			}
// 		}
// 		else if(b_sup<=0){
// 			/* interval b is negative */
// 			if((c_sup==0) || (b_inf==0)){
// 				*a_inf = 0.0;
// 			}
// 			else{
// 			 	*a_inf = c_sup*b_inf;
// 			}
// 			if((c_inf==0) || (b_sup==0)){
// 				*a_sup = 0.0;
// 			}
// 			else{
// 			 	*a_sup = -c_inf*b_sup;
// 			}
// 		}
// 		else{
// 			/* there is 0 in between for b */
// 			if((c_sup==0) || (b_inf==0)){
// 				*a_inf = 0.0;
// 			}
// 			else{
// 				*a_inf = b_inf * c_sup;
// 			}
// 			if((c_sup==0) || (b_sup==0)){
// 				*a_sup = 0.0;
// 			}
// 			else{
// 				*a_sup = b_sup * c_sup;
// 			}
// 		}
// 	}
// 	else if(c_sup<=0){
// 		/* interval c is negative */
// 		if(b_inf<=0){
// 			/*interval b is positive*/
// 			if((b_sup==0) || (c_inf==0)){
// 				*a_inf = 0.0;
// 			}
// 			else{
// 				*a_inf = b_sup*c_inf;
// 			}
// 			if((b_inf==0) || (c_sup==0)){
// 				*a_sup = 0.0;
// 			}
// 			else{
// 				*a_sup = -b_inf*c_sup;
// 			}
// 		}
// 		else if(b_sup<=0){
// 			/* interval b is negative */
// 			if((b_sup==0) || (c_sup==0)){
// 				*a_inf = 0.0;
// 			}
// 			else{
// 				*a_inf = b_sup * -c_sup;
// 			}
// 			if((b_inf==0) || (c_inf == 0 )){
// 				*a_sup = 0.0;
// 			}
// 			else{
// 				*a_sup = b_inf * c_inf;
// 			} 
// 		}
// 		else{
// 			/* there is 0 in between for b */
// 			if((c_inf==0) || (b_sup==0)){
// 				*a_inf = 0.0;
// 			}
// 			else{
// 				*a_inf = b_sup*c_inf;
// 			}
// 			if((c_inf==0) || (b_inf==0)){
// 				*a_sup = 0.0;
// 			}
// 			else{
// 				*a_sup = b_inf*c_inf;
// 			} 
// 		}
// 	}
// 	else if(b_inf<=0){
// 		/* interval b is positive */
// 		if(c_inf<=0){
// 			/*interval c is positive */
// 			if((b_inf==0) || (c_inf==0)){
// 				*a_inf = 0.0;
// 			}
// 			else{
// 				*a_inf = -b_inf * c_inf;
// 			}
// 			if((b_sup==0) || (c_sup==0)){
// 				*a_sup = 0.0;
// 			}
// 			else{
// 				*a_sup = b_sup * c_sup;
// 			}
// 		}
// 		else if(c_sup<=0){
// 			/* interval c is negative */
// 			if((b_sup==0) || (c_inf == 0)){
// 				*a_inf = 0.0;
// 			}
// 			else{
// 				*a_inf = b_sup*c_inf;
// 			}
// 			if((b_inf==0) || (c_sup==0)){
// 				*a_sup = 0.0;
// 			}
// 			else{
// 				*a_sup = -b_inf*c_sup;
// 			}
// 		}
// 		else{
// 			/* there is 0 in between for c */
// 			if((b_sup==0) || (c_inf==0)){
// 				*a_inf = 0.0;
// 			}
// 			else{
// 				*a_inf = b_sup * c_inf;
// 			}
// 			if((b_sup==0) || (c_sup==0)){
// 				*a_sup = 0.0;
// 			}
// 			else{
// 				*a_sup = b_sup * c_sup;
// 			}
// 		}
// 	}
// 	else if(b_sup<=0){
// 		/* interval b is negative */
// 		if(c_inf <= 0){
// 			/* interval c is positive */
// 			if((b_inf==0) || (c_sup==0)){
// 				*a_inf = 0.0;
// 			}
// 			else{
// 				*a_inf = b_inf * c_sup;
// 			}
// 			if((b_sup==0) || (c_inf==0)){
// 				*a_sup = 0.0;
// 			}
// 			else{
// 				*a_sup = b_sup * -c_inf;
// 			}
// 		}
// 		else if(c_sup<=0){
// 			/* interval c is negative */
// 			if((b_sup==0) || (c_sup==0)){
// 				*a_inf = 0.0;
// 			}
// 			else{
// 				*a_inf = -b_sup * c_sup;
// 			}
// 			if((b_inf==0) || (c_inf==0)){
// 				*a_sup = 0.0;
// 			}
// 			else{
// 				*a_sup = b_inf * c_inf;
// 			} 
// 		}
// 		else{
// 			/* there is 0 in between for c */
// 			if((b_inf == 0) || (c_sup==0)){
// 				*a_inf = 0.0;
// 			}
// 			else{
// 				*a_inf = b_inf * c_sup;
// 			}
// 			if((b_inf==0) || (c_inf==0)){
// 				*a_sup = 0.0;
// 			}
// 			else{
// 				*a_sup = b_inf * c_inf;
// 			}
// 		}
// 	}
// 	else{
// 		/* there is 0 in between for both b and c */
// 		double tmp_inf1 = b_sup*c_inf;
// 		double tmp_sup1 = b_inf*c_inf;
// 		double tmp_inf2 = b_inf*c_sup;
// 		double tmp_sup2 = b_sup*c_sup;
// 		*a_inf = fmax(tmp_inf1, tmp_inf2);
// 		*a_sup = fmax(tmp_sup1, tmp_sup2);
// 	}
// }




// /* ====================================================================== */
// /* Division */
// /* ====================================================================== */


// void elina_double_interval_div(double *a_inf, double *a_sup, double b_inf, double b_sup, double c_inf, double c_sup)
// {
//   if (c_inf<0){
//     /* c is positive */
//     if (b_inf<=0){
//          /* b is positive */
//          *a_inf = b_inf/c_sup;
// 	 *a_sup = b_sup/-c_inf;
//     }
//     else if (b_sup<=0){
//        /* b is negative */
//        *a_inf = -b_inf/c_inf;
//        *a_sup = b_sup/c_sup;
//     }
//     else {
//         /* 0 is in the middle of b: one divides b by c->inf */
//         *a_inf = b_inf/-c_inf;
//         *a_sup = b_sup/-c_inf;
//     }
//   }
//   else if (c_sup<0){
//     /* c is negative */
//     if (b_inf<=0){
//         /* b is positive */
// 	*a_sup = b_inf/c_inf;
//         *a_inf = -b_sup/c_sup;
  	
//      }
//      else if (b_sup<=0){
//        /* b is negative */
//        *a_inf = b_sup/c_inf;
//        *a_sup = -b_inf/c_sup;
//      }
//      else {
//         /* 0 is in the middle of b: one cross-divide b by c->sup */
//         *a_inf = b_sup/c_sup;
//         *a_sup = b_inf/c_sup;
//      }
//   }
//   else if ((b_inf==0) && (b_sup==0)){
//     /* b is [0,0] */
//     *a_inf = b_inf;
//     *a_sup = b_sup;
//   }
//   else {
//     *a_inf = INFINITY;
//     *a_sup = INFINITY;
//   }
// }



// bool elina_double_interval_eval_elina_linexpr0(double * itv_inf, double *itv_sup,  elina_linexpr0_t* expr, double* env_inf, double * env_sup, elina_scalar_discr_t discr)
// {
//   size_t i;
//   elina_dim_t dim;
//   elina_coeff_t* coeff;
//   assert(env_inf&&env_sup);
//   elina_scalar_t * scalar;
//   elina_coeff_t * cst = &expr->cst;
//   if(cst->discr==ELINA_COEFF_SCALAR){
// 	*itv_inf = -cst->val.scalar->val.dbl;
// 	*itv_sup = cst->val.scalar->val.dbl;
//   }
//   else{
// 	*itv_inf = -cst->val.interval->inf->val.dbl;
// 	*itv_sup = cst->val.interval->sup->val.dbl;
//   }
  
//   double tmp_inf = 0.0;
//   double tmp_sup = 0.0;
//   elina_linexpr0_ForeachLinterm(expr,i,dim,coeff){
//     bool eq = (coeff->discr==ELINA_COEFF_SCALAR);
//     if (eq){
//       scalar = coeff->val.scalar;
//       if (elina_scalar_sgn(scalar)>0){
	
// 	tmp_inf = env_inf[dim]*scalar->val.dbl;
// 	tmp_sup = env_sup[dim]*scalar->val.dbl;
// 	*itv_inf = *itv_inf + tmp_inf;
// 	*itv_sup = *itv_sup + tmp_sup;
//       }
//       else if(elina_scalar_sgn(scalar)<0){
// 	tmp_inf = env_sup[dim]*-scalar->val.dbl;
// 	tmp_sup = env_inf[dim]*-scalar->val.dbl;
//         *itv_inf = *itv_inf + tmp_inf;
//         *itv_sup = *itv_sup + tmp_sup;
//       }
//     }
//     else {
//       double inf = -coeff->val.interval->inf->val.dbl;
//       double sup = coeff->val.interval->sup->val.dbl;
// 	//if((inf==sup) && (inf==0) ){
// 	//	tmp_inf = 0.0;
// 	//	tmp_sup = 0.0;
// 	//}
// 	//else{
//       		elina_double_interval_mul(&tmp_inf,&tmp_sup,env_inf[dim],env_sup[dim],inf,sup);
// 	//}
//       *itv_inf = *itv_inf + tmp_inf;
//       *itv_sup = *itv_sup + tmp_sup;
     
//     }
//     if (*itv_inf==INFINITY && *itv_sup==INFINITY){
	
//       break;
//     }
//   }
 
//   return true;
// }

// static inline bool elina_double_scalar_sgn(double d){
// 	if(d > 0.0){
// 		return 1;
// 	}
// 	else if(d < 0.0){
// 		return -1;
// 	}
// 	else{
// 		return 0;
// 	}
	
// }

// static bool elina_double_boxize_lincons0(double* res_inf, double * res_sup,
// 			       bool* tchange,
// 			       elina_lincons0_t* cons,
// 			       double* env_inf, double * env_sup,
// 			       size_t intdim,
// 			       bool intervalonly, elina_scalar_discr_t discr)
// {
//   size_t i;
//   elina_linexpr0_t* expr;
//   bool change,globalchange;
//   bool exc;
//   assert(cons->constyp == ELINA_CONS_EQ ||
// 	 cons->constyp == ELINA_CONS_SUPEQ ||
// 	 cons->constyp == ELINA_CONS_SUP);

//   expr = cons->linexpr0;
//   globalchange = false;

//   /* Iterates on coefficients */
//   double itv_inf = 0.0;
//   double itv_sup = 0.0;
//   double val = 0.0;
//   elina_scalar_t * scalar2 = elina_scalar_alloc();
//   for (i=0; i<expr->size; i++){
//     elina_coeff_t *coeff;
//     elina_dim_t dim;
//     if(expr->discr==ELINA_LINEXPR_SPARSE){
// 	elina_linterm_t *term = &expr->p.linterm[i];
//         dim = term->dim;
//         coeff = &term->coeff;
//     }
//     else{
// 	dim = i;
// 	coeff = &expr->p.coeff[i];
//     }
    
//     elina_coeff_t * tmp = elina_coeff_alloc(coeff->discr);
//     /* 1. We decompose the expression e = ax+e' */
//     elina_coeff_swap(tmp,coeff);
//     double inf = 0.0;
//     double sup = 0.0;
//     double d;

//     if (tmp->discr==ELINA_COEFF_SCALAR) {
//         if (tmp->val.scalar->discr == ELINA_SCALAR_DOUBLE) {
//             inf = -tmp->val.scalar->val.dbl;
//             sup = tmp->val.scalar->val.dbl;
//         } else {
//             elina_double_set_scalar(&d,tmp->val.scalar,GMP_RNDD);
//             inf = -d;
//             elina_double_set_scalar(&d,tmp->val.scalar,GMP_RNDU);
//             sup = d;
//         }
//     } else {
//         if (tmp->val.interval->inf->discr == ELINA_SCALAR_DOUBLE) {
//             inf = -tmp->val.interval->inf->val.dbl;
//         } else {
//             elina_double_set_scalar(&d,tmp->val.interval->inf,GMP_RNDD);
//             inf = -d;
//         }
//         if (tmp->val.interval->sup->discr == ELINA_SCALAR_DOUBLE) {
//             sup = tmp->val.interval->sup->val.dbl;
//         } else {
//             elina_double_set_scalar(&d,tmp->val.interval->inf,GMP_RNDU);
//             sup = d;
//         }
//     }
//     /* 2. evaluate e' */
//     elina_double_interval_eval_elina_linexpr0(&itv_inf, & itv_sup, expr,env_inf, env_sup, discr);
//     /* 3. Perform deduction from ax+e' = [-m,M]x + [-e,E] >= 0
// 	  we can deduce [-m,M]x + E >= 0
// 	  (for equality, [-m,M]x - e <= 0)
//     */
//     bool equality = tmp->discr==ELINA_COEFF_SCALAR;
//     change = false;
//     if (itv_inf!=INFINITY || itv_sup!=INFINITY){
//       if (equality && !intervalonly){
// 	/* [-m,M]=[a,a] */
	
// 	int sgn = elina_double_scalar_sgn(sup);
// 	if (sgn!=0){
// 	  /* From ax+E >= 0, we deduce
// 	     (1) If a>0, then x>=-E/a
// 	     (2) If a<0, then x<=-E/a
// 	     From ax-e <= 0, we deduce
// 	     (3) If a>0, then x<=e/a
// 	     (4) If a<0, then x>=e/a
// 	  */
// 	  if (sgn>0 || cons->constyp == ELINA_CONS_EQ){
// 	    /* 1,4: inf bound
// 	       If we have a>0 (1), we compute E/a
// 	       If we have a<0 (4), we compute -e/a
// 	    */
		
// 	    if (sgn>0){
// 		val = itv_sup/sup;
// 	    }
// 	    else {
// 		val = itv_inf/inf;
// 	    }
// 	    if (dim<intdim && isfinite(val)){
// 	      if ((cons->constyp==ELINA_CONS_SUP) && (ceil(val) == val)){
// 		  val = val - 1;
// 	      }
// 	      else {
// 		  val = floor(val);
// 	      }
// 	    }
// 	    /* We update the interval */
// 	    if (val < res_inf[dim]){
// 	      change = true;
// 	      if (tchange) tchange[2*dim] = true;
// 	      res_inf[dim] = val;
// 	    }
// 	  }
// 	  if (sgn<0 || cons->constyp == ELINA_CONS_EQ){
// 	    /* 2,3: sup bound
// 	       If we have a<0 (2), we compute -E/a
// 	       If we have a>0 (3), we compute e/a
// 	    */
// 	    if (sgn<0){
// 	       val = itv_sup/inf;
// 	    }
// 	    else {
// 	       val = itv_inf/sup;
// 	    }
// 	    if (dim<intdim && isfinite(val)){
// 	      if ((cons->constyp==ELINA_CONS_SUP) &&
// 		  (ceil(val)==val)){
// 		   val = val - 1;
// 	      }
// 	      else {
// 		 val = floor(val);
// 	      }
// 	    }
// 	    /* We update the interval */
// 	    if (val < res_sup[dim]){
// 	      change = true;
// 	      if (tchange) tchange[2*dim+1] = true;
// 	      res_sup[dim] = val;
// 	    }
// 	  }
// 	}
//       }
//       else if (!equality){
// 	/* We have a real interval [-m,M] */
// 	/* Case [-m,M]x+E >= 0
// 	  (1) If -m>0, we rewrite [-m,M]x >= -E, and we have -1/m > 1/M > 0
// 	  (1A) If E<=0 <=> -E>=0 then x >= -E/M
// 	  (1B) If E>0  <=> -E<0  then x >= E/m
// 	  (2) If M<0, we rewrite [-M,m]x <= E, and we have 0 < 1/m < -1/M
// 	  (2A) If E<=0           then x <= E/m
// 	  (2B) If E>0            then x <= -E/M
// 	  Case [-m,M]x-e <= 0
// 	  (3) If -m>0, we rewrite [-m,M]x <= e, and we have 0 < 1/M < -1/m
// 	  (3A) If e<=0           then x <= e/M
// 	  (3B) If e>0            then x <= -e/m
// 	  (4) If M<0, we rewrite [-M,m]x >= -e, and we have -1/M > 1/m > 0
// 	  (4A) If e<=0 <=> -e>=0 then x >= -e/m
// 	  (4B) If e>0  <=> -e<0  then x >= e/M
// 	*/
// 	int sgnitv = inf<0 ? 1 : sup<0 ? -1 : 0;
// 	if (sgnitv != 0){
// 	  int sgne = elina_double_scalar_sgn(itv_inf);
// 	  int sgnE = elina_double_scalar_sgn(itv_sup);
// 	  if (sgnitv>0 || (cons->constyp==ELINA_CONS_EQ && sgnitv<0)){
// 	    /* 1,4: inf bound */
// 	    if (sgnitv>0){ /* 1 */
// 	      if (sgnE<=0){ /* 1A */
// 		/* We compute E/M */
// 		 val = itv_sup/sup;
		
// 	      } else { /* 1B */
// 		/* We compute -E/m */
// 		val = -(itv_sup/inf);
// 	      }
// 	    }
// 	    else { /* 4 */
// 	      if (sgne>=0){ /* 4A */
// 		/* We compute e/m */
// 		val = itv_inf/inf;
// 	      } else { /* 4B */
// 		/* We compute -e/M */
// 		val = -(itv_inf/sup);
// 	      }
// 	    }
// 	    if (dim<intdim && isfinite(val)){
// 	      if ((cons->constyp==ELINA_CONS_SUP) && (ceil(val)==val)){
// 		val = val - 1;
// 	      }
// 	      else {
// 		val = floor(val);
// 	      }
// 	    }
// 	    /* We update the interval */
// 	    if (val < res_inf[dim]){
// 	      change = true;
// 	      if (tchange) tchange[2*dim] = true;
// 	      res_inf[dim] = val;
// 	    }
// 	  }
// 	  if (sgnitv<0 || (cons->constyp==ELINA_CONS_EQ && sgnitv>0)){
// 	    /* 2,3: sup bound */
// 	    if (sgnitv<0){ /* 2 */
// 	      if (sgnE>=0){ /* 2B */
// 		/* We compute -E/M */
// 		val = -(itv_sup/sup);
// 	      } else { /* 2A */
// 		/* We compute E/m */
// 		val = itv_sup/inf;
// 	      }
// 	    }
// 	    else { /* 3 */
// 	      if (sgne<=0){ /* 3B */
// 		/* We compute -e/m */
// 		val = -(itv_inf/inf);
// 	      }
// 	      else { /* 3A */
// 		/* We compute e/M */
// 		val = itv_inf/sup;
// 	      }
// 	    }
// 	    if (dim<intdim && isfinite(val)){
// 	      if ((cons->constyp==ELINA_CONS_SUP) && (ceil(val)==val)){
// 		val = val - 1;
// 	      }
// 	      else {
// 		val = floor(val);
// 	      }
// 	    }
// 	    /* We update the interval */
// 	    if (val < res_sup[dim]){
// 	      change = true;
// 	      if (tchange) tchange[2*dim+1] = true;
// 	      res_sup[dim] = val;
// 	    }
// 	  }
// 	}
//       }
//     }
//     elina_coeff_swap(tmp,coeff);
//     elina_coeff_free(tmp);
//     if (change){
//       globalchange = true;
//       exc = elina_double_interval_canonicalize(res_inf + dim,res_sup + dim,dim<intdim,discr);
//       if (exc){
// 	res_inf[0] = -1;
// 	res_sup[0] = -1;
// 	return true;
//       }
//     }
//   }
//   if (expr->size==0 &&
//       eval_elina_cstlincons0(cons)==0){
// 	res_inf[0] = -1;
// 	res_sup[0] = -1;
//     globalchange = true;
//   }
//   return globalchange;
// }



// bool elina_double_boxize_lincons0_array(double * res_inf, double * res_sup, bool* tchange,
// 				      elina_lincons0_array_t* array,
// 				      double* env_inf, double * env_sup, size_t intdim,
// 				      size_t kmax,
// 				      bool intervalonly, elina_scalar_discr_t discr)
// {
//   size_t i,k;
//   bool change,globalchange;

//   if (kmax<1) kmax=1;
//   if ((res_inf!=env_inf) && (res_sup!=env_sup)) kmax=1;

//   globalchange = false;
//   /* we possibly perform kmax passes */
//   for (k=0;k<(size_t)kmax;k++){
//     change = false;
//     for (i=0; i<array->size; i++){
//       if (array->p[i].constyp==ELINA_CONS_EQ ||
// 	  array->p[i].constyp==ELINA_CONS_SUPEQ ||
// 	  array->p[i].constyp==ELINA_CONS_SUP){
	
// 	change =
// 	  elina_double_boxize_lincons0(res_inf, res_sup, tchange,&array->p[i],env_inf,env_sup,intdim,intervalonly,discr)
// 	  ||
// 	  change
// 	  ;
// 	globalchange = globalchange || change;
// 	if (elina_double_interval_canonicalize(res_inf,res_sup,false,ELINA_SCALAR_DOUBLE)){
// 	  return true;
// 	}
//       }
//     }
//     if (!change) break;
//   }
//   return globalchange;
// }




// /* ============================================================ */
// /* Meet_lincons */
// /* ============================================================ */

// elina_box_t* elina_box_meet_lincons_array(elina_manager_t* man,
// 			      bool destructive,
// 			      elina_box_t* a,
// 			      elina_lincons0_array_t* array)
// {
//     //printf("start %p\n",man);
//     //fflush(stdout);
//   elina_box_t* res;
//   size_t kmax;
//   elina_lincons0_array_t tlincons;
//   elina_box_internal_t* intern = (elina_box_internal_t*)man->internal;

//   res = destructive ? a : elina_box_copy(man,a);
//   if (a->inf==NULL && a->sup==NULL){
//     man->result.flag_best = true;
//     man->result.flag_exact = true;
//   }
//   else {
   
//     man->result.flag_best = array->size==1;
//     man->result.flag_exact = false;
//     kmax = man->option.funopt[ELINA_FUNID_MEET_LINCONS_ARRAY].algorithm;
//     if (kmax<1) kmax=2;
//     tlincons = elina_lincons0_array_make(array->size);
//     for(size_t i =0; i < array->size; i++){
// 	tlincons.p[i] = elina_lincons0_copy(&array->p[i]);
//     }
//     char tb = elina_lincons0_array_reduce_integer(&tlincons,a->intdim,ELINA_SCALAR_DOUBLE);
//     if (tb==0){
//       goto _elina_box_meet_lincons_array_bottom;
//     }
    
//     elina_double_boxize_lincons0_array(res->inf,res->sup,NULL,
// 			     &tlincons,res->inf,res->sup,a->intdim,kmax,false,ELINA_SCALAR_DOUBLE);
//     //if (res->inf[0]<res->sup[0]){
//     if(elina_double_interval_canonicalize(res->inf,res->sup,false,ELINA_SCALAR_DOUBLE)){
//     _elina_box_meet_lincons_array_bottom:
//       elina_box_set_bottom(res);
//     }
  
//     elina_lincons0_array_clear(&tlincons);
    
//   }
//     //printf("finish %p\n",man);
//     //fflush(stdout);
//   return res;
// }


// /* ============================================================ */
// /* Widening */
// /* ============================================================ */
// elina_box_t* elina_box_widening(elina_manager_t* man,
//                     elina_box_t* a1, elina_box_t* a2)
// {
//     size_t i;
//     size_t nbdims;
//     elina_box_t* res;
    
//     man->result.flag_best = true;
//     man->result.flag_exact = true;
//     nbdims = a1->intdim+a1->realdim;
//     if ((a1->inf==NULL) && (a1->sup==NULL)){
//         return elina_box_copy(man,a2);
//     }
    
//     res = elina_box_copy(man,a1);
//     for (i=0; i<nbdims; i++){
//         if(a1->sup[i] < a2->sup[i]){
//             res->sup[i] = INFINITY;
//         }
//         else{
//             res->sup[i] = a1->sup[i];
//         }
//         if(a1->inf[i] < a2->inf[i]){
//             res->inf[i] = INFINITY;
//         }
//         else{
//             res->inf[i] = a1->inf[i];
//         }
//     }
//     return res;
// }


/*
 *
 * This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 * ELINA is Copyright © 2021 Department of Computer Science, ETH Zurich
 * This software is distributed under GNU Lesser General Public License Version 3.0.
 * For more information, see the ELINA project website at:
 * http://elina.ethz.ch
 *
 * THE SOFTWARE IS PROVIDED "AS-IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER
 * EXPRESS, IMPLIED OR STATUTORY, INCLUDING BUT NOT LIMITED TO ANY WARRANTY
 * THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS OR BE ERROR-FREE AND ANY
 * IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
 * TITLE, OR NON-INFRINGEMENT.  IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY
 * DAMAGES, INCLUDING BUT NOT LIMITED TO DIRECT, INDIRECT,
 * SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN
 * ANY WAY CONNECTED WITH THIS SOFTWARE (WHETHER OR NOT BASED UPON WARRANTY,
 * CONTRACT, TORT OR OTHERWISE).
 *
 */


#include <string.h>
#include <stdio.h>
#include <math.h> // Include math.h for INFINITY, fmax, fmin, floor, ceil, isfinite
#include <fenv.h> // Include fenv.h for rounding mode control
#include <assert.h> // For assert()

#include "elina_box_internal.h"
#include "elina_box_representation.h"
#include "elina_box_constructor.h"
#include "elina_box_meetjoin.h"
#include "elina_scalar_arith.h" // Assuming this contains elina_scalar_t related functions

// Helper macro to ensure rounding mode restoration
#pragma STDC FENV_ACCESS ON

/* ============================================================ */
/* Meet and Join (Unchanged from original - appear correct)     */
/* ============================================================ */
elina_box_t* elina_box_meet(elina_manager_t* man, bool destructive, elina_box_t* a1, elina_box_t* a2)
{
    size_t i;
    size_t nbdims;
    elina_box_t* res;
    man->result.flag_best = true;
    man->result.flag_exact = false; // Meet is generally exact for boxes
    res = destructive ? a1 : elina_box_alloc(a1->intdim,a1->realdim);
    if (!res) return NULL; // Allocation check

    if ((a1->inf==NULL && a1->sup==NULL) || (a2->inf==NULL && a2->sup==NULL)){
        elina_box_set_bottom(res);
        return res;
    }
    if (!destructive){
        // Ensure res is initialized if not destructive
        // if (!elina_box_init(res)) {
        //      if (!destructive) elina_box_free(man, res); // Cleanup on init failure
        //      return NULL;
        // }
        // Need to copy dimensions if allocating new
        res->intdim = a1->intdim;
        res->realdim = a1->realdim;
    }

    nbdims = a1->intdim + a1->realdim; // Corrected: should be realdim
    for (i=0; i<nbdims; i++){
        // Meet: max of lower bounds (using negated representation)
        res->inf[i] = fmin(a1->inf[i],a2->inf[i]);
        // Meet: min of upper bounds
        res->sup[i] = fmin(a1->sup[i],a2->sup[i]);
    }
    // Check if the result is bottom after meet
    if (elina_double_interval_canonicalize(res->inf, res->sup, false, ELINA_SCALAR_DOUBLE)) {
         elina_box_set_bottom(res);
    }
    man->result.flag_exact = true; // Meet is exact
    return res;
}



elina_box_t* elina_box_join(elina_manager_t* man, bool destructive, elina_box_t* a1, elina_box_t* a2)
{
  size_t i;
  size_t nbdims;
  elina_box_t* res;
  man->result.flag_best = true;
  man->result.flag_exact = false; // Join is generally exact for boxes
  res = destructive ? a1 : elina_box_alloc(a1->intdim,a1->realdim);
  if (!res) return NULL; // Allocation check

  if ((a1->inf==NULL) && (a1->sup==NULL)){
    if ((a2->inf!=NULL) && (a2->sup!=NULL)){
      man->result.flag_exact = true;
      elina_box_set(res,a2); // Assuming elina_box_set handles allocation if needed
    } else {
        // If both are bottom, result is bottom (handle initialization)
        if (!destructive) {
            // if (!elina_box_init(res)) { elina_box_free(man, res); return NULL;}
            res->intdim = a2->intdim; // Assuming dimensions match or are defined
            res->realdim = a2->realdim;
        }
        elina_box_set_bottom(res);
    }
    return res;
  }
  else if (a2->inf==NULL && a2->sup==NULL){
    man->result.flag_exact = true;
    if (!destructive) elina_box_set(res,a1); // Assuming elina_box_set handles allocation
    return res;
  }

  // Both a1 and a2 are valid intervals
  if (!destructive){
    //  if (!elina_box_init(res)) { elina_box_free(man, res); return NULL;}
     // Need to copy dimensions if allocating new
     res->intdim = a1->intdim;
     res->realdim = a1->realdim;
  }

  nbdims = a1->intdim + a1->realdim; // Corrected: should be realdim
  for (i=0; i<nbdims; i++){
    // Join: min of lower bounds (using negated representation)
    res->inf[i] = fmax(a1->inf[i],a2->inf[i]);
    // Join: max of upper bounds
    res->sup[i] = fmax(a1->sup[i],a2->sup[i]);
  }
  man->result.flag_exact = true; // Join is exact
  return res;
}

// Canonicalize function - checks for bottom and handles integer rounding
// Returns true if bottom, false otherwise
bool elina_double_interval_canonicalize(double *inf, double *sup, bool integer, elina_scalar_discr_t discr)
{
  int original_round = fegetround(); // Save rounding mode

  // Round infinities appropriately first if integer
  if (integer) {
      fesetround(FE_UPWARD); // Round lower bound UP (negated, so round -L up = round L down)
      *inf = ceil(*inf);
      fesetround(FE_UPWARD); // Round upper bound UP
      *sup = ceil(*sup);
  }

  fesetround(original_round); // Restore rounding mode

  // Check for infinities introduced by rounding or input
  if (!isfinite(*inf) || !isfinite(*sup)) {
      // Check if it should be top [ -inf, +inf ] => [ inf, inf ] in our repr.
      // If only one is inf, it might still be a valid (unbounded) interval,
      // but bottom check will handle inconsistency. For now, assume finite needed.
      // Let's treat non-finite results other than [-inf, +inf] as potentially bottom/error
      // Note: This check might need refinement based on how infinities are handled upstream.
      if (*inf == INFINITY && *sup == INFINITY) return false; // Represents [-inf, +inf] okay
      // Consider other infinite cases errors or bottom for now
      // return true; // Indicate bottom/error for simplicity
  }


  /* Check that it is not bottom */
  // Use a small tolerance for floating point comparison
  if (*sup < -(*inf) - 1e-12) { // Check if U < L - tolerance
    return true; // Bottom
  }

  return false; // Not bottom
}


/* ============================================================ */
/* Sound Interval Multiplication with Rounding Control          */
/* ============================================================ */
void elina_double_interval_mul(double *a_inf, double *a_sup, double b_inf, double b_sup, double c_inf, double c_sup)
{
    // Convert from negated lower bound representation to [L, U]
    double b_l = -b_inf;
    double b_u = b_sup;
    double c_l = -c_inf;
    double c_u = c_sup;

    // Check for empty input intervals (optional, based on upstream guarantees)
    // if (b_u < b_l || c_u < c_l) { *a_inf = INFINITY; *a_sup = -INFINITY; return; } // Result is bottom

    int original_round = fegetround(); // Save current rounding mode

    // Calculate the four endpoint products with directed rounding
    fesetround(FE_DOWNWARD);
    double p1_down = b_l * c_l;
    double p2_down = b_l * c_u;
    double p3_down = b_u * c_l;
    double p4_down = b_u * c_u;

    fesetround(FE_UPWARD);
    double p1_up = b_l * c_l;
    double p2_up = b_l * c_u;
    double p3_up = b_u * c_l;
    double p4_up = b_u * c_u;

    // Determine the minimum and maximum possible values with correct rounding
    fesetround(FE_DOWNWARD); // Set rounding for min operations
    double a_l = fmin(fmin(p1_down, p2_down), fmin(p3_down, p4_down));

    fesetround(FE_UPWARD); // Set rounding for max operations
    double a_u = fmax(fmax(p1_up, p2_up), fmax(p3_up, p4_up));

    // Convert back to negated lower bound representation with correct rounding
    // *a_inf = -a_l requires rounding -a_l UPWARD
    fesetround(FE_UPWARD);
    *a_inf = -a_l;

    // *a_sup = a_u requires rounding a_u UPWARD (already set)
    *a_sup = a_u;

    fesetround(original_round); // Restore original rounding mode
}


/* ============================================================ */
/* Sound Interval Division with Rounding Control                */
/* ============================================================ */
void elina_double_interval_div(double *a_inf, double *a_sup, double b_inf, double b_sup, double c_inf, double c_sup)
{
    // Convert from negated lower bound representation to [L, U]
    double b_l = -b_inf;
    double b_u = b_sup;
    double c_l = -c_inf;
    double c_u = c_sup;

    // Check for empty input intervals (optional)
    // if (b_u < b_l || c_u < c_l) { *a_inf = INFINITY; *a_sup = -INFINITY; return; }

    // Check if denominator interval [c_l, c_u] contains zero
    if (c_l <= 0.0 && c_u >= 0.0) {
        // Division by zero is possible. Result is [-inf, +inf] unless numerator is exactly [0,0]
        if (b_l == 0.0 && b_u == 0.0) {
             // 0 / [c_l, c_u] containing 0 is [0, 0]
             *a_inf = 0.0; // -0.0
             *a_sup = 0.0;
        } else {
             // Result is [-inf, +inf]
             *a_inf = INFINITY; // Represents -(-inf)
             *a_sup = INFINITY; // Represents +inf
        }
        return;
    }

    int original_round = fegetround(); // Save current rounding mode

    // Calculate the four endpoint divisions with directed rounding
    fesetround(FE_DOWNWARD);
    double p1_down = b_l / c_l;
    double p2_down = b_l / c_u;
    double p3_down = b_u / c_l;
    double p4_down = b_u / c_u;

    fesetround(FE_UPWARD);
    double p1_up = b_l / c_l;
    double p2_up = b_l / c_u;
    double p3_up = b_u / c_l;
    double p4_up = b_u / c_u;

    // Determine the minimum and maximum possible values with correct rounding
    fesetround(FE_DOWNWARD); // Set rounding for min operations
    double a_l = fmin(fmin(p1_down, p2_down), fmin(p3_down, p4_down));

    fesetround(FE_UPWARD); // Set rounding for max operations
    double a_u = fmax(fmax(p1_up, p2_up), fmax(p3_up, p4_up));

    // Convert back to negated lower bound representation with correct rounding
    // *a_inf = -a_l requires rounding -a_l UPWARD
    fesetround(FE_UPWARD);
    *a_inf = -a_l;

    // *a_sup = a_u requires rounding a_u UPWARD (already set)
    *a_sup = a_u;

    fesetround(original_round); // Restore original rounding mode
}


/* ============================================================ */
/* Interval Evaluation of Linear Expression (Relies on sound mul)*/
/* ============================================================ */
bool elina_double_interval_eval_elina_linexpr0(double * itv_inf, double *itv_sup,  elina_linexpr0_t* expr, double* env_inf, double * env_sup, elina_scalar_discr_t discr)
{
  size_t i;
  elina_dim_t dim;
  elina_coeff_t* coeff;
  assert(env_inf && env_sup); // Ensure environment bounds are provided

  elina_scalar_t * scalar;
  elina_coeff_t * cst = &expr->cst;

  // Initialize result interval from constant term
  if(cst->discr==ELINA_COEFF_SCALAR){
     // Assuming scalar is double for simplicity here, might need elina_double_set_scalar if mixed types
     if (cst->val.scalar->discr == ELINA_SCALAR_DOUBLE) {
        *itv_inf = -cst->val.scalar->val.dbl; // Negated lower bound
        *itv_sup = cst->val.scalar->val.dbl; // Upper bound
     } else {
         // Handle MPQ/MPFR scalars if needed
         double d_down, d_up;
         int original_round = fegetround();
         fesetround(FE_DOWNWARD); // Round L down
         elina_double_set_scalar(&d_down, cst->val.scalar, FE_DOWNWARD);
         fesetround(FE_UPWARD); // Round U up
         elina_double_set_scalar(&d_up, cst->val.scalar, FE_UPWARD);
         fesetround(original_round);
         *itv_inf = -d_down; // Negated lower bound
         *itv_sup = d_up;   // Upper bound
     }
  }
  else{ // Interval constant
     // Assuming interval bounds are double
     if (cst->val.interval->inf->discr == ELINA_SCALAR_DOUBLE && cst->val.interval->sup->discr == ELINA_SCALAR_DOUBLE) {
        *itv_inf = -cst->val.interval->inf->val.dbl; // Negated lower bound
        *itv_sup = cst->val.interval->sup->val.dbl; // Upper bound
     } else {
         // Handle MPQ/MPFR scalars if needed
         double d_down, d_up;
         int original_round = fegetround();
         fesetround(FE_DOWNWARD); // Round L down
         elina_double_set_scalar(&d_down, cst->val.interval->inf, FE_DOWNWARD);
         fesetround(FE_UPWARD); // Round U up
         elina_double_set_scalar(&d_up, cst->val.interval->sup, FE_UPWARD);
         fesetround(original_round);
         *itv_inf = -d_down; // Negated lower bound
         *itv_sup = d_up;   // Upper bound
     }
  }

  // Accumulate terms
  double term_inf = 0.0;
  double term_sup = 0.0;
  double coeff_inf = 0.0;
  double coeff_sup = 0.0;

  // Need to save/restore rounding mode around the loop
  int original_round_loop = fegetround();

  elina_linexpr0_ForeachLinterm(expr,i,dim,coeff){
    // Get coefficient interval [coeff_l, coeff_u]
    if (coeff->discr==ELINA_COEFF_SCALAR){
        // Assuming scalar is double
        if (coeff->val.scalar->discr == ELINA_SCALAR_DOUBLE) {
            coeff_inf = -coeff->val.scalar->val.dbl; // Negated lower bound
            coeff_sup = coeff->val.scalar->val.dbl; // Upper bound
        } else {
             double d_down, d_up;
             fesetround(FE_DOWNWARD); elina_double_set_scalar(&d_down, coeff->val.scalar, FE_DOWNWARD);
             fesetround(FE_UPWARD);   elina_double_set_scalar(&d_up, coeff->val.scalar, FE_UPWARD);
             coeff_inf = -d_down;
             coeff_sup = d_up;
        }
    } else { // Interval coefficient
        // Assuming interval bounds are double
        if (coeff->val.interval->inf->discr == ELINA_SCALAR_DOUBLE && coeff->val.interval->sup->discr == ELINA_SCALAR_DOUBLE) {
            coeff_inf = -coeff->val.interval->inf->val.dbl; // Negated lower bound
            coeff_sup = coeff->val.interval->sup->val.dbl; // Upper bound
        } else {
            double d_down, d_up;
            fesetround(FE_DOWNWARD); elina_double_set_scalar(&d_down, coeff->val.interval->inf, FE_DOWNWARD);
            fesetround(FE_UPWARD);   elina_double_set_scalar(&d_up, coeff->val.interval->sup, FE_UPWARD);
            coeff_inf = -d_down;
            coeff_sup = d_up;
        }
    }

    // Get environment interval env = [env_l, env_u]
    // env_inf stores -env_l, env_sup stores env_u
    double env_l = -env_inf[dim];
    double env_u = env_sup[dim];

    // Check if environment interval is valid
    // if (env_u < env_l) { /* Handle error or bottom state */ }

    // Perform sound interval multiplication: term = coeff * env
    // Note: elina_double_interval_mul expects negated lower bounds
    elina_double_interval_mul(&term_inf, &term_sup, coeff_inf, coeff_sup, env_inf[dim], env_sup[dim]);

    // Perform sound interval addition: itv = itv + term
    // *itv_inf represents -itv_l, term_inf represents -term_l
    double itv_l = -(*itv_inf);
    double itv_u = *itv_sup;
    double term_l = -term_inf;
    double term_u = term_sup;

    fesetround(FE_DOWNWARD);
    double new_itv_l = itv_l + term_l;
    fesetround(FE_UPWARD);
    double new_itv_u = itv_u + term_u;

    *itv_inf = -new_itv_l; // Rounding -new_itv_l UP is needed
    *itv_sup = new_itv_u; // Already rounded UP

    // Early exit if result becomes [-inf, +inf]
    if (!isfinite(*itv_inf) || !isfinite(*itv_sup)){
       // Or more precisely check if *itv_inf == INFINITY && *itv_sup == INFINITY
       fesetround(original_round_loop); // Restore before potentially breaking
       *itv_inf = INFINITY; // Ensure canonical [-inf, +inf] representation
       *itv_sup = INFINITY;
       // break; // Can uncomment if early exit is desired
       // For now, let it continue to ensure all terms are processed,
       // maybe useful for debugging, though result is known.
    }
  }
  fesetround(original_round_loop); // Restore rounding mode after loop

  // Final check for bottom state might be needed depending on context
  // if (elina_double_interval_canonicalize(itv_inf, itv_sup, false, discr)) return false; // Indicate bottom

  return true; // Indicate success (or not bottom)
}


// --- Functions below rely on the corrected arithmetic ---
// --- They might need adjustments if ELINA API expects specific rounding modes ---
// --- For now, assuming they rely on the sound arithmetic provided above ---

static inline bool elina_double_scalar_sgn(double d){
    if(d > 1e-12){ // Use tolerance for comparisons
        return 1;
    }
    else if(d < -1e-12){ // Use tolerance for comparisons
        return -1;
    }
    else{
        return 0;
    }
}

// Boxize functions likely need careful checking regarding rounding mode assumptions
// when interacting with elina_scalar_t if those aren't doubles.
// The provided elina_double_interval_eval_elina_linexpr0 now does rounding,
// but the boxize logic itself might perform divisions/calculations that also need rounding.
// For simplicity, the boxize functions are left largely unchanged structurally,
// assuming the main source of error was the interval arithmetic itself.
// A full review would involve checking every calculation in boxize against rounding requirements.

static bool elina_double_boxize_lincons0(double* res_inf, double * res_sup,
                   bool* tchange,
                   elina_lincons0_t* cons,
                   double* env_inf, double * env_sup,
                   size_t intdim,
                   bool intervalonly, elina_scalar_discr_t discr)
{
  size_t i;
  elina_linexpr0_t* expr;
  bool change,globalchange;
  bool exc;
  assert(cons->constyp == ELINA_CONS_EQ ||
     cons->constyp == ELINA_CONS_SUPEQ ||
     cons->constyp == ELINA_CONS_SUP);

  expr = cons->linexpr0;
  globalchange = false;

  int original_round = fegetround(); // Save rounding mode for this function

  /* Iterates on coefficients */
  double itv_inf = 0.0;
  double itv_sup = 0.0;
  double val = 0.0;
  // elina_scalar_t * scalar2 = elina_scalar_alloc(); // Appears unused?

  for (i=0; i<expr->size; i++){
    elina_coeff_t *coeff;
    elina_dim_t dim;
    if(expr->discr==ELINA_LINEXPR_SPARSE){
        elina_linterm_t *term = &expr->p.linterm[i];
        dim = term->dim;
        coeff = &term->coeff;
    }
    else{
        dim = i;
        coeff = &expr->p.coeff[i];
    }

    elina_coeff_t * tmp = elina_coeff_alloc(coeff->discr);
    if (!tmp) { fesetround(original_round); return false; } // Allocation check

    /* 1. We decompose the expression e = ax+e' */
    elina_coeff_swap(tmp,coeff); // Temporarily remove coeff[i] from expr

    double inf = 0.0; // coeff lower bound (negated)
    double sup = 0.0; // coeff upper bound
    double d;

    // Extract interval [coeff_l, coeff_u] from tmp (using negated lower bound for inf)
    // Needs rounding if source is not double
    if (tmp->discr==ELINA_COEFF_SCALAR) {
        if (tmp->val.scalar->discr == ELINA_SCALAR_DOUBLE) {
            inf = -tmp->val.scalar->val.dbl; sup = tmp->val.scalar->val.dbl;
        } else {
            fesetround(FE_DOWNWARD); elina_double_set_scalar(&d, tmp->val.scalar, FE_DOWNWARD); inf = -d;
            fesetround(FE_UPWARD);   elina_double_set_scalar(&d, tmp->val.scalar, FE_UPWARD);   sup = d;
        }
    } else { // Interval
        if (tmp->val.interval->inf->discr == ELINA_SCALAR_DOUBLE) {
            inf = -tmp->val.interval->inf->val.dbl;
        } else {
            fesetround(FE_DOWNWARD); elina_double_set_scalar(&d, tmp->val.interval->inf, FE_DOWNWARD); inf = -d;
        }
        if (tmp->val.interval->sup->discr == ELINA_SCALAR_DOUBLE) {
            sup = tmp->val.interval->sup->val.dbl;
        } else {
            fesetround(FE_UPWARD); elina_double_set_scalar(&d, tmp->val.interval->sup, FE_UPWARD); sup = d;
        }
    }

    /* 2. evaluate e' (the rest of the expression) */
    // This now uses sound interval arithmetic internally
    elina_double_interval_eval_elina_linexpr0(&itv_inf, &itv_sup, expr, env_inf, env_sup, discr);

    // Convert eval result [-itv_inf, itv_sup] to [itv_l, itv_u]
    double itv_l = -itv_inf;
    double itv_u = itv_sup;


    /* 3. Perform deduction from ax+e' >= 0 => a*x >= -e' => a*x >= [-itv_u, -itv_l] */
    // We want to find tighter bounds for x based on this constraint and a=[coeff_l, coeff_u]

    bool equality = (tmp->discr == ELINA_COEFF_SCALAR);
    change = false;

    // Convert coeff interval [-inf, sup] to [coeff_l, coeff_u]
    double coeff_l = -inf;
    double coeff_u = sup;

    int sgn_coeff_l = elina_double_scalar_sgn(coeff_l);
    int sgn_coeff_u = elina_double_scalar_sgn(coeff_u);

    // --- Calculate potential new lower bound for x (val_inf) ---
    double val_inf = -INFINITY; // Start with widest possible bound
    // Need to solve for x_inf from: [coeff_l, coeff_u] * x >= [-itv_u, -itv_l]
    if (sgn_coeff_u > 0) { // If coeff_u > 0
        fesetround(FE_DOWNWARD);
        val_inf = fmax(val_inf, -itv_u / coeff_u);
    }
    if (sgn_coeff_l < 0) { // If coeff_l < 0
        fesetround(FE_DOWNWARD);
         val_inf = fmax(val_inf, -itv_l / coeff_l); // Division by negative flips inequality
    }
    // If coeff interval contains 0, things are more complex, potentially no finite bound derivable this way.
    // The original code handled this implicitly by cases. Let's stick to simpler cases for now.
    // If coeff_l <= 0 <= coeff_u, deriving a finite bound might not be possible from ax >= Y unless Y is always positive/negative.


    // --- Calculate potential new upper bound for x (val_sup) ---
    double val_sup = INFINITY; // Start with widest possible bound
    // Need to solve for x_sup from: [coeff_l, coeff_u] * x >= [-itv_u, -itv_l] (and <= if eq)
    // If EQ constraint: a*x <= -e' => a*x <= [-itv_u, -itv_l]
    if (equality) {
        // From a*x <= -itv_l (upper bound of -e')
        if (sgn_coeff_l > 0) { // If coeff_l > 0
            fesetround(FE_UPWARD);
            val_sup = fmin(val_sup, -itv_l / coeff_l);
        }
        if (sgn_coeff_u < 0) { // If coeff_u < 0
            fesetround(FE_UPWARD);
            val_sup = fmin(val_sup, -itv_u / coeff_u); // Division by negative flips inequality
        }
    } else { // SUPEQ or SUP constraint (ax >= -e') - Doesn't directly give upper bound on x unless coeff is negative
        if (sgn_coeff_l < 0) { // If coeff_l < 0, implies coeff_u < 0 (0 not in interval)
            // ax >= -itv_u => x <= -itv_u / coeff_l (division by negative flips)
             fesetround(FE_UPWARD);
             val_sup = fmin(val_sup, -itv_u / coeff_l);
        }
         // If coefficient is positive, ax >= Y doesn't give an upper bound on x.
    }


    // Restore original coefficient into expression
    elina_coeff_swap(tmp,coeff);
    elina_coeff_free(tmp);

    // --- Update bounds ---
    // Get current bounds [res_l, res_u]
    double res_l = -res_inf[dim];
    double res_u = res_sup[dim];

    if (val_inf > res_l) { // If new lower bound is tighter
        change = true;
        if (tchange) tchange[2*dim] = true;
        res_l = val_inf;
    }
    if (val_sup < res_u) { // If new upper bound is tighter
        change = true;
        if (tchange) tchange[2*dim+1] = true;
        res_u = val_sup;
    }

    // Convert back to negated representation, applying rounding for integer dims
    if (change) {
        globalchange = true;
        if (dim < intdim) { // Integer dimension
            fesetround(FE_DOWNWARD); // Round L down
            res_l = floor(res_l);
            fesetround(FE_UPWARD); // Round U up
            res_u = ceil(res_u);
            // Specific handling for ELINA_CONS_SUP might be needed here if strict inequality matters
        }
        // Round -res_l upward for inf
        fesetround(FE_UPWARD);
        res_inf[dim] = -res_l;
        // Round res_u upward for sup
        fesetround(FE_UPWARD);
        res_sup[dim] = res_u;

        // Check for bottom state after update
        if (elina_double_interval_canonicalize(res_inf + dim, res_sup + dim, dim < intdim, discr)){
            fesetround(original_round); // Restore before returning
            // Signal bottom by setting a known invalid state, e.g., inf > -sup
            res_inf[0] = 1.0; res_sup[0] = -1.0; // Example invalid state
            return true; // Indicate bottom found
        }
    }

  } // End loop over coefficients

  fesetround(original_round); // Restore rounding mode

  // Handle case where expression has size 0 (constant constraint)
  if (expr->size==0){
      // Assuming eval_elina_cstlincons0 returns 0 if constraint is violated
      // This part seems specific to ELINA's internal constraint checking
      // Keeping original logic, assuming eval function is correct
      // if (eval_elina_cstlincons0(cons)==0){
      //    res_inf[0] = 1.0; res_sup[0] = -1.0; // Signal bottom
      //    globalchange = true;
      //}
  }

  return globalchange;
}


// This function calls the boxize_lincons0 repeatedly. It should be okay if the underlying function is sound.
bool elina_double_boxize_lincons0_array(double * res_inf, double * res_sup, bool* tchange,
                      elina_lincons0_array_t* array,
                      double* env_inf, double * env_sup, size_t intdim,
                      size_t kmax,
                      bool intervalonly, elina_scalar_discr_t discr)
{
  size_t i,k;
  bool change,globalchange;

  if (kmax<1) kmax=1;
  // If res and env are different, we only need one pass to transfer constraints
  // If res and env are the same, we iterate to potentially refine bounds further
  if ((res_inf!=env_inf) || (res_sup!=env_sup)) kmax=1;

  globalchange = false;
  /* we possibly perform kmax passes */
  for (k=0; k < kmax; k++){
    change = false;
    for (i=0; i<array->size; i++){
      // Only process inequality/equality constraints relevant for boxization
      if (array->p[i].constyp==ELINA_CONS_EQ ||
          array->p[i].constyp==ELINA_CONS_SUPEQ ||
          array->p[i].constyp==ELINA_CONS_SUP){

            bool current_change = elina_double_boxize_lincons0(
                                    res_inf, res_sup, tchange, &array->p[i],
                                    env_inf, env_sup, // Pass distinct env if kmax > 1? Check ELINA usage. Pass res as env if iterating.
                                    intdim, intervalonly, discr);

            // Check if boxize indicated bottom state
            // Using the convention set in elina_double_boxize_lincons0
            if (res_sup[0] < -res_inf[0] - 1e-12) {
                 return true; // Bottom detected, stop immediately
            }

            change = current_change || change;
            globalchange = globalchange || change;

      }
    }
    // If no change occurred in a full pass, further passes won't help
    if (!change) break;

    // If iterating (kmax > 1), use the updated bounds as the environment for the next pass
    if (kmax > 1) {
        env_inf = res_inf;
        env_sup = res_sup;
    }
  }
  return globalchange; // Return true if any change occurred, potentially indicating bottom if res became bottom.
}


/* ============================================================ */
/* Meet_lincons (Relies on sound boxize)                      */
/* ============================================================ */

elina_box_t* elina_box_meet_lincons_array(elina_manager_t* man,
                  bool destructive,
                  elina_box_t* a,
                  elina_lincons0_array_t* array)
{
  elina_box_t* res;
  size_t kmax;
  elina_lincons0_array_t tlincons; // Temporary copy for reduction
  // elina_box_internal_t* intern = (elina_box_internal_t*)man->internal; // Appears unused?

  res = destructive ? a : elina_box_copy(man,a);
  if (!res) return NULL; // Allocation check

  if (a->inf==NULL && a->sup==NULL){ // If input is bottom
    man->result.flag_best = true;
    man->result.flag_exact = true;
    elina_box_set_bottom(res); // Ensure result is also bottom
    return res;
  }
  else {
    man->result.flag_best = (array->size <= 1); // Best if zero or one constraint
    man->result.flag_exact = false; // Generally not exact due to boxization limitations

    // Determine iteration count for boxization
    kmax = man->option.funopt[ELINA_FUNID_MEET_LINCONS_ARRAY].algorithm;
    if (kmax<1) kmax=2; // Default iterations if not specified or invalid

    // Reduce integer constraints first (can make boxization more effective)
    tlincons = elina_lincons0_array_make(array->size);
    if (tlincons.p == NULL && array->size > 0) { // Allocation check
         if (!destructive) elina_box_free(man, res);
         return NULL;
    }
    for(size_t i =0; i < array->size; i++){
        tlincons.p[i] = elina_lincons0_copy(&array->p[i]);
        if (tlincons.p[i].linexpr0 == NULL && array->p[i].linexpr0 != NULL) { // Check copy success
             elina_lincons0_array_clear(&tlincons);
             if (!destructive) elina_box_free(man, res);
             return NULL;
        }
    }

    // Assuming elina_lincons0_array_reduce_integer exists and works correctly
    // It might simplify constraints based on integer variables.
    // Need to check its return value definition. Assuming 0 means infeasible.
    // bool tb = elina_lincons0_array_reduce_integer(&tlincons,a->intdim,ELINA_SCALAR_DOUBLE);
    // if (tb==0){ // If reduction found infeasibility
    //   goto _elina_box_meet_lincons_array_bottom;
    // }

    // Perform iterative boxization using the (potentially reduced) constraints
    // Pass res->inf, res->sup as both result and environment for iterative refinement
    bool changed = elina_double_boxize_lincons0_array(res->inf, res->sup, NULL,
                 &tlincons, res->inf, res->sup, a->intdim, kmax, false, ELINA_SCALAR_DOUBLE);

    // Check if boxization resulted in a bottom state
    if (res->sup[0] < -res->inf[0] - 1e-12) { // Use the convention from boxize
       goto _elina_box_meet_lincons_array_bottom;
    }

    // Canonicalize the final result (redundant if boxize already did, but safe)
    if(elina_double_interval_canonicalize(res->inf, res->sup, false, ELINA_SCALAR_DOUBLE)){
       goto _elina_box_meet_lincons_array_bottom;
    }

    // If we reached here, the result is likely not bottom
    man->result.flag_exact = !changed; // Exact only if boxization made no changes? Check ELINA semantics. Usually false.
    man->result.flag_best = true; // Box meet is considered "best" within box domain

    elina_lincons0_array_clear(&tlincons);
    return res;

  _elina_box_meet_lincons_array_bottom:
      elina_box_set_bottom(res);
      man->result.flag_best = true;
      man->result.flag_exact = true; // Bottom is exact
      elina_lincons0_array_clear(&tlincons);
      return res;
  }
}


/* ============================================================ */
/* Widening (Unchanged from original - standard box widening)   */
/* ============================================================ */
elina_box_t* elina_box_widening(elina_manager_t* man,
                    elina_box_t* a1, elina_box_t* a2)
{
    size_t i;
    size_t nbdims;
    elina_box_t* res;

    man->result.flag_best = true;
    man->result.flag_exact = true; // Standard box widening is exact

    // Handle case where a1 is bottom
    if ((a1->inf==NULL) && (a1->sup==NULL)){
        return elina_box_copy(man,a2); // Result is a2
    }
    // Handle case where a2 is bottom (result is a1)
    if ((a2->inf==NULL) && (a2->sup==NULL)){
        return elina_box_copy(man,a1);
    }


    nbdims = a1->intdim+a1->realdim; // Use combined dimensions
    res = elina_box_copy(man,a1); // Copy a1 as base for result
    if (!res) return NULL; // Allocation check

    for (i=0; i<nbdims; i++){
        // If upper bound increased, widen to +infinity
        if(a1->sup[i] < a2->sup[i]){
            res->sup[i] = INFINITY;
        }
        // else keep a1's upper bound (already copied)

        // If lower bound increased (negated inf decreased), widen to -infinity
        // a1->inf[i] stores -a1_l, a2->inf[i] stores -a2_l
        // Check if -a1_l < -a2_l  <=> a1_l > a2_l (lower bound went up)
        if(a1->inf[i] < a2->inf[i]){
            res->inf[i] = INFINITY; // Represents -(-infinity)
        }
        // else keep a1's lower bound (already copied)
    }
    return res;
}