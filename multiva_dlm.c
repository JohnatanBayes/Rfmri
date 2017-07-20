#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics.h>

void multiDlm(double *Vm0, double *Vc0, int *Npar, double *Vbeta,
              double *series, double *dist, int *dimSer, int *dimDis, 
              int *Nneigh, int *Stop, double *X, int *dimX1, double *n01,
              double *Cbeta1, double *VecProb, int *CONDITIONALMOD){
  
  int p = *Npar, k33 = dimX1[1], s;
  double VarC0 = *Vc0, Vm01 = *Vm0,  Vbeta1 = *Vbeta, Cbeta = *Cbeta1, n0 = *n01;
  
  int k1 = dimSer[0], k11 = dimDis[0], k22 = dimDis[1], Nneighbors = *Nneigh, stop1 = *Stop;
   
  
  
  
  
  
  /*SETING VIEW MATRICES OF THE SERIES AND VOXELS POSITIONS IN THE BRAING IMAGE*/
  gsl_matrix_view DIST = gsl_matrix_view_array(dist, k11, k22);
  gsl_matrix_view SERI = gsl_matrix_view_array(series, k1, k22);
  gsl_matrix_view X1  = gsl_matrix_view_array(X, k1, k33);
  gsl_matrix_view Proba = gsl_matrix_view_array(VecProb, p, k22);
  /*ALLOCATING AUXILIAR VECTORS THAT WILL HELP ON THE NEIGHBORS SELECTION PROCESS*/
  gsl_vector * V1 = gsl_vector_alloc(k1);
  gsl_vector * V2 = gsl_vector_alloc(k1);
  gsl_vector * pos1 = gsl_vector_alloc(k11);
  gsl_vector * pos2 = gsl_vector_alloc(k11);
  gsl_vector * pos2Aux = gsl_vector_alloc(k11);
  gsl_matrix * Y1  = gsl_matrix_alloc(k1, Nneighbors);
  for(int jjj=0; jjj<k22; jjj++){
  //int jjj = jjj1[0];
  int conta = 0;
  /*THE PURPOUS OF THIS BLOCK OF CODE IS TO IDENTIFY CLUSTERS OF RELATED VOXELS - BLOCK 1 */
  gsl_matrix_get_col(V1, &SERI.matrix, jjj);
  gsl_matrix_get_col(pos1, &DIST.matrix, jjj);
  gsl_matrix_set_col(Y1, conta, V1);
  for(int jj=0; jj<k22; jj++){// FOR 2
    
    gsl_matrix_get_col(pos2, &DIST.matrix, jj);
    gsl_vector_memcpy(pos2Aux, pos1);
    gsl_vector_sub(pos2Aux, pos2);
    double dist_aux = gsl_blas_dnrm2(pos2Aux);
    //ONLY TAKE A CANDIDATE VOXEL IF THE DISTANCE IS EQUAL OR LESS THAN 1 OR WHATEVER THE USER WANTED TO FIX
    if(dist_aux==1){
      gsl_matrix_get_col(V2, &SERI.matrix, jj);
      // THE VOXEL AS PART OF THE NEIGHBORHOOD 
      conta = conta + 1;
      gsl_matrix_set_col(Y1, conta, V2);
    }
    
  }//END FOR 2
  //END BLOCK 1
  
  conta = conta + 1;
  
  if(conta==Nneighbors){
      /*ON THIS BLOCK THE DYNAMIC MATRICES m0, c0, S0 ARE CONSTRUCTED ALONG THE 
     THE COVARINCE MATRIX S0, WHO DEFINES THE CORRELATION AMONG THE NEIGHBORIN VOXELS - BLOCK 2*/ 
     int q = conta;
     gsl_vector * y1 = gsl_vector_alloc(q); 
     gsl_matrix * y11 = gsl_matrix_calloc(1, q); 
     gsl_matrix * c0 = gsl_matrix_calloc(p, p);
     gsl_matrix_set_identity(c0);
     gsl_matrix_scale(c0, VarC0);
     gsl_matrix * beta  = gsl_matrix_calloc(p, p);
     gsl_matrix_set_identity(beta);
     gsl_matrix_scale(beta, Vbeta1);
     //gsl_matrix_set(beta, 0, 0, 1);
     gsl_matrix * R1 = gsl_matrix_calloc(p, p);
     gsl_matrix * R1aux = gsl_matrix_calloc(p, p);
     gsl_matrix * m0 = gsl_matrix_calloc(p, q);
     gsl_matrix * m0Aux = gsl_matrix_calloc(p, q);
     gsl_matrix_add_constant(m0, Vm01);    
     gsl_matrix * f1 = gsl_matrix_calloc(1, q);
     gsl_matrix * F1t = gsl_matrix_calloc(p, 1);
     gsl_matrix * F1aux = gsl_matrix_alloc(p, 1);    
     gsl_matrix * Q1 = gsl_matrix_calloc(1, 1);    
     gsl_matrix * e1 = gsl_matrix_calloc(1, q);     
     gsl_matrix * Qinv = gsl_matrix_alloc(1, 1);
     gsl_matrix * A1 = gsl_matrix_calloc(p, 1);
     gsl_matrix * A1t = gsl_matrix_alloc(1, p);
     gsl_matrix * Aaux = gsl_matrix_calloc(p, 1);
     gsl_matrix * S0 = gsl_matrix_calloc(q, q);
     gsl_matrix_set_identity(S0); 
     gsl_matrix * e1t = gsl_matrix_calloc(q, 1);
     gsl_matrix * e2 = gsl_matrix_calloc(q, q);
     gsl_matrix * C1aux = gsl_matrix_calloc(p, p);
     //END BLOCK 2
     
     // Only for testing
     //int w = *w1;
     //gsl_matrix_view Y11t = gsl_matrix_view_array(y11t, 1, 1);     

    //NEW MATRIX WITH THE VOXELS (OR SERIES) SELECTED ON THE NEIGHBORHOOD
    gsl_matrix_view Y11 = gsl_matrix_submatrix(Y1, 0, 0, k1, q);
    // fortesting gsl_matrix_memcpy(&Y11t.matrix, &Y11.matrix);
    
    for(int i=0; i<stop1; i++){
       //fortesting 
       
      gsl_matrix_get_row(y1, &Y11.matrix, i);
      gsl_matrix_set_row(y11, 0, y1);
      //fortesting  gsl_matrix_memcpy(&Y11t.matrix, y11);
      
      // R1 = B*C0*B
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                     1.0, c0, beta,
                     0.0, R1aux);
      
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                     1.0, beta, R1aux,
                     
                     0.0, R1);
      
      //fortesting gsl_matrix_memcpy(&Y11t.matrix, R1);
    
      
      //test****************************
      //gsl_matrix_memcpy (&TEST.matrix, R);
      //endtest************************
      
      //f1 <- F1*m0
      gsl_matrix_view F1 = gsl_matrix_submatrix(&X1.matrix, i, 0, 1, p);
      //fortesting  gsl_matrix_memcpy(&Y11t.matrix, &F1.matrix);
      //test****************************
      //gsl_matrix_memcpy (&TEST.matrix, &F1.matrix);
      //endtest************************
      
      
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                     1.0, &F1.matrix, m0,
                     0.0, f1);
      
      //fortesting gsl_matrix_memcpy(&Y11t.matrix, f1);
    
      //Q1 = F1*C0*t(F1) + 1
      gsl_matrix_transpose_memcpy(F1t, &F1.matrix);
      //test****************************
      //gsl_matrix_memcpy (&TEST.matrix, F1t);
      //endtest************************
      
      
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                     1.0, c0, F1t,
                     0.0, F1aux);
      
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                     1.0, &F1.matrix, F1aux,
                     0.0, Q1);
      
      gsl_matrix_add_constant(Q1, 1);
      //fortesting gsl_matrix_memcpy(&Y11t.matrix, Q1);
      
      //test****************************
      //gsl_matrix_memcpy (&TEST.matrix, Q1);
      //endtest************************
      
      
      //e1 <- y[i] - f1
      //gsl_matrix_view y1 = gsl_matrix_submatrix(&Y1.matrix, i, 0, 1, 1);
      gsl_matrix_sub(y11, f1);
      
      gsl_matrix_memcpy(e1, y11);
      //fortesting gsl_matrix_memcpy(&Y11t.matrix, e1);
      
      //test****************************
      //gsl_matrix_memcpy (&TEST.matrix, e1);
      //endtest************************
      
      
      
      //Q inverse**************
      gsl_permutation * p = gsl_permutation_alloc (1);
      gsl_linalg_LU_decomp (Q1, p, &s);    
      gsl_linalg_LU_invert (Q1, p, Qinv);
      //fortesting gsl_matrix_memcpy(&Y11t.matrix, Qinv);
      //Q inverse**************
      
      //test****************************
      //gsl_matrix_memcpy (&TEST.matrix, Qinv);
      //endtest************************
      
      //A1 <- R1%*%t(F1)%*%solve(Q1)//
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                     1.0, F1t, Qinv,
                     0.0, Aaux);
       
       //gsl_matrix_memcpy(&Y11t.matrix, Aaux);
      //test****************************
      //gsl_matrix_memcpy (&TEST.matrix, Aaux);
      //endtest************************
      
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                     1.0, R1, Aaux,
                     0.0, A1);
      
    //  gsl_matrix_memcpy(&Y11t.matrix, A1);
      
      //test****************************
      //gsl_matrix_memcpy (&TEST.matrix, A1);
      //endtest************************
      
      
      //n1 <- n0 + 1 
      double n1 = Cbeta*n0 + 1;
             
            //fortesting  gsl_matrix_set(&Y11t.matrix, 0,0, n1);
              
      //S1 <- (1/n1)*(beta1*n0*S0 + t(e1)%*%e1/as.numeric(Q1))
      //gsl_matrix_memcpy(S0Aux, S0);
      gsl_matrix_scale(S0, Cbeta*n0);
      gsl_matrix_transpose_memcpy(e1t, e1);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                     1.0, e1t, e1,
                     0.0, e2);
      gsl_matrix_scale(e2, 1/gsl_matrix_get(Q1, 0, 0));
      gsl_matrix_add(S0, e2);
      gsl_matrix_scale(S0, 1/n1);
      n0 = n1;
      //fortesting 
      
      
           
      
      //C1 <- (R1 - A1%*%t(A1)*as.numeric(Q1))
      
      gsl_matrix_transpose_memcpy(A1t, A1);
      
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                     gsl_matrix_get(Q1, 0, 0), A1, A1t,
                     0.0, C1aux);
      
      gsl_matrix_sub (R1 , C1aux);
      // Copying the matrix R on c0
      gsl_matrix_memcpy (c0, R1);
      
     //fortesting gsl_matrix_memcpy(&Y11t.matrix, c0);
      
      
      
      //m1 <- m0 + A1*e1
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                     1.0, A1, e1,
                     0.0, m0Aux);
       
      gsl_matrix_add(m0, m0Aux);
      
      //fortesting gsl_matrix_memcpy(&Y11t.matrix, m0);

      }
    
    //SCALING S0 BY C0[1,1]
    gsl_matrix_scale(S0, gsl_matrix_get(c0, 0, 0));
    // CREATING AUXILIAR MATRICES
    gsl_matrix * SUBMAINV = gsl_matrix_calloc(q-1, q-1); 
    gsl_matrix * SUBMATRIX2T = gsl_matrix_calloc(q-1, 1);
    gsl_matrix * SUBAUX = gsl_matrix_calloc(q-1, 1);
    gsl_matrix * SUBAUX2 = gsl_matrix_calloc(1, 1);
    
    if(CONDITIONALMOD[0]==1){
    // TAKING SUBMATRICES FROM S0*CO[1,1]
    gsl_matrix_view SUBMATRIX = gsl_matrix_submatrix(S0, 1, 1, q-1, q-1);
    gsl_matrix_view SUBMATRIX2 = gsl_matrix_submatrix(S0, 0, 1, 1, q-1);
    gsl_matrix_transpose_memcpy(SUBMATRIX2T, &SUBMATRIX2.matrix);
    // INVERTING THE SUBMATRIX
    gsl_permutation * p = gsl_permutation_alloc (q-1);
    gsl_linalg_LU_decomp (&SUBMATRIX.matrix, p, &s);    
    gsl_linalg_LU_invert (&SUBMATRIX.matrix, p, SUBMAINV);
    
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                   1.0, SUBMAINV, SUBMATRIX2T,
                   0.0, SUBAUX);
    
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                   1.0, &SUBMATRIX2.matrix, SUBAUX,
                   0.0, SUBAUX2);
    
    double SIGMA = gsl_matrix_get(S0, 0, 0) - gsl_matrix_get(SUBAUX2, 0, 0);
    gsl_permutation_free (p);
    
    gsl_matrix_set(&Proba.matrix, 0, jjj, gsl_matrix_get(m0, 0, 0));
    gsl_matrix_set(&Proba.matrix, 1, jjj, SIGMA);
    
    
    }else{
      gsl_matrix_set(&Proba.matrix, 0, jjj, gsl_matrix_get(m0, 0, 0));
      gsl_matrix_set(&Proba.matrix, 1, jjj, gsl_matrix_get(S0, 0, 0));
     }
    
    
    gsl_matrix_free (SUBMAINV);
    gsl_matrix_free (SUBAUX);
    gsl_matrix_free (SUBAUX2);
    gsl_vector_free (y1);
    gsl_matrix_free (y11);
    gsl_matrix_free (c0);
    gsl_matrix_free (beta);
    gsl_matrix_free (R1);
    gsl_matrix_free (R1aux);
    gsl_matrix_free (m0);
    gsl_matrix_free (m0Aux);
    gsl_matrix_free (f1);
    gsl_matrix_free (F1t);
    gsl_matrix_free (F1aux);
    gsl_matrix_free (Q1);
    gsl_matrix_free (A1);
    gsl_matrix_free (A1t);
    gsl_matrix_free (Aaux);
    gsl_matrix_free (S0);
    gsl_matrix_free (e1t);
    gsl_matrix_free (e2);
    gsl_matrix_free (C1aux);
    
   
  }else{
    for(int ii=0; ii<p; ii++){gsl_matrix_set(&Proba.matrix, ii, jjj, 1000);}
  } 
   
    
  }
  
   

  gsl_vector_free (V1);
  gsl_vector_free (V2);
  gsl_vector_free (pos1);
  gsl_vector_free (pos2);
  gsl_vector_free (pos2Aux);
  gsl_matrix_free (Y1);
  


}
