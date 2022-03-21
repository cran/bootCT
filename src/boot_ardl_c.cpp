//' @name boot_ardl_c
//' @title Code to generate data having a VECM-ARDL structure
//' @param r_in input residuals
//' @param GAMMAX short run parameter matrices column bound, first row in conditional form
//' @param A long-run parameter matrix, first row in conditional form
//' @param start_z data matrix of starting point
//' @param omegat parameter vector of the unlagged differences
//' @param interc vecm intercept, first element in conditional form
//' @param trend vecm trend, first element in conditional form
//' @return A list containing the generated dataset
//' @keywords internal

#include "RcppArmadillo.h"
#include "Rcpp.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List boot_ardl_c(arma::mat r_in, //input residuals
                      arma::mat GAMMAX,
                      arma::mat A,
                      arma::mat start_z,
                      arma::vec omegat,
                      arma::vec interc,//estimated intercept
                      arma::vec trend){  //estimated trend

  Rcpp::Environment base("package:base");
  Rcpp::Function intersect = base["intersect"];

  arma::mat ut = r_in;

  //PI matrix
  arma::mat PIM = -A;

  int nvdi= r_in.n_cols;
  int nss = r_in.n_rows;
  int MAXl = GAMMAX.n_cols/A.n_cols;

  //PSI matrix
  arma::rowvec psi = GAMMAX.row(0);

  //ay.x
  arma::rowvec ayx = PIM.row(0);

  //error conditional
  arma::vec e_cond = r_in.col(0);

  int nstart = start_z.n_rows;

  arma::mat df_oss (nss, (2 * nvdi + nvdi * (MAXl+1) + 1),fill::zeros);
  df_oss.rows(0,nstart-1) = start_z;

  //0_1 level or difference
  Rcpp::Range tv = Rcpp::seq(0,1);
  //1_2_..._(d-1)
  Rcpp::Range nind = Rcpp::seq(1,nvdi-1);
  //0_1_..._(maxlag)
  Rcpp::Range nlag = Rcpp::seq(0,MAXl);

  int dimtot = tv.size()*nlag.size()*nind.size();

  //level patterns
  arma::mat lev_patterns(dimtot,4,fill::zeros);
  //diff patterns
  arma::mat diff_patterns(dimtot,4,fill::zeros);
  int t = 0;
  for(int i=0; i < nlag.size(); ++i){ // lagdiff loop
    for(int j=1; j <= nind.size(); ++j){ //indep loop
      for(int k=0; k < tv.size(); ++k){ //level or difference

        //index vector
        arma::rowvec idx(4,fill::zeros);
        //create pattern combinations e.g. tv = 0, ind = 1, lag = 0
        idx(1) = k;
        idx(2) = j;
        idx(3) = i;
        if(!(k==0 && j>1)){
          //need all the lags
          lev_patterns.row(t) = idx;
          diff_patterns.row(t) = idx;
          t++;
        }
      }
    }
  }

  lev_patterns.shed_rows(t,lev_patterns.n_rows-1);
  arma::mat lev_patternsx=lev_patterns;

  diff_patterns.shed_rows(PIM.n_cols*(MAXl+1),diff_patterns.n_rows-1);
  if(MAXl>1){
    lev_patterns.shed_rows(PIM.n_cols*2,lev_patterns.n_rows-1);
  }

  arma::colvec fldiff(diff_patterns.n_rows,fill::ones);
  diff_patterns.col(0) = fldiff;
  arma::rowvec constantx(4);
  constantx = constantx - 1;
  arma::mat patterns = join_cols(lev_patterns,diff_patterns);
  patterns = join_cols(patterns,constantx);

//patterns matrix:
//first column: 0 = level; 1 = first difference
//second column: 0 = dep; 1 = ind
//third column: from 1 to J; 1 only for dep, 1:J for ind
//fourth column: from 0 to I lags,
//               lag0 and lag1 for levels, lag_i 0<i<I for diffs
  arma::uvec criteria_in;
  arma::uvec criteria_out;
  arma::uvec criteria1 = find(patterns.col(0) == 1);
  arma::uvec criteria0 = find(patterns.col(0) == 0);
  arma::uvec criteria0_out;
  arma::uvec criteria0_in;
  arma::uvec criteria1_out;
  arma::uvec criteria1_in;
   arma::uvec row_in(1);
   arma::uvec row_out(1);

   arma::uvec criteriad = find(patterns.col(0) == 1);
   arma::uvec criterialvl = find(patterns.col(0) == 0);
   arma::uvec criterial = find(patterns.col(3) > 0);
   arma::uvec criteria_dl = Rcpp::as<arma::uvec>(intersect(criteriad,criterial));
   arma::uvec criteria_lvll = Rcpp::as<arma::uvec>(intersect(criterialvl,criterial));
   arma::uvec criteriaind = find(patterns.col(1) == 1);
   arma::uvec criteriainst = find(patterns.col(3) == 0);
   arma::uvec criteria_dindinst = Rcpp::as<arma::uvec>(intersect(intersect(criteriaind,criteriainst),criteriad));

   // index vectors will help in reconstructing the matrix
   for(int w = nstart ; w < nss; ++w){
       if(w > 1){
          //shift at previous loop step moving observations in lagged rows
          for(int j = 0; j < MAXl; j++){
             for(int n = 1; n <= MAXl;n++){
               criteria_out = find(patterns.col(3) == j);
               criteria_in = find(patterns.col(3) == (j+n));
               criteria0_out = Rcpp::as<arma::uvec>(intersect(criteria0,criteria_out));
               criteria0_in = Rcpp::as<arma::uvec>(intersect(criteria0,criteria_in));
               criteria1_out = Rcpp::as<arma::uvec>(intersect(criteria1,criteria_out));
               criteria1_in = Rcpp::as<arma::uvec>(intersect(criteria1,criteria_in));
               row_in(0) = w;
               row_out(0) = (w-n);
               if(criteria1_in.n_rows>0){df_oss.submat(row_in,criteria1_in)=df_oss.submat(row_out,criteria1_out);}
               if(criteria0_in.n_rows>0){df_oss.submat(row_in,criteria0_in)=df_oss.submat(row_out,criteria0_out);}
             }
          }
          }
   arma::uvec row_sel(1);

//marginal VECM model for X
   row_sel(0) = w;
//selecting Z variables difflag
   arma::rowvec v_dlag = df_oss.submat(row_sel,criteria_dl);
   //selecting Z variables lagged levels
   arma::rowvec v_lag = df_oss.submat(row_sel,criteria_lvll);
   //selecting X variables unlagged diffs
   arma::rowvec idx0(nvdi);
   for(int j = 0; j < nvdi; ++j){
      idx0(j) = w;
    }
   //filling marginal VECM model for X
   df_oss.submat(row_sel,criteria_dindinst)=
     ut.submat(w,1,w,(nvdi-1)) +
     interc.subvec(1,(nvdi-1)).t() +
     (trend.subvec(1,(nvdi-1)).t())%idx0.subvec(1,(nvdi-1))+ //error + trend + intercept X
     v_dlag*(GAMMAX.rows(1,GAMMAX.n_rows-1)).t()+ //difflag Z
     v_lag*(PIM.rows(1,PIM.n_rows-1)).t(); //lag levels X

   //levels X
   df_oss.submat(w,1,w,(nvdi-1)) =
     df_oss.submat(w-1,1,w-1,(nvdi-1))+ //previous row
     df_oss.submat(row_sel,criteria_dindinst); //current row

   //selecting X variables unlagged diffs
   arma::rowvec x_diff = df_oss(row_sel,criteria_dindinst);

   //filling conditional ARDL model for y
   vec newy =
     e_cond(w) + //conditional error
     interc.at(0) + //conditional intercept
     trend.at(0) + //conditional trend
     (v_lag)*ayx.t() + //conditional long-run
     x_diff*omegat + //conditional unlagged diffs
     v_dlag*psi.t(); //conditional lagged diffs

   df_oss.at(w,(2*nvdi)) = newy(0);

   //levels y
   df_oss.at(w,0) =
     df_oss.at(w-1,0) + //previous row
     df_oss.at(w,(2*nvdi)); //current row
}

   return Rcpp::List::create(
     Rcpp::Named("df") = df_oss
   );
}

