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
  
  //0_1
  Rcpp::Range tv = Rcpp::seq(0,1);
  //1_2
  Rcpp::Range nind = Rcpp::seq(1,nvdi-1);
  //0_1_2
  Rcpp::Range nlag = Rcpp::seq(0,MAXl);
  
  //2*3*2 = 12
  int dimtot = tv.size()*nlag.size()*nind.size();

  //12*4 patterns in livello
  arma::mat lev_patterns(dimtot,4,fill::zeros);
  //12*4 patterns in differenza
  arma::mat diff_patterns(dimtot,4,fill::zeros);
  int t = 0;
  for(int i=0; i < nlag.size(); ++i){ // ciclo 0_1_2 lags
    for(int j=1; j <= nind.size(); ++j){ //ciclo 1_2 indip
      for(int k=0; k < tv.size(); ++k){ //ciclo 0_1 tv
        
        //creo vettore indice
        arma::rowvec idx(4,fill::zeros);
        //contiene il pattern es. tv = 0, ind = 1, lag = 0
        idx(1) = k;
        idx(2) = j;
        idx(3) = i;
        if(!(k==0 && j>1)){
          //mi servono tutti i lag
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
//last row = constant (-1)
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

   for(int w = nstart ; w < nss; ++w){
       if(w > 1){
          //shift delle osservazioni al passo precedente
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
 
//riempimento colonne per modello marginale x di ARDL
   row_sel(0) = w;
//seleziono variabili diff in lag per z
   arma::rowvec v_dlag = df_oss.submat(row_sel,criteria_dl);
//seleziono variabili lagged levels per z
   arma::rowvec v_lag = df_oss.submat(row_sel,criteria_lvll);
//diff istantanee delle x
   arma::rowvec idx0(nvdi);
   for(int j = 0; j < nvdi; ++j){
      idx0(j) = w;
    }
   df_oss.submat(row_sel,criteria_dindinst)=
     ut.submat(w,1,w,(nvdi-1)) +
     interc.subvec(1,(nvdi-1)).t() +
     (trend.subvec(1,(nvdi-1)).t())%idx0.subvec(1,(nvdi-1))+ //errore + trend + intercetta per x
     v_dlag*(GAMMAX.rows(1,GAMMAX.n_rows-1)).t()+ //diff in lag per z
     v_lag*(PIM.rows(1,PIM.n_rows-1)).t(); //lag levels per z lungo periodo
   
   //levels per x
   df_oss.submat(w,1,w,(nvdi-1)) =
     df_oss.submat(w-1,1,w-1,(nvdi-1))+ //riga precedente
     df_oss.submat(row_sel,criteria_dindinst); //riga appena creata

   //seleziono variabili diff istantanee per x
   arma::rowvec x_diff = df_oss(row_sel,criteria_dindinst);

   //riempimento colonne per modello CONDIZIONATO y di ARDL
   vec newy =
     e_cond(w) + //errore in formulazione condizionata
     interc.at(0) + // intercept conditional
     trend.at(0) + //trend conditional
     (v_lag)*ayx.t() +  //relazione di lungo periodo in formulazione condizionata
     x_diff*omegat + //differenze istantanee x in formulazione condizionata
     v_dlag*psi.t(); //diff in lag per z in formulazione condizionata
   df_oss.at(w,(2*nvdi)) = newy(0);

   //levels per y
   df_oss.at(w,0) =
     df_oss.at(w-1,0) + //riga precedente
     df_oss.at(w,(2*nvdi)); //riga appena creata
}

   return Rcpp::List::create(
     Rcpp::Named("df") = df_oss
   );
}

