#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::rowvec projection_vector_directions(arma::vec xxT, arma::mat directions){
  
   int nrow = sqrt(xxT.n_elem);  
   arma::mat mxxT;
   mxxT.insert_cols(0, xxT);
   mxxT.reshape(nrow, nrow);
   
   //int p = directions.n_rows;
   int ndir = directions.n_cols;
   
   arma::rowvec aprojection(ndir);

   for (int i=0; i < ndir; i++) {
        arma::mat A = directions.col(i).t()*mxxT*directions.col(i);
        double aa = A(0, 0);  
        aprojection(i) = aa ;
    }  
  return(aprojection);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat projection_directions(arma::mat XXT, arma::mat directions){
  
   //int nrow = XXT.n_rows;
   int ncol = XXT.n_cols;
   int ndir = directions.n_cols;
   
   arma::mat allprojection(ncol, ndir);
   for (int i=0; i < ncol; i++) {
      arma::vec acov = XXT.col(i);
      allprojection.row(i) = projection_vector_directions(acov, directions);
   }
   return(allprojection.t());
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::rowvec projection_a_direction(arma::mat XXT, arma::vec adirection){
  
   int nrow = XXT.n_rows;
   int ncol = XXT.n_cols;
   int p = sqrt(nrow);
  
   arma::rowvec projection(ncol);

   for (int i=0; i < ncol; i++) {
      arma::mat mxxT;
      mxxT.insert_cols(0, XXT.col(i));
      mxxT.reshape(p, p);
      arma::mat A = adirection.t()*mxxT*adirection;
      double aa = A(0, 0);  
      projection(i) = aa ;
   }
   return(projection);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List find_depth_on_u_move(arma::vec u_move, arma::mat joined_Minter){
  
   //int nrow = joined_Minter.n_rows;
   int ncol = joined_Minter.n_cols;
        
   arma::rowvec project_to_u_move = projection_a_direction(joined_Minter, u_move);
   arma::uvec sorted = sort_index(sort_index(project_to_u_move));
   int nelem = sorted.n_elem;
   arma::urowvec rank_u_move(nelem);
   for(int i=0; i < nelem; i++){
     rank_u_move(i) = sorted(i);
   }
   arma::urowvec depth_on_u_move(ncol); 
   //Rcpp::Rcout << "sort 1" << std::endl;    
   for (int i=0; i < ncol; i++) {
      if(rank_u_move(i) > ncol/2){
         depth_on_u_move(i) = ncol - rank_u_move(i);
      } else {
         depth_on_u_move(i) = rank_u_move(i);
      }
    }
    
   return Rcpp::List::create(Rcpp::Named("depth_on_u_move") = depth_on_u_move,
                             Rcpp::Named("u_move") = u_move,
                             Rcpp::Named("project_to_u_move") = project_to_u_move); 
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List find_depth_on_multiple_u_move(arma::mat u_min, arma::mat joined_Minter){
       
   // int nrow = joined_Minter.n_rows;
   // int ncol = joined_Minter.n_cols;
   
   // Rcpp::Rcout << "Mark 2" << std::endl;
    arma::uvec aa = arma::randi<arma::uvec>(1, arma::distr_param(0, u_min.n_cols-1));
    arma::vec u_move = vectorise(u_min.col(aa(0)));  
    arma::uvec all(u_min.n_cols);    
    for (int i=0; i < u_min.n_cols; i++) {
        all(i) = i;
    }
        
    arma::uvec others = all.elem(arma::find(all != aa(0)));
    arma::mat backup = u_min.cols(others);
        
    Rcpp::List res = find_depth_on_u_move(u_move, joined_Minter); 
    res["backup"] = backup;
    return(res);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List find_depth_on_all_directions(arma::mat M, arma::mat projectto, arma::mat directions){
  
    arma::rowvec Minter_projection = projection_vector_directions(M, directions);
    arma::mat joined = join_rows(projectto, Minter_projection.t());
    
    int nrow = joined.n_rows;
    int ncol = joined.n_cols;
    
   // Rcpp::Rcout << "Mark 1" << std::endl;
    arma::umat projectMrank(nrow, ncol);
    for (int i=0; i < nrow; i++) {
        projectMrank.row(i) = trans(sort_index(sort_index(joined.row(i))));    
    }
    arma::uvec Minter_depth(ncol);  
    arma::uvec Minter_rank(ncol); 
   // Rcpp::Rcout << "Mark 2" << std::endl;
    Minter_depth = Minter_rank = projectMrank.col(ncol-1);
    for (int i=0; i < nrow; i++) {
      if(Minter_rank(i) > ncol/2){
          Minter_depth(i) = ncol - Minter_rank(i);
      }  
    }  
    // Rcpp::Rcout << "Mark 3" << std::endl;
    return Rcpp::List::create(Rcpp::Named("Minter_depth") = Minter_depth); 
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]


Rcpp::List move_a_step(arma::uvec depth_on_u_move, int Minter_depth_u_move, arma::vec u_move, arma::mat Minter, arma::rowvec project_to_u_move){
  
    int ind = project_to_u_move.n_elem;
    int towards = depth_on_u_move.elem(arma::find(depth_on_u_move > Minter_depth_u_move)).min();
    arma::vec dis = project_to_u_move(ind-1) - project_to_u_move(arma::find(depth_on_u_move == towards));
    arma::vec step = dis(arma::find(abs(dis) == abs(dis).min())); 
    
    //Rcpp::Rcout << "step = " << step << std::endl; 
    if (std::abs(step(0)) < 0.0000001){
       towards = towards + 1;
       dis = project_to_u_move(ind-1) - project_to_u_move(arma::find(depth_on_u_move == towards));
       step = dis(arma::find(abs(dis) == abs(dis).min()));  
    }
    arma::mat M = Minter - arma::vectorise(step(0)*u_move*u_move.t());
    
    return Rcpp::List::create(Rcpp::Named("step") = step, Rcpp::Named("M") = M);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List descent_algorithm_matrix_depth(arma::mat XXT, arma::vec Minitial, arma::mat directions, int ntry, int step_check, bool print){
  
   
  
   arma::mat projectto = projection_directions(XXT, directions);
   int ncol = projectto.n_cols + 1;
   arma::mat Minter = Minitial;
   arma::mat best = Minitial;
   
   
   Rcpp::List depth_on_all_directions = find_depth_on_all_directions(Minter, projectto, directions);
     
  // Rcpp::Rcout << "test 1" << std::endl;
     
   arma::uvec Minter_depth = depth_on_all_directions["Minter_depth"]; 
   int current_depth = Minter_depth.min();
   int best_depth = current_depth;  
   int previous_depth = current_depth;
   int counter = 0;
   int Minter_depth_u_move = 1;
   arma::mat backup = directions.cols(0, 9);
   arma::mat record_u_min1 = directions.cols(0, 9);
   arma::mat record_u_min2 = directions.cols(0, 9);
   arma::mat record_u_min3 = directions.cols(0, 9);
   arma::mat record_u_min4 = directions.cols(0, 9);
   arma::mat record_u_min5 = directions.cols(0, 9);
   arma::mat u_min = directions.cols(0, 9);
   arma::mat u_move = directions.col(0);

   arma::vec print_depth(ntry);
   
   int len = 1;
   int active_backup = 0;
   arma::uvec record_len(5);
   record_len.zeros();
     
   for (int i=0; i < ntry; i++) {
      
      if(counter < step_check){
        
          arma::uvec ind = find(Minter_depth == current_depth); 
          len = ind.n_elem;      
          if(len ==  1) {
            u_min.col(0) = directions.cols(ind);
            active_backup = 0;
          }
          if(len > 1 && len < 10){
            u_min.cols(0, len-1) = directions.cols(ind);
            active_backup = len - 1;
            if(previous_depth != current_depth) {
              record_u_min5 = record_u_min4;
              record_len(4) = record_len(3);
              record_u_min4 = record_u_min3;
              record_len(3) = record_len(2);
              record_u_min3 = record_u_min2;
              record_len(2) = record_len(1);
              record_u_min2 = record_u_min1;
              record_len(1) = record_len(0);
              record_u_min1 = u_min;
              record_len(0) = len;
            }
          }
          if(len >= 10){
            u_min.cols(0, 9) = directions.cols(ind.subvec(0, 9));
            len = 10;
            active_backup = len - 1;
            if(previous_depth != current_depth) {
              record_u_min5 = record_u_min4;
              record_len(4) = record_len(3);
              record_u_min4 = record_u_min3;
              record_len(3) = record_len(2);
              record_u_min3 = record_u_min2;
              record_len(2) = record_len(1);
              record_u_min2 = record_u_min1;
              record_len(1) = record_len(0);
              record_u_min1 = u_min;
              record_len(0) = len;
            }
          }       
      }  
      
      previous_depth = current_depth;
      
      if(len > 1){
       // Rcpp::Rcout << "Mark 4" << std::endl;
        arma::uvec aa = arma::randi<arma::uvec>(1, arma::distr_param(0, len-1));
        u_move = vectorise(u_min.col(aa(0)));  
        arma::uvec all(len);    
        for (int i=0; i < len; i++) {
           all(i) = i;
        }    
        arma::uvec others = all.elem(arma::find(all != aa(0)));
        backup.cols(0, len-2) = u_min.cols(others);   
        active_backup = len - 1;
      } else {
        u_move = u_min.col(0);
        active_backup = 0;
      }
      
     // Rcpp::Rcout << "len" << len << std::endl; 
     // Rcpp::Rcout << counter << std::endl;
      
      arma::mat joined_Minter = join_rows(XXT, Minter);
      
      Rcpp::List on_u_move = find_depth_on_u_move(u_move, joined_Minter);
      arma::uvec depth_on_u_move = on_u_move["depth_on_u_move"];   
      arma::rowvec project_to_u_move = on_u_move["project_to_u_move"];
      Minter_depth_u_move = depth_on_u_move(depth_on_u_move.n_elem-1);
             
     // Rcpp::Rcout << Minter(0) << std::endl;       
     // Rcpp::Rcout << Minter_depth_u_move << std::endl;
     // Rcpp::Rcout << project_to_u_move(0) << std::endl;  
    
      if(counter >= step_check && active_backup > 1) {
        List on_multiple_u_move = find_depth_on_multiple_u_move(backup.cols(0, active_backup-1), joined_Minter);     
        arma::mat new_u_move = on_multiple_u_move["u_move"];
        u_move = new_u_move;
        arma::uvec new_depth_on_u_move = on_multiple_u_move["depth_on_u_move"];
        depth_on_u_move = new_depth_on_u_move;
        Minter_depth_u_move = depth_on_u_move(depth_on_u_move.n_elem-1);
        arma::mat new_backup = on_multiple_u_move["backup"];
        backup.cols(0, active_backup-2) = new_backup;
        len = active_backup;
        active_backup = active_backup - 1;    
        counter = 0;
      } else if(counter >= step_check && active_backup <= 1) {
        
      // Rcpp::Rcout << "mark 1" << std::endl;
         
        arma::uvec large_ind = find(record_len > 1); 
        int ind_len = large_ind.n_elem; 
        arma::uvec bb = arma::randi<arma::uvec>(1, arma::distr_param(0, ind_len-1));
        
       // Rcpp::Rcout << "bb " << bb(0) << std::endl;
        
        if(bb(0) == 0){
          backup.cols(0, record_len(0)-1) = record_u_min1.cols(0, record_len(0)-1);
          active_backup = record_len(0);
        } 
        if(bb(0) == 1){
          backup.cols(0, record_len(1)-1) = record_u_min2.cols(0, record_len(1)-1);
          active_backup = record_len(1);
        }
        if(bb(0) == 2){
          backup.cols(0, record_len(2)-1) = record_u_min3.cols(0, record_len(2)-1);
          active_backup = record_len(2);
        }
        if(bb(0) == 3){
          backup.cols(0, record_len(3)-1) = record_u_min4.cols(0, record_len(3)-1);
          active_backup = record_len(3);
        }
        if(bb(0) == 4){
          backup.cols(0, record_len(4)-1) = record_u_min5.cols(0, record_len(4)-1);
          active_backup = record_len(4);
        }
        
        List on_multiple_u_move = find_depth_on_multiple_u_move(backup.cols(0, active_backup-1), joined_Minter);     
        arma::mat new_u_move = on_multiple_u_move["u_move"];
        u_move = new_u_move;
        arma::uvec new_depth_on_u_move = on_multiple_u_move["depth_on_u_move"];
        depth_on_u_move = new_depth_on_u_move;
        Minter_depth_u_move = depth_on_u_move(depth_on_u_move.n_elem-1);
        arma::mat new_backup = on_multiple_u_move["backup"];
        backup.cols(0, active_backup-2) = new_backup;
        len = active_backup;
        active_backup = active_backup - 1;
        counter = 0;
     } 
        
      while(Minter_depth_u_move >= ncol/2-1){     
        if(active_backup > 1){
          List on_multiple_u_move = find_depth_on_multiple_u_move(backup.cols(0, active_backup-1), joined_Minter);     
          arma::mat new_u_move = on_multiple_u_move["u_move"];
          u_move = new_u_move;
          arma::uvec new_depth_on_u_move = on_multiple_u_move["depth_on_u_move"];
          depth_on_u_move = new_depth_on_u_move;
          Minter_depth_u_move = depth_on_u_move(depth_on_u_move.n_elem-1);
          arma::mat new_backup = on_multiple_u_move["backup"]; 
          backup.cols(0, active_backup-2) = new_backup;
          active_backup = active_backup - 1;  
        } else {         
        arma::uvec large_ind = find(record_len > 1); 
        int ind_len = large_ind.n_elem; 
        arma::uvec bb = arma::randi<arma::uvec>(1, arma::distr_param(0, ind_len-1));
        
        if(bb(0) == 0){
          backup.cols(0, record_len(0)-1) = record_u_min1.cols(0, record_len(0)-1);
          active_backup = record_len(0);
        } 
        if(bb(0) == 1){
          backup.cols(0, record_len(1)-1) = record_u_min2.cols(0, record_len(1)-1);
          active_backup = record_len(1);
        }
        if(bb(0) == 2){
          backup.cols(0, record_len(2)-1) = record_u_min3.cols(0, record_len(2)-1);
          active_backup = record_len(2);
        }
        if(bb(0) == 3){
          backup.cols(0, record_len(3)-1) = record_u_min4.cols(0, record_len(3)-1);
          active_backup = record_len(3);
        }
        if(bb(0) == 4){
          backup.cols(0, record_len(4)-1) = record_u_min5.cols(0, record_len(4)-1);
          active_backup = record_len(4);
        }    
          List on_multiple_u_move = find_depth_on_multiple_u_move(backup.cols(0, active_backup-1), joined_Minter);     
          arma::mat new_u_move = on_multiple_u_move["u_move"];
          u_move = new_u_move;
          arma::uvec new_depth_on_u_move = on_multiple_u_move["depth_on_u_move"];
          depth_on_u_move = new_depth_on_u_move;
          Minter_depth_u_move = depth_on_u_move(depth_on_u_move.n_elem-1);
          arma::mat new_backup = on_multiple_u_move["backup"];
          backup.cols(0, active_backup-2) = new_backup;
          len = active_backup;
          active_backup = active_backup - 1;
        }
      }
      
      if(Minter_depth_u_move >= ncol/2-2) {
          Minitial = Minter; 
      } else { 
          Rcpp::List M_move_a_step = move_a_step(depth_on_u_move, Minter_depth_u_move, u_move, Minter, project_to_u_move);
          arma::mat Minitial_new = M_move_a_step["M"];
          Minitial = Minitial_new;
      }
      
      List depth_on_all_directions = find_depth_on_all_directions(Minitial, projectto, directions);
      arma::uvec Minter_depth_new = depth_on_all_directions["Minter_depth"]; 
      Minter_depth = Minter_depth_new;
      current_depth = Minter_depth.min();
      
      //Rcpp::Rcout << "current" << current_depth << std::endl;
      if(current_depth > best_depth){
         best_depth = current_depth;
         best = Minitial;
         counter = 0;
       }
     
     if(print == TRUE){
       print_depth(i) = current_depth;
     }
     
      counter = counter + 1;
      Minter = Minitial;
     // Rcpp::Rcout << "best" << best_depth << std::endl;
      if(best_depth >= ncol/2) break;
   } 
   if(print == TRUE){
     return Rcpp::List::create(Rcpp::Named("best") = best, Rcpp::Named("best_depth") = best_depth, Rcpp::Named("print_depth") = print_depth);  
   } else {
     return Rcpp::List::create(Rcpp::Named("best") = best, Rcpp::Named("best_depth") = best_depth);  
   }
}   

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat banding_a_matrix(arma::mat M, int k){
  
  int ncol = M.n_cols;
  for (int i = 0; i < ncol-k+1; i++) {
    for(int j = k + i; j < ncol; j++){
      M(i, j) = 0; 
    }
  }
  
  for (int i = k; i < ncol; i++) {
    for(int j = 0; j < i-k+1; j++){
      M(i, j) = 0; 
    }
  }
  
  return(M);
}
