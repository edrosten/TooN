
/*                       
	 Copyright (C) 2005 Tom Drummond

     This library is free software; you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public
     License as published by the Free Software Foundation; either
     version 2.1 of the License, or (at your option) any later version.

     This library is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     Lesser General Public License for more details.

     You should have received a copy of the GNU Lesser General Public
     License along with this library; if not, write to the Free Software
     Foundation, Inc.
     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
*/
#ifndef __BLASOPERATORS_H
#define __BLASOPERATORS_H

extern "C" {
  void dgemm_(const char* TRANSA, const char* TRANSB, int* width, int* height, int* inter,
	      double* alpha, const double* A, int *lda, const double* B, int *ldb,
	      double* beta, double* C, int *ldc);  
}




template <class RET, class LHS, class RHS>
inline void blasmmmult(RET& ret, const LHS& lhs, const RHS& rhs){

  // This awkwardness is the only way I can figure out to get at
  // the is_rowmajor() static function in the layout classes
  typedef typename RET::layout retlayout;
  typedef typename LHS::layout lhslayout;
  typedef typename RHS::layout rhslayout;

  const bool CM = !retlayout::is_rowmajor();
  const bool LCM = !lhslayout::is_rowmajor();
  const bool RCM = !rhslayout::is_rowmajor();

  int width = CM ? ret.num_cols():ret.num_rows();
  int height = CM ? ret.num_rows():ret.num_cols();
  int inter = lhs.num_cols();
  double alpha = 1;
  int lda = CM ? lhs.num_skip() : rhs.num_skip();
  int ldb = CM ? rhs.num_skip() : lhs.num_skip();
  double beta = 0;
  int ldc = ret.num_skip();
  const char* TrA = CM ? (LCM ? "N" : "T") : (RCM ? "T" : "N");
  const char* TrB = CM ? (RCM ? "N" : "T") : (LCM ? "T" : "N");
  
  dgemm_(TrA,TrB,&width,&height,&inter,&alpha,
	 CM ? lhs.get_data_ptr() : rhs.get_data_ptr(),&lda,
	 CM ? rhs.get_data_ptr() : lhs.get_data_ptr(),&ldb,
	 &beta, ret.get_data_ptr(), &ldc);
}
  


#endif
