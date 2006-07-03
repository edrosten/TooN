// Generated for J*C*J^T, C symmetric
template <class Accessor1, class Accessor2> inline Matrix<2> transformCovariance(const FixedMatrix<2,2,Accessor1>& A, const FixedMatrix<2,2,Accessor2>& B)
{
	Matrix<2> M;
	M[0][0] = A[0][0]*(2*Dot<1,1>::eval(A[0], B[0]) + A[0][0]*B[0][0])
		+ A[0][1]*A[0][1]*B[1][1];
	M[0][1] = M[1][0] = A[0][0]*(B[0]*A[1])
		+ A[0][1]*(B[1]*A[1]);
	M[1][1] = A[1][0]*(2*Dot<1,1>::eval(A[1], B[0]) + A[1][0]*B[0][0])
		+ A[1][1]*A[1][1]*B[1][1];
	return M;
}

// Generated for J*C*J^T, C symmetric
template <class Accessor1, class Accessor2> inline Matrix<2> transformCovariance(const FixedMatrix<2,3,Accessor1>& A, const FixedMatrix<3,3,Accessor2>& B)
{
	Matrix<2> M;
	M[0][0] = A[0][0]*(2*Dot<1,2>::eval(A[0], B[0]) + A[0][0]*B[0][0])
		+ A[0][1]*(2*Dot<2,2>::eval(A[0], B[1]) + A[0][1]*B[1][1])
		+ A[0][2]*A[0][2]*B[2][2];
	M[0][1] = M[1][0] = A[0][0]*(B[0]*A[1])
		+ A[0][1]*(B[1]*A[1])
		+ A[0][2]*(B[2]*A[1]);
	M[1][1] = A[1][0]*(2*Dot<1,2>::eval(A[1], B[0]) + A[1][0]*B[0][0])
		+ A[1][1]*(2*Dot<2,2>::eval(A[1], B[1]) + A[1][1]*B[1][1])
		+ A[1][2]*A[1][2]*B[2][2];
	return M;
}

// Generated for J*C*J^T, C symmetric
template <class Accessor1, class Accessor2> inline Matrix<3> transformCovariance(const FixedMatrix<3,2,Accessor1>& A, const FixedMatrix<2,2,Accessor2>& B)
{
	Matrix<3> M;
	M[0][0] = A[0][0]*(2*Dot<1,1>::eval(A[0], B[0]) + A[0][0]*B[0][0])
		+ A[0][1]*A[0][1]*B[1][1];
	M[0][1] = M[1][0] = A[0][0]*(B[0]*A[1])
		+ A[0][1]*(B[1]*A[1]);
	M[0][2] = M[2][0] = A[0][0]*(B[0]*A[2])
		+ A[0][1]*(B[1]*A[2]);
	M[1][1] = A[1][0]*(2*Dot<1,1>::eval(A[1], B[0]) + A[1][0]*B[0][0])
		+ A[1][1]*A[1][1]*B[1][1];
	M[1][2] = M[2][1] = A[1][0]*(B[0]*A[2])
		+ A[1][1]*(B[1]*A[2]);
	M[2][2] = A[2][0]*(2*Dot<1,1>::eval(A[2], B[0]) + A[2][0]*B[0][0])
		+ A[2][1]*A[2][1]*B[1][1];
	return M;
}

// Generated for J*C*J^T, C symmetric
template <class Accessor1, class Accessor2> inline Matrix<2> transformCovariance(const FixedMatrix<2,6,Accessor1>& A, const FixedMatrix<6,6,Accessor2>& B)
{
	Matrix<2> M;
	M[0][0] = A[0][0]*(2*Dot<1,5>::eval(A[0], B[0]) + A[0][0]*B[0][0])
		+ A[0][1]*(2*Dot<2,5>::eval(A[0], B[1]) + A[0][1]*B[1][1])
		+ A[0][2]*(2*Dot<3,5>::eval(A[0], B[2]) + A[0][2]*B[2][2])
		+ A[0][3]*(2*Dot<4,5>::eval(A[0], B[3]) + A[0][3]*B[3][3])
		+ A[0][4]*(2*Dot<5,5>::eval(A[0], B[4]) + A[0][4]*B[4][4])
		+ A[0][5]*A[0][5]*B[5][5];
	M[0][1] = M[1][0] = A[0][0]*(B[0]*A[1])
		+ A[0][1]*(B[1]*A[1])
		+ A[0][2]*(B[2]*A[1])
		+ A[0][3]*(B[3]*A[1])
		+ A[0][4]*(B[4]*A[1])
		+ A[0][5]*(B[5]*A[1]);
	M[1][1] = A[1][0]*(2*Dot<1,5>::eval(A[1], B[0]) + A[1][0]*B[0][0])
		+ A[1][1]*(2*Dot<2,5>::eval(A[1], B[1]) + A[1][1]*B[1][1])
		+ A[1][2]*(2*Dot<3,5>::eval(A[1], B[2]) + A[1][2]*B[2][2])
		+ A[1][3]*(2*Dot<4,5>::eval(A[1], B[3]) + A[1][3]*B[3][3])
		+ A[1][4]*(2*Dot<5,5>::eval(A[1], B[4]) + A[1][4]*B[4][4])
		+ A[1][5]*A[1][5]*B[5][5];
	return M;
}

// Generated for J*C*J^T, C symmetric
template <class Accessor1, class Accessor2> inline Matrix<6> transformCovariance(const FixedMatrix<6,2,Accessor1>& A, const FixedMatrix<2,2,Accessor2>& B)
{
	Matrix<6> M;
	M[0][0] = A[0][0]*(2*Dot<1,1>::eval(A[0], B[0]) + A[0][0]*B[0][0])
		+ A[0][1]*A[0][1]*B[1][1];
	M[0][1] = M[1][0] = A[0][0]*(B[0]*A[1])
		+ A[0][1]*(B[1]*A[1]);
	M[0][2] = M[2][0] = A[0][0]*(B[0]*A[2])
		+ A[0][1]*(B[1]*A[2]);
	M[0][3] = M[3][0] = A[0][0]*(B[0]*A[3])
		+ A[0][1]*(B[1]*A[3]);
	M[0][4] = M[4][0] = A[0][0]*(B[0]*A[4])
		+ A[0][1]*(B[1]*A[4]);
	M[0][5] = M[5][0] = A[0][0]*(B[0]*A[5])
		+ A[0][1]*(B[1]*A[5]);
	M[1][1] = A[1][0]*(2*Dot<1,1>::eval(A[1], B[0]) + A[1][0]*B[0][0])
		+ A[1][1]*A[1][1]*B[1][1];
	M[1][2] = M[2][1] = A[1][0]*(B[0]*A[2])
		+ A[1][1]*(B[1]*A[2]);
	M[1][3] = M[3][1] = A[1][0]*(B[0]*A[3])
		+ A[1][1]*(B[1]*A[3]);
	M[1][4] = M[4][1] = A[1][0]*(B[0]*A[4])
		+ A[1][1]*(B[1]*A[4]);
	M[1][5] = M[5][1] = A[1][0]*(B[0]*A[5])
		+ A[1][1]*(B[1]*A[5]);
	M[2][2] = A[2][0]*(2*Dot<1,1>::eval(A[2], B[0]) + A[2][0]*B[0][0])
		+ A[2][1]*A[2][1]*B[1][1];
	M[2][3] = M[3][2] = A[2][0]*(B[0]*A[3])
		+ A[2][1]*(B[1]*A[3]);
	M[2][4] = M[4][2] = A[2][0]*(B[0]*A[4])
		+ A[2][1]*(B[1]*A[4]);
	M[2][5] = M[5][2] = A[2][0]*(B[0]*A[5])
		+ A[2][1]*(B[1]*A[5]);
	M[3][3] = A[3][0]*(2*Dot<1,1>::eval(A[3], B[0]) + A[3][0]*B[0][0])
		+ A[3][1]*A[3][1]*B[1][1];
	M[3][4] = M[4][3] = A[3][0]*(B[0]*A[4])
		+ A[3][1]*(B[1]*A[4]);
	M[3][5] = M[5][3] = A[3][0]*(B[0]*A[5])
		+ A[3][1]*(B[1]*A[5]);
	M[4][4] = A[4][0]*(2*Dot<1,1>::eval(A[4], B[0]) + A[4][0]*B[0][0])
		+ A[4][1]*A[4][1]*B[1][1];
	M[4][5] = M[5][4] = A[4][0]*(B[0]*A[5])
		+ A[4][1]*(B[1]*A[5]);
	M[5][5] = A[5][0]*(2*Dot<1,1>::eval(A[5], B[0]) + A[5][0]*B[0][0])
		+ A[5][1]*A[5][1]*B[1][1];
	return M;
}

// Generated for lower triangular L, inverse diagonal invdiag
template <class A1, class A2, class A3, class A4> inline void cholesky_backsub(const FixedMatrix<6,6,A1>& L, const FixedVector<6,A2>& invdiag, const FixedVector<6,A3>& v, FixedVector<6,A4>& x)
{
	Vector<6> t;
	//forward substitution
	t[0] = invdiag[0]*(v[0]);
	t[1] = invdiag[1]*(v[1] - L[1][0]*t[0]);
	t[2] = invdiag[2]*(v[2] - L[2][0]*t[0] - L[2][1]*t[1]);
	t[3] = invdiag[3]*(v[3] - L[3][0]*t[0] - L[3][1]*t[1] - L[3][2]*t[2]);
	t[4] = invdiag[4]*(v[4] - L[4][0]*t[0] - L[4][1]*t[1] - L[4][2]*t[2] - L[4][3]*t[3]);
	t[5] = invdiag[5]*(v[5] - L[5][0]*t[0] - L[5][1]*t[1] - L[5][2]*t[2] - L[5][3]*t[3] - L[5][4]*t[4]);
	//backward substitution
	x[5] = invdiag[5]*(t[5]);
	x[4] = invdiag[4]*(t[4] - L[5][4]*x[5]);
	x[3] = invdiag[3]*(t[3] - L[4][3]*x[4] - L[5][3]*x[5]);
	x[2] = invdiag[2]*(t[2] - L[3][2]*x[3] - L[4][2]*x[4] - L[5][2]*x[5]);
	x[1] = invdiag[1]*(t[1] - L[2][1]*x[2] - L[3][1]*x[3] - L[4][1]*x[4] - L[5][1]*x[5]);
	x[0] = invdiag[0]*(t[0] - L[1][0]*x[1] - L[2][0]*x[2] - L[3][0]*x[3] - L[4][0]*x[4] - L[5][0]*x[5]);
}
// Generated for postive-definite symmetric M, returns rank
template <class A1, class A2, class A3> inline int cholesky_compute(const FixedMatrix<6,6,A1>& M, FixedMatrix<6,6,A2>& L, FixedVector<6,A3>& invdiag)
{
	int rank = 6;
	double a;
	const double eps = std::numeric_limits<double>::min();
	a = M[0][0];
	if (eps < a) {
		const double factor = invdiag[0] = 1.0 / (L[0][0] = sqrt(a));
		L[0][1]=0;
		L[1][0] = factor * (M[0][1]);
		L[0][2]=0;
		L[2][0] = factor * (M[0][2]);
		L[0][3]=0;
		L[3][0] = factor * (M[0][3]);
		L[0][4]=0;
		L[4][0] = factor * (M[0][4]);
		L[0][5]=0;
		L[5][0] = factor * (M[0][5]);
	} else {
		--rank;
		L[0][0] = invdiag[0] = 0;
		L[0][1] = L[1][0] = 0;
		L[0][2] = L[2][0] = 0;
		L[0][3] = L[3][0] = 0;
		L[0][4] = L[4][0] = 0;
		L[0][5] = L[5][0] = 0;
	}
	a = M[1][1] - L[1][0]*L[1][0];
	if (eps < a) {
		const double factor = invdiag[1] = 1.0 / (L[1][1] = sqrt(a));
		L[1][2]=0;
		L[2][1] = factor * (M[1][2] - L[1][0]*L[2][0]);
		L[1][3]=0;
		L[3][1] = factor * (M[1][3] - L[1][0]*L[3][0]);
		L[1][4]=0;
		L[4][1] = factor * (M[1][4] - L[1][0]*L[4][0]);
		L[1][5]=0;
		L[5][1] = factor * (M[1][5] - L[1][0]*L[5][0]);
	} else {
		--rank;
		L[1][1] = invdiag[1] = 0;
		L[1][2] = L[2][1] = 0;
		L[1][3] = L[3][1] = 0;
		L[1][4] = L[4][1] = 0;
		L[1][5] = L[5][1] = 0;
	}
	a = M[2][2] - L[2][0]*L[2][0] - L[2][1]*L[2][1];
	if (eps < a) {
		const double factor = invdiag[2] = 1.0 / (L[2][2] = sqrt(a));
		L[2][3]=0;
		L[3][2] = factor * (M[2][3] - L[2][0]*L[3][0] - L[2][1]*L[3][1]);
		L[2][4]=0;
		L[4][2] = factor * (M[2][4] - L[2][0]*L[4][0] - L[2][1]*L[4][1]);
		L[2][5]=0;
		L[5][2] = factor * (M[2][5] - L[2][0]*L[5][0] - L[2][1]*L[5][1]);
	} else {
		--rank;
		L[2][2] = invdiag[2] = 0;
		L[2][3] = L[3][2] = 0;
		L[2][4] = L[4][2] = 0;
		L[2][5] = L[5][2] = 0;
	}
	a = M[3][3] - L[3][0]*L[3][0] - L[3][1]*L[3][1] - L[3][2]*L[3][2];
	if (eps < a) {
		const double factor = invdiag[3] = 1.0 / (L[3][3] = sqrt(a));
		L[3][4]=0;
		L[4][3] = factor * (M[3][4] - L[3][0]*L[4][0] - L[3][1]*L[4][1] - L[3][2]*L[4][2]);
		L[3][5]=0;
		L[5][3] = factor * (M[3][5] - L[3][0]*L[5][0] - L[3][1]*L[5][1] - L[3][2]*L[5][2]);
	} else {
		--rank;
		L[3][3] = invdiag[3] = 0;
		L[3][4] = L[4][3] = 0;
		L[3][5] = L[5][3] = 0;
	}
	a = M[4][4] - L[4][0]*L[4][0] - L[4][1]*L[4][1] - L[4][2]*L[4][2] - L[4][3]*L[4][3];
	if (eps < a) {
		const double factor = invdiag[4] = 1.0 / (L[4][4] = sqrt(a));
		L[4][5]=0;
		L[5][4] = factor * (M[4][5] - L[4][0]*L[5][0] - L[4][1]*L[5][1] - L[4][2]*L[5][2] - L[4][3]*L[5][3]);
	} else {
		--rank;
		L[4][4] = invdiag[4] = 0;
		L[4][5] = L[5][4] = 0;
	}
	a = M[5][5] - L[5][0]*L[5][0] - L[5][1]*L[5][1] - L[5][2]*L[5][2] - L[5][3]*L[5][3] - L[5][4]*L[5][4];
	if (eps < a) {
		invdiag[5] = 1.0 / (L[5][5] = sqrt(a));
	} else {
		--rank;
		L[5][5] = invdiag[5] = 0;
	}
	return rank;
}
// Generated for lower triangular L, inverse diagonal invdiag
template <class A1, class A2, class A3> inline void cholesky_inverse(const FixedMatrix<6,6,A1>& L, const FixedVector<6,A2>& invdiag, FixedMatrix<6,6,A3>& I)
{
	{ // column 0
		//forward substitution
		const double t0 = invdiag[0];
		const double t1 = -invdiag[1]*(L[1][0]*t0);
		const double t2 = -invdiag[2]*(L[2][0]*t0 + L[2][1]*t1);
		const double t3 = -invdiag[3]*(L[3][0]*t0 + L[3][1]*t1 + L[3][2]*t2);
		const double t4 = -invdiag[4]*(L[4][0]*t0 + L[4][1]*t1 + L[4][2]*t2 + L[4][3]*t3);
		const double t5 = -invdiag[5]*(L[5][0]*t0 + L[5][1]*t1 + L[5][2]*t2 + L[5][3]*t3 + L[5][4]*t4);
		//backward substitution
		const double x5 = invdiag[5]*(t5);
		const double x4 = invdiag[4]*(t4-L[5][4]*x5);
		const double x3 = invdiag[3]*(t3-L[4][3]*x4-L[5][3]*x5);
		const double x2 = invdiag[2]*(t2-L[3][2]*x3-L[4][2]*x4-L[5][2]*x5);
		const double x1 = invdiag[1]*(t1-L[2][1]*x2-L[3][1]*x3-L[4][1]*x4-L[5][1]*x5);
		const double x0 = invdiag[0]*(t0-L[1][0]*x1-L[2][0]*x2-L[3][0]*x3-L[4][0]*x4-L[5][0]*x5);
		I[5][0] = I[0][5] = x5;
		I[4][0] = I[0][4] = x4;
		I[3][0] = I[0][3] = x3;
		I[2][0] = I[0][2] = x2;
		I[1][0] = I[0][1] = x1;
		I[0][0] = x0;
	}
	{ // column 1
		//forward substitution
		const double t1 = invdiag[1];
		const double t2 = -invdiag[2]*(L[2][1]*t1);
		const double t3 = -invdiag[3]*(L[3][1]*t1 + L[3][2]*t2);
		const double t4 = -invdiag[4]*(L[4][1]*t1 + L[4][2]*t2 + L[4][3]*t3);
		const double t5 = -invdiag[5]*(L[5][1]*t1 + L[5][2]*t2 + L[5][3]*t3 + L[5][4]*t4);
		//backward substitution
		const double x5 = invdiag[5]*(t5);
		const double x4 = invdiag[4]*(t4-L[5][4]*x5);
		const double x3 = invdiag[3]*(t3-L[4][3]*x4-L[5][3]*x5);
		const double x2 = invdiag[2]*(t2-L[3][2]*x3-L[4][2]*x4-L[5][2]*x5);
		const double x1 = invdiag[1]*(t1-L[2][1]*x2-L[3][1]*x3-L[4][1]*x4-L[5][1]*x5);
		I[5][1] = I[1][5] = x5;
		I[4][1] = I[1][4] = x4;
		I[3][1] = I[1][3] = x3;
		I[2][1] = I[1][2] = x2;
		I[1][1] = x1;
	}
	{ // column 2
		//forward substitution
		const double t2 = invdiag[2];
		const double t3 = -invdiag[3]*(L[3][2]*t2);
		const double t4 = -invdiag[4]*(L[4][2]*t2 + L[4][3]*t3);
		const double t5 = -invdiag[5]*(L[5][2]*t2 + L[5][3]*t3 + L[5][4]*t4);
		//backward substitution
		const double x5 = invdiag[5]*(t5);
		const double x4 = invdiag[4]*(t4-L[5][4]*x5);
		const double x3 = invdiag[3]*(t3-L[4][3]*x4-L[5][3]*x5);
		const double x2 = invdiag[2]*(t2-L[3][2]*x3-L[4][2]*x4-L[5][2]*x5);
		I[5][2] = I[2][5] = x5;
		I[4][2] = I[2][4] = x4;
		I[3][2] = I[2][3] = x3;
		I[2][2] = x2;
	}
	{ // column 3
		//forward substitution
		const double t3 = invdiag[3];
		const double t4 = -invdiag[4]*(L[4][3]*t3);
		const double t5 = -invdiag[5]*(L[5][3]*t3 + L[5][4]*t4);
		//backward substitution
		const double x5 = invdiag[5]*(t5);
		const double x4 = invdiag[4]*(t4-L[5][4]*x5);
		const double x3 = invdiag[3]*(t3-L[4][3]*x4-L[5][3]*x5);
		I[5][3] = I[3][5] = x5;
		I[4][3] = I[3][4] = x4;
		I[3][3] = x3;
	}
	{ // column 4
		//forward substitution
		const double t4 = invdiag[4];
		const double t5 = -invdiag[5]*(L[5][4]*t4);
		//backward substitution
		const double x5 = invdiag[5]*(t5);
		const double x4 = invdiag[4]*(t4-L[5][4]*x5);
		I[5][4] = I[4][5] = x5;
		I[4][4] = x4;
	}
	{ // column 5
		//forward substitution
		const double t5 = invdiag[5];
		//backward substitution
		const double x5 = invdiag[5]*(t5);
		I[5][5] = x5;
	}
}
