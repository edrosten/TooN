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
