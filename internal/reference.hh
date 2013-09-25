// -*- c++ -*-

// Copyright (C) 2009 Tom Drummond (twd20@cam.ac.uk),
// Ed Rosten (er258@cam.ac.uk)

//All rights reserved.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions
//are met:
//1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//2. Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
//THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND OTHER CONTRIBUTORS ``AS IS''
//AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR OTHER CONTRIBUTORS BE
//LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//POSSIBILITY OF SUCH DAMAGE.

namespace TooN {

////////////////////////////////////////////////////////////////////////////////
//
// Helper classes for matrices constructed as references to foreign data
//

struct Reference
{

	template<int Size, typename Precision>
	struct VLayout
		: public Internal::GenericVBase<Size, Precision, 1, Internal::VectorSlice<Size, Precision> >
	{

		VLayout(Precision* p, int sz=0)
			: Internal::GenericVBase<Size, Precision, 1, Internal::VectorSlice<Size, Precision> >(p, sz, 0)
		{}
	};


	struct RowMajor
	{
		template<int Rows, int Cols, class Precision>
		struct MLayout
			: public Internal::GenericMBase<Rows, Cols, Precision, (Cols==-1?-2:Cols), 1, Internal::MatrixSlice<Rows, Cols, Precision> >
		{

			MLayout(Precision* p)
				: Internal::GenericMBase<Rows,Cols,Precision, (Cols==-1?-2:Cols), 1, Internal::MatrixSlice<Rows, Cols, Precision> > (p)
			{}
			MLayout(Precision* p, int r, int c)
				: Internal::GenericMBase<Rows,Cols,Precision, (Cols==-1?-2:Cols), 1, Internal::MatrixSlice<Rows, Cols, Precision> > (p, r, c, 0, 0)
			{}
		};
	};

	struct ColMajor
	{
		template<int Rows, int Cols, class Precision> struct MLayout: public Internal::GenericMBase<Rows, Cols, Precision, 1, (Rows==-1?-2:Rows), Internal::MatrixSlice<Rows, Cols, Precision> >
		{
			MLayout(Precision* p)
				: Internal::GenericMBase<Rows, Cols, Precision, 1, (Rows==-1?-2:Rows), Internal::MatrixSlice<Rows, Cols, Precision> >(p)
			{}
			MLayout(Precision* p, int r, int c)
				: Internal::GenericMBase<Rows, Cols, Precision, 1, (Rows==-1?-2:Rows), Internal::MatrixSlice<Rows, Cols, Precision> >(p, r, c, 0, 0)
			{}
		};
	};
};


                                    inline       Vector<Dynamic, double,          Reference> wrapVector(double* data, int size)          { return Vector<Dynamic, double,          Reference>(data, size); }
                                    inline const Vector<Dynamic, const double,    Reference> wrapVector(const double* data, int size)    { return Vector<Dynamic, const double,    Reference>(data, size); }
template<int Size>                  inline       Vector<Size,    double,          Reference> wrapVector(double* data)                    { return Vector<Size,    double,          Reference>(data); }
template<int Size>                  inline const Vector<Size,    const double,    Reference> wrapVector(const double* data)              { return Vector<Size,    const double,    Reference>(data); }
template<class Precision>           inline       Vector<Dynamic, Precision,       Reference> wrapVector(Precision* data, int size)       { return Vector<Dynamic, Precision,       Reference>(data, size); }
template<class Precision>           inline const Vector<Dynamic, const Precision, Reference> wrapVector(const Precision* data, int size) { return Vector<Dynamic, const Precision, Reference>(data, size); }
template<int Size, class Precision> inline       Vector<Size,    Precision,       Reference> wrapVector(Precision* data)                 { return Vector<Size,    Precision,       Reference>(data); }
template<int Size, class Precision> inline const Vector<Size,    const Precision, Reference> wrapVector(const Precision* data)           { return Vector<Size,    const Precision, Reference>(data); }

///Wrap external data as a \link TooN::Matrix Matrix \endlink
///As usual, if template sizes are provided, then the run-time size is only
///used if the template size is not Dynamic.
///@ingroup gLinAlg
///@{
//Fully static matrices, ie no size parameters
template<int Rows, int Cols>                  inline       Matrix<Rows, Cols,       double,          Reference::RowMajor> wrapMatrix(double*    data)                           { return Matrix<Rows, Cols,       double,          Reference::RowMajor>(data);}
template<int Rows, int Cols>                  inline const Matrix<Rows, Cols,       const double,    Reference::RowMajor> wrapMatrix(const double*    data)                     { return Matrix<Rows, Cols,       const double,    Reference::RowMajor>(data);}
template<int Rows, int Cols, class Precision> inline       Matrix<Rows, Cols,       Precision,       Reference::RowMajor> wrapMatrix(Precision* data)                           { return Matrix<Rows, Cols,       Precision,       Reference::RowMajor>(data);}
template<int Rows, int Cols, class Precision> inline const Matrix<Rows, Cols,       const Precision, Reference::RowMajor> wrapMatrix(const Precision* data)                     { return Matrix<Rows, Cols,       const Precision, Reference::RowMajor>(data);}
//Static sizes with size parameters (useful for half-dynamic matrices)
template<int Rows, int Cols>                  inline       Matrix<Rows, Cols,       double,          Reference::RowMajor> wrapMatrix(double*    data, int rows, int cols)       { return Matrix<Rows, Cols,       double,          Reference::RowMajor>(data, rows, cols);}
template<int Rows, int Cols>                  inline const Matrix<Rows, Cols,       const double,    Reference::RowMajor> wrapMatrix(const double*    data, int rows, int cols) { return Matrix<Rows, Cols,       const double,    Reference::RowMajor>(data, rows, cols);}
template<int Rows, int Cols, class Precision> inline       Matrix<Rows, Cols,       Precision,       Reference::RowMajor> wrapMatrix(Precision* data, int rows, int cols)       { return Matrix<Rows, Cols,       Precision,       Reference::RowMajor>(data, rows, cols);}
template<int Rows, int Cols, class Precision> inline const Matrix<Rows, Cols,       const Precision, Reference::RowMajor> wrapMatrix(const Precision* data, int rows, int cols) { return Matrix<Rows, Cols,       const Precision, Reference::RowMajor>(data, rows, cols);}
//Fully dynamic
                                              inline       Matrix<Dynamic, Dynamic, double,          Reference::RowMajor> wrapMatrix(double*          data, int rows, int cols) { return Matrix<Dynamic, Dynamic, double,          Reference::RowMajor>(data, rows, cols);}
                                              inline const Matrix<Dynamic, Dynamic, const double,    Reference::RowMajor> wrapMatrix(const double*    data, int rows, int cols) { return Matrix<Dynamic, Dynamic, const double,    Reference::RowMajor>(data, rows, cols);}
template<class Precision>                     inline       Matrix<Dynamic, Dynamic, Precision,       Reference::RowMajor> wrapMatrix(Precision* data, int rows, int cols)       { return Matrix<Dynamic, Dynamic, Precision,       Reference::RowMajor>(data, rows, cols);}
template<class Precision>                     inline const Matrix<Dynamic, Dynamic, const Precision, Reference::RowMajor> wrapMatrix(const Precision* data, int rows, int cols) { return Matrix<Dynamic, Dynamic, const Precision, Reference::RowMajor>(data, rows, cols);}
///@}
}
