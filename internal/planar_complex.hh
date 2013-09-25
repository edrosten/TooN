//Copyright (C) Edward Rosten 2010

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

namespace TooN{
	namespace Internal{
		
		///@internal
		///@brief A simple smart pointer type representing planar complex data.
		///Not resizable. Also, only the minimum number of smart pointer
		///functionality has been implemented to make Vector work. The class returns
		///reference types, so it can represent mutable data.
		///@ingroup gInternal
		template<class Precision> struct PointerToPlanarComplex;
		template<class Precision> struct PointerToPlanarComplex<std::complex<Precision> >
		{
			const Precision* const re;
			const Precision* const im;
		
			PointerToPlanarComplex(std::pair<Precision*, Precision*> d)
			:re(d.first),im(d.second)
			{}

			PointerToPlanarComplex(std::pair<const Precision*, const Precision*> d)
			:re(d.first),im(d.second)
			{}

			PointerToPlanarComplex<std::complex<Precision> > operator+(int i) const
			{
				return PointerToPlanarComplex<std::complex<Precision> >(std::make_pair(re+i, im+i));
			}

			const std::complex<Precision> operator[](int i) const
			{
				return std::complex<Precision>(re[i], im[i]);
			}

		};
	}

	
	struct ReferencePlanarComplex
	{
		template<int Size, typename Precision>
		struct VLayout;

		template<int Size, typename Precision>
		struct VLayout<Size, std::complex<Precision> >: 
		  public Internal::GenericVBase<Size, std::complex<Precision>, 1, Internal::VectorSlice<Size, std::complex<Precision>, Internal::PointerToPlanarComplex<std::complex<double> >,
	                                                    Internal::PointerToPlanarComplex<std::complex<double> >,
	                                                    const std::complex<double>,
	                                                    const std::complex<double> > >
		{
			VLayout(Internal::PointerToPlanarComplex<std::complex<Precision> > p, int sz=0)
			:Internal::GenericVBase<Size, std::complex<Precision>, 1, Internal::VectorSlice<Size, std::complex<Precision>, Internal::PointerToPlanarComplex<std::complex<double> >,
	                                                    Internal::PointerToPlanarComplex<std::complex<double> >,
	                                                    const std::complex<double>,
	                                                    const std::complex<double> > >(p, sz, 1)
			{}

			
		};

	};
}

