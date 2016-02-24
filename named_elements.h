#ifndef TOON_INCLUDE_NAMED_ELEMENTS_H
#define TOON_INCLUDE_NAMED_ELEMENTS_H

#include <TooN/TooN.h>
using namespace TooN;
using namespace std;

namespace TooN
{
	namespace Internal
	{
		//All the functionality we need for a statically allocated base class
		//i.e. not much.
		template<int Size, class Precision> class StaticAllocFunctions
		{
			public:
				StaticAllocFunctions() = default;
				StaticAllocFunctions(int) {};
				StaticAllocFunctions(const StaticAllocFunctions&)=default;

				///Construction from an Operator. See Operator::size().
				template<class Op>
					StaticAllocFunctions(const Operator<Op>&) {}

				using PointerType = Precision*;
				using ConstPointerType = const Precision*;
				using ReferenceType = Precision&;
				using ConstReferenceType = const Precision&;

				///Return the size of the vector.
				int size() const {
					return Size;
				}

			protected:
				void try_destructive_resize(int)
				{}

				template<class Op> void try_destructive_resize(const Operator<Op>&) 
				{}
		};



		//Vector base class to use a supplied templated allocator.
		template<template<int,typename> class Alloc>
			struct StaticBase
			{
				// this class is really just a typedef
				template<int Size, class Precision>
					struct VLayout : public Internal::GenericVBase<Size, Precision, 1, Alloc<Size, Precision> > {

						VLayout(){}
						static_assert(sizeof(Alloc<Size, Precision>) == Size * sizeof(Precision), "Weird Padding");

						VLayout(VLayout&&) = default;
						VLayout(const VLayout&) = default;

						VLayout(int s)
							:Internal::GenericVBase<Size, Precision, 1, Alloc<Size, Precision> >(s)
						{}

						template<class Op>
							VLayout(const Operator<Op>& op)
							:Internal::GenericVBase<Size, Precision, 1, Alloc<Size, Precision> >(op) {}
					};
			};
	}
}

//Now some variadic macros for generating the correct classes.
#define TOON_MAKE_NAMED_ELEMENT_VECTOR(Name, ...)\
namespace TooN{namespace Internal{\
	struct Name##Dummy\
	{\
		char __VA_ARGS__;\
	};\
	template<int S, class Precision> struct Name##Alloc\
	{\
		static_assert(S == sizeof(Name##Dummy), #Name " Is the wrong size.");\
	};\
	\
	template<class Precision> struct Name##Alloc<sizeof(Name##Dummy), Precision>: public StaticAllocFunctions<sizeof(Name##Dummy), Precision>\
	{\
		using StaticAllocFunctions<sizeof(Name##Dummy), Precision>::StaticAllocFunctions;\
		union\
		{\
			struct\
			{\
				Precision __VA_ARGS__;\
			};\
			Precision my_data[sizeof(Name##Dummy)];\
		};\
		Precision *data()\
		{\
			return my_data;\
		}\
		const Precision *data() const\
		{\
			return my_data;\
		}\
	};\
}}\
\
template<class Precision=DefaultPrecision> using Name = Vector<sizeof(TooN::Internal::Name##Dummy), Precision, TooN::Internal::StaticBase<Internal::Name##Alloc>>;

#endif
