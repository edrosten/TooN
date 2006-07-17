
namespace VectorMagic {
  template <int N=1>
  struct ComponentPlaceHolder {
  };

  struct InsertionStyle {};
  struct CommaStyle {};


  template<bool Good> struct sentinel
  {
	typedef int too_many_elements_inserted;
  };

  template<> struct sentinel<false>
  {
  };

  template <int Index, int Limit, class Vec, class Style> struct VectorFiller;
  template <int Index, int Limit, class Vec> struct VectorFiller<Index,Limit,Vec,CommaStyle> {


    // a static assertion
    typedef typename sentinel<(Limit >= Index)>::too_many_elements_inserted dummy; 
  
    Vec& v;
    bool final_initializer_but_Vector_not_filled;
    inline VectorFiller(Vec& vec) : v(vec), final_initializer_but_Vector_not_filled(Index!=Limit) {}

    inline VectorFiller<Index+1,Limit,Vec,CommaStyle> operator,(double t) {
      v[Index] = t;
      final_initializer_but_Vector_not_filled = false;
      return VectorFiller<Index+1,Limit,Vec,CommaStyle>(v);
    }
    template <int N> inline VectorFiller<Index+N,Limit,Vec,CommaStyle> operator,(const ComponentPlaceHolder<N>& ph) {
      final_initializer_but_Vector_not_filled = false;
      return (VectorFiller<Index+1,Limit,Vec,CommaStyle>(v), ComponentPlaceHolder<N-1>());
    }
    inline VectorFiller<Index+1,Limit,Vec,CommaStyle> operator,(const ComponentPlaceHolder<1>& ph) {
      final_initializer_but_Vector_not_filled = false;
      return VectorFiller<Index+1,Limit,Vec,CommaStyle>(v);
    }
    inline ~VectorFiller() {
      assert(!final_initializer_but_Vector_not_filled);
    }
    inline operator Vec () const { return v; }
  };
  
  template <int Index, int Limit, class Vec> struct VectorFiller<Index,Limit,Vec,InsertionStyle> {
	typedef typename sentinel<(Limit >= Index)>::too_many_elements_inserted dummy; 
    Vec& v;
    inline VectorFiller(Vec& vec) : v(vec){}

    template <int N> inline VectorFiller<Index+N,Limit,Vec,InsertionStyle> operator<<(const Vector<N>& t) {
      v.template slice<Index,N>() = t;
      return VectorFiller<Index+N,Limit,Vec,InsertionStyle>(v);
    }

    inline VectorFiller<Index+1,Limit,Vec,InsertionStyle> operator<<(double t) {
      v[Index] = t;
      return VectorFiller<Index+1,Limit,Vec,InsertionStyle>(v);
    }
    template <int N> inline VectorFiller<Index+N,Limit,Vec,InsertionStyle> operator<<(const ComponentPlaceHolder<N>& ph) {
      return (VectorFiller<Index+1,Limit,Vec,InsertionStyle>(v) << ComponentPlaceHolder<N-1>());
    }
    inline VectorFiller<Index+1,Limit,Vec,InsertionStyle> operator<<(const ComponentPlaceHolder<1>& ph) {
      return VectorFiller<Index+1,Limit,Vec,InsertionStyle>(v);
    }
    inline operator Vec () const { return v; }
  };


  template <class Left, int Size> struct VectorCreator
  {
    const Left& left;
    double val;
    VectorCreator(const Left& l, double v) : left(l), val(v) { }
    VectorCreator<VectorCreator<Left,Size>, Size+1> operator,(double t) const {
      return VectorCreator<VectorCreator<Left,Size>, Size+1>(*this, t);
    }
    template <class V> void assign(V& v) const { v[Size-1] = val; left.assign(v); }
    operator Vector<Size> () const {
      Vector<Size> v;
      assign(v);
      return v;
    }
  };

  struct BaseVectorCreator
  {
    inline VectorCreator<BaseVectorCreator, 1> operator,(double t) const {
      return VectorCreator<BaseVectorCreator, 1>(*this, t);
    }
    template <class V> inline void assign(V& ) const {}
  };
}

static VectorMagic::BaseVectorCreator make_Vector;

namespace VectorMagic 
{
  inline void dummy_make_Vector_user() { int i; make_Vector.assign(i); }
}

template <int N> VectorMagic::ComponentPlaceHolder<N> no_change() { return VectorMagic::ComponentPlaceHolder<N>(); }
inline VectorMagic::ComponentPlaceHolder<1> no_change() { return VectorMagic::ComponentPlaceHolder<1>(); }
