AC_INIT(TooN, version-1.0.1)
AC_COPYRIGHT([Copyright T. Drummond, E. Rosten 2006, 2007, 2008])

toon2=no

AC_ARG_ENABLE(toon2, [AS_HELP_STRING([--enable-toon2],[Enable TooN 2.0. DO NOT USE!])], [toon2=yes])

if test "$toon2" != yes
then
cat<<FOO
You're on the completely broken development branch of TooN 2.0, and it WILL NOT
WORK. If you want a working TooN version, then run, get it using:

cvs -z3 -d:pserver:anoncvs@cvs.savannah.nongnu.org:/cvsroot/toon co -r Maintenance_Branch_1_x TooN

or update what you currently have using:

cvs up -r Maintenance_Branch_1_x

or head over to:

http://mi.eng.cam.ac.uk/~er258/cvd/

and get the latest snapshot. If you are trying to hack the latest version of
TooN in to a working form, then run:

./configure --enable-toon2

FOO

AC_MSG_ERROR([Please don't use this version of TooN])

fi

AC_LANG(C++)
AC_PROG_CXX

typeof=check
AC_ARG_ENABLE(typeof, [AS_HELP_STRING([--enable-typeof=X],[Force typeof to be decltype, __typeof__, typeof, boost or none])], [typeof=$enableval])

define(TYPEOF_TEST,[
	if test "$typeof" == check
	then
		AC_MSG_CHECKING([for $1])
		AC_COMPILE_IFELSE([
		$3
		int main()
		{
			typedef $2(int()*double()) foo;
		}
		], [ AC_MSG_RESULT(yes);typeof=$1], AC_MSG_RESULT(no))
	fi
])

TYPEOF_TEST(decltype, decltype)
TYPEOF_TEST(typeof, typeof)
TYPEOF_TEST(__typeof__, __typeof__)
TYPEOF_TEST(boost, BOOST_TYPEOF, [#include<boost/typeof/typeof.hpp>])

if test "$typeof" == typeof
then
	AC_DEFINE(TOON_TYPEOF_TYPEOF, 1)
elif test "$typeof" == __typeof__
then
	AC_DEFINE(TOON_TYPEOF___TYPEOF__, 1)
elif test "$typeof" == decltype
then
	AC_DEFINE(TOON_TYPEOF_DECLTYPE, 1)
elif test "$typeof" == boost
then
	AC_DEFINE(TOON_TYPEOF_BOOST, 1)
fi


AC_CONFIG_HEADERS(internal/config.hh)
AC_OUTPUT(Makefile)

