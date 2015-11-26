AC_INIT(TooN, version-3.0.0-beta1)
AC_COPYRIGHT([Copyright T. Drummond, E. Rosten, G. Reitmayr 2006--2015])
AC_LANG(C++)
AC_PROG_CXX


m4_include([m4/ax_cxx_compile_stdcxx_14.m4])


################################################################################
#
# Useful macros
#

dnl APPEND(var, value)
dnl This appends vale to a shell variable var
define(APPEND, [$1="$$1 $2"])



dnl TEST_AND_SET_CXXFLAG(flag, [program], [run])
dnl
dnl This attempts to compile a and run program with a certain compiler flag.
dnl If no program is given, then the minimal C++ program is compiled, and 
dnl this tests just the validity of the compiler flag. 
dnl
define([TEST_AND_SET_CXXFLAG],[
	if test "$3" == ""
	then
		AC_MSG_CHECKING([if compiler flag $1 works])	
	else
		AC_MSG_CHECKING([$3])	
	fi
	save_CXXFLAGS="$CXXFLAGS"
	APPEND(CXXFLAGS, [$1])

	m4_if([$2],[],[pushdef(prog, [int main(){}])], [pushdef(prog, [$2])])

	m4_if([$4],[run],
		[AC_RUN_IFELSE([AC_LANG_SOURCE([prog])], [cvd_conf_test=1],[cvd_conf_test=0], [cvd_conf_test=0])],
		[AC_COMPILE_IFELSE([AC_LANG_SOURCE([prog])], [cvd_conf_test=1],[cvd_conf_test=0])]
	)
	
	
	popdef([prog])

	if test $cvd_conf_test = 1
	then
		AC_MSG_RESULT(yes)
		ts_success=yes
	else
		AC_MSG_RESULT(no)
		CXXFLAGS="$save_CXXFLAGS"
		ts_success=no
	fi
])






AX_CXX_COMPILE_STDCXX_14


#
# Test for Package Config
#
AC_PATH_PROG(PKG_CONFIG, pkg-config)
if test "x$PKG_CONFIG" = "x"; then
	AC_MSG_WARN([Could not find pkg-config, will not create pc file.])
else
	# we need sed to find the pkg-config lib directory
	AC_CHECK_PROG(SED,sed,sed,AC_MSG_ERROR([You Must install sed]))
	AC_MSG_CHECKING([for pkg-config library dir])
	PKGCONFIG_LIBDIR="`echo $PKG_CONFIG | $SED -e 's~.*/bin/pkg-config$~~'`${libdir}/pkgconfig"

	AC_MSG_RESULT($PKGCONFIG_LIBDIR)
	AC_SUBST(PKGCONFIG_LIBDIR)
	AC_CONFIG_FILES([ TooN.pc ])

	# This will be put into the pc file
	VERSION=AC_PACKAGE_VERSION()
	AC_SUBST(VERSION)
fi


# See if we can find LAPACK
#ACX_LAPACK(
#			 LAPACK_LINKER_FLAGS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
#			 AC_SUBST(LAPACK_LINKER_FLAGS)
#		   , AC_ERROR([Could not find LAPACK])
#		   )
#

typeof=check
AC_ARG_ENABLE(lapack, [AS_HELP_STRING([--enable-lapack],[Use LAPACK where optional])], [lapack=$enableval])
AC_ARG_WITH(default_precision, [AS_HELP_STRING([--with-default_precision=X],[Override default precision from double to X])], [default_precision="$withval"])

if test "$default_precision" != ""
then
	AC_DEFINE_UNQUOTED(TOON_DEFAULT_PRECISION, $default_precision)
fi

if test "$lapack" == "" 
then
	AC_CHECK_LIB(lapack, dgesvd_)
	if test $ac_cv_lib_lapack_dgesvd_ = yes
	then
		lapack=yes
	else
		lapack=no
	fi
fi

if test "$lapack" != "no" 
then
	AC_DEFINE(TOON_USE_LAPACK, 1)
	AC_SUBST(use_lapack, yes)
fi



TEST_AND_SET_CXXFLAG(-Wall)
TEST_AND_SET_CXXFLAG(-Wextra)
TEST_AND_SET_CXXFLAG(-W)
TEST_AND_SET_CXXFLAG(-O3)
TEST_AND_SET_CXXFLAG(-Werror)

AC_CONFIG_HEADERS(internal/config.hh)
AC_OUTPUT(Makefile)

