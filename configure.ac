AC_INIT(TooN, version-1.0.1)

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

AC_OUTPUT(Makefile)
