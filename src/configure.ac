dnl :    edit this file and run 'autoconf' to generate a new 'configure' file

AC_INIT(linear_wave.c)
AC_PROG_CC
AC_PROG_INSTALL
AC_OUTPUT(Makefile athena.def)



dnl --- set DIM ---

AC_ARG_WITH(dim, [ --with-dim=DIM     Dimensionality of program],
		with_dim=$withval, with_dim=1)
	

AC_SUBST(XXX_D)
if test "$with_dim" == "1"; then
  XXX_D="ONE_D"
elif test "$with_dim" == "2"; then
  XXX_D="TWO_D"
elif test "$with_dim" == "3"; then
  XXX_D="THREE_D"
fi
	
