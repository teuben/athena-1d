dnl :    edit this file and run 'autoconf' to generate a new 'configure' file

AC_INIT(linear_wave.c)
AC_PROG_CC
AC_PROG_INSTALL


dnl --- set DIM ---

AC_SUBST(DIM)
AC_SUBST(XXX_D)

AC_ARG_WITH(dim, [ --with-dim=DIM     Dimensionality of program],
		with_dim=$withval, with_dim=1)
	
DIM=$withval
if test "$with_dim" == "1"; then
  XXX_D="ONE_D"
elif test "$with_dim" == "2"; then
  XXX_D="TWO_D"
elif test "$with_dim" == "3"; then
  XXX_D="THREE_D"
else
  XXX_D="INVALID_D"
fi
echo Dimensionality of ATHENA = $XXX_D	


dnl --- write final files --

AC_CONFIG_FILES(Makefile athena.def)
AC_OUTPUT

