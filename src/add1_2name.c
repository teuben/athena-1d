/*==============================================================================
 *                        Function ADD1_2NAME 
 * Adds one to numerical counter embedded in ATHENA filenames.  Filename is a
 * character string of length 9 with format XXXNNNXX\n, where NNN=counter
 *============================================================================*/
void add1_2name(char file[9])
{
  if (*(file+5) < '9') {
    ++*(file+5);
  }
  else if (*(file+4) < '9') {
    ++*(file+4);
    *(file+5)='0';
  }
  else if (*(file+3) < '9') {
    ++*(file+3);
    *(file+4)='0';
    *(file+5)='0';
  }
  else {
    *(file+3)='X';
    *(file+4)='X';
    *(file+5)='X';
  }
}
