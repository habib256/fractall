/*
 *     Complexmath.c
 *     Complex functions
 *     based on Akihiro Sato work (27/Oct/1994)
 * 
 *     Rewritten by VERHILLE A. (GPL 2002)
 *     Modif: A new coding style and a better Exception handler
 */

#include <math.h>
#include "complexmath.h"


/*  Set of Complex functions
 *****************************/

complex MakeComplex(double real, double imag)  /* rz = real + i * imag */
{
  complex rz;
  rz.x = real;
  rz.y = imag;
  return rz;
}  /* end. MakeComplex */


complex ZeroSetofComplex()  /* rz = 0.0 + i * 0.0 */
{
  complex rz;
  rz.x = 0;
  rz.y = 0;
  return rz;
}  /* end. ZeroSetofComplex */


complex ESetofComplex()  /* rz = 1.0 + i * 0.0 */
{
  complex rz;
  rz.x = 1.0;
  rz.y = 0;
  return rz;
}  /* end. ESetofComplex */


complex ISetofComplex()  /* rz = 0.0 + i * 1.0 */
{
  complex rz;
  rz.x = 0;
  rz.y = 1.0;
  return rz;
}  /* end. ISetofComplex */


complex InfSetofComplex()  /* rz = Infinity */
{
  complex rz;
  rz.x = INFINITY;
  rz.y = 0.0;
  return rz;
}  /* end. InfSetofComplex */

complex NaNSetofComplex()  /* rz = NaN + i * NaN */
{
  complex rz;
  rz.x = NOT_A_NUMBER;
  rz.y = NOT_A_NUMBER;
  return rz;
}  /* end. NaNSetofComplex */


complex SetofComplex(complex z)  /* rz = z */
{
  return z;
}  /* end. SetofComplex */


complex ScalarTimesofComplex(double c, complex z)  /* rz = c * z */
{
  complex rz;
  rz.x = c * z.x;
  rz.y = c * z.y;
  return rz;
}  /* end. ScalarTimesofComplex */


/* Neg and Read Complex functions
************************************/ 

complex Negz(complex z)   /* rz = -z */
{
  complex rz;
  if( z.x != 0.0) rz.x = -z.x;
  if( z.y != 0.0) rz.y = -z.y;
  return rz;
}   /* end. Negz */ 


double Rez(complex z) 
{
  return(z.x);
}   /* end. Rez */


double Imz(complex z)
{
  return(z.y);
}   /* end. Imz */


double Magz(complex z)  /* Magz(z)=|rz| */
{
  
  double abx,aby,temp;
  
  if( z.x == INFINITY || z.y == INFINITY) return(INFINITY);

  abx = fabs( z.x); aby = fabs( z.y);

  if( z.x == 0.0 && z.y == 0.0) return( 0.0); 
  if( z.x == 0.0 && z.y !=0.0) return( aby); 
  if( z.x != 0.0 && z.y == 0.0) return( abx); 
  if( abx!= 0.0 && aby != 0.0){ 
    if( abx >= aby){ 
      temp = aby / abx;
      return( abx*sqrt( 1.0 + temp*temp));
    } 
    else{ 
      temp = abx / aby; 
      return( aby*sqrt( temp*temp + 1.0));
    } 
  } 
  return NOT_A_NUMBER;
} /* end. Magz */


double Magz2(complex z)
{
  
  double abx,aby,temp; 
  if( z.x == INFINITY || z.y == INFINITY) return( INFINITY);

  abx = fabs( z.x); aby = fabs( z.y);

  if( abx == 0.0 && aby == 0.0) return( 0.0);
  if( abx == 0.0 && aby != 0.0) return( z.y*z.y); 
  if( abx != 0.0 && aby == 0.0) return( z.x*z.x);

  if( abx != 0.0 && aby != 0.0){ 
    if( abx >= aby){ 
      temp = z.y / z.x;
      return( z.x*z.x*( 1.0 + temp*temp));
    } 
    else{ 
      temp = z.x / z.y; 
      return( z.y*z.y*( temp*temp + 1.0));
    } 
  }
  return NOT_A_NUMBER;
} /* end. Magz2 */

double Argz(complex z) 
{ 
  return( atan2(z.y, z.x));
} /* end. Argz */


complex Addz(complex z1, complex z2) /* rz = z1 + z2 */ 
{ 
  complex rz;
  rz.x = z1.x + z2.x;
  rz.y = z1.y + z2.y; 
  return rz;
} /* end. Addz */


complex Subz(complex z1, complex z2) /* rz = z1 - z2 */
{ 
  complex rz;
  rz.x =z1.x - z2.x; 
  rz.y = z1.y - z2.y;
  return rz;
} /* end. Subz */


complex Mulz(complex z1, complex z2) /* rz = z1 * z2 */  
{ 
  complex rz;
  double abz1,abz2;

  abz1 = Magz(z1); abz2 = Magz(z2);

  if( abz1 == INFINITY){ 
    if( abz2 == INFINITY) rz = InfSetofComplex(); /* Infinity * Infinity */ 
    if( abz2 == 0.0) rz = NaNSetofComplex(); /* Infinity * 0 */ 
  } 
  if( abz1 == 0.0){ 
    if( abz2 == INFINITY) rz = NaNSetofComplex(); /* 0 * Infinity */
    if( abz2 == 0.0) rz = ZeroSetofComplex(); /* 0 * 0 */ 
  } 
  else{ 
    rz.x = (z1.x) * (z2.x) - (z1.y) * (z2.y);
    rz.y = (z1.x) * (z2.y) + (z1.y) * (z2.x);
  }
  return rz;
} /*end. Mulz */


complex Divz(complex z1, complex z2) /* rz = z1 / z2 */  
{ 
  complex rz;
  double abx2,aby2,abz1,abz2,s,temp;

  abz1 = Magz(z1);  abz2 = Magz(z2);
  abx2 = fabs(z2.x);  aby2 = fabs(z2.y); 
  if( abz2 == 0.0){ 
    if( abz1 == 0.0 || abz1 == INFINITY){ 
      rz = NaNSetofComplex();
      THROW (EX_DIVZ);
      return rz;
    } /* in case of 0/0 or Inf/0 */
    else{ 
      rz = InfSetofComplex(); 
      return rz;
    } /* in case of z/0 */
  }

  if( abz2 == INFINITY){ 
    if( abz1 == 0.0 || abz1 == INFINITY){ 
      rz = NaNSetofComplex(); 
      THROW (EX_DIVZ);
      return rz;
    } /* in case of 0/Inf or Inf/Inf */
    else{ 
      rz = ZeroSetofComplex(); 
      return rz;
    } /* in case of 0/z */
  }

  if( z1.x != 0.0 && z1.y == 0.0 && z2.x != 0.0 && z2.y == 0.0){
    rz = MakeComplex( z1.x / z2.x, 0.0); 
    return rz;
  }

  if( z1.x == 0.0 && z1.y != 0.0 && z2.x == 0.0 && z2.y != 0.0){
    rz = MakeComplex( z1.y / z2.y, 0.0);
    return rz;
  }
  else{
    if( abx2 >= aby2){ 
      temp = z2.y / z2.x; 
      s = z2.x + z2.y * temp;
      rz.x = ( z1.x + z1.y * temp ) / s;
      rz.y = ( -z1.x * temp + z1.y ) / s;
    }
    else{
      temp = z2.x / z2.y;
      s = z2.x * temp + z2.y;
      rz.x = ( z1.x * temp + z1.y ) / s; 
      rz.y = ( -z1.x + z1.y * temp ) / s;
    } 
  } return rz;

} /* end. Divz */


double sign(double x) /* sign returns plus plus, minus or zero. */ 
{ 
  if( x == 0.0) return( 0.0);
  if( x>0 ) return( 1.0);
  if( x<0 ) return( -1.0);
  return (0.0);
} /* end. sign */


complex Sqrtz(complex z) /* rz = sqrt( z) */ 
{
  complex rz;
  double r,th,abz;
  
  abz = Magz( z);

  if( abz == INFINITY){
    rz = InfSetofComplex(); 
    return rz;
  }
  if( z.x == 0.0 && z.y == 0.0){
    rz = ZeroSetofComplex();
    return rz;
  }
 
  if( z.y == 0.0){
    if( z.x > 0.0) rz = MakeComplex( sqrt( z.x), 0.0);
    else rz = MakeComplex( 0.0, sqrt( fabs( z.x)));
  } else{
    r = Magz( z);
    th = Argz( z) / 2;

    rz.x = sign( cos( th)) * sqrt( ( r + z.x)/2);
    rz.y = sign( sin( th)) * sqrt( ( r - z.x)/2);
  }
  return rz;
}    /* end. Sqrtz */


complex powzn(complex z1,complex z2, int n)    /* rz = pow( z1, z2, n) */
{
  complex rz = NaNSetofComplex();
  double zz;
  double r,rr,lnr,arg;
  double th,e,sith,coth;

  zz = Magz( z1);

  if( zz == INFINITY){
    rz = InfSetofComplex();
    return rz;
  }

  if( zz == 0.0){
    if( z2.x == 0.0 && z2.y == 0.0){
      rz = NaNSetofComplex();
      THROW (EX_POWZN);
      return rz;
    }
    if( z2.y != 0.0 ){
      rz = NaNSetofComplex();
      THROW (EX_POWZN);
      return rz;
    }
    if( z2.x != 0.0 && z2.y == 0.0){
      rz = ZeroSetofComplex();
      return rz;
    }
  }
  else{
    if( z2.x == 0.0 && z2.y == 0.0){
      rz = ESetofComplex();
      return rz;
    } else{

      r = zz;
      rr = pow( r, z2.x);
      lnr = log( r);
      arg = Argz( z1) + 2 * n * M_PI;

      e = exp( - z2.y * arg );
      th = z2.y * lnr + z2.x * arg;

      coth = cos( th);
      sith = sin( th);

      if( rr * e == INFINITY) rz = InfSetofComplex();
      else{
	rz.x = rr * e *coth;
	rz.y = rr * e *sith;
       }
      return rz;
    }
  }
  return rz;
}  /* end. powzn */


complex Powz(complex z1, complex z2)   /* rz = pow( z1, z2) */
{
  return (powzn( z1, z2, 0));
}  /* end. powz */


complex sinz(complex z)    /* rz = sin( z) */
{
  complex rz;
  double ch, sh;

  ch = cosh( z.y);  sh = sinh( z.y);

  if( ch == INFINITY || sh == INFINITY) rz = InfSetofComplex();
  else{
    rz.x = sin( z.x) * cosh( z.y);
    rz.y = cos( z.x) * sinh( z.y);
  }
  return rz;
}  /* end. sinz */


complex cosz(complex z)     /* rz = cos( z) */
{
  complex rz;
  double ch,sh;

  ch = cosh( z.y); sh = sinh( z.y);

  if( ch == INFINITY || sh == INFINITY) rz = InfSetofComplex();
  else{
    rz.x = cos( z.x) * cosh( z.y);
    rz.y = -sin( z.x) * sinh( z.y);
  }
  return rz;
}  /* end. cosz */


complex tanz(complex z)     /* rz = tan( z) */ 
{
  complex rz;
  complex s,c;

  if( z.x == 0.0 && z.y == 0.0){
    rz = ZeroSetofComplex();
    return rz;
  }
  if( z.x != 0.0 && z.y == 0.0){
    rz = MakeComplex( tan( z.x), 0.0);
    return rz;
  }
  if( z.x == 0.0 && z.y != 0.0){
    rz = MakeComplex( 0.0, tanh( z.y));
    return rz;
  }
  if( z.x != 0.0 && z.y != 0.0){

    s = sinz( z);
    c = cosz( z);
    
    TRY {
    rz = Divz( s, c);
    return rz;
    } 
    CATCH (EX_DIVZ) {
    rz = NaNSetofComplex();
    THROW (EX_TANZ);
    return rz;
    }
    ENDTRY;
  }
} 

complex expz(complex z)     /* rz = exp( z) */ 
{
  complex rz;
  double ex;

  ex = exp( z.x);
  rz.x = ex*cos( z.y);
  rz.y = ex*sin( z.y);

  return rz;
}  /* end. expz */ 


complex sinhz(complex z)  /* rz = sinhz( z) */
{
  complex rz;
  double sh, ch;

  sh = sinh( z.x); ch = cosh( z.x);
  if( sh == INFINITY || ch == INFINITY) rz = InfSetofComplex();
  else{
    rz.x = sh * cos( z.y);
    rz.y = ch * sin( z.y);
  }
  return rz;
}  /* end. sinhz */


complex coshz(complex z)  /* rz = cosh( z) */
{
  complex rz;
  double sh, ch;

  sh = sinh( z.x); ch = cosh( z.x);
  if( sh == INFINITY || ch == INFINITY) rz = InfSetofComplex();
  else{
    rz.x = ch * cos( z.y);
    rz.y = sh * sin( z.y);
  }
  return rz;
}  /* end. coshz */


complex tanhz(complex z)  /* rz = tanh( z) */
{
  complex rz;
  complex sh,ch;

  if( z.x == 0.0 && z.y == 0.0){
    rz = ZeroSetofComplex();
    return rz;
  }
  if( z.x != 0.0 && z.y == 0.0){
    rz = MakeComplex( tanh( z.x), 0.0);
    return rz;
  }
  if( z.x == 0.0 && z.y != 0.0){
    rz = MakeComplex( 0.0, tan( z.y));
    return rz;
  }
  if( z.x != 0.0 && z.y != 0.0){
 
    sh = sinhz( z);
    ch = coshz( z);

    TRY {
    rz = Divz( sh, ch);
    return rz;
    }
    CATCH (EX_DIVZ) {
      rz = NaNSetofComplex();
      THROW (EX_TANHZ);
      return rz;
    }
    ENDTRY;
  }
  return rz;
}  /* end. tanhz */

complex Logz(complex z)  /* rz = log( z) */
{
  complex rz;
  double zz,th;

  zz = Magz(z);

  if( zz == 0.0){
    rz = NaNSetofComplex();
    THROW (EX_LOGZ);
    return rz;
  }
  else{
    rz.x = log( zz);
    th = Argz( z);
    rz.y = th;
    return rz;
  }

}  /* end. Logz */ 


complex Arcsinz(complex z) /* rz=Arcsin(z)=-i*Log{i*z+sqrt(1-z*z)} */
{
  complex rz;
  complex tempz1,tempz2;
  complex e,i;

  e = ESetofComplex();
  i = ISetofComplex();
  
  tempz1 = Mulz( z, z);
  tempz1 = Subz( e, tempz1);
  tempz1 = Sqrtz( tempz1);

  tempz2 = Mulz( i, z);
  tempz1 = Addz( tempz1, tempz2);
  
  TRY {
  tempz1 = Logz( tempz1);
  tempz1 = Mulz( i, tempz1);
  rz =  Negz( tempz1);
  return rz;
  }
  CATCH (EX_LOGZ) {
    rz = NaNSetofComplex ();
    THROW (EX_ARCSINZ);
    return rz;
  }
  ENDTRY;
}  

complex Arccosz(complex z)  /* rz=Arccos(z)=-i*Log{z+sqrt(z*z-1)} */  
{
  complex rz;
  complex temp;
  complex e,i;

  e = ESetofComplex();
  i = ISetofComplex();

  temp = Mulz( z, z);
  temp = Subz( temp, e);
  temp = Sqrtz( temp);

  temp = Addz( z, temp);
  
  TRY {
  temp = Logz( temp);
  temp = Mulz( i, temp); 
  rz = Negz( temp);
  return rz;
  }
  CATCH (EX_LOGZ) {
    rz = NaNSetofComplex ();
    THROW (EX_ARCCOSZ);
    return rz;
  }
  ENDTRY;
}   


complex Arctanz(complex z)   /* rz=Arctan(z)=i/2*Log{(1-i*z)/(1+i*z)} */
{
  complex rz;
  complex temp0,temp1,temp2,temp3;
 
  complex e,i,halfi;

  if( z.x == 0.0 && z.y == 0.0){
    rz = ZeroSetofComplex();
    return rz;
  }
  if( z.x != 0.0 && z.y == 0.0){
    rz = MakeComplex( atan( z.x), 0.0);
    return rz;
  }
  if( z.x == 0.0 && z.y != 0.0){

    temp0.x = Imz( z);  temp0.y = 0.0;
    i = ISetofComplex();
    temp0 = Arctanz( temp0);

    rz = Mulz( i, temp0);
    return rz;
  } 
  if( z.x != 0.0 && z.y != 0.0){
    
    e = ESetofComplex();
    i = ISetofComplex();
    halfi = MakeComplex( 0.0, 0.5);

    temp0 = Mulz( i, z );
    temp1 = Subz( e, temp0);
    temp2 = Addz( e, temp0);

    TRY {
      temp3 = Divz (temp1, temp2);
      temp3 = Logz (temp3);
      rz = Mulz (halfi, temp3);
    }
    CATCH (EX_DIVZ) {
      rz = NaNSetofComplex();
      THROW (EX_ARCTANZ);
      return rz;
    }
    CATCH (EX_LOGZ) {
      rz = NaNSetofComplex();
      THROW (EX_ARCTANZ);
      return rz;
    }
    ENDTRY;
  }
  return rz;
}   /* end. Arctanz */


complex Arcsinhz( complex z)    /* rz=Arcsinh(z)=Log{z+sqrt(z*z+1)} */
{
  complex rz;
  complex temp;
  complex e;

  e = ESetofComplex();

  temp = Mulz( z, z);
  temp = Addz( temp, e);
  temp = Sqrtz( temp);
  temp = Addz( z, temp);

  TRY {
    rz = Logz (temp);
    return rz;
  }
  CATCH (EX_LOGZ) {
    rz = NaNSetofComplex ();
    THROW (EX_ARCSINHZ);
    return rz;
  }
  ENDTRY;
}  


complex Arccoshz(complex z)   /* rz=Arccosh(z)=Log(z+sqrt(z*z-1)} */
{
  complex rz;
  complex tempz;
  complex e;

  e = ESetofComplex();

  tempz = Mulz( z, z);
  tempz = Subz( tempz, e);
  tempz = Sqrtz( tempz);
  tempz = Addz( z, tempz);

  TRY {
    rz = Logz (tempz);
    return rz;
  }
  CATCH (EX_LOGZ) {
    rz = NaNSetofComplex();
    THROW (EX_ARCCOSHZ);
    return rz;
  }
  ENDTRY;
}   


 /*complex Arctanhz(complex z)    // rz=Arctanh(z)=1/2*Log{(1+z)/(1-z)}
{
  complex rz;
  complex temp0,temp1,temp2;
  complex e,halfe;
 
  if( z.x == 0.0){
    rz = MakeComplex( 0.0, atan( z.y));
    return rz;
  }
  else{

    if( fabs( z.x) <= 1.0 && z.y == 0.0){
      rz = MakeComplex( atanh( z.x), 0.0);
      return rz;
    } 
    else{

      e = ESetofComplex();
      halfe = MakeComplex( 0.5, 0.0);

      temp0 = Addz( e, z);
      temp1 = Subz( e, z);

      TRY {
	temp2 = Divz (temp0, temp1);
	temp2 = Logz (temp2);
	rz = Mulz (halfe, temp2);
	return rz;
      }
      CATCH (EX_DIVZ) {
	rz = NaNSetofComplex();
	THROW (EX_ARCTANHZ);
	return rz;
      }
      CATCH (EX_LOGZ) {
	rz = NaNSetofComplex();
	THROW (EX_ARCTANHZ);
	return rz;
      }
      ENDTRY;
    }
  }
}    /*end. Arctanhz */



/* complexmath.c EOF */
