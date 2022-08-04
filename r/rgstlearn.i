%module(directors="1") gstlearn // TODO : configure this using CMake configure_file

// Note : Keep order in this file!

// https://stackoverflow.com/a/26035360/3952924
#%import "doc/documentation.i"

//////////////////////////////////////////////////////////////
//       Specific typemaps and fragments for R language     //
//////////////////////////////////////////////////////////////

%fragment("ToCpp", "header")
{

  template <typename Type> int convertToCpp(SEXP obj, Type& value);
  
  template <> int convertToCpp(SEXP obj, int& value)
  {
    // TODO : Handle undefined or NA values
    return SWIG_AsVal_int(obj, &value);
  }
  template <> int convertToCpp(SEXP obj, double& value)
  {
    // TODO : Handle undefined or NA values
    return SWIG_AsVal_double(obj, &value);
  }
  template <> int convertToCpp(SEXP obj, String& value)
  {
    // TODO : Handle undefined or NA values
    return SWIG_AsVal_std_string(obj, &value);
  }
  template <> int convertToCpp(SEXP obj, float& value)
  {
    // TODO : Handle undefined or NA values
    return SWIG_AsVal_float(obj, &value);
  }
  template <> int convertToCpp(SEXP obj, UChar& value)
  {
    // TODO : Handle undefined or NA values
    int v = 0;
    int myres = SWIG_AsVal_int(obj, &v);
    if (v < 0 || v > 255)
      myres = SWIG_TypeError;
    else
      value = static_cast<UChar>(v);
    return myres;
  }
  template <> int convertToCpp(SEXP obj, bool& value)
  {
    int v = 0;
    int myres = SWIG_AsVal_int(obj, &v);
    if (v == 0)
      value = false;
    else
      value = true;
    return myres;
  }
  
  // Certainly not the most efficient way to convert vectors.
  // But at least, I can test each value for particular NAs
  SEXP getElem(SEXP obj, int i)
  {
    if (Rf_isInteger(obj))      return Rf_ScalarInteger(INTEGER(obj)[i]);
    if (Rf_isReal(obj))         return Rf_ScalarReal(REAL(obj)[i]);
    if (Rf_isString(obj))       return Rf_ScalarString(STRING_ELT(obj, i));
    if (Rf_isLogical(obj))      return Rf_ScalarLogical(LOGICAL(obj)[i]);
    if (TYPEOF(obj) == VECSXP)  return VECTOR_ELT(obj, i);
    return SEXP();
  }
  
  template <typename Vector>
  int vectorToCpp(SEXP obj, Vector& vec)
  {
    // Type definitions
    using ValueType = typename Vector::value_type;
    
    // Conversion
    vec.clear();
    int myres = SWIG_OK;
    int size = (int)Rf_length(obj);
    if (size < 0)
    {
      // Not a vector
      ValueType value;
      // Try to convert
      myres = convertToCpp(obj, value);
      if (SWIG_IsOK(myres))
        vec.push_back(value);
    }
    else if (size > 0)
    {
      // Real vector
      vec.reserve(size);
      for (int i = 0; i < size && SWIG_IsOK(myres); i++)
      {
        SEXP item = getElem(obj,i);
        ValueType value;
        myres = convertToCpp(item, value);
        if (SWIG_IsOK(myres))
          vec.push_back(value);
      }
    }
    return myres;
  }
  
  template <typename VectorVector>
  int vectorVectorToCpp(SEXP obj, VectorVector& vvec)
  {
    // Type definitions
    using InputVector = typename VectorVector::value_type;
    
    // Conversion
    int myres = SWIG_OK;
    int size = (int)Rf_length(obj);
    if (size <= 1)
    {
      // Not a vector (or a single value)
      InputVector vec;
      // Try to convert
      myres = vectorToCpp(obj, vec);
      if (SWIG_IsOK(myres))
        vvec.push_back(vec);
    }
    else if (size > 1)
    {
      for (int i = 0; i < size && SWIG_IsOK(myres); i++)
      {
        SEXP item = getElem(obj,i);
        InputVector vec;
        myres = vectorToCpp(item, vec);
        if (SWIG_IsOK(myres))
          vvec.push_back(vec);
      }
    }
    return myres;
  }
}

// Add typecheck typemaps for dispatching functions
// TODO : rtypecheck doesn't take into account precedence. Vector with one item will be seen as a single scalar
%typemap(rtypecheck, noblock=1) const int&, int                     { length($arg) == 1 && (is.integer($arg) || is.numeric($arg)) }
%typemap(rtypecheck, noblock=1) const double&, double               { length($arg) == 1 &&  is.numeric($arg) }
%typemap(rtypecheck, noblock=1) const String&, String               { length($arg) == 1 &&  is.character($arg) }
%typemap(rtypecheck, noblock=1) const float&, float                 { length($arg) == 1 &&  is.character($arg) }
%typemap(rtypecheck, noblock=1) const UChar&, UChar { length($arg) == 1 && (is.integer($arg) || is.numeric($arg)) }
%typemap(rtypecheck, noblock=1) const bool&, bool                   { length($arg) == 1 &&  is.logical($arg) }
%typemap(rtypecheck, noblock=1) const VectorInt&, VectorInt         { length($arg)  > 1 && (is.integer($arg) || is.numeric($arg)) }
%typemap(rtypecheck, noblock=1) const VectorDouble&, VectorDouble   { length($arg)  > 1 &&  is.numeric($arg) }
%typemap(rtypecheck, noblock=1) const VectorString&, VectorString   { length($arg)  > 1 &&  is.character($arg) }
%typemap(rtypecheck, noblock=1) const VectorFloat&, VectorFloat     { length($arg)  > 1 &&  is.numeric($arg) }
%typemap(rtypecheck, noblock=1) const VectorUChar&, VectorUChar     { length($arg)  > 1 && (is.integer($arg) || is.numeric($arg)) }
%typemap(rtypecheck, noblock=1) const VectorBool&, VectorBool       { length($arg)  > 1 && is.logical($arg) }

%fragment("FromCpp", "header")
{
  template <typename Vector>
  int vectorFromCpp(SEXP* obj, const Vector& vec)
  {
    *obj = swig::from(vec.getVector());
    return (*obj) == NULL ? -1 : 0;
  }

  template <typename VectorVector>
  int vectorVectorFromCpp(SEXP* obj, const VectorVector& vec)
  {
    int myres = SWIG_TypeError;
    // https://cpp.hotexamples.com/examples/-/-/Rf_allocVector/cpp-rf_allocvector-function-examples.html
    const unsigned int size = vec.size();
    PROTECT(*obj = Rf_allocVector(VECSXP, size));
    if(*obj != NULL)
    {
      myres = SWIG_OK;
      for(unsigned int i = 0; i < size && SWIG_IsOK(myres); i++)
      {
        SEXP rvec;
        myres = vectorFromCpp(&rvec, vec.at(i));
        if (SWIG_IsOK(myres))
          SET_VECTOR_ELT(*obj, i, rvec);
      }
    }
    UNPROTECT(1);
    return myres;
  }
}

%typemap(scoerceout) int,    int*,    int&,
                     double, double*, double&,
                     float,  float*,  float&,
                     UChar,  UChar*,  UChar&,
                     bool,   bool*,   bool&
 %{    %}

%typemap(scoerceout) VectorInt,    VectorInt*,    VectorInt&,
                     VectorDouble, VectorDouble*, VectorDouble&,
                     VectorString, VectorString*, VectorString&,
                     VectorFloat,  VectorFloat*,  VectorFloat&,
                     VectorUChar,  VectorUChar*,  VectorUChar&,
                     VectorBool,   VectorBool*,   VectorBool&
 %{    %}

%typemap(scoerceout) VectorVectorInt,    VectorVectorInt*,    VectorVectorInt&,
                     VectorVectorDouble, VectorVectorDouble*, VectorVectorDouble&,
                     VectorVectorFloat,  VectorVectorFloat*,  VectorVectorFloat&
 %{    %}

//////////////////////////////////////////////////////////////
//         C++ library SWIG includes and typemaps           //
//////////////////////////////////////////////////////////////

%include ../swig/swig_inc.i

//////////////////////////////////////////////////////////////
//                C++ library SWIG interface                //
//////////////////////////////////////////////////////////////

%include ../swig/swig_exp.i

//////////////////////////////////////////////////////////////
//                    Add C++ extension                     //
//////////////////////////////////////////////////////////////

// TODO: to be kept ?
%rename(__getitem__) Db::operator[];

%{
  #include <stdio.h>
  
  #include <R_ext/Print.h>
  #include <R_ext/Error.h>
  
  void R_Write(const char *string)
  {
    if (string == NULL) return;
    int length = strlen(string);
    if (length > 0) Rprintf(string);
  }
  
  void R_Warning(const char *string)
  {
    if (string == NULL) return;
    int length = strlen(string);
    if (length > 0) Rprintf(string);
  }
  #ifndef _WIN32
  #define R_INTERFACE_PTRS 1
  #include <Rinterface.h>
  void R_Read(const char *prompt,char *answer)
  {
  #define LNG 100
    char reponse[LNG];
    (void) strcpy(reponse,"");
    (void) strcpy(answer ,"");
    ptr_R_ReadConsole(prompt,(unsigned char*) reponse,LNG,0);
    int longueur = strlen(reponse);
    reponse[longueur-1] = '\0';
    if (strlen(reponse) > 0) (void) strcpy(answer,reponse);
  }
  #endif // Not _WIN32
  
  void R_Exit(void)
  {
    Rf_error("Abort caught by C++ and return to R");
  }
%}

%init %{
  redefine_message(R_Write);
  redefine_error(R_Warning);
#ifndef _WIN32
  redefine_read(R_Read);
#endif // Not _WIN32
  redefine_exit(R_Exit);
%}

//////////////////////////////////////////////////////////////
//       Add target language additional features below      //
//////////////////////////////////////////////////////////////

%insert(s)
%{

## Add operator [] to VectorXXX R class [1-based index] ##
## ---------------------------------------------------- ##

"getitem" <-
function(x, i)
{
  idx = as.integer(i)
  sapply(idx, function(n) {
    if (n < 1 || n > x$length())
      stop("Index out of range")
    x$get(n-1)
  })
}
"setitem" <-
function(x, i, value)
{
  idx = as.integer(i)
  sapply(1:length(i), function(n) {
    if (i[n] < 1 || i[n] > x$length())
      stop("Index out of range")
    x$set(i[n]-1, value[n])
  })
  x
}

setMethod('[',   '_p_VectorTT_int_t',                  getitem)
setMethod('[<-', '_p_VectorTT_int_t',                  setitem)
setMethod('[',   '_p_VectorTT_double_t',               getitem)
setMethod('[<-', '_p_VectorTT_double_t',               setitem)
setMethod('[',   '_p_VectorTT_String_t',               getitem) # TODO : Different from myfibo and don't know why (_p_VectorTT_std__string_t)
setMethod('[<-', '_p_VectorTT_String_t',               setitem) # TODO : Different from myfibo and don't know why (_p_VectorTT_std__string_t)
setMethod('[',   '_p_VectorNumTT_int_t',               getitem)
setMethod('[<-', '_p_VectorNumTT_int_t',               setitem)
setMethod('[',   '_p_VectorNumTT_double_t',            getitem)
setMethod('[<-', '_p_VectorNumTT_double_t',            setitem)
#setMethod('[',   '_p_VectorTT_VectorNumTT_int_t_t',    getitem) # TODO : VectorVectorXXX getitem doesn't work yet
#setMethod('[<-', '_p_VectorTT_VectorNumTT_int_t_t',    setitem) # TODO : VectorVectorXXX setitem doesn't work yet
#setMethod('[',   '_p_VectorTT_VectorNumTT_double_t_t', getitem) # TODO : VectorVectorXXX getitem doesn't work yet
#setMethod('[<-', '_p_VectorTT_VectorNumTT_double_t_t', setitem) # TODO : VectorVectorXXX setitem doesn't work yet

%}