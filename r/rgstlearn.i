%module(directors="1") gstlearn // TODO : configure this using CMake configure_file

// Note : Keep order in this file!

// https://stackoverflow.com/a/26035360/3952924
#%import "doc/documentation.i"

//////////////////////////////////////////////////////////////
//       Specific typemaps and fragments for R language     //
//////////////////////////////////////////////////////////////

%begin %{
// For isnan
#include <cmath>
// For isinf
#include <math.h>
// For cout
#include <iostream>
%}

%fragment("ToCpp", "header")
{
  template <typename Type> int convertToCpp(SEXP obj, Type& value);
  
  template <> int convertToCpp(SEXP obj, int& value)
  {
    // Test argument
    if (obj == NULL) return SWIG_TypeError;

    int myres = SWIG_TypeError;
    if (Rf_length(obj) > 0) // Prevent NULL value from becoming NA
    {
      myres = SWIG_AsVal_int(obj, &value);
      //std::cout << "convertToCpp(int): value=" << value << std::endl;
      if (SWIG_IsOK(myres) && value == R_NaInt) // NA, NaN, Inf or out of bounds value becomes NA
        value = getNA<int>();
    }
    return myres;
  }
  template <> int convertToCpp(SEXP obj, double& value)
  {
    // Test argument
    if (obj == NULL) return SWIG_TypeError;
      
    int myres = SWIG_TypeError;
    if (Rf_length(obj) > 0) // Prevent NULL value from becoming NA
    {
       myres = SWIG_AsVal_double(obj, &value);
      //std::cout << "convertToCpp(double): value=" << value << std::endl;
      if (SWIG_IsOK(myres) && !R_finite(value)) // NA, NaN, Inf or out of bounds value becomes NA
        value = getNA<double>();
    }
    return myres; 
  }
  template <> int convertToCpp(SEXP obj, String& value)
  {
    // Test argument
    if (obj == NULL) return SWIG_TypeError;
      
    int myres = SWIG_TypeError;
    if (Rf_length(obj) > 0) // Prevent NULL value from being accepted
    {
      myres = SWIG_AsVal_std_string(obj, &value);
      //std::cout << "convertToCpp(String): value=" << value << std::endl;
      // No undefined
    }
    return myres;
  }
  template <> int convertToCpp(SEXP obj, float& value)
  {
    // Test argument
    if (obj == NULL) return SWIG_TypeError;
    
    int myres = SWIG_TypeError;
    if (Rf_length(obj) > 0) // Prevent NULL value from becoming NA
    {
       myres = SWIG_AsVal_float(obj, &value);
      //std::cout << "convertToCpp(float): value=" << value << std::endl;
      if (SWIG_IsOK(myres) && !R_finite(value)) // NA, NaN, Inf or out of bounds value becomes NA
        value = getNA<float>();
    }
    return myres; 
  }
  template <> int convertToCpp(SEXP obj, UChar& value)
  {
    // Test argument
    if (obj == NULL) return SWIG_TypeError;
    
    int myres = SWIG_TypeError;
    if (Rf_length(obj) > 0) // Prevent NULL value from becoming NA
    {
      int v = 0;
      myres = SWIG_AsVal_int(obj, &v);
      //std::cout << "convertToCpp(int): value=" << v << std::endl;
      if (myres == SWIG_OverflowError || 
          v < std::numeric_limits<UChar>::min() ||
          v > std::numeric_limits<UChar>::max()) // Out of bound value is error (no NA for UChar)
      {
        myres = SWIG_OverflowError;
      }
      else if (!SWIG_IsOK(myres) || v == R_NaInt) // NA, NaN or Inf is error (no NA for UChar)
      {
        myres = SWIG_TypeError;
      }
      else
      {
        value = static_cast<UChar>(v);
      }
    }
    return myres;
  }
  template <> int convertToCpp(SEXP obj, bool& value)
  {
    // Test argument
    if (obj == NULL) return SWIG_TypeError;
    
    int v = 0;
    int myres = SWIG_AsVal_int(obj, &v);
    //std::cout << "convertToCpp(bool): value=" << v << std::endl;
    if (v == 0)
      value = false;
    else
      value = true;
    return myres;
  }
  
  // Certainly not the most efficient way to convert vectors,
  // but at least, I can test each value for particular NAs
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
    vec.clear();
    
    // Test argument
    if (obj == NULL) return SWIG_TypeError;
    if (TYPEOF(obj) == EXTPTRSXP) return SWIG_TypeError;

    // Conversion
    int myres = SWIG_OK;
    int size = (int)Rf_length(obj);
    if (size > 0)
    {
      // Real vector
      vec.reserve(size);
      for (int i = 0; i < size && SWIG_IsOK(myres); i++)
      {
        SEXP item = getElem(obj,i); // item could be NULL or NIL
        if (TYPEOF(item) == NILSXP) continue; // If NIL, no error
        ValueType value;
        myres = convertToCpp(item, value);
        if (SWIG_IsOK(myres))
          vec.push_back(value);
      }
    }
    // else length = 0, empty vector
    return myres;
  }

  template <typename VectorVector>
  int vectorVectorToCpp(SEXP obj, VectorVector& vvec)
  {
    // Type definitions
    using InputVector = typename VectorVector::value_type;
    vvec.clear();
    
    // Test argument
    if (obj == NULL) return SWIG_TypeError;
    if (TYPEOF(obj) == EXTPTRSXP) return SWIG_TypeError;

    // Conversion
    int myres = SWIG_OK;
    int size = (int)Rf_length(obj);
    if (size == 1)
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
    // else length = 0, empty vector
    return myres;
  }
}

// Add typecheck typemaps for dispatching functions
%typemap(rtypecheck, noblock=1) const int&, int                               { length($arg) == 1 && (is.integer(unlist($arg)) || is.numeric(unlist($arg))) }
%typemap(rtypecheck, noblock=1) const double&, double                         { length($arg) == 1 &&  is.numeric(unlist($arg)) }
%typemap(rtypecheck, noblock=1) const String&, String                         { length($arg) == 1 &&  is.character(unlist($arg)) }
%typemap(rtypecheck, noblock=1) const float&, float                           { length($arg) == 1 &&  is.character(unlist($arg)) }
%typemap(rtypecheck, noblock=1) const UChar&, UChar                           { length($arg) == 1 && (is.integer(unlist($arg)) || is.numeric(unlist($arg))) }
%typemap(rtypecheck, noblock=1) const bool&, bool                             { length($arg) == 1 &&  is.logical(unlist($arg)) }
%typemap(rtypecheck, noblock=1) const VectorInt&, VectorInt                   { length($arg) == 0 || (length($arg) > 0 && (is.integer(unlist($arg)) || is.numeric(unlist($arg)))) }
%typemap(rtypecheck, noblock=1) const VectorDouble&, VectorDouble             { length($arg) == 0 || (length($arg) > 0 &&  is.numeric(unlist($arg))) }
%typemap(rtypecheck, noblock=1) const VectorString&, VectorString             { length($arg) == 0 || (length($arg) > 0 &&  is.character(unlist($arg))) }
%typemap(rtypecheck, noblock=1) const VectorFloat&, VectorFloat               { length($arg) == 0 || (length($arg) > 0 &&  is.numeric(unlist($arg))) }
%typemap(rtypecheck, noblock=1) const VectorUChar&, VectorUChar               { length($arg) == 0 || (length($arg) > 0 && (is.integer(unlist($arg)) || is.numeric(unlist($arg)))) }
%typemap(rtypecheck, noblock=1) const VectorBool&, VectorBool                 { length($arg) == 0 || (length($arg) > 0 &&  is.logical(unlist($arg))) }
%typemap(rtypecheck, noblock=1) const VectorVectorInt&, VectorVectorInt       { length($arg) == 0 || (length($arg) > 0 && 
                                                                               (length($arg[[1]]) == 0 || (length($arg[[1]]) > 0 && (is.integer(unlist($arg[[1]])) || is.numeric(unlist($arg[[1]])))))) }
%typemap(rtypecheck, noblock=1) const VectorVectorDouble&, VectorVectorDouble { length($arg) == 0 || (length($arg) > 0 && 
                                                                               (length($arg[[1]]) == 0 || (length($arg[[1]]) > 0 && is.numeric(unlist($arg[[1]]))))) }

%fragment("FromCpp", "header")
{  
  template <typename InputType> struct OutTraits;
  template <> struct OutTraits<int>     { using OutputType = int; };
  template <> struct OutTraits<double>  { using OutputType = double; };
  template <> struct OutTraits<String>  { using OutputType = String; };
  template <> struct OutTraits<float>   { using OutputType = float; };
  template <> struct OutTraits<UChar>   { using OutputType = UChar; };
  template <> struct OutTraits<bool>    { using OutputType = bool; };
  
  template <typename Type> typename OutTraits<Type>::OutputType convertFromCpp(const Type& value);
  template <> int convertFromCpp(const int& value)
  {
    //std::cout << "convertFromCpp(int): value=" << value << std::endl;
    if (isNA<int>(value))
      return R_NaInt;
    return value;
  }
  template <> double convertFromCpp(const double& value)
  {
    //std::cout << "convertFromCpp(double): value=" << value << std::endl;
    if (isNA<double>(value))
      return R_NaReal;
    return value;
  }
  template <> String convertFromCpp(const String& value)
  {
    //std::cout << "convertFromCpp(String): value=" << value << std::endl;
    return value; // No special conversion provided
  }
  template <> float convertFromCpp(const float& value)
  {
    //std::cout << "convertFromCpp(float): value=" << value << std::endl;
    if (isNA<float>(value))
      return R_NaReal;
    return value;
  }
  template <> UChar convertFromCpp(const UChar& value)
  {
    //std::cout << "convertFromCpp(UChar): value=" << value << std::endl;
    return value; // No special conversion provided
  }
    template <> bool convertFromCpp(const bool& value)
  {
    //std::cout << "convertFromCpp(bool): value=" << value << std::endl;
    return value; // No special conversion provided
  }
  
  template <typename Type> SEXP objectFromCpp(const Type& value);
  template <> SEXP objectFromCpp(const int& value)
  {
    return Rf_ScalarInteger(convertFromCpp(value));
  }
  template <> SEXP objectFromCpp(const double& value)
  {
    return Rf_ScalarReal(convertFromCpp(value));
  }
  template <> SEXP objectFromCpp(const String& value)
  {
    return Rf_ScalarString(Rf_mkChar(convertFromCpp(value).c_str()));
  }
  template <> SEXP objectFromCpp(const float& value)
  {
    return Rf_ScalarReal(static_cast<double>(convertFromCpp(value)));
  }
  template <> SEXP objectFromCpp(const UChar& value)
  {
    return Rf_ScalarInteger(static_cast<int>(convertFromCpp(value)));
  }
  template <> SEXP objectFromCpp(const bool& value)
  {
    return Rf_ScalarLogical(static_cast<int>(convertFromCpp(value)));
  }
  
  template <typename Vector>
  int vectorFromCpp(SEXP* obj, const Vector& vec)
  {
    // Type definitions
    int myres = SWIG_TypeError;
    using SizeType = typename Vector::size_type;
 
    // Test NA values
    auto vec2 = vec.getVector();
    SizeType size = vec2.size();
    for(SizeType i = 0; i < size; i++)
    {
      vec2[i] = convertFromCpp(vec2[i]);
    }
    // Convert to R vector
    *obj = swig::from(vec2);
    myres = (*obj) == NULL ? SWIG_TypeError : SWIG_OK;
    return myres;
  }

  template <typename VectorVector>
  int vectorVectorFromCpp(SEXP* obj, const VectorVector& vec)
  {
    // Type definitions
    int myres = SWIG_TypeError;
    using SizeType = typename VectorVector::size_type;

    // https://cpp.hotexamples.com/examples/-/-/Rf_allocVector/cpp-rf_allocvector-function-examples.html
    SizeType size = vec.size();
    PROTECT(*obj = Rf_allocVector(VECSXP, size));
    if(*obj != NULL)
    {
      myres = SWIG_OK;
      for(SizeType i = 0; i < size && SWIG_IsOK(myres); i++)
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

// This for automatically convert R lists to externalptr
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

# Add automatic display for all AStringable objects (see onAttach comment below)
setMethod(f = "show", signature = "_p_AStringable", definition = function(object){ AStringable_display(object) })

##
## Add operator [] to VectorXXX R class [1-based index] ##
## ---------------------------------------------------- ##

"getVitem" <-
function(x, i)
{
  idx = as.integer(i)
  if (length(idx) > 1) {
    sapply(idx, function(n) {
      if (n < 1 || n > x$length())
        stop("Index out of range")
      x$getAt(n-1)
    }) 
  }
  else {
    x$getAt(idx-1)
  }
}

"setVitem" <-
function(x, i, value)
{
  idx = as.integer(i)
  if (length(idx) > 1) {
    sapply(1:length(i), function(n) {
      if (i[n] < 1 || i[n] > x$length())
        stop("Index out of range")
      x$setAt(i[n]-1, value[n])
    })
  }
  else {
    x$setAt(idx-1, value)
  }
  x
}

setMethod('[',    '_p_VectorTT_int_t',                  getVitem)
setMethod('[<-',  '_p_VectorTT_int_t',                  setVitem)
setMethod('[',    '_p_VectorTT_double_t',               getVitem)
setMethod('[<-',  '_p_VectorTT_double_t',               setVitem)
setMethod('[',    '_p_VectorTT_String_t',               getVitem) # TODO : Different from swigex and don't know why (_p_VectorTT_std__string_t)
setMethod('[<-',  '_p_VectorTT_String_t',               setVitem) # TODO : Different from swigex and don't know why (_p_VectorTT_std__string_t)
setMethod('[',    '_p_VectorTT_float_t',                getVitem)
setMethod('[<-',  '_p_VectorTT_float_t',                setVitem)
setMethod('[',    '_p_VectorTT_UChar_t',                getVitem)
setMethod('[<-',  '_p_VectorTT_UChar_t',                setVitem)
setMethod('[',    '_p_VectorNumTT_int_t',               getVitem)
setMethod('[<-',  '_p_VectorNumTT_int_t',               setVitem)
setMethod('[',    '_p_VectorNumTT_double_t',            getVitem)
setMethod('[<-',  '_p_VectorNumTT_double_t',            setVitem)
setMethod('[',    '_p_VectorNumTT_float_t',             getVitem)
setMethod('[<-',  '_p_VectorNumTT_float_t',             setVitem)
setMethod('[',    '_p_VectorNumTT_UChar_t',             getVitem)
setMethod('[<-',  '_p_VectorNumTT_UChar_t',             setVitem)
setMethod('[[',   '_p_VectorTT_VectorNumTT_int_t_t',    getVitem)
setMethod('[[<-', '_p_VectorTT_VectorNumTT_int_t_t',    setVitem)
setMethod('[[',   '_p_VectorTT_VectorNumTT_double_t_t', getVitem)
setMethod('[[<-', '_p_VectorTT_VectorNumTT_double_t_t', setVitem)
setMethod('[[',   '_p_VectorTT_VectorNumTT_float_t_t',  getVitem)
setMethod('[[<-', '_p_VectorTT_VectorNumTT_float_t_t',  setVitem)

##
## Add operator [] to Db R class ##
## ----------------------------- ##

"is.undef" <- function(x)
{
  if (length(x) <= 0) return(TRUE)
  if (length(x) > 1) return(FALSE)
  if (is.na(x)) return(TRUE)
  return(FALSE) 
}

"getDbitem" <-
function (x,i,j,...,drop=TRUE)
{
  db   <- x
  irow <- NA
  icol <- NA
  flag_i_defined = (deparse(substitute(i)) != "")
  flag_j_defined = (deparse(substitute(j)) != "")
  nargs = 0
  if (flag_i_defined) nargs = nargs + 1
  if (flag_j_defined) nargs = nargs + 1
  
  nech = db$getSampleNumber()
  ncol = db$getColumnNumber()

  if (nargs == 2) {
  
    # Case of both arguments are defined
    
    if (is.logical(i) & length(i)==nech) i=(1:nech)[i]
    if (is.numeric(i)) irow <- i
  
    if (is.numeric(j)) {
      icol <- j
    } else {
      icol <- grep(paste(i,collapse="|"),db$getAllNames())
      if (length(icol) > 0) icol = icol - 1
    }
  } else if (nargs == 1) {
  
    # Case where only one argument in defined

    if (flag_i_defined) k = i
    if (flag_j_defined) k = j
    if (is.numeric(k)) {
      icol <- k
    } else {
      icol <- grep(paste(k,collapse="|"),db$getAllNames())
      if (length(icol) > 0) icol = icol - 1
    }
  }

  icol_dim = icol
  irow_dim = irow
  if (is.undef(irow)) irow_dim = seq(0, nech - 1)
  if (is.undef(icol)) icol_dim = seq(0, ncol - 1)
  col_names = db$getNamesByColIdx(icol_dim)
  row_names = irow_dim
  ncol_dim = length(icol_dim)
  nrow_dim = length(irow_dim)
    
  if (length(icol) <= 0)
  {
    messerr("The variable does not exist")
    messerr("This is not authorized in this function")
    stop()
  } else {
  if (! is.undef(irow) && ! is.undef(icol))
     res <- db$getValuesByColIdx(irow_dim,icol_dim)
   if (! is.undef(irow) &&   is.undef(icol))
     res <- db$getArrayBySample(irow_dim)
   if (  is.undef(irow) && ! is.undef(icol))
     res <- db$getColumnsByColIdx(icol_dim)
  if (  is.undef(irow) &&   is.undef(icol))
     res <- db$getAllColumns()
  if (nrow_dim > 1 && ncol_dim > 1)
  {
    res <- as.data.frame(matrix(res, nrow=nrow_dim, ncol=ncol_dim))
    names(res) = col_names
    row.names(res) = row_names
  }
  }
  res
}

"setDbitem" <-
  function (x,i,j,...,value)
{
  db   <- x
  irow <- NA
  icol <- NA
  flag_i_defined = (deparse(substitute(i)) != "")
  flag_j_defined = (deparse(substitute(j)) != "")
  nargs = 0
  if (flag_i_defined) nargs = nargs + 1
  if (flag_j_defined) nargs = nargs + 1

  nech = db$getSampleNumber()
  ncol = db$getColumnNumber()
  value = as.numeric(unlist(value))

  new_names = "New"
  if (nargs == 2) {
  
    # Both arguments are defined
    
    if (is.numeric(i) || is.logical(i)) irow <- i
  
    if (is.numeric(j) || is.logical(j)) {
      icol <- j
    } else {
      icol <- grep(paste(i,collapse="|"),db$getAllNames())
      if (length(icol) > 0) icol = icol - 1
      new_names = j
    }
  } else if (nargs == 1) {
  
    # Only one argument is defined: it corresponds to the column

    if (flag_i_defined) k = i
    if (flag_j_defined) k = j
    if (is.numeric(k)) {
      icol <- k
    } else {
      icol <- grep(paste(k,collapse="|"),db$getAllNames())
      if (length(icol) > 0) icol = icol - 1
      new_names = k
    }
  }

  icol_dim = icol
  irow_dim = irow
  if (is.undef(irow)) irow_dim = seq(0, nech - 1)
  if (is.undef(icol)) icol_dim = seq(0, ncol - 1)
  col_names = db$getNamesByColIdx(icol_dim)
  row_names = irow_dim
  ncol_dim = length(icol_dim)
  nrow_dim = length(irow_dim)

  if (length(icol) <= 0)
  {
  
   # Case of a new variable
  
    icol = db$addColumns(value, new_names)
  }
  else
  {

    # Case of already an existing variable: replacement
    
  if (! is.undef(irow) && ! is.undef(icol))
    db$setValuesByColIdx(irow_dim,icol_dim,value)
   if (! is.undef(irow) &&   is.undef(icol))
    db$setArrayBySample(irow_dim,value)
   if (  is.undef(irow) && ! is.undef(icol))
    db$setColumnsByColIdx(value,icol_dim)
   if (  is.undef(irow) &&   is.undef(icol))
    db$setAllColumns(value)
  }
  db
}

setMethod('[',    '_p_Db',               getDbitem)
setMethod('[<-',  '_p_Db',               setDbitem)
setMethod('[',    '_p_DbGrid',           getDbitem)
setMethod('[<-',  '_p_DbGrid',           setDbitem)
setMethod('[',    '_p_VectorTT_int_t',                  getVitem)
setMethod('[<-',  '_p_VectorTT_int_t',                  setVitem)
setMethod('[',    '_p_VectorTT_double_t',               getVitem)
setMethod('[<-',  '_p_VectorTT_double_t',               setVitem)
setMethod('[',    '_p_VectorTT_String_t',               getVitem) # TODO : Different from myfibo and don't know why (_p_VectorTT_std__string_t)
setMethod('[<-',  '_p_VectorTT_String_t',               setVitem) # TODO : Different from myfibo and don't know why (_p_VectorTT_std__string_t)
setMethod('[',    '_p_VectorTT_float_t',                getVitem)
setMethod('[<-',  '_p_VectorTT_float_t',                setVitem)
setMethod('[',    '_p_VectorTT_UChar_t',                getVitem)
setMethod('[<-',  '_p_VectorTT_UChar_t',                setVitem)
setMethod('[',    '_p_VectorNumTT_int_t',               getVitem)
setMethod('[<-',  '_p_VectorNumTT_int_t',               setVitem)
setMethod('[',    '_p_VectorNumTT_double_t',            getVitem)
setMethod('[<-',  '_p_VectorNumTT_double_t',            setVitem)
setMethod('[',    '_p_VectorNumTT_float_t',             getVitem)
setMethod('[<-',  '_p_VectorNumTT_float_t',             setVitem)
setMethod('[',    '_p_VectorNumTT_UChar_t',             getVitem)
setMethod('[<-',  '_p_VectorNumTT_UChar_t',             setVitem)
setMethod('[[',   '_p_VectorTT_VectorNumTT_int_t_t',    getVitem)
setMethod('[[<-', '_p_VectorTT_VectorNumTT_int_t_t',    setVitem)
setMethod('[[',   '_p_VectorTT_VectorNumTT_double_t_t', getVitem)
setMethod('[[<-', '_p_VectorTT_VectorNumTT_double_t_t', setVitem)
setMethod('[[',   '_p_VectorTT_VectorNumTT_float_t_t',  getVitem)
setMethod('[[<-', '_p_VectorTT_VectorNumTT_float_t_t',  setVitem)

##
## Add operator [] to Db R class ##
## ----------------------------- ##

"is.undef" <- function(x)
{
  if (length(x) <= 0) return(TRUE)
  if (length(x) > 1) return(FALSE)
  if (is.na(x)) return(TRUE)
  return(FALSE) 
}

"getDbitem" <-
function (x,i,j,...,drop=TRUE)
{
  db   <- x
  irow <- NA
  icol <- NA
  flag_i_defined = (deparse(substitute(i)) != "")
  flag_j_defined = (deparse(substitute(j)) != "")
  nargs = 0
  if (flag_i_defined) nargs = nargs + 1
  if (flag_j_defined) nargs = nargs + 1
  
  nech = db$getSampleNumber()
  ncol = db$getColumnNumber()

  if (nargs == 2) {
  
    # Case of both arguments are defined
    
    if (is.logical(i) & length(i)==nech) i=(1:nech)[i]
    if (is.numeric(i)) irow <- i
  
    if (is.numeric(j)) {
      icol <- j
    } else {
      icol <- grep(paste(i,collapse="|"),db$getAllNames())
      if (length(icol) > 0) icol = icol - 1
    }
  } else if (nargs == 1) {
  
    # Case where only one argument in defined

    if (flag_i_defined) k = i
    if (flag_j_defined) k = j
    if (is.numeric(k)) {
      icol <- k
    } else {
      icol <- grep(paste(k,collapse="|"),db$getAllNames())
      if (length(icol) > 0) icol = icol - 1
    }
  }

  icol_dim = icol
  irow_dim = irow
  if (is.undef(irow)) irow_dim = seq(0, nech - 1)
  if (is.undef(icol)) icol_dim = seq(0, ncol - 1)
  col_names = db$getNamesByColIdx(icol_dim)
  row_names = irow_dim
  ncol_dim = length(icol_dim)
  nrow_dim = length(irow_dim)
    
  if (length(icol) <= 0) {
    messerr("The variable does not exist")
    messerr("This is not authorized in this function")
    stop()
  } else {
    if (! is.undef(irow) && ! is.undef(icol))
      res <- db$getValuesByColIdx(irow_dim,icol_dim)
    if (! is.undef(irow) &&   is.undef(icol))
      res <- db$getArrayBySample(irow_dim)
    if (  is.undef(irow) && ! is.undef(icol))
      res <- db$getColumnsByColIdx(icol_dim)
    if (  is.undef(irow) &&   is.undef(icol))
      res <- db$getAllColumns()
    if (nrow_dim > 1 && ncol_dim > 1)
    {
      res <- as.data.frame(matrix(res, nrow=nrow_dim, ncol=ncol_dim))
      names(res) = col_names
      row.names(res) = row_names
    }
  }
  res
}

"setDbitem" <-
  function (x,i,j,...,value)
{
  db   <- x
  irow <- NA
  icol <- NA
  flag_i_defined = (deparse(substitute(i)) != "")
  flag_j_defined = (deparse(substitute(j)) != "")
  nargs = 0
  if (flag_i_defined) nargs = nargs + 1
  if (flag_j_defined) nargs = nargs + 1

  nech = db$getSampleNumber()
  ncol = db$getColumnNumber()
  value = as.numeric(unlist(value))

  new_names = "New"
  if (nargs == 2) {
  
    # Both arguments are defined
    
    if (is.numeric(i) || is.logical(i)) irow <- i
  
    if (is.numeric(j) || is.logical(j)) {
      icol <- j
    } else {
      icol <- grep(paste(i,collapse="|"),db$getAllNames())
      if (length(icol) > 0) icol = icol - 1
      new_names = j
    }
  } else if (nargs == 1) {
  
    # Only one argument is defined: it corresponds to the column

    if (flag_i_defined) k = i
    if (flag_j_defined) k = j
    if (is.numeric(k)) {
      icol <- k
    } else {
      icol <- grep(paste(k,collapse="|"),db$getAllNames())
      if (length(icol) > 0) icol = icol - 1
      new_names = k
    }
  }

  icol_dim = icol
  irow_dim = irow
  if (is.undef(irow)) irow_dim = seq(0, nech - 1)
  if (is.undef(icol)) icol_dim = seq(0, ncol - 1)
  col_names = db$getNamesByColIdx(icol_dim)
  row_names = irow_dim
  ncol_dim = length(icol_dim)
  nrow_dim = length(irow_dim)

  if (length(icol) <= 0)
  {
    # Case of a new variable
    icol = db$addColumns(value, new_names)
  }
  else
  {
    # Case of already an existing variable: replacement
    
    if (! is.undef(irow) && ! is.undef(icol))
      db$setValuesByColIdx(irow_dim,icol_dim,value)
    if (! is.undef(irow) &&   is.undef(icol))
      db$setArrayBySample(irow_dim,value)
    if (  is.undef(irow) && ! is.undef(icol))
      db$setColumnsByColIdx(value,icol_dim)
    if (  is.undef(irow) &&   is.undef(icol))
      db$setAllColumns(value)
  }
  db
}

setMethod('[',    '_p_Db',               getDbitem)
setMethod('[<-',  '_p_Db',               setDbitem)
setMethod('[',    '_p_DbGrid',           getDbitem)
setMethod('[<-',  '_p_DbGrid',           setDbitem)

%}