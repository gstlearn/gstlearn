// TODO: test directors feature (inheritance in target language)
%module(directors="1") gstlearn // TODO : configure this using CMake configure_file

// https://stackoverflow.com/a/26035360/3952924
#%import "doc/documentation.i"

// TODO: to be kept ?
%rename(__getitem__) Db::operator[];

// Include C++ library SWIG interface (Keep Order !!!!)
%include ../swig/swig_inc.i
%include ../swig/swig_exp.i

// For suppressing SWIG warning due to -keyword option
#pragma SWIG nowarn=511
#pragma SWIG nowarn=506

/// TODO, include numpy.i ?

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

