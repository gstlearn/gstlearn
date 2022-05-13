// TODO: restore directors feature
%module(directors="1") gstlearn // TODO : configure this using CMake configure_file

// https://stackoverflow.com/a/26035360/3952924
#%import "doc/documentation.i"

%rename(__getitem__) Db::operator[];

// Include C++ library SWIG interface (Keep Order !!!!)
%include ../swig/swig_inc.i
%include ../swig/swig_exp.i

// For suppressing SWIG warning due to -keyword option
#pragma SWIG nowarn=511
#pragma SWIG nowarn=506

/// TODO, include numpy.i ?

%{
#include <string>
#include <sstream>

#define R_INTERFACE_PTRS 1
#include <Rembedded.h>
#include <Rinterface.h>

// TODO : Clean the code below

// Look at rgeostats_loader.R for knowing how to display R object directly
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

void R_Exit(void)
{
  Rf_error("Abort caught by C and return to R");
}

%}

%init %{
  redefine_message(R_Write);
  redefine_error(R_Warning);
  redefine_read(R_Read);
  redefine_exit(R_Exit);
%}

