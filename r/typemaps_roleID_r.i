/* =========================================================================
 * Typemaps for gstlrn::RoleID (R Interface)
 * Fixes OBJSXP / ENVSXP unwrapping inside R lists
 * ========================================================================= */

// -------------------------------------------------------------------------
// Helper function to extract underlying SEXP pointer from SWIG R wrapper
// Marked inline and [[maybe_unused]] to prevent unused-function warnings.
// -------------------------------------------------------------------------
%{
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#if defined(__GNUG__) || defined(__clang__)
  #define SWIG_MAYBE_UNUSED __attribute__((unused))
#else
  #define SWIG_MAYBE_UNUSED
#endif

static inline SEXP unwrap_swig_r_obj(SEXP obj) SWIG_MAYBE_UNUSED;

static inline SEXP unwrap_swig_r_obj(SEXP obj) {
    if (obj == R_NilValue) return obj;

    if (TYPEOF(obj) == ENVSXP) {
        // SWIG S3 environment wrapper: extract .S3Type or ref pointer
        SEXP sym = Rf_install(".S3Type");
        SEXP ref = Rf_findVarInFrame(obj, sym);
        if (ref != R_UnboundValue) return ref;
    }
    else if (TYPEOF(obj) == OBJSXP) {
        // SWIG S4 object wrapper: extract "ref" slot
        SEXP sym = Rf_install("ref");
        if (R_has_slot(obj, sym)) {
            return R_do_slot(obj, sym);
        }
    }
    return obj;
}
%}
