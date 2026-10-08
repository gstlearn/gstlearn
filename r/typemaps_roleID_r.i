/* =========================================================================
 * Typemaps for gstlrn::RoleID (R Interface)
 * Fixes OBJSXP / ENVSXP unwrapping inside R lists for overload dispatch.
 * Prevents memory buffer overread when converting ERole vs RoleID pointers.
 * ========================================================================= */

%{
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#if defined(__GNUG__) || defined(__clang__)
  #define SWIG_MAYBE_UNUSED __attribute__((unused))
#else
  #define SWIG_MAYBE_UNUSED
#endif

// Declare C symbol explicitly to prevent C++ name mangling issues with R's C-API
extern "C" {
    SEXP Rf_findVar(SEXP, SEXP);
}

// -------------------------------------------------------------------------
// Helper function to extract underlying SEXP pointer from SWIG R wrapper
// Marked inline and SWIG_MAYBE_UNUSED to prevent unused-function warnings.
// -------------------------------------------------------------------------
static inline SEXP unwrap_swig_r_obj(SEXP obj) SWIG_MAYBE_UNUSED;

static inline SEXP unwrap_swig_r_obj(SEXP obj) {
    if (obj == R_NilValue) return obj;

    int type = TYPEOF(obj);

    if (type == ENVSXP) {
        // SWIG S3 environment wrapper: search for .S3Type in the environment
        SEXP sym = Rf_install(".S3Type");
        SEXP ref = Rf_findVar(sym, obj);
        if (ref != R_UnboundValue) return ref;
    }
    else if (type == OBJSXP) {
        // SWIG S4 object wrapper: extract "ref" slot
        SEXP sym = Rf_install("ref");
        if (R_has_slot(obj, sym)) {
            return R_do_slot(obj, sym);
        }
    }
    return obj;
}
%}

// -------------------------------------------------------------------------
// Typecheck Typemap
// -------------------------------------------------------------------------
%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER) const gstlrn::RoleID&
{
    SEXP obj = $input;
    int type = TYPEOF(obj);

    if (type == VECSXP) {
        // Valid if it's an R list of length 2: list(RoleID|ERole, index)
        $1 = (Rf_length(obj) == 2) ? 1 : 0;
    } else if (type == ENVSXP || type == OBJSXP || type == EXTPTRSXP) {
        $1 = 1;
    } else {
        $1 = 0;
    }
}

// -------------------------------------------------------------------------
// In Typemap
// -------------------------------------------------------------------------
%typemap(in) const gstlrn::RoleID& (gstlrn::RoleID temp)
{
    SEXP obj = $input;
    gstlrn::Id index = 0;
    bool has_custom_index = false;

    // 1. Unpack list input: list(RoleID|ERole, index)
    if (TYPEOF(obj) == VECSXP)
    {
        if (Rf_length(obj) != 2)
        {
            Rf_error("RoleID list must contain exactly two elements: list(RoleID|ERole, index)");
        }

        SEXP first  = VECTOR_ELT(obj, 0);
        SEXP second = VECTOR_ELT(obj, 1);

        if (!Rf_isInteger(second) && !Rf_isReal(second))
        {
            Rf_error("Second list element must be a numeric integer index");
        }

        index = static_cast<gstlrn::Id>(Rf_asInteger(second));
        has_custom_index = true;

        // Extract underlying native pointer (OBJSXP/ENVSXP -> EXTPTRSXP)
        obj = unwrap_swig_r_obj(first);
    }
    else
    {
        obj = unwrap_swig_r_obj(obj);
    }

    // 2. Query SWIG type descriptors
    static swig_type_info *type_RoleID = SWIG_TypeQuery("gstlrn::RoleID *");
    static swig_type_info *type_ERole  = SWIG_TypeQuery("gstlrn::ERole *");

    gstlrn::RoleID *roleID = nullptr;
    gstlrn::ERole  *role   = nullptr;

    // 3. MUST check ERole* FIRST!
    // Attempting SWIG_ConvertPtr to RoleID* on an ERole* object causes a
    // memory buffer overread (reading 48 bytes out of a 40-byte ERole object).
    int res_role = SWIG_ConvertPtr(
        obj,
        reinterpret_cast<void**>(&role),
        type_ERole ? type_ERole : SWIGTYPE_p_gstlrn__ERole,
        0);

    if (SWIG_IsOK(res_role) && role != nullptr)
    {
        temp = gstlrn::RoleID(*role, index);
        $1 = &temp;
    }
    else
    {
        // 4. If not an ERole*, check for gstlrn::RoleID*
        int res_roleID = SWIG_ConvertPtr(
            obj,
            reinterpret_cast<void**>(&roleID),
            type_RoleID ? type_RoleID : SWIGTYPE_p_gstlrn__RoleID,
            0);

        if (SWIG_IsOK(res_roleID) && roleID != nullptr)
        {
            temp = *roleID;
            if (has_custom_index)
            {
                temp.setIndex(index);
            }
            $1 = &temp;
        }
        else
        {
            Rf_error("Expected RoleID, ERole, list(RoleID, index) or list(ERole, index)");
        }
    }
}

// -------------------------------------------------------------------------
// Apply typemaps to value signatures (after definitions to avoid Warning 453)
// -------------------------------------------------------------------------
%apply const gstlrn::RoleID& { gstlrn::RoleID };
