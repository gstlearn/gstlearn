/***************************************************************************/
/*                                                                         */
/*  R typemaps for gstlrn::ColID                                           */
/*                                                                         */
/*  Dual-binding support for both:                                         */
/*    1. const gstlrn::ColID&  (Preferred API facade in Db.hpp)            */
/*    2. gstlrn::ColID&&       (Legacy / Low-level engine in DbData.hpp)   */
/*                                                                         */
/*  This typemap transparently converts high-level R arguments into        */
/*  a C++ ColID temporary instance.                                        */
/*                                                                         */
/*  Accepted R inputs:                                                     */
/*                                                                         */
/*    - "name"                  -> ColID("name")                           */
/*    - "name", version         -> ColID("name", version)                  */
/*    - icol                    -> ColID(icol)                             */
/*    - icol, version           -> ColID(icol, version)                    */
/*    - RoleID                  -> ColID(RoleID)                           */
/*    - RoleID, version         -> ColID(RoleID, version)                  */
/*    - ERole                   -> ColID(ERole)                            */
/*    - ERole, version          -> ColID(ERole, version)                   */
/*                                                                         */
/*  Memory Safety:                                                         */
/*    A heap-allocated ColID object is constructed in %typemap(in) and     */
/*    automatically reclaimed in %typemap(freearg) after the C++ method    */
/*    call returns.                                                        */
/*                                                                         */
/***************************************************************************/

%{
#include <R.h>
#include <Rinternals.h>
%}

/* ========================================================================= */
/* 1. BYPASS SWIG R AUTOMATIC S4 COERCION & AS() CONVERSIONS                 */
/* ========================================================================= */

%typemap(rtype)      gstlrn::ColID, const gstlrn::ColID&, gstlrn::ColID&& "ANY";
%typemap(scoercein)  gstlrn::ColID, const gstlrn::ColID&, gstlrn::ColID&& "$input";
%typemap(sclass)     gstlrn::ColID, const gstlrn::ColID&, gstlrn::ColID&& "ANY";

/* ========================================================================= */
/* 2. UNIVERSAL TYPECHECK FOR SWIG OVERLOAD RESOLUTION                       */
/* ========================================================================= */

%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER)
    gstlrn::ColID,
    const gstlrn::ColID&,
    gstlrn::ColID&&
{
    SEXP obj = $input;
    if (TYPEOF(obj) == VECSXP && LENGTH(obj) == 2)
    {
        obj = VECTOR_ELT(obj, 0);
    }

    const SEXPTYPE t = TYPEOF(obj);
    if (t == STRSXP || t == INTSXP || t == REALSXP || t == EXTPTRSXP || Rf_isS4(obj))
    {
        $1 = 1;
    }
    else
    {
        $1 = 0;
    }
}

/* ========================================================================= */
/* 3. INPUT CONVERSION TYPEMAP FOR REFERENCES & VALUES ONLY                  */
/* ========================================================================= */

%typemap(in, typemap="in")
    gstlrn::ColID                (gstlrn::ColID *temp_colid = nullptr),
    const gstlrn::ColID&         (gstlrn::ColID *temp_colid = nullptr),
    gstlrn::ColID&&              (gstlrn::ColID *temp_colid = nullptr)
{
    SEXP obj = $input;
    SEXP target = obj;
    gstlrn::Id version = 0;

    /* ----------------------------------------------------------------- */
    /* Optional two-element list: (object, version)                      */
    /* ----------------------------------------------------------------- */

    if (TYPEOF(obj) == VECSXP && LENGTH(obj) == 2)
    {
        if (LENGTH(obj) != 2)
        {
            Rf_error("[ColID Typemap Error] ColID list must contain exactly two elements");
        }

        target = VECTOR_ELT(obj, 0);
        SEXP param = VECTOR_ELT(obj, 1);
        const SEXPTYPE param_type = TYPEOF(param);

        if (param_type == INTSXP && LENGTH(param) > 0)
        {
            version = static_cast<gstlrn::Id>(INTEGER(param)[0]);
        }
        else if (param_type == REALSXP && LENGTH(param) > 0)
        {
            version = static_cast<gstlrn::Id>(REAL(param)[0]);
        }
        else
        {
            Rf_error("[ColID Typemap Error] Second list element must be an integer version");
        }
    }

    const SEXPTYPE target_type = TYPEOF(target);

    /* ----------------------------------------------------------------- */
    /* 1. Column name (STRSXP)                                           */
    /* ----------------------------------------------------------------- */

    if (target_type == STRSXP && LENGTH(target) > 0)
    {
        const char *name = CHAR(STRING_ELT(target, 0));
        temp_colid = new gstlrn::ColID(std::string(name), version);
    }

    /* ----------------------------------------------------------------- */
    /* 2. Raw index (INTSXP / REALSXP)                                   */
    /* ----------------------------------------------------------------- */

    else if ((target_type == INTSXP || target_type == REALSXP) && LENGTH(target) > 0)
    {
        gstlrn::Id icol = (target_type == INTSXP)
            ? static_cast<gstlrn::Id>(INTEGER(target)[0])
            : static_cast<gstlrn::Id>(REAL(target)[0]);

        temp_colid = new gstlrn::ColID(icol, version);
    }

    /* ----------------------------------------------------------------- */
    /* 3. Direct external pointer (EXTPTRSXP)                            */
    /* ----------------------------------------------------------------- */

    else if (target_type == EXTPTRSXP)
    {
        static swig_type_info *type_ColID  = SWIG_TypeQuery("gstlrn::ColID *");
        static swig_type_info *type_RoleID = SWIG_TypeQuery("gstlrn::RoleID *");
        static swig_type_info *type_ERole  = SWIG_TypeQuery("gstlrn::ERole *");

        SEXP tag = R_ExternalPtrTag(target);
        swig_type_info *tag_type = nullptr;

        if (TYPEOF(tag) == EXTPTRSXP)
        {
            tag_type = reinterpret_cast<swig_type_info *>(R_ExternalPtrAddr(tag));
        }

        swig_type_info *expected_ColID  = type_ColID  ? type_ColID  : SWIGTYPE_p_gstlrn__ColID;
        swig_type_info *expected_RoleID = type_RoleID ? type_RoleID : SWIGTYPE_p_gstlrn__RoleID;
        swig_type_info *expected_ERole  = type_ERole  ? type_ERole  : SWIGTYPE_p_gstlrn__ERole;

        if (tag_type == expected_ColID)
        {
            gstlrn::ColID *col_ptr = nullptr;
            if (SWIG_IsOK(SWIG_ConvertPtr(target, reinterpret_cast<void **>(&col_ptr), expected_ColID, 0)) && col_ptr != nullptr)
            {
                temp_colid = (version != 0) ? new gstlrn::ColID(gstlrn::ColID::create(*col_ptr, version))
                                            : new gstlrn::ColID(*col_ptr);
            }
            else
            {
                Rf_error("[ColID Typemap Error] Cannot extract C++ pointer for ColID");
            }
        }
        else if (tag_type == expected_RoleID)
        {
            gstlrn::RoleID *roleid_ptr = nullptr;
            if (SWIG_IsOK(SWIG_ConvertPtr(target, reinterpret_cast<void **>(&roleid_ptr), expected_RoleID, 0)) && roleid_ptr != nullptr)
            {
                temp_colid = new gstlrn::ColID(*roleid_ptr, version);
            }
            else
            {
                Rf_error("[ColID Typemap Error] Cannot extract C++ pointer for RoleID");
            }
        }
        else if (tag_type == expected_ERole)
        {
            gstlrn::ERole *role_ptr = nullptr;
            if (SWIG_IsOK(SWIG_ConvertPtr(target, reinterpret_cast<void **>(&role_ptr), expected_ERole, 0)) && role_ptr != nullptr)
            {
                temp_colid = new gstlrn::ColID(*role_ptr, version);
            }
            else
            {
                Rf_error("[ColID Typemap Error] Cannot extract C++ pointer for ERole");
            }
        }
        else
        {
            Rf_error("[ColID Typemap Error] Unknown SWIG external pointer");
        }
    }

    /* ----------------------------------------------------------------- */
    /* 4. S4 ERole object                                                */
    /* ----------------------------------------------------------------- */

    else if (Rf_inherits(target, "_p_gstlrn__ERole"))
    {
        static SEXP sym_ref = Rf_install("ref");
        static swig_type_info *type_ERole = SWIG_TypeQuery("gstlrn::ERole *");

        SEXP ref = R_do_slot(target, sym_ref);
        gstlrn::ERole *role_ptr = nullptr;
        swig_type_info *expected = type_ERole ? type_ERole : SWIGTYPE_p_gstlrn__ERole;

        if (SWIG_IsOK(SWIG_ConvertPtr(ref, reinterpret_cast<void **>(&role_ptr), expected, 0)) && role_ptr != nullptr)
        {
            temp_colid = new gstlrn::ColID(*role_ptr, version);
        }
        else
        {
            Rf_error("[ColID Typemap Error] Cannot extract C++ pointer for ERole");
        }
    }

    /* ----------------------------------------------------------------- */
    /* 5. S4 ColID object                                                */
    /* ----------------------------------------------------------------- */

    else if (Rf_inherits(target, "_p_gstlrn__ColID"))
    {
        static SEXP sym_ref = Rf_install("ref");
        static swig_type_info *type_ColID = SWIG_TypeQuery("gstlrn::ColID *");

        SEXP ref = R_do_slot(target, sym_ref);
        gstlrn::ColID *col_ptr = nullptr;
        swig_type_info *expected = type_ColID ? type_ColID : SWIGTYPE_p_gstlrn__ColID;

        if (SWIG_IsOK(SWIG_ConvertPtr(ref, reinterpret_cast<void **>(&col_ptr), expected, 0)) && col_ptr != nullptr)
        {
            temp_colid = (version != 0) ? new gstlrn::ColID(gstlrn::ColID::create(*col_ptr, version))
                                        : new gstlrn::ColID(*col_ptr);
        }
        else
        {
            Rf_error("[ColID Typemap Error] Cannot extract C++ pointer for ColID");
        }
    }

    /* ----------------------------------------------------------------- */
    /* 6. S4 RoleID object                                               */
    /* ----------------------------------------------------------------- */

    else if (Rf_inherits(target, "_p_gstlrn__RoleID"))
    {
        static SEXP sym_ref = Rf_install("ref");
        static swig_type_info *type_RoleID = SWIG_TypeQuery("gstlrn::RoleID *");

        SEXP ref = R_do_slot(target, sym_ref);
        gstlrn::RoleID *roleid_ptr = nullptr;
        swig_type_info *expected = type_RoleID ? type_RoleID : SWIGTYPE_p_gstlrn__RoleID;

        if (SWIG_IsOK(SWIG_ConvertPtr(ref, reinterpret_cast<void **>(&roleid_ptr), expected, 0)) && roleid_ptr != nullptr)
        {
            temp_colid = new gstlrn::ColID(*roleid_ptr, version);
        }
        else
        {
            Rf_error("[ColID Typemap Error] Cannot extract C++ pointer for RoleID");
        }
    }

    /* ----------------------------------------------------------------- */
    /* 7. Failure fallback                                               */
    /* ----------------------------------------------------------------- */

    else
    {
        Rf_error("[ColID Typemap Error] Cannot convert R object to ColID");
    }

    $1 = temp_colid;
}

/* ========================================================================= */
/* 4. FREE ARG TYPEMAP (CLEANUP TEMPORARY INSTANCES)                         */
/* ========================================================================= */

%typemap(freearg)
    gstlrn::ColID,
    const gstlrn::ColID&,
    gstlrn::ColID&&
{
    if ($1)
        delete $1;
}

/* Strict neutralization for raw pointers (prevents double delete) */
%typemap(freearg) gstlrn::ColID * "";
