/***************************************************************************/
/*                                                                         */
/*  R typemaps for gstlrn::ColID (Version #10 - In-file Toggle)            */
/*                                                                         */
/***************************************************************************/

%{
#include <R.h>
#include <Rinternals.h>

/* Set to false to disable all debug prints for non-regression tests */
static const bool COLID_DEBUG_ENABLED = false;

#define COLID_DEBUG_LOG(...) \
    do { \
        if (COLID_DEBUG_ENABLED) { \
            Rprintf(__VA_ARGS__); \
        } \
    } while(0)
%}

/* ========================================================================= */
/* 1. BYPASS SWIG R AUTOMATIC S4 COERCION & AS() CONVERSIONS                 */
/* ========================================================================= */

%typemap(rtype)      gstlrn::ColID, const gstlrn::ColID&, gstlrn::ColID&&, gstlrn::ColID* "ANY";
%typemap(scoercein)  gstlrn::ColID, const gstlrn::ColID&, gstlrn::ColID&&, gstlrn::ColID* "$input";
%typemap(sclass)     gstlrn::ColID, const gstlrn::ColID&, gstlrn::ColID&&, gstlrn::ColID* "ANY";

/* ========================================================================= */
/* 2. UNIVERSAL TYPECHECK FOR SWIG OVERLOAD RESOLUTION                       */
/* ========================================================================= */

%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER)
    gstlrn::ColID,
    const gstlrn::ColID&,
    gstlrn::ColID&&,
    gstlrn::ColID*
{
    SEXP obj = $input;
    if (TYPEOF(obj) == VECSXP && LENGTH(obj) == 2)
    {
        obj = VECTOR_ELT(obj, 0);
    }

    const SEXPTYPE t = TYPEOF(obj);
    COLID_DEBUG_LOG("[ColID Typemap V10] TYPECHECK: SEXPTYPE = %d\n", static_cast<int>(t));

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
/* 3. INPUT CONVERSION TYPEMAP WITH EXPLICIT TEMPORARY DECLARATION          */
/* ========================================================================= */

%typemap(in, typemap="in")
    gstlrn::ColID                (gstlrn::ColID *temp_colid = nullptr),
    const gstlrn::ColID&         (gstlrn::ColID *temp_colid = nullptr),
    gstlrn::ColID&&              (gstlrn::ColID *temp_colid = nullptr),
    gstlrn::ColID*               (gstlrn::ColID *temp_colid = nullptr)
{
    SEXP obj = $input;
    SEXP target = obj;

    gstlrn::Id index = 0;
    gstlrn::Id version = 0;

    COLID_DEBUG_LOG("[ColID Typemap V10] IN: Entering conversion. TYPEOF = %d\n", static_cast<int>(TYPEOF(obj)));

    /* ----------------------------------------------------------------- */
    /* 1. Two-element list                                               */
    /* ----------------------------------------------------------------- */

    if (TYPEOF(obj) == VECSXP && LENGTH(obj) == 2)
    {
        COLID_DEBUG_LOG("[ColID Typemap V10] Matched: 2-element VECSXP (list)\n");
        target = VECTOR_ELT(obj, 0);
        SEXP param = VECTOR_ELT(obj, 1);

        const SEXPTYPE param_type = TYPEOF(param);
        gstlrn::Id val = 0;

        if (param_type == INTSXP && LENGTH(param) > 0)
        {
            val = static_cast<gstlrn::Id>(INTEGER(param)[0]);
        }
        else if (param_type == REALSXP && LENGTH(param) > 0)
        {
            val = static_cast<gstlrn::Id>(REAL(param)[0]);
        }

        if (Rf_inherits(target, "_p_gstlrn__ERole"))
        {
            index = val;
        }
        else
        {
            version = val;
        }
    }

    /* ----------------------------------------------------------------- */
    /* 2. Column name (STRSXP)                                           */
    /* ----------------------------------------------------------------- */

    const SEXPTYPE target_type = TYPEOF(target);

    if (target_type == STRSXP && LENGTH(target) > 0)
    {
        const char *name = CHAR(STRING_ELT(target, 0));
        COLID_DEBUG_LOG("[ColID Typemap V10] Matched: STRSXP (String) -> '%s' (version=%d)\n", name, static_cast<int>(version));

        temp_colid = new gstlrn::ColID(
            std::string(name),
            version);
    }

    /* ----------------------------------------------------------------- */
    /* 3. Direct external pointer                                        */
    /* ----------------------------------------------------------------- */

    else if (target_type == EXTPTRSXP)
    {
        COLID_DEBUG_LOG("[ColID Typemap V10] Matched: EXTPTRSXP\n");

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
            int res = SWIG_ConvertPtr(target, reinterpret_cast<void **>(&col_ptr), expected_ColID, 0);

            if (SWIG_IsOK(res) && col_ptr != nullptr)
            {
                temp_colid = new gstlrn::ColID(
                    gstlrn::ColID::create(*col_ptr, version));
            }
            else
            {
                Rf_error("[ColID Typemap Error] Cannot extract C++ pointer for ColID");
            }
        }
        else if (tag_type == expected_RoleID)
        {
            gstlrn::RoleID *roleid_ptr = nullptr;
            int res = SWIG_ConvertPtr(target, reinterpret_cast<void **>(&roleid_ptr), expected_RoleID, 0);

            if (SWIG_IsOK(res) && roleid_ptr != nullptr)
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
            int res = SWIG_ConvertPtr(target, reinterpret_cast<void **>(&role_ptr), expected_ERole, 0);

            if (SWIG_IsOK(res) && role_ptr != nullptr)
            {
                temp_colid = new gstlrn::ColID(*role_ptr, index, version);
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
    /* 4. S4 ERole                                                       */
    /* ----------------------------------------------------------------- */

    else if (Rf_inherits(target, "_p_gstlrn__ERole"))
    {
        COLID_DEBUG_LOG("[ColID Typemap V10] Matched: S4 ERole\n");
        static SEXP sym_ref = Rf_install("ref");
        static swig_type_info *type_ERole = SWIG_TypeQuery("gstlrn::ERole *");

        SEXP ref = R_do_slot(target, sym_ref);
        gstlrn::ERole *role_ptr = nullptr;

        swig_type_info *expected = type_ERole ? type_ERole : SWIGTYPE_p_gstlrn__ERole;
        int res = SWIG_ConvertPtr(ref, reinterpret_cast<void **>(&role_ptr), expected, 0);

        if (SWIG_IsOK(res) && role_ptr != nullptr)
        {
            temp_colid = new gstlrn::ColID(*role_ptr, index, version);
        }
        else
        {
            Rf_error("[ColID Typemap Error] Cannot extract C++ pointer for ERole");
        }
    }

    /* ----------------------------------------------------------------- */
    /* 5. S4 ColID                                                       */
    /* ----------------------------------------------------------------- */

    else if (Rf_inherits(target, "_p_gstlrn__ColID"))
    {
        COLID_DEBUG_LOG("[ColID Typemap V10] Matched: S4 ColID\n");
        static SEXP sym_ref = Rf_install("ref");
        static swig_type_info *type_ColID = SWIG_TypeQuery("gstlrn::ColID *");

        SEXP ref = R_do_slot(target, sym_ref);
        gstlrn::ColID *col_ptr = nullptr;

        swig_type_info *expected = type_ColID ? type_ColID : SWIGTYPE_p_gstlrn__ColID;
        int res = SWIG_ConvertPtr(ref, reinterpret_cast<void **>(&col_ptr), expected, 0);

        if (SWIG_IsOK(res) && col_ptr != nullptr)
        {
            temp_colid = new gstlrn::ColID(
                gstlrn::ColID::create(*col_ptr, version));
        }
        else
        {
            Rf_error("[ColID Typemap Error] Cannot extract C++ pointer for ColID");
        }
    }

    /* ----------------------------------------------------------------- */
    /* 6. S4 RoleID                                                      */
    /* ----------------------------------------------------------------- */

    else if (Rf_inherits(target, "_p_gstlrn__RoleID"))
    {
        COLID_DEBUG_LOG("[ColID Typemap V10] Matched: S4 RoleID\n");
        static SEXP sym_ref = Rf_install("ref");
        static swig_type_info *type_RoleID = SWIG_TypeQuery("gstlrn::RoleID *");

        SEXP ref = R_do_slot(target, sym_ref);
        gstlrn::RoleID *roleid_ptr = nullptr;

        swig_type_info *expected = type_RoleID ? type_RoleID : SWIGTYPE_p_gstlrn__RoleID;
        int res = SWIG_ConvertPtr(ref, reinterpret_cast<void **>(&roleid_ptr), expected, 0);

        if (SWIG_IsOK(res) && roleid_ptr != nullptr)
        {
            temp_colid = new gstlrn::ColID(*roleid_ptr, version);
        }
        else
        {
            Rf_error("[ColID Typemap Error] Cannot extract C++ pointer for RoleID");
        }
    }

    /* ----------------------------------------------------------------- */
    /* 7. Raw index (INTSXP / REALSXP)                                   */
    /* ----------------------------------------------------------------- */

    else if ((target_type == INTSXP || target_type == REALSXP) && LENGTH(target) > 0)
    {
        gstlrn::Id icol = (target_type == INTSXP)
            ? static_cast<gstlrn::Id>(INTEGER(target)[0])
            : static_cast<gstlrn::Id>(REAL(target)[0]);

        COLID_DEBUG_LOG("[ColID Typemap V10] Matched: Integer/Real index -> %d (version=%d)\n", static_cast<int>(icol), static_cast<int>(version));

        temp_colid = new gstlrn::ColID(icol, version);
    }

    /* ----------------------------------------------------------------- */
    /* 8. Failure                                                        */
    /* ----------------------------------------------------------------- */

    else
    {
        COLID_DEBUG_LOG("[ColID Typemap V10] Conversion FAILED for target_type = %d\n", static_cast<int>(target_type));
        Rf_error("[ColID Typemap Error] Cannot convert R object to ColID");
    }

    $1 = temp_colid;
}

/* ========================================================================= */
/* 4. FREE ARG TYPEMAP                                                       */
/* ========================================================================= */

%typemap(freearg)
    gstlrn::ColID,
    const gstlrn::ColID&,
    gstlrn::ColID&&,
    gstlrn::ColID*
{
    COLID_DEBUG_LOG("[ColID Typemap V10] FREEARG: Cleaning up temporary ColID\n");
    if ($1)
        delete $1;
}
