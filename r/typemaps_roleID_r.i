/***************************************************************************/
/*                                                                         */
/*  R typemap for const gstlrn::RoleID&                                    */
/*                                                                         */
/*  Accepted R objects:                                                    */
/*                                                                         */
/*      RoleID(...)               -> RoleID instance                       */
/*      list(RoleID(...), index)  -> RoleID(RoleID(...), index)            */
/*                                                                         */
/*      ERole                     -> RoleID(ERole, index=0)                */
/*      list(ERole, index)        -> RoleID(ERole, index)                  */
/*                                                                         */
/***************************************************************************/

%typemap(in) const gstlrn::RoleID& (gstlrn::RoleID temp_roleid)
{
    SEXP obj = $input;
    SEXP target = obj;

    gstlrn::Id index = 0;

    /* ================================================================== */
    /* 1. Two-element list: list(object, index)                           */
    /* ================================================================== */

    if (TYPEOF(obj) == VECSXP && LENGTH(obj) == 2)
    {
        target = VECTOR_ELT(obj, 0);
        SEXP param = VECTOR_ELT(obj, 1);

        const SEXPTYPE param_type = TYPEOF(param);

        if (param_type == INTSXP && LENGTH(param) > 0)
        {
            index = static_cast<gstlrn::Id>(INTEGER(param)[0]);
        }
        else if (param_type == REALSXP && LENGTH(param) > 0)
        {
            index = static_cast<gstlrn::Id>(REAL(param)[0]);
        }
        else
        {
            Rf_error("Second list element must be an integer index");
        }
    }

    /* ================================================================== */
    /* 2. Direct external pointer                                         */
    /* ================================================================== */

    const SEXPTYPE target_type = TYPEOF(target);

    if (target_type == EXTPTRSXP)
    {
        // Cache SWIG type descriptors to avoid repeated lookups
        static swig_type_info *type_RoleID = SWIG_TypeQuery("gstlrn::RoleID *");
        static swig_type_info *type_ERole  = SWIG_TypeQuery("gstlrn::ERole *");

        SEXP tag = R_ExternalPtrTag(target);
        swig_type_info *tag_type = nullptr;

        if (TYPEOF(tag) == EXTPTRSXP)
        {
            tag_type = reinterpret_cast<swig_type_info *>(R_ExternalPtrAddr(tag));
        }

        swig_type_info *expected_RoleID = type_RoleID ? type_RoleID : SWIGTYPE_p_gstlrn__RoleID;
        swig_type_info *expected_ERole  = type_ERole  ? type_ERole  : SWIGTYPE_p_gstlrn__ERole;

        /* -------------------------------------------------------------- */
        /* 2a. RoleID                                                     */
        /* -------------------------------------------------------------- */

        if (tag_type == expected_RoleID)
        {
            gstlrn::RoleID *roleid_ptr = nullptr;
            int res = SWIG_ConvertPtr(target, reinterpret_cast<void **>(&roleid_ptr), expected_RoleID, 0);

            if (SWIG_IsOK(res) && roleid_ptr != nullptr)
            {
                temp_roleid = gstlrn::RoleID(*roleid_ptr, index);
                $1 = &temp_roleid;
            }
            else
            {
                Rf_error("Cannot extract C++ pointer for RoleID");
            }
        }

        /* -------------------------------------------------------------- */
        /* 2b. ERole                                                      */
        /* -------------------------------------------------------------- */

        else if (tag_type == expected_ERole)
        {
            gstlrn::ERole *role_ptr = nullptr;
            int res = SWIG_ConvertPtr(target, reinterpret_cast<void **>(&role_ptr), expected_ERole, 0);

            if (SWIG_IsOK(res) && role_ptr != nullptr)
            {
                temp_roleid = gstlrn::RoleID(*role_ptr, index);
                $1 = &temp_roleid;
            }
            else
            {
                Rf_error("Cannot extract C++ pointer for ERole");
            }
        }

        /* -------------------------------------------------------------- */
        /* 2c. Unknown tag                                                */
        /* -------------------------------------------------------------- */

        else
        {
            Rf_error("Unknown SWIG external pointer: unable to determine RoleID or ERole");
        }
    }

    /* ================================================================== */
    /* 3. S4 ERole                                                       */
    /* ================================================================== */

    else if (Rf_inherits(target, "_p_gstlrn__ERole"))
    {
        static SEXP sym_ref = Rf_install("ref");
        static swig_type_info *type_ERole = SWIG_TypeQuery("gstlrn::ERole *");

        SEXP ref = R_do_slot(target, sym_ref);
        gstlrn::ERole *role_ptr = nullptr;

        swig_type_info *expected = type_ERole ? type_ERole : SWIGTYPE_p_gstlrn__ERole;
        int res = SWIG_ConvertPtr(ref, reinterpret_cast<void **>(&role_ptr), expected, 0);

        if (SWIG_IsOK(res) && role_ptr != nullptr)
        {
            temp_roleid = gstlrn::RoleID(*role_ptr, index);
            $1 = &temp_roleid;
        }
        else
        {
            Rf_error("Cannot extract C++ pointer for ERole");
        }
    }

    /* ================================================================== */
    /* 4. S4 RoleID                                                      */
    /* ================================================================== */

    else if (Rf_inherits(target, "_p_gstlrn__RoleID"))
    {
        static SEXP sym_ref = Rf_install("ref");
        static swig_type_info *type_RoleID = SWIG_TypeQuery("gstlrn::RoleID *");

        SEXP ref = R_do_slot(target, sym_ref);
        gstlrn::RoleID *roleid_ptr = nullptr;

        swig_type_info *expected = type_RoleID ? type_RoleID : SWIGTYPE_p_gstlrn__RoleID;
        int res = SWIG_ConvertPtr(ref, reinterpret_cast<void **>(&roleid_ptr), expected, 0);

        if (SWIG_IsOK(res) && roleid_ptr != nullptr)
        {
            temp_roleid = gstlrn::RoleID(*roleid_ptr, index);
            $1 = &temp_roleid;
        }
        else
        {
            Rf_error("Cannot extract C++ pointer for RoleID");
        }
    }

    /* ================================================================== */
    /* 5. Failure                                                         */
    /* ================================================================== */

    else
    {
        Rf_error("Cannot convert R object to RoleID");
    }
}
