/***************************************************************************/
/*                                                                         */
/*  Typemap SWIG simple pour gstlrn::ColID&&                               */
/*                                                                         */
/*  VERSION V17 - OPTIMIZED                                                */
/*                                                                         */
/***************************************************************************/

%typemap(in) gstlrn::ColID&& (gstlrn::ColID *temp_colid = nullptr)
{
    SEXP obj = $input;
    SEXP target = obj;

    gstlrn::Id index = 0;
    gstlrn::Id version = 0;

    /* ================================================================== */
    /* 1. Liste de deux éléments                                          */
    /* ================================================================== */

    if (TYPEOF(obj) == VECSXP && LENGTH(obj) == 2)
    {
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

        /*
         * list(ERole, Id) représente un RoleID :
         *     list(ERole_Z(), 0) -> RoleID(ERole_Z(), 0) -> ColID(RoleID, version=0)
         */
        if (Rf_inherits(target, "_p_gstlrn__ERole"))
        {
            index = val;
        }
        else
        {
            /* Pour ColID ou RoleID, le second élément représente la version */
            version = val;
        }
    }

    /* ================================================================== */
    /* 2. Diagnostic général (ignoré en release)                        */
    /* ================================================================== */

    /* ================================================================== */
    /* 3. Nom de colonne                                                 */
    /* ================================================================== */

    const SEXPTYPE target_type = TYPEOF(target);

    if (target_type == STRSXP && LENGTH(target) > 0)
    {
        const char *name = CHAR(STRING_ELT(target, 0));

        temp_colid = new gstlrn::ColID(
            std::string(name),
            version);
    }

    /* ================================================================== */
    /* 4. External pointer direct                                        */
    /* ================================================================== */

    else if (target_type == EXTPTRSXP)
    {
        // Cache SWIG type descriptors to avoid repeated lookups
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

        /* -------------------------------------------------------------- */
        /* 4a. ColID                                                      */
        /* -------------------------------------------------------------- */

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
                Rf_error("Impossible de récupérer le pointeur C++ de ColID");
            }
        }

        /* -------------------------------------------------------------- */
        /* 4b. RoleID                                                     */
        /* -------------------------------------------------------------- */

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
                Rf_error("Impossible de récupérer le pointeur C++ de RoleID");
            }
        }

        /* -------------------------------------------------------------- */
        /* 4c. ERole                                                      */
        /* -------------------------------------------------------------- */

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
                Rf_error("Impossible de récupérer le pointeur C++ de ERole");
            }
        }

        /* -------------------------------------------------------------- */
        /* 4d. Tag inconnu                                                */
        /* -------------------------------------------------------------- */

        else
        {
            Rf_error("External pointer SWIG inconnu : impossible de déterminer ColID, RoleID ou ERole");
        }
    }

    /* ================================================================== */
    /* 5. ERole S4                                                       */
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
            temp_colid = new gstlrn::ColID(*role_ptr, index, version);
        }
        else
        {
            Rf_error("Impossible de récupérer le pointeur C++ de ERole");
        }
    }

    /* ================================================================== */
    /* 6. ColID S4                                                       */
    /* ================================================================== */

    else if (Rf_inherits(target, "_p_gstlrn__ColID"))
    {
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
            Rf_error("Impossible de récupérer le pointeur C++ de ColID");
        }
    }

    /* ================================================================== */
    /* 7. RoleID S4                                                      */
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
            temp_colid = new gstlrn::ColID(*roleid_ptr, version);
        }
        else
        {
            Rf_error("Impossible de récupérer le pointeur C++ de RoleID");
        }
    }

    /* ================================================================== */
    /* 8. Index brut                                                     */
    /* ================================================================== */

    else if ((target_type == INTSXP || target_type == REALSXP) && LENGTH(target) > 0)
    {
        gstlrn::Id icol = (target_type == INTSXP)
            ? static_cast<gstlrn::Id>(INTEGER(target)[0])
            : static_cast<gstlrn::Id>(REAL(target)[0]);

        temp_colid = new gstlrn::ColID(icol, version);
    }

    /* ================================================================== */
    /* 9. Échec                                                          */
    /* ================================================================== */

    else
    {
        Rf_error("Impossible de convertir l'objet R en ColID");
    }

    $1 = temp_colid;
}

%typemap(freearg) gstlrn::ColID&&
{
    if ($1)
        delete $1;
}
