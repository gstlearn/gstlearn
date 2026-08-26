/***************************************************************************/
/*                                                                         */
/*  Typemap SWIG simple pour gstlrn::ColID&&                               */
/*                                                                         */
/*  VERSION V17                                                            */
/*                                                                         */
/***************************************************************************/

%typemap(in) gstlrn::ColID&& (gstlrn::ColID *temp_colid = nullptr)
{
    bool debug = false;

    SEXP obj = $input;
    SEXP target = obj;

    gstlrn::Id index = 0;
    gstlrn::Id version = 0;

    if (debug)
        Rprintf("\n[DEBUG Typemap ColID V17] --- Début conversion ---\n");

    /* ================================================================== */
    /* 1. Liste de deux éléments                                          */
    /* ================================================================== */

    if (TYPEOF(obj) == VECSXP && Rf_length(obj) == 2)
    {
        target = VECTOR_ELT(obj, 0);
        SEXP param = VECTOR_ELT(obj, 1);

        gstlrn::Id val =
            Rf_isInteger(param)
                ? static_cast<gstlrn::Id>(INTEGER(param)[0])
                : static_cast<gstlrn::Id>(REAL(param)[0]);

        /*
         * list(ERole, Id) représente un RoleID :
         *
         *     list(ERole_Z(), 0)
         *             -> RoleID(ERole_Z(), 0)
         *             -> ColID(RoleID, version=0)
         */

        if (Rf_inherits(target, "_p_gstlrn__ERole"))
        {
            index = val;

            if (debug)
                Rprintf(
                    "[DEBUG Typemap ColID V17] "
                    "Liste : ERole -> RoleID "
                    "(index=%d, version=%d)\n",
                    (int) index,
                    (int) version);
        }
        else
        {
            /*
             * Pour ColID ou RoleID, le second élément représente
             * la version.
             */
            version = val;

            if (debug)
                Rprintf(
                    "[DEBUG Typemap ColID V17] "
                    "Liste : autre objet (version=%d)\n",
                    (int) version);
        }
    }

    /* ================================================================== */
    /* 2. Diagnostic général                                             */
    /* ================================================================== */

    if (debug)
    {
        Rprintf(
            "[DEBUG Typemap ColID V17] target TYPEOF = %d\n",
            TYPEOF(target));

        Rprintf(
            "[DEBUG Typemap ColID V17] target isS4 = %d\n",
            Rf_isS4(target));

        Rprintf(
            "[DEBUG Typemap ColID V17] target is.object = %d\n",
            Rf_isObject(target));

        Rprintf(
            "[DEBUG Typemap ColID V17] inherits _p_gstlrn__ERole = %d\n",
            Rf_inherits(target, "_p_gstlrn__ERole"));

        Rprintf(
            "[DEBUG Typemap ColID V17] inherits _p_gstlrn__ColID = %d\n",
            Rf_inherits(target, "_p_gstlrn__ColID"));

        Rprintf(
            "[DEBUG Typemap ColID V17] inherits _p_gstlrn__RoleID = %d\n",
            Rf_inherits(target, "_p_gstlrn__RoleID"));
    }

    /* ================================================================== */
    /* 3. Nom de colonne                                                 */
    /* ================================================================== */

    if (Rf_isString(target) && Rf_length(target) > 0)
    {
        const char *name =
            CHAR(STRING_ELT(target, 0));

        if (debug)
            Rprintf(
                "[DEBUG Typemap ColID V17] "
                "Succès -> Par nom: '%s' (version=%d)\n",
                name,
                (int) version);

        temp_colid =
            new gstlrn::ColID(
                std::string(name),
                version);
    }

    /* ================================================================== */
    /* 4. External pointer direct                                        */
    /* ================================================================== */

    else if (TYPEOF(target) == EXTPTRSXP)
    {
        if (debug)
            Rprintf(
                "[DEBUG Typemap ColID V17] "
                ">>> BRANCHE EXTPTR DIRECTE <<<\n");

        SEXP tag = R_ExternalPtrTag(target);
        swig_type_info *tag_type = nullptr;

        if (TYPEOF(tag) == EXTPTRSXP)
            tag_type =
                (swig_type_info *) R_ExternalPtrAddr(tag);

        if (debug)
        {
            Rprintf(
                "[DEBUG Typemap ColID V17] "
                "R_ExternalPtrTag = %p\n",
                (void *) tag);

            Rprintf(
                "[DEBUG Typemap ColID V17] "
                "tag TYPEOF = %d\n",
                TYPEOF(tag));

            Rprintf(
                "[DEBUG Typemap ColID V17] "
                "tag_type = %p\n",
                (void *) tag_type);

            Rprintf(
                "[DEBUG Typemap ColID V17] "
                "SWIGTYPE_p_gstlrn__ERole = %p\n",
                (void *) SWIGTYPE_p_gstlrn__ERole);

            Rprintf(
                "[DEBUG Typemap ColID V17] "
                "SWIGTYPE_p_gstlrn__ColID = %p\n",
                (void *) SWIGTYPE_p_gstlrn__ColID);

            Rprintf(
                "[DEBUG Typemap ColID V17] "
                "SWIGTYPE_p_gstlrn__RoleID = %p\n",
                (void *) SWIGTYPE_p_gstlrn__RoleID);
        }

        /* -------------------------------------------------------------- */
        /* 4a. ColID                                                      */
        /* -------------------------------------------------------------- */

        if (tag_type == SWIGTYPE_p_gstlrn__ColID)
        {
            gstlrn::ColID *col_ptr = nullptr;

            int res =
                SWIG_ConvertPtr(
                    target,
                    (void **) &col_ptr,
                    SWIGTYPE_p_gstlrn__ColID,
                    0);

            if (debug)
                Rprintf(
                    "[DEBUG Typemap ColID V17] "
                    "TAG = ColID ; SWIG_ConvertPtr res=%d ptr=%p\n",
                    res,
                    (void *) col_ptr);

            if (SWIG_IsOK(res) && col_ptr != nullptr)
            {
                temp_colid =
                    new gstlrn::ColID(
                        gstlrn::ColID::create(
                            *col_ptr,
                            version));

                if (debug)
                    Rprintf(
                        "[DEBUG Typemap ColID V17] "
                        "Succès -> ColID direct (version=%d)\n",
                        (int) version);
            }
            else
            {
                Rf_error(
                    "Impossible de récupérer le pointeur C++ de ColID");
            }
        }

        /* -------------------------------------------------------------- */
        /* 4b. RoleID                                                     */
        /* -------------------------------------------------------------- */

        else if (tag_type == SWIGTYPE_p_gstlrn__RoleID)
        {
            gstlrn::RoleID *roleid_ptr = nullptr;

            int res =
                SWIG_ConvertPtr(
                    target,
                    (void **) &roleid_ptr,
                    SWIGTYPE_p_gstlrn__RoleID,
                    0);

            if (debug)
                Rprintf(
                    "[DEBUG Typemap ColID V17] "
                    "TAG = RoleID ; SWIG_ConvertPtr res=%d ptr=%p\n",
                    res,
                    (void *) roleid_ptr);

            if (SWIG_IsOK(res) && roleid_ptr != nullptr)
            {
                temp_colid =
                    new gstlrn::ColID(
                        *roleid_ptr,
                        version);

                if (debug)
                    Rprintf(
                        "[DEBUG Typemap ColID V17] "
                        "Succès -> RoleID direct (version=%d)\n",
                        (int) version);
            }
            else
            {
                Rf_error(
                    "Impossible de récupérer le pointeur C++ de RoleID");
            }
        }

        /* -------------------------------------------------------------- */
        /* 4c. ERole                                                      */
        /* -------------------------------------------------------------- */

        else if (tag_type == SWIGTYPE_p_gstlrn__ERole)
        {
            gstlrn::ERole *role_ptr = nullptr;

            int res =
                SWIG_ConvertPtr(
                    target,
                    (void **) &role_ptr,
                    SWIGTYPE_p_gstlrn__ERole,
                    0);

            if (debug)
                Rprintf(
                    "[DEBUG Typemap ColID V17] "
                    "TAG = ERole ; SWIG_ConvertPtr res=%d ptr=%p\n",
                    res,
                    (void *) role_ptr);

            if (SWIG_IsOK(res) && role_ptr != nullptr)
            {
                temp_colid =
                    new gstlrn::ColID(
                        *role_ptr,
                        index,
                        version);

                if (debug)
                    Rprintf(
                        "[DEBUG Typemap ColID V17] "
                        "Succès -> ERole direct "
                        "(index=%d version=%d)\n",
                        (int) index,
                        (int) version);
            }
            else
            {
                Rf_error(
                    "Impossible de récupérer le pointeur C++ de ERole");
            }
        }

        /* -------------------------------------------------------------- */
        /* 4d. Tag inconnu                                                */
        /* -------------------------------------------------------------- */

        else
        {
            Rf_error(
                "External pointer SWIG inconnu : "
                "impossible de déterminer ColID, RoleID ou ERole");
        }
    }

    /* ================================================================== */
    /* 5. ERole S4                                                       */
    /* ================================================================== */

    else if (Rf_inherits(target, "_p_gstlrn__ERole"))
    {
        if (debug)
            Rprintf(
                "[DEBUG Typemap ColID V17] "
                ">>> BRANCHE S4 ERole <<<\n");

        SEXP ref =
            R_do_slot(
                target,
                Rf_install("ref"));

        gstlrn::ERole *role_ptr = nullptr;

        int res =
            SWIG_ConvertPtr(
                ref,
                (void **) &role_ptr,
                SWIGTYPE_p_gstlrn__ERole,
                0);

        if (debug)
            Rprintf(
                "[DEBUG Typemap ColID V17] "
                "SWIG_ConvertPtr ERole = %d ptr=%p\n",
                res,
                (void *) role_ptr);

        if (SWIG_IsOK(res) && role_ptr != nullptr)
        {
            temp_colid =
                new gstlrn::ColID(
                    *role_ptr,
                    index,
                    version);

            if (debug)
                Rprintf(
                    "[DEBUG Typemap ColID V17] "
                    "Succès -> ERole "
                    "(index=%d version=%d)\n",
                    (int) index,
                    (int) version);
        }
        else
        {
            Rf_error(
                "Impossible de récupérer le pointeur C++ de ERole");
        }
    }

    /* ================================================================== */
    /* 6. ColID S4                                                       */
    /* ================================================================== */

    else if (Rf_inherits(target, "_p_gstlrn__ColID"))
    {
        if (debug)
            Rprintf(
                "[DEBUG Typemap ColID V17] "
                ">>> BRANCHE S4 ColID <<<\n");

        SEXP ref =
            R_do_slot(
                target,
                Rf_install("ref"));

        gstlrn::ColID *col_ptr = nullptr;

        int res =
            SWIG_ConvertPtr(
                ref,
                (void **) &col_ptr,
                SWIGTYPE_p_gstlrn__ColID,
                0);

        if (debug)
            Rprintf(
                "[DEBUG Typemap ColID V17] "
                "SWIG_ConvertPtr ColID = %d ptr=%p\n",
                res,
                (void *) col_ptr);

        if (SWIG_IsOK(res) && col_ptr != nullptr)
        {
            temp_colid =
                new gstlrn::ColID(
                    gstlrn::ColID::create(
                        *col_ptr,
                        version));
        }
        else
        {
            Rf_error(
                "Impossible de récupérer le pointeur C++ de ColID");
        }
    }

    /* ================================================================== */
    /* 7. RoleID S4                                                      */
    /* ================================================================== */

    else if (Rf_inherits(target, "_p_gstlrn__RoleID"))
    {
        if (debug)
            Rprintf(
                "[DEBUG Typemap ColID V17] "
                ">>> BRANCHE S4 RoleID <<<\n");

        SEXP ref =
            R_do_slot(
                target,
                Rf_install("ref"));

        gstlrn::RoleID *roleid_ptr = nullptr;

        int res =
            SWIG_ConvertPtr(
                ref,
                (void **) &roleid_ptr,
                SWIGTYPE_p_gstlrn__RoleID,
                0);

        if (debug)
            Rprintf(
                "[DEBUG Typemap ColID V17] "
                "SWIG_ConvertPtr RoleID = %d ptr=%p\n",
                res,
                (void *) roleid_ptr);

        if (SWIG_IsOK(res) && roleid_ptr != nullptr)
        {
            temp_colid =
                new gstlrn::ColID(
                    *roleid_ptr,
                    version);

            if (debug)
                Rprintf(
                    "[DEBUG Typemap ColID V17] "
                    "Succès -> RoleID (version=%d)\n",
                    (int) version);
        }
        else
        {
            Rf_error(
                "Impossible de récupérer le pointeur C++ de RoleID");
        }
    }

    /* ================================================================== */
    /* 8. Index brut                                                     */
    /* ================================================================== */

    else if (Rf_isInteger(target) || Rf_isReal(target))
    {
        gstlrn::Id icol =
            Rf_isInteger(target)
                ? static_cast<gstlrn::Id>(INTEGER(target)[0])
                : static_cast<gstlrn::Id>(REAL(target)[0]);

        temp_colid =
            new gstlrn::ColID(
                icol,
                version);

        if (debug)
            Rprintf(
                "[DEBUG Typemap ColID V17] "
                "Par index de colonne: %d (version=%d)\n",
                (int) icol,
                (int) version);
    }

    /* ================================================================== */
    /* 9. Échec                                                          */
    /* ================================================================== */

    else
    {
        Rf_error(
            "Impossible de convertir l'objet R en ColID");
    }

    $1 = temp_colid;
}

%typemap(freearg) gstlrn::ColID&&
{
    if ($1)
        delete $1;
}
