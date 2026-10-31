/***************************************************************************/
/*                                                                         */
/*  R typemap for gstlrn::MultID&& and const gstlrn::MultID&               */
/*                                                                         */
/*  Accepted R objects:                                                    */
/*                                                                         */
/*      "name"                    -> MultID("name")                        */
/*      c("name1", "name2")       -> MultID(VectorString)                  */
/*      list("name1", "name2")    -> MultID(VectorString)                  */
/*                                                                         */
/*      icol                      -> MultID(Id icol)                       */
/*      c(icol1, icol2)           -> MultID(VectorInt)                     */
/*      list(icol1, icol2)        -> MultID(VectorInt)                     */
/*                                                                         */
/*      ERole                     -> MultID(ERole)                         */
/*      MultID(...)               -> copy constructor                      */
/*                                                                         */
/***************************************************************************/

%typemap(in) gstlrn::MultID&& (gstlrn::MultID *temp_multid = nullptr), const gstlrn::MultID& (gstlrn::MultID *temp_multid = nullptr)
{
    SEXP obj = $input;
    const SEXPTYPE obj_type = TYPEOF(obj);

    /* ================================================================== */
    /* 1. Vector of character strings: STRSXP                             */
    /* ================================================================== */

    if (obj_type == STRSXP && LENGTH(obj) > 0)
    {
        R_xlen_t len = LENGTH(obj);
        if (len == 1)
        {
            const char *str = CHAR(STRING_ELT(obj, 0));
            temp_multid = new gstlrn::MultID(std::string(str));
        }
        else
        {
            gstlrn::VectorString names;
            names.reserve(static_cast<size_t>(len));
            for (R_xlen_t i = 0; i < len; ++i)
            {
                names.push_back(std::string(CHAR(STRING_ELT(obj, i))));
            }
            temp_multid = new gstlrn::MultID(names);
        }
    }

    /* ================================================================== */
    /* 2. Numeric vector (Integers or Reals/Doubles): INTSXP / REALSXP   */
    /* ================================================================== */

    else if ((obj_type == INTSXP || obj_type == REALSXP) && LENGTH(obj) > 0)
    {
        R_xlen_t len = LENGTH(obj);

        if (len == 1)
        {
            gstlrn::Id icol = (obj_type == INTSXP)
                ? static_cast<gstlrn::Id>(INTEGER(obj)[0])
                : static_cast<gstlrn::Id>(REAL(obj)[0]);
            temp_multid = new gstlrn::MultID(icol);
        }
        else
        {
            gstlrn::VectorInt icols;
            icols.reserve(static_cast<size_t>(len));

            if (obj_type == INTSXP)
            {
                int *data = INTEGER(obj);
                for (R_xlen_t i = 0; i < len; ++i)
                {
                    icols.push_back(static_cast<gstlrn::Id>(data[i]));
                }
            }
            else
            {
                double *data = REAL(obj);
                for (R_xlen_t i = 0; i < len; ++i)
                {
                    icols.push_back(static_cast<gstlrn::Id>(data[i]));
                }
            }

            temp_multid = new gstlrn::MultID(icols);
        }
    }

    /* ================================================================== */
    /* 3. List R: VECSXP                                                  */
    /* ================================================================== */

    else if (obj_type == VECSXP && LENGTH(obj) > 0)
    {
        R_xlen_t len = LENGTH(obj);
        SEXP first_elt = VECTOR_ELT(obj, 0);
        SEXPTYPE first_type = TYPEOF(first_elt);

        if (first_type == STRSXP)
        {
            gstlrn::VectorString names;
            names.reserve(static_cast<size_t>(len));
            for (R_xlen_t i = 0; i < len; ++i)
            {
                SEXP item = VECTOR_ELT(obj, i);
                if (TYPEOF(item) != STRSXP || LENGTH(item) == 0)
                {
                    Rf_error("All list elements for MultID names must be non-empty character strings");
                }
                names.push_back(std::string(CHAR(STRING_ELT(item, 0))));
            }
            temp_multid = new gstlrn::MultID(names);
        }
        else if (first_type == INTSXP || first_type == REALSXP)
        {
            gstlrn::VectorInt icols;
            icols.reserve(static_cast<size_t>(len));
            for (R_xlen_t i = 0; i < len; ++i)
            {
                SEXP item = VECTOR_ELT(obj, i);
                SEXPTYPE item_type = TYPEOF(item);
                if ((item_type != INTSXP && item_type != REALSXP) || LENGTH(item) == 0)
                {
                    Rf_error("All list elements for MultID column indices must be numeric");
                }
                gstlrn::Id val = (item_type == INTSXP)
                    ? static_cast<gstlrn::Id>(INTEGER(item)[0])
                    : static_cast<gstlrn::Id>(REAL(item)[0]);
                icols.push_back(val);
            }
            temp_multid = new gstlrn::MultID(icols);
        }
        else
        {
            Rf_error("List for MultID must contain either strings or numeric indices");
        }
    }

    /* ================================================================== */
    /* 4. Direct External Pointer                                         */
    /* ================================================================== */

    else if (obj_type == EXTPTRSXP)
    {
        static swig_type_info *type_MultID = SWIG_TypeQuery("gstlrn::MultID *");
        static swig_type_info *type_ERole  = SWIG_TypeQuery("gstlrn::ERole *");

        SEXP tag = R_ExternalPtrTag(obj);
        swig_type_info *tag_type = nullptr;

        if (TYPEOF(tag) == EXTPTRSXP)
        {
            tag_type = reinterpret_cast<swig_type_info *>(R_ExternalPtrAddr(tag));
        }

        swig_type_info *expected_MultID = type_MultID ? type_MultID : SWIGTYPE_p_gstlrn__MultID;
        swig_type_info *expected_ERole  = type_ERole  ? type_ERole  : SWIGTYPE_p_gstlrn__ERole;

        if (tag_type == expected_MultID)
        {
            gstlrn::MultID *mult_ptr = nullptr;
            int res = SWIG_ConvertPtr(obj, reinterpret_cast<void **>(&mult_ptr), expected_MultID, 0);

            if (SWIG_IsOK(res) && mult_ptr != nullptr)
            {
                temp_multid = new gstlrn::MultID(*mult_ptr);
            }
            else
            {
                Rf_error("Cannot extract C++ pointer for MultID");
            }
        }
        else if (tag_type == expected_ERole)
        {
            gstlrn::ERole *role_ptr = nullptr;
            int res = SWIG_ConvertPtr(obj, reinterpret_cast<void **>(&role_ptr), expected_ERole, 0);

            if (SWIG_IsOK(res) && role_ptr != nullptr)
            {
                temp_multid = new gstlrn::MultID(*role_ptr);
            }
            else
            {
                Rf_error("Cannot extract C++ pointer for ERole");
            }
        }
        else
        {
            Rf_error("Unknown SWIG external pointer: unable to determine MultID or ERole");
        }
    }

    /* ================================================================== */
    /* 5. S4 ERole Object                                                */
    /* ================================================================== */

    else if (Rf_inherits(obj, "_p_gstlrn__ERole"))
    {
        static SEXP sym_ref = Rf_install("ref");
        static swig_type_info *type_ERole = SWIG_TypeQuery("gstlrn::ERole *");

        SEXP ref = R_do_slot(obj, sym_ref);
        gstlrn::ERole *role_ptr = nullptr;

        swig_type_info *expected = type_ERole ? type_ERole : SWIGTYPE_p_gstlrn__ERole;
        int res = SWIG_ConvertPtr(ref, reinterpret_cast<void **>(&role_ptr), expected, 0);

        if (SWIG_IsOK(res) && role_ptr != nullptr)
        {
            temp_multid = new gstlrn::MultID(*role_ptr);
        }
        else
        {
            Rf_error("Cannot extract C++ pointer for ERole");
        }
    }

    /* ================================================================== */
    /* 6. S4 MultID Object                                              */
    /* ================================================================== */

    else if (Rf_inherits(obj, "_p_gstlrn__MultID"))
    {
        static SEXP sym_ref = Rf_install("ref");
        static swig_type_info *type_MultID = SWIG_TypeQuery("gstlrn::MultID *");

        SEXP ref = R_do_slot(obj, sym_ref);
        gstlrn::MultID *mult_ptr = nullptr;

        swig_type_info *expected = type_MultID ? type_MultID : SWIGTYPE_p_gstlrn__MultID;
        int res = SWIG_ConvertPtr(ref, reinterpret_cast<void **>(&mult_ptr), expected, 0);

        if (SWIG_IsOK(res) && mult_ptr != nullptr)
        {
            temp_multid = new gstlrn::MultID(*mult_ptr);
        }
        else
        {
            Rf_error("Cannot extract C++ pointer for MultID");
        }
    }

    /* ================================================================== */
    /* 7. Failure                                                         */
    /* ================================================================== */

    else
    {
        Rf_error("Cannot convert R object to MultID");
    }

    $1 = temp_multid;
}

%typemap(freearg) gstlrn::MultID&&, const gstlrn::MultID&
{
    if ($1)
        delete $1;
}
