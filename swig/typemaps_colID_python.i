/***************************************************************************/
/*                                                                         */
/*  Python typemap for gstlrn::ColID&&                                     */
/*                                                                         */
/*  Accepted Python objects:                                               */
/*                                                                         */
/*      "name"                  -> ColID("name")                           */
/*      ("name",version)        -> ColID("name",version)                   */
/*                                                                         */
/*      icol                    -> ColID(icol)                             */
/*      (icol,version)          -> ColID(icol,version)                     */
/*                                                                         */
/*      RoleID(role,index)      -> ColID(RoleID(...))                      */
/*      (RoleID(...),version)   -> ColID(RoleID(...),version)              */
/*                                                                         */
/*      ERole.X                 -> ColID(ERole::X)                         */
/*                                                                         */
/*      ColID(...)              -> copy constructor                        */
/*                                                                         */
/***************************************************************************/
%typemap(in) gstlrn::ColID&&
{
    PyObject *obj = $input;
    gstlrn::Id version = 0;

    /**********************************************************************/
    /* Optional tuple: (object,version)                                   */
    /**********************************************************************/

    if (PyTuple_Check(obj))
    {
        if (PyTuple_Size(obj) != 2)
        {
            SWIG_exception_fail(
                SWIG_TypeError,
                "ColID tuple must contain exactly two elements");
        }

        PyObject *first  = PyTuple_GET_ITEM(obj,0);
        PyObject *second = PyTuple_GET_ITEM(obj,1);

        if (!PyLong_Check(second))
        {
            SWIG_exception_fail(
                SWIG_TypeError,
                "Second tuple element must be an integer version");
        }

        version =
            static_cast<gstlrn::Id>(PyLong_AsLong(second));

        obj = first;
    }


    /**********************************************************************/
    /* Existing ColID                                                     */
    /**********************************************************************/

    gstlrn::ColID *colid = nullptr;

    int res = SWIG_ConvertPtr(
        obj,
        (void**)&colid,
        SWIGTYPE_p_gstlrn__ColID,
        0);

    if (SWIG_IsOK(res) && colid != nullptr)
    {
        $1 = new gstlrn::ColID(
            gstlrn::ColID::create(*colid, version));
    }


    /**********************************************************************/
    /* String                                                             */
    /**********************************************************************/

    else if (PyUnicode_Check(obj))
    {
        Py_ssize_t size = 0;

        const char *str =
            PyUnicode_AsUTF8AndSize(obj, &size);

        if (str == nullptr)
        {
            SWIG_exception_fail(
                SWIG_RuntimeError,
                "Cannot convert Python string");
        }

        $1 = new gstlrn::ColID(
            gstlrn::ColID::create(
                std::string_view(str, size),
                version));
    }


    /**********************************************************************/
    /* Integer                                                            */
    /**********************************************************************/

    else if (PyLong_Check(obj))
    {
        $1 = new gstlrn::ColID(
            gstlrn::ColID::create(
                static_cast<gstlrn::Id>(PyLong_AsLong(obj)),
                version));
    }


    /**********************************************************************/
    /* RoleID                                                             */
    /**********************************************************************/

    else
    {
        gstlrn::RoleID *roleID = nullptr;

        res = SWIG_ConvertPtr(
            obj,
            (void**)&roleID,
            SWIGTYPE_p_gstlrn__RoleID,
            0);

        if (SWIG_IsOK(res) && roleID != nullptr)
        {
            $1 = new gstlrn::ColID(
                gstlrn::ColID::create(*roleID, version));
        }


        /******************************************************************/
        /* ERole                                                           */
        /******************************************************************/

        else
        {
            gstlrn::ERole *role = nullptr;

            res = SWIG_ConvertPtr(
                obj,
                (void**)&role,
                SWIGTYPE_p_gstlrn__ERole,
                0);

            if (SWIG_IsOK(res) && role != nullptr)
            {
                /*
                 * Tuple (ERole,version) deliberately not supported:
                 * ambiguous with (ERole,index)
                 */
                if (version != 0)
                {
                    SWIG_exception_fail(
                        SWIG_TypeError,
                        "(ERole,version) syntax is not supported");
                }

                $1 = new gstlrn::ColID(
                    gstlrn::ColID::create(*role));
            }

            else
            {
                SWIG_exception_fail(
                    SWIG_TypeError,
                    "Expected ColID, string, integer, RoleID, ERole, "
                    "(string,version), (integer,version) or "
                    "(RoleID,version)");
            }
        }
    }
}


%typemap(freearg) gstlrn::ColID&&
{
    delete $1;
}
