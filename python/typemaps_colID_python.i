/***************************************************************************/
/*                                                                         */
/*  Python typemaps for gstlrn::ColID                                      */
/*                                                                         */
/*  Dual-binding support for both:                                         */
/*    1. const gstlrn::ColID&  (Preferred API facade in Db.hpp)            */
/*    2. gstlrn::ColID&&       (Legacy / Low-level engine in DbData.hpp)   */
/*                                                                         */
/*  This typemap transparently converts high-level Python arguments into   */
/*  a C++ ColID temporary instance.                                        */
/*                                                                         */
/*  Accepted Python inputs:                                                */
/*                                                                         */
/*    - String:              "name"                -> ColID("name")        */
/*    - String Tuple:        ("name", version)     -> ColID("name", v)     */
/*    - Integer:             icol                  -> ColID(icol)          */
/*    - Integer Tuple:       (icol, version)       -> ColID(icol, v)       */
/*    - RoleID object:       RoleID(role, index)   -> ColID(RoleID)        */
/*    - RoleID Tuple:        (RoleID, version)     -> ColID(RoleID, v)     */
/*    - ERole enum:          ERole.X               -> ColID(ERole::X)      */
/*    - ColID object:        ColID(...)            -> Copy constructor     */
/*    - ColID Tuple:         (ColID, version)      -> ColID(ColID, v)      */
/*                                                                         */
/*  Memory Safety:                                                         */
/*    A heap-allocated ColID object is constructed in %typemap(in) and     */
/*    automatically reclaimed in %typemap(freearg) after the C++ method    */
/*    call returns.                                                        */
/*                                                                         */
/***************************************************************************/

%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER)
    const gstlrn::ColID&,
    gstlrn::ColID&&
{
    PyObject *obj = $input;
    if (PyTuple_Check(obj) && PyTuple_GET_SIZE(obj) == 2)
    {
        obj = PyTuple_GET_ITEM(obj, 0);
    }

    if (PyUnicode_Check(obj) || PyLong_Check(obj))
    {
        $1 = 1;
    }
    else
    {
        static swig_type_info *type_ColID  = SWIG_TypeQuery("gstlrn::ColID *");
        static swig_type_info *type_RoleID = SWIG_TypeQuery("gstlrn::RoleID *");
        static swig_type_info *type_ERole  = SWIG_TypeQuery("gstlrn::ERole *");

        void *ptr = nullptr;
        if (SWIG_IsOK(SWIG_ConvertPtr(obj, &ptr, type_ColID ? type_ColID : SWIGTYPE_p_gstlrn__ColID, 0)) && ptr)
        {
            $1 = 1;
        }
        else if (SWIG_IsOK(SWIG_ConvertPtr(obj, &ptr, type_ERole ? type_ERole : SWIGTYPE_p_gstlrn__ERole, 0)) && ptr)
        {
            $1 = 1;
        }
        else if (SWIG_IsOK(SWIG_ConvertPtr(obj, &ptr, type_RoleID ? type_RoleID : SWIGTYPE_p_gstlrn__RoleID, 0)) && ptr)
        {
            $1 = 1;
        }
        else
        {
            $1 = 0;
        }
    }
}


%typemap(in)
    const gstlrn::ColID&,
    gstlrn::ColID&&
{
    PyObject *obj = $input;
    gstlrn::Id version = 0;

    /**********************************************************************/
    /* Optional tuple: (object, version)                                  */
    /**********************************************************************/

    if (PyTuple_Check(obj))
    {
        if (PyTuple_GET_SIZE(obj) != 2)
        {
            SWIG_exception_fail(
                SWIG_TypeError,
                "ColID tuple must contain exactly two elements");
        }

        PyObject *first  = PyTuple_GET_ITEM(obj, 0);
        PyObject *second = PyTuple_GET_ITEM(obj, 1);

        if (!PyLong_Check(second))
        {
            SWIG_exception_fail(
                SWIG_TypeError,
                "Second tuple element must be an integer version");
        }

        version = static_cast<gstlrn::Id>(PyLong_AsLong(second));
        obj = first;
    }


    /**********************************************************************/
    /* Fast path: Native Python types (String & Integer)                  */
    /**********************************************************************/

    if (PyUnicode_Check(obj))
    {
        Py_ssize_t size = 0;
        const char *str = PyUnicode_AsUTF8AndSize(obj, &size);

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

    else if (PyLong_Check(obj))
    {
        $1 = new gstlrn::ColID(
            gstlrn::ColID::create(
                static_cast<gstlrn::Id>(PyLong_AsLong(obj)),
                version));
    }


    /**********************************************************************/
    /* C++ Wrapped Objects (ColID, RoleID, ERole) with Type Query Cache   */
    /**********************************************************************/

    else
    {
        static swig_type_info *type_ColID  = SWIG_TypeQuery("gstlrn::ColID *");
        static swig_type_info *type_RoleID = SWIG_TypeQuery("gstlrn::RoleID *");
        static swig_type_info *type_ERole  = SWIG_TypeQuery("gstlrn::ERole *");

        gstlrn::ColID *colid = nullptr;

        int res = SWIG_ConvertPtr(
            obj,
            (void**)&colid,
            type_ColID ? type_ColID : SWIGTYPE_p_gstlrn__ColID,
            0);

        if (SWIG_IsOK(res) && colid != nullptr)
        {
            $1 = new gstlrn::ColID(
                gstlrn::ColID::create(*colid, version));
        }
        else
        {
            // Note: ERole checked BEFORE RoleID to avoid SWIG buffer overread / type confusion
            gstlrn::ERole *role = nullptr;

            res = SWIG_ConvertPtr(
                obj,
                (void**)&role,
                type_ERole ? type_ERole : SWIGTYPE_p_gstlrn__ERole,
                0);

            if (SWIG_IsOK(res) && role != nullptr)
            {
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
                gstlrn::RoleID *roleID = nullptr;

                res = SWIG_ConvertPtr(
                    obj,
                    (void**)&roleID,
                    type_RoleID ? type_RoleID : SWIGTYPE_p_gstlrn__RoleID,
                    0);

                if (SWIG_IsOK(res) && roleID != nullptr)
                {
                    $1 = new gstlrn::ColID(
                        gstlrn::ColID::create(*roleID, version));
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
}


%typemap(freearg)
    const gstlrn::ColID&,
    gstlrn::ColID&&
{
    delete $1;
}
