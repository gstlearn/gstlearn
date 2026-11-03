/***************************************************************************/
/*                                                                         */
/*  Python typemap for gstlrn::MultID&& and const gstlrn::MultID&          */
/*                                                                         */
/*  Accepted Python objects:                                               */
/*                                                                         */
/*      "name"                  -> MultID("name")                          */
/*      ["name1", "name2"]      -> MultID(VectorString)                    */
/*      ("name1", "name2"]      -> MultID(VectorString)                    */
/*      np.array([...])         -> MultID(VectorString or VectorInt)       */
/*                                                                         */
/*      icol                    -> MultID(Id icol)                         */
/*      [icol1, icol2]          -> MultID(VectorInt)                       */
/*      (icol1, icol2)          -> MultID(VectorInt)                       */
/*                                                                         */
/*      ERole.X                 -> MultID(ERole::X)                        */
/*      MultID(...)             -> copy constructor                        */
/*                                                                         */
/***************************************************************************/
%typemap(in) gstlrn::MultID&&, const gstlrn::MultID&
{
    PyObject *obj = $input;

    /**********************************************************************/
    /* Helper lambda/logic to extract a vector of strings from sequence   */
    /**********************************************************************/
    auto extract_vector_string = [](PyObject *seq, gstlrn::VectorString &out_vec) -> bool {
        Py_ssize_t size = PySequence_Size(seq);
        for (Py_ssize_t i = 0; i < size; ++i)
        {
            PyObject *item = PySequence_GetItem(seq, i);
            if (!PyUnicode_Check(item))
            {
                Py_DECREF(item);
                return false;
            }
            Py_ssize_t len = 0;
            const char *str = PyUnicode_AsUTF8AndSize(item, &len);
            if (str == nullptr)
            {
                PyErr_Clear();
                Py_DECREF(item);
                return false;
            }
            out_vec.push_back(std::string(str, len));
            Py_DECREF(item);
        }
        return true;
    };

    /**********************************************************************/
    /* Helper lambda/logic to extract a vector of integers from sequence  */
    /**********************************************************************/
    auto extract_vector_int = [](PyObject *seq, gstlrn::VectorInt &out_vec) -> bool {
        Py_ssize_t size = PySequence_Size(seq);
        for (Py_ssize_t i = 0; i < size; ++i)
        {
            PyObject *item = PySequence_GetItem(seq, i);
            if (!PyLong_Check(item) && !PyNumber_Check(item))
            {
                Py_DECREF(item);
                return false;
            }
            long val = PyLong_AsLong(item);
            if (PyErr_Occurred())
            {
                PyErr_Clear();
                Py_DECREF(item);
                return false;
            }
            out_vec.push_back(static_cast<gstlrn::Id>(val));
            Py_DECREF(item);
        }
        return true;
    };

    /**********************************************************************/
    /* 1. Fast path: Native Python String -> MultID(const String&)        */
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

        $1 = new gstlrn::MultID(std::string(str, size));
    }
    /**********************************************************************/
    /* 2. Fast path: Native Python Integer -> MultID(Id icol)             */
    /**********************************************************************/
    else if (PyLong_Check(obj))
    {
        long val = PyLong_AsLong(obj);
        if (PyErr_Occurred())
        {
            PyErr_Clear();
            SWIG_exception_fail(SWIG_TypeError, "Cannot convert number to integer ID");
        }
        $1 = new gstlrn::MultID(static_cast<gstlrn::Id>(val));
    }
    /**********************************************************************/
    /* 3. General Sequences (List, Tuple, NumPy Array, etc.)              */
    /**********************************************************************/
    else if (PySequence_Check(obj) && !PyUnicode_Check(obj) && !PyBytes_Check(obj))
    {
        Py_ssize_t size = PySequence_Size(obj);
        if (size == 0)
        {
            SWIG_exception_fail(
                SWIG_TypeError,
                "Sequence cannot be empty for MultID conversion");
        }

        // Check the type of the first element to decide strictly between strings and numbers
        PyObject *first_item = PySequence_GetItem(obj, 0);
        bool is_string_seq = PyUnicode_Check(first_item);
        Py_DECREF(first_item);

        if (is_string_seq)
        {
            gstlrn::VectorString names;
            if (!extract_vector_string(obj, names))
            {
                SWIG_exception_fail(
                    SWIG_TypeError,
                    "All elements in sequence must be strings");
            }
            $1 = new gstlrn::MultID(names);
        }
        else
        {
            gstlrn::VectorInt icols;
            if (!extract_vector_int(obj, icols))
            {
                SWIG_exception_fail(
                    SWIG_TypeError,
                    "All elements in sequence must be integers");
            }
            $1 = new gstlrn::MultID(icols);
        }
    }
    /**********************************************************************/
    /* 4. C++ Wrapped Objects (MultID, ERole)                             */
    /**********************************************************************/
    else
    {
        static swig_type_info *type_MultID = SWIG_TypeQuery("gstlrn::MultID *");
        static swig_type_info *type_ERole  = SWIG_TypeQuery("gstlrn::ERole *");

        gstlrn::MultID *multid = nullptr;
        int res = SWIG_ConvertPtr(
            obj,
            (void**)&multid,
            type_MultID ? type_MultID : SWIGTYPE_p_gstlrn__MultID,
            0);

        if (SWIG_IsOK(res) && multid != nullptr)
        {
            $1 = new gstlrn::MultID(*multid);
        }
        else
        {
            gstlrn::ERole *role = nullptr;
            res = SWIG_ConvertPtr(
                obj,
                (void**)&role,
                type_ERole ? type_ERole : SWIGTYPE_p_gstlrn__ERole,
                0);

            if (SWIG_IsOK(res) && role != nullptr)
            {
                $1 = new gstlrn::MultID(*role);
            }
            else
            {
                SWIG_exception_fail(
                    SWIG_TypeError,
                    "Expected string, sequence of strings, number, "
                    "sequence of numbers, MultID, or ERole");
            }
        }
    }
}

%typemap(freearg) gstlrn::MultID&&, const gstlrn::MultID&
{
    delete $1;
}
