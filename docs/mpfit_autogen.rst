mpfit overview
==============

``mpfit`` is a wrapper around the `cMPFIT
library <https://www.physics.wisc.edu/~craigm/idl/cmpfit.html>`__ for
Nim.

Usage of the library is centered around a single ``fit`` procedure, with
signature:

.. code:: nim

   proc fit*[T](f: FuncProto[T], pS: openArray[T], x, y, ey: openArray[T]): (seq[T], mp_result) =

the first argument is a user defined function (see below), the following
arguments are:

-  ``pS``: the first guess for the parameters
-  ``x``: data for x
-  ``y``: data for y
-  ``ey``: errors for y

The ``FuncProto[T]`` type is simply a specific signature for a proc,
which takes a sequence of parameters and a value. It is the user defined
function to be fitted.

Internally the function and the data is stored in a ``VarStruct``
object:

.. code:: nim

   type
     VarStruct*[T] = ref object
       x: seq[T]
       y: seq[T]
       ey: seq[T]
       f: FuncProto[T]

This is done, because the C library accepts a function to be fitted of
the form:

.. code:: nim

   mp_func* = proc (m: cint; n: cint; x: ptr cdouble; fvec: ptr cdouble;
                    dvec: ptr ptr cdouble; private_data: var pointer): cint {.cdecl.}

where ``m`` is the number of data points, ``n`` the number of
parameters, ``x`` a pointer to the sequence of parameters, ``fvec`` a
pointer to the sequence of deviates and ``dvec`` allows for user
computed derivates (which are currently not implemented).

The last parameter ``private_data`` is an opaque pointer and the reason
for the ``VarStruct`` object. We hand the user data and user function to
the ``mp_func`` via a cast to a pointer. A ``funcImpl`` wrapper proc of
the signature ``mp_func`` casts this opaque pointer back to
``VarStruct[float]`` (yes, currently it explicitly casts back to float,
so the allowed signature is actually only ``FuncProto[float]`` for the
time being), and calls the user function with the data computing the
deviates:

.. code:: nim

   for i in 0 ..< m:
     f = ff(pCall, x[i])
     dy[i] = (y[i] - f) / ey[i]

where ``ff`` is the user function and ``pCall`` the parameters as a
``seq[float]`` (since they have to be handed as a ``ptr cdouble`` we
have to convert to ``seq[float]`` first).

For an example see the `file:../README.org <../README.org>`__.
