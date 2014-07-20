#cython: c_string_encoding=ascii  # for cython>=0.19
from  libcpp.string  cimport string as libcpp_string
from  libcpp.set     cimport set as libcpp_set
from  libcpp.vector  cimport vector as libcpp_vector
from  libcpp.pair    cimport pair as libcpp_pair
from  libcpp.map     cimport map  as libcpp_map
from  smart_ptr cimport shared_ptr
from  AutowrapRefHolder cimport AutowrapRefHolder
from  libcpp cimport bool
from  libc.string cimport const_char
from cython.operator cimport dereference as deref, preincrement as inc, address as address
from Distance cimport Distance as _Distance
cdef extern from "autowrap_tools.hpp":
    char * _cast_const_away(char *) 

cdef class Distance:

    cdef shared_ptr[_Distance] inst

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self, bytes filename , bytes file_format , bytes datatype , bytes model_name ,  interleaved ):
        assert isinstance(filename, bytes), 'arg filename wrong type'
        assert isinstance(file_format, bytes), 'arg file_format wrong type'
        assert isinstance(datatype, bytes), 'arg datatype wrong type'
        assert isinstance(model_name, bytes), 'arg model_name wrong type'
        assert isinstance(interleaved, (int, long)), 'arg interleaved wrong type'
    
    
    
    
    
        self.inst = shared_ptr[_Distance](new _Distance((<libcpp_string>filename), (<libcpp_string>file_format), (<libcpp_string>datatype), (<libcpp_string>model_name), (<bool>interleaved)))
    
    def set_rates(self, list in_0 , bytes order ):
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, float) for elemt_rec in in_0), 'arg in_0 wrong type'
        assert isinstance(order, bytes), 'arg order wrong type'
        cdef libcpp_vector[double] v0 = in_0
    
        self.inst.get().set_rates(v0, (<libcpp_string>order))
        
    
    def is_dna(self):
        cdef bool _r = self.inst.get().is_dna()
        py_result = <bool>_r
        return py_result
    
    def set_frequencies(self, list in_0 ):
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, float) for elemt_rec in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[double] v0 = in_0
        self.inst.get().set_frequencies(v0)
        
    
    def is_protein(self):
        cdef bool _r = self.inst.get().is_protein()
        py_result = <bool>_r
        return py_result
    
    def compute_distances(self):
        self.inst.get().compute_distances()
    
    def set_model(self, bytes model_name ):
        assert isinstance(model_name, bytes), 'arg model_name wrong type'
    
        self.inst.get().set_model((<libcpp_string>model_name))
    
    def read_alignment(self, bytes filename , bytes file_format , bytes datatype ,  interleaved ):
        assert isinstance(filename, bytes), 'arg filename wrong type'
        assert isinstance(file_format, bytes), 'arg file_format wrong type'
        assert isinstance(datatype, bytes), 'arg datatype wrong type'
        assert isinstance(interleaved, (int, long)), 'arg interleaved wrong type'
    
    
    
    
        self.inst.get().read_alignment((<libcpp_string>filename), (<libcpp_string>file_format), (<libcpp_string>datatype), (<bool>interleaved))
    
    def set_alpha(self,  ncat , double alpha ):
        assert isinstance(ncat, (int, long)), 'arg ncat wrong type'
        assert isinstance(alpha, float), 'arg alpha wrong type'
    
    
        self.inst.get().set_alpha((<int>ncat), (<double>alpha))
    
    def get_distances(self):
        _r = self.inst.get().get_distances()
        cdef list py_result = _r
        return py_result
    
    def get_names(self):
        _r = self.inst.get().get_names()
        cdef list py_result = _r
        return py_result 
