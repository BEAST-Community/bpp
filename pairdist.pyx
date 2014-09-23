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
from pairdist_h cimport Distance as _Distance
from pairdist_h cimport Simulator as _Simulator
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
    
    def get_tree(self):
        cdef libcpp_string _r = self.inst.get().get_tree()
        py_result = <libcpp_string>_r
        return py_result
    
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
    
    def get_likelihood(self):
        cdef double _r = self.inst.get().get_likelihood()
        py_result = <double>_r
        return py_result
    
    def initialise_likelihood(self, bytes tree ):
        assert isinstance(tree, bytes), 'arg tree wrong type'
    
        self.inst.get().initialise_likelihood((<libcpp_string>tree))
    
    def get_names(self):
        _r = self.inst.get().get_names()
        cdef list py_result = _r
        return py_result
    
    def optimise_parameters(self,  fix_branch_lengths ):
        assert isinstance(fix_branch_lengths, (int, long)), 'arg fix_branch_lengths wrong type'
    
        self.inst.get().optimise_parameters((<bool>fix_branch_lengths)) 

cdef class Simulator:

    cdef shared_ptr[_Simulator] inst

    def __dealloc__(self):
         self.inst.reset()

    
    def set_rates(self, list in_0 , bytes order ):
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, float) for elemt_rec in in_0), 'arg in_0 wrong type'
        assert isinstance(order, bytes), 'arg order wrong type'
        cdef libcpp_vector[double] v0 = in_0
    
        self.inst.get().set_rates(v0, (<libcpp_string>order))
        
    
    def __init__(self, bytes model_name ):
        assert isinstance(model_name, bytes), 'arg model_name wrong type'
    
        self.inst = shared_ptr[_Simulator](new _Simulator((<libcpp_string>model_name)))
    
    def set_alpha(self,  ncat , double alpha ):
        assert isinstance(ncat, (int, long)), 'arg ncat wrong type'
        assert isinstance(alpha, float), 'arg alpha wrong type'
    
    
        self.inst.get().set_alpha((<int>ncat), (<double>alpha))
    
    def get_frequencies(self):
        _r = self.inst.get().get_frequencies()
        cdef list py_result = _r
        return py_result
    
    def is_protein(self):
        cdef bool _r = self.inst.get().is_protein()
        py_result = <bool>_r
        return py_result
    
    def set_model(self, bytes model_name ):
        assert isinstance(model_name, bytes), 'arg model_name wrong type'
    
        self.inst.get().set_model((<libcpp_string>model_name))
    
    def is_dna(self):
        cdef bool _r = self.inst.get().is_dna()
        py_result = <bool>_r
        return py_result
    
    def set_frequencies(self, list in_0 ):
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, float) for elemt_rec in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[double] v0 = in_0
        self.inst.get().set_frequencies(v0)
        
    
    def _simulate_0(self,  nsites , bytes tree ):
        assert isinstance(nsites, (int, long)), 'arg nsites wrong type'
        assert isinstance(tree, bytes), 'arg tree wrong type'
    
    
        _r = self.inst.get().simulate((<unsigned int>nsites), (<libcpp_string>tree))
        cdef list py_result = _r
        return py_result
    
    def _simulate_1(self,  nsites ):
        assert isinstance(nsites, (int, long)), 'arg nsites wrong type'
    
        _r = self.inst.get().simulate((<unsigned int>nsites))
        cdef list py_result = _r
        return py_result
    
    def simulate(self, *args):
        if (len(args)==2) and (isinstance(args[0], (int, long))) and (isinstance(args[1], bytes)):
            return self._simulate_0(*args)
        elif (len(args)==1) and (isinstance(args[0], (int, long))):
            return self._simulate_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def write_alignment(self, bytes filename , bytes file_format ,  interleaved ):
        assert isinstance(filename, bytes), 'arg filename wrong type'
        assert isinstance(file_format, bytes), 'arg file_format wrong type'
        assert isinstance(interleaved, (int, long)), 'arg interleaved wrong type'
    
    
    
        self.inst.get().write_alignment((<libcpp_string>filename), (<libcpp_string>file_format), (<bool>interleaved))
    
    def set_simulator(self, bytes tree ):
        assert isinstance(tree, bytes), 'arg tree wrong type'
    
        self.inst.get().set_simulator((<libcpp_string>tree))
    
    def get_rates(self, bytes order ):
        assert isinstance(order, bytes), 'arg order wrong type'
    
        _r = self.inst.get().get_rates((<libcpp_string>order))
        cdef list py_result = _r
        return py_result
    
    def get_alpha(self):
        cdef double _r = self.inst.get().get_alpha()
        py_result = <double>_r
        return py_result
    
    def get_sequences(self):
        _r = self.inst.get().get_sequences()
        cdef list py_result = _r
        return py_result 
