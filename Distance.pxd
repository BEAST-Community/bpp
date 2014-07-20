from  libcpp.string  cimport string as libcpp_string
from  libcpp.vector  cimport vector as libcpp_vector
from  libcpp cimport bool
cdef extern from "src/Distance.h":
    cdef cppclass Distance:
        Distance(libcpp_string filename, libcpp_string file_format, libcpp_string datatype, libcpp_string model_name, bool interleaved) except +
        void read_alignment(libcpp_string filename, libcpp_string file_format, libcpp_string datatype, bool interleaved) except +
        void set_model(libcpp_string model_name) except +
        void set_alpha(int ncat, double alpha) except +
        void set_rates(libcpp_vector[double], libcpp_string order) except +
        void set_frequencies(libcpp_vector[double]) except +
        libcpp_vector[libcpp_vector[double]] get_distances() except +
        libcpp_vector[libcpp_string] get_names() except +
        void compute_distances() except +
        bool is_dna() except +
        bool is_protein() except +
