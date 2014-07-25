from  libcpp.string  cimport string as libcpp_string
from  libcpp.vector  cimport vector as libcpp_vector
from  libcpp.pair    cimport pair   as libcpp_pair
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


cdef extern from "src/Simulator.h":
    cdef cppclass Simulator:
        Simulator(libcpp_string model_name) except +
        void set_model(libcpp_string model_name) except +
        void set_alpha(int ncat, double alpha) except +
        void set_rates(libcpp_vector[double], libcpp_string order) except +
        void set_frequencies(libcpp_vector[double]) except +
        double get_alpha() except +
        libcpp_vector[double] get_rates(libcpp_string order) except +
        libcpp_vector[double] get_frequencies() except +
        bool is_dna() except +
        bool is_protein() except +
        void set_simulator(libcpp_string tree) except +
        libcpp_vector[libcpp_pair[libcpp_string, libcpp_string]] simulate(unsigned int nsites, libcpp_string tree) except +
        libcpp_vector[libcpp_pair[libcpp_string, libcpp_string]] simulate(unsigned int nsites) except +
        void write_alignment(libcpp_string filename, libcpp_string file_format, bool interleaved) except +
        libcpp_vector[libcpp_pair[libcpp_string, libcpp_string]] get_sequences() except +
