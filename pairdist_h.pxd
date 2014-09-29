from  libcpp.string  cimport string as libcpp_string
from  libcpp.vector  cimport vector as libcpp_vector
from  libcpp.pair    cimport pair   as libcpp_pair
from  libcpp cimport bool
cdef extern from "src/Alignment.h":
    cdef cppclass Alignment:
        Alignment(libcpp_vector[libcpp_pair[libcpp_string, libcpp_string]], libcpp_string datatype) except +
        Alignment(libcpp_string filename, libcpp_string file_format, libcpp_string datatype, bool interleaved) except +
        Alignment(libcpp_string filename, libcpp_string file_format, libcpp_string datatype, libcpp_string model_name, bool interleaved) except +
        void read_alignment(libcpp_string filename, libcpp_string file_format, libcpp_string datatype, bool interleaved) except +
        void set_model(libcpp_string model_name) except +
        void set_alpha(double alpha) except +
        void set_gamma(int ncat, double alpha) except +
        void set_rates(libcpp_vector[double], libcpp_string order) except +
        void set_frequencies(libcpp_vector[double]) except +
        double get_alpha() except +
        libcpp_vector[double] get_rates(libcpp_string order) except +
        libcpp_vector[double] get_frequencies() except +
        libcpp_vector[libcpp_string] get_names() except +
        libcpp_vector[libcpp_vector[double]] get_exchangeabilities() except +
        libcpp_string get_model() except +
        bool is_dna() except +
        bool is_protein() except +

        # Distance
        void compute_distances() except +
        void fast_compute_distances() except +
        void set_distance_matrix(libcpp_vector[libcpp_vector[double]] matrix) except +
        libcpp_string get_nj_tree() except +
        libcpp_string get_nj_tree(libcpp_vector[libcpp_vector[double]] matrix) except +
        libcpp_vector[libcpp_vector[double]] get_distances() except +
        libcpp_vector[libcpp_vector[double]] get_variances() except +
        libcpp_vector[libcpp_vector[double]] get_distance_variance_matrix() except +

        # Likelihood
        void initialise_likelihood(libcpp_string tree) except +
        void optimise_parameters(bool fix_branch_lengths) except +
        double get_likelihood() except +
        libcpp_string get_tree() except +

        # Simulator
        void write_alignment(libcpp_string filename, libcpp_string file_format, bool interleaved) except +
        void set_simulator(libcpp_string tree) except +
        libcpp_vector[libcpp_pair[libcpp_string, libcpp_string]] simulate(unsigned int nsites, libcpp_string tree) except +
        libcpp_vector[libcpp_pair[libcpp_string, libcpp_string]] simulate(unsigned int nsites) except +
        libcpp_vector[libcpp_pair[libcpp_string, libcpp_string]] get_simulated_sequences() except +
