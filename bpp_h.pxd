from  libcpp.string  cimport string as libcpp_string
from  libcpp.vector  cimport vector as libcpp_vector
from  libcpp.pair    cimport pair   as libcpp_pair
from  libcpp cimport bool
cdef extern from "src/Alignment.h":
    cdef cppclass Alignment:
        Alignment() except +
        Alignment(libcpp_vector[Alignment] alignments) except +
        Alignment(libcpp_vector[libcpp_pair[libcpp_string, libcpp_string]], libcpp_string datatype) except +
        Alignment(libcpp_string filename, libcpp_string file_format, bool interleaved) except +
        Alignment(libcpp_string filename, libcpp_string file_format, libcpp_string datatype, bool interleaved) except +
        Alignment(libcpp_string filename, libcpp_string file_format, libcpp_string datatype, libcpp_string model_name, bool interleaved) except +
        void read_alignment(libcpp_string filename, libcpp_string file_format, bool interleaved) except +
        void read_alignment(libcpp_string filename, libcpp_string file_format, libcpp_string datatype, bool interleaved) except +
        void sort_alignment(bool ascending) except +
        void write_alignment(libcpp_string filename, libcpp_string file_format, bool interleaved) except +
        void set_substitution_model(libcpp_string model_name) except +
        void set_gamma_rate_model(size_t ncat, double alpha) except +
        void set_constant_rate_model() except +
        void set_alpha(double alpha) except +
        void set_number_of_gamma_categories(size_t ncat) except +
        void set_rates(libcpp_vector[double], libcpp_string order) except +
        void set_frequencies(libcpp_vector[double]) except +
        void set_namespace(libcpp_string name) except +
        libcpp_vector[libcpp_pair[libcpp_string, libcpp_string]] get_sequences() except +
        double get_alpha() except +
        size_t get_number_of_gamma_categories() except +
        libcpp_vector[double] get_rates(libcpp_string order) except +
        libcpp_vector[double] get_rate_model_categories() except +
        libcpp_vector[double] get_frequencies() except +
        libcpp_vector[libcpp_string] get_names() except +
        libcpp_vector[libcpp_vector[double]] get_exchangeabilities() except +
        libcpp_string get_namespace() except +
        libcpp_string get_substitution_model() except +
        size_t get_number_of_sequences() except +
        libcpp_vector[libcpp_string] get_sites() except +
        size_t get_number_of_sites() except +
        libcpp_vector[libcpp_string] get_informative_sites(bool exclude_gaps) except +
        size_t get_number_of_informative_sites(bool exclude_gaps) except +
        size_t get_number_of_distinct_sites() except +
        size_t get_number_of_free_parameters() except +
        bool is_dna() except +
        bool is_protein() except +
        void _print_params() except +

        # Distance
        void compute_distances() except +
        void fast_compute_distances() except +
        void set_distance_matrix(libcpp_vector[libcpp_vector[double]] matrix) except +
        libcpp_string get_bionj_tree() except +
        libcpp_string get_bionj_tree(libcpp_vector[libcpp_vector[double]] matrix) except +
        libcpp_vector[libcpp_vector[double]] get_distances() except +
        libcpp_vector[libcpp_vector[double]] get_variances() except +
        libcpp_vector[libcpp_vector[double]] get_distance_variance_matrix() except +

        # Likelihood
        void initialise_likelihood() except +
        void initialise_likelihood(libcpp_string tree) except +
        void optimise_parameters(bool fix_branch_lengths) except +
        void optimise_topology(bool fix_model_params) except +
        double get_likelihood() except +
        libcpp_string get_tree() except +

        # Simulator
        void write_simulation(size_t nsites, libcpp_string filename, libcpp_string file_format, bool interleaved) except +
        void set_simulator(libcpp_string tree) except +
        libcpp_vector[libcpp_pair[libcpp_string, libcpp_string]] simulate(size_t nsites, libcpp_string tree) except +
        libcpp_vector[libcpp_pair[libcpp_string, libcpp_string]] simulate(size_t nsites) except +
        libcpp_vector[libcpp_pair[libcpp_string, libcpp_string]] get_simulated_sequences() except +

        # Bootstrap
        libcpp_vector[libcpp_pair[libcpp_string, libcpp_string]] get_bootstrapped_sequences() except +

        # Misc
        libcpp_string get_mrp_supertree(libcpp_vector[libcpp_string]) except +

        # Test
        void chkdst() except +