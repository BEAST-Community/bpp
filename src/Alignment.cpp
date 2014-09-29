/*
 * Alignment.cpp
 *
 *  Created on: Jul 20, 2014
 *      Author: kgori
 */

#include "Alignment.h"
#include "SiteContainerBuilder.h"
#include "ModelFactory.h"

#include "Bpp/Numeric/Prob/GammaDiscreteDistribution.h"
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include "Bpp/Seq/SymbolListTools.h"
#include "Bpp/Phyl/Distance/DistanceEstimation.h"
#include "Bpp/Phyl/Distance/BioNJ.h"
#include <Bpp/Phyl/Io/Newick.h>
#include "Bpp/Phyl/OptimizationTools.h"
#include <Bpp/Phyl/Simulation/HomogeneousSequenceSimulator.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Io/Phylip.h>

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>

using namespace bpp;
using namespace std;

double VARMIN = 0.000001;
double DISTMAX = 10000;

Alignment::Alignment(vector<pair<string, string>> headers_sequences, string datatype) {
    sequences = SiteContainerBuilder::construct_alignment_from_strings(headers_sequences, datatype);
    _set_datatype();
    set_gamma();
}

Alignment::Alignment(string filename, string file_format, string datatype, bool interleaved) {
    read_alignment(filename, file_format, datatype, interleaved);
    _set_datatype();
    set_gamma();
}

Alignment::Alignment(string filename, string file_format, string datatype, string model_name, bool interleaved) {
    read_alignment(filename, file_format, datatype, interleaved);
    set_gamma();
    try {
        _check_compatible_model(datatype, model_name);
        set_model(model_name);
    }
    catch (Exception& e) {
        cerr << "Alignment was initialised, but the model wasn't valid: " <<e.what() << endl;
    }
}

void Alignment::read_alignment(string filename, string file_format, string datatype, bool interleaved) {
    sequences = SiteContainerBuilder::read_alignment(filename, file_format, datatype, interleaved);
    _set_datatype();
    _clear_distances();
    _clear_likelihood();
}

void Alignment::set_model(string model_name) {
    string datatype = is_protein() ? "protein" : "dna";
    _check_compatible_model(datatype, model_name);
    model = ModelFactory::create(model_name);
    _model = model_name;
    _clear_distances();
    _clear_likelihood();
}

void Alignment::set_alpha(double alpha) {
    if(!rates) throw Exception("No rate model is set");
    rates->setParameterValue("alpha", alpha);
}

void Alignment::set_gamma(int ncat, double alpha) {
    rates = make_shared<GammaDiscreteDistribution>(ncat, alpha, alpha);
    rates->aliasParameters("alpha", "beta");
    _clear_distances();
    _clear_likelihood();
}

void Alignment::set_rates(vector<double> rates, string order) {
    if (is_dna()) {
        if (order == "acgt" || order == "ACGT") {
            double normaliser = rates[1];
            model->setParameterValue("a", rates[4] / normaliser);
            model->setParameterValue("b", rates[2] / normaliser);
            model->setParameterValue("c", rates[5] / normaliser);
            model->setParameterValue("d", rates[0] / normaliser);
            model->setParameterValue("e", rates[3] / normaliser);
        }
        else if (order == "tcag" || order == "TCAG") {
            model->setParameterValue("a", rates[0]);
            model->setParameterValue("b", rates[1]);
            model->setParameterValue("c", rates[2]);
            model->setParameterValue("d", rates[3]);
            model->setParameterValue("e", rates[4]);
        }
        else {
            throw Exception("Unrecognised order for rates: " + order);
        }
        _clear_distances();
        _clear_likelihood();
    }
}

void Alignment::set_frequencies(vector<double> freqs) {
    size_t reqd = is_dna() ? 4 : 20;
    if (freqs.size() != reqd) {
        throw Exception("Frequencies vector is the wrong length (dna: 4; aa: 20)");
    }
    map<int, double> m = _vector_to_map(freqs);
    model->setFreq(m);
    _clear_distances();
    _clear_likelihood();
}

double Alignment::get_alpha() {
    if (rates) {
        return rates->getParameterValue("alpha");
    }
    else throw Exception("Rate model not set");
}

vector<double> Alignment::get_rates(string order) {
    if (is_dna()) {
        vector<double> rates_vec;
        if (order == "acgt" || order == "ACGT") { //{a-c, a-g, a-t, c-g, c-t, g-t=1}
            double normaliser = model->getParameterValue("c");
            rates_vec.push_back(model->getParameterValue("d") / normaliser);
            rates_vec.push_back(1.0 / normaliser);
            rates_vec.push_back(model->getParameterValue("b") / normaliser);
            rates_vec.push_back(model->getParameterValue("e") / normaliser);
            rates_vec.push_back(model->getParameterValue("a") / normaliser);
            rates_vec.push_back(1.0);
        }
        else if (order == "tcag" || order == "TCAG") { //{a=t-c, b=t-a, c=t-g, d=c-a, e=c-g, f=a-g=1}
            rates_vec.push_back(model->getParameterValue("a"));
            rates_vec.push_back(model->getParameterValue("b"));
            rates_vec.push_back(model->getParameterValue("c"));
            rates_vec.push_back(model->getParameterValue("d"));
            rates_vec.push_back(model->getParameterValue("e"));
            rates_vec.push_back(1.0);
        }
        else {
            cerr << "Unknown order: " << order << ". Accepted orders are {tcag, acgt}" << endl;
            throw Exception("Unknown order error");
        }
        return rates_vec;
    }
    else {
        cerr << "Getting and setting rates is not implemented for protein models" << endl;
        throw Exception("Protein model disallowed error");
    }
}

vector<double> Alignment::get_frequencies() {
    return model->getFrequencies();
}

vector<string> Alignment::get_names() {
    return sequences->getSequencesNames();
}

string Alignment::get_model() {
    if (_model.empty()) throw Exception("No model name is set");
    return _model;
}

vector<vector<double>> Alignment::get_exchangeabilities() {
    if(!model) throw Exception("No model has been set.");
    RowMatrix<double> exch = model->getExchangeabilityMatrix();
    vector<vector<double>> matrix;
    for (size_t i = 0; i < exch.getNumberOfRows(); ++i) {
        vector<double> row;
        for (size_t j = 0; j < exch.getNumberOfColumns(); ++j) {
            row.push_back(exch(i, j));
        }
        matrix.push_back(row);
    }
    return matrix;
}

bool Alignment::is_dna() {
    return dna && !protein;
}

bool Alignment::is_protein() {
    return protein && !dna;
}

// Distance
void Alignment::compute_distances() {
    if (!model) throw Exception("No model of evolution available");
    VectorSiteContainer* sites_ = sequences->clone();
    SiteContainerTools::changeGapsToUnknownCharacters(*sites_);
    size_t n = sites_->getNumberOfSequences();
    vector<string> names = get_names();
    double var;

    _clear_distances();
    distances = make_shared<DistanceMatrix>(names);
    variances = make_shared<DistanceMatrix>(names);
    for (size_t i = 0; i < n; i++) {
        (*distances)(i, i) = 0;
        for (size_t j = i + 1; j < n; j++) {
            auto lik = make_shared<TwoTreeLikelihood>(names[i], names[j], *sites_, model.get(), rates.get(), false);
            lik->initialize();
            lik->enableDerivatives(true);
            size_t d = SymbolListTools::getNumberOfDistinctPositions(sites_->getSequence(i), sites_->getSequence(j));
            size_t g = SymbolListTools::getNumberOfPositionsWithoutGap(sites_->getSequence(i), sites_->getSequence(j));
            lik->setParameterValue("BrLen", g == 0 ? lik->getMinimumBranchLength() : std::max(lik->getMinimumBranchLength(), static_cast<double>(d) / static_cast<double>(g)));
            // Optimization:
            ParameterList params = lik->getBranchLengthsParameters();
            OptimizationTools::optimizeNumericalParameters(lik.get(), params, 0, 1, 0.000001, 1000000, NULL, NULL, false, 0, OptimizationTools::OPTIMIZATION_NEWTON, OptimizationTools::OPTIMIZATION_BRENT);
            // Store results:
            (*distances)(i, j) = (*distances)(j, i) = lik->getParameterValue("BrLen");
            var = 1.0 / lik->d2f("BrLen", params);
            (*variances)(i, j) = (*variances)(j, i) = var > VARMIN ? var : VARMIN;
        }
    }
    delete sites_;
}

void Alignment::fast_compute_distances() {
    unsigned int s;
    if (is_dna()) {
        s = 4;
    }
    if (is_protein()) {
        s = 20;
    }
    size_t n = sequences->getNumberOfSequences();
    vector<string> names = sequences->getSequencesNames();
    if (distances) distances.reset();
    distances = make_shared<DistanceMatrix>(names);
    variances = make_shared<DistanceMatrix>(names);
    for (size_t i = 0; i < n; i++) {
        (*distances)(i, i) = 0;
        for (size_t j=i+1; j < n; j++) {
            size_t d = SymbolListTools::getNumberOfDistinctPositions(sequences->getSequence(i), sequences->getSequence(j));
            size_t g = SymbolListTools::getNumberOfPositionsWithoutGap(sequences->getSequence(i), sequences->getSequence(j));
            double dist = _jcdist(d, g, s);
            double var = _jcvar(d, g, s);
            (*distances)(i, j) = (*distances)(j, i) = dist;
            (*variances)(i, j) = (*variances)(j, i) = var;
        }
    }
}

void Alignment::set_distance_matrix(vector<vector<double>> matrix) {
    distances = _create_distance_matrix(matrix);
}

string Alignment::get_nj_tree() {
    if (!distances) throw Exception("No distances have been calculated yet");
    BioNJ bionj(*distances, false, true, false); // rooted=false, positiveLengths=true, verbose=false
    auto njtree = bionj.getTree();
    stringstream ss;
    Newick treeWriter;
    treeWriter.write(*njtree, ss);
    delete njtree;
    string s{ss.str()};
    s.erase(s.find_last_not_of(" \n\r\t")+1);
    return s;
}

string Alignment::get_nj_tree(vector<vector<double>> matrix) {
    auto dm = _create_distance_matrix(matrix);
    BioNJ bionj(*dm, false, true, false); // rooted=false, positiveLengths=true, verbose=false
    auto njtree = bionj.getTree();
    stringstream ss;
    Newick treeWriter;
    treeWriter.write(*njtree, ss);
    delete njtree;
    string s{ss.str()};
    s.erase(s.find_last_not_of(" \n\r\t")+1);
    return s;
}

vector<vector<double>> Alignment::get_distances() {
    if(!distances) throw Exception("No distances have been calculated yet");
    vector<vector<double>> vec;
    size_t nrow = distances->getNumberOfRows();
    for (size_t i = 0; i < nrow; ++i) {
        vec.push_back(distances->row(i));
    }
    return vec;
}

vector<vector<double>> Alignment::get_variances() {
    if(!variances) throw Exception("No distances have been calculated yet");
    vector<vector<double>> vec;
    size_t nrow = variances->getNumberOfRows();
    for (size_t i = 0; i < nrow; ++i) {
        vec.push_back(variances->row(i));
    }
    return vec;
}

vector<vector<double>> Alignment::get_distance_variance_matrix() {
    if(!variances || !distances) throw Exception("No distances have been calculated yet");
    vector<vector<double>> vec;
    size_t nrow = variances->getNumberOfRows();
    for (size_t i = 0; i < nrow; ++i) {
        vector<double> row;
        for (size_t j = 0; j < nrow; ++j) {
            if (j < i) row.push_back((*variances)(i, j));
            else row.push_back((*distances)(i, j));
        }
        vec.push_back(row);
    }
    return vec;
}

// Likelihood
void Alignment::initialise_likelihood(string tree) {
    if (!model) {
        cerr << "Model not set" << endl;
        throw Exception("Model not set error");
    }
    if (!rates) {
        cerr << "Rates not set" << endl;
        throw Exception("Rates not set error");
    }
    Tree * liktree;
    auto reader = make_shared<Newick>(false);
    if (_is_file(tree)) {
        liktree = reader->read(tree);
    }
    else if (_is_tree_string(tree)) {
        stringstream ss{tree};
        liktree = reader->read(ss);
    }
    else {
        cerr << "Couldn\'t understand this tree: " << tree << endl;
        throw Exception("Tree error");
    }
    VectorSiteContainer* sites_ = sequences->clone();
    SiteContainerTools::changeGapsToUnknownCharacters(*sites_);
    likelihood = make_shared<RHomogeneousTreeLikelihood>(*liktree, *sites_, model.get(), rates.get(), true, false, true);
    likelihood->initialize();
    delete liktree;
}

void Alignment::optimise_parameters(bool fix_branch_lengths) {
    if (!likelihood) {
        cerr << "Likelihood calculator not set - call initialise_likelihood" << endl;
        throw Exception("Uninitialised likelihood error");
    }
    ParameterList pl;
    if (fix_branch_lengths) {
        pl = likelihood->getSubstitutionModelParameters();
        pl.addParameters(likelihood->getRateDistributionParameters());
    }
    else {
        pl = likelihood->getParameters();
    }
    OptimizationTools::optimizeNumericalParameters2(likelihood.get(), pl, 0, 0.0001, 1000000, NULL, NULL, false, false, 0);
}

double Alignment::get_likelihood() {
    if (!likelihood) {
        cerr << "Likelihood calculator not set - call initialise_likelihood" << endl;
        throw Exception("Uninitialised likelihood error");
    }
    return likelihood->getLogLikelihood();
}

string Alignment::get_tree() {
    if (!likelihood) {
        cerr << "Likelihood calculator not set - call initialise_likelihood" << endl;
        throw Exception("Uninitialised likelihood error");
    }
    auto *tree = likelihood->getTree().clone();
    stringstream ss;
    Newick treeWriter;
    treeWriter.write(*tree, ss);
    delete tree;
    string s{ss.str()};
    s.erase(s.find_last_not_of(" \n\r\t")+1);
    return s;
}

// Simulator
void Alignment::write_alignment(string filename, string file_format, bool interleaved) {
    if (file_format == "fas" || file_format == "fasta") {
        _write_fasta(filename);
    }
    else if (file_format == "phy" || file_format == "phylip") {
        _write_phylip(filename, interleaved);
    }
    else {
        cerr << "Unrecognised file format: " << file_format << endl;
        throw exception();
    }
}

void Alignment::set_simulator(string tree) {
    if (!model) {
        cerr << "Model not set" << endl;
        throw exception();
    }
    if (!rates) {
        cerr << "Rates not set" << endl;
        throw exception();
    }
    Tree * simtree;
    Newick * reader = new Newick(false);
    if (_is_file(tree)) {
        simtree = reader->read(tree);
        delete reader;
    }
    else if (_is_tree_string(tree)) {
        stringstream ss{tree};
        simtree = reader->read(ss);
        delete reader;
    }
    else {
        cerr << "Couldn\'t understand this tree: " << tree << endl;
        delete reader;
        throw exception();
    }
    simulator = make_shared<HomogeneousSequenceSimulator>(model.get(), rates.get(), simtree);
    delete simtree;
}

vector<pair<string, string>> Alignment::simulate(unsigned int nsites, string tree) {
    set_simulator(tree);
    return simulate(nsites);
}

vector<pair<string, string>> Alignment::simulate(unsigned int nsites) {
    if (!simulator) {
        cout << "Tried to simulate without a simulator" << endl;
        throw exception();
    }
    size_t nsites_{nsites};
    SiteContainer * tmp = simulator->simulate(nsites_);
    simulated_sequences = make_shared<VectorSiteContainer>(*tmp);
    delete tmp;
    return get_simulated_sequences();
}

vector<pair<string, string>> Alignment::get_simulated_sequences() {
    vector<pair<string, string>> ret;
    if (!simulated_sequences) {
        cerr << "No sequences exist yet" << endl;
        throw exception();
    }
    for (size_t i = 0; i < simulated_sequences->getNumberOfSequences(); ++i) {
        BasicSequence seq = simulated_sequences->getSequence(i);
        ret.push_back(make_pair(seq.getName(), seq.toString()));
    }
    return ret;
}

// Private methods
void Alignment::_set_dna() {
    dna = true;
    protein = false;
}

void Alignment::_set_protein() {
    dna = false;
    protein = true;
}

void Alignment::_set_datatype() {
    string type = sequences->getAlphabet()->getAlphabetType();
    if (type == "DNA alphabet") {
        _set_dna();
    }
    else if (type == "Proteic alphabet") {
        _set_protein();
    }
    else {
        cerr << "Type = " << type << endl;
        throw Exception("Unrecognised data type");
    }
}

void Alignment::_write_fasta(string filename) {
    if (!simulated_sequences) {
        cerr << "No sequences exist yet" << endl;
        throw exception();
    }
    Fasta writer;
    writer.writeAlignment(filename, *simulated_sequences);
}

void Alignment::_write_phylip(string filename, bool interleaved) {
    if (!simulated_sequences) {
        cerr << "No sequences exist yet" << endl;
        throw exception();
    }
    Phylip writer{true, !interleaved, 100, true, "  "};
    writer.writeAlignment(filename, *simulated_sequences, true);
}

map<int, double> Alignment::_vector_to_map(vector<double> vec) {
    map<int, double> m;
    size_t l = vec.size();
    for (size_t i = 0; i < l; ++i) {
        m[i] = vec[i];
    }
    return m;
}

void Alignment::_check_compatible_model(string datatype, string model) {
    bool incompat = false;
    if (((datatype == "dna") | (datatype == "nt")) & ((model == "JTT92") | (model == "JCprot") | (model == "DSO78") | (model == "WAG01") | (model == "LG08"))) {
        incompat = true;
    }
    else if (((datatype == "aa") | (datatype == "protein")) & ((model == "JCnuc") | (model == "JC69") | (model == "K80") | (model == "HKY85") | (model == "TN93" ) | (model == "GTR") | (model == "T92") | (model == "F84"))) {
        incompat = true;
    }
    if (incompat) {
        cerr << "Incompatible model (" << model << ") and datatype (" << datatype << ")" << endl;
        throw Exception("Incompatible model and datatype error");
    }
}

void Alignment::_clear_distances() {
    if (distances) {
        distances.reset();
    }
    if (variances) {
        variances.reset();
    }
}

void Alignment::_clear_likelihood() {
    if (likelihood) {
        likelihood.reset();
    }
}

/*
Calculate Jukes Cantor distance between sequences:
d = number of differences positions
g = ungapped length
s = number of states in the alphabet
NB - This is not defined when the log term is <= 0, so we set to DISTMAX in this case.
*/
double Alignment::_jcdist(double d, double g, double s) {
    double p;
    double dist;
    p = g > 0 ? d / g : 0;
    dist = (1 - (s/(s-1) * p)) > 0 ? - ((s-1)/s) * log(1 - (s/(s-1) * p)) : DISTMAX;
    return dist;
}

double Alignment::_jcvar(double d, double g, double s) {
    double p;
    double var;
    p = g > 0 ? d / g : 0;
    var = (p * (1 - p)) / (g * (pow(1 - ((s/(s-1)) * p), 2)));
    if (var < VARMIN) var = VARMIN;
    return var;
}

shared_ptr<DistanceMatrix> Alignment::_create_distance_matrix(vector<vector<double>> matrix) {
    size_t n = sequences->getNumberOfSequences();
    if (matrix.size() != n) throw Exception("Matrix wrong size error");
    if (distances) distances.reset();
    vector<string> names = sequences->getSequencesNames();
    shared_ptr<DistanceMatrix> dm = make_shared<DistanceMatrix>(names);
    for (vector<double> &row : matrix) {
        dm->addRow(row);
    }
    return dm;
}

bool Alignment::_is_file(string filename) {
    ifstream fl(filename.c_str());
    bool result = true;
    if (!fl) {
        result = false;
    }
    fl.close();
    return result;
}

bool Alignment::_is_tree_string(string tree_string) {
    size_t l = tree_string.length();
    return (tree_string[0]=='(' && tree_string[l-1]==';');
}
