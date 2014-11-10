/*
 * Alignment.cpp
 *
 *  Created on: Jul 20, 2014
 *      Author: kgori
 */

#include "Alignment.h"
#include "SiteContainerBuilder.h"
#include "ModelFactory.h"

#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Phyl/Model/FrequenciesSet/NucleotideFrequenciesSet.h>
#include <Bpp/Phyl/Model/FrequenciesSet/ProteinFrequenciesSet.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Simulation/HomogeneousSequenceSimulator.h>
#include <Bpp/Phyl/Likelihood/NNIHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/TreeTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/CompressedVectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerIterator.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Io/Phylip.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/SymbolListTools.h>

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>

using namespace bpp;
using namespace std;

#define DISTMIN  0.000001
#define VARMIN  0.000001
#define DISTMAX  10000
#define MIN_BRANCH_LENGTH 0.000001

size_t getNumberOfDistinctPositionsWithoutGap(const SymbolList& l1, const SymbolList& l2) {
      if (l1.getAlphabet()->getAlphabetType() != l2.getAlphabet()->getAlphabetType()) throw AlphabetMismatchException("SymbolListTools::getNumberOfDistinctPositions.", l1.getAlphabet(), l2.getAlphabet());
      const Alphabet* alpha = l1.getAlphabet();
      int gapCode = alpha->getGapCharacterCode();
      size_t n = min(l1.size(), l2.size());
      size_t count = 0;
      for (size_t i = 0; i < n; i++) {
          int x = l1[i];
          int y = l2[i];
          if (alpha->isUnresolved(x)) x = gapCode;
          if (alpha->isUnresolved(x)) y = gapCode;
          if (x != gapCode && y != gapCode && x != y) count++;
      }
      return count;
}

void ensure_minval_and_sum(std::vector<double>& v, double minval) {
    double added = 0;
    double diff = 0;
    bool make_adjustment = false;
    std::vector<size_t> changes;

    for (size_t i = 0; i < v.size(); ++i) {
        if (v[i] < minval) {
            make_adjustment = true;
            diff = minval - v[i];
            v[i] = minval;
            added += diff;
        }
        else {
            changes.push_back(i);
        }
    }
    if (make_adjustment) {
    double to_subtract = added / changes.size();
    for (size_t i = 0; i < changes.size(); ++i) {
            v[changes[i]] -= to_subtract;
        }
    }
}

Alignment::Alignment() {}

Alignment::Alignment(vector<Alignment> alignments) {
    vector<shared_ptr<VectorSiteContainer>> vec_of_vsc;
    vec_of_vsc.reserve(alignments.size());
    for (auto al : alignments) {
        if (!al.sequences) throw Exception("At least one alignment has no sequences");
        vec_of_vsc.push_back(al.sequences);
    }
    sequences = SiteContainerBuilder::concatenate_alignments(vec_of_vsc);
}

Alignment::Alignment(vector<pair<string, string>> headers_sequences, string datatype) {
    sequences = SiteContainerBuilder::construct_alignment_from_strings(headers_sequences, datatype);
}

Alignment::Alignment(string filename, string file_format, bool interleaved) {
    read_alignment(filename, file_format, interleaved);
}

Alignment::Alignment(string filename, string file_format, string datatype, bool interleaved) {
    read_alignment(filename, file_format, datatype, interleaved);
}

Alignment::Alignment(string filename, string file_format, string datatype, string model_name, bool interleaved) {
    read_alignment(filename, file_format, datatype, interleaved);
    set_gamma_rate_model();
    try {
        _check_compatible_model(model_name);
        set_substitution_model(model_name);
    }
    catch (Exception& e) {
        cerr << "Alignment was initialised, but the model wasn't valid: " <<e.what() << endl;
    }
}

void Alignment::read_alignment(string filename, string file_format, bool interleaved) {
    sequences = SiteContainerBuilder::read_alignment(filename, file_format, interleaved);
    _clear_distances();
    _clear_likelihood();
}

void Alignment::read_alignment(string filename, string file_format, string datatype, bool interleaved) {
    sequences = SiteContainerBuilder::read_alignment(filename, file_format, datatype, interleaved);
    _clear_distances();
    _clear_likelihood();
}

void Alignment::sort_alignment(bool ascending) {
    if (!sequences) throw Exception("No sequences to sort");
    sequences = SiteContainerBuilder::construct_sorted_alignment(sequences.get(), ascending);
}

void Alignment::write_alignment(string filename, string file_format, bool interleaved) {
    if (file_format == "fas" || file_format == "fasta") {
        _write_fasta(sequences, filename);
    }
    else if (file_format == "phy" || file_format == "phylip") {
        _write_phylip(sequences, filename, interleaved);
    }
    else {
        cerr << "Unrecognised file format: " << file_format << endl;
        throw exception();
    }
}

void Alignment::set_substitution_model(string model_name) {
    _check_compatible_model(model_name);
    model = ModelFactory::create(model_name);
    _clear_likelihood();
}

void Alignment::set_gamma_rate_model(size_t ncat, double alpha) {
    if (ncat == 1) {
        cerr << "A discrete gamma distribution with 1 category is a constant distribution, so that's what I'm setting." << endl;
        set_constant_rate_model();
    }
    else {
        rates = make_shared<GammaDiscreteDistribution>(ncat, alpha, alpha);
        rates->aliasParameters("alpha", "beta");
    }
    _clear_likelihood();
}

void Alignment::set_constant_rate_model() {
    rates = make_shared<ConstantDistribution>(1.0);
}

void Alignment::set_alpha(double alpha) {
    if(!rates) throw Exception("No rate model is set");
    rates->setParameterValue("alpha", alpha);
    _clear_likelihood();
}

void Alignment::set_number_of_gamma_categories(size_t ncat) {
    if (!rates) throw Exception("No rate model is set");
    rates->setNumberOfCategories(ncat);
    _clear_likelihood();
}

void Alignment::set_rates(vector<double> rates, string order) {
    if (!model) throw Exception("Model not set");
    if (!is_dna() || model->getName() != "GTR") throw Exception("Setting rates is only implemented for DNA GTR model.");
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
        else throw Exception("Unrecognised order for rates: " + order);
        _clear_likelihood();
    }
}

void Alignment::set_frequencies(vector<double> freqs) {
    if (!model) throw Exception("Model not set");
    size_t reqd = is_dna() ? 4 : 20;
    if (freqs.size() != reqd) throw Exception("Frequencies vector is the wrong length (dna: 4; aa: 20)");
    ensure_minval_and_sum(freqs, 1.1e-6);
    map<int, double> m = _vector_to_map(freqs);
    if (is_protein()) {
        model = ModelFactory::create(model->getName(), freqs);
    }
    model->setFreq(m);
    _clear_likelihood();
}

void Alignment::set_namespace(string name) {
    if ((!rates) | (!model)) throw Exception("Substitution and rate models need to be fully set before adding a namespace");
    rates->setNamespace(name);
    model->setNamespace(name);
    //_clear_likelihood();
}

double Alignment::get_alpha() {
    if (rates) return rates->getParameterValue("alpha");
    else throw Exception("Gamma distributed rate model not set");
}

size_t Alignment::get_number_of_gamma_categories() {
    if (rates) return rates->getNumberOfCategories();
    else throw Exception("Rate model not set");
}

vector<double> Alignment::get_rates(string order) {
    if (is_dna()) {
        vector<double> rates_vec;
        RowMatrix<double> exch = model->getExchangeabilityMatrix();
        if (order == "acgt" || order == "ACGT") { //{a-c, a-g, a-t, c-g, c-t, g-t=1}
            double normaliser = exch(2,3);
            rates_vec.push_back(exch(0,1) / normaliser);
            rates_vec.push_back(exch(0,2) / normaliser);
            rates_vec.push_back(exch(0,3) / normaliser);
            rates_vec.push_back(exch(1,2) / normaliser);
            rates_vec.push_back(exch(1,3) / normaliser);
            rates_vec.push_back(1.0);
        }
        else if (order == "tcag" || order == "TCAG") { //{a=t-c, b=t-a, c=t-g, d=c-a, e=c-g, f=a-g=1}
            double normaliser = exch(0,2);
            rates_vec.push_back(exch(1,3) / normaliser);
            rates_vec.push_back(exch(0,3) / normaliser);
            rates_vec.push_back(exch(2,3) / normaliser);
            rates_vec.push_back(exch(0,1) / normaliser);
            rates_vec.push_back(exch(1,2) / normaliser);
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
    if (!model) throw Exception("Substitution model not set");
    return model->getFrequencies();
}

vector<double> Alignment::get_rate_model_categories() {
    if (!rates) throw Exception("Rate model not set");
    return rates->getCategories();
}

vector<string> Alignment::get_names() {
    if (!sequences) throw Exception("This instance has no sequences");
    return sequences->getSequencesNames();
}

size_t Alignment::get_number_of_sequences() {
    if (!sequences) throw Exception("This instance has no sequences");
    return sequences->getNumberOfSequences();
}

size_t Alignment::get_number_of_sites() {
    if (!sequences) throw Exception("This instance has no sequences");
    return sequences->getNumberOfSites();
}

size_t Alignment::get_number_of_distinct_sites() {
    if (!sequences) throw Exception("This instance has no sequences");
    SitePatterns sp{sequences.get()};
    auto unique_sites = sp.getSites();
    size_t distinct_sites = unique_sites->getNumberOfSites();
    delete unique_sites;
    return distinct_sites;
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

string Alignment::get_substitution_model() {
    if (!model) throw Exception("No model has been set");
    return model->getName();
}

string Alignment::get_namespace() {
    if (_name.empty()) throw Exception("No namespace is set");
    return _name;
}

vector<string> Alignment::get_sites() {
    if (!sequences) throw Exception("No sequences present.");
    vector<string> sites;
    auto si = new SimpleSiteContainerIterator(*sequences);
    while(si->hasMoreSites()) {
        sites.push_back(si->nextSite()->toString());
    }
    delete si;
    return sites;
}

vector<string> Alignment::get_informative_sites(bool exclude_gaps) {
    if (!sequences) throw Exception("No sequences present.");
    vector<string> inf_sites;
    ConstSiteIterator* si = nullptr;
    if (exclude_gaps) si = new CompleteSiteContainerIterator(*sequences);
    else si = new SimpleSiteContainerIterator(*sequences);
    const Site* site = 0;
    while (si->hasMoreSites()) {
        site = si->nextSite();
        if (SiteTools::isParsimonyInformativeSite(*site)) inf_sites.push_back(site->toString());
    }
    delete si;
    return inf_sites;
}

size_t Alignment::get_number_of_informative_sites(bool exclude_gaps) {
    return get_informative_sites(exclude_gaps).size();
}

size_t Alignment::get_number_of_free_parameters() {
    if (!likelihood) throw Exception("Likelihood model not initialised");
    ParameterList pl = likelihood->getBranchLengthsParameters();
    pl.addParameters(model->getIndependentParameters());
    if (rates->getName() == "Gamma") pl.addParameters(rates->getIndependentParameters());
    return pl.size();
}

void Alignment::_print_params() {
    if (likelihood) {
            ParameterList pl = likelihood->getParameters();
            pl.printParameters(cout);
    }
    else if (rates && model) {
         ParameterList pl = rates->getIndependentParameters();
         pl.addParameters(model->getIndependentParameters());
         pl.printParameters(cout);
         cout << "----------" << endl;
    }
}

bool Alignment::is_dna() {
    return _get_datatype() == "DNA alphabet";
}

bool Alignment::is_protein() {
    return _get_datatype() == "Proteic alphabet";
}

// Distance
void Alignment::compute_distances() {
    if (!sequences) throw Exception("This instance has no sequences");
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
    if (!sequences) throw Exception("This instance has no sequences");
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
            size_t d = getNumberOfDistinctPositionsWithoutGap(sequences->getSequence(i), sequences->getSequence(j));
            size_t g = SymbolListTools::getNumberOfPositionsWithoutGap(sequences->getSequence(i), sequences->getSequence(j));
//            cout << sequences->getSequence(i).toString() << endl << sequences->getSequence(j).toString() << endl;
            double dist = _jcdist(d, g, s);
            double var = _jcvar(d, g, s);
            (*distances)(i, j) = (*distances)(j, i) = dist;
            (*variances)(i, j) = (*variances)(j, i) = var;
        }
    }
}

void Alignment::set_distance_matrix(vector<vector<double>> matrix) {
    try {
        distances = _create_distance_matrix(matrix);
    }
    catch (Exception &e) {
        cout << e.what() << endl;
        throw Exception("Error setting distance matrix");
    }
}

void Alignment::chkdst() {
    if (!distances) throw Exception("No distances have been calculated yet");
    cout << "Dims = (" << distances->getNumberOfRows() << ", " << distances->getNumberOfColumns() << ")" << endl;
    for (size_t i=0; i < distances->getNumberOfRows(); ++i) {
        for (size_t j=0; j < distances->getNumberOfColumns(); ++j) {
            cout << (*distances)(i,j) << " ";
        }
        cout <<endl;
    }
    cout << endl;
}

string Alignment::get_bionj_tree() {
    if (!distances) throw Exception("No distances have been calculated yet");
    if (!variances) return _computeTree(*distances, *distances);
    return _computeTree(*distances, *variances);
}

string Alignment::get_bionj_tree(vector<vector<double>> matrix) {
    shared_ptr<DistanceMatrix> dm = _create_distance_matrix(matrix);
    return _computeTree(*dm, *dm);
}

vector<vector<double>> Alignment::get_distances() {
    if(!distances) throw Exception("No distances have been calculated yet");
    vector<vector<double>> vec;
    vector<string> names = sequences->getSequencesNames();
    size_t nrow = distances->getNumberOfRows();
    for (size_t i = 0; i < nrow; ++i) {
        vector<double> row;
        for (size_t j = 0; j < nrow; ++j) {
            row.push_back((*distances)(names[i], names[j]));
        }
        vec.push_back(row);
    }
    return vec;
}

vector<vector<double>> Alignment::get_variances() {
    if(!variances) throw Exception("No distances have been calculated yet");
    vector<vector<double>> vec;
    vector<string> names = sequences->getSequencesNames();
    size_t nrow = variances->getNumberOfRows();
    for (size_t i = 0; i < nrow; ++i) {
        vector<double> row;
        for (size_t j = 0; j < nrow; ++j) {
            row.push_back((*variances)(names[i], names[j]));
        }
        vec.push_back(row);
    }
    return vec;
}

vector<vector<double>> Alignment::get_distance_variance_matrix() {
    if(!variances || !distances) throw Exception("No distances have been calculated yet");
    vector<vector<double>> vec;
    vector<string> names = sequences->getSequencesNames();
    size_t nrow = variances->getNumberOfRows();
    for (size_t i = 0; i < nrow; ++i) {
        vector<double> row;
        for (size_t j = 0; j < nrow; ++j) {
            if (j < i) row.push_back((*variances)(names[i], names[j]));
            else row.push_back((*distances)(names[i], names[j]));
        }
        vec.push_back(row);
    }
    return vec;
}

// Likelihood
void Alignment::initialise_likelihood() {
    if (!distances) fast_compute_distances();
    try {
        initialise_likelihood(get_bionj_tree());
    }
    catch (Exception& e) {
        cerr << e.what();
        _clear_distances();
    }
}

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
    auto sites_ = new CompressedVectorSiteContainer(*sequences);
    SiteContainerTools::changeGapsToUnknownCharacters(*sites_);
    likelihood = make_shared<NNIHomogeneousTreeLikelihood>(*liktree, *sites_, model.get(), rates.get(), true, false);
    likelihood->initialize();
    delete liktree;
    delete sites_;
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
    OptimizationTools::optimizeNumericalParameters2(likelihood.get(), pl, 0, 0.001, 1000000, NULL, NULL, false, false, 10);
}

void Alignment::optimise_topology(bool fix_model_params) {
    if (!likelihood) {
        cerr << "Likelihood calculator not set - call initialise_likelihood" << endl;
        throw Exception("Uninitialised likelihood error");
    }
    ParameterList pl = likelihood->getBranchLengthsParameters();
    if (!fix_model_params) {
        pl.addParameters(model->getIndependentParameters());
        if (rates->getName() == "Gamma") pl.addParameters(rates->getIndependentParameters());
    }
    likelihood = make_shared<NNIHomogeneousTreeLikelihood>(*OptimizationTools::optimizeTreeNNI2(likelihood.get(), pl, true, 0.001, 0.1, 1000000, 1, NULL, NULL, false, 10));
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
        throw Exception("Likelihood calculator not set - call initialise_likelihood");
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
void Alignment::write_simulation(size_t nsites, string filename, string file_format, bool interleaved) {
    simulate(nsites);

    if (file_format == "fas" || file_format == "fasta") {
        _write_fasta(simulated_sequences, filename);
    }
    else if (file_format == "phy" || file_format == "phylip") {
        _write_phylip(simulated_sequences, filename, interleaved);
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

vector<pair<string, string>> Alignment::simulate(size_t nsites, string tree) {
    set_simulator(tree);
    return simulate(nsites);
}

vector<pair<string, string>> Alignment::simulate(size_t nsites) {
    if (!simulator) {
        cout << "Tried to simulate without a simulator" << endl;
        throw exception();
    }
    SiteContainer * tmp = simulator->simulate(nsites);
    simulated_sequences = make_shared<VectorSiteContainer>(*tmp);
    /* For future use: to mask gaps in a simulated alignment:

    auto alphabet = sequences->getAlphabet();
    int gapCode = alphabet->getGapCharacterCode();
    size_t numseq = this->get_number_of_sequences();
    for (unsigned int i = 0; i < numseq; i++) {
        BasicSequence oldseq = sequences->getSequence(i);
        BasicSequence newseq = sim_sites->getSequence(i);
        for (unsigned int j = 0; j < number_of_sites; j++) {
            if (alphabet->isGap(oldseq[j])) {
                newseq.setElement(j, gapCode);
            }
        }
        string name = newseq.getName();
        sim_sites->setSequence(name, newseq, true);
    }

    */
    delete tmp;
    return get_simulated_sequences();
}

vector<pair<string, string>> Alignment::get_sequences() {
    if (!sequences) throw Exception("No sequences to return");
    return _get_sequences(sequences.get());
}

vector<pair<string, string>> Alignment::get_simulated_sequences() {
    if (!simulated_sequences) throw Exception("No simulated sequences to return");
    return _get_sequences(simulated_sequences.get());
}

// Bootstrap
vector<pair<string, string>> Alignment::get_bootstrapped_sequences() {
    if (!sequences) throw Exception("No sequences to bootstrap.");
    VectorSiteContainer *tmp = SiteContainerTools::bootstrapSites(*sequences);
    auto ret = _get_sequences(tmp);
    delete tmp;
    return ret;
}

// Misc
string Alignment::get_mrp_supertree(vector<string> trees) {
    vector<Tree*> input_trees;
    stringstream ss;
    unique_ptr<Newick> newickIO(new Newick(false));

    for (string tree : trees) {
        ss.str(tree);
        ss.clear();
        input_trees.push_back(newickIO->read(ss));
    }

    ss.str(string());
    ss.clear();

    auto mrptree = TreeTools::MRP(input_trees);
    newickIO->write(*mrptree, ss);
    delete mrptree;

    string s{ss.str()};
    s.erase(s.find_last_not_of(" \n\r\t")+1);

    return s;
}

// Private methods
string Alignment::_get_datatype() {
    if (!sequences) throw Exception("This instance has no sequences");
    return sequences->getAlphabet()->getAlphabetType();
}

vector<pair<string, string>> Alignment::_get_sequences(VectorSiteContainer *seqs) {
    vector<pair<string, string>> ret;
    if (!seqs) {
        cerr << "Empty sequences pointer" << endl;
        throw exception();
    }
    for (size_t i = 0; i < seqs->getNumberOfSequences(); ++i) {
        BasicSequence seq = seqs->getSequence(i);
        ret.push_back(make_pair(seq.getName(), seq.toString()));
    }
    return ret;
}

void Alignment::_write_fasta(shared_ptr<VectorSiteContainer> seqs, string filename) {
    Fasta writer;
    writer.writeAlignment(filename, *seqs);
}

void Alignment::_write_phylip(shared_ptr<VectorSiteContainer> seqs, string filename, bool interleaved) {
    Phylip writer{true, !interleaved, 100, true, "  "};
    writer.writeAlignment(filename, *seqs, true);
}

map<int, double> Alignment::_vector_to_map(vector<double> vec) {
    map<int, double> m;
    size_t l = vec.size();
    for (size_t i = 0; i < l; ++i) {
        m[i] = vec[i];
    }
    return m;
}

void Alignment::_check_compatible_model(string model) {
    bool incompat = false;
    if (is_dna() & ((model == "JTT92") | (model == "JCprot") | (model == "DSO78") | (model == "WAG01") | (model == "LG08"))) {
        incompat = true;
    }
    else if (is_dna() & ((model == "JTT92+F") | (model == "JCprot+F") | (model == "DSO78+F") | (model == "WAG01+F") | (model == "LG08+F"))) {
        incompat = true;
    }
    else if (is_protein() & ((model == "JCnuc") | (model == "JC69") | (model == "K80") | (model == "HKY85") | (model == "TN93" ) | (model == "GTR") | (model == "T92") | (model == "F84"))) {
        incompat = true;
    }
    if (incompat) {
        string datatype;
        if (is_dna()) datatype = "dna";
        else datatype = "protein";
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
//    cout << "d,g,s,p,dist : " << d << " " << g << " " << s << " " << p << " " << dist << endl;
    return dist > DISTMIN ? dist : DISTMIN;
}

double Alignment::_jcvar(double d, double g, double s) {
    double p;
    double var;
    p = g > 0 ? d / g : 0;
    var = (p * (1 - p)) / (g * (pow(1 - ((s/(s-1)) * p), 2)));
    return var > VARMIN ? var : VARMIN;
}

shared_ptr<DistanceMatrix> Alignment::_create_distance_matrix(vector<vector<double>> matrix) {
    if (!sequences) throw Exception("This instance has no sequences");
    size_t n = sequences->getNumberOfSequences();
    if (matrix.size() != n) throw Exception("Matrix wrong size error");
    vector<string> names = sequences->getSequencesNames();
    shared_ptr<DistanceMatrix> dm = make_shared<DistanceMatrix>(names);
    for (size_t i=0; i < matrix.size(); ++i) {
        auto row = matrix[i];
        for (size_t j=i+1; j < row.size(); ++j) {
            (*dm)(i, j) = (*dm)(j, i) = row[j];
        }
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

string Alignment::_computeTree(DistanceMatrix dists, DistanceMatrix vars) throw (Exception) {
    // Initialization:
    std::map<size_t, Node*> currentNodes_;
    std::vector<double> sumDist_(dists.size());
    double lambda_;

    for (size_t i = 0; i < dists.size(); i++) {
        currentNodes_[i] = new Node(static_cast<int>(i), dists.getName(i));
    }
    int idNextNode = dists.size();
    vector<double> newDist(dists.size());
    vector<double> newVar(dists.size());

    // Build tree:
    while (currentNodes_.size() > 3) {
        // get best pair
        for (std::map<size_t, Node*>::iterator i = currentNodes_.begin(); i != currentNodes_.end(); i++) {
            size_t id = i->first;
            sumDist_[id] = 0;
            for (map<size_t, Node*>::iterator j = currentNodes_.begin(); j != currentNodes_.end(); j++) {
                size_t jd = j->first;
                sumDist_[id] += dists(id, jd);
            }
        }
        vector<size_t> bestPair(2);
        double critMax = std::log(0.);
        for (map<size_t, Node*>::iterator i = currentNodes_.begin(); i != currentNodes_.end(); i++) {
            size_t id = i->first;
            map<size_t, Node*>::iterator j = i;
            j++;
            for ( ; j != currentNodes_.end(); j++) {
                size_t jd = j->first;
                double crit = sumDist_[id] + sumDist_[jd] - static_cast<double>(currentNodes_.size() - 2) * dists(id, jd);
                // cout << "\t" << id << "\t" << jd << "\t" << crit << endl;
                if (crit > critMax) {
                    critMax = crit;
                    bestPair[0] = id;
                    bestPair[1] = jd;
                }
            }
        }
        if (critMax == std::log(0.)) throw Exception("Unexpected error: no maximum criterium found.");

        // get branch lengths for pair
        double ratio = (sumDist_[bestPair[0]] - sumDist_[bestPair[1]]) / static_cast<double>(currentNodes_.size() - 2);
        vector<double> d(2);

        d[0] = std::max(.5 * (dists(bestPair[0], bestPair[1]) + ratio), MIN_BRANCH_LENGTH);
        d[1] = std::max(.5 * (dists(bestPair[0], bestPair[1]) - ratio), MIN_BRANCH_LENGTH);

        Node* best1 = currentNodes_[bestPair[0]];
        Node* best2 = currentNodes_[bestPair[1]];

        // Distances may be used by getParentNodes (PGMA for instance).
        best1->setDistanceToFather(d[0]);
        best2->setDistanceToFather(d[1]);
        Node* parent = new Node(idNextNode++);
        parent->addSon(best1);
        parent->addSon(best2);

        // compute lambda
        lambda_ = 0;
        if (vars(bestPair[0], bestPair[1]) == 0) lambda_ = .5;
        else {
            for (map<size_t, Node*>::iterator i = currentNodes_.begin(); i != currentNodes_.end(); i++) {
                size_t id = i->first;
                if (id != bestPair[0] && id != bestPair[1]) lambda_ += (vars(bestPair[1], id) - vars(bestPair[0], id));
            }
            double div = 2 * static_cast<double>(currentNodes_.size() - 2) * vars(bestPair[0], bestPair[1]);
            lambda_ /= div;
            lambda_ += .5;
        }
        if (lambda_ < 0.) lambda_ = 0.;
        if (lambda_ > 1.) lambda_ = 1.;

        for (map<size_t, Node*>::iterator i = currentNodes_.begin(); i != currentNodes_.end(); i++) {
            size_t id = i->first;
            if (id != bestPair[0] && id != bestPair[1]) {
                newDist[id] = std::max(lambda_ * (dists(bestPair[0], id) - d[0]) + (1 - lambda_) * (dists(bestPair[1], id) - d[1]), 0.);
                newVar[id] = lambda_ * vars(bestPair[0], id) + (1 - lambda_) * vars(bestPair[1], id) - lambda_ * (1 - lambda_) * vars(bestPair[0], bestPair[1]);
            }
          else newDist[id] = 0;
        }
        // Actualize currentNodes_:
        currentNodes_[bestPair[0]] = parent;
        currentNodes_.erase(bestPair[1]);
        for (map<size_t, Node*>::iterator i = currentNodes_.begin(); i != currentNodes_.end(); i++) {
            size_t id = i->first;
            dists(bestPair[0], id) = dists(id, bestPair[0]) = newDist[id];
            vars(bestPair[0], id) =  vars(id, bestPair[0]) = newVar[id];
        }
    }
    // final step
    Node* root = new Node(idNextNode);
    map<size_t, Node* >::iterator it = currentNodes_.begin();
    size_t i1 = it->first;
    Node* n1       = it->second;
    it++;
    size_t i2 = it->first;
    Node* n2       = it->second;
    if (currentNodes_.size() == 2) {
        // Rooted
        double d = dists(i1, i2) / 2;
        root->addSon(n1);
        root->addSon(n2);
        n1->setDistanceToFather(d);
        n2->setDistanceToFather(d);
    }
    else {
        // Unrooted
        it++;
        size_t i3 = it->first;
        Node* n3       = it->second;
        double d1 = std::max(dists(i1, i2) + dists(i1, i3) - dists(i2, i3), MIN_BRANCH_LENGTH);
        double d2 = std::max(dists(i2, i1) + dists(i2, i3) - dists(i1, i3), MIN_BRANCH_LENGTH);
        double d3 = std::max(dists(i3, i1) + dists(i3, i2) - dists(i1, i2), MIN_BRANCH_LENGTH);
        root->addSon(n1);
        root->addSon(n2);
        root->addSon(n3);
        n1->setDistanceToFather(d1 / 2.);
        n2->setDistanceToFather(d2 / 2.);
        n3->setDistanceToFather(d3 / 2.);
    }
    Tree *tree_ = new TreeTemplate<Node>(root);
    stringstream ss;
    Newick treeWriter;
    if (!tree_) throw Exception("The tree is empty");
    treeWriter.write(*tree_, ss);
    delete tree_;
    string s{ss.str()};
    s.erase(s.find_last_not_of(" \n\r\t")+1);
    return s;
}
