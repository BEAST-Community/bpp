#include "AlignmentList.h"

using namespace std;
using namespace bpp;

AlignmentList::AlignmentList() {}

AlignmentList::AlignmentList(vector<Alignment> alignments) {
    _alignments = vector<Alignment>(alignments);
}

Alignment AlignmentList::operator[](const size_t idx) {return _alignments[idx];}

void Alignment::initialise_likelihood(string tree) {

    VectorSiteContainer* sites_ = sequences->clone();
    SiteContainerTools::changeGapsToUnknownCharacters(*sites_);
    likelihood = make_shared<NNIHomogeneousTreeLikelihood>(*liktree, *sites_, model.get(), rates.get(), true, true);
    likelihood->initialize();
    delete liktree;
}