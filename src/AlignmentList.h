#ifndef _ALIGNMENT_LIST_H_
#define _ALIGNMENT_LIST_H_

#include "Alignment.h"

using namespace std;
using namespace bpp;

class AlignmentList {
    public :
        AlignmentList();
        AlignmentList(vector<Alignment>);
        virtual ~AlignmentList();
        void initialise_likelihood();
        void optimise_parameters(bool fix_branch_lengths);
        void optimise_topology(bool fix_model_params);
        Alignment& operator[](size_t idx);

    private :
        vector<Alignment> _alignments;
        shared_ptr<NNIHomogeneousTreeLikelihood> likelihood = nullptr;
};


#endif // _ALIGNMENT_LIST_H_