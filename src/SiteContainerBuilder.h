/*
 * SiteContainerBuilder.h
 *
 *  Created on: Jul 20, 2014
 *      Author: kgori
 */

#ifndef SITECONTAINERBUILDER_H_
#define SITECONTAINERBUILDER_H_

#include <memory>
#include <string>
#include <vector>
#include <Bpp/Exceptions.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/VectorSequenceContainer.h>
using namespace bpp;
using namespace std;

class SiteContainerBuilder {
public:
    SiteContainerBuilder();
    virtual ~SiteContainerBuilder();
    static shared_ptr<VectorSiteContainer> read_alignment(string filename, string file_format, bool interleaved=true) throw (Exception);
    static shared_ptr<VectorSiteContainer> read_alignment(string filename, string file_format, string datatype, bool interleaved=true) throw (Exception);
    static shared_ptr<VectorSiteContainer> construct_alignment_from_strings(vector<pair<string, string>> headers_sequences, string datatype) throw (Exception);
    static shared_ptr<VectorSiteContainer> construct_sorted_alignment(VectorSiteContainer *sites, bool ascending);
    static shared_ptr<VectorSiteContainer> concatenate_alignments(vector<shared_ptr<VectorSiteContainer>> vec_of_vsc);
private:
    static bool asking_for_fasta(string file_format);
    static bool asking_for_phylip(string file_format);
    static bool asking_for_dna(string datatype);
    static bool asking_for_protein(string datatype);
    static shared_ptr<VectorSiteContainer> read_fasta_dna_file(string filename);
    static shared_ptr<VectorSiteContainer> read_fasta_protein_file(string filename);
    static shared_ptr<VectorSiteContainer> read_phylip_dna_file(string filename, bool interleaved);
    static shared_ptr<VectorSiteContainer> read_phylip_protein_file(string filename, bool interleaved);
    static shared_ptr<BasicSequence> _convert_pair_to_dna_sequence(pair<string, string>);
    static shared_ptr<BasicSequence> _convert_pair_to_protein_sequence(pair<string, string>);
    static shared_ptr<VectorSequenceContainer> _convert_vector_to_dna_sequence_container(vector<pair<string, string>> v);
    static shared_ptr<VectorSequenceContainer> _convert_vector_to_protein_sequence_container(vector<pair<string, string>> v);
};

#endif /* SITECONTAINERBUILDER_H_ */
