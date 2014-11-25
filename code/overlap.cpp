#include "Alignment.h"
#include "Needleman_Wunsch.h"
#include "Smith_Waterman.h"
#include "End1_Beginning2.h"
#include "Local1_Global2.h"
#include "End1_local2.h"
#include "String.h"
#include "common.h"
#include "SeqIORead_fasta.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <map>
#include <vector>

using namespace std;

double VERSION = 1.01;

namespace Bio {
	typedef Alignment<Bio::DNAString> DNA_Alignment;
	typedef Needleman_Wunsch<Bio::DNAString> DNA_NW;
	typedef Smith_Waterman<Bio::DNAString> DNA_SW;

/***********************************************************************************************************************************/
class DNAScoring : public Bio::Alignment<Bio::DNAString>::Scoring {
public:
	DNAScoring(double _gap=-2, double _match=2, double _mismatch=-3) : gap(_gap), match(_match), mismatch(_mismatch)	{} 
	double score(char ch1, char ch2) const 
	{
		if((ch1 == '-') || (ch2 == '-'))
			return gap;
		else if((ch1 == 'N') || (ch2 == 'N') || (ch1 != ch2))
			return mismatch;
		// else if(ch1 == ch2)
		return match;
	}
	void	set_gap(double _gap)			{gap = _gap;}
	double	get_gap() const				{return gap;}
	void	set_mismatch(double _mismatch)		{mismatch = _mismatch;}
	double	get_mismatch() const			{return mismatch;}
	void	set_match(double _match)		{match = _match;}
	double	get_match() const			{return match;}
protected:
	double	gap, match, mismatch;
};

}

/***********************************************************************************************************************************/
bool read_params(int argc, const char* argv[]);
void usage(const char* prog_name);
void write_alignment(const Bio::Alignment<Bio::DNAString>& aligner, const Bio::DNASequence& seq1, Bio::DNASequence& seq2, bool seq2_strand_is_plus);

/***********************************************************************************************************************************/
const char* seq_file = NULL;
const char* m8_blast_file = NULL;
Bio::Alignment<Bio::DNAString>* aligner = NULL;
enum Alignment_type {global=1111, local=2222} alignment_type = global;
Bio::DNAScoring score(-5, 1, -3);
size_t	overlap_lenght_threshold = 500;
size_t	seed_lenght_threshold = 300;
double	pidentity_threshold = 99;

typedef enum{Forward, Reverse} Direction;

/***********************************************************************************************************************************/
struct Pair_info {
	Pair_info() : seq1(""), seq2(""), s1(0), e1(0), s2(0), e2(0), pidentity(0), dir(Forward) {}
	Pair_info(string _seq1, string _seq2, size_t _s1, size_t _e1, size_t _s2, size_t _e2, double _pidentity) : seq1(_seq1), seq2(_seq2), 
		s1(_s1), e1(_e1), s2(_s2), e2(_e2), pidentity(_pidentity), dir((_s2<_e2)? Forward : Reverse) 
	{if(s2 > e2) {s2=_e2; e2=_s2;}}
	string		seq1;
	string		seq2;
	size_t		s1, e1, s2, e2;
	double		pidentity;
	Direction	dir;
};

/***********************************************************************************************************************************/
int main(int argc, const char* argv[])
{
	cerr << argv[0] << " v" << VERSION << endl;
	if(!read_params(argc, argv))
		return 1;

	/*------------- Step 1: extract all pairs to check from the m8 file -------------*/
	cerr << time_since_epoch() << " Reading " << m8_blast_file << endl;
	map<string, Bio::DNASequence*> sequences;
	vector<Pair_info>	pairs;
	FILE* fp = fopen(m8_blast_file, "r");
	if(fp == NULL) {
		cerr << endl << "Could not read file " << m8_blast_file << endl << endl;
		return 1;
	}

	string	line;
	map<string, int> seen_pair;

	while(read_line(line, fp)) {
		// mol-32-15fa-028585	mol-32-15fa-049805	98.18	7857	141	2	4122	11977	570	8425	0.0	1.403e+04
		char	seq1[1024], seq2[1024];
		double	pidentity=0;
		size_t	nidentical=0, nmismatches=0, ngaps=0, s1=0, e1=0, s2=0, e2=0;
		if(sscanf(line.c_str(), "%s %s %lf %ld %ld %ld %ld %ld %ld %ld", seq1, seq2, &pidentity, &nidentical, &nmismatches, &ngaps, &s1, &e1, &s2, &e2) != 10) {
			cerr << endl << "Unexpected line in an m8 blast file:" << endl << line << endl << endl;
			cerr << seq1 << '\t' << seq2 << '\t' << pidentity << '\t' << nidentical << '\t' << nmismatches << '\t' << ngaps << '\t' << s1 
			     << '\t' << e1 << '\t' << s2 << '\t' << e2 << endl;
			return 1;
		}

		// Each qualifying pair will be seen twice - we arbitrarily pick one
		
//		if(strcmp(seq1, seq2) >= 0)
//			continue;
		if((pidentity < pidentity_threshold) || (abs(e1-s1)+1 < seed_lenght_threshold) || (abs(e2-s2)+1 < seed_lenght_threshold)) 
			continue;
		string id = string(seq1) + string(":") + string(seq2);
		if(seen_pair.find(id) != seen_pair.end()) 
			continue;
	
		sequences[seq1] = NULL;
		sequences[seq2] = NULL;
		pairs.push_back(Pair_info(seq1, seq2, s1, e1, s2, e2, pidentity));
		seen_pair.insert(pair<string, int>(id, 0));
	}
	fclose(fp);
	cerr << time_since_epoch() << " ok, identified " << pairs.size() << " pairs to check over " << sequences.size() << " sequences"  << endl;

	/*------------- Step 2: Read all sequences, keep those that are needed -------------*/
	cerr << time_since_epoch() << " Reading sequence file" << endl;
	Bio::SeqIORead_fasta<Bio::DNASequence>	reader(seq_file);
	Bio::DNASequence*	seq_obj;
	while(seq_obj = reader.next_seq()) {
		map<string, Bio::DNASequence*>::iterator it = sequences.find(seq_obj->display_id());
		if(it != sequences.end()) {
			it->second = seq_obj;
		}
	}
	cerr << time_since_epoch() << " ok" << endl;

	/*------------- Step 3: go over all pairs -------------*/
	cerr << time_since_epoch() << " Analyzing pairs" << endl;
	Bio::Local1_Global2<Bio::DNAString>* 	l1g2 = new Bio::Local1_Global2<Bio::DNAString>();
	Bio::End1_Beginning2<Bio::DNAString>*	edges = new Bio::End1_Beginning2<Bio::DNAString>(overlap_lenght_threshold);
	Bio::End1_Local2<Bio::DNAString>*	edge_local = new Bio::End1_Local2<Bio::DNAString>();
	Bio::Smith_Waterman<Bio::DNAString>*	local =  new Bio::Smith_Waterman<Bio::DNAString>();
	int count = 0;
	cerr << time_since_epoch() << " 0%" << endl;
	for(vector<Pair_info>::const_iterator it=pairs.begin(); it!=pairs.end(); it++) {
		if(int(1000.0*(it-pairs.begin())/pairs.size()) > count) {
			count++;
			cerr << time_since_epoch() << " " << int(1000*(it-pairs.begin())/pairs.size())/10.0 << "%" << endl;
		}

		map<string, Bio::DNASequence*>::const_iterator its = sequences.find(it->seq1);
		if(its == sequences.end()) {
			cerr << "Error: sequence " << it->seq1 << " appears in " << m8_blast_file << " but not in " << seq_file << endl << endl;
			return 1;
		}
		Bio::DNASequence* seq1_obj = its->second;
		its = sequences.find(it->seq2);
		if(its == sequences.end()) {
			cerr << "Error: sequence " << it->seq2 << " appears in " << m8_blast_file << " but not in " << seq_file << endl << endl;
			return 1;
		}

		Bio::DNASequence* seq2_obj = its->second;
		// First predict the region that needs to be aligned
		size_t align_s1, align_e1;
		size_t align_s2, align_e2;

		size_t gap1, gap2;

		// align_s1
		gap1 = it->s1-1;
		gap2 = (it->dir == Forward)? (it->s2-1) : (seq2_obj->seq().size()-it->e2);
		align_s1 = (gap1 < gap2)? 1 : (it->s1-gap2);

		// align_s2
		gap2 = it->s2-1;
		gap1 = (it->dir == Forward)? (it->s1-1) : (seq1_obj->seq().size()-it->e1);
		align_s2 = (gap2 < gap1)? 1 : (it->s2-gap1);

		// align_e1
		gap1 = seq1_obj->seq().size()-it->e1;
		gap2 = (it->dir == Forward)? (seq2_obj->seq().size()-it->e2) : (it->s2-1);
		align_e1 = (gap1 < gap2)? seq1_obj->seq().size() : (it->e1+gap2);

		// align_e2
		gap2 = seq2_obj->seq().size()-it->e2;
		gap1 = (it->dir == Forward)? (seq1_obj->seq().size()-it->e1) : (it->s1-1);
		align_e2 = (gap2 < gap1)? seq2_obj->seq().size() : (it->e2+gap1);

		size_t grace_bp = (align_e1-align_s1)/100;	// 11/12/13: change this to allow much less freedom
		size_t cs1 = (align_s1<grace_bp)? 0 : (align_s1-grace_bp);
		size_t ce1 = (align_e1+grace_bp >= seq1_obj->seq().size())? (seq1_obj->seq().size()-1) : (align_e1+grace_bp);
		size_t cs2 = (align_s2<grace_bp)? 0 : (align_s2-grace_bp);
		size_t ce2 = (align_e2+grace_bp >= seq2_obj->seq().size())? (seq2_obj->seq().size()-1) : (align_e2+grace_bp);

		double pidentity1=0;
		size_t s1, e1, s2, e2, alignment_size;

		// First, try to see if seq2 is contained within seq1 or vice-versa (or both)
		// If the following is true then seq2 can be contained within seq1
		if((cs2 == 0) && (ce2 == (seq2_obj->seq().size()-1))) {
			// If the blast result already indicate that seq2 is contained in seq1 - no need to align
			if((it->s2 == 1) && (it->e2 == seq2_obj->seq().size())) {
				s1 = it->s1;
				e1 = it->e1;
				s2 = (it->dir == Forward)? it->s2 : it->e2;
				e2 = (it->dir == Forward)? it->e2 : it->s2;
				pidentity1 = it->pidentity;
				alignment_size = e1-s1+1;
			}
			else {
				if(it->dir == Forward) {
					l1g2->align(seq1_obj->seq().subseq(cs1, ce1), seq2_obj->seq(), score);
				}
				else {
					l1g2->align(seq1_obj->seq().subseq(cs1, ce1), seq2_obj->seq().reverse_complement(), score);
				}
				pidentity1 = 100.0*l1g2->num_identical()/l1g2->alignment_size();
	
				if((pidentity1 >= pidentity_threshold) && (l1g2->alignment_size() >= overlap_lenght_threshold)) {
					s1 = l1g2->seq1_start()+cs1+1;
					e1 = l1g2->seq1_end()+cs1+1;
					s2 = (it->dir == Forward)? (l1g2->seq2_start()+1) : (seq2_obj->seq().size()-l1g2->seq2_start());
					e2 = (it->dir == Forward)? (l1g2->seq2_end()+1) : (seq2_obj->seq().size()-l1g2->seq2_end());
					alignment_size = l1g2->alignment_size();
				}
			}
		}

		// Next, check the other way around
		double pidentity2 = 0;
		size_t os1, oe1, os2, oe2, oalignment_size;

		if((cs1 == 0) && (ce1 == (seq1_obj->seq().size()-1))) {
			// If the blast result already indicate that seq2 is contained in seq1 - no need to align
			if((it->s1 == 1) && (it->e1 == seq1_obj->seq().size())) {
				os1 = it->s1;
				oe1 = it->e1;
				os2 = (it->dir == Forward)? it->s2 : it->e2;
				oe2 = (it->dir == Forward)? it->e2 : it->s2;
				pidentity2 = it->pidentity;
				oalignment_size = oe1-os1+1;
			}
			else {
				if(it->dir == Forward) {
					l1g2->align(seq2_obj->seq().subseq(cs2, ce2), seq1_obj->seq(), score);
				}
				else {
					l1g2->align(seq2_obj->seq().subseq(cs2, ce2), seq1_obj->seq().reverse_complement(), score);
				}
				pidentity2 = 100.0*l1g2->num_identical()/l1g2->alignment_size();
				if((pidentity2 >= pidentity_threshold) && (l1g2->alignment_size() >= overlap_lenght_threshold)) {
					os1 = (it->dir == Forward)? (l1g2->seq2_start()+1) : (seq1_obj->seq().size()-l1g2->seq2_start());
					oe1 = (it->dir == Forward)? (l1g2->seq2_end()+1) : (seq1_obj->seq().size()-l1g2->seq2_end());
					os2 = (l1g2->seq1_start()+cs2+1);
					oe2 = (l1g2->seq1_end()+cs2+1);
					oalignment_size = l1g2->alignment_size();
				}
			}
		}
		if((pidentity1 >= pidentity_threshold) && (alignment_size >= overlap_lenght_threshold)) {
		 	const char* type = ((pidentity2 >= pidentity_threshold) && (oalignment_size >= overlap_lenght_threshold))? 
					"IDENTICAL" : "CONTAINS";
			// CONTAINS	mol-32-15fa-031894	mol-32-15fa-009682$s2	4137	722	4858	4137	4137	9296	99
			cout 	<< type << '\t' << seq1_obj->display_id() << '\t' << seq2_obj->display_id() << '\t' << s1 << '\t' << e1 
				<< '\t' << s2 << '\t' << e2 << '\t' << seq1_obj->seq().size() << '\t' << seq2_obj->seq().size() << '\t' << pidentity1 << endl;
			continue; 
		}
		else if((pidentity2 >= pidentity_threshold) && (oalignment_size >= overlap_lenght_threshold)) {
			cout 	<< "CONTAINS" << '\t' << seq2_obj->display_id() << '\t' << seq1_obj->display_id() << '\t' << os2 << '\t' 
				<< oe2 << '\t' << os1 << '\t' << oe1 << '\t' << seq2_obj->seq().size() << '\t' << seq1_obj->seq().size() 
				<< '\t' << pidentity2 << endl;
			continue; 
		}

		// If we are here then no containment was found. Look for edge-edge overlap

		// First: seq1_obj:3' ...
		pidentity1 = 0;

		if((it->e1 == seq1_obj->seq().size()) && (((it->s2 == 1) && (it->dir == Forward)) || ((it->e2 == seq2_obj->seq().size()) && (it->dir == Reverse)))) {
			// Nothing to alignment - blast alignment is good enough
			s1 = it->s1;
			e1 = it->e1;
			s2 = (it->dir == Forward)? it->s2 : it->e2;
			e2 = (it->dir == Forward)? it->e2 : it->s2;
			pidentity1 = it->pidentity;
			alignment_size = e1-s1+1;
		} 
		else if((it->dir == Forward) && (cs2 == 0) && (ce1 == (seq1_obj->seq().size()-1))) {
			edges->align(seq1_obj->seq().subseq(cs1, (seq1_obj->seq().size()-1)), seq2_obj->seq().subseq(0, ce2), score);
			pidentity1 = 100.0*edges->num_identical()/edges->alignment_size();
			s1 = edges->seq1_start()+1+cs1;
			e1 = edges->seq1_end()+1+cs1;
			s2 = edges->seq2_start()+1;
			e2 = edges->seq2_end()+1;
			alignment_size = edges->alignment_size();
		}
		else if((it->dir == Reverse) && (ce2 == (seq2_obj->seq().size()-1)) && (ce1 == (seq1_obj->seq().size()-1))) {
			edges->align(seq1_obj->seq().subseq(cs1, (seq1_obj->seq().size()-1)), seq2_obj->seq().subseq(cs2, (seq2_obj->seq().size()-1)).reverse_complement(), score);
			pidentity1 = 100.0*edges->num_identical()/edges->alignment_size();
			s1 = edges->seq1_start()+1+cs1;
			e1 = edges->seq1_end()+1+cs1;
			s2 = seq2_obj->seq().size()-edges->seq2_start();
			e2 = seq2_obj->seq().size()-edges->seq2_end();
			alignment_size = edges->alignment_size();
		}
		if((pidentity1 >= pidentity_threshold) && (alignment_size >= overlap_lenght_threshold)) {
			// CONNECTED	mol-32-15fa-013605	3	mol-32-15fa-045455	5	7297	8388	1	1092	8388	8885	100
			int side2 = (it->dir == Forward)? 5 : 3;
			
			cout 	<< "CONNECTED" << '\t' << seq1_obj->display_id() << "\t3\t" << seq2_obj->display_id() << '\t' << side2 << '\t'
				<< s1 << '\t' << e1 << '\t' << s2 << '\t' << e2 << '\t' << seq1_obj->seq().size() << '\t' << seq2_obj->seq().size() 
				<< '\t' << pidentity1 << endl;
			continue;
		}

		// Next: seq1_obj:5'
		pidentity1 = 0;
		if((it->s1 == 0) && (((it->s2 == 0) && (it->dir == Reverse)) || ((it->e2 == seq2_obj->seq().size()) && (it->dir == Forward)))) {
			// Nothing to alignment - blast alignment is good enough
			s1 = it->s1;
			e1 = it->e1;
			s2 = (it->dir == Forward)? it->s2 : it->e2;
			e2 = (it->dir == Forward)? it->e2 : it->s2;
			pidentity1 = it->pidentity;
			alignment_size = e1-s1+1;
		} 
		else if((it->dir == Reverse) && (cs2 == 0) && (cs1 == 0)) {
			edges->align(seq1_obj->seq().subseq(0, ce1).reverse_complement(), seq2_obj->seq().subseq(0, ce2), score);
			pidentity1 = 100.0*edges->num_identical()/edges->alignment_size();
			s1 = ce1-edges->seq1_start()+1;
			e1 = ce1-edges->seq1_end()+1;
			s2 = edges->seq2_start()+1;
			e2 = edges->seq2_end()+1;
			alignment_size = edges->alignment_size();
		}
		else if((it->dir == Forward) && (ce2 == (seq2_obj->seq().size()-1)) && (cs1 == 0)) {
			edges->align(seq1_obj->seq().subseq(0, ce1).reverse_complement(), seq2_obj->seq().subseq(cs2, (seq2_obj->seq().size()-1)).reverse_complement(), score);
			pidentity1 = 100.0*edges->num_identical()/edges->alignment_size();
			s1 = ce1-edges->seq1_start()+1;
			e1 = ce1-edges->seq1_end()+1;
			s2 = seq2_obj->seq().size()-edges->seq2_start();
			e2 = seq2_obj->seq().size()-edges->seq2_end();
			alignment_size = edges->alignment_size();
		}
		if((pidentity1 >= pidentity_threshold) && (alignment_size >= overlap_lenght_threshold)) {
			// CONNECTED	mol-32-15fa-013605	3	mol-32-15fa-045455	5	7297	8388	1	1092	8388	8885	100
			int side2 = (it->dir == Forward)? 3 : 5;
			
			cout 	<< "CONNECTED" << '\t' << seq1_obj->display_id() << "\t5\t" << seq2_obj->display_id() << '\t' << side2 << '\t'
				<< s1 << '\t' << e1 << '\t' << s2 << '\t' << e2 << '\t' << seq1_obj->seq().size() << '\t' 
				<< seq2_obj->seq().size() << '\t' << pidentity1 << endl;
			continue;
		}

		// If we are here then we need to check the option of end-middle connection
		// First: seq1_obj:3' ...
		pidentity1 = 0;

		if(it->e1 == seq1_obj->seq().size()) {
			// Nothing to alignment - blast alignment is good enough
			s1 = it->s1;
			e1 = it->e1;
			s2 = (it->dir == Forward)? it->s2 : it->e2;
			e2 = (it->dir == Forward)? it->e2 : it->s2;
			pidentity1 = it->pidentity;
			alignment_size = e1-s1+1;
		} 
		else if(ce1 == (seq1_obj->seq().size()-1)) {
			edge_local->align(seq1_obj->seq(), (it->dir == Forward)? seq2_obj->seq() : seq2_obj->seq().reverse_complement(), score);
			pidentity1 = 100.0*edge_local->num_identical()/edge_local->alignment_size();
			s1 = edge_local->seq1_start()+1;
			e1 = edge_local->seq1_end()+1;
			s2 = (it->dir == Forward)? edge_local->seq2_start()+1 : seq2_obj->seq().size()-edge_local->seq2_start();
			e2 = (it->dir == Forward)? edge_local->seq2_end()+1 : seq2_obj->seq().size()-edge_local->seq2_end();
			alignment_size = edge_local->alignment_size();
		}
		if((pidentity1 >= pidentity_threshold) && (alignment_size >= overlap_lenght_threshold)) {
			// CONNECTED	mol-32-15fa-013605	3	mol-32-15fa-045455	5	7297	8388	1	1092	8388	8885	100
			int side2 = (it->dir == Forward)? 5 : 3;
			cout 	<< "CONNECTED" << '\t' << seq1_obj->display_id() << "\t3\t" << seq2_obj->display_id() << '\t' << "Middle" << '\t'
				<< s1 << '\t' << e1 << '\t' << s2 << '\t' << e2 << '\t' << seq1_obj->seq().size() << '\t' << seq2_obj->seq().size() 
				<< '\t' << pidentity1 << endl;
			continue;
		}

		// seq1:5
		pidentity1 = 0;

		if(it->s1 == 1) {
			// Nothing to alignment - blast alignment is good enough
			s1 = it->s1;
			e1 = it->e1;
			s2 = (it->dir == Forward)? it->s2 : it->e2;
			e2 = (it->dir == Forward)? it->e2 : it->s2;
			pidentity1 = it->pidentity;
			alignment_size = e1-s1+1;
		} 
		else if(cs1 == 0) {
			edge_local->align(seq1_obj->seq().reverse_complement(), (it->dir == Forward)? seq2_obj->seq() : seq2_obj->seq().reverse_complement(), score);
			pidentity1 = 100.0*edge_local->num_identical()/edge_local->alignment_size();
			s1 = seq1_obj->seq().size()-edge_local->seq1_start();
			e1 = seq1_obj->seq().size()-edge_local->seq1_end();
			s2 = (it->dir == Forward)? edge_local->seq2_start()+1 : seq2_obj->seq().size()-edge_local->seq2_start();
			e2 = (it->dir == Forward)? edge_local->seq2_end()+1 : seq2_obj->seq().size()-edge_local->seq2_end();
			alignment_size = edge_local->alignment_size();
		}
		if((pidentity1 >= pidentity_threshold) && (alignment_size >= overlap_lenght_threshold)) {
			// CONNECTED	mol-32-15fa-013605	3	mol-32-15fa-045455	5	7297	8388	1	1092	8388	8885	100
			int side2 = (it->dir == Forward)? 5 : 3;
			
			cout 	<< "CONNECTED" << '\t' << seq1_obj->display_id() << "\t5\t" << seq2_obj->display_id() << '\t' << "Middle" << '\t'
				<< s1 << '\t' << e1 << '\t' << s2 << '\t' << e2 << '\t' << seq1_obj->seq().size() << '\t' << seq2_obj->seq().size() 
				<< '\t' << pidentity1 << endl;
			continue;
		}

		// Seq2:3'
		if(it->e2 == seq2_obj->seq().size()) {
			// Nothing to alignment - blast alignment is good enough
			s1 = it->s1;
			e1 = it->e1;
			s2 = (it->dir == Forward)? it->s2 : it->e2;
			e2 = (it->dir == Forward)? it->e2 : it->s2;
			pidentity1 = it->pidentity;
			alignment_size = e1-s1+1;
		} 
		else if(ce2 == (seq2_obj->seq().size()-1)) {
			edge_local->align(seq2_obj->seq(), (it->dir == Forward)? seq1_obj->seq() : seq1_obj->seq().reverse_complement(), score);
			pidentity1 = 100.0*edge_local->num_identical()/edge_local->alignment_size();
			s2 = edge_local->seq1_start()+1;
			e2 = edge_local->seq1_end()+1;
			s1 = (it->dir == Forward)? edge_local->seq2_start()+1 : seq1_obj->seq().size()-edge_local->seq2_start();
			e1 = (it->dir == Forward)? edge_local->seq2_end()+1 : seq1_obj->seq().size()-edge_local->seq2_end();
			alignment_size = edge_local->alignment_size();
		}
		if((pidentity1 >= pidentity_threshold) && (alignment_size >= overlap_lenght_threshold)) {
			// CONNECTED	mol-32-15fa-013605	3	mol-32-15fa-045455	5	7297	8388	1	1092	8388	8885	100
			int side2 = (it->dir == Forward)? 5 : 3;
			
			cout 	<< "CONNECTED" << '\t' << seq2_obj->display_id() << "\t3\t" << seq1_obj->display_id() << '\t' << "Middle" << '\t'
				<< s2 << '\t' << e2 << '\t' << s1 << '\t' << e1 << '\t' << seq2_obj->seq().size() << '\t' << seq1_obj->seq().size() 
				<< '\t' << pidentity1 << endl;
			continue;
		}

		// seq2:5
		pidentity1 = 0;
		if(it->s2 == 1) {
			// Nothing to alignment - blast alignment is good enough
			s1 = it->s1;
			e1 = it->e1;
			s2 = (it->dir == Forward)? it->s2 : it->e2;
			e2 = (it->dir == Forward)? it->e2 : it->s2;
			pidentity1 = it->pidentity;
			alignment_size = e1-s1+1;
		} 
		else if(cs2 == 0) {
			edge_local->align(seq2_obj->seq().reverse_complement(), (it->dir == Forward)? seq1_obj->seq() : seq1_obj->seq().reverse_complement(), score);
			pidentity1 = 100.0*edge_local->num_identical()/edge_local->alignment_size();
			s2 = seq2_obj->seq().size()-edge_local->seq1_start();
			e2 = seq2_obj->seq().size()-edge_local->seq1_end();
			s1 = (it->dir == Forward)? edge_local->seq2_start()+1 : seq1_obj->seq().size()-edge_local->seq2_start();
			e1 = (it->dir == Forward)? edge_local->seq2_end()+1 : seq1_obj->seq().size()-edge_local->seq2_end();
			alignment_size = edge_local->alignment_size();
		}
		if((pidentity1 >= pidentity_threshold) && (alignment_size >= overlap_lenght_threshold)) {
			// CONNECTED	mol-32-15fa-013605	3	mol-32-15fa-045455	5	7297	8388	1	1092	8388	8885	100
			int side2 = (it->dir == Forward)? 5 : 3;
			
			cout 	<< "CONNECTED" << '\t' << seq2_obj->display_id() << "\t5\t" << seq1_obj->display_id() << '\t' << "Middle" << '\t'
				<< s2 << '\t' << e2 << '\t' << s1 << '\t' << e1 << '\t' << seq2_obj->seq().size() << '\t' << seq1_obj->seq().size() 
				<< '\t' << pidentity1 << endl;
			continue;
		}

		// Finally - if all failed then maybe this is a shared region?
		local->align(seq1_obj->seq().subseq(cs1, ce1), (it->dir == Forward)? seq2_obj->seq().subseq(cs2, ce2) : seq2_obj->seq().subseq(cs2, ce2).reverse_complement(), score);
		pidentity1 = 100.0*local->num_identical()/local->alignment_size();
		if((pidentity1 >= pidentity_threshold) && (local->alignment_size() >= overlap_lenght_threshold)) {
			s1 = cs1+local->seq1_start()+1;
			e1 = cs1+local->seq1_end()+1;
			s2 = (it->dir == Forward)? (cs2+local->seq2_start()+1) : (seq2_obj->seq().size()-((seq2_obj->seq().size()-ce2)+local->seq2_start()));
			e2 = (it->dir == Forward)? (cs2+local->seq2_end()+1) : (seq2_obj->seq().size()-((seq2_obj->seq().size()-ce2)+local->seq2_end()));

			cout 	<< "SHARED" << '\t' << seq1_obj->display_id() << "\tMiddle\t" << seq2_obj->display_id() << '\t' << "Middle" << '\t'
				<< s1 << '\t' << e1 << '\t' << s2 << '\t' << e2 << '\t' << seq1_obj->seq().size() << '\t' << seq2_obj->seq().size() 
				<< '\t' << pidentity1 << endl;
			continue;
		}
	}

	cerr << time_since_epoch() << "Finished successfully" << endl << endl;
	delete(l1g2);
	delete(edges);
	delete(edge_local);
	delete(local);
	return 0;
}

/***********************************************************************************************************************************/
bool read_params(int argc, const char* argv[])
{
	if(argc < 3) {
		usage(argv[0]);
		return false;
	}
	m8_blast_file = argv[--argc];
	seq_file = argv[--argc];
	if(!file_exists(seq_file)) {
		cerr << "Error: file " << seq_file << " does not exist" << endl;
		return false;
	}
	if(!file_exists(m8_blast_file)) {
		cerr << "Error: file " << m8_blast_file << " does not exist" << endl;
		return false;
	}
	for(int i=1; i<argc; i++) {
		if(!strcmp(argv[i], "-d")) {
			if(i == (argc-1)) {
				cerr << "Error: seed length threshold (-d) must be specified" << endl << endl;
				return false; 
			}
			if(atoi(argv[++i]) <= 0) {
				cerr << "Error: seed length threshold (-d) must be a positive integer" << endl;
				return false; 
			}
			seed_lenght_threshold = atoi(argv[i]);
		}
		else if(!strcmp(argv[i], "-G")) {
			if(i == (argc-1)) {
				cerr << "Error: gap cost (-G) must be specified" << endl << endl;
				return false; 
			}
			if(atof(argv[++i]) >= 0) {
				cerr << "Error: gap cost (-G) must be a negative integer" << endl;
				return false; 
			}
			score.set_gap(atof(argv[i]));
		}
		else if(!strcmp(argv[i], "-os")) {
			if(i == (argc-1)) {
				cerr << "Error: overlap threshold (-os) must be specified" << endl << endl;
				return false; 
			}
			if(atoi(argv[++i]) <= 0) {
				cerr << "Error: overlap threshold (-os) must be a positive integer" << endl;
				return false; 
			}
			overlap_lenght_threshold = atoi(argv[i]);
		}
		else if(!strcmp(argv[i], "-p")) {
			if(i == (argc-1)) {
				cerr << "Error: \% identity threshold (-p) must be specified" << endl << endl;
				return false; 
			}
			++i;
			if((atof(argv[i]) < 0) || (atof(argv[i]) > 100)) {
				cerr << "Error: \% identity threshold (-p) must be in the range [0, 100]" << endl;
				return false; 
			}
			pidentity_threshold = atof(argv[i]);
		}
		else if(!strcmp(argv[i], "-q")) {
			if(i == (argc-1)) {
				cerr << "Error: mismatch penalty (-q) must be specified" << endl << endl;
				return false; 
			}
			if(atof(argv[++i]) >= 0) {
				cerr << "Error: mismatch penalty (-q) must be a negative integer" << endl;
				return false; 
			}
			score.set_mismatch(atof(argv[i]));
		}
		else if(!strcmp(argv[i], "-r")) {
			if(i == (argc-1)) {
				cerr << "Error: match reward (-r) must be specified" << endl << endl;
				return false; 
			}
			if(atof(argv[++i]) <= 0) {
				cerr << "Error: match reward (-r) must be a positive integer" << endl;
				return false; 
			}
			score.set_match(atof(argv[i]));
		}
		else {
			cerr << "Unknown option, " << argv[i] << endl << endl;
			usage(argv[0]);
			return false;
		}
	}
	return true;
}

/***********************************************************************************************************************************/
void usage(const char* prog_name)
{
	cerr << endl << "Usage: overlap [-os <min-overlap-size>] [-G <gap-penalty>] [-r <match-reward>] [-d <seed-dize>]" << endl;
	cerr <<         "               [-q <mismatch-penalty>] [-p <\% identity>] <seq-file> <m8-blast-file>" << endl << endl;
	cerr << "Where" << endl;
	cerr << " -os      | minimum overlap required for determining connection (default: 500 bp)" << endl;
	cerr << " -d       | minimum alignment size required for a couple to be considered in the m8 blast file (default: 300)" << endl;
	cerr << " -p       | \% identity threshold for a connection (default: 99%)" << endl;
	cerr << " -G       | Cost to open a gap (defualt: -5)" << endl;
	cerr << " -r       | Reward for a nucleotide match (default: 1)" << endl;
	cerr << " -q       | Penalty for a nucleotide mismatch (default: -3)" << endl;
	cerr << " m8-blast-file is a self-blast report (recommended aruments: -F F -m 8 -r 1 -q -3)" << endl;
	cerr << " seq-file is a FASTA file of the analyzed sequences" << endl; 
	cerr << endl;
}
