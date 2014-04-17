/*
 * End1_Local2.h
 *
 *  Created on: 11/Oct/2013
 *      Author: itai Sharon, itai.sharon@gmail.com
 */

#ifndef END1_LOCAL2_H
#define END1_LOCAL2_H

#include "Alignment.h"

#include <iostream>
#include <stdlib.h>

using namespace std;

namespace Bio {

/*****************************************************************************************************************************************
 * class End1_Local2
 * Looks for the best alignment of the first k or more bps of seq2 with the last k or more bps of seq1
 *****************************************************************************************************************************************/
template <class T>
class End1_Local2 : public Alignment<T> {
public:
	End1_Local2(size_t min_ovlp=23) : Alignment<T>(), min_overlap(min_ovlp)	{}
	End1_Local2<T>*	duplicate() const					{return new End1_Local2(*this);}
protected:
        void		set_alignment();
	void		prepare_matrix();
	size_t		min_overlap;
};

/*****************************************************************************************************************************************/
template <class T>
void End1_Local2<T>::set_alignment()
{
	this->endc = this->col_seq.size();
	this->endr = this->row_seq.size();
	this->alignment_score = this->matrix[this->endr][this->endc].score;
	// The path must end at the last position of row_seq (row_seq.size()) and not before column min_overlap 
	for(size_t c = this->col_seq.size(); c >= min_overlap; c--) {
		if(this->matrix[this->endr][c].score > this->alignment_score) {
			this->alignment_score = this->matrix[this->endr][c].score;
			this->endc = c;
		}
	}

	this->startr = this->endr;
	this->startc = this->endc;

	this->row_seq_line.resize(0);
	this->col_seq_line.resize(0);
	this->consensus_line.resize(0);

	while(this->elongate_alignment())
		;

        this->endr--;
        this->endc--;
}

/*****************************************************************************************************************************************/
template <class T>
void End1_Local2<T>::prepare_matrix()
{
        for(size_t r=0; r<=this->row_seq.size(); r++)
        {
		for(size_t c=0; c<=this->col_seq.size(); c++)
                {
			// We want the path to begin at position 0 to col_seq, and not after position row_seq.size()-min_overlap of row_seq  
			if((c == 0) && (r > this->row_seq.size() - min_overlap)) 
			{
				// Setting the score in these entries to something that will never lead to a path starting from them
				this->matrix[r][c].score = -((this->col_seq.size()+this->row_seq.size())*100000.0);
				this->matrix[r][c].dir = Alignment<T>::None;
			}
			else {
				this->fill_matrix_entry(r, c);
				if(this->matrix[r][c].score <= 0) {
					this->matrix[r][c].score = 0;
					this->matrix[r][c].dir = Alignment<T>::None;
				}
			}
                }
        }
}

}

#endif // END1_LOCAL2_H
