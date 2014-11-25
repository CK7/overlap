/*
 * Needleman_Wunsch.h
 *
 *  Created on: 07/Mar/2013
 *      Author: itai Sharon, itai.sharon@gmail.com
 */

#ifndef LOCAL1_GLOBAL2_H
#define LOCAL1_GLOBAL2_H

#include "Alignment.h"

namespace Bio {

/*****************************************************************************************************************************************
 * class Local1_Global2
 * Aligns seq2 globally and seq1 locally
 *****************************************************************************************************************************************/
template <class T>
class Local1_Global2 : public Alignment<T> {
public:
	Local1_Global2() : Alignment<T>()	{}
	Local1_Global2<T>*	duplicate() const			{return new Local1_Global2<T>(*this);}
protected:
        void			set_alignment();
	void			prepare_matrix();
};

/*****************************************************************************************************************************************/
template <class T>
void Local1_Global2<T>::set_alignment()
{
	size_t row_seq_size = this->row_seq.size();

	this->endc = this->col_seq.size();
	this->endr = this->row_seq.size();
	this->alignment_score = this->matrix[this->endr][this->endc].score;
	for(size_t r=0; r<=row_seq_size; r++) {
		if(this->matrix[r][this->endc].score > this->alignment_score) {
			this->alignment_score = this->matrix[r][this->endc].score;
			this->endr = r;
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
void Local1_Global2<T>::prepare_matrix()
{
	size_t row_seq_size = this->row_seq.size();
	size_t col_seq_size = this->col_seq.size();
        for(size_t r=0; r<=row_seq_size; r++)
        {
		for(size_t c=0; c<=col_seq_size; c++)
                {
			this->fill_matrix_entry(r, c);
			if((c == 0) && (this->matrix[r][c].score < 0)) {
				this->matrix[r][c].score = 0;
				this->matrix[r][c].dir = Alignment<T>::None;
			}
                }
        }
}

}

#endif // LOCAL1_GLOBAL2_H
