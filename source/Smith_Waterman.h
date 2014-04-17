/*
 * Smith_Waterman.h
 *
 *  Created on: 06/Mar/2013
 *      Author: itai Sharon, itai.sharon@gmail.com
 */

#ifndef SMITH_WATERMAN_H
#define SMITH_WATERMAN_H

#include "Alignment.h"

namespace Bio {

/*****************************************************************************************************************************************
 * class Smith_Waterman
 * Base class for alignment algorithms such as Needleman-Wunsch or Smith-Waterman. T parameter is the type of string (DNAString, Protein).
 *****************************************************************************************************************************************/
template <class T>
class Smith_Waterman : public Alignment<T> {
public:
	Smith_Waterman<T>*	duplicate() const			{return new Smith_Waterman(*this);}
protected:
        void		set_alignment();
	void		prepare_matrix();
};

/*****************************************************************************************************************************************/
template <class T>
void Smith_Waterman<T>::set_alignment()
{
	this->endr = this->endc = 0;
	this->alignment_score = 0;
	for(size_t r=1; r<=this->row_seq.size(); r++) { 
		for(size_t c=1; c<=this->col_seq.size(); c++) { 
			if(this->matrix[r][c].score > this->alignment_score) {
				this->endr = r;
				this->endc = c;
				this->alignment_score = this->matrix[r][c].score;
			}
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
void Smith_Waterman<T>::prepare_matrix()
{
        for(size_t r=0; r<=this->row_seq.size(); r++)
        {
		for(size_t c=0; c<=this->col_seq.size(); c++)
                {
			this->fill_matrix_entry(r, c);
			if(this->matrix[r][c].score <= 0) {
				this->matrix[r][c].score = 0;
				this->matrix[r][c].dir = Alignment<T>::None;
			}
                }
        }
}

}

#endif // SMITH_WATERMAN_H
