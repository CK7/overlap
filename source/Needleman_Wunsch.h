/*
 * Needleman_Wunsch.h
 *
 *  Created on: 07/Mar/2013
 *      Author: itai Sharon, itai.sharon@gmail.com
 */

#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

#include "Alignment.h"

namespace Bio {

/*****************************************************************************************************************************************
 * class Needleman_Wunsch
 * Implements the global alignment algorithm  Needleman-Wunsch. T parameter is the type of string (DNAString, Protein).
 *****************************************************************************************************************************************/
template <class T>
class Needleman_Wunsch : public Alignment<T> {
public:
	Needleman_Wunsch() : Alignment<T>()			{}
	Needleman_Wunsch<T>*   duplicate() const 		{return new Needleman_Wunsch(*this);}  
protected:
        void		set_alignment();
	void		prepare_matrix();
};

/*****************************************************************************************************************************************/
template <class T>
void Needleman_Wunsch<T>::set_alignment()
{
	this->endr = this->row_seq.size();
	this->endc = this->col_seq.size();
	this->alignment_score = this->matrix[this->endr][this->endc].score;
	this->nmismatches = this->npositives = this->nidentical = 0;
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
void Needleman_Wunsch<T>::prepare_matrix()
{
        for(size_t r=0; r<=this->row_seq.size(); r++)
        {
		for(size_t c=0; c<=this->col_seq.size(); c++)
                {
			this->fill_matrix_entry(r, c);
                }
        }
}

}

#endif // NEEDLEMAN_WUNSCH_H
