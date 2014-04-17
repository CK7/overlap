/*
 * Alignment.h
 *
 *  Created on: 05/Mar/2013
 *      Author: itai Sharon, itai.sharon@gmail.com
 */

#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

namespace Bio {

/*****************************************************************************************************************************************
 * class Alignment
 * Base class for alignment algorithms such as Needleman-Wunsch or Smith-Waterman. T parameter is the type of string (DNAString, Protein).
 *****************************************************************************************************************************************/
template <class T>
class Alignment {
public:
	class Scoring {
	public:
		Scoring()		{}
		virtual double score(char ch1, char ch2) const = 0;
	};
public:
	Alignment() : scoring(NULL), matrix(0), matrix_ncells(0), matrix_nrows(0), startr(0), endr(0), startc(0), endc(0), nidentical(0), 
		npositives(0), nmismatches(0), ngaps(0), alignment_score(0)	{}
	virtual ~Alignment()						{if(matrix) {free(matrix[0]); free(matrix);}}
	virtual Alignment<T>*	duplicate() const = 0;
	void			align(const T& seq1, const T& seq2, const Scoring& scoring);
	const string&		seq1_alignment_line() const		{return row_seq_line;}
	const string&		seq2_alignment_line() const		{return col_seq_line;}
	const string&		consensus_alignment_line() const	{return consensus_line;}
	size_t			seq1_start() const			{return startr;}
	size_t			seq1_end() const			{return endr;}
	size_t			seq2_start() const			{return startc;}
	size_t			seq2_end() const			{return endc;}
	size_t			num_identical() const			{return nidentical;}
	size_t			num_positives() const			{return npositives;}
	size_t			num_mismatches() const			{return nmismatches;}
	size_t			num_gaps() const			{return ngaps;}
	double			score() const				{return alignment_score;}
	size_t			alignment_size() const			{return consensus_line.size();}
protected:
	virtual	void		set_alignment() = 0;
	virtual void		prepare_matrix() = 0;
	void 			fill_matrix_entry(size_t r, size_t c);
	bool 			elongate_alignment();
protected:
	typedef enum {None, Row, Column, Diagonal} Direction;
	struct MatrixEntry {
		MatrixEntry() : dir(None), score(-1)	{}
		Direction	dir;
		double		score;
	};
protected:
	T			row_seq, col_seq;
	const Scoring*		scoring;
	MatrixEntry**		matrix;
	size_t			matrix_ncells, matrix_nrows;
	size_t			startr, endr;
	size_t			startc, endc;
	string			row_seq_line;
	string			col_seq_line;
	string			consensus_line;
	size_t			nidentical;
	size_t			npositives;
	size_t			nmismatches;
	size_t			ngaps;
	double			alignment_score;
};

/*****************************************************************************************************************************************/
template <class T>
void Alignment<T>::align(const T& seq1, const T& seq2, const Scoring& _scoring)
{
	row_seq = seq1;
	col_seq = seq2;
	scoring = &_scoring;
	startr = 0; 
	endr = 0;
	startc = 0;
	endc = 0;
	row_seq_line = "";
	col_seq_line = "";
	consensus_line = "";
	nmismatches = npositives = nidentical = 0;
	alignment_score = 0;

	// The following code is supposed to save time when this is ran many times. Instead of allocating and freeing memory again and again this will
	// just keep the largest chunk of memory allocated
	if(matrix) {
		if(row_seq.size()+1 > matrix_nrows) {
			MatrixEntry* p = matrix[0];
			free(matrix);
			matrix = ((MatrixEntry**)calloc(row_seq.size()+1, sizeof(MatrixEntry*)));
			matrix[0] = p;
			matrix_nrows = row_seq.size()+1;		
		}
		if((row_seq.size()+1)*(col_seq.size()+1) > matrix_ncells) {
			free(matrix[0]);
			matrix[0] = (MatrixEntry*)calloc((row_seq.size()+1)*(col_seq.size()+1), sizeof(MatrixEntry));
			matrix_ncells = (row_seq.size()+1)*(col_seq.size()+1);
		}
	}
	else {
		matrix = ((MatrixEntry**)calloc(row_seq.size()+1, sizeof(MatrixEntry*)));
		matrix[0] = (MatrixEntry*)calloc((row_seq.size()+1)*(col_seq.size()+1), sizeof(MatrixEntry));
		matrix_nrows = row_seq.size()+1;
		matrix_ncells = (row_seq.size()+1)*(col_seq.size()+1);
	}


	for(size_t i=1; i<=row_seq.size(); i++)
		matrix[i] = matrix[0]+i*(col_seq.size()+1);

	this->prepare_matrix();
	// set_alignment should set the alignment strings as well as start and end coordinated of both sequences
	this->set_alignment();

        this->nidentical = this->npositives = this->nmismatches = 0;
        for(string::const_iterator it = this->consensus_line.begin(); it!=this->consensus_line.end(); it++) {
                if(*it == ' ')
                        this->nmismatches++;
                else if(*it == '+')
                        this->npositives++;
                else
                        this->nidentical++;
        }
        this->npositives += this->nidentical;

        this->ngaps = 0;
        for(string::const_iterator it = this->row_seq_line.begin(); it!=this->row_seq_line.end(); it++)
                this->ngaps += (*it == '-');
        for(string::const_iterator it = this->col_seq_line.begin(); it!=this->col_seq_line.end(); it++)
                this->ngaps += (*it == '-');
}

/*****************************************************************************************************************************************/
template <class T>
void Alignment<T>::fill_matrix_entry(size_t r, size_t c)
{
	MatrixEntry& entry = matrix[r][c];
	if((r == 0) && (c == 0)) {
		entry.score = 0;
		entry.dir = None;
	}
	else if(c == 0) {
		entry.score = matrix[r-1][c].score+scoring->score(row_seq.get(r-1), '-');
		entry.dir = Row;
	}
	else if(r == 0) {
		entry.score = matrix[r][c-1].score+scoring->score('-', col_seq.get(c-1));
		entry.dir = Column;
	}
	else {
		entry.score = matrix[r-1][c-1].score+scoring->score(row_seq.get(r-1), col_seq.get(c-1));
		entry.dir = Diagonal;
		double s = matrix[r-1][c].score+scoring->score(row_seq.get(r-1), '-');
		if(s > entry.score) {
			entry.score = s;
			entry.dir = Row;
		}
		s = matrix[r][c-1].score+scoring->score('-', col_seq.get(c-1));
		if(s > entry.score) {
			entry.score = s;
			entry.dir = Column;
		}
	}
}

/*****************************************************************************************************************************************/
template <class T>
bool Alignment<T>::elongate_alignment()
{
	if((startr == 0) && (startc == 0) || (matrix[startr][startc].dir == None))
		return false;

	if(matrix[startr][startc].dir == Alignment::Diagonal) {
		char	chc = col_seq.get(startc-1), chr = row_seq.get(startr-1);
		row_seq_line.insert(0, 1, chr);		
		col_seq_line.insert(0, 1, chc);
		if(chc==chr) {
			consensus_line.insert(0, 1, '|');
			nidentical++;
		}
		else if(scoring->score(chc, chr) > 0)
			consensus_line.insert(0, 1, '+');
		else
			consensus_line.insert(0, 1, ' ');
		startr--;
		startc--;
	}
	else if(matrix[startr][startc].dir == Alignment::Row) {
		row_seq_line.insert(0, 1, row_seq.get(startr-1));		
		col_seq_line.insert(0, 1, '-');
		consensus_line.insert(0, 1, ' ');
		startr--;
	}
	else if(matrix[startr][startc].dir == Alignment::Column) {
		row_seq_line.insert(0, 1, '-');
		col_seq_line.insert(0, 1, col_seq.get(startc-1));
		consensus_line.insert(0, 1, ' ');
		startc--;
	}
	return true;
}

}
#endif // ALIGNMENT_H
