/* junctionFinder.hpp
 * Declarations for functions for finding the telomere junction location in a telomeric read.
 * Author: Seth Baunach
 * Date: 6/30/2020
 *
 */

#ifndef JUNCTIONFINDER_HPP
#define JUNCTIONFINDER_HPP
#include <stdexcept>
#include <limits>
#include "core/telomere_core.hpp"

namespace core = telomere_core;

namespace junctionFinder {

static int windowLen1; // window length for first pass
static int stepLen1; // step length for first pass
static double cutoff1; // cutoff for first pass
static int shortWindowLen1; // cutoff for second pass; to match the bad part.
static double shortCutoff1; // cutoff for second pass; to match the bad part.
static int startTries1; // amount of times to try skipping when first window fails
static int startSkipLen1; // amount of bp to skip and try looking for a window if first window fails completely

/* Returns sum of insertions and deletions in first n positions of optimal WF edit path.
 *  Deletions are +1, insertions are -1.
 *  Uses traceback matrix
 */
int indelSum(vector<vector<pair<int, int> > > tb_matrix, int n);

/* Uses a sliding window approach to determine where the telomere junction ends.
 * Does two passes; one with a larger window size, one with a smaller size.
 *  0. Determine offset of initial sequence by matching first 2 repeats, aligning, picking best offset
 *  1. Run W-F algorithm on the actual sequence and expected sequence based on offset
 *		If penalty exceeds cutoff, go to 2
 *  	Otherwise, adjust offset (based on stepLen and also if there was an insertion or deletion in the step)
 *		Adjust window and repeat
 *	2. Take a window of 2 telomere repeats and read until quality drops off
 *		When quality drops off, look for next place where offset is 0 and return as junction
 *
 *	Returns 0 if no telomeric sequence was matched
 *  Returns INF ( > seq.size()) if the whole sequence was telomeric
 */
size_t findTelomereJunction(string seq, core::TelRepeatInfo telRepeatInfo,
                            size_t windowLen = windowLen1, size_t stepLen = stepLen1,
                            double cutoff = cutoff1, double shortWindowLen = shortWindowLen1,
                            double shortCutoff = shortCutoff1, int startTries = startTries1,
                            int startSkipLen = startSkipLen1);

// includes information for last window matched for first and second pass
size_t findTelomereJunction_i(string seq, core::TelRepeatInfo telRepeatInfo,
                              size_t & firstWindowStart, size_t & secondWindowStart,
                              size_t windowLen = windowLen1, size_t stepLen = stepLen1,
                              double cutoff = cutoff1, double shortWindowLen = shortWindowLen1,
                              double shortCutoff = shortCutoff1, int startTries = startTries1,
                              int startSkipLen = startSkipLen1);

} // end namespace junctionFinder

#endif // JUNCTIONFINDER_HPP
