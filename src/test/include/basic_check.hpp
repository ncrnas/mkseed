#include <seqan/seq_io.h>

int gtest_compare_two_files(seqan::CharString const &f1, seqan::CharString const &f2);
int gtest_read_file(seqan::CharString const &fname, seqan::StringSet<seqan::CharString> &lines);